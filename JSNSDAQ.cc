#include <RAT/JSNSDAQ.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DB.hh>
#include <RAT/PMTWaveform.hh>
#include <RAT/Digitizer.hh>
#include <RAT/DS/DigitizedWaveform.hh>
#include <RAT/FEE.hh>
#include <RAT/Trigger.hh>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <float.h>

namespace RAT {

JSNSDAQProc::JSNSDAQProc() : Processor("jsnsdaq") {

    DBLinkPtr ldaq = DB::Get()->GetLink("DAQ");
    double min_input = ldaq->GetD("caen_min_input");
    double range = ldaq->GetD("caen_dynamic_range");
    unsigned int nbits = ldaq->GetI("caen_nbits");
    fNSamples = ldaq->GetI("caen_nsamples");
    fSamplingTime = ldaq->GetD("caen_sampling_time");
    CAENDigitizer = Digitizer(min_input, range, nbits);    

    fthreshold = ldaq->GetD("FADC_self_trigger_threshold");
    fVetothreshold = ldaq->GetD("FADC_veto_trigger_threshold");
    ftrigsampTime = 1.0;//add to .ratdb
    fchunksize = ldaq->GetI("FEE_nchannels");
    fVetochunksize = ldaq->GetI("FEE_nVetochannels");
    fscale = ldaq->GetI("attenuator_scale");

    fPEThreshold = 10;
    ftimeBlock = 10.0;

    std::vector<double> freq_filter = ldaq->GetDArray("fee_bandwidth_freqs");
    std::vector<double> low_gain_bandwidth = ldaq->GetDArray("fee_bandwidth_lowgain");
    std::vector<double> high_gain_bandwidth = ldaq->GetDArray("fee_bandwidth_highgain");
    fee.SetBandwidthCurve(freq_filter, low_gain_bandwidth, high_gain_bandwidth);


    std::vector<double> gain_x = ldaq->GetDArray("fee_gain_curve_high_x");
    std::vector<double> gain_y = ldaq->GetDArray("fee_gain_curve_high_y");
    fee.SetGainCurve(gain_x, gain_y, FEE::HIGHGAIN);
    gain_x =  ldaq->GetDArray("fee_gain_curve_low_x");
    gain_y =  ldaq->GetDArray("fee_gain_curve_low_y");
    fee.SetGainCurve(gain_x, gain_y, FEE::LOWGAIN);

    fee.fHighGainDeltaT = ldaq->GetD("fee_high_gain_time_delay");
    fee.fLowGainDeltaT = ldaq->GetD("fee_low_gain_time_delay");
    fee.fAnalogSumDeltaT = ldaq->GetD("fee_analog_sum_time_delay");
    fee.fSaturationValue = ldaq->GetD("fee_saturation_voltage");
    fee.fHighGainNoise = ldaq->GetD("fee_high_gain_noise");
    fee.fLowGainNoise = ldaq->GetD("fee_low_gain_noise");
    fee.fAnalogSumNoise = ldaq->GetD("fee_analog_sum_noise");

    PMTPulse::Initialize();

    fEventCounter = 0;
}

Processor::Result JSNSDAQProc::DSEvent(DS::Root *ds) {
  DS::MC *mc = ds->GetMC();
  DS::Run* run = DS::RunStore::GetRun(ds->GetRunID());
  RAT::DS::PMTInfo *pmtinfo = run->GetPMTInfo();
  // there is already a EV branch present, rmove it, otherwise we'll have
  // multiple detector events. TODO, add warning output
  if(ds->ExistEV()) {
    ds->PruneEV();
  }
  // create an associative array between PMT ids and the signals they produce.
  std::map<int, PMTWaveform> pmt_wfs;

    // First cycle through all PMTs, then all hits, and
    // then all PMT responses to generate waveforms.
  for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
      DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
      int pmtid = mcpmt->GetID();
      for(int imcphoton=0; imcphoton < mcpmt->GetMCPhotonCount(); imcphoton++) {
          DS::MCPhoton* mcphoton = mcpmt->GetMCPhoton(imcphoton);
          for(int ipmtresponse=0; ipmtresponse < mcphoton->GetNPMTResponse(); ipmtresponse++) {
            double time = mcphoton->GetFrontEndTimes()[ipmtresponse];
            double charge = mcphoton->GetCharges()[ipmtresponse];

            // Produce PMT waveforms
            PMTPulse thisPulse(time, charge);
            pmt_wfs[pmtid].AddPulse(thisPulse);
          }
      }
      fee.PlugInSignal(pmtid, &pmt_wfs[pmtid]);
  }
  // Now add waveforms for the PMTs that did not get hit since they get readout too.
  // The PMTInfo class has PMTs in order of ID, i.e. the 0th pmt in the PMTInfo
  // has ID 0, the first has ID 1, etc.
  PMTPulse emptyPulse = PMTPulse::EmptyPulse();
  for(int iPMT=0; iPMT < pmtinfo->GetPMTCount(); iPMT++) {
      if(pmt_wfs.count(iPMT) > 0) {
          continue;
      }
      fee.PlugInSignal(iPMT, &emptyPulse);
  }
  
  //Initialize triggers
  double dtime = fSamplingTime*fNSamples;
  TargetTrigger = Trigger(fthreshold, fscale, fchunksize, dtime, ftrigsampTime);
  TopVetoTrigger = Trigger(fVetothreshold, fscale, fVetochunksize, dtime, ftrigsampTime);
  BotVetoTrigger = Trigger(fVetothreshold, fscale, fVetochunksize, dtime, ftrigsampTime);

  //Finds the first and last light entering the PMTs
  //Not tecknically how the trigger will work. The amount of time is reduced to 
  //save time in the computation
  double min_time = DBL_MAX;
  double max_time = DBL_MIN;

  for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) {
     DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
     for(int imcphoton=0; imcphoton < mcpmt->GetMCPhotonCount(); imcphoton++) {
	DS::MCPhoton* mcphoton = mcpmt->GetMCPhoton(imcphoton);
	double itime = mcphoton->GetFrontEndTime();

	if(min_time > itime){
	  min_time = itime-20; //-20 fudge factor probably not needed
	}

        if(max_time < itime){
          max_time = itime+20; //+20 fudge factor probably not needed
        }
     }
  }

  //Sets the computation threshold for charge by finding the PMT with the lowest value of 
  //spe_charge in the detector. I have set it to the lowest for saftey
  double chargeScale = DBL_MAX;
  int numModels = pmtinfo->GetModelCount();

  for (int i = 0; i < numModels; i++) {
    std::string modelName = pmtinfo->GetModelName(i);
    DBLinkPtr model = DB::Get()->GetLink("PMTRESPONSE", modelName);
    double spe_charge = model->GetD("spe_charge");
    if(spe_charge<chargeScale){
      chargeScale=spe_charge;
    }
  }
  fChargeThreshold = fPEThreshold*chargeScale;

  //Determines the number of target and veto PMTs and chunks
  for( int ipmt=0; ipmt<pmtinfo->GetPMTCount(); ipmt++){
    if(pmtinfo->GetType(ipmt) == 1){
      TargetTrigger.AddPMT(ipmt);
      continue;
    }
    if(pmtinfo->GetPosition(ipmt)[2]>0){
      TopVetoTrigger.AddPMT(ipmt);
      continue;
    }
    BotVetoTrigger.AddPMT(ipmt);
  }
  //Computational Trigger: not part of the actual DAQ trigger, implemented only to save computing time.
  //This trigger counts PE and computes the total charge in some time block and if either exceeds a threshold
  //the program proceeds to record samples for that time block.
  //Note: timeBlock must be bigger than trigSampTime or program will freak out.
  for(int itime=min_time; itime<max_time; itime+=ftimeBlock){
    //Counts PE and charge in this time block
    int nPE = 0;
    double charge = 0.0;
    //double pmtcharge = 0.0;
    for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++){
      DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
      for(int imcphoton=0; imcphoton < mcpmt->GetMCPhotonCount(); imcphoton++){
        DS::MCPhoton* mcphoton = mcpmt->GetMCPhoton(imcphoton);
        double jtime = mcphoton->GetFrontEndTime();
        if(itime<=jtime && jtime<=itime+ftimeBlock){
          nPE++;
          for(int ipmtresponse=0; ipmtresponse < mcphoton->GetNPMTResponse(); ipmtresponse++) {
            charge += mcphoton->GetCharges()[ipmtresponse];
          }
        }
      }
    }

    //Check if enough samples occured in this time block to record
    if(nPE<fPEThreshold && charge<fChargeThreshold){
    continue;
    }
    TargetTrigger.AddSamples(itime, ftimeBlock, &fee);
    TopVetoTrigger.AddSamples(itime, ftimeBlock, &fee);
    BotVetoTrigger.AddSamples(itime, ftimeBlock, &fee);
  }  
  
  //This returns a vector of all the times at which an event was triggered
				std::cout<<"Target Signals:  "<<std::endl;
  std::vector<double> TargetTrigTimes = TargetTrigger.GetTrigTimes();
				std::cout<<"Top Veto Signals:  "<<std::endl;
  std::vector<double> TopVetoTrigTimes = TopVetoTrigger.GetTrigTimes();
				std::cout<<"Bottom Veto Signals:  "<<std::endl;
  std::vector<double> BotVetoTrigTimes = BotVetoTrigger.GetTrigTimes();

  //Loop through all triggered events, add new EV event and the digitize the signals
  int numtriggers = TargetTrigTimes.size();
  int numTopVetotriggers = TopVetoTrigTimes.size();
  int numBotVetotriggers = BotVetoTrigTimes.size();
                                std::cout<<"Target Trigger Times:  "<<std::endl;
                                for(int i=0; i<numtriggers; i++){std::cout<<TargetTrigTimes[i]<<std::endl;}
                                std::cout<<"Top Veto Trigger Times:  "<<std::endl;
                                for(int i=0; i<numTopVetotriggers; i++){std::cout<<TopVetoTrigTimes[i]<<std::endl;}
                                std::cout<<"Bottom Veto Trigger Times:  "<<std::endl;
                                for(int i=0; i<numBotVetotriggers; i++){std::cout<<BotVetoTrigTimes[i]<<std::endl;}



  for(int i=0; i<numtriggers; i++){  

    //Adds new EV event
    DS::EV *ev = ds->AddNewEV();
    ev->SetID(fEventCounter);
    fEventCounter++;
  
    // And finally digitze all the signals
    for(int iPMT=0; iPMT < pmtinfo->GetPMTCount(); iPMT++) {

        fee.OutputSelect(iPMT, FEE::HIGHGAIN);
        std::vector<uint16_t> high_gain_digitized_wf = CAENDigitizer.Digitize(&fee,
                                                                            TargetTrigTimes.at(i),
                                                                            fNSamples,
                                                                            fSamplingTime);
        fee.OutputSelect(iPMT, FEE::LOWGAIN);
        std::vector<uint16_t> low_gain_digitized_wf = CAENDigitizer.Digitize(&fee,
                                                                           TargetTrigTimes.at(i),
                                                                           fNSamples,
                                                                           fSamplingTime);
        DS::PMT* pmt = ev->AddNewPMT();

        DS::DigitizedWaveform ds_wf;
        ds_wf.SetSamples(high_gain_digitized_wf);
        ds_wf.SetSamplingTime(fSamplingTime);
        pmt->SetHighGainWaveform(ds_wf);

        ds_wf.SetSamples(low_gain_digitized_wf);
        ds_wf.SetSamplingTime(fSamplingTime);
        pmt->SetLowGainWaveform(ds_wf);

        pmt->SetID(iPMT);
    }

  }
  return Processor::OK;
}

} // namesnssor::OK; RAT

