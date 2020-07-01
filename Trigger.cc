#include <RAT/Trigger.hh>
#include <vector>
#include <map>
#include <iostream>

#include <RAT/AnalogSignal.hh>
#include <RAT/FEE.hh>

namespace RAT{
  
  Trigger::Trigger(double _threshold, double _scale, int _chunksize, double _dtime, double _sampgapTime){ 
    fthreshold = _threshold;
    fscale = _scale;
    fchunksize = _chunksize;
    fdtime = _dtime;
    fsampgapTime = _sampgapTime;
  }

  void Trigger::AddSample(int chunk, double time, double signal){
    fallSamples[time][chunk] = signal;
  }

  void Trigger::AddPMT(int ipmt){
    fPMTs.push_back(ipmt);
    fnumPMTs = fPMTs.size();
  } 
  
  void Trigger::AddPMTs(std::vector<int> ipmts){
    fPMTs.insert(fPMTs.end(), ipmts.begin(), ipmts.end());
    fnumPMTs = fPMTs.size();
  }


  void Trigger::AddSamples(double time, double timeblock, FEE* fee){
    int numsamp = timeblock/fsampgapTime;
    int nchunks = fnumPMTs/fchunksize;
    int nremPMTs = fnumPMTs%fchunksize;
    //adds PMTS that are part of a complete chunk of size fchunksize
    for(int chunk=0; chunk<nchunks; chunk++){
      std::vector<int> ichannels;
      for(int j=0; j<fchunksize; j++){
        ichannels.push_back(fPMTs[j + (chunk*fchunksize)]);
      }
      fee->OutputSelect(ichannels);
      std::vector<double> sample = fee->Sample(time, fsampgapTime, numsamp);
      for(int i=0; i<numsamp; i++){
        double time_of_sample = time+(i*fsampgapTime);
        fallSamples[time_of_sample][chunk] = sample[i];
      }
    }

    //Adds any PMTs that were not part of a complete set size fchunksize
    std::vector<int> ichannels;
    for(int j=0; j<nremPMTs; j++){
      ichannels.push_back(fPMTs[j + (nchunks*fchunksize)]);
    }
    fee->OutputSelect(ichannels);
    std::vector<double> sample = fee->Sample(time, fsampgapTime, numsamp);
    for(int i=0; i<numsamp; i++){
      double time_of_sample = time+(i*fsampgapTime);
      fallSamples[time_of_sample][nchunks] = sample[i];
    }
  }


  std::vector<double> Trigger::GetTrigTimes(){
    std::vector<double> trigtimes;
				double max_signal=0; 
    //cycles through all times fed to the trigger
    auto itime = fallSamples.begin();
    while(itime != fallSamples.end()){
      double signal = 0.0;
      //looks at each seperate chunk and adds them to the current signal. Also divides by the attenuator scale
      for(auto ichunk = itime->second.begin(); ichunk != itime->second.end(); ichunk++){
	signal += (ichunk->second)/fscale;
      }
				//std::cout<<signal<<std::endl;
				if(signal<max_signal){max_signal=signal;}
      if(signal>fthreshold){
        itime++;
        continue;
      }

      trigtimes.push_back(itime->first);
      
      //Finds next time past the digitalizer window
      auto jtime = itime;
      do{
        jtime++;
        if(jtime == fallSamples.end()){
          break;
        }
      }while(itime->first + fdtime > jtime->first);
      itime=jtime;
    }
				std::cout<<max_signal<<std::endl;
  return trigtimes;
  }

}
