#ifndef __RAT_Trigger__
#define __RAT_Trigger__

#include<vector>
#include<map>

namespace RAT {

  class AnalogSignal;
  class FEE;

class Trigger {
public:
  Trigger(double _threshold=0, double _scale=0, int _chunksize=0, double _dtime=0, double sampgapTime=0);
  void AddPMT(int ipmt);
  void AddPMTs(std::vector<int> ipmts);
  void AddSample(int chunk, double time, double signal);
  void AddSamples(double time, double timeblock, FEE* fee);
  std::vector<double> GetTrigTimes();

protected:
  std::map<double, std::map<int,double>> fallSamples; //stores samples with corresponding pmt chuck
  double fthreshold; //trigger threshold in mV
  double fscale; // amount the attenuator in the trigger scales the samples by
  int fchunksize; //Number of PMTs read at a time
  double fdtime;  // amount of time it take the digitizer to run. Serves as the gap in time after the trigger is tripped before it reads again
  double fsampgapTime; //gap between each sample the trigger reads
  std::vector<int> fPMTs; //Stores the pmts that are being read by this trigger
  int fnumPMTs; //Total number of PMTs being read by the trigger
};

}

#endif

