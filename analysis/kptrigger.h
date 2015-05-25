#ifndef __KPTRIGGER__
#define __KPTRIGGER__

#include <iostream>
#include <string>
#include <vector>

#include "RAT/DS/MC.hh"

class TRandom3;

class KPPulse {
  
 public:
  
  KPPulse();
  ~KPPulse();

  double tstart;
  double tend;
  double tpeak;
  double peakamp;

  typedef enum { kRising, kFalling, kDefined } Status;
  Status fStatus;

  double last_max;
  int nfallingbins;
  double pe;
  double z;
  int hits_assigned;
  double pe_dark;
  double pe_adjusted;
  
};

typedef std::vector< KPPulse* > KPPulseList;
typedef std::vector< KPPulse* >::iterator KPPulseListIter;

int find_trigger( RAT::DS::MC* mc, 
		  double threshold, double window_ns, double darkrate_hz,
		  bool hoop_cut, double min_hoop, double max_hoop,
		  bool time_cut, double min_time, double max_time,
		  int n_decay_constants, double decay_weights[], double decay_constants_ns[],
		  KPPulseList& pulses, int first_od_sipmid, bool veto, std::vector<double>& wfm, int version=2 );
void assign_pulse_charge( RAT::DS::MC* mc, std::string pmtinfo, KPPulseList& pulselist, 
			  double darkrate_hz,
			  bool hoop_cut, double min_hoop, double max_hoop,
			  double decay_const, int first_od_sipmid, bool veto, int version=2);

#endif
