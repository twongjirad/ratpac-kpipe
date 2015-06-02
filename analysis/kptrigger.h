#ifndef __KPTRIGGER__
#define __KPTRIGGER__

#include <iostream>
#include <string>
#include <vector>

#include "RAT/DS/MC.hh"
#include "kpdaq.h"

class TRandom3;

class KPPulse {
  
 public:
  
  KPPulse();
  ~KPPulse();

  double tstart;
  double tend;
  double tpeak;
  double peakamp;

  typedef enum { kRising, kFalling, kDefined, kRejected } Status;
  Status fStatus;

  double last_max;
  int nfallingbins;
  double pe;
  double petrig;
  double z;
  int ihoop;
  int hits_assigned;
  double pe_dark;
  double pe_adjusted;
};

typedef std::vector< KPPulse* > KPPulseList;
typedef std::vector< KPPulse* >::iterator KPPulseListIter;
void free_pulse_list( KPPulseList& );
bool is_pulse_in_list( KPPulseList&, KPPulse& );

int find_trigger( RAT::DS::MC* mc, 
		  double threshold, double window_ns, double darkrate_hz,
		  bool hoop_cut, double min_hoop, double max_hoop,
		  bool time_cut, double min_time, double max_time,
		  int n_decay_constants, double decay_weights[], double decay_constants_ns[],
		  KPPulseList& pulses, int first_od_sipmid, bool veto, std::vector<double>& wfm, int version );
int find_trigger2( const KPDAQ& daq, 
		   double threshold, double window_ns, double darkrate_hz,
		   int chstart, int chend,
		   bool time_cut, double min_time, double max_time,
		   int n_decay_constants, double decay_weights[], double decay_constants_ns[],
		   KPPulseList& pulses, int first_od_sipmid, bool veto, std::vector<double>& wfm, int version );
int find_trigger3( std::vector<double>& tbins,
		   double threshold, double window_ns, double darkrate_hz,
		   bool time_cut, double min_time, double max_time,
		   int n_decay_constants, double decay_weights[], double decay_constants_ns[], 
		   KPPulseList& pulses, int first_od_sipmid, bool veto, int version );

int find_trigger4( std::vector<double>& tbins,
		   double window_ns, double sigfactor, double min_pulse_width,
		   double orig_window_ns, double darkrate_hz, int nhoops,
		   int n_decay_constants, double decay_weights[], double decay_constants_ns[], 
		   KPPulseList& pulses );


void assign_pulse_charge( RAT::DS::MC* mc, std::string pmtinfo, KPPulseList& pulselist, 
			  double darkrate_hz,
			  bool hoop_cut, double min_hoop, double max_hoop,
			  double decay_const, int first_od_sipmid, int nod_sipms_per_hoop_side, int nod_sipms_per_hoop_endcap,
			  bool veto, int version);

#endif
