#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "kptrigger.h"

// int find_trigger4( std::vector<double>& tbins,
//                    double window_ns, double sigfactor, double min_pulse_width,
//                    double orig_window_ns, double darkrate_hz, int nhoops,
//                    int n_decay_constants, double decay_weights[], double decay_constants_ns[],
//                    KPPulseList& pulses );

int main(int narg, char** argv) {

  int nevents = 250000;
  const int NBINS = 10000;

  KPPulseList results;
  double window_ns_veto = 20;
  double ODsigma_threshold = 3.0;
  double min_pulse_width = 25.0;
  double sipm_darkrate_hz = 1.6e6;
  int nsipms = 100;

  int n_decay_constants_veto = 3;
  double decay_weights_veto[3] = { 0.5, 0.2, 0.3 };
  double decay_constants_ns_veto[3] = { 30.0, 90.0, 400.0 };

  KPPulseList pulselist_veto;

  //int trig_version = 4;
  //int nod_sipms_per_hoop_endcap = 100;

  // ------------------------------------------------------------
  double mean_events_fullwindow_per_sipm = sipm_darkrate_hz*1.0e-9*double(NBINS);
  
  // ------------------------------------------------------------

  TFile* out = new TFile("output_darknoise_study.root", "RECREATE");
  TTree* tree = new TTree("vetodark","Dark Rate study");

  TRandom3 rand(time(NULL));
  std::vector<double> wfm(NBINS,0.0);

  double tstart1;
  tree->Branch("tstart1", &tstart1,"tstart1/D");
  //tree->Branch("wfm",&wfm);

  int nevents_wpulse = 0;
  int nevents_wpulse_125ns = 0;
  
  for (int ievent=0; ievent<nevents; ievent++) {
    // zero vector
    std::fill( wfm.begin(), wfm.end(), 0.0 );

    // generate noise
    int tothits = 0;
    for (int i=0; i<nsipms; i++) {
      int ndarkhits = rand.Poisson( mean_events_fullwindow_per_sipm );
      for (int n=0; n<ndarkhits; n++) {
	double t = rand.Uniform()*double(NBINS);
	wfm.at( (int)t ) += 1.0;
	tothits++;
      }
    }

    if ( ievent%1000==0 ) {
      std::cout << "Event " << ievent << ": " << tothits << std::endl;
    }



    // run trigger
    int npulses_veto = find_trigger4( wfm,
				      0.5*window_ns_veto, ODsigma_threshold, min_pulse_width,
				      window_ns_veto, sipm_darkrate_hz, nsipms,
				      n_decay_constants_veto, decay_weights_veto, decay_constants_ns_veto,
				      pulselist_veto );
    tstart1 = -1;
    if (npulses_veto>0 && pulselist_veto.at(0)->tstart<NBINS-0.5*window_ns_veto) {
      nevents_wpulse++;
      tstart1 = pulselist_veto.at(0)->tstart;
      
      for ( KPPulseListIter it=pulselist_veto.begin(); it!=pulselist_veto.end(); it++ ) {
	if ( (*it)->tstart<250.0 ) { // ns
	  nevents_wpulse_125ns +=1;
	  break;
	}
      }
    }

    tree->Fill();
    free_pulse_list( pulselist_veto );
  }//event loop

  out->Write();
  std::cout << "Dark rate: " << nevents_wpulse << " (" << nevents_wpulse_125ns << ", 125 ns window) of " << nevents << std::endl;

}
