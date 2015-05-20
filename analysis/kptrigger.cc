#include "kptrigger.h"
#include <iostream>
#include <ctime>
#include "pmtinfo.hh"

KPPulse::KPPulse() {
  nfallingbins = 0;
  last_max = 0.0;
  tpeak = tstart = tend = peakamp = 0.0;
  pe = 0.0;
  z = 0.0;
  hits_assigned = 0;
  fStatus = kRising;
}

KPPulse::~KPPulse() {}; 

int find_trigger( RAT::DS::MC* mc, 
		  double threshold, double window_ns, double darkrate_hz,
		  bool hoop_cut, double min_hoop, double max_hoop, 
		  int n_decay_constants, double decay_weights[], double decay_constants_ns[], 
		  KPPulseList& pulses, int first_od_sipmid, bool veto, std::vector<double>& hitwfm ) {

  // (1) bin hits out to 20 microseconds.
  // (2) scan until a bin over threshold
  // (3) scan until max found: using averaging (-n,+n) bins
  // (4) adjust threshold level using maxamp*exp(-t/(t0))
  hitwfm.clear();
  if ( mc->GetNumPE()==0 )
    return 0;

  // ------------------------------------------------
  // INTERNAL PARAMETERS
  const int nbins = 10000;
  const double nspertic = 1.0;
  const int NFALLING = 5;
  double tbins[nbins];
  memset( tbins, 0, sizeof(double)*nbins );
  int windowbins = (int)window_ns/nspertic;
  if ( windowbins==0 ) windowbins++;
  double darkrate_window = window_ns*(1.0e-9*darkrate_hz)*(90000); // 1 MHz dark rate * Window  * NSiPMs
  if ( hoop_cut )
    darkrate_window = window_ns*(1.0e-9*darkrate_hz)*(max_hoop-min_hoop); // 1 MHz dark rate * Window  * NSiPMs
  hitwfm.resize( nbins );

  // ------------------------------------------------
  // Fill time bins
  int npmts = mc->GetMCPMTCount();
  for ( int ipmt=0; ipmt<npmts; ipmt++ ) {

    RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
    int nhits = pmt->GetMCPhotonCount();
    int pmtid = pmt->GetID();
    int hoopid = pmtid/100;
    if ( hoop_cut && ( hoopid<min_hoop || hoopid>max_hoop ) )
      continue;

    if ( veto && pmtid<first_od_sipmid ) 
      continue;
    else if ( !veto && pmtid>=first_od_sipmid )
      continue;

    for (int ihit=0; ihit<nhits; ihit++) {
      RAT::DS::MCPhoton* hit = pmt->GetMCPhoton( ihit );
      double thit = hit->GetHitTime();
      int ibin = (int)thit/nspertic;
      if ( ibin>=0 && ibin<nbins )
	tbins[ibin]++;
    }
  }//end of pmt loop

  // ------------------------------------------------
  // Find peaks by scanning
  int npulses = 0;
  for (int ibin=windowbins; ibin<nbins; ibin++ ) {
    
    // count number active pulses
    int nactive = 0;
    for ( KPPulseListIter it=pulses.begin(); it!=pulses.end(); it++ ) {
      if ( (*it)->fStatus!=KPPulse::kDefined )
	nactive++;
    }

    // calc triggering quantities
    int nhits_window = 0;
    for (int i=ibin-windowbins; i<ibin; i++) {
      nhits_window += tbins[i];
    }

    // deprecated
//     double ave_window = 0;
//     int nave = 0;
//     int ave_start = ibin-tave_bins;
//     int ave_end   = ibin+tave_bins;
//     if ( ave_start<0 ) ave_start = 0;
//     if ( ave_end>nbins ) ave_end = nbins;
//     for (int i=ave_start; i<ave_end;  i++ ) {
//       ave_window += tbins[i];
//       nave++;
//     }
//     ave_window /= double(nave);

    // Check for new pulse
    if ( nactive==0 ) {
      // no active pulses. do search based on hits ver threshold of moving window
      if ( nhits_window>(int)threshold ) {
	// make a new pulse!
	KPPulse* apulse = new KPPulse;
	apulse->tstart = (ibin-windowbins)*nspertic;
	apulse->fStatus = KPPulse::kRising; // start of rising edge (until we find max)
	pulses.push_back( apulse );
	npulses++;
      }
    }
    else {
      // we have active pulses. we look for a second peak with a trigger algorithm that accounts for scintillator decay time

      // find modified threshold
      double modthresh = threshold;
      bool allfalling = true;
      for ( KPPulseListIter it=pulses.begin(); it!=pulses.end(); it++ ) {
	if ( (*it)->fStatus==KPPulse::kRising ) {
	  // have a rising peak. will block creation of new pulse.
	  allfalling = false;
	  modthresh += 2.0*nhits_window;
	}
	else {
	  // for pulses considered falling, we modify the threshold to be 3 sigma (roughly) above
	  //double expectation = ((*it)->peakamp)*exp( -( ibin*nspertic - (*it)->tpeak )/decay_constant ); // later can expand to multiple components
	  double arg = 0.0;
	  for (int idcy=0; idcy<n_decay_constants; idcy++) {
	    arg += decay_weights[idcy]*( ibin*nspertic - (*it)->tpeak )/decay_constants_ns[idcy];
	  }
	  double expectation = ((*it)->peakamp)*exp( -arg );
	  //double expectation = ((*it)->peakamp)*exp( -( ibin*nspertic - (*it)->tpeak )/decay_constant ); // later can expand to multiple components
	  modthresh += ( expectation + 3.0*sqrt(expectation) )*windowbins;
	}
      }//end of loop over pulses

      //std::cout << " ibin=" << ibin << " hits_window=" << nhits_window << " modthresh=" << modthresh << std::endl;

      // apply threshold
      if ( allfalling && nhits_window>modthresh ) {
	// new pulse! (on top of other pulse)
	KPPulse* apulse = new KPPulse;
        apulse->tstart = ibin*nspertic;
        apulse->fStatus = KPPulse::kRising; // start of rising edge (until we find max) 
	pulses.push_back( apulse );
        npulses++;
	//std::cout << "  new pulse!" << std::endl;
      }

      // now we find max of rising pulses and end of falling pulses
      for ( KPPulseListIter it=pulses.begin(); it!=pulses.end(); it++ ) {
	KPPulse* apulse = *it;
	if ( apulse->fStatus==KPPulse::kDefined )
	  continue;

	if ( apulse->fStatus==KPPulse::kRising ) {
	  if ( (nhits_window/double(windowbins)) < apulse->last_max )
	    apulse->nfallingbins += 1;
	  else {
	    apulse->nfallingbins = 0;
	    apulse->last_max = nhits_window/double(windowbins);
	  }

	  if ( apulse->nfallingbins>NFALLING ) {
	    // found our max
	    apulse->fStatus = KPPulse::kFalling;
	    apulse->tpeak = ibin*nspertic - apulse->nfallingbins;
	    apulse->peakamp = apulse->last_max;
	  }
	}//end of if rising
	else if ( apulse->fStatus==KPPulse::kFalling ) {
	  double decay_constant = 0.0;
	  for (int idcy=0; idcy<n_decay_constants; idcy++)
	    decay_constant += decay_weights[idcy]*decay_constants_ns[idcy];
	  if ( ibin*nspertic > apulse->tpeak + 8*decay_constant ) {
	    apulse->tend = ibin*nspertic;
	    apulse->fStatus=KPPulse::kDefined;
	  }
	}
      }//end of loop over pulses
    }//end of active pulse condition
    
  }//end of scan over timing histogram

  // save hit histogram
  for (int ibin=0; ibin<nbins; ibin++)
    hitwfm.at(ibin) = tbins[ibin];

  return npulses;
}

void assign_pulse_charge( RAT::DS::MC* mc, std::string pmtinfofile, KPPulseList& pulselist, 
			  double darkrate_hz,
			  bool hoop_cut, double min_hoop, double max_hoop,
			  double decay_const_ns, 
			  int first_od_sipmid, bool veto ) {
  // get pulse info: assigned charge
  if ( pulselist.size()==0 )
    return;

  double nsipms = (double)first_od_sipmid;
  if ( hoop_cut )
    nsipms = (double)(max_hoop-min_hoop)*100.0; // 1 MHz dark rate * Window  * NSiPMs

  PMTinfo* pmtinfo = PMTinfo::GetPMTinfo( pmtinfofile );

  const int npulses = pulselist.size();
  double heights[ npulses ];
  double height_tot;

  // clear pe, z
  for ( int ipulse=0; ipulse<npulses; ipulse++ ) {
    KPPulse* apulse = pulselist.at(ipulse);
    apulse->pe = 0.0;
    apulse->z  = 0.0;
    apulse->hits_assigned = 0;
  }
  
  // Fill time bins
  int npmts = mc->GetMCPMTCount();

  for ( int ipmt=0; ipmt<npmts; ipmt++ ) {
    RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
    int nhits = pmt->GetMCPhotonCount();
    int pmtid = pmt->GetID();
    int hoopid = pmtid/100;
    if ( hoop_cut && ( hoopid<min_hoop || hoopid>max_hoop ) )
      continue;    
    if ( veto && pmtid<first_od_sipmid )
      continue;
    else if ( !veto && pmtid>=first_od_sipmid )
      continue;

    float pmtpos[3];
    pmtinfo->getposition( pmtid, pmtpos );

    for (int ihit=0; ihit<nhits; ihit++) {
      RAT::DS::MCPhoton* hit = pmt->GetMCPhoton( ihit );
      double thit = hit->GetHitTime();

      height_tot = 0;
      memset( heights, 0, sizeof(double)*npulses );

      // make claims on hits based on pulse height
      for ( int ipulse=0; ipulse<npulses; ipulse++ ) {
	KPPulse* apulse = pulselist.at(ipulse);
	double ndarkrate = (darkrate_hz*1.0e-9)*nsipms*( apulse->tend-apulse->tstart );
	apulse->pe_dark = ndarkrate;
	
	if ( apulse->tstart<= thit && thit<= apulse->tend ) {
	  double pheight;
	  if ( thit<apulse->tpeak ) {
	    pheight = (apulse->peakamp)/( apulse->tpeak-apulse->tstart )*( thit-apulse->tstart );
	  }
	  else {
	    double dt = thit-apulse->tpeak;
	    pheight = apulse->peakamp*exp( -dt/decay_const_ns );
	  }
	  heights[ipulse] = pheight;
	  height_tot += pheight;
	}
	else {
	  heights[ipulse] = 0.0;
	}
      }//end of kppulse list
      
      // assign
      if ( height_tot>0 ) {
	for ( int ipulse=0; ipulse<npulses; ipulse++ ) {
	  KPPulse* apulse = pulselist.at(ipulse);
	  apulse->pe += heights[ipulse]/height_tot;
	  if (  heights[ipulse]>0.0 ) {
	    apulse->hits_assigned++;
	    apulse->z += (double)pmtpos[2];
// 	    if ( ipulse==0 ) {
// 	      std::cout << " pulse 1: pmtid=" << pmtid << "  pe=" << heights[ipulse] << ", z=" << pmtpos[2] << std::endl;
// 	    }
	  }
	}
      }
    }
  }//end of pmt loop

  // take care of averages and subtract dark noise
  for ( int ipulse=0; ipulse<npulses; ipulse++ ) {
    KPPulse* apulse = pulselist.at(ipulse);
    if ( apulse->hits_assigned>0 ) 
      apulse->z /= double(apulse->hits_assigned);
    apulse->pe_adjusted = apulse->pe - apulse->pe_dark;
  }

}
