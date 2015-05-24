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
		  bool time_cut, double min_time, double max_time,
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
  const int NFALLING = 3;
  bool use_ave = true;
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

    // VETO OR INNER
    if ( veto && pmtid<first_od_sipmid ) 
      continue;
    else if ( !veto && pmtid>=first_od_sipmid )
      continue;

    int hoopid = pmtid/100;
    if ( !veto ) {
      // inner pipe hoop structure
    }
    else {
      // version 1
      // veto has string structure, which in retrospect, was dumb, we're interested in z
      // this uses info only in my head. correlates with build_geometry.py
//       if ( pmtid<first_od_sipmid+1000 )
// 	hoopid = 900+pmtid%100; // 0-999 inclusive are radial
//       else {
// 	if (pmtid<first_od_sipmid+1100)
// 	  hoopid = 1000;     // endcap
// 	else
// 	  hoopid = 1001;     // endcap
//       }

      //version 2, this is fixed
      if ( pmtid<first_od_sipmid+1000 )
	hoopid = 900 + (pmtid-first_od_sipmid)/10; // 0-999 inclusive are radial
      else {
	if (pmtid<first_od_sipmid+1100)
	  hoopid = 1000;     // endcap
	else
	  hoopid = 1001;     // endcap
      }

    }

    if ( hoop_cut && ( hoopid<min_hoop || hoopid>max_hoop ) )
      continue;

    for (int ihit=0; ihit<nhits; ihit++) {
      RAT::DS::MCPhoton* hit = pmt->GetMCPhoton( ihit );
      double thit = hit->GetHitTime();
      if ( time_cut && (thit<min_time || thit>max_time) )
	continue;
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
    double ave_window = 0.;
    for (int i=ibin-windowbins; i<ibin; i++) {
      nhits_window += tbins[i];
    }
    ave_window = double(nhits_window)/double(windowbins);

    // expectations
    std::vector<double> pulse_expectation;
    pulse_expectation.resize( pulses.size() );

    // Check for new pulse
    if ( nactive==0 ) {
      // no active pulses. do search based on hits ver threshold of moving window
      if ( (!use_ave && nhits_window>(int)threshold ) || (use_ave && ave_window>threshold/windowbins) ) {
	// make a new pulse!
	KPPulse* apulse = new KPPulse;
	apulse->tstart = (ibin-windowbins)*nspertic;
	apulse->fStatus = KPPulse::kRising; // start of rising edge (until we find max)
	apulse->last_max = ave_window;
	pulses.push_back( apulse );
	npulses++;
// 	if ( veto )
// 	  std::cout << " ibin=" << ibin << " hits_window=" << nhits_window << " ave_window=" << ave_window << " thresh=" << threshold << " (" << threshold/float(windowbins) << ")" << std::endl;
      }
    }
    else {
      // we have active pulses. we look for a second peak with a trigger algorithm that accounts for scintillator decay time

      // find modified threshold
      double modthresh = threshold;
      if ( use_ave )
	modthresh = threshold/double(windowbins);
      bool allfalling = true;
      int ipulse = 0;
      for ( KPPulseListIter it=pulses.begin(); it!=pulses.end(); it++ ) {
	if ( (*it)->fStatus==KPPulse::kRising ) {
	  // have a rising peak. will block creation of new pulse.
	  allfalling = false;
	  if ( !use_ave )
	    modthresh += 2.0*nhits_window;
	  else
	    modthresh += 2.0*ave_window;
	  pulse_expectation.at(ipulse) = modthresh;
	}
	else {
	  // for pulses considered falling, we modify the threshold to be 3 sigma (roughly) above
	  //double expectation = ((*it)->peakamp)*exp( -( ibin*nspertic - (*it)->tpeak )/decay_constant ); // later can expand to multiple components
	  double arg = 0.0;
	  for (int idcy=0; idcy<n_decay_constants; idcy++) {
	    arg += decay_weights[idcy]*( (ibin-windowbins)*nspertic - (*it)->tpeak )/decay_constants_ns[idcy];
	  }
	  // calc expectation for bin: note 'peakamp' is always ave of hits in window
	  double expectation = 0;
	  if ( !use_ave ) 
	    expectation = ((*it)->peakamp)*windowbins*exp( -arg ); // window sum
	  else
	    expectation = ((*it)->peakamp)*exp( -arg );
	  modthresh += ( expectation + 3.0*sqrt(expectation) );
	  pulse_expectation.at(ipulse) = expectation;
	}
	ipulse++;
      }//end of loop over pulses

//       if ( veto )
// 	std::cout << " ibin=" << ibin << " hits_window=" << nhits_window << " ave_window=" << ave_window << " modthresh=" << modthresh << std::endl;

      // apply threshold
      if ( allfalling && ( (!use_ave && nhits_window>modthresh) || (use_ave && ave_window>modthresh) ) ) {
	// new pulse! (on top of other pulse)
	KPPulse* apulse = new KPPulse;
        apulse->tstart = (ibin-windowbins)*nspertic;
        apulse->fStatus = KPPulse::kRising; // start of rising edge (until we find max) 
	apulse->last_max = ave_window;
	pulses.push_back( apulse );
        npulses++;
	//std::cout << "  new pulse!" << std::endl;
      }

      // now we find max of rising pulses and end of falling pulses
      ipulse = 0;
      for ( KPPulseListIter it=pulses.begin(); it!=pulses.end(); it++ ) {
	KPPulse* apulse = *it;
	if ( apulse->fStatus==KPPulse::kDefined )
	  continue;

	if ( apulse->fStatus==KPPulse::kRising ) {
	  if ( ave_window < apulse->last_max )
	    apulse->nfallingbins += 1;
	  else {
	    apulse->nfallingbins = 0;
	    if ( ave_window>apulse->last_max )
	      apulse->last_max = ave_window;
	  }

	  if ( apulse->nfallingbins>NFALLING ) {
	    // found our max
	    apulse->fStatus = KPPulse::kFalling;
	    apulse->tpeak = (ibin-windowbins-apulse->nfallingbins)*nspertic;
	    apulse->peakamp = apulse->last_max; // note, always uses averages
	  }
	}//end of if rising
	else if ( apulse->fStatus==KPPulse::kFalling ) {
	  double decay_constant = 0.0;
	  for (int idcy=0; idcy<n_decay_constants; idcy++)
	    decay_constant += decay_weights[idcy]*decay_constants_ns[idcy];
	  //if (  pulse_expectation.at(ipulse) <=threshold*0.5 || ((ibin-windowbins)*nspertic > apulse->tpeak + 10*decay_constant) ) {
	  if ( ((ibin-windowbins)*nspertic > apulse->tpeak + 10*decay_constant) ) {
	    apulse->tend = (ibin-windowbins)*nspertic;
	    apulse->fStatus=KPPulse::kDefined;
	  }
	}
	ipulse++;
      }//end of loop over pulses
    }//end of active pulse condition
    
  }//end of scan over timing histogram

  // save unclosed pulses
  for ( KPPulseListIter it=pulses.begin(); it!=pulses.end(); it++ ) {
    KPPulse* apulse = *it;
    if ( apulse->fStatus!=KPPulse::kDefined ) {
      if ( apulse->fStatus==KPPulse::kRising ){
	apulse->peakamp = apulse->last_max;
	apulse->tpeak = ( nbins-windowbins-1 )*nspertic;
	apulse->tend  = ( nbins-windowbins )*nspertic;
      }
      apulse->fStatus = KPPulse::kDefined;
    }
  }

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
  if ( hoop_cut ) {
    // inner pipe
    if ( !veto )
      nsipms = (double)(max_hoop-min_hoop)*100.0; // 1 MHz dark rate * Window  * NSiPMs
    // outer pipe
    if ( veto ) {
      nsipms = 0;
      for (int ihoop=min_hoop; ihoop<=max_hoop; ihoop++) {
	if ( min_hoop<1000 )
	  nsipms += 10;  // sides
	else
	  nsipms += 100; // end caps
      }
    }
  }
  
  PMTinfo* pmtinfo = PMTinfo::GetPMTinfo( pmtinfofile );

  const int npulses = pulselist.size();
  double heights[ npulses ];
  double height_tot;

  // clear pe, z
  for ( int ipulse=0; ipulse<npulses; ipulse++ ) {
    KPPulse* apulse = pulselist.at(ipulse);
    apulse->pe = 0.0;
    apulse->pe_dark = 0.0;
    apulse->pe_adjusted = 0.0;
    apulse->z  = 0.0;
    apulse->hits_assigned = 0;
  }
  
  // Fill time bins
  int npmts = mc->GetMCPMTCount();
  int ntothits = 0;
  int ncuthits = 0;

  for ( int ipmt=0; ipmt<npmts; ipmt++ ) {
    RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
    int nhits = pmt->GetMCPhotonCount();
    int pmtid = pmt->GetID();

    // VETO OR INNER
    if ( veto && pmtid<first_od_sipmid )  {
      continue;
    }
    else if ( !veto && pmtid>=first_od_sipmid ) {
      continue;
    }

    int hoopid = pmtid/100;
    if ( !veto ) {
      // inner pipe hoop structure
    }
    else {
      // version 1
      // veto has string structure, which in retrospect, was dumb, we're interested in z
      // this uses info only in my head. correlates with build_geometry.py
//       if ( pmtid<first_od_sipmid+1000 )
// 	hoopid = 900+pmtid%100; // 0-999 inclusive are radial
//       else {
// 	if (pmtid<first_od_sipmid+1100)
// 	  hoopid = 1000;     // endcap
// 	else
// 	  hoopid = 1001;     // endcap
//       }
      // version 2
      if ( pmtid<first_od_sipmid+1000 )
	hoopid = 900 + (pmtid-first_od_sipmid)/10; // 0-999 inclusive are radial
      else {
	if (pmtid<first_od_sipmid+1100)
	  hoopid = 1000;     // endcap
	else
	  hoopid = 1001;     // endcap
      }
    }
    //std::cout << " assign pmtid=" << pmtid << " hoopid=" << hoopid << std::endl;
    if ( hoop_cut && ( hoopid<min_hoop || hoopid>max_hoop ) ) {
      ncuthits+=nhits;
      continue;
    }

    float pmtpos[3];
    pmtinfo->getposition( pmtid, pmtpos );
    //std::cout << " assign pmtid=" << pmtid << " hoopid=" << hoopid << " z=" << pmtpos[2] << std::endl;

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
//   	    if ( ipulse==0 )
// 	      std::cout << " pulse 1: pmtid=" << pmtid 
// 			<< "  height[ipulse]=" 
// 			<< heights[ipulse] 
// 			<< ", z=" << pmtpos[2] 
// 			<< " frac of pulse: " << heights[ipulse]/height_tot 
// 			<< " total: " << apulse->pe
// 			<< std::endl;
	  }
	} // loop over pulses
      }//end of if height>0
      else {
	ncuthits++;
      }
      ntothits++;
    }//end of loop over hits
  }//end of pmt loop
  //std::cout << "  total hits=" << ntothits << " ncut hits=" << ncuthits << std::endl;

  // take care of averages and subtract dark noise
  for ( int ipulse=0; ipulse<npulses; ipulse++ ) {
    KPPulse* apulse = pulselist.at(ipulse);
    if ( apulse->hits_assigned>0 ) 
      apulse->z /= double(apulse->hits_assigned);
    apulse->pe_adjusted = apulse->pe - apulse->pe_dark;
  }

}
