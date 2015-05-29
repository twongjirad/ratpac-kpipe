#include "pmtinfo.hh"
#include "gen_dark_noise.hh"
#include <iostream>
#include "TRandom3.h"
#include "RAT/DS/MCPhoton.hh"
#include "RAT/DS/MCPMT.hh"

void gen_dark_noise( RAT::DS::MC* mc, std::string pmtinfofile, double sipm_darkrate, int nodpmts, double window_ns ) {

  TRandom3* fRand = new TRandom3( time(NULL) );
  //const int NPMTS = 90000 + 2000; // ID + OD tubes
  const int NPMTS = 90000 + nodpmts;
  std::vector<int> hashit( NPMTS, 0 );

  double lambda = window_ns*(sipm_darkrate*1.0e-9);
  PMTinfo* pmtinfo = PMTinfo::GetPMTinfo( pmtinfofile );
  std::cout << "load pmt info: " << pmtinfofile << std::endl;

  // we have to go through list of SiPMS with hits first
  // then go through ones without
  int npe = 0;
  int ndark = 0;
  for (int ipmt=0; ipmt<mc->GetMCPMTCount(); ipmt++) {
    RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
    if ( pmt==NULL )
      continue;
    int pmtid = pmt->GetID();
    hashit.at(pmtid) = 1;
    npe += pmt->GetMCPhotonCount();

    // generate photons
    double ndarkhits = fRand->Poisson( lambda );
    npe += ndarkhits;
    ndark += ndarkhits;
    //std::cout << "pmt=" << pmtid << ", dark hits=" << ndarkhits << std::endl;
    for (int idark=0; idark<ndarkhits; idark++) {
      RAT::DS::MCPhoton* aphoton = pmt->AddNewMCPhoton();
      double hittime = fRand->Uniform( 10000.0 );
      aphoton->SetHitTime( hittime );
      aphoton->SetFrontEndTime( hittime );
      float pmtpos[3];
      pmtinfo->getposition( pmtid, pmtpos );
      aphoton->SetPosition( TVector3( pmtpos[0], pmtpos[1], pmtpos[2] ) );
      aphoton->SetCharge(1.0);
      aphoton->SetDarkHit(true);
      aphoton->SetOriginFlag( 4 );
      aphoton->SetTrackID( -1 );
    }
  }

  // Now add dark hits to pmts without hits
  int totpmts = 0;
  for (int ipmt=0; ipmt<NPMTS; ipmt++) {
    totpmts++;
    if ( hashit.at(ipmt)==1 ) {
      continue;
    }
    
    RAT::DS::MCPMT* pmt = mc->AddNewMCPMT();
    pmt->SetID( ipmt );
    
    // generate photons
    double ndarkhits = fRand->Poisson( lambda );
    npe += ndarkhits;
    ndark += ndarkhits;
//     if ( ipmt>=90000 ) 
//       std::cout << "pmt=" << ipmt << ", dark hits=" << ndarkhits <<  std::endl;
    for (int idark=0; idark<ndarkhits; idark++) {
      RAT::DS::MCPhoton* aphoton = pmt->AddNewMCPhoton();
      double hittime = fRand->Uniform( 10000.0 );
      aphoton->SetHitTime( hittime );
      aphoton->SetFrontEndTime( hittime );
      float pmtpos[3];
      pmtinfo->getposition( ipmt, pmtpos );
      aphoton->SetPosition( TVector3( pmtpos[0], pmtpos[1], pmtpos[2] ) );
      aphoton->SetCharge(1.0);
      aphoton->SetDarkHit(true);
      aphoton->SetOriginFlag( 4 );
      aphoton->SetTrackID( -1 );
    }
  }

  std::cout << "dark nose applied to " << totpmts << " pmts." << std::endl;

  mc->SetNumPE( npe );
  mc->SetNumDark( ndark );

  delete fRand;
}
