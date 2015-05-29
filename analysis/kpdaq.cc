#include "kpdaq.h"
#include "pmtinfo.hh"
#include "RAT/DS/MCPMT.hh"
#include "RAT/DS/MCPhoton.hh"
#include "TRandom3.h"
#include <ctime>
#include <assert.h>
#include <iostream>

KPDAQ::KPDAQ( int nidsipms_perhoop, int nodsipms_perhoop, int nidhoops, int nodhoops, int nidchperhoop, int version, std::string pmtinfofilename )
  : fNIDSiPMs_perhoop( nidsipms_perhoop ), fNODSiPMS_perhoop( nodsipms_perhoop), 
    fNIDhoops( nidhoops ), fNODhoops( nodhoops ), fNIDchperhoop( nidchperhoop ), fVersion( version ), pmtinfo_fname( pmtinfofilename )
{

  fNIDSiPMs = fNIDSiPMs_perhoop*fNIDhoops;
  fNODSiPMs = fNODSiPMS_perhoop*(fNODhoops-2) + 2*100;

  fNIDChannels = fNIDchperhoop*fNIDhoops;
  fNODChannels = fNODhoops;

  fNbins = 10000; // 10 microseconds
  
  fpmtinfo = PMTinfo::GetPMTinfo( pmtinfo_fname );

  // reserve the daq space
  fID_wfm.reserve( fNIDChannels );
  for (int i=0; i<fNIDChannels; i++ ) {
    std::vector<double> wfm( fNbins, 0.0 );
    //fID_wfm.emplace_back( wfm );
    fID_wfm.push_back( wfm );
  }

  fOD_wfm.reserve( fNODChannels );
  for (int i=0; i<fNODChannels; i++ ) {
    std::vector<double> wfm( fNbins, 0.0 );
    fOD_wfm.push_back( wfm );
  }

  fThinPE = false;

  frand = new TRandom3( time(NULL) );

  std::cout << "============================================================" << std::endl;
  std::cout << "Configured KPIPE DAQ" << std::endl;
  std::cout << "  NIDSiPMs: " << fNIDSiPMs << std::endl;
  std::cout << "  NODSiPMs: " << fNODSiPMs << std::endl;
  std::cout << "  NIDChannels: " << fNIDChannels << std::endl;
  std::cout << "  NODChannels: " << fNODChannels << std::endl;
  std::cout << "============================================================" << std::endl;
}


KPDAQ::~KPDAQ() 
{
  delete frand;
}

int KPDAQ::getHoopIndexFromPMTID( int pmtid ) {
  int hoopid = 0;

  if ( pmtid<fNIDSiPMs ) {
    hoopid = pmtid/fNIDSiPMs_perhoop;
  }
  else {
    
    // version 1
    if ( fVersion==1 ) {
      // veto has string structure, which in retrospect, was dumb, we're interested in z
      // this uses info only in my head. correlates with build_geometry.py
      if ( pmtid<fNIDSiPMs+1000 )
	hoopid = 900+pmtid%100; // 0-999 inclusive are radial
      else {
	if (pmtid<fNIDSiPMs+1100)
	    hoopid = 1000;     // endcap
	else
	  hoopid = 1001;     // endcap
      }
    }
    else if ( fVersion==2 ) {
      //version 2, this is fixed
      if ( pmtid<fNIDSiPMs+1000 )
	hoopid = 900 + (pmtid-fNIDSiPMs)/10; // 0-999 inclusive are radial
      else {
	if (pmtid<fNIDSiPMs+1100)
	  hoopid = 1000;     // endcap
	else
	  hoopid = 1001;     // endcap
      }
    }
    else if ( fVersion==3 ) {
      if ( pmtid<fNIDSiPMs+5000 )
	  hoopid = 900 + (pmtid-fNIDSiPMs)/50; // 0-999 inclusive are radial
      else {
	if (pmtid<fNIDSiPMs+5100)
	  hoopid = 1000;     // endcap
	else
	  hoopid = 1001;     // endcap
      }	
    }
    else if ( fVersion==4 ) {
      if ( pmtid<fNIDSiPMs+(fNODSiPMs-(2)*fNODSiPMS_perhoop))
	hoopid = pmtid/100; // 0-999 inclusive are radial
      else {
	if (pmtid<fNIDSiPMs+(fNODSiPMs-fNODSiPMS_perhoop))
	  hoopid = 1000;     // endcap
	else
	  hoopid = 1001;     // endcap
      }	
    }
    else 
      assert(false);
  }//end of od
  return hoopid;
}

int KPDAQ::getChannelIndex( int ihoop, int pmtid ) {
  int chnum = 0;
  if ( pmtid<fNIDSiPMs ) {
    // ID
    int hoop_ch = (pmtid%fNIDSiPMs_perhoop)/( (fNIDSiPMs_perhoop/fNIDchperhoop) );
    chnum = ihoop*fNIDchperhoop + hoop_ch;
  }
  else {
    // OD
    chnum = ihoop; 
  }
  return chnum;
}


void KPDAQ::processEvent(  RAT::DS::MC& mc ) {

  // clear daq
  for ( int i=0; i<fNIDChannels; i++ )
    fID_wfm.at(i).assign( fNbins, 0.0 );
  for (int j=0; j<fNODChannels; j++ )
    fOD_wfm.at(j).assign( fNbins, 0.0 );

  // fill wfms
  int removedpe = 0;
  int npmts = mc.GetMCPMTCount();
  //npmts = (fNIDSiPMs+fNODSiPMs);
  std::cout << "npmts: " << npmts << std::endl;
  for ( int ipmt=0; ipmt<npmts; ipmt++ ) {

    RAT::DS::MCPMT* pmt = mc.GetMCPMT( ipmt );

    if ( pmt==NULL )
      continue;

    int nhits = 0;
    int pmtid = 0;
    int hoopid = 0;
    int chnum = 0;

    nhits = pmt->GetMCPhotonCount();
    pmtid = pmt->GetID();
    hoopid = getHoopIndexFromPMTID( pmtid );
    chnum = getChannelIndex( hoopid, pmtid );
  
    for (int ihit=0; ihit<nhits; ihit++) {
      RAT::DS::MCPhoton* hit = pmt->GetMCPhoton( ihit );
      if ( fThinPE && !hit->IsDarkHit() ) {
	// reduce scintillator light yield in this way
	if ( frand->Uniform()>fThinPEfactor ) {
	  removedpe++;
	  continue;
	}
      }
      
      double thit = hit->GetHitTime();
      int ibin = (int)thit;
      if ( ibin>=0 && ibin<fNbins ) {

	try {
	  if ( chnum < fNIDChannels )
	    fID_wfm.at( chnum ).at( ibin ) += 1.0;
	  else
	    fOD_wfm.at( chnum-fNIDChannels ).at( ibin ) += 1.0;
	}
	catch (...) {
	  std::cout << "oops. " << chnum << " " << hoopid << " " << pmtid << std::endl;
	  assert(false);
	}
      }
    }//end of hit loop
  }//end of pmt loop


//   for (int i=0; i<fOD_wfm.size(); i++) {
//     double vetope = 0.0;
//     for (int j=0; j<10000; j++) 
//       vetope += fOD_wfm.at(i).at(j);
//     std::cout << "veto daq " << i+fNIDChannels << ": " << vetope << std::endl;
//   }
}

double KPDAQ::getWindowSum( int ch, double start, double end ) {
  double sum = 0;
  int istart = (int)start;
  int iend   = (int)end;

  if ( istart<0 ) istart = 0;
  if ( iend>=fNbins ) iend = fNbins-1;

  if ( ch<fNIDChannels )
    for ( int i=istart; i<=iend; i++ ) 
      sum += fID_wfm.at(ch).at(i);
  else
    for ( int i=istart; i<=iend; i++ ) 
      sum += fOD_wfm.at(ch-fNIDChannels).at(i);
    
  return sum;
}

void KPDAQ::copyWaveforms( std::vector<double>& copy, int chstart, int chend ) const {
  if ( copy.size()!=fNbins )
    copy.resize(fNbins);
  copy.assign(fNbins,0.0);

  int istart = chstart;
  int iend = chend;
  if ( istart<0 ) istart = 0;
  if ( iend>=(fNIDChannels+fNODChannels) ) iend = fNIDChannels+fNODChannels-1;

  for (int ich=istart; ich<=iend; ich++) {
    for ( int ibin=0; ibin<fNbins; ibin++ ) {
      if ( ich<fNIDChannels )
	copy.at(ibin) +=  fID_wfm.at(ich).at(ibin);
      else
	copy.at(ibin) += fOD_wfm.at(ich-fNIDChannels).at(ibin);
    }
  }
}

void KPDAQ::addWaveform( std::vector<double>& wfm, int ich, double tstart, double tend, double offset ) const {
  int istart = (int)tstart;
  int iend   = (int)tend;
  if ( iend>=fNbins ) iend = fNbins-1;
  for ( int ibin=istart; ibin<=iend; ibin++ ) {
    if ( ich<fNIDChannels )
      wfm.at(ibin) +=  fID_wfm.at(ich).at(ibin)-offset;
    else
      wfm.at(ibin) += fOD_wfm.at(ich-fNIDChannels).at(ibin)-offset;
  }  
}

void KPDAQ::getChannelPos( int ich, float* pos ) {
  int pmtid;
  if ( ich<fNIDChannels ) {
    pmtid = ich*fNIDSiPMs_perhoop + int(0.5*fNIDSiPMs_perhoop);
    fpmtinfo->getposition( pmtid, pos );
  }

}
