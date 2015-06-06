#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "RAT/DSReader.hh"
#include "RAT/DS/MC.hh"

#include "tree3.h"
#include "PMTinfo.hh"

int main( int nargs, char** argv ) {

  // ------------------------------------------------------------------
  // Purpose of this code is to reach back into the truth information
  // and record some stuff.

  // ------------------------------------------------------------------
  // ARGUMENTS
  if ( nargs!=6 ) {
    std::cout << "usage: scrape_data <input RAT root file> <input KDAR source file> <pmt info file> <output rootfile> <offset>" << std::endl;
    return 0;
  }

  std::string pmtinfofile = "../data/kpipe/PMTINFO.root";

  std::string input_rat = argv[1];
  std::string input_kdar = argv[2];
  pmtinfofile = argv[3];
  std::string output = argv[4];
  int offset = atoi(argv[5]);

  PMTinfo* pmtinfo = PMTinfo::GetPMTinfo( pmtinfofile );
    
  // --------------------------------
  // CONSTANTS/PARAMETERS
  double detR = 150; // cm
  double detZ = 6000; // cm
  double posv[3];

  // --------------------------------
  // LOAD RAT FILE

  RAT::DSReader* ds = new RAT::DSReader( input_rat.c_str() ); 

  // --------------------------------
  // LOAD KDAR FILE
  TFile* finput = new TFile( input_kdar.c_str(), "open");
  TTree* tevents = (TTree*)finput->Get( "treeout" );
  tree data( tevents );
  

  // -----------------------------------------------------------------------------------------
  // DEFINE OUTPUT

  TFile* out = new TFile(output.c_str(), "RECREATE" );

  // variables we want
  int enu;  // starting neutrino energy
  int cc;   // 1 if CC
  int ccqe; // 1 if CCQE
  int kdar; // kdar neutrino event
  double evisv;  // visible energy
  
  TTree* outtree = new TTree( "mcdata", "MC Data" );
  // recon
  outtree->Branch( "cc", &cc, "cc/I" );
  outtree->Branch( "ccqe", &ccqe, "ccqe/I" );
  outtree->Branch( "kdar", &kdar, "kdar/I" );
  outtree->Branch( "enu", &enu, "enu/D" );

  int ievent = 20;
  int nevents = ds->GetTotal();
  nevents = 30;

  std::cout << "Number of events: " << nevents << std::endl;

  // --------------------------------
  // Event Loop
  while (ievent<nevents) {
    RAT::DS::Root* root = ds->GetEvent(ievent);

    // --------------------------------
    // Clear Variables

    //summary vars
    enu = 0;
    cc = 0;
    ccqe = 0;
    kdar = 0;
    evisv = 0.0;

    // --------------------------------

    if ( ievent%1000==0 )
      std::cout << "Event " << ievent << std::endl;

    std::cout << "------------------------------------------" << std::endl;
    std::cout << "EVENT " << ievent << std::endl;

    // --------------------------------
    // PROCESS CRY FILE


    // --------------------------------
    // GET RAT MC OBJECT
    RAT::DS::MC* mc = root->GetMC();
    if ( mc==NULL )
      break;

    if ( mc->GetMCParticleCount()==0 ) {
      outtree->Fill();
      ievent++;
      continue;
    }

    // ==================================================================================================================

    // posv
    posv[0] = mc->GetMCParticle(0)->GetPosition().X()/10.0; //change to cm
    posv[1] = mc->GetMCParticle(0)->GetPosition().Y()/10.0; //change to cm
    posv[2] = mc->GetMCParticle(0)->GetPosition().Z()/10.0; //change to cm

    // enu
    enu = data.in_t[0];

    // KDAR: Josh's esoteric criteron
    if ( enu>235 && enu<237 )
      kdar = 1;

    // mode flags
    if (  data.flag_cc ) {
      cc = 1;
      if ( data.flag_qel )
	ccqe = 1;
      else
	ccqe = 0;
    }
    else
      cc = 0;

//     // visible energy: muon
//     if ( abs(data.post_pdg[0])==13 ) {
//       double norm = sqrt( data.post_x[0]*data.post_x[0] + data.post_y[0]*data.post_y[0] + data.post_z[0]*data.post_z[0] );
//       double v[3] = { post_x[0]/norm, post_y[0]/norm, post_z[0]/norm };

//       // we need distance along direction to edge of det
//       // endcaps
//       double sz = 0;
//       if ( v[2]>0 )
// 	sz = detZ - posv[2];
//       else
// 	sz = -detZ - posv[2];
//       sz /= v[2];

//       double endxy[2] = { sz*v[0] + posv[0], sz*v[1] + posv[1] };
//       double endr = sqrt( endxy[0]*endxy[0] + endxy[1]*endxy[1] );

//       // sides: circle line intersection
//       double x1[2] = { posv[0], posv[1] };
//       double x2[2] = { posv[0]+v[0], posv[2]+v[1] };
//       double signvy = v[1]/fabs(v[1]);
//       double dr = sqrt(v[0]*v[0] + v[1]*v[1] );
//       double D = x1[0]*x2[1] - x2[0]*x1[1];

//       double discr = detR*detR*dr*dr - D*D;
//       if ( discr >= 0 ) {
// 	// intersection occurs
	

//       }

//     }//visible energy calc for muon


    // ==================================================================================================================

    ievent++;
   
    outtree->Fill();
    //std::cin.get();
  }//end of while loop

  std::cout << "write." << std::endl;

  out->cd();
  outtree->Write();

  std::cout << "finished." << std::endl;

  return 0;
}
