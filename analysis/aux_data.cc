#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "RAT/DSReader.hh"
#include "RAT/DS/MC.hh"

#include "tree3.h"
#include "pmtinfo.hh"

#include "geotools.h"

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
  float posv[3];

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
  int cc;   // 1 if CC
  int ccqe; // 1 if CCQE
  int kdar; // kdar neutrino event
  double enu;  // starting neutrino energy
  double evisv;  // visible energy
  geotools::GeoDef geom;
  geom.radius = 1500;
  geom.height = 2.0*45000;
  
  TTree* outtree = new TTree( "mcdata", "MC Data" );
  // recon
  outtree->Branch( "cc", &cc, "cc/I" );
  outtree->Branch( "ccqe", &ccqe, "ccqe/I" );
  outtree->Branch( "kdar", &kdar, "kdar/I" );
  outtree->Branch( "enu", &enu, "enu/D" );
  outtree->Branch( "posv", posv, "posv[3]/D" );
  outtree->Branch( "evisv", &evisv, "evisv/D" );

  int ievent = 0;
  int nevents = ds->GetTotal();
  //nevents = 10;

  std::cout << "Number of events: " << nevents << std::endl;

  // --------------------------------
  // Event Loop
  while (ievent<nevents) {
    RAT::DS::Root* root = ds->GetEvent(ievent);
    data.GetEntry( offset + ievent );

    // --------------------------------
    // Clear Variables

    //summary vars
    enu = 0;
    cc = 0;
    ccqe = 0;
    kdar = 0;
    evisv = 0.0;
    posv[0] = posv[1] = posv[2] = 0.0;

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
    double rv = sqrt( posv[0]*posv[0] + posv[1]*posv[1] );
    double zv = posv[2];

    // enu
    enu = data.in_t[0];

    // KDAR: Josh's esoteric criteron
    if ( enu>235.0 && enu<237.0 )
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

    evisv = 0;
    if ( rv<geom.radius && fabs(zv)<0.5*geom.height ) {
      // only if in target volume
      for (int ipart=0; ipart<data.post_; ipart++) {
	// visible energy: muon
	if ( abs(data.post_pdg[ipart])==13 ) {
	  double norm = sqrt( data.post_x[ipart]*data.post_x[ipart] + data.post_y[ipart]*data.post_y[ipart] + data.post_z[ipart]*data.post_z[ipart] );
	  float dir[3] = { data.post_x[ipart]/norm, data.post_y[ipart]/norm, data.post_z[ipart]/norm }; // dir
	  double muon_ke_mev = data.post_t[ipart]-105.7; // Energy
	  
	  // we want to find intersection against the wall
	  geotools::Intersection ans;
	  geotools::RayCylinderIntersection( posv, dir, geom, ans );
	  if ( ans.intersectionType==geotools::kForBack ) {
	    double dist2wall = ans.distToNearFor;
	    double est_range = muon_ke_mev/(1.7*0.86); // stopping power, density, distance
	    if ( dist2wall<est_range )
	      evisv+=muon_ke_mev;
	    else
	      evisv+=muon_ke_mev*(est_range/dist2wall);
	  }
	}
	else if ( abs(data.post_pdg[ipart])==2212 ) {
	  // proton
	  // assume range is short
	  evisv += data.post_t[ipart];
	}
      }
    }//if true vertex
      



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
