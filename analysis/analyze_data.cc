#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "RAT/DSReader.hh"
#include "RAT/DS/MC.hh"

#include "kptrigger.h"
#include "gen_dark_noise.hh"
#include "prefitz.hh"

int main( int nargs, char** argv ) {

  if ( nargs<4 ) {
    std::cout << "usage: scrape_data <input RAT root file> <output rootfile> <pmt info file>" << std::endl;
    return 0;
  }

  std::string inputfile = argv[1];
  std::string outfile = argv[2];
  std::string pmtinfofile = "../data/kpipe/PMTINFO.root";
  pmtinfofile = argv[3];

  RAT::DSReader* ds = new RAT::DSReader( inputfile.c_str() ); 
  int first_od_sipmid = 90000;
  int n_decay_constants = 2;
  double decay_weights[2] = { 0.6, 0.4 };
  double decay_constants_ns[2] = { 45.0, 67.6 };

  TFile* out = new TFile(outfile.c_str(), "RECREATE" );
  // variables we want
  int npe, idpe, odpe;
  int npmts, idpmts, odpmts;
  int nhoops;
  double posv[3]; // true vertex
  double mudirv[3]; // true direction
  double muendr;
  double muendv[3]; // projected muon endpoint
  double mumomv;  // truth momentum (0 if no muon! NC event)
  double totkeprotonv;  // truth momentum (0 if no muon! NC event)
  double rv, zv;  // from truth
  // trigger info
  int npulses = 0;
  double prefit_z_cm = 0;
  std::vector<double> ttrig;
  std::vector<double> tpeak;
  std::vector<double> peakamp;
  std::vector<double> tend;
  std::vector<double> pulsepe;
  std::vector<double> pulsez;
  std::vector<double> twfm;
//   double ttrig[10];
//   double tpeak[10];
//   double peakamp[10];
//   double tend[10];
//   double pulsepe[10];
  
  TTree* tree = new TTree( "mcdata", "MC Data" );
  // recon
  tree->Branch( "npe", &npe, "npe/I" );
  tree->Branch( "idpe", &idpe, "idpe/I" );
  tree->Branch( "odpe", &odpe, "odpe/I" );
  tree->Branch( "npmts", &npmts, "npmts/I" );
  tree->Branch( "idpmts", &idpmts, "idpmts/I" );
  tree->Branch( "odpmts", &odpmts, "odpmts/I" );
  // truth
  tree->Branch( "posv", posv, "posv[3]/D" );
  tree->Branch( "mudirv", mudirv, "mudirv[3]/D" );
  tree->Branch( "muendv", muendv, "muendv[3]/D" );
  tree->Branch( "muendr", &muendr, "muendr/D" );
  tree->Branch( "mumomv", &mumomv, "mumomv/D" );
  tree->Branch( "totkeprotonv", &totkeprotonv, "totkeprotonv/D" );
  tree->Branch( "rv", &rv, "rv/D" );  
  tree->Branch( "zv", &zv, "zv/D" );
  // trigger vars
  tree->Branch( "npulses", &npulses, "npulses/I" );
  tree->Branch( "prefit_z_cm", &prefit_z_cm, "prefit_z_cm/D" );
  tree->Branch( "ttrig",  &ttrig );
  tree->Branch( "tpeak",  &tpeak );
  tree->Branch( "tend",  &tend );
  tree->Branch( "peakamp",  &peakamp );
  tree->Branch( "pulsepe",  &pulsepe );
  tree->Branch( "pulsez",  &pulsez );
  tree->Branch( "twfm", &twfm );

  int ievent = 0;
  int nevents = ds->GetTotal();
  nevents = 30;

  KPPulseList pulselist;

  std::cout << "Number of events: " << nevents << std::endl;
  
  while (ievent<nevents) {
    RAT::DS::Root* root = ds->NextEvent();

    // --------------------------------
    // Clear Variables
    npe = idpe = odpe = 0;
    npmts = idpmts = odpmts = 0;
    for (int i=0; i<3; i++)
      posv[i] = 0.0;
    rv = zv = 0.0;
    mumomv = 0;
    totkeprotonv = 0.;
    mudirv[0] = mudirv[1] = mudirv[2] = 0.0;
    muendv[0] = muendv[1] = muendv[2] = 0.0;
    muendr = 0;
    npulses = 0;
    prefit_z_cm = 0;
    ttrig.clear();
    tpeak.clear();
    tend.clear();
    peakamp.clear();
    pulsepe.clear();
    pulsez.clear();
    twfm.clear();
    // --------------------------------

    if ( ievent%1000==0 )
      std::cout << "Event " << ievent << std::endl;

    RAT::DS::MC* mc = root->GetMC();

    std::cout << "event " << ievent << ", npe=" << mc->GetNumPE() << " in " << mc->GetMCPMTCount() << std::endl;
    gen_dark_noise( mc, pmtinfofile, 1.0e6, 10000 );
    std::cout << "  finished dark noise gen, npe=" << mc->GetNumPE() << " in " << mc->GetMCPMTCount() << std::endl;

    if ( mc==NULL )
      break;
    npe = mc->GetNumPE();
    npmts = mc->GetMCPMTCount();

    if ( mc->GetMCParticleCount()==0 ) {
      tree->Fill();
      ievent++;
      continue;
    }

    // true vertex
    //std::cout << "mc part: " << mc->GetMCParticleCount() << " " <<  mc->GetMCParticle(0) << std::endl;
    posv[0] = mc->GetMCParticle(0)->GetPosition().X()/10.0; //change to cm
    posv[1] = mc->GetMCParticle(0)->GetPosition().Y()/10.0; //change to cm
    posv[2] = mc->GetMCParticle(0)->GetPosition().Z()/10.0; //change to cm
    rv = sqrt(posv[0]*posv[0] + posv[1]*posv[1]);
    zv = posv[2];

    // true muon momentum
    for (int ipart=0; ipart<mc->GetMCParticleCount(); ipart++) {
      if ( mc->GetMCParticle(ipart)->GetPDGCode()==13 ) {
	TVector3 mom( mc->GetMCParticle(ipart)->GetMomentum() );
	mumomv = sqrt( mom.X()*mom.X() + mom.Y()*mom.Y() + mom.Z()*mom.Z() );
	mudirv[0] = mom.X()/mumomv;
	mudirv[1] = mom.Y()/mumomv;
	mudirv[2] = mom.Z()/mumomv;
	double muke = mc->GetMCParticle(ipart)->GetKE();
	
	for (int i=0; i<3; i++) {
	  muendv[i] = posv[i] + mudirv[i]*(muke/1.7);
	}
	muendr = sqrt( muendv[0]*muendv[0] + muendv[1]*muendv[1] );
      }
      else if ( mc->GetMCParticle(ipart)->GetPDGCode()==2212 ) {
	totkeprotonv += mc->GetMCParticle(ipart)->GetKE();
      }
    }

    // count stuff
    npmts = mc->GetMCPMTCount();
    idpmts = odpmts = 0;
    npe = mc->GetNumPE();
    idpe = odpe = 0;
    for ( int ipmt=0; ipmt<npmts; ipmt++ ) {
      RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
      int nhits = pmt->GetMCPhotonCount();
      int pmtid = pmt->GetID();

      if ( pmtid<first_od_sipmid ) {
	idpe += nhits;
	idpmts++;
      }
      else {
	odpe += nhits;
	odpmts++;
      }
    }

    std::cout << "------------------------------------------" << std::endl;
    std::cout << "EVENT " << ievent << std::endl;
    std::cout << "  ID PEs: " << idpe << " PMTs: " << idpmts << std::endl;
    std::cout << "  OD PEs: " << odpe << " PMTs: " << odpmts << std::endl;

    int maxhoop = 0;
    prefit_z_cm = calc_prefitz( mc, pmtinfofile, 1.0e6, 500.0, 90000, 2000, maxhoop );
    std::cout << "  prefit z: " << prefit_z_cm << " maxhoop=" << maxhoop << std::endl;
    int min_hoopid = maxhoop-50;
    int max_hoopid = maxhoop+50;
    if ( min_hoopid<0 )
      min_hoopid = 0;
    if ( max_hoopid>=900 )
      max_hoopid = 900;

    // TRIGGER
    npulses = find_trigger( mc, 500.0+50.0, 5.0, 10.0, 
			    true, min_hoopid, max_hoopid,
			    n_decay_constants, decay_weights, decay_constants_ns,
			    pulselist, 90000, false, twfm );

    assign_pulse_charge( mc, pmtinfofile, pulselist, 45.0, 90000, false );
    std::cout << "  posv: " << posv[0] << ", " << posv[1] << ", " << posv[2] << std::endl;
    std::cout << "  npulses=" << npulses << std::endl;
    for ( KPPulseListIter it=pulselist.begin(); it!=pulselist.end(); it++ )
      std::cout << "    - tstart=" << (*it)->tstart << " tpeak=" << (*it)->tpeak << " pe=" << (*it)->pe << " z=" << (*it)->z << std::endl;

    for ( KPPulseListIter it=pulselist.begin(); it!=pulselist.end(); it++ ) {
      ttrig.push_back( (*it)->tstart );
      tpeak.push_back( (*it)->tpeak );
      tend.push_back( (*it)->tend );
      peakamp.push_back( (*it)->peakamp );
      pulsepe.push_back( (*it)->pe ); // not yet calculated
      pulsez.push_back( (*it)->z ); // not yet calculated
      delete *it;
      *it = NULL;
    }

    pulselist.clear();
    ievent++;
   
    tree->Fill();
    //std::cin.get();
  }//end of while loop

  std::cout << "write." << std::endl;

  out->cd();
  tree->Write();

  std::cout << "finished." << std::endl;

  return 0;
}
