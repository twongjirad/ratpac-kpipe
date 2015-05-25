#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "RAT/DSReader.hh"
#include "RAT/DS/MC.hh"

#include "kptrigger.h"
#include "gen_dark_noise.hh"

int main( int nargs, char** argv ) {

  if ( nargs!=5 ) {
    std::cout << "usage: scrape_data <input RAT root file> <input cryfile> <output rootfile> <pmt info file>" << std::endl;
    return 0;
  }

  std::string inputfile = argv[1];
  std::string cryfile = argv[2];
  std::string outfile = argv[3];
  std::string pmtinfofile = "../data/kpipe/PMTINFO.root";
  pmtinfofile = argv[4];

  RAT::DSReader* ds = new RAT::DSReader( inputfile.c_str() ); 
  int first_od_sipmid = 90000;
  int n_decay_constants = 2;
  double decay_weights[2] = { 0.6, 0.4 };
  double decay_constants_ns[2] = { 45.0, 67.6 };



  // -----------------------------------------------------------------------------------------
  // CRY OUTPUT FILE. CONTAINS TRUTH ABOUT EVENT
  TChain* crytree = new TChain("crytree");
  crytree->Add( cryfile.c_str() );

  int nparticles;
  std::vector< int >* status = 0;
  std::vector< int >* pdg = 0;
  std::vector< double >* momx_gev = 0;
  std::vector< double >* momy_gev = 0;
  std::vector< double >* momz_gev = 0;
  std::vector< double >* mass_gev = 0;
  std::vector< double >* posx_mm = 0;
  std::vector< double >* posy_mm = 0;
  std::vector< double >* posz_mm = 0;
  std::vector< double >* hitx_mm = 0;
  std::vector< double >* hity_mm = 0;
  std::vector< double >* hitz_mm = 0;
  std::vector< double >* telapsed_sec = 0;
  std::vector< double >* delta_time_ns = 0;

  TBranch* b_status = 0;
  TBranch* b_pdg = 0;
  TBranch* b_momx_gev = 0;
  TBranch* b_momy_gev = 0;
  TBranch* b_momz_gev = 0;
  TBranch* b_mass_gev = 0;
  TBranch* b_posx_mm = 0;
  TBranch* b_posy_mm = 0;
  TBranch* b_posz_mm = 0;
  TBranch* b_hitx_mm = 0;
  TBranch* b_hity_mm = 0;
  TBranch* b_hitz_mm = 0;
  TBranch* b_telapsed_sec = 0;
  TBranch* b_delta_time_ns = 0;

  crytree->SetBranchAddress( "nparticles", &nparticles );
  crytree->SetBranchAddress( "pdg", &pdg, &b_pdg );
  crytree->SetBranchAddress( "momx_gev", &momx_gev, &b_momx_gev );
  crytree->SetBranchAddress( "momy_gev", &momy_gev, &b_momy_gev );
  crytree->SetBranchAddress( "momz_gev", &momz_gev, &b_momz_gev );
  crytree->SetBranchAddress( "mass_gev", &mass_gev, &b_mass_gev );
  crytree->SetBranchAddress( "posx_mm", &posx_mm, &b_posx_mm );
  crytree->SetBranchAddress( "posy_mm", &posy_mm, &b_posy_mm );
  crytree->SetBranchAddress( "posz_mm", &posz_mm, &b_posz_mm );
  crytree->SetBranchAddress( "hitx_mm", &hitx_mm, &b_hitx_mm );
  crytree->SetBranchAddress( "hity_mm", &hity_mm, &b_hity_mm );
  crytree->SetBranchAddress( "hitz_mm", &hitz_mm, &b_hitz_mm );
  crytree->SetBranchAddress( "telapsed_sec", &telapsed_sec, &b_telapsed_sec );
  crytree->SetBranchAddress( "delta_time_ns", &delta_time_ns, &b_delta_time_ns );

  // -----------------------------------------------------------------------------------------


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
  int npulses_veto = 0;
  std::vector<double> ttrig;
  std::vector<double> tpeak;
  std::vector<double> peakamp;
  std::vector<double> tend;
  std::vector<double> pulsepe;
  std::vector<double> pulsez;
  std::vector<double> twfm;
  std::vector<double> twfm_veto;
  // cosmic info
  int ncr_photons;
  int ncr_electrons;
  int ncr_muons;
  int ncr_neutrons;
  int ncr_other;
  std::vector<double> ke_crmuons;
  std::vector<double> ke_crphotons;
  std::vector<double> ke_crelectrons;
  std::vector<double> ke_crneutrons;
  std::vector<double> ke_crother;
  // dark rate
  
  // output tree
  TTree* tree = new TTree( "mcdata", "MC Data (Cosmic ray version)" );
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
  tree->Branch( "ttrig",  &ttrig );
  tree->Branch( "tpeak",  &tpeak );
  tree->Branch( "tend",  &tend );
  tree->Branch( "peakamp",  &peakamp );
  tree->Branch( "pulsepe",  &pulsepe );
  tree->Branch( "pulsez",  &pulsez );
  tree->Branch( "twfm",  &twfm );
  // cosmic info
  tree->Branch( "ncr_photons",  &ncr_photons, "ncr_photons/I" );
  tree->Branch( "ncr_electrons",  &ncr_electrons, "ncr_electrons/I" );
  tree->Branch( "ncr_muons",  &ncr_muons, "ncr_muons/I" );
  tree->Branch( "ncr_neutrons",  &ncr_neutrons, "ncr_neutrons/I" );
  tree->Branch( "ncr_other",  &ncr_other, "ncr_other/I" );
  tree->Branch( "ke_crphotons",  &ke_crphotons );
  tree->Branch( "ke_crelectrons",  &ke_crelectrons );
  tree->Branch( "ke_crmuons",  &ke_crmuons );
  tree->Branch( "ke_crneutrons",  &ke_crneutrons );
  tree->Branch( "ke_crother",  &ke_crother );


  int ievent = 0;
  int nevents = ds->GetTotal();

  KPPulseList pulselist;

  std::cout << "Number of events: " << nevents << std::endl;
  
  while (ievent<nevents) {
    RAT::DS::Root* root = ds->NextEvent();
    crytree->GetEntry(ievent);

    // --------------------------------
    // Clear Variables
    npe = idpe = odpe = 0;
    ncr_photons = ncr_electrons = ncr_muons = ncr_neutrons;
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
    ttrig.clear();
    tpeak.clear();
    tend.clear();
    peakamp.clear();
    pulsepe.clear();
    pulsez.clear();
    ke_crmuons.clear();
    ke_crphotons.clear();
    ke_crelectrons.clear();
    ke_crneutrons.clear();
    // --------------------------------

    if ( ievent%1000==0 )
      std::cout << "Event " << ievent << std::endl;

    // --------------------------------
    // PROCESS CRY FILE

    for (int icr=0; icr<nparticles; icr++) {
      if ( hitx_mm->at(icr)==0 && hity_mm->at(icr)==0 && hitz_mm->at(icr)==0 )
	continue;

      double crke = sqrt( momx_gev->at(icr)*momx_gev->at(icr) + momy_gev->at(icr)*momy_gev->at(icr) + momz_gev->at(icr)*momz_gev->at(icr) );

      if ( pdg->at(icr)==13 || pdg->at(icr)==-13) {
	ncr_muons++;
	ke_crmuons.push_back( crke );
      }
      else if ( pdg->at(icr)==11 || pdg->at(icr)==-11) {
	ncr_electrons++;
	ke_crelectrons.push_back( crke );
      }
      else if ( pdg->at(icr)==22 ) {
	ncr_photons++;
	ke_crphotons.push_back( crke );
      }
      else if ( pdg->at(icr)==2112 ) {
	ncr_neutrons++;
	ke_crneutrons.push_back( crke );
      }
      else {
	ncr_other++;
	ke_crother.push_back( crke );
      }

    }

    // --------------------------------
    // GET RAT MC OBJECT
    RAT::DS::MC* mc = root->GetMC();

    // --------------------------------
    // GENERATE DARK NOISE
    gen_dark_noise( mc, pmtinfofile, 1.0e6, 10000 );
    
    // --------------------------------
    // PROCESS RAT FILE
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

    // TRIGGER
    npulses = find_trigger( mc, 5.0, 5.0, 10.0, 
			    false, 0, 0,
			    false, 0, 0,
			    n_decay_constants, decay_weights, decay_constants_ns,
			    pulselist, 90000, false, twfm );

    assign_pulse_charge( mc, pmtinfofile, pulselist,
			 1.0e6,
			 false, 0, 0,
			 45.0, 90000, false );
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
