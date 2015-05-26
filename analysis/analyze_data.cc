#include <iostream>
#include <assert.h>

#include "TFile.h"
#include "TTree.h"

#include "RAT/DSReader.hh"
#include "RAT/DS/MC.hh"

#include "kptrigger.h"
#include "gen_dark_noise.hh"
#include "prefitz.hh"
#include "pmtinfo.hh"


int main( int nargs, char** argv ) {

  // --------------------------------
  // ARGUMENTS
  if ( nargs!=4 and nargs!=6) {
    std::cout << "usage: scrape_data <input RAT root file> <output rootfile> <pmt info file> [optional: <crymode=1> <cry file>]" << std::endl;
    return 0;
  }

  std::string inputfile = argv[1];
  std::string outfile = argv[2];
  std::string pmtinfofile = "../data/kpipe/PMTINFO.root";
  pmtinfofile = argv[3];
  std::cout << "pmt info file: " << pmtinfofile << std::endl;
  PMTinfo* pmtinfo = PMTinfo::GetPMTinfo( pmtinfofile );
  bool cry_mode  = false;
  std::string cryfile = "tears.root";
  if (nargs==6) {
    cry_mode = true;
    cryfile = argv[5];
  }
    
  // --------------------------------
  // CONSTANTS/PARAMETERS

  RAT::DSReader* ds = new RAT::DSReader( inputfile.c_str() ); 
  int first_od_sipmid = 90000;
  int n_decay_constants = 2;
  double window_ns = 40.0;
  double window_ns_veto = 10.0;
  double decay_weights[2] = { 0.6, 0.4 };
  double decay_constants_ns[2] = { 45.0, 67.6 };
  int n_decay_constants_veto = 1;
  double decay_weights_veto[1] = { 1.0 };
  double decay_constants_ns_veto[1] = { 50.0 };
  int trig_version = 3;
  int nod_sipms_per_hoop = 50;
  int nod_sipms_per_hoop_endcap = 100;
  int nodpmts[4] = {0,1200,1200,5200}; 
  int nodhoops[4] = { 0,102,102,102 };
  const int NPMTS = 90000+nodpmts[trig_version];

  double sipm_darkrate_hz = 10.0e6;
  double threshold = 500.0;
//   double sipm_darkrate_hz = 0.0;
//   double threshold = 10.0;

  // --------------------------------
  // INPUT CRY VARS

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
  // DEFINE OUTPUT

  TFile* out = new TFile(outfile.c_str(), "RECREATE" );
  // variables we want
  int npe, idpe, odpe, predark_idpe, predark_odpe;
  int npmts, idpmts, odpmts;
  int nhoops;
  // truth
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
  double prefit_z_cm = 0;
  double pulse_totodpe = 0.0;
  // inner pipe varso
  std::vector<double> ttrig;
  std::vector<double> tpeak;
  std::vector<double> peakamp;
  std::vector<double> tend;
  std::vector<double> pulsepe;
  std::vector<double> pulsepedark;
  std::vector<double> pulseperaw;
  std::vector<double> pulsez;
  // outer pipe vars
  std::vector<double> pulsez_veto;
  std::vector<double> pulsepe_veto;
  std::vector<double> ttrig_veto;
  std::vector<double> tend_veto;
  std::vector<double> twfm;
  std::vector<double> twfm_veto;
  // cry cosmic info
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


  TTree* tree = new TTree( "mcdata", "MC Data" );
  // recon
  tree->Branch( "npe", &npe, "npe/I" );
  tree->Branch( "idpe", &idpe, "idpe/I" );
  tree->Branch( "odpe", &odpe, "odpe/I" );
  tree->Branch( "predark_idpe", &predark_idpe, "predark_idpe/I" );
  tree->Branch( "predark_odpe", &predark_odpe, "predark_odpe/I" );
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
  tree->Branch( "npulses_veto", &npulses_veto, "npulses_veto/I" );
  tree->Branch( "pulse_totodpe", &pulse_totodpe, "pulse_totodpe/D" );
  tree->Branch( "prefit_z_cm", &prefit_z_cm, "prefit_z_cm/D" );
  // inner pipe pusles
  tree->Branch( "ttrig",  &ttrig );
  tree->Branch( "tpeak",  &tpeak );
  tree->Branch( "tend",  &tend );
  tree->Branch( "peakamp",  &peakamp );
  tree->Branch( "pulsepe",  &pulsepe );
  tree->Branch( "pulsepedark",  &pulsepedark );
  tree->Branch( "pulseperaw",  &pulseperaw );
  tree->Branch( "pulsez",  &pulsez );
  // outer pipe pusles
  tree->Branch( "pulsepe_veto",  &pulsepe_veto );
  tree->Branch( "pulsez_veto",  &pulsez_veto );
  tree->Branch( "ttrig_veto",  &ttrig_veto );
  tree->Branch( "tend_veto",  &tend_veto );
//   tree->Branch( "twfm", &twfm );
  //tree->Branch( "twfm_veto", &twfm_veto );
  // cosmic truth
  if ( cry_mode ) {
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
  }


  int ievent = 0;
  int nevents = ds->GetTotal();
  nevents = 50;

  KPPulseList pulselist;
  KPPulseList pulselist_veto;

  std::cout << "Number of events: " << nevents << std::endl;
  
  while (ievent<nevents) {
    RAT::DS::Root* root = ds->GetEvent(ievent);
    if ( cry_mode )
      crytree->GetEntry( ievent );

    // --------------------------------
    // Clear Variables

    //summary vars
    npe = idpe = odpe = predark_idpe = predark_odpe = 0;
    npmts = idpmts = odpmts = 0;
    // truth
    for (int i=0; i<3; i++)
      posv[i] = 0.0;
    rv = zv = 0.0;
    mumomv = 0;
    totkeprotonv = 0.;
    mudirv[0] = mudirv[1] = mudirv[2] = 0.0;
    muendv[0] = muendv[1] = muendv[2] = 0.0;
    muendr = 0;
    npulses = 0;
    npulses_veto = 0;
    prefit_z_cm = 0;
    pulse_totodpe = 0;
    // inner pipe pulses
    ttrig.clear();
    tpeak.clear();
    tend.clear();
    peakamp.clear();
    pulsepe.clear();
    pulsepedark.clear();
    pulseperaw.clear();
    pulsez.clear();
    twfm.clear();
    twfm_veto.clear();
    // veto pulses
    pulsepe_veto.clear();
    pulsez_veto.clear();
    ttrig_veto.clear();
    tend_veto.clear();
    // cry/cosmics
    if ( cry_mode ) {
      ncr_photons = ncr_electrons = ncr_muons = ncr_neutrons = 0;
      ke_crmuons.clear();
      ke_crphotons.clear();
      ke_crelectrons.clear();
      ke_crneutrons.clear();
    }

    // --------------------------------

    if ( ievent%1000==0 )
      std::cout << "Event " << ievent << std::endl;

    std::cout << "------------------------------------------" << std::endl;
    std::cout << "EVENT " << ievent << std::endl;

    // --------------------------------
    // PROCESS CRY FILE

    if ( cry_mode ) {
      for (int icr=0; icr<nparticles; icr++) {
	if ( hitx_mm->at(icr)==0 && hity_mm->at(icr)==0 && hitz_mm->at(icr)==0 )
	  continue;

	double crke = sqrt( momx_gev->at(icr)*momx_gev->at(icr) 
			    + momy_gev->at(icr)*momy_gev->at(icr) 
			    + momz_gev->at(icr)*momz_gev->at(icr) );

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
      std::cout << "  ncr_muons=" << ncr_muons << ", ncr_photos=" << ncr_photons << " neutrons=" << ncr_neutrons << " ncr_other=" << ncr_neutrons << std::endl;
    }//end of CRY mode

    // --------------------------------
    // GET RAT MC OBJECT
    RAT::DS::MC* mc = root->GetMC();
    if ( mc==NULL )
      break;
    npe = mc->GetNumPE();
    npmts = mc->GetMCPMTCount();

    if ( mc->GetMCParticleCount()==0 ) {
      tree->Fill();
      ievent++;
      continue;
    }

    // --------------------------------
    // GENERATE DARK NOISE
    if ( sipm_darkrate_hz>0 )  {

      npmts = mc->GetMCPMTCount();
      idpmts = odpmts = 0;
      npe = mc->GetNumPE();
      predark_idpe = predark_odpe = 0;
      int num = npmts;
      for ( int ipmt=0; ipmt<num; ipmt++ ) {
	RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
	int nhits = pmt->GetMCPhotonCount();
	int pmtid = pmt->GetID();
	
	if ( pmtid<first_od_sipmid ) {
	  predark_idpe += nhits;
	  idpmts++;
	}
	else {
	  predark_odpe += nhits;
	  odpmts++;
	}
      }
      std::cout << "  pre-dark noise ID PEs: " << predark_idpe << " PMTs: " << idpmts << std::endl;
      std::cout << "  pre-dark noise OD PEs: " << predark_odpe << " PMTs: " << odpmts << std::endl;      
      gen_dark_noise( mc, pmtinfofile, sipm_darkrate_hz, nodpmts[trig_version], 10000 );
    }

    // --------------------------------
    // PROCESS RAT FILE

    // true vertex
    std::cout << "mc part: " << mc->GetMCParticleCount() << " " <<  mc->GetMCParticle(0) << std::endl;
    posv[0] = mc->GetMCParticle(0)->GetPosition().X()/10.0; //change to cm
    posv[1] = mc->GetMCParticle(0)->GetPosition().Y()/10.0; //change to cm
    posv[2] = mc->GetMCParticle(0)->GetPosition().Z()/10.0; //change to cm
    rv = sqrt(posv[0]*posv[0] + posv[1]*posv[1]);
    zv = posv[2];
    std::cout << "  posv: " << posv[0] << ", " << posv[1] << ", " << posv[2] << " rv=" << rv << std::endl;

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
      std::cout << "  particle=" << ipart << " pdg=" << mc->GetMCParticle(ipart)->GetPDGCode() << " ke=" << mc->GetMCParticle(ipart)->GetKE() << std::endl;
    }

    // count stuff (post dark noise)
    npmts = mc->GetMCPMTCount();
    idpmts = odpmts = 0;
    npe = mc->GetNumPE();
    idpe = odpe = 0;
    int num = npmts;
    if ( sipm_darkrate_hz>0 )
      num = NPMTS;
    for ( int ipmt=0; ipmt<num; ipmt++ ) {
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

    std::cout << "  ID PEs: " << idpe << " PMTs: " << idpmts << std::endl;
    std::cout << "  OD PEs: " << odpe << " PMTs: " << odpmts << std::endl;

    // --------------------------------
    // Find Z of the first pulse
    int maxhoop = 0;
    prefit_z_cm = calc_prefitz( mc, pmtinfofile, sipm_darkrate_hz, 500.0, 90000, 2000, maxhoop );
    std::cout << "  prefit z: " << prefit_z_cm << " maxhoop=" << maxhoop << std::endl;
    int min_hoopid = maxhoop-75;
    int max_hoopid = maxhoop+75;
    if ( min_hoopid<0 )
      min_hoopid = 0;
    if ( max_hoopid>=900 )
      max_hoopid = 900;

    double expected_darkrate = (double( max_hoopid - min_hoopid )*100)*(sipm_darkrate_hz*1.0e-9)*window_ns;
    double sig_darkrate = sqrt( expected_darkrate );
    threshold = expected_darkrate + 5.0*sig_darkrate;

    // --------------------------------
    // TRIGGER
    std::cout << "  -- pulse finder -- " << std::endl;

    // INNER PIPE
    npulses = find_trigger( mc, 
			    threshold, window_ns, sipm_darkrate_hz,
			    true, min_hoopid, max_hoopid,
			    false, 0, 0,
			    n_decay_constants, decay_weights, decay_constants_ns,
			    pulselist, 90000, false, twfm );

    assign_pulse_charge( mc, pmtinfofile, pulselist, 
			 sipm_darkrate_hz,
			 true, min_hoopid, max_hoopid,
			 60.0, 90000, nod_sipms_per_hoop, nod_sipms_per_hoop_endcap,
			 false );
    std::cout << "  ID npulses=" << npulses << "  with threshold=" <<  threshold << std::endl;
    for ( KPPulseListIter it=pulselist.begin(); it!=pulselist.end(); it++ )
      std::cout << "    - tstart=" << (*it)->tstart 
		<< " tpeak=" << (*it)->tpeak 
		<< " tend=" << (*it)->tend
		<< " pe=" << (*it)->pe << " (dark=" << (*it)->pe_dark << ", adjusted=" << (*it)->pe_adjusted << ") z=" << (*it)->z << std::endl;

    for ( KPPulseListIter it=pulselist.begin(); it!=pulselist.end(); it++ ) {
      ttrig.push_back( (*it)->tstart );
      tpeak.push_back( (*it)->tpeak );
      tend.push_back( (*it)->tend );
      peakamp.push_back( (*it)->peakamp );
      pulseperaw.push_back( (*it)->pe ); 
      pulsepedark.push_back( (*it)->pe_dark ); 
      pulsepe.push_back( (*it)->pe_adjusted ); 
      pulsez.push_back( (*it)->z ); 
      delete *it;
      *it = NULL;
    }
    pulselist.clear();

    // VETO: LOOKING FOR SMALLER PULSE, SEARCH, HOOP BY HOOP
    double vetothreshold = 10;
    double vetoerr = 1.;

    std::cout << "  OD pulses (expected dark rate=" << nod_sipms_per_hoop*(sipm_darkrate_hz*1.0e-9)*window_ns_veto << ")" << std::endl;

    std::vector<double> temp_twfm_veto;
    for (int ihoop=900; ihoop<900+102; ihoop++) {

      int ihoop_end = ihoop;
      if ( ihoop<1000 ) {
	if ( ihoop_end>=1000 )
	  ihoop_end = ihoop_end=999;
      }
      else {
	ihoop_end = ihoop;
      }

      if ( ihoop<1000) {
	vetothreshold = nod_sipms_per_hoop*((ihoop_end-ihoop+1))*(sipm_darkrate_hz*1.0e-9)*window_ns_veto; // side
      }
      else {
	vetothreshold = 100*(sipm_darkrate_hz*1.0e-9)*window_ns_veto; // endcap
	//std::cout << "vetothresh: " << vetothreshold << std::endl;
      }

      if ( vetothreshold<1.0 ) {
	if ( trig_version<=2 )
	  vetoerr = 5.0;
	else
	  vetoerr = 7.0;
      }
      else if ( vetothreshold<5.0 )
	vetoerr = 8.0*vetothreshold;
      else {
	// big threshold
	if ( trig_version==3 ) {
	  if ( ihoop<1000 )
	    vetoerr = 17.0;
	  else
	    vetoerr = 25.0;
	}
	else {
	  vetoerr = 5.0*sqrt(vetothreshold); // normal
	}
      }

//       if ( ihoop>=1000 )
// 	std::cout << "vetothresh: " << vetothreshold << " +/- " << vetoerr << ": " << (sipm_darkrate_hz*1.0e-9)*window_ns_veto << std::endl;
      vetothreshold += vetoerr;

      if ( sipm_darkrate_hz<=0 )
	vetothreshold = 0.5;

      int npulses_vetohoop = find_trigger( mc, 
					   vetothreshold, window_ns_veto, sipm_darkrate_hz,
					   true, ihoop, ihoop_end,
					   false, 0, 0.0,
					   n_decay_constants_veto, decay_weights_veto, decay_constants_ns_veto,
					   pulselist_veto, 90000, true, temp_twfm_veto, trig_version );
      npulses_veto += npulses_vetohoop;
      
      if ( npulses_vetohoop>0 ) {
	assign_pulse_charge( mc, pmtinfofile, pulselist_veto,
			     sipm_darkrate_hz,
			     true, ihoop, ihoop,
			     50.0, 90000, nod_sipms_per_hoop, nod_sipms_per_hoop_endcap,
			     true, trig_version );
      }
      double odintegral=0.0;
      for ( std::vector<double>::iterator od_it=temp_twfm_veto.begin(); od_it!=temp_twfm_veto.end(); od_it++)
	odintegral += *od_it; 
      float odpos[3];
      //pmtinfo->getposition( 90000 + (ihoop-900), odpos );
      int odpmt_fromhoop;
      if ( trig_version<=2 ) {
	odpmt_fromhoop = 90000 + (ihoop-900)*10;
      }
      else {
	if ( ihoop<1000 )
	  odpmt_fromhoop = 90000 + (ihoop-900)*50;
	else if ( ihoop==1000 )
	  odpmt_fromhoop = 95000;
	else if ( ihoop==1001 )
	  odpmt_fromhoop = 95100;
	else
	  assert(false);
      }

      pmtinfo->getposition( odpmt_fromhoop, odpos );
      if ( npulses_vetohoop>0 || ( sipm_darkrate_hz==0 && odintegral>0 ) ) {
	std::cout << "    -- od hoopid: [" << ihoop << "," << ihoop_end << "]"
		  << ": z=" << odpos[2] 
		  << " integral=" << odintegral 
		  << " dark rate in window=" << nod_sipms_per_hoop*((ihoop_end-ihoop+1))*(sipm_darkrate_hz*1.0e-9)*window_ns_veto << " +/- " << vetoerr
		  << " vetothreshold=" << vetothreshold << " (ave=" << vetothreshold/window_ns_veto << ")"
		  << " npulses=" << npulses_vetohoop
		  << std::endl;
	for (int iod=0; iod<npulses_vetohoop; iod++) {
	  std::cout << "     odpulse " << iod << ": "
		    << " pe=" << pulselist_veto.at(iod)->pe_adjusted
		    << " t=(" << pulselist_veto.at(iod)->tstart << ", " << pulselist_veto.at(0)->tpeak << ", " << pulselist_veto.at(0)->tend << ")"
		    << " z=" << pulselist_veto.at(iod)->z
		    << std::endl;
	  pulse_totodpe += pulselist_veto.at(iod)->pe_adjusted;
	  pulsepe_veto.push_back( pulselist_veto.at(iod)->pe_adjusted );
	  pulsez_veto.push_back(  pulselist_veto.at(iod)->z );
	  ttrig_veto.push_back(  pulselist_veto.at(iod)->tstart );
	  tend_veto.push_back(  pulselist_veto.at(iod)->tend );
	}
      }
      if ( ihoop>=1000 && npulses_vetohoop>0 ) {
	twfm_veto = temp_twfm_veto;
	std::cout << "save veto wfm: " << ihoop << std::endl;
      }
      pulselist_veto.clear();
    }//hoop loop
    std::cout << "  Total OD pulse pe: " << pulse_totodpe << " (vs. pre pe " << predark_odpe <<")" << std::endl;

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
