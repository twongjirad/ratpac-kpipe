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
#include "kpdaq.h"

//#define __CH_VERBOSE__
//#define __VETO_VERBOSE__

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
  double window_ns = 50.0;
  double window_ns_veto = 40.0;
  double decay_weights[2] = { 0.6, 0.4 };
  double decay_constants_ns[2] = { 45.0, 67.6 };
  int n_decay_constants_veto = 1;
  double decay_weights_veto[1] = { 1.0 };
  double decay_constants_ns_veto[1] = { 50.0 };
  int trig_version = 4;
  int nod_sipms_per_hoop_endcap = 100;
  int nodpmts[5] = {0,1200,1200,5200,10200}; 
  int nodhoops[5] = { 0,102,102,102,102 };
  int nodsipms_per_hoop[5] = { 0, 10, 10, 50, 100 };
  const int NPMTS = 90000+nodpmts[trig_version];

  int num_id_hoop_neighbors = 4;
  int num_id_hoop_sum = 100;

  double sipm_darkrate_hz = 1.6e6;
  double threshold = 500.0;
  //double sipm_darkrate_hz = 0.0;
  //double threshold = 10.0;


  // --------------------------------
  // Build DAQ
  KPDAQ daq( 100, nodsipms_per_hoop[trig_version], 900, nodhoops[trig_version], 1, trig_version, pmtinfofile );

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
  double twfm_integral = 0;
  double twfm_veto_integral = 0;
  // inner pipe vars
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
  tree->Branch( "twfm_integral", &twfm_integral, "twfm_integral/D" );
  tree->Branch( "twfm_veto_integral", &twfm_veto_integral, "twfm_veto_integral/D" );
  //tree->Branch( "twfm", &twfm );
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


  int ievent = 8;
  int nevents = ds->GetTotal();
  nevents = 11;

  KPPulseList pulselist;
  KPPulseList pulselist_veto;

  std::cout << "Number of events: " << nevents << std::endl;

  twfm.reserve(10000);
  twfm_veto.reserve(10000);
  
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
    twfm_integral = 0;
    twfm_veto_integral = 0;
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
    // LOAD DAQ
    daq.processEvent( *mc );

    // --------------------------------
    // TRIGGER
    std::cout << "  -- pulse finder -- " << std::endl;
    int ncoinhoops = 2;
    int nhoops_group = 2*int(ncoinhoops/2)+1;
    double expected_darkrate = 100*(sipm_darkrate_hz*1.0e-9)*window_ns*( nhoops_group ); // in channel
    double sig_darkrate = sqrt( expected_darkrate );
    //if ( expected_darkrate>10.0 )
    threshold = expected_darkrate + 5.0*sig_darkrate;
//     else
//       threshold = 10;
    if ( sipm_darkrate_hz== 0 )
      threshold = 0.5;
    std::cout << "Total hoops in a group that are searched: " << 
    std::cout << "Expected Dark Rate: " << expected_darkrate << " +/-  " << sig_darkrate << std::endl;
    std::cout << "Threshold set to " << threshold << std::endl;

    // INNER PIPE
    //  1) find pulses for all channels
    //  2) assemble pulses into one master waveform.  Run pulse finder again.
    std::map< int, KPPulseList* > idch_pulse_list;
    std::vector<double> tmp_twfm_id(10000, 0.0);

    // assemble list of pulses from channels
    int maxhoop = -1;
    double maxhoop_pe = 0;
    int npre_pulses = 0;
    for (int ich=1; ich<daq.getNIDChannels()-1; ich ++) {
      idch_pulse_list[ich] = new KPPulseList;
      int up_ch = ich-ncoinhoops/2;
      int ds_ch = ich+ncoinhoops/2;
      if ( up_ch<0 )
	up_ch = 0;
      if ( ds_ch >=daq.getNIDChannels() )
	ds_ch = daq.getNIDChannels()-1;
      int ch_npulses = find_trigger2( daq, 
				      threshold, window_ns, sipm_darkrate_hz,
				      up_ch, ds_ch,
				      false, 0, 0,
				      n_decay_constants, decay_weights, decay_constants_ns,
				      *(idch_pulse_list[ich]), 90000, false, tmp_twfm_id, trig_version );
#ifdef __CH_VERBOSE__
      std::cout << "[ channel " << ich << "]: ";
      std::cout << " npulses=" << ch_npulses;
#endif
      double chpe = 0.0;
      int ich_p = 1;
      for ( KPPulseListIter itp=idch_pulse_list[ich]->begin(); itp!=idch_pulse_list[ich]->end(); itp++ ) {
	double ppe = 0;
	double ppe_dark = 0;
#ifdef __CH_VERBOSE__
	std::cout << " (" << ich_p << ") ";
#endif
	int iend = (*itp)->tstart+window_ns;
	if ( iend>=10000 )
	  iend = 9999;
	for (int ibin= (*itp)->tstart; ibin< iend; ibin++ ) {
	  ppe += tmp_twfm_id.at( ibin ); 
	  ppe_dark += (sipm_darkrate_hz*1.0e-9)*100*ncoinhoops;
	}
	chpe += ppe-ppe_dark;
#ifdef __CH_VERBOSE__
	std::cout << " tpeak=" << (*itp)->tpeak << " pe=" << ppe-ppe_dark  << " (dark=" << ppe_dark << "+/-" << sqrt(ppe_dark) << "), ";
#endif
	ich_p++;
      }

      npre_pulses += ch_npulses;
#ifdef __CH_VERBOSE__ 
      std::cout << ": total channel pe=" << chpe << std::endl;
#endif
    }
    std::cout << "Number of ID channel pulses: " << npre_pulses << std::endl;

    // merge regions of interested based on pulses into one waveform
    twfm.assign(10000, 0.0);
    for ( std::map< int, KPPulseList* >::iterator it=idch_pulse_list.begin(); it!=idch_pulse_list.end(); it++ ) {
      int ich = (*it).first;
      double chpe = 0.0;
      if ( (*it).second->size()>0 ) {
	for ( KPPulseListIter pit=(*it).second->begin(); pit!=(*it).second->end(); pit++ ) {

	  bool fill = false;
	  // require neighboor cooincidence
// 	  int jch = ich-1;
// 	  if ( ich==0 )
// 	    jch = ich+1;
// 	  KPPulseList* us_pulses = idch_pulse_list[jch];
// 	  for (int jpulse=0; jpulse<us_pulses->size(); jpulse++) {
// 	    if ( fabs( us_pulses->at(jpulse)->tstart-(*pit)->tstart )<5.0 ) {
// 	      fill = true;
// 	      break;
// 	    }
// 	  }
	  fill = true;
	  
	  if ( fill ) {
	    daq.addWaveform( twfm, ich, (*pit)->tstart, (*pit)->tend, sipm_darkrate_hz*1.0e-9*100  );
	    if ( (*pit)->peakamp > maxhoop_pe ) {
	      maxhoop_pe = (*pit)->peakamp;
	      maxhoop = ich;
	    }
	  }
	}//end of loop over channel pulses
      }//if pulses found
    }//end of loop over channels and pulselist

    // look for final pulses
    npulses = find_trigger3( twfm,
			     threshold, window_ns_veto, sipm_darkrate_hz,
			     false, 0, 0,
			     n_decay_constants, decay_weights, decay_constants_ns,
			     pulselist, 90000, false, trig_version );
    for ( std::vector<double>::iterator itwfm=twfm.begin(); itwfm!=twfm.end(); itwfm++ )
      twfm_integral += *itwfm;

    if ( npulses>0 ) {
      int min_hoopid = maxhoop - num_id_hoop_sum;
      int max_hoopid = maxhoop + num_id_hoop_sum;
      if ( min_hoopid<0 ) min_hoopid = 0;
      if ( max_hoopid>=daq.getNIDChannels() ) max_hoopid = daq.getNIDChannels()-1;

      assign_pulse_charge( mc, pmtinfofile, pulselist, 
			   sipm_darkrate_hz,
			   true, min_hoopid, max_hoopid,
			   60.0, 90000, nodsipms_per_hoop[trig_version], nod_sipms_per_hoop_endcap,
			   false, trig_version );
    }
    float maxhoop_pos[3] = { 0 };
    if ( npulses>0 )
      daq.getChannelPos( maxhoop, &maxhoop_pos[0] );
    std::cout << "  ID npulses=" << npulses << "  with threshold=" <<  threshold << " ( exp. dark rate=" << expected_darkrate << " +/- " << sig_darkrate << ") maxhoop=" << maxhoop << " maxhoop[z]=" << maxhoop_pos[2] << std::endl;
    
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

    // ==================================================================================================================
    // VETO: SAME STRATEGY AS SIGAL REGION
    int veto_ncoinhoops = 2; // we analyze hoop with neighboring hoops
    int veto_nhoops_group_sides = 2*int(ncoinhoops/2)+1; 
    int veto_nhoops_group_caps = 1; // except for caps, these analyzed alone
    double veto_expected_darkrate_sides = nodsipms_per_hoop[trig_version]*(sipm_darkrate_hz*1.0e-9)*window_ns_veto*( veto_nhoops_group_sides );
    double veto_sig_darkrate_sides = sqrt( veto_expected_darkrate_sides );
    double veto_threshold_sides = veto_expected_darkrate_sides + 4.0*veto_sig_darkrate_sides;
    double veto_expected_darkrate_caps = nod_sipms_per_hoop_endcap*(sipm_darkrate_hz*1.0e-9)*window_ns_veto*veto_nhoops_group_caps;
    double veto_sig_darkrate_caps = sqrt( veto_expected_darkrate_caps );
    double veto_threshold_caps = veto_expected_darkrate_caps + 4.0*veto_sig_darkrate_caps;

    std::cout << "  OD pulses [thresh: SIDES=" << veto_threshold_sides << ", CAPS=" << veto_threshold_caps << "] "
	      << " (expected dark rate:"
	      << " SIDES=" << veto_expected_darkrate_sides << " +/- " << veto_sig_darkrate_sides 
	      << " CAPS=" << veto_expected_darkrate_caps << " +/- " << veto_sig_darkrate_caps << " )" << std::endl;

    std::map< int, KPPulseList* > hoop_pulse_list;
    std::vector<double> temp_twfm_veto(10000, 0.0);

    // assemble list of pulses from channels
    for (int ich=veto_ncoinhoops/2; ich<daq.getNODChannels(); ich ++) {

      
      hoop_pulse_list[ich] = new KPPulseList;
      double vetothreshold;
      double veto_nhoops_group;
      int up_ch;
      int ds_ch;
      if ( ich<daq.getNODChannels()-2 ) {
	// sides
	up_ch = ich-int(veto_ncoinhoops/2);
	ds_ch = ich+int(veto_ncoinhoops/2);
	if ( ds_ch>=daq.getNODChannels()-2 ){
	  ds_ch = daq.getNODChannels()-3;
	  if ( ds_ch<ich ) ds_ch = ich;
	}
	vetothreshold = veto_threshold_sides;
	veto_nhoops_group = veto_nhoops_group_sides;
      }
      else {
	// caps
	up_ch = ds_ch = ich;
	vetothreshold = veto_threshold_caps;
	veto_nhoops_group = veto_nhoops_group_caps;
      }

      if ( sipm_darkrate_hz==0 )
	vetothreshold = 0.5;
      
      KPPulseList* mypulselist = hoop_pulse_list[ich];
      int npulses_vetohoop = find_trigger2( daq,
					    vetothreshold, window_ns_veto, sipm_darkrate_hz,
					    900+up_ch, 900+ds_ch,
					    false, 0, 0.0,
					    n_decay_constants_veto, decay_weights_veto, decay_constants_ns_veto,
					    *mypulselist, 90000, true, temp_twfm_veto, trig_version );

#ifdef __VETO_VERBOSE__
      std::cout << "[ veto channel " << ich << "]: (" << up_ch << ", " << ds_ch << ", " << vetothreshold << ") ";
#endif

      double wfm_integral = 0;
      for ( std::vector< double >::iterator itwfm=temp_twfm_veto.begin(); itwfm!=temp_twfm_veto.end(); itwfm++ )
	wfm_integral += (*itwfm);
#ifdef __VETO_VERBOSE__
      std::cout << " integral=" << wfm_integral << " ";
      
      std::cout << " npulses=" << npulses_vetohoop;      
#endif
      double chpe = 0.0;
      int ich_p = 1;
      for ( KPPulseListIter itp=mypulselist->begin(); itp!=mypulselist->end(); itp++ ) {
	double ppe = 0;
	double ppe_dark = 0;
#ifdef __VETO_VERBOSE__
	std::cout << " (" << ich_p << ") ";
#endif
	int iend = (*itp)->tstart+window_ns;
	if ( iend>=10000 )
	  iend = 9999;
	for (int ibin= (*itp)->tstart; ibin< iend; ibin++ ) {
	  ppe += temp_twfm_veto.at( ibin ); 
	  ppe_dark += (sipm_darkrate_hz*1.0e-9)*100*veto_nhoops_group;
	}
	chpe += ppe-ppe_dark;
#ifdef __VETO_VERBOSE__
	std::cout << " tpeak=" << (*itp)->tpeak << " pe=" << ppe-ppe_dark  << " (dark=" << ppe_dark << "+/-" << sqrt(ppe_dark) << "), ";
#endif
	ich_p++;
      }
#ifdef __VETO_VERBOSE__
      std::cout << ": total channel pe=" << chpe << std::endl;
#endif
    }//loop over veto channels

    // merge regions of interested based on pulses into one waveform
    twfm_veto.assign(10000, 0.0);
    double veto_maxhoop_pe = 0.0;
    int veto_maxhoop = -1;
    for ( std::map< int, KPPulseList* >::iterator it=hoop_pulse_list.begin(); it!=hoop_pulse_list.end(); it++ ) {
      int ich = (*it).first;

      double nvetosipms = nodsipms_per_hoop[trig_version];
      if ( ich>=daq.getNODChannels()-2 )
	nvetosipms = nod_sipms_per_hoop_endcap;

      double chpe = 0.0;
      if ( (*it).second->size()>0 ) {
	for ( KPPulseListIter pit=(*it).second->begin(); pit!=(*it).second->end(); pit++ ) {

	  bool fill = false;
	  // require neighboor cooincidence
// 	  int jch = ich-1;
// 	  if ( ich==0 )
// 	    jch = ich+1;
// 	  KPPulseList* us_pulses = idch_pulse_list[jch];
// 	  for (int jpulse=0; jpulse<us_pulses->size(); jpulse++) {
// 	    if ( fabs( us_pulses->at(jpulse)->tstart-(*pit)->tstart )<5.0 ) {
// 	      fill = true;
// 	      break;
// 	    }
// 	  }
	  fill = true;
	  
	  if ( fill ) {
	    daq.addWaveform( twfm_veto, 900+ich, (*pit)->tstart, (*pit)->tend, sipm_darkrate_hz*1.0e-9*nvetosipms  );
	    if ( (*pit)->peakamp > veto_maxhoop_pe ) {
	      veto_maxhoop_pe = (*pit)->peakamp;
	      veto_maxhoop = 900+ich;
	    }
	  }
	}//end of loop over channel pulses
      }//if pulses found
    }//end of loop over channels and pulselist
    
    // look for final pulses
    std::cout << "  post veto pulse finding: " << veto_threshold_sides << ", " << window_ns_veto << std::endl;
    npulses_veto = find_trigger3( twfm_veto,
				  veto_threshold_sides, window_ns_veto, sipm_darkrate_hz,
				  false, 0, 0,
				  n_decay_constants_veto, decay_weights_veto, decay_constants_ns_veto,
				  pulselist_veto, 90000, true, trig_version );
    for ( std::vector<double>::iterator itwfm=twfm_veto.begin(); itwfm!=twfm_veto.end(); itwfm++ )
      twfm_veto_integral += *itwfm;

    if ( npulses_veto>0 ) {
      int min_hoopid = veto_maxhoop - 5;
      int max_hoopid = veto_maxhoop + 5;
      if ( min_hoopid<0 ) min_hoopid = 0;
      if ( max_hoopid-daq.getNIDChannels()>=daq.getNODChannels() ) max_hoopid = daq.getNODChannels()-1;
      
      assign_pulse_charge( mc, pmtinfofile, pulselist_veto, 
			   sipm_darkrate_hz,
			   true, min_hoopid, max_hoopid,
			   60.0, 90000, nodsipms_per_hoop[trig_version], nod_sipms_per_hoop_endcap,
			   true, trig_version );
    }
//     float maxhoop_pos[3] = { 0 };
//     if ( npulses>0 )
//       daq.getChannelPos( maxhoop, &maxhoop_pos[0] );
    std::cout << "  OD npulses=" << npulses_veto << std::endl;
    
    pulse_totodpe = 0.0;
    for ( KPPulseListIter it=pulselist_veto.begin(); it!=pulselist_veto.end(); it++ ) {
      std::cout << "    - tstart=" << (*it)->tstart 
		<< " tpeak=" << (*it)->tpeak 
		<< " tend=" << (*it)->tend
		<< " pe=" << (*it)->pe << " (dark=" << (*it)->pe_dark << ", adjusted=" << (*it)->pe_adjusted << ") z=" << (*it)->z << std::endl;
      pulse_totodpe += (*it)->pe_adjusted;
      pulsepe_veto.push_back( (*it)->pe_adjusted );
      pulsez_veto.push_back(  (*it)->z );
      ttrig_veto.push_back(  (*it)->tstart );
      tend_veto.push_back(  (*it)->tend );
      delete *it;
      *it = NULL;
    }
    
    //std::cout << "  Total OD pulse pe: " << pulse_totodpe << " (vs. pre pe " << predark_odpe <<")" << std::endl;

    // ==================================================================================================================
    // CLEAN UP PULSE LISTS

    pulselist.clear();
    for ( std::map< int, KPPulseList* >::iterator it=idch_pulse_list.begin(); it!=idch_pulse_list.end(); it++ ) {
      free_pulse_list( *(*it).second );
    }    

    //free_pulse_list( pulselist_veto );
    pulselist_veto.clear(); // does not own
    for ( std::map< int, KPPulseList* >::iterator it=hoop_pulse_list.begin(); it!=hoop_pulse_list.end(); it++ ) {
      free_pulse_list( *(*it).second );
    }


    // ==================================================================================================================

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
