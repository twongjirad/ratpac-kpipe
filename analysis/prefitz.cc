#include "prefitz.hh"
#include "pmtinfo.hh"
#include <iostream>

// Goal is to find z of pulse.  Then we can limit a range of the detector to look over to reduce dark noise.
double calc_prefitz( RAT::DS::MC* mc, std::string pmtinfofile, double darkrate_hz, double window_ns, int n_id_sipms, int n_od_sipms, int& maxhoop ) {

  PMTinfo* pmtinfo = PMTinfo::GetPMTinfo( pmtinfofile );

  int nhoops = 1000;
  double hoop_totals[1000] = { 0 };
  double hoop_totals_adjusted[1000] = { 0 };
  double hoop_z[1000] = { 0.0 };
  double ndark_expt = darkrate_hz*1.0e-9*window_ns;

  for (int ipmt=0; ipmt<mc->GetMCPMTCount(); ipmt++) {

    RAT::DS::MCPMT* pmt = mc->GetMCPMT( ipmt );
    int pmtid = pmt->GetID();

    if ( pmtid>=n_id_sipms )
      continue;

    int hoopid = pmtid/100;

    double npe = 0.0;
    for (int iphoton=0; iphoton<pmt->GetMCPhotonCount(); iphoton++) {
      RAT::DS::MCPhoton* photon = pmt->GetMCPhoton( iphoton );
      if ( photon->GetHitTime()<window_ns )
	npe += photon->GetCharge();
    }
      
    float pmtpos[3];
    pmtinfo->getposition( pmtid, pmtpos );
    
    hoop_z[hoopid] = pmtpos[2]/100.0;
    hoop_totals[hoopid] += npe;
    hoop_totals_adjusted[hoopid] += npe - (ndark_expt);

  }

  //std::cout << "HOOP TOTALS" << std::endl;
  double z_hoop = 0.0;
  double hooppe = 0.0;
  maxhoop = -1;
  double maxhoop_pe = 0.;
  for (int ihoop=0; ihoop<nhoops; ihoop++) {
    //std::cout << " [hoop " << ihoop << ", z=" << hoop_z[ihoop] << "] " << hoop_totals[ihoop] << ", " << hoop_totals_adjusted[ihoop] << std::endl;
    z_hoop += hoop_z[ihoop]*(hoop_totals[ihoop]*hoop_totals[ihoop]);
    hooppe += hoop_totals[ihoop]*hoop_totals[ihoop];

    // find max hoop
    if ( maxhoop_pe<hoop_totals[ihoop] ) {
      maxhoop = ihoop;
      maxhoop_pe = hoop_totals[ihoop];
    }
  }
    
  z_hoop /= hooppe;

  //std::cout << "Max hoop: " << maxhoop << " z=" << hoop_z[maxhoop] << std::endl;
  //std::cout << "PE-weighted hoop z: " << z_hoop << std::endl;

  if ( hooppe>0.0 ) {
    return hoop_z[maxhoop];
  }
  else {
    return 0.0;
  }
  
}
