#ifndef __KPDAQ__
#define __KPDAQ__

#include <string>
#include <vector>
#include "RAT/DS/MC.hh"

class PMTinfo;
class TRandom3;

class KPDAQ {

 public:
  
  KPDAQ( int nidsipms_perhoop, int nodsipms_perhoop, int nidhoops, int nodhoops, int nidchperhoop, int version, std::string pmtinfofilename  );
  ~KPDAQ();

  int fNIDSiPMs;
  int fNODSiPMs;
  int fNIDSiPMs_perhoop;
  int fNODSiPMS_perhoop;
  int fNIDhoops;
  int fNODhoops;
  int fNIDchperhoop;
  int fVersion;
  int fNbins;
  std::string pmtinfo_fname;

  int getNIDChannels() { return fNIDChannels; };
  int getNODChannels() { return fNODChannels; };
  //void getChannelPos( double* pos );

  const std::vector< double >& getIDwfm( int ch ) const { return fID_wfm.at(ch); };
  const std::vector< double >& getODwfm( int ch ) const { return fOD_wfm.at(ch); };
  double getWindowSum( int ch, double start, double end );
  void copyWaveforms( std::vector<double>& copy, int chstart, int chend ) const;

  void processEvent( RAT::DS::MC& mc );
  
  void setThinPE( bool set ) { fThinPE = set; };
  void setThinPEfactor( double factor ) { fThinPEfactor = factor; if ( fThinPEfactor>1.0 ) fThinPEfactor = 1.0; };

 protected:

  int getHoopIndexFromPMTID( int pmtid );
  int getChannelIndex( int ihoop, int pmtid );

  int fNIDChannels;
  int fNODChannels;

  std::vector< std::vector<double> > fID_wfm;
  std::vector< std::vector<double> > fOD_wfm;

  bool fThinPE;
  double fThinPEfactor;

  PMTinfo* fpmtinfo;
  TRandom3* frand;
};


#endif
