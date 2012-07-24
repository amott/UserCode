//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef GenZjets_h
#define GenZjets_h

#include "GenVecbos.h"
#include "Vecbos.hh"

using namespace std;

class GenZjets : public GenVecbos {
public:

  GenZjets(TTree *tree=0, double xsec=100., double lumi=100.); /// Class Constructor
  virtual ~GenZjets();     /// Class Destructor
  void Loop(string outFileName);
  void SetEtaMax(double max);
  void SetPtZMin(double max);

private:
  int numJetsAbovePtAndEta(std::vector<Jet> & theJets, double thePtCut, double theEtaCut);
  double _weight;
  double _EtaMax;
  double _PtZMin;
};
#endif
