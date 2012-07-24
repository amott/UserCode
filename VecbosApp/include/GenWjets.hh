//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef GenWjets_h
#define GenWjets_h

#include "GenVecbos.h"
#include "Vecbos.hh"

using namespace std;

class GenWjets : public GenVecbos {
public:

  GenWjets(TTree *tree=0, double xsec=100., double lumi=100.); /// Class Constructor
  virtual ~GenWjets();     /// Class Destructor
  void Loop(string outFileName);
  void SetEtaMax(double max);
  void SetPtWMin(double max);

private:
  int numJetsAbovePtAndEta(std::vector<Jet> & theJets, double thePtCut, double theEtaCut);
  double _weight;
  double _EtaMax;
  double _PtWMin;

};
#endif
