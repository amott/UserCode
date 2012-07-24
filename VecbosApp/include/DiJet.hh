//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef DiJet_h
#define DiJet_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class DiJet : public Vecbos{
public:

  DiJet(TTree *tree=0); /// Class Constructor
  DiJet(TTree *tree=0, bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~DiJet();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);

private:
  int HighestPtJet(vector<TLorentzVector> Jet, int firstJet);
  vector<TLorentzVector> JetRecovery(vector<TLorentzVector> Jet, int iJ1, int iJ2, double dRRatio, double pTthreshold);
  bool _isData;
  bool _goodRunLS;

};
#endif
