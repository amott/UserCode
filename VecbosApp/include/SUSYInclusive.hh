//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef SUSYInclusive_h
#define SUSYInclusive_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class SUSYInclusive : public Vecbos{
public:

  SUSYInclusive(TTree *tree=0); /// Class Constructor
  SUSYInclusive(TTree *tree=0, bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~SUSYInclusive();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);

private:
  double SumPt(int iMu);
  int EvCategory();
  bool HTtrigger(double HTmin);
  bool isTightMuon(int iMu);
  vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets);
  int HighestPtJet(vector<TLorentzVector> Jet, int firstJet);
  vector<TLorentzVector> JetRecovery(vector<TLorentzVector> Jet, int iJ1, int iJ2, int dRRatio, double pTthreshold);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;

};
#endif
