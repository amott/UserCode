//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef SUSYTau_h
#define SUSYTau_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class SUSYTau : public Vecbos{
public:

  SUSYTau(TTree *tree=0); /// Class Constructor
  SUSYTau(TTree *tree=0, string json="none", bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~SUSYTau();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);


private:
  vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets); //, TLorentzVector DiTau);
  int ZWtype(int);
  int TTtype(int);
  // TLorentzVector CalcMHT(vector<TLorentzVector> myjets);
  // double CalcHT(vector<TLorentzVector> myjets);
  // int FindHighestTaNC(vector<int>FakeTaus);
  // int FindHighestPt(vector<int>GoodElectrons);
  // int FindHighestPtTau(vector<int>GoodTaus);
  // int FindMostIsolatedTau(vector<int>GoodTaus);
  // int HighestPtJet(vector<TLorentzVector> Jet, int firstJet);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;

};
#endif
