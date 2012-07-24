 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorBoostedTop_h
#define RazorBoostedTop_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class RazorBoostedTop : public Vecbos{
public:

  RazorBoostedTop(TTree *tree=0); /// Class Constructor
  RazorBoostedTop(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false,int mod=-1); /// Class Constructor
  virtual ~RazorBoostedTop();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetFirstTop(TLorentzVector myTop, int merged); 
  void SetSecondTop(TLorentzVector myTop, int merged); 

  void PARSE_EVENT();
private:

  vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

  // Top variables 
  double PFT1Pt;
  double PFT1Eta;
  double PFT1Phi;
  double PFT1Mass;
  int MergedT1;

  double PFT2Pt;
  double PFT2Eta;
  double PFT2Phi;
  double PFT2Mass;
  int MergedT2;

  int MODEL;
  double m0;
  double m12;
  double tanb;
  double A0;
  double mu;
};
#endif
