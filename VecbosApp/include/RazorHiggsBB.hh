//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorHiggsBB_h
#define RazorHiggsBB_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class RazorHiggsBB : public Vecbos{
public:

  RazorHiggsBB(TTree *tree=0); /// Class Constructor
  RazorHiggsBB(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorHiggsBB();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  vector<TLorentzVector> SetFirstHiggs(TLorentzVector myHiggs, int merged); 
  vector<TLorentzVector> SetSecondHiggs(TLorentzVector myHiggs, int merged); 

private:

  vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

  // Higgs variables 
  double PFH1Pt;
  double PFH1Eta;
  double PFH1Phi;
  double PFH1Mass;
  int MergedH1;

  double PFH2Pt;
  double PFH2Eta;
  double PFH2Phi;
  double PFH2Mass;
  int MergedH2;

};
#endif
