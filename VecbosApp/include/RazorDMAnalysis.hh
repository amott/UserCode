//-------------------------------------------------------
// Description:
//    Razor analysis for 2b2mu+MET final state
// Authors:
//
//-------------------------------------------------------

#ifndef RazorDMAnalysis_h
#define RazorDMAnalysis_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class RazorDMAnalysis : public Vecbos{
public:

  RazorDMAnalysis(TTree *tree=0); /// Class Constructor
  RazorDMAnalysis(TTree *tree=0, string jsonFile=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorDMAnalysis();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetWeight(double);
  double _weight;

private:
  int HighestPt(vector<TLorentzVector> p, int iHIGHEST);
  bool _isSMS;
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  
};
#endif
