//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef SUSYMultiB_h
#define SUSYMultiB_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class SUSYMultiB : public Vecbos{
public:

  SUSYMultiB(TTree *tree=0); /// Class Constructor
  SUSYMultiB(TTree *tree=0, bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~SUSYMultiB();     /// Class Destructor
  int HighestPtJet(vector<TLorentzVector> Jet, int firstJet);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  bool isTightMuon(int);
  bool isGlobalMuon(int);
  bool isLooseMuon(int);
  bool is80Electron(int);
  bool is95Electron(int);
  void SetWeight(double);
  vector<TLorentzVector> CombineJets(vector<TLorentzVector>);
  void BubbleSort(vector<float> &);
  double _weight;

private:
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  
};
#endif
