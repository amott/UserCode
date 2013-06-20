 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
// Alex Mott (Caltech)
//
//-------------------------------------------------------

#ifndef VecbosExample_h
#define VecbosExample_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class VecbosExample : public Vecbos{
public:

  VecbosExample(TTree *tree=0); /// Class Constructor
  VecbosExample(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false,int mod=-1); /// Class Constructor
  virtual ~VecbosExample();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
private:

  vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

};
#endif
