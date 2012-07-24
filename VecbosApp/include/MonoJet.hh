//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef MonoJet_h
#define MonoJet_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class MonoJet : public Vecbos{
public:

  MonoJet(TTree *tree=0); /// Class Constructor
  MonoJet(TTree *tree=0, string jsonFile=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~MonoJet();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetWeight(double);
  double _weight;

private:
  bool _isSMS;
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  int HighestPt(vector<TLorentzVector> v, int ignore);
  
};
#endif
