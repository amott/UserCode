//-------------------------------------------------------
// Description:
//    Class to test VBTF cuts efficiency
// Authors: Raffaele Tito D'Agnolo (SNS)
//
//-------------------------------------------------------

#ifndef VBTFLeptEff_h
#define VBTFLeptEff_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class VBTFLeptEff : public Vecbos{
public:

  VBTFLeptEff(TTree *tree=0); /// Class Constructor
  VBTFLeptEff(TTree *tree=0, bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~VBTFLeptEff();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  bool isGlobalMuon(int iMu);
  bool MatchEle(int iEle);

private:
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;

};
#endif
