//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef Razor_h
#define Razor_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class Razor : public Vecbos{
public:

  Razor(TTree *tree=0); /// Class Constructor
  Razor(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~Razor();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);

private:

  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

};
#endif
