//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef GammaPlusJet_h
#define GammaPlusJet_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class GammaPlusJet : public Vecbos{
public:

  GammaPlusJet(TTree *tree=0); /// Class Constructor
  GammaPlusJet(TTree *tree=0, bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~GammaPlusJet();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetLuminosity(double lumi);
  void SetXsection(double xsec);

private:

  bool isElectron(int iPh);
  int HighestPtSC(vector<int> indPhotons);
  TTree* _treeCond;
  bool _isData;
  bool _goodRunLS;
  double _Lumi;
  double _xsec;

};
#endif
