//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorDiPhoton_h
#define RazorDiPhoton_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class RazorDiPhoton : public Vecbos{
public:

  RazorDiPhoton(TTree *tree=0); /// Class Constructor
  RazorDiPhoton(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorDiPhoton();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);

private:
  bool PhotonIdBarrel(int iPh, TString Selector);
  bool PhotonIdBarrelisEM(int i);
  bool PhotonIdBarrelisLoose(int i);
  bool PhotonIdBarrelisTight(int i);
  int  AntiPhotonIdBarrel(int iPh);
  bool IsPhotonBarrel(int iPh);

  vector<TLorentzVector> CombineJets_R(vector<TLorentzVector> myjets);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

};
#endif
