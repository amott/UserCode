//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorLeptons_h
#define RazorLeptons_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class RazorLeptons : public Vecbos{
public:

  RazorLeptons(TTree *tree=0); /// Class Constructor
  RazorLeptons(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorLeptons();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);

private:

  vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets);
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

  void FirstLepton(TLorentzVector lep, int id);
  void SecondLepton(TLorentzVector lep, int id);

  // Lepton Block
  double pTLep1;
  double etaLep1;
  double phiLep1;
  int idLep1;

  double pTLep2;
  double etaLep2;
  double phiLep2;
  int idLep2;

};
#endif
