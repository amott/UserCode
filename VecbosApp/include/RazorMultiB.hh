//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorMultiB_h
#define RazorMultiB_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

using namespace std;

class RazorMultiB : public Vecbos{
public:
  RazorMultiB(TTree *tree=0); /// Class Constructor
  RazorMultiB(TTree *tree=0, string jsonFile=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  RazorMultiB(TTree *tree=0, string jsonFile=string("none"), bool goodRunLS=false, bool isData=false, string smsName=string("none")); /// Class Constructor
  virtual ~RazorMultiB();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetWeight(double);
  double _weight;
  vector<Jet> FastJetAlgorithmForceTwo(vector<TLorentzVector> InputCollection, double Rparam=0.5, double thePtMin=1.);
  vector<TLorentzVector> CombineJetsTs(vector<TLorentzVector> myjets, vector<TLorentzVector> Ts);
  double shat3D(TLorentzVector TOT, TVector3 P);
  TVector3 Boost_type1(TLorentzVector H1, TLorentzVector H2);
  TVector3 BetaTR(TLorentzVector TOT, TVector3 MET);
  double shatR(TLorentzVector TOT, TVector3 MET);
  vector<float> ParseEvent();
  bool goodJetID(int i);
  bool FailFilters();
  bool isFlagged();
  void InitEventFlag(char *s_Event);


private:
  bool _isSMS;
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  struct JetConfig;
  struct EventIndex;
  map<EventIndex, int> EventCounts;
  JetConfig *theJetConfig;
  
};
#endif
