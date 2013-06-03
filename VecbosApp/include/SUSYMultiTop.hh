//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef SUSYMultiTop_h
#define SUSYMultiTop_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "CommonTools/include/LumiReWeightingStandAlone.h"

using namespace std;

class SUSYMultiTop : public Vecbos{
public:

  SUSYMultiTop(TTree *tree=0, double weight=1.0, int tWfile=0); /// Class Constructor
  SUSYMultiTop(TTree *tree=0, bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~SUSYMultiTop();     /// Class Destructor
  int HighestPtJet(vector<TLorentzVector> Jet, int firstJet);
  int NewTag(bool isBTagged, int pdgIdPart, double Btag_SF, double Btag_eff, double Bmistag_SF, double Bmistag_eff);
  int NewTagPOG(bool isBTagged, int pdgIdPart, double Btag_SF, double Btag_eff, double Bmistag_SF, double Bmistag_eff);
  bool SwitchTag(bool isBTagged, double Btag_SF, double Btag_eff);
  bool SwitchTagPOG(bool isBTagged, double Btag_SF, double Btag_eff);
  int NearestInt(double num);
  double FindWeight(string, double);
  bool isLooseJetMva(float pt, float eta, float id);
  std::string FindSample(string fileWeight, double weight);
  vector <double> ReadBins(string _file, bool isVerbose);
  vector <double> GetSFLight(string _file, int Error);
  vector <double> GetSF(string _file, int Error, bool isC);
  vector <double> GetEff(string _file);
  void SetTags(int JetFlavor, double JetpT, double JetEta, int JetTag, int& JetTagNew, int& JetTagNewBUp, int& JetTagNewBDown, int& JetTagNewLUp, int& JetTagNewLDown, bool isPOG);  
  vector<int> SetTagsSMS(int JetFlavor, double JetpT, double JetEta, int JetTag, bool isPOG);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  //  bool isTightMuon(int);
  // bool isGlobalMuon(int);
  //bool isLooseMuon(int);
  //bool is80Electron(int);
  //bool is95Electron(int);
  void SetWeight(double);
  int JetFlavorFull(TLorentzVector);
  int JetFlavorEasy(TLorentzVector);
  vector<TLorentzVector> CombineJets(vector<TLorentzVector>);
  void BubbleSort(vector<float> &);
  double _weight;
  string _sample;
  vector <double> pTB;
  vector <double> pTL;
  vector <double> EtaB;
  vector <double> EtaL;
  vector <double> pTBFast;
  vector <double> pTLFast;
  vector <double> EtaBFast;
  vector <double> EtaLFast;
  vector <double> SFBv;
  vector <double> SFCvUp;
  vector <double> SFCvDown;
  vector <double> SFLv;
  vector <double> SFBvUp;
  vector <double> SFBvDown;
  vector <double> SFLvUp;
  vector <double> SFLvDown;
  vector <double> SFBvFast;
  vector <double> SFCvFast;
  vector <double> SFCvUpFast;
  vector <double> SFCvDownFast;
  vector <double> SFLvFast;
  vector <double> SFBvUpFast;
  vector <double> SFBvDownFast;
  vector <double> SFLvUpFast;
  vector <double> SFLvDownFast;
  vector <double> EffL;
  vector <double> EffB;
  vector <double> EffC;
  vector <double> EffLFast;
  vector <double> EffBFast;
  vector <double> EffCFast;

private:
  bool _isData;
  bool isT1tttt;
  bool isEffOnly;
  bool _goodRunLS;
  TTree* _treeCond;
  
};
#endif
