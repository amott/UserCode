// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
// #include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "SUSYTau.hh"

const int numTrees = 7;

SUSYTau::SUSYTau(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
}

SUSYTau::SUSYTau(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

SUSYTau::~SUSYTau() {}

void SUSYTau::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void SUSYTau::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  // prescaled triggers for May10 data
  int HLT_IsoMu17;
  int HLT_Ele8_CaloIdL_CaloIsoVL;
  int HLT_Ele8_CaloIdL_TrkIdVL;
  int HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT;
  int HLT_Ele17_CaloIdL_CaloIsoVL;
  int HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;

  // prescaled Razor Triggers
  int HLT_R014_MR150;
  int HLT_R020_MR150;
  int HLT_R025_MR150;

  // hadronic razor triggers
  int HLT_R020_MR500;
  int HLT_R020_MR550;
  int HLT_R025_MR400;
  int HLT_R025_MR450;
  int HLT_R033_MR300;
  int HLT_R033_MR350;
  int HLT_R038_MR200;
  int HLT_R038_MR250;

  // Muon Razor Triggers
  int HLT_Mu8_R005_MR200;
  int HLT_Mu8_R020_MR200;
  int HLT_Mu8_R025_MR200;
  int HLT_Mu8_R029_MR200;
  int HLT_Mu10_R005_MR200; // these three are new!
  int HLT_Mu10_R025_MR200;
  int HLT_Mu10_R029_MR200;
  int HLT_Mu10_R033_MR200;
  
  // Ele Razor Triggers
  int HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200;
  int HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200;
  int HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200;
  int HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200;
  int HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200;
  int HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200;
  int HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200;
  int HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200;

  // general event info
  int run;
  long evNum;
  double bx;
  int ls;
  double orbit;
  int Nevents;
  int nPv;

  // PF block
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double PFR, PFRsq;
  double PFMR;

  // MC bkg flags
  int ZWlep;
  int TTlep;

  // Tau block
  double pTTau;
  double etaTau;
  double phiTau;
  int chargeTau;

  // Electron Block
  double pTEle;
  double etaEle;
  double phiEle;
  int chargeEl;

  // Muon Block
  double pTMu;
  double etaMu;
  double phiMu;
  int chargeMu;

  // b-tag
  int hasbtag;

  // // MET filters
  // int ECALTPFilterFlag;
  // int drBoundary;
  // int drDead;
  // int CSCHaloFilterFlag;
  // int trackerFailureFilterFlag;
  // int BEECALflag;

  // // Tracker failures flags
  // int tooManySeeds;
  // int tooManyClusters;

  // prepare the output tree
  TTree* outTree[numTrees]; 
  outTree[0] = new TTree("outTreeHad", "outTreeHad");
  outTree[1] = new TTree("outTreeEle", "outTreeEle");
  outTree[2] = new TTree("outTreeEleTauSS", "outTreeEleTauSS");
  outTree[3] = new TTree("outTreeEleTauOS", "outTreeEleTauOS");
  outTree[4] = new TTree("outTreeMu", "outTreeMu");
  outTree[5] = new TTree("outTreeMuTauSS", "outTreeMuTauSS");
  outTree[6] = new TTree("outTreeMuTauOS", "outTreeMuTauOS");

  for(int i=0; i<numTrees; i++) {
    outTree[i]->Branch("HLT_IsoMu17", &HLT_IsoMu17, "HLT_IsoMu17/I");
    outTree[i]->Branch("HLT_Ele8_CaloIdL_CaloIsoVL", &HLT_Ele8_CaloIdL_CaloIsoVL, "HLT_Ele8_CaloIdL_CaloIsoVL/I");
    outTree[i]->Branch("HLT_Ele8_CaloIdL_TrkIdVL", &HLT_Ele8_CaloIdL_TrkIdVL, "HLT_Ele8_CaloIdL_TrkIdVL/I");
    outTree[i]->Branch("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT", &HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT/I");
    outTree[i]->Branch("HLT_Ele17_CaloIdL_CaloIsoVL", &HLT_Ele17_CaloIdL_CaloIsoVL, "HLT_Ele17_CaloIdL_CaloIsoVL/I");
    outTree[i]->Branch("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", &HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL/I");
    outTree[i]->Branch("HLT_R014_MR150", &HLT_R014_MR150, "HLT_R014_MR150/I");
    outTree[i]->Branch("HLT_R020_MR150", &HLT_R020_MR150, "HLT_R020_MR150/I");
    outTree[i]->Branch("HLT_R025_MR150", &HLT_R025_MR150, "HLT_R025_MR150/I");
    outTree[i]->Branch("HLT_R020_MR500", &HLT_R020_MR500, "HLT_R020_MR500/I");
    outTree[i]->Branch("HLT_R020_MR550", &HLT_R020_MR550, "HLT_R020_MR550/I");
    outTree[i]->Branch("HLT_R025_MR400", &HLT_R025_MR400, "HLT_R025_MR400/I");
    outTree[i]->Branch("HLT_R025_MR450", &HLT_R025_MR450, "HLT_R025_MR450/I");
    outTree[i]->Branch("HLT_R033_MR300", &HLT_R033_MR300, "HLT_R033_MR300/I");
    outTree[i]->Branch("HLT_R033_MR350", &HLT_R033_MR350, "HLT_R033_MR350/I");
    outTree[i]->Branch("HLT_R038_MR200", &HLT_R038_MR200, "HLT_R038_MR200/I");
    outTree[i]->Branch("HLT_R038_MR250", &HLT_R038_MR250, "HLT_R038_MR250/I");
    outTree[i]->Branch("HLT_Mu8_R005_MR200", &HLT_Mu8_R005_MR200, "HLT_Mu8_R005_MR200/I");
    outTree[i]->Branch("HLT_Mu8_R020_MR200", &HLT_Mu8_R020_MR200, "HLT_Mu8_R020_MR200/I");
    outTree[i]->Branch("HLT_Mu8_R025_MR200", &HLT_Mu8_R025_MR200, "HLT_Mu8_R025_MR200/I");
    // outTree[i]->Branch("HLT_Mu8_R029_MR200", &HLT_Mu8_R029_MR200, "HLT_Mu8_R029_MR200/I"); // do I have this?
    outTree[i]->Branch("HLT_Mu10_R005_MR200", &HLT_Mu10_R005_MR200, "HLT_Mu10_R005_MR200/I");
    outTree[i]->Branch("HLT_Mu10_R025_MR200", &HLT_Mu10_R025_MR200, "HLT_Mu10_R025_MR200/I");
    outTree[i]->Branch("HLT_Mu10_R029_MR200", &HLT_Mu10_R029_MR200, "HLT_Mu10_R029_MR200/I");
    outTree[i]->Branch("HLT_Mu10_R033_MR200", &HLT_Mu10_R033_MR200, "HLT_Mu10_R033_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200", &HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200, "HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200", &HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200, "HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200", &HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200, "HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200", &HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200, "HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200", &HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200, "HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200/I");
    outTree[i]->Branch("HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200", &HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200, "HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200/I");
    outTree[i]->Branch("HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200", &HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200, "HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200/I");
    outTree[i]->Branch("HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200", &HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200, "HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200/I");
    outTree[i]->Branch("run", &run, "run/I");
    outTree[i]->Branch("evNum", &evNum, "evNum/L");
    outTree[i]->Branch("bx", &bx, "bx/D");
    outTree[i]->Branch("ls", &ls, "ls/I");
    outTree[i]->Branch("orbit", &orbit, "orbit/D");
    outTree[i]->Branch("Nevents", &Nevents, "Nevents/I");
    outTree[i]->Branch("nPV", &nPv, "nPV/I");

    // PF block
    outTree[i]->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
    outTree[i]->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
    outTree[i]->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
    outTree[i]->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
    outTree[i]->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
    outTree[i]->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
    outTree[i]->Branch("PFR", &PFR, "PFR/D");
    outTree[i]->Branch("PFRsq", &PFRsq, "PFRsq/D");
    outTree[i]->Branch("PFMR", &PFMR, "PFMR/D");
    outTree[i]->Branch("ZWlep", &ZWlep, "ZWlep/I");
    outTree[i]->Branch("TTlep", &TTlep, "TTlep/I");

    // b-tag
    outTree[i]->Branch("hasbtag", &hasbtag, "hasbtag/I");

    // // Filters' block
    // outTree[i]->Branch("ECALTPFilterFlag", &ECALTPFilterFlag, "ECALTPFilterFlag/I");
    // outTree[i]->Branch("drBoundary", &drBoundary, "drBoundary/I");
    // outTree[i]->Branch("drDead", &drDead, "drDead/I");
    // outTree[i]->Branch("CSCHaloFilterFlag", &CSCHaloFilterFlag, "CSCHaloFilterFlag/I");
    // outTree[i]->Branch("trackerFailureFilterFlag", &trackerFailureFilterFlag, "trackerFailureFilterFlag/I");
    // outTree[i]->Branch("BEECALflag", &BEECALflag, "BEECALflag/I");
    // outTree[i]->Branch("tooManySeeds", &tooManySeeds, "tooManySeeds/I");
    // outTree[i]->Branch("tooManyClusters", &tooManyClusters, "tooManyClusters/I");

    // Tau block
    if (i!=0 && i!=1 && i!=4) {
      outTree[i]->Branch("pTTau", &pTTau, "pTTau/D");
      outTree[i]->Branch("etaTau", &etaTau, "etaTau/D");
      outTree[i]->Branch("phiTau", &phiTau, "phiTau/D");
      outTree[i]->Branch("chargeTau", &chargeTau, "chargeTau/I");
    }
    // Ele block
    if (1==i || 2==i || 3==i) {
      outTree[i]->Branch("pTEle", &pTEle, "pTEle/D");
      outTree[i]->Branch("etaEle", &etaEle, "etaEle/D");
      outTree[i]->Branch("phiEle", &phiEle, "phiEle/D");
      outTree[i]->Branch("chargeEle", &chargeEl, "chargeEle/I");
    }
    // Muon block
    if (4==i || 5==i || 6==i) {
      outTree[i]->Branch("pTMu", &pTMu, "pTMu/D");
      outTree[i]->Branch("etaMu", &etaMu, "etaMu/D");
      outTree[i]->Branch("phiMu", &phiMu, "phiMu/D");
      outTree[i]->Branch("chargeMu", &chargeMu, "chargeMu/I");
    }
  }

  // double pT;
  // double eta;
  // double HPSloose;
  // double HPSMu;
  // double HPSEle;

  // TTree* tauTree = new TTree("tauTree", "tauTreeHad");
  // tauTree->Branch("pT", &pT, "pT/D");
  // tauTree->Branch("eta", &eta, "eta/D");
  // tauTree->Branch("HPSloose", &HPSloose, "HPSloose/D");
  // tauTree->Branch("HPSMu", &HPSMu, "HPSMu/D");
  // tauTree->Branch("HPSEle", &HPSEle, "HPSEle/D");

  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // prescaled Triggers for May10 data
  vector<string> maskHLT_IsoMu17; maskHLT_IsoMu17.push_back("HLT_IsoMu17_v");
  vector<string> maskHLT_Ele8_CaloIdL_CaloIsoVL; maskHLT_Ele8_CaloIdL_CaloIsoVL.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v");
  vector<string> maskHLT_Ele8_CaloIdL_TrkIdVL; maskHLT_Ele8_CaloIdL_TrkIdVL.push_back("HLT_Ele8_CaloIdL_TrkIdVL_v");
  vector<string> maskHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT; maskHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT.push_back("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v");
  vector<string> maskHLT_Ele17_CaloIdL_CaloIsoVL; maskHLT_Ele17_CaloIdL_CaloIsoVL.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v");
  vector<string> maskHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL; maskHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");

  // prescaled Razor Triggers
  vector<string> maskHLT_R014_MR150; maskHLT_R014_MR150.push_back("HLT_R014_MR150_v");
  vector<string> maskHLT_R020_MR150; maskHLT_R020_MR150.push_back("HLT_R020_MR150_v");
  vector<string> maskHLT_R025_MR150; maskHLT_R025_MR150.push_back("HLT_R025_MR150_v");

  // hadronic razor triggers
  vector<string> maskHLT_R020_MR500; maskHLT_R020_MR500.push_back("HLT_R020_MR500_v");
  vector<string> maskHLT_R020_MR550; maskHLT_R020_MR550.push_back("HLT_R020_MR550_v");
  vector<string> maskHLT_R025_MR400; maskHLT_R025_MR400.push_back("HLT_R025_MR400_v");
  vector<string> maskHLT_R025_MR450; maskHLT_R025_MR450.push_back("HLT_R025_MR450_v");
  vector<string> maskHLT_R033_MR300; maskHLT_R033_MR300.push_back("HLT_R033_MR300_v");
  vector<string> maskHLT_R033_MR350; maskHLT_R033_MR350.push_back("HLT_R033_MR350_v");
  vector<string> maskHLT_R038_MR200; maskHLT_R038_MR200.push_back("HLT_R038_MR200_v");
  vector<string> maskHLT_R038_MR250; maskHLT_R038_MR250.push_back("HLT_R038_MR250_v");

  // Ele Razor Triggers
  vector<string> maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200; maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200_v");
  vector<string> maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200; maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200_v");
  vector<string> maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200; maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200_v");
  vector<string> maskHLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200; maskHLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200.push_back("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200_v");
  vector<string> maskHLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200; maskHLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200.push_back("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200_v");
  vector<string> maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200; maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200.push_back("HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200_v");
  vector<string> maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200; maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200.push_back("HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200_v");
  vector<string> maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200; maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200.push_back("HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200_v");

  // Muon Razor Triggers
  vector<string> maskHLT_Mu8_R005_MR200; maskHLT_Mu8_R005_MR200.push_back("HLT_Mu8_R005_MR200_v");
  vector<string> maskHLT_Mu8_R020_MR200; maskHLT_Mu8_R020_MR200.push_back("HLT_Mu8_R020_MR200_v");
  vector<string> maskHLT_Mu8_R025_MR200; maskHLT_Mu8_R025_MR200.push_back("HLT_Mu8_R025_MR200_v");
  vector<string> maskHLT_Mu10_R005_MR200; maskHLT_Mu10_R005_MR200.push_back("HLT_Mu10_R005_MR200_v");
  vector<string> maskHLT_Mu10_R025_MR200; maskHLT_Mu10_R025_MR200.push_back("HLT_Mu10_R025_MR200_v");
  vector<string> maskHLT_Mu10_R029_MR200; maskHLT_Mu10_R029_MR200.push_back("HLT_Mu10_R029_MR200_v");
  vector<string> maskHLT_Mu10_R033_MR200; maskHLT_Mu10_R033_MR200.push_back("HLT_Mu10_R033_MR200_v");

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  Nevents = stop;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      // prescaled Triggers for May10 data
      setRequiredTriggers(maskHLT_IsoMu17); reloadTriggerMask(true); HLT_IsoMu17 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele8_CaloIdL_CaloIsoVL); reloadTriggerMask(true); HLT_Ele8_CaloIdL_CaloIsoVL = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele8_CaloIdL_TrkIdVL); reloadTriggerMask(true); HLT_Ele8_CaloIdL_TrkIdVL = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT); reloadTriggerMask(true); HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele17_CaloIdL_CaloIsoVL); reloadTriggerMask(true); HLT_Ele17_CaloIdL_CaloIsoVL = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL); reloadTriggerMask(true); HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL = hasPassedHLT();
      // prescaled Razor Triggers
      setRequiredTriggers(maskHLT_R014_MR150); reloadTriggerMask(true); HLT_R014_MR150 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R020_MR150); reloadTriggerMask(true); HLT_R020_MR150 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R025_MR150); reloadTriggerMask(true); HLT_R025_MR150 = hasPassedHLT();
      // hadronic razor triggers
      setRequiredTriggers(maskHLT_R020_MR500); reloadTriggerMask(true); HLT_R020_MR500 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R020_MR550); reloadTriggerMask(true); HLT_R020_MR550 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R025_MR400); reloadTriggerMask(true); HLT_R025_MR400 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R025_MR450); reloadTriggerMask(true); HLT_R025_MR450 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R033_MR300); reloadTriggerMask(true); HLT_R033_MR300 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R033_MR350); reloadTriggerMask(true); HLT_R033_MR350 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R038_MR200); reloadTriggerMask(true); HLT_R038_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R038_MR250); reloadTriggerMask(true); HLT_R038_MR250 = hasPassedHLT();      
      // Muon Razor Triggers
      setRequiredTriggers(maskHLT_Mu8_R005_MR200); reloadTriggerMask(true); HLT_Mu8_R005_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Mu8_R020_MR200); reloadTriggerMask(true); HLT_Mu8_R020_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Mu8_R025_MR200); reloadTriggerMask(true); HLT_Mu8_R025_MR200 = hasPassedHLT();   
      setRequiredTriggers(maskHLT_Mu10_R005_MR200); reloadTriggerMask(true); HLT_Mu10_R005_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Mu10_R025_MR200); reloadTriggerMask(true); HLT_Mu10_R025_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Mu10_R029_MR200); reloadTriggerMask(true); HLT_Mu10_R029_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Mu10_R033_MR200); reloadTriggerMask(true); HLT_Mu10_R033_MR200 = hasPassedHLT();
      // Ele Razor Triggers
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200); reloadTriggerMask(true); HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200); reloadTriggerMask(true); HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200); reloadTriggerMask(true); HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200); reloadTriggerMask(true); HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200); reloadTriggerMask(true); HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200); reloadTriggerMask(true); HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200); reloadTriggerMask(true); HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200); reloadTriggerMask(true); HLT_Ele12_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R033_MR200 = hasPassedHLT();
    }
    
    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    // // HCAL noise filter
    if(!eventPassHcalFilter()) continue;

    if (_isData) {
      // MET filters
      if ( (METFlags >> 0)%2 == 0) continue;
      if ( (METFlags >> 1)%2 == 0) continue;
      if ( (METFlags >> 2)%2 == 0) continue;
      if ( (METFlags >> 3)%2 == 0) continue;
      if ( (METFlags >> 4)%2 == 0) continue;
      if ( (METFlags >> 5)%2 == 0) continue;
    
      // Tracker failures
      if ( (tooManyTrackerFailures >> 0)%2 == 1) continue;
      if ( (tooManyTrackerFailures >> 1)%2 == 1) continue;
    }

    if (!_isData && nMc>0) {
    // ------------------------
    // W -> enu, mu nu, tau nu
    // Z -> ee, nu nu, tau tau
    // tt -> WW
    //--------------------------
      ZWlep = ZWtype(nMc);
      TTlep = TTtype(nMc);
    }

    // Find the best Tight Muon
    int nmuTight = 0;
    TLorentzVector myMuTight(0.1, 0., 0., 0.1);
    int iMuTight = -99;
    for(int i=0; i<nMuon; i++) {
      if (muonPassTight(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 14. && fabs(thisMu.Eta())<2.1) {
	  nmuTight++;
	  if(thisMu.Pt() > myMuTight.Pt()) {
	    myMuTight = thisMu;
	    iMuTight = i;
	  }
	}
      }
    }

    // Find the best WP80 Electron
    int neleWP80 = 0;
    TLorentzVector myEleWP80(0.1, 0., 0., 0.1);
    int iEleWP80 = -99;
    for(int i=0; i<nEle; i++) {
      if(electronPassWP80(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 14. && fabs(thisEle.Eta()<3.0)) {
	  neleWP80++;
	  if(thisEle.Pt() > myEleWP80.Pt()) {
	    myEleWP80 = thisEle;
	    iEleWP80 = i;
	  }
	}
      }
    }

    hasbtag = 0;
    bool badjet = false;

    // Jet selection : these go into the hemispheres
    vector<TLorentzVector> PFPUcorrJet;    
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      if (pfJetPassTCHEM(i)==true)
	hasbtag = 1;

      // double scale = 1.;
      // if (_isData)
      // 	scale = L2L3CorrEnergyAK5PFPUcorrJet[i]/energyAK5PFPUcorrJet[i];
      double EU = uncorrEnergyAK5PFPUcorrJet[i];
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);
      // Apply jet correction first
      bool good_jet = false;
      if (myJet.Pt() > 40.0 && fabs(myJet.Eta()) < 3.0) {
	double fHAD = (neutralHadronEnergyAK5PFPUcorrJet[i]+chargedHadronEnergyAK5PFPUcorrJet[i])/EU;
	if(fHAD > 0.99) {
	  badjet = true;
	  break;
	}
	if (fabs(myJet.Eta()) < 2.4 && chargedEmEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (fabs(myJet.Eta()) < 2.4 && chargedHadronEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (muonEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (fabs(myJet.Eta()) >= 2.4 && neutralHadronEnergyAK5PFPUcorrJet[i]/EU < 0.99 && neutralEmEnergyAK5PFPUcorrJet[i]/EU < 0.99) good_jet = true;
	if (good_jet==false) {
	  badjet = true;
	  break;
	}
      }

      if (nmuTight>=1 && 0==neleWP80) {
       	if (myJet.Pt()>40. && fabs(myJet.Eta())< 2.4 && myJet.DeltaR(myMuTight)>0.3 )
       	  PFPUcorrJet.push_back(myJet);
      } else if (myJet.Pt()>40. && fabs(myJet.Eta())< 2.4)
       	PFPUcorrJet.push_back(myJet);
    }

    // Look for two 60 GeV/c jets
    vector<TLorentzVector> HighPtPFJet;
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      double EU = uncorrEnergyAK5PFPUcorrJet[i];
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);
      // Apply jet correction first
      bool good_jet = false;
      if (myJet.Pt() > 40.0 && fabs(myJet.Eta()) < 3.0) {
	double fHAD = (neutralHadronEnergyAK5PFPUcorrJet[i]+chargedHadronEnergyAK5PFPUcorrJet[i])/EU;
	if(fHAD > 0.99) {
	  badjet = true;
	  break;
	}
	if (fabs(myJet.Eta()) < 2.4 && chargedEmEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (fabs(myJet.Eta()) < 2.4 && chargedHadronEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (muonEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (fabs(myJet.Eta()) >= 2.4 && neutralHadronEnergyAK5PFPUcorrJet[i]/EU < 0.99 && neutralEmEnergyAK5PFPUcorrJet[i]/EU < 0.99) good_jet = true;
	if (good_jet==false) {
	  badjet = true;
	  break;
	}
      }

      if (nmuTight>=1 && 0==neleWP80) {
       	if (myJet.Pt()>60. && fabs(myJet.Eta())< 2.4 && myJet.DeltaR(myMuTight)>0.3 )
       	  HighPtPFJet.push_back(myJet);
      } else if (myJet.Pt()>60. && fabs(myJet.Eta())< 2.4)
       	HighPtPFJet.push_back(myJet);
    }

    if (badjet == true) continue;
    if ((int)HighPtPFJet.size()<2) continue;

    // Hemispheres
    pTPFHem1 = -9999;
    etaPFHem1 = -9999;
    phiPFHem1 = -9999;
    pTPFHem2 = -9999;
    etaPFHem2 = -9999;
    phiPFHem2 = -9999;
    PFR = -99999.;
    PFRsq = -99999.;
    PFMR = -99999.;
    TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);

    vector<TLorentzVector> tmpJet = CombineJets(PFPUcorrJet); //HighPtPFJet);
    if(tmpJet.size() >= 2) {
      TLorentzVector PFHem1 = tmpJet[0];
      TLorentzVector PFHem2 = tmpJet[1];      
      // compute boost
      // double num = PFHem1.P()-PFHem2.P();
      // double den = PFHem1.Pz()-PFHem2.Pz();            
      double MTR = CalcMTR(PFHem1, PFHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcGammaMRstar(PFHem1, PFHem2);
      if(variable >0) Rvariable = MTR/variable;      
      // fill the tree
      pTPFHem1 = PFHem1.Pt();
      etaPFHem1 = PFHem1.Eta();
      phiPFHem1 = PFHem1.Phi();
      pTPFHem2 = PFHem2.Pt();
      etaPFHem2 = PFHem2.Eta();
      phiPFHem2 = PFHem2.Phi();
      PFR = Rvariable;
      PFRsq = Rvariable*Rvariable;
      PFMR = variable;    
    }

    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;
    nPv = nPV;

    // // Find the best Tight Muon
    // int nmuTight = 0;
    // TLorentzVector myMuTight(0.1, 0., 0., 0.1);
    // int iMuTight = -99;
    // for(int i=0; i<nMuon; i++) {
    //   if(muonPassTight(i)) {
    // 	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
    // 	if(thisMu.Pt() > 10.) {
    // 	  nmuTight++;
    // 	  if(thisMu.Pt() > myMuTight.Pt()) {
    // 	    myMuTight = thisMu;
    // 	    iMuTight = i;
    // 	  }
    // 	}
    //   }
    // }

    // // Find the best WP80 Electron
    // int neleWP80 = 0;
    // TLorentzVector myEleWP80(0.1, 0., 0., 0.1);
    // int iEleWP80 = -99;
    // for(int i=0; i<nEle; i++) {
    //   if(electronPassWP80(i)) {
    // 	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
    // 	if(thisEle.Pt() > 12.) {
    // 	  neleWP80++;
    // 	  if(thisEle.Pt() > myEleWP80.Pt()) {
    // 	    myEleWP80 = thisEle;
    // 	    iEleWP80 = i;
    // 	  }
    // 	}
    //   }
    // }

    // Find the best HPS Tau
    int ntauHPS = 0;
    TLorentzVector myTauHPS(0.1, 0., 0., 0.1);
    int iTauHPS = -99;
    for (int i=0; i<nPFTau; i++) {
      if (abs(thehpsTauDiscrByTightElectronRejectionPFTau[i]-1.)<0.001 && abs(thehpsTauDiscrByTightMuonRejectionPFTau[i]-1.)<0.0001 &&
     	  abs(thehpsTauDiscrByLooseIsolationPFTau[i]-1.)<0.0001) {
     	TLorentzVector thisTau(pxPFTau[i], pyPFTau[i], pzPFTau[i], energyPFTau[i]);
     	if (thisTau.Pt() > 15. && fabs(thisTau.Eta())<2.3) {
     	  ntauHPS++;
     	  if (thisTau.Pt() > myTauHPS.Pt()) {
     	    myTauHPS = thisTau;
     	    iTauHPS = i;
     	  }
     	}
      }
    }

    // // Find the fake HPS Tau in the fit region -- to be commented when reconstructing taus
    // int ntauHPS = 0;
    // TLorentzVector myTauHPS(0.1, 0., 0., 0.1);
    // int iTauHPS = -99;
    // for (int i=0; i<nPFTau; i++) {
    //   if (abs(thehpsTauDiscrByTightElectronRejectionPFTau[i]-1.)<0.001 && abs(thehpsTauDiscrByTightMuonRejectionPFTau[i]-1.)<0.0001 &&
    // 	  abs(thehpsTauDiscrByLooseIsolationPFTau[i]-0.)<0.0001) {
    // 	TLorentzVector thisTau(pxPFTau[i], pyPFTau[i], pzPFTau[i], energyPFTau[i]);
    // 	if (thisTau.Pt() > 15. && PFMR>300 &&
    // 	    ( (PFMR<650 && PFRsq>=0.1 && PFRsq<0.2) ||
    // 	      (PFMR<400 && PFRsq>=0.2 && PFRsq<0.3) ||
    // 	      (PFMR<400 && PFRsq>=0.3 && PFRsq<0.45) ||
    // 	      (PFMR<350 && PFRsq>=0.45 && PFRsq<0.5) )
    // 	    ) {
    // 	  ntauHPS++;
    // 	  if (thisTau.Pt() > myTauHPS.Pt()) {
    // 	    myTauHPS = thisTau;
    // 	    iTauHPS = i;
    // 	  }
    // 	}
    //   }
    // }

    // // Fill the tauTree for the HPS efficiency
    // for (int i=0; i<nPFTau; i++) {
    //   pT = ptPFTau[i];      
    //   eta = etaPFTau[i];
    //   HPSloose = thehpsTauDiscrByLooseIsolationPFTau[i];
    //   HPSMu = thehpsTauDiscrByTightMuonRejectionPFTau[i];
    //   HPSEle = thehpsTauDiscrByLooseIsolationPFTau[i];
    //   if (pT>=15. && fabs(eta)<2.3 && nmuTight>=1 && 0==neleWP80)
    // 	tauTree->Fill();
    // }



    // Fill the tree per box
    // NOTE: WE HAVE NO MU-ELE BOX
    pTTau = -999.;
    etaTau = -999.;
    phiTau = -999.;
    chargeTau = -999;    
    pTMu = -999.;
    etaMu = -999.;
    phiMu = -999.;
    chargeMu = -999;
    pTEle = -999.;
    etaEle = -999.;
    phiEle = -999.;
    chargeEl = -999;
    if (neleWP80>=1 && 0==nmuTight) {
      pTEle = myEleWP80.Pt();
      etaEle = myEleWP80.Eta();
      phiEle = myEleWP80.Phi();
      chargeEl = chargeEle[iEleWP80];
      if (ntauHPS>=1) {
	pTTau = myTauHPS.Pt();
	etaTau = myTauHPS.Eta();
	phiTau = myTauHPS.Phi();
	chargeTau = chargePFTau[iTauHPS];
	if (chargeTau==chargeEl)
	  outTree[2]->Fill(); // EleTauSS box
	else
	  outTree[3]->Fill(); // EleTauOS box
      } else
	outTree[1]->Fill(); // Ele box
    } else if (nmuTight>=1 && 0==neleWP80) {
      pTMu = myMuTight.Pt();
      etaMu = myMuTight.Eta();
      phiMu = myMuTight.Phi();
      chargeMu = chargeMuon[iMuTight];
      if (ntauHPS>=1) {
	pTTau = myTauHPS.Pt();
	etaTau = myTauHPS.Eta();
	phiTau = myTauHPS.Phi();
	chargeTau = chargePFTau[iTauHPS];
	if (chargeTau==chargeMu)
	  outTree[5]->Fill(); // MuTauSS box
	else
	  outTree[6]->Fill(); // MuTauOS box
      } else
	outTree[4]->Fill(); // Mu box
    } else if (0==neleWP80 && 0==nmuTight)
      outTree[0]->Fill(); // Had box

  }

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  for(int i=0; i<numTrees; i++) 
    outTree[i]->Write();
  // tauTree->Write();
  file->Close();
}

vector<TLorentzVector> SUSYTau::CombineJets(vector<TLorentzVector> myjets)
{  
  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;
  bool foundGood = false;
  
  int N_comb = 1;
  for(int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }
  
  double M_min = 9999999999.0;
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
	j_temp1 += myjets[count];
      } else {
	j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }    
    double M_temp = j_temp1.M2()+j_temp2.M2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  // set masses to 0
  j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;  
}

int SUSYTau::TTtype(int nMc)
{
  Int_t Tfirst=0;
  Int_t Tsecond=0;
  Int_t TTlep=0;
  Int_t Tindex=0;
  Int_t firstMoth=0;

 for (int i=0; i<nMc; i++) {
   Int_t mothIndex = mothMc[i];
   if (mothIndex<=nMc && mothIndex>=0) {
     if ((abs(idMc[i])==11 || abs(idMc[i])==12) && abs(idMc[mothIndex])==24) {
       Tfirst=-1;
       Tindex=i;
       firstMoth=mothIndex;
       break;
     }
     if ((abs(idMc[i])==13 || abs(idMc[i])==14) && abs(idMc[mothIndex])==24) {
       Tfirst=-2;
       Tindex=i;
       firstMoth=mothIndex;
       break;
     }
     if ((abs(idMc[i])==15 || abs(idMc[i])==16) && abs(idMc[mothIndex])==24) {
       Tfirst=-3;
       Tindex=i;
       firstMoth=mothIndex;
       break;
     }
   }
 }

 if (Tfirst!=0) {
    for (int j=Tindex+1; j<nMc; j++) {
      Int_t mothIndex = mothMc[j];

      if (mothIndex<=nMc && mothIndex>=0) {
	if ((abs(idMc[j])==11 || abs(idMc[j])==12) && abs(idMc[mothIndex])==24 && mothIndex!=firstMoth) {
	  Tsecond=-1;
	  break;
	}
	if ((abs(idMc[j])==13 || abs(idMc[j])==14) && abs(idMc[mothIndex])==24 && mothIndex!=firstMoth) {
	  Tsecond=-2;
	  break;
	}
	if ((abs(idMc[j])==15 || abs(idMc[j])==16) && abs(idMc[mothIndex])==24 && mothIndex!=firstMoth) {
	  Tsecond=-3;
	  break;
	}
      }
    }
  }

 if (Tfirst==-1 && Tsecond==-1) // EE
    TTlep=-1;
  if (Tfirst==-2 && Tsecond==-2) // MM
    TTlep=-2;
  if (Tfirst==-3 && Tsecond==-3) // TT
    TTlep=-3;
  if ((Tfirst==-1 && Tsecond==-2) || (Tfirst==-2 && Tsecond==-1)) // EM
    TTlep=-4;
  if ((Tfirst==-1 && Tsecond==-3) || (Tfirst==-3 && Tsecond==-1)) // TE
    TTlep=-5;
  if ((Tfirst==-2 && Tsecond==-3) || (Tfirst==-3 && Tsecond==-2)) // TM
    TTlep=-6;

  if (Tfirst==-1 && Tsecond==0) // Single E
    TTlep=+1;
  if (Tfirst==-2 && Tsecond==0) // Single M
    TTlep=+2;
  if (Tfirst==-3 && Tsecond==0) // Single T
    TTlep=+3;

  return TTlep;
}

int SUSYTau::ZWtype(int nMc)
{
  int ZWlep=0;
  for(Int_t i=0; i<nMc; i++) {
    Int_t mothIndex = mothMc[i];
    if (mothIndex<=nMc && mothIndex>=0)
      {            
	if (abs(idMc[i])==11 && idMc[mothIndex]==23)
	  ZWlep=1;
	else if (abs(idMc[i])==13 && idMc[mothIndex]==23)
	  ZWlep=2;
	else if (abs(idMc[i])==15 && idMc[mothIndex]==23)
	  ZWlep=3;
	else if ((abs(idMc[i])==11 || abs(idMc[i])==12) && abs(idMc[mothIndex])==24)
	  ZWlep=-1;
	else if ((abs(idMc[i])==13 || abs(idMc[i])==14) && abs(idMc[mothIndex])==24)
	  ZWlep=-2;
	else if ((abs(idMc[i])==15 || abs(idMc[i])==16) && abs(idMc[mothIndex])==24)
	  ZWlep=-3;
      }
  }
  return ZWlep;
}
