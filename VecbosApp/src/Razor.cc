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
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "Razor.hh"

Razor::Razor(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
}

Razor::Razor(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

Razor::~Razor() {}

void Razor::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void Razor::SetWeight(double weight) {
  _weight = weight;
}

void Razor::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  cout << "ciaoo" << endl;

  // Prescaled Jet Triggers
  int HLT_DiJetAve30;
  int HLT_DiJetAve60;
  int HLT_Jet30;
  int HLT_Jet60;

  // other hadronic triggers
  int HLT_HT300_MHT75;
  int HLT_Meff640;
  int HLT_HT250_AlphaT0p55;
  int HLT_HT300_AlphaT0p52;

  // Razor BJet Triggers
  int HLT_R014_MR150_CentralJet40_BTagIP;
  int HLT_R014_MR450_CentralJet40_BTagIP;
  int HLT_R020_MR350_CentralJet40_BTagIP;
  int HLT_R025_MR250_CentralJet40_BTagIP;

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
  
  // Ele Razor Triggers
  int HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200;
  int HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200;
  int HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200;

  // PF  block
  int    passedPF;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double PFR;
  double PFRsq;
  double PFMR;

  // calo block
  int    passedCalo;
  double pTCaloHem1;
  double etaCaloHem1;
  double phiCaloHem1;
  double pTCaloHem2;
  double etaCaloHem2;
  double phiCaloHem2;
  double CaloR;
  double CaloRsq;
  double CaloMR;

  // general event info
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double W;

  // Muon Block
  double pTMu1;
  double etaMu1;
  double phiMu1;
  double chargeMu1;
  double pTMu2;
  double etaMu2;
  double phiMu2;
  double chargeMu2;

  // Electron Block
  double pTEle1;
  double etaEle1;
  double phiEle1;
  double chargeEle1;
  double pTEle2;
  double etaEle2;
  double phiEle2;
  double chargeEle2;

  // prepare the output tree
  TTree* outTree[6]; 
  outTree[0] = new TTree("outTreeHad", "outTreeHad");
  outTree[1] = new TTree("outTreeEle", "outTreeEle");
  outTree[2] = new TTree("outTreeMu", "outTreeMu");
  outTree[3] = new TTree("outTreeMuMu", "outTreeMuMu");
  outTree[4] = new TTree("outTreeMuEle", "outTreeMuEle");
  outTree[5] = new TTree("outTreeEleEle", "outTreeEleEle");

  for(int i=0; i<6; i++) {

    outTree[i]->Branch("HLT_DiJetAve30", &HLT_DiJetAve30, "HLT_DiJetAve30/I");
    outTree[i]->Branch("HLT_DiJetAve60", &HLT_DiJetAve60, "HLT_DiJetAve60/I");
    outTree[i]->Branch("HLT_Jet30", &HLT_Jet30, "HLT_Jet30/I");
    outTree[i]->Branch("HLT_Jet60", &HLT_Jet60, "HLT_Jet60/I");
    outTree[i]->Branch("HLT_HT300_MHT75", &HLT_HT300_MHT75, "HLT_HT300_MHT75/I");
    outTree[i]->Branch("HLT_Meff640", &HLT_Meff640, "HLT_Meff640/I");
    outTree[i]->Branch("HLT_HT250_AlphaT0p55", &HLT_HT250_AlphaT0p55, "HLT_HT250_AlphaT0p55/I");
    outTree[i]->Branch("HLT_HT300_AlphaT0p52", &HLT_HT300_AlphaT0p52, "HLT_HT300_AlphaT0p52/I");
    outTree[i]->Branch("HLT_R014_MR150_CentralJet40_BTagIP", &HLT_R014_MR150_CentralJet40_BTagIP, "HLT_R014_MR150_CentralJet40_BTagIP/I");
    outTree[i]->Branch("HLT_R014_MR450_CentralJet40_BTagIP", &HLT_R014_MR450_CentralJet40_BTagIP, "HLT_R014_MR450_CentralJet40_BTagIP/I");
    outTree[i]->Branch("HLT_R020_MR350_CentralJet40_BTagIP", &HLT_R020_MR350_CentralJet40_BTagIP, "HLT_R020_MR350_CentralJet40_BTagIP/I");
    outTree[i]->Branch("HLT_R025_MR250_CentralJet40_BTagIP", &HLT_R025_MR250_CentralJet40_BTagIP, "HLT_R025_MR250_CentralJet40_BTagIP/I");
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
    outTree[i]->Branch("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200", &HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200, "HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200", &HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200, "HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200/I");
    outTree[i]->Branch("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200", &HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200, "HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200/I");

    outTree[i]->Branch("run", &run, "run/D");
    outTree[i]->Branch("evNum", &evNum, "evNum/D");
    outTree[i]->Branch("bx", &bx, "bx/D");
    outTree[i]->Branch("ls", &ls, "ls/D");
    outTree[i]->Branch("orbit", &orbit, "orbit/D");
    outTree[i]->Branch("W", &W, "W/D");
    
    // PF block
    outTree[i]->Branch("passedPF", &passedPF, "passedPF/I");
    outTree[i]->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
    outTree[i]->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
    outTree[i]->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
    outTree[i]->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
    outTree[i]->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
    outTree[i]->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
    outTree[i]->Branch("PFR", &PFR, "PFR/D");
    outTree[i]->Branch("PFRsq", &PFRsq, "PFRsq/D");
    outTree[i]->Branch("PFMR", &PFMR, "PFMR/D");
    
    // Calo block
    outTree[i]->Branch("passedCalo", &passedCalo, "passedCalo/I");
    outTree[i]->Branch("pTCaloHem1", &pTCaloHem1, "pTCaloHem1/D");
    outTree[i]->Branch("etaCaloHem1", &etaCaloHem1, "etaCaloHem1/D");
    outTree[i]->Branch("phiCaloHem1", &phiCaloHem1, "phiCaloHem1/D");
    outTree[i]->Branch("pTCaloHem2", &pTCaloHem2, "pTCaloHem2/D");
    outTree[i]->Branch("etaCaloHem2", &etaCaloHem2, "etaCaloHem2/D");
    outTree[i]->Branch("phiCaloHem2", &phiCaloHem2, "phiCaloHem2/D");
    outTree[i]->Branch("CaloR", &CaloR, "CaloR/D");
    outTree[i]->Branch("CaloRsq", &CaloRsq, "CaloRsq/D");
    outTree[i]->Branch("CaloMR", &CaloMR, "CaloMR/D");
    
    // First muon for Mu,  DiMu and MuEle
    if(i ==2 || i == 3 || i == 4) {
      outTree[i]->Branch("pTMu1", &pTMu1, "pTMu1/D");
      outTree[i]->Branch("etaMu1", &etaMu1, "etaMu1/D");
      outTree[i]->Branch("phiMu1", &phiMu1, "phiMu1/D");
      outTree[i]->Branch("chargeMu1", &chargeMu1, "chargeMu11/D");
    }
    // Second MU
    if(i == 3) {
      outTree[i]->Branch("pTMu2", &pTMu2, "pTMu2/D");
      outTree[i]->Branch("etaMu2", &etaMu2, "etaMu2/D");
      outTree[i]->Branch("phiMu2", &phiMu2, "phiMu2/D");
      outTree[i]->Branch("chargeMu2", &chargeMu2, "chargeMu21/D");
    }
    // First Ele for Ele, DiEle and MuEle
    if(i ==1 || i == 4 || i == 5) {
      outTree[i]->Branch("pTEle1", &pTEle1, "pTEle1/D");
      outTree[i]->Branch("etaEle1", &etaEle1, "etaEle1/D");
      outTree[i]->Branch("phiEle1", &phiEle1, "phiEle1/D");
      outTree[i]->Branch("chargeEle1", &chargeEle1, "chargeEle11/D");
    }
    // Second Ele
    if(i == 5) {
      outTree[i]->Branch("pTEle2", &pTEle2, "pTEle2/D");
      outTree[i]->Branch("etaEle2", &etaEle2, "etaEle2/D");
      outTree[i]->Branch("phiEle2", &phiEle2, "phiEle2/D");
      outTree[i]->Branch("chargeEle2", &chargeEle2, "chargeEle21/D");
    }
  }

  // prepare vectors for efficiency
  double Npassed_In = 0;
  double Npassed_HLT = 0; 
  double Npassed_PV = 0;
  double Npassed_PFHem = 0;
  double Npassed_CaloHem = 0;

  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // Prescaled Jet Triggers
  std::vector<std::string> maskHLT_DiJetAve30; maskHLT_DiJetAve30.push_back("HLT_DiJetAve30_v");
  std::vector<std::string> maskHLT_DiJetAve60; maskHLT_DiJetAve60.push_back("HLT_DiJetAve60_v");
  std::vector<std::string> maskHLT_Jet30; maskHLT_Jet30.push_back("HLT_Jet30_v");
  std::vector<std::string> maskHLT_Jet60; maskHLT_Jet60.push_back("HLT_Jet60_v");

  // other hadronic triggers
  std::vector<std::string> maskHLT_HT300_MHT75; maskHLT_HT300_MHT75.push_back("HLT_HT300_MHT75_v");
  std::vector<std::string> maskHLT_Meff640; maskHLT_Meff640.push_back("HLT_Meff640_v");
  std::vector<std::string> maskHLT_HT250_AlphaT0p55; maskHLT_HT250_AlphaT0p55.push_back("HLT_HT250_AlphaT0p55_v");
  std::vector<std::string> maskHLT_HT300_AlphaT0p52; maskHLT_HT300_AlphaT0p52.push_back("HLT_HT300_AlphaT0p52_v");

  // Razor BJet Triggers
  std::vector<std::string> maskHLT_R014_MR150_CentralJet40_BTagIP; maskHLT_R014_MR150_CentralJet40_BTagIP.push_back("HLT_R014_MR150_CentralJet40_BTagIP_v");
  std::vector<std::string> maskHLT_R014_MR450_CentralJet40_BTagIP; maskHLT_R014_MR450_CentralJet40_BTagIP.push_back("HLT_R014_MR450_CentralJet40_BTagIP_v");
  std::vector<std::string> maskHLT_R020_MR350_CentralJet40_BTagIP; maskHLT_R020_MR350_CentralJet40_BTagIP.push_back("HLT_R020_MR350_CentralJet40_BTagIP_v");
  std::vector<std::string> maskHLT_R025_MR250_CentralJet40_BTagIP; maskHLT_R025_MR250_CentralJet40_BTagIP.push_back("HLT_R025_MR250_CentralJet40_BTagIP_v");

  // prescaled Razor Triggers
  std::vector<std::string> maskHLT_R014_MR150; maskHLT_R014_MR150.push_back("HLT_R014_MR150_v");
  std::vector<std::string> maskHLT_R020_MR150; maskHLT_R020_MR150.push_back("HLT_R020_MR150_v");
  std::vector<std::string> maskHLT_R025_MR150; maskHLT_R025_MR150.push_back("HLT_R025_MR150_v");

  // hadronic razor triggers
  std::vector<std::string> maskHLT_R020_MR500; maskHLT_R020_MR500.push_back("HLT_R020_MR500_v");
  std::vector<std::string> maskHLT_R020_MR550; maskHLT_R020_MR550.push_back("HLT_R020_MR550_v");
  std::vector<std::string> maskHLT_R025_MR400; maskHLT_R025_MR400.push_back("HLT_R025_MR400_v");
  std::vector<std::string> maskHLT_R025_MR450; maskHLT_R025_MR450.push_back("HLT_R025_MR450_v");
  std::vector<std::string> maskHLT_R033_MR300; maskHLT_R033_MR300.push_back("HLT_R033_MR300_v");
  std::vector<std::string> maskHLT_R033_MR350; maskHLT_R033_MR350.push_back("HLT_R033_MR350_v");
  std::vector<std::string> maskHLT_R038_MR200; maskHLT_R038_MR200.push_back("HLT_R038_MR200_v");
  std::vector<std::string> maskHLT_R038_MR250; maskHLT_R038_MR250.push_back("HLT_R038_MR250_v");

  // Muon Razor Triggers
  std::vector<std::string> maskHLT_Mu8_R005_MR200; maskHLT_Mu8_R005_MR200.push_back("HLT_Mu8_R005_MR200_v");
  std::vector<std::string> maskHLT_Mu8_R020_MR200; maskHLT_Mu8_R020_MR200.push_back("HLT_Mu8_R020_MR200_v");
  std::vector<std::string> maskHLT_Mu8_R025_MR200; maskHLT_Mu8_R025_MR200.push_back("HLT_Mu8_R025_MR200_v");
  
  // Ele Razor Triggers
  std::vector<std::string> maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200; 
  maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200_v");
  std::vector<std::string> maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200; 
  maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200_v");
  std::vector<std::string> maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200; 
  maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200.push_back("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200_v");

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      
      // Prescaled Jet Triggers
      setRequiredTriggers(maskHLT_DiJetAve30); reloadTriggerMask(true); HLT_DiJetAve30 = hasPassedHLT();
      setRequiredTriggers(maskHLT_DiJetAve60); reloadTriggerMask(true); HLT_DiJetAve60 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Jet30); reloadTriggerMask(true); HLT_Jet30 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Jet60); reloadTriggerMask(true); HLT_Jet60 = hasPassedHLT();
      
      // other hadronic triggers
      setRequiredTriggers(maskHLT_HT300_MHT75); reloadTriggerMask(true); HLT_HT300_MHT75 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Meff640); reloadTriggerMask(true); HLT_Meff640 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT250_AlphaT0p55); reloadTriggerMask(true); HLT_HT250_AlphaT0p55 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT300_AlphaT0p52); reloadTriggerMask(true); HLT_HT300_AlphaT0p52 = hasPassedHLT();

      // Razor BJet Triggers
      setRequiredTriggers(maskHLT_R014_MR150_CentralJet40_BTagIP); reloadTriggerMask(true); HLT_R014_MR150_CentralJet40_BTagIP = hasPassedHLT();
      setRequiredTriggers(maskHLT_R014_MR450_CentralJet40_BTagIP); reloadTriggerMask(true); HLT_R014_MR450_CentralJet40_BTagIP = hasPassedHLT();
      setRequiredTriggers(maskHLT_R020_MR350_CentralJet40_BTagIP); reloadTriggerMask(true); HLT_R020_MR350_CentralJet40_BTagIP = hasPassedHLT();
      setRequiredTriggers(maskHLT_R025_MR250_CentralJet40_BTagIP); reloadTriggerMask(true); HLT_R025_MR250_CentralJet40_BTagIP = hasPassedHLT();
      
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
      
      // Ele Razor Triggers
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200); reloadTriggerMask(true); 
      HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200); reloadTriggerMask(true); 
      HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200); reloadTriggerMask(true); 
      HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200 = hasPassedHLT();
      
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
    
    Npassed_In += _weight;

    // HLT 
    //   vector<int> myHLT = getHLTOutput();
    //    HLTbit = 0;
    //    if(myHLT.size() == 1) HLTbit = myHLT[0];
    //    HLTbitEmulator = HTtrigger(150.);
    Npassed_HLT += _weight;

    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    Npassed_PV += _weight;

    int goodPFevent = true;
    int goodCaloevent = true;

    // HCAL FLAGS
    if(!eventPassHcalFilter()) continue;

    // Jet selection 
    vector<TLorentzVector> PFPUcorrJet;    
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);   
      if(myJet.Pt()>40. && fabs(myJet.Eta())< 2.4) 
	//	if(myJet.DeltaR(PFPhoton[iPh1]) > 0.35 && myJet.DeltaR(PFPhoton[iPh2]) > 0.35)
	PFPUcorrJet.push_back(myJet);
    }
    
    vector<TLorentzVector> CaloJet;    
    for(int i = 0; i < nAK5Jet; i++){
      TLorentzVector jet;
      double px = pxAK5Jet[i];
      double py = pyAK5Jet[i];
      double pz = pzAK5Jet[i];
      double E = sqrt(px*px+py*py+pz*pz);
      jet.SetPxPyPzE(px,py,pz,E);
      
      if(jet.Pt() < 40. || fabs(jet.Eta()) >= 3.0) continue; //jet kinematic cuts
      CaloJet.push_back(jet);
      
      // if(nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98){ //fails jet id
      // 	goodCaloevent = false;
      // }
      // if(fabs(jet.Eta()) < 2.55 && emFracAK5Jet[i] <= 0.01){ //fails jet id
      // 	goodCaloevent = false;
      // }
      // if(fabs(jet.Eta()) >= 2.55){ //fails jet ID
      // 	if(jet.Pt() > 80. && emFracAK5Jet[i] >= 1.){
      // 	  goodCaloevent = false;	
      // 	}
      // 	if(emFracAK5Jet[i] <= -0.9){
      // 	  goodCaloevent = false;
      // 	}
      //      }
    }

    // use PFMET
    TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);

    // dummy values
    passedPF = 0;
    pTPFHem1 = -9999;
    etaPFHem1 = -9999;
    phiPFHem1 = -9999;
    pTPFHem2 = -9999;
    etaPFHem2 = -9999;
    phiPFHem2 = -9999;
    PFR = -99999.;
    PFRsq = -99999.;
    PFMR = -99999.;

    // hemispheres
    vector<TLorentzVector> tmpJet = CombineJets(PFPUcorrJet);
    if(tmpJet.size() >= 2) {
      Npassed_PFHem += _weight;
      
      TLorentzVector PFHem1 = tmpJet[0];
      TLorentzVector PFHem2 = tmpJet[1];
      
      // compute boost
      double num = PFHem1.P()-PFHem2.P();
      double den = PFHem1.Pz()-PFHem2.Pz();      
      
      double MT = CalcMTR(PFHem1, PFHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcGammaMRstar(PFHem1, PFHem2);
      if(variable >0) Rvariable = MT/variable;
      
      // fill the tree
      passedPF = 1;
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

    // dummy values
    passedCalo = 0;
    pTCaloHem1 = -9999;
    etaCaloHem1 = -9999;
    phiCaloHem1 = -9999;
    pTCaloHem2 = -9999;
    etaCaloHem2 = -9999;
    phiCaloHem2 = -9999;
    CaloR = -99999.;
    CaloRsq = -99999.;
    CaloMR = -99999.;

    // hemispheres
    tmpJet = CombineJets(CaloJet);
    if(tmpJet.size() >= 2) {
      Npassed_CaloHem += _weight;
    
      TLorentzVector CaloHem1 = tmpJet[0];
      TLorentzVector CaloHem2 = tmpJet[1];
      
      // compute boost
      double num = CaloHem1.P()-CaloHem2.P();
      double den = CaloHem1.Pz()-CaloHem2.Pz();      
      double beta = num/den;
      
      double MT = CalcMTR(CaloHem1, CaloHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcGammaMRstar(CaloHem1, CaloHem2);
      if(variable >0) Rvariable = MT/variable;
      
      // fill the tree
      passedCalo = 1;
      pTCaloHem1 = CaloHem1.Pt();
      etaCaloHem1 = CaloHem1.Eta();
      phiCaloHem1 = CaloHem1.Phi();
      pTCaloHem2 = CaloHem2.Pt();
      etaCaloHem2 = CaloHem2.Eta();
      phiCaloHem2 = CaloHem2.Phi();
      CaloR = Rvariable;
      CaloRsq = Rvariable*Rvariable;
      CaloMR = variable;    
    }

    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;
    W = _weight;

    int nmuLoose = 0;
    int nmuTight = 0;
    TLorentzVector myMuLoose(0.1, 0., 0., 0.1);
    TLorentzVector myMuTight(0.1, 0., 0., 0.1);
    int iMuLoose = -99;
    int iMuTight = -99;

    for(int i=0; i<nMuon; i++) {
      if(muonPassTight(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 10.) {
	  nmuTight++;
	  if(thisMu.Pt() > myMuTight.Pt()) {
	    myMuTight = thisMu;
	    iMuTight = i;
	  }
	}
      }
    }

    for(int i=0; i<nMuon; i++) {
      if(muonPassLoose(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 10. && iMuTight != i) {
	  // we don't count the best Tight mu among the Loose muons
	  nmuLoose++;
	  if(thisMu.Pt() > myMuLoose.Pt()) {
	    myMuLoose = thisMu;
	    iMuLoose = i;
	  }
	}
      }
    }

    int neleWP95 = 0;
    int neleWP80 = 0;
    TLorentzVector myEleWP95(0.1, 0., 0., 0.1);
    TLorentzVector myEleWP80(0.1, 0., 0., 0.1);
    int iEleWP95 = -99;
    int iEleWP80 = -99;

    for(int i=0; i<nEle; i++) {
      if(electronPassWP80(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 12.) {
	  neleWP80++;
	  if(thisEle.Pt() > myEleWP80.Pt()) {
	    myEleWP80 = thisEle;
	    iEleWP80 = i;
	  }
	}
      }
    }

    for(int i=0; i<nEle; i++) {
      if(electronPassWP95(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 12. && i != iEleWP80) {
	  // we don't count the best WP80 Ele among the WP95 muons
	  neleWP95++;
	  if(thisEle.Pt() > myEleWP95.Pt()) {
	    myEleWP95 = thisEle;
	    iEleWP95 = i;
	  }
	}
      }
    }

    // Fill the tree per box
    pTMu1 = -999.;
    etaMu1 = -999.;
    phiMu1 = -999.;
    chargeMu1 = -999.;
    pTMu2 = -999.;
    etaMu2 = -999.;
    phiMu2 = -999.;
    chargeMu2 = -999.;    
    pTEle1 = -999.;
    etaEle1 = -999.;
    phiEle1 = -999.;
    chargeEle1 = -999.;
    pTEle2 = -999.;
    etaEle2 = -999.;
    phiEle2 = -999.;
    chargeEle2 = -999.;

    if(nmuTight>=1) { 
      // tight Mu leg
      pTMu1 = myMuTight.Pt();
      etaMu1 = myMuTight.Eta();
      phiMu1 = myMuTight.Phi();
      chargeMu1 = chargeMuon[iMuTight];
      if(nmuLoose>=1) { // DiMuon Box
	// loose Mu leg
	pTMu2 = myMuLoose.Pt();
	etaMu2 = myMuLoose.Eta();
	phiMu2 = myMuLoose.Phi();
	chargeMu2 = chargeMuon[iMuLoose];	  
	outTree[3]->Fill();
      } else if(neleWP80>=1) { // MuEle Box
	// tight Ele leg
	pTEle1 = myEleWP80.Pt();
	etaEle1 = myEleWP80.Eta();
	phiEle1 = myEleWP80.Phi();
	chargeEle1 = chargeEle[iEleWP80];
	outTree[4]->Fill();
      } else { // Mu Box
	outTree[2]->Fill();
      }
    } else if(neleWP80>=1) {
      // tight Ele leg
      pTEle1 = myEleWP80.Pt();
      etaEle1 = myEleWP80.Eta();
      phiEle1 = myEleWP80.Phi();
      chargeEle1 = chargeEle[iEleWP80];
      if(neleWP95>=1) { // DiEleBox
	// loose Ele leg
	pTEle2 = myEleWP95.Pt();
	etaEle2 = myEleWP95.Eta();
	phiEle2 = myEleWP95.Phi();
	chargeEle2 = chargeEle[iEleWP95];
	outTree[5]->Fill();
      } else { // Ele Box
	outTree[1]->Fill();
      }
    } else { // Had box
      outTree[0]->Fill();
    }

  }

  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
    
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_HLT",      &Npassed_HLT,      "Npassed_HLT/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("Npassed_PFHem",      &Npassed_PFHem,      "Npassed_PFHem/D");
  effTree->Branch("Npassed_CaloHem",      &Npassed_CaloHem,      "Npassed_CaloHem/D");
  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  effTree->Write();
  for(int i=0; i<6; i++) 
    outTree[i]->Write();
  file->Close();
}
  

