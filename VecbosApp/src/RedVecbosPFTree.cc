#include "../include/RedVecbosPFTree.h"
#include <assert.h>

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

// Root
#include "TFile.h"
#include "TTree.h"

using namespace std;

RedVecbosPFTree::RedVecbosPFTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","vecbos analysis tree for Z studies");
  
  myTree->Branch("run",                 &myRun,                 "run/I");
  myTree->Branch("lumi",                &myLumi,                "lumi/I");
  myTree->Branch("event",               &myEvent,               "event/I");
  myTree->Branch("nJets",               &myNJets,               "nJets/I");  
  myTree->Branch("nPFJets",             &myNPFJets,             "nPFJets/I");  
  myTree->Branch("nGenJets",            &myNGenJets,            "nGenJets/I");
  myTree->Branch("calomt",              &myCaloTransvMass,      "calomt/F");  
  myTree->Branch("tcmt",                &myTCTransvMass,        "tcmt/F");  
  myTree->Branch("pfmt",                &myPFTransvMass,        "pfmt/F");  
  myTree->Branch("calomet",             &myCaloMet,             "calomet/F");  
  myTree->Branch("tcmet",               &myTcMet,               "tcmet/F");  
  myTree->Branch("pfmet",               &myPFMet,               "pfmet/F");  
  myTree->Branch("caloJetSelected",     &myCaloJetSelected,     "caloJetSelected/I");    
  myTree->Branch("PFJetSelected",       &myPFJetSelected,       "PFJetSelected/I");    
  // myTree->Branch("genWpt",              &myGenWpt,              "genWpt/F");
}

RedVecbosPFTree::~RedVecbosPFTree() {delete myFile;}

void RedVecbosPFTree::addMcTruthInfos() {
  myTree->Branch("promptDecay",    &myPromptDecay, "promptDecay/I");
}

void RedVecbosPFTree::addElectronPFInfos() {

  myTree->Branch("PFPt",          myPFPt,          "myPFPt[2]/F");
  myTree->Branch("PFEta",         myPFEta,         "myPFEta[2]/F");
  myTree->Branch("PFMva",         myPFMva,         "myPFMva[2]/F");
  myTree->Branch("PFChargedIso",  myPFChargedIso,  "myPFChargedIso[2]/F");
  myTree->Branch("PFNeutralIso",  myPFNeutralIso,  "myPFNeutralIso[2]/F");
  myTree->Branch("PFPhotonsIso",  myPFPhotonsIso,  "myPFPhotonsIso[2]/F");
  myTree->Branch("PFCombinedIso", myPFCombinedIso, "myPFCombinedIso[2]/F");
  myTree->Branch("PFCharge",      myPFCharge,      "myPFCharge[2]/I");
}

void RedVecbosPFTree::store() { myTree->Fill(); }

void RedVecbosPFTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedVecbosPFTree::fillJetMultiplicities(int nJets, int nPFJets, int nGenJets, int caloJetSelected, int PFJetSelected) {
  myNJets              = nJets;
  myNPFJets            = nPFJets;
  myNGenJets           = nGenJets;
  myCaloJetSelected    = caloJetSelected;
  myPFJetSelected      = PFJetSelected;
}

void RedVecbosPFTree::fillElectronsPF(float pt[2], float eta[2], float mva[2], float ci[2], float ni[2], float pi[2], float combIso[2], int charge[2]) {
  
  for(int i=0; i<2; i++) {
    myPFPt[i]          = pt[i];
    myPFEta[i]         = eta[i];
    myPFMva[i]         = mva[i];
    myPFChargedIso[i]  = ci[i];
    myPFNeutralIso[i]  = ni[i];
    myPFPhotonsIso[i]  = pi[i];
    myPFCombinedIso[i] = combIso[i];
    myPFCharge[i]      = charge[i];
  }
}

void RedVecbosPFTree::fillKinematicsPF(float calotmass, float tctmass, float pftmass, float calomet, float tcmet, float pfmet) {

  myCaloTransvMass = calotmass;
  myTCTransvMass   = tctmass;
  myPFTransvMass   = pftmass;
  myCaloMet        = calomet;
  myTcMet          = tcmet;
  myPFMet          = pfmet;
}

void RedVecbosPFTree::fillMcTruth(int promptdecay) {
  myPromptDecay = promptdecay;
}

void RedVecbosPFTree::fillRunInfo(int run, int lumi, int event) {
  myRun   = run;
  myLumi  = lumi;
  myEvent = event;
}

// void RedVecbosPFTree::fillArcs(float gpt) {
//  myGenWpt = gpt;
//}
