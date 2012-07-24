#include "../include/RedTopTree.h"
#include <assert.h>

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <assert.h>

// Root
#include "TFile.h"
#include "TTree.h"

using namespace std;

RedTopTree::RedTopTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","top studies");

  myTree->Branch("run",   &myRun,   "run/I");
  myTree->Branch("lumi",  &myLumi,  "lumi/I");
  myTree->Branch("event", &myEvent, "event/I");

  myTree->Branch("channelEE", &myChannelEE, "channelEE/I");  
  myTree->Branch("channelMM", &myChannelMM, "channelMM/I");  
  myTree->Branch("channelEM", &myChannelEM, "channelEM/I");  
  myTree->Branch("pt1",       &myPt1,       "pt1/F");
  myTree->Branch("pt2",       &myPt2,       "pt2/F");
  myTree->Branch("eta1",      &myEta1,      "eta1/F");
  myTree->Branch("eta2",      &myEta2,      "eta2/F");
  myTree->Branch("charge1",   &myCharge1,   "charge1/I");
  myTree->Branch("charge2",   &myCharge2,   "charge2/I");
  myTree->Branch("invMass",   &myInvMass,   "invMass/F");
  myTree->Branch("nPFJets",   &myNPFJets,   "nPFJets/I");
  myTree->Branch("met",       &myMet,       "met/F");  

  myTree->Branch("dzEle",     &myDzEle,     "dzEle/F");
  myTree->Branch("dzMu",      &myDzMu,      "dzMu/F");

  myTree->Branch("nBTagJets", myNJetsBTagged, "nBTagJets[5]/I");

  myTree->Branch("etJet1",      &myEtJet1,     "etJet1/F");
  myTree->Branch("etJet2",      &myEtJet2,     "etJet2/F");
  myTree->Branch("foundB",      &myFoundB,     "foundB/I");
  myTree->Branch("foundBhad",   &myFoundBhad,  "foundBhad/I");
  myTree->Branch("matchedAll",  &myMatchedAll, "matchedAll/I");
  myTree->Branch("matchedJet1", &myMatchedJet1,"matchedJet1/I");
  myTree->Branch("matchedJet2", &myMatchedJet2,"matchedJet2/I");
  myTree->Branch("taggedJet1",  &myTaggedJet1, "taggedJet1/I");
  myTree->Branch("taggedJet2",  &myTaggedJet2, "taggedJet2/I");
}

RedTopTree::~RedTopTree() {delete myFile;}

void RedTopTree::store() { myTree->Fill(); }

void RedTopTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedTopTree::fillGeneral(int cEE, int cMM, int cEM, float p1, float p2, float e1, float e2, int c1, int c2, float im, int npfj, float met, float ml1, float ml2) {
  myChannelEE = cEE;
  myChannelMM = cMM;
  myChannelEM = cEM;
  myPt1 = p1;
  myPt2 = p2;
  myEta1 = e1; 
  myEta2 = e2; 
  myCharge1 = c1; 
  myCharge2 = c2;
  myInvMass = im; 
  myNPFJets = npfj; 
  myMet = met; 
  myMcLeptDr1 = ml1;
  myMcLeptDr2 = ml2;
}

void RedTopTree::fillVertexComp(float de, float dm) {
  
  myDzEle = de;
  myDzMu  = dm;
}

void RedTopTree::fillBTagJetMultiplicities(int nJetsBTagged[5]) {
  for(int i=0; i<5; i++) myNJetsBTagged[i] = nJetsBTagged[i];
}

void RedTopTree::fillBTag(float etj1, float etj2, int fb, int fbh, int ma,int mj1, int mj2, int tj1, int tj2) {

  myEtJet1      = etj1;
  myEtJet2      = etj2;
  myFoundB      = fb;
  myFoundBhad   = fbh;
  myMatchedAll  = ma;
  myMatchedJet1 = mj1;
  myMatchedJet2 = mj2;
  myTaggedJet1  = tj1;
  myTaggedJet2  = tj2;
}

void RedTopTree::fillRunInfo(int run, int lumi, int event) {
  myRun   = run;
  myLumi  = lumi;
  myEvent = event;
}
