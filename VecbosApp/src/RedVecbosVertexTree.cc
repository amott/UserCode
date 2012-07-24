#include "../include/RedVecbosVertexTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedVecbosVertexTree::RedVecbosVertexTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","vecbos analysis tree for Z studies");

  myTree->Branch("nJets",         &myNJets,         "nJets/I");  
  myTree->Branch("trackDxyPV",    &myTrackDxyPV,    "trackDxyPV/F");  
  myTree->Branch("trackDxyOr",    &myTrackDxyOr,    "trackDxyOr/F");  
  myTree->Branch("trackDxyError", &myTrackDxyError, "trackDxyError/F");  
  myTree->Branch("trackDszPV",    &myTrackDszPV,    "trackDszPV/F");  
  myTree->Branch("trackDszOr",    &myTrackDszOr,    "trackDszOr/F");  
  myTree->Branch("trackDszError", &myTrackDszError, "trackDszError/F");  
  myTree->Branch("trackVertexX",  &myTrackVertexX,  "trackVertexX/F");  
  myTree->Branch("trackVertexY",  &myTrackVertexY,  "trackVertexY/F");  
  myTree->Branch("trackVertexZ",  &myTrackVertexZ,  "trackVertexZ/F");  
  myTree->Branch("PVx",           &myPVx,           "PVx/F");  
  myTree->Branch("PVy",           &myPVy,           "PVy/F");  
  myTree->Branch("PVz",           &myPVz,           "PVz/F");  
  myTree->Branch("trackDzError",  &myTrackDzError,  "trackDzError/F");  
  myTree->Branch("trackerIso",    &myTrackerIso,    "trackerIso/F");  
  myTree->Branch("ecalIso",       &myEcalIso,       "ecalIso/F");  
  myTree->Branch("hcalIso",       &myHcalIso,       "hcalIso/F");  
  myTree->Branch("wTransvMass",   &myWTransvMass,   "wTransvMass/F");  
  myTree->Branch("met",           &myMet,           "met/F");  
  myTree->Branch("chargeEle",     &myChargeEle,     "chargeEle/I");  
}

RedVecbosVertexTree::~RedVecbosVertexTree() {delete myFile;}

void RedVecbosVertexTree::addCSA07Infos() {

  myTree->Branch("CSA07weight",     &myWeight,     "CSA07weight/D");
  myTree->Branch("CSA07processId",  &myProcesId,   "CSA07processId/D");
  myTree->Branch("CSA07lumi",       &myLumi,       "CSA07lumi/F");
}

void RedVecbosVertexTree::addMcTruthInfos() {
  myTree->Branch("WToENuDecay",    &myWToENuDecay, "WToENuDecay/I");
}

void RedVecbosVertexTree::store() { myTree->Fill(); }

void RedVecbosVertexTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedVecbosVertexTree::fillAll(int nj, float dxyPv, float dxyOr, float dxyErr, float dszPv, float dszOr, float dszErr, float vtxX, float vtxY, float vtxZ, float pvx, float pvy, float pvz, float dzErr, float trackerIso, float ecalIso, float hcalIso, float wtmass, float met, int chargeEle){ 
  
  myNJets         = nj;
  myTrackDxyPV    = dxyPv;
  myTrackDxyOr    = dxyOr;
  myTrackDxyError = dxyErr;
  myTrackDszPV    = dszPv;
  myTrackDszOr    = dszOr;
  myTrackDszError = dszErr;
  myTrackVertexX  = vtxX;
  myTrackVertexY  = vtxY;
  myTrackVertexZ  = vtxZ;
  myPVx           = pvx;
  myPVy           = pvy;
  myPVz           = pvz;
  myTrackDzError  = dzErr;
  myTrackerIso    = trackerIso;
  myEcalIso       = ecalIso;
  myHcalIso       = hcalIso;
  myWTransvMass   = wtmass;
  myMet           = met;
  myChargeEle     = chargeEle;
}

void RedVecbosVertexTree::fillCSA07(double weight, double processId, float lumi) {
  myWeight = weight;
  myProcesId = processId;
  myLumi = lumi;
}

void RedVecbosVertexTree::fillMcTruth(int wtoenudecay) {
  myWToENuDecay = wtoenudecay;
}
