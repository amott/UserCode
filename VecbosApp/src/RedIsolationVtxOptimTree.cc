#include "include/RedIsolationVtxOptimTree.hh"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedIsolationVtxOptimTree::RedIsolationVtxOptimTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("promptDecay",  &myPromptDecay,  "promptDecay/I");
  myTree->Branch("trackerIsol",  &myTrackerIsol,  "trackerIsol/F");  
  myTree->Branch("hcalIsol",     &myHcalIsol,     "hcalIsol/F");  
  myTree->Branch("ecalJurIsol",  &myEcalJurIsol,  "ecalJurIsol/F");  
  myTree->Branch("ecalGTIsol",   &myEcalGTIsol,   "ecalGTIsol/F");  
  myTree->Branch("combinedIsol", &myCombinedIsol, "combinedIsol/F");  
  myTree->Branch("dzVtx",        &myDzVtx,        "dzVtx/F");
  myTree->Branch("dxyVtx",       &myDxyVtx,       "dxyVtx/F");
  myTree->Branch("dxyErrVtx",    &myDxyErrVtx,    "dxyErrVtx/F");

}

RedIsolationVtxOptimTree::~RedIsolationVtxOptimTree() 
{
  delete myFile;
}

void RedIsolationVtxOptimTree::addCSA07Infos() {

  myTree->Branch("CSA07weight",     &myWeight,     "CSA07weight/D");
  myTree->Branch("CSA07processId",  &myProcesId,   "CSA07processId/D");
  myTree->Branch("CSA07lumi",       &myLumi,       "CSA07lumi/F");

}

void RedIsolationVtxOptimTree::addKinematicsInfos() {
  myTree->Branch("pt",  &myPt, "pt/F");
  myTree->Branch("eta", &myEta, "eta/F");
}

void RedIsolationVtxOptimTree::store()
{
  myTree->Fill();
}


void RedIsolationVtxOptimTree::save() 
{
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedIsolationVtxOptimTree::fillAll(int promptDecay, float ftracker, float fhcal, float fecalJ, float fecalGT, float fcombined, float dzvtx, float dxyvtx, float dxyerrvtx)
{
  myPromptDecay = promptDecay;
  myTrackerIsol = ftracker;
  myHcalIsol    = fhcal;
  myEcalJurIsol = fecalJ;
  myEcalGTIsol  = fecalGT;
  myDzVtx       = dzvtx;
  myDxyVtx      = dxyvtx;
  myDxyErrVtx   = dxyerrvtx;
  myCombinedIsol = fcombined;
}

void RedIsolationVtxOptimTree::fillKinematics(float pt, float eta) {
  myPt = pt;
  myEta = eta;
}

void RedIsolationVtxOptimTree::fillCSA07(double weight, double processId, float lumi) {
  myWeight = weight;
  myProcesId = processId;
  myLumi = lumi;
}
