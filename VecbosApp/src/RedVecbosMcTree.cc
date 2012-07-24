#include "../include/RedVecbosMcTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedVecbosMcTree::RedVecbosMcTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","vecbos analysis tree for Mc studies");
  myTree->Branch("processId",  &myProcessId, "processId/D");
  myTree->Branch("lept1Id",    &myLept1Id,   "lept1Id/F");  
  myTree->Branch("lept1Eta",   &myLept1Eta,  "lept1Eta/F");  
  myTree->Branch("lept1Phi",   &myLept1Phi,  "lept1Phi/F");  
  myTree->Branch("lept1Ene",   &myLept1Ene,  "lept1Ene/F");  
  myTree->Branch("lept2Id",    &myLept2Id,   "lept2Id/F");  
  myTree->Branch("lept2Eta",   &myLept2Eta,  "lept2Eta/F");  
  myTree->Branch("lept2Phi",   &myLept2Phi,  "lept2Phi/F");  
  myTree->Branch("lept2Ene",   &myLept2Ene,  "lept2Ene/F");  
}

RedVecbosMcTree::~RedVecbosMcTree() {delete myFile;}

void RedVecbosMcTree::store() { myTree->Fill(); }


void RedVecbosMcTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedVecbosMcTree::fillAll(int procId, float id1, float eta1, float phi1, float ene1, float id2, float eta2, float phi2, float ene2) { 
  myProcessId = procId;
  myLept1Id   = id1;
  myLept1Eta  = eta1;
  myLept1Phi  = phi1;
  myLept1Ene  = ene1;
  myLept2Id   = id2;
  myLept2Eta  = eta2;
  myLept2Phi  = phi2;
  myLept2Ene  = ene2;
}

