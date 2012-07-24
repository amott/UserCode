#include "include/RedEleIDOptimTree.hh"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedEleIDOptimTree::RedEleIDOptimTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("class",    &myClassification, "class/I");
  myTree->Branch("recoFlag", &myRecoFlag,  "recoFlag/I");
  myTree->Branch("eta",      &myEta,       "eta/F");
  myTree->Branch("pt",       &myPt,        "pt/F");
  myTree->Branch("deltaEta", &myDeltaEta,  "deltaEta/F");  
  myTree->Branch("deltaPhi", &myDeltaPhi,  "deltaPhi/F");  
  myTree->Branch("hoe",      &myHoe,       "hoe/F");  
  myTree->Branch("s9s25",    &myS9s25,     "s9s25/F");  
  myTree->Branch("see",      &mySee,       "see/F");  
  myTree->Branch("eopOut",   &myEopOut,    "eopOut/F");  
  myTree->Branch("fbrem",    &myFBrem,     "fbrem/F");  

}

RedEleIDOptimTree::~RedEleIDOptimTree()
{ 
  delete myFile;
}

void RedEleIDOptimTree::store()
{ 
  myTree->Fill();
}


void RedEleIDOptimTree::save()
{ 
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedEleIDOptimTree::fillAll(int classification, int recoflag, float eta, float pt,
                                float deta, float dphi, float hoe, float s9s25, float see, float eopout, float fbrem)
{ 
  myClassification = classification;
  myRecoFlag = recoflag;
  myEta = eta;
  myPt = pt;
  myDeltaEta = deta;
  myDeltaPhi = dphi;
  myHoe      = hoe;
  myS9s25    = s9s25;
  mySee      = see;
  myEopOut   = eopout;
  myFBrem    = fbrem;
}
