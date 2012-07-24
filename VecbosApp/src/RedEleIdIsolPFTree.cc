#include "include/RedEleIdIsolPFTree.hh"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedEleIdIsolPFTree::RedEleIdIsolPFTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","eleID tree");

  myTree->Branch("eta",          &myEta,          "eta/F");
  myTree->Branch("pt",           &myPt,           "pt/F");
  myTree->Branch("mva",          &myMva,          "mva/F");
  myTree->Branch("hoe",          &myHoE,          "hoe/F");
  myTree->Branch("deltaEta",     &myDeltaEta,     "deltaEta/F");
  myTree->Branch("deltaPhi",     &myDeltaPhi,     "deltaPhi/F");
  myTree->Branch("eop",          &myEoP,          "eop/F");
  myTree->Branch("eopout",       &myEoPout,       "eopout/F");
  myTree->Branch("sigmaIeIe",    &mySigmaIeIe,    "sigmaIeIe/F");
  myTree->Branch("charged",      &myCharged,      "charged/F");
  myTree->Branch("neutral",      &myNeutral,      "neutral/F");
  myTree->Branch("photon",       &myPhoton,       "photon/F");
  myTree->Branch("combined",     &myCombined,     "combined/F");
  myTree->Branch("innerLayers",  &myInnerLayers,  "innerLayers/I");
  myTree->Branch("tranImpPar",   &myTranImpPar,   "tranImpPar/F");

  myTree->Branch("trackerIsol",  &myTrackerIsol,  "trackerIsol/F");
  myTree->Branch("ecalIsol",     &myEcalIsol,     "ecalIsol/F");
  myTree->Branch("hcalIsol",     &myHcalIsol,     "hcalIsol/F");
  
  myTree->Branch("passedEleId",  &myPassedEleId,  "passedEleId/I");
  myTree->Branch("passedIsol",   &myPassedIsol,   "passedIsol/I");
  myTree->Branch("passedConv",   &myPassedConv,   "passedConv/I");
  
  myTree->Branch("eleChargedIso03nV", &myEleChargedIso03nV, "eleChargedIso03nV/F");
  myTree->Branch("eleChargedIso04nV", &myEleChargedIso04nV, "eleChargedIso04nV/F");
  myTree->Branch("eleChargedIso05nV", &myEleChargedIso05nV, "eleChargedIso05nV/F");
  myTree->Branch("eleChargedIso03v",  &myEleChargedIso03v,  "eleChargedIso03v/F");
  myTree->Branch("eleChargedIso04v",  &myEleChargedIso04v,  "eleChargedIso04v/F");
  myTree->Branch("eleChargedIso05v",  &myEleChargedIso05v,  "eleChargedIso05v/F");

  myTree->Branch("eleChargedIsoNvc03nV", &myEleChargedIsoNvc03nV, "eleChargedIsoNvc03nV/F");
  myTree->Branch("eleChargedIsoNvc04nV", &myEleChargedIsoNvc04nV, "eleChargedIsoNvc04nV/F");
  myTree->Branch("eleChargedIsoNvc05nV", &myEleChargedIsoNvc05nV, "eleChargedIsoNvc05nV/F");
  myTree->Branch("eleChargedIsoNvc03v",  &myEleChargedIsoNvc03v,  "eleChargedIsoNvc03v/F");
  myTree->Branch("eleChargedIsoNvc04v",  &myEleChargedIsoNvc04v,  "eleChargedIsoNvc04v/F");
  myTree->Branch("eleChargedIsoNvc05v",  &myEleChargedIsoNvc05v,  "eleChargedIsoNvc05v/F");

  myTree->Branch("eleNeutralIso03nV", &myEleNeutralIso03nV, "eleNeutralIso03nV/F");
  myTree->Branch("eleNeutralIso04nV", &myEleNeutralIso04nV, "eleNeutralIso04nV/F");
  myTree->Branch("eleNeutralIso05nV", &myEleNeutralIso05nV, "eleNeutralIso05nV/F");
  myTree->Branch("eleNeutralIso03v",  &myEleNeutralIso03v,  "eleNeutralIso03v/F");
  myTree->Branch("eleNeutralIso04v",  &myEleNeutralIso04v,  "eleNeutralIso04v/F");
  myTree->Branch("eleNeutralIso05v",  &myEleNeutralIso05v,  "eleNeutralIso05v/F");

  myTree->Branch("elePhotonsIso03nV", &myElePhotonsIso03nV, "elePhotonsIso03nV/F");
  myTree->Branch("elePhotonsIso04nV", &myElePhotonsIso04nV, "elePhotonsIso04nV/F");
  myTree->Branch("elePhotonsIso05nV", &myElePhotonsIso05nV, "elePhotonsIso05nV/F");
  myTree->Branch("elePhotonsIso03v",  &myElePhotonsIso03v,  "elePhotonsIso03v/F");
  myTree->Branch("elePhotonsIso04v",  &myElePhotonsIso04v,  "elePhotonsIso04v/F");
  myTree->Branch("elePhotonsIso05v",  &myElePhotonsIso05v,  "elePhotonsIso05v/F");
}

RedEleIdIsolPFTree::~RedEleIdIsolPFTree() { 

  delete myFile;
}

void RedEleIdIsolPFTree::store() { 

  myTree->Fill();
}


void RedEleIdIsolPFTree::save() { 

  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedEleIdIsolPFTree::fillAll(float pt, float eta, float mva, float hoe, float deta, float dphi, float eop, float eopo, float see, float ch, float nh, float ph, float comb, int il, float tip, float tisol, float eisol, float hisol) {
  
  myPt          = pt;
  myEta         = eta;
  myMva         = mva;
  myHoE         = hoe;
  myDeltaEta    = deta;
  myDeltaPhi    = dphi;
  myEoP         = eop;
  myEoPout      = eopo;
  mySigmaIeIe   = see;
  myCharged     = ch;
  myNeutral     = nh;
  myPhoton      = ph;
  myCombined    = comb;
  myInnerLayers = il;
  myTranImpPar  = tip;
  myTrackerIsol = tisol;
  myEcalIsol    = eisol;
  myHcalIsol    = hisol;
}

void RedEleIdIsolPFTree::fillSetpsAfter(int eleId, int isol, int conv) {

  myPassedEleId = eleId;
  myPassedIsol  = isol;
  myPassedConv  = conv;
}

void RedEleIdIsolPFTree::fillIsol(
float cnv3, float cnv4, float cnv5, float cv3, float cv4, float cv5, float nvcnv3, float nvcnv4, float nvcnv5, float nvcv3, float nvcv4, float nvcv5, float nnv3, float nnv4, float nnv5, float nv3, float nv4, float nv5, float pnv3, float pnv4, float pnv5, float pv3, float pv4, float pv5){

  myEleChargedIso03nV = cnv3;
  myEleChargedIso04nV = cnv4;
  myEleChargedIso05nV = cnv5;
  myEleChargedIso03v  = cv3;
  myEleChargedIso04v  = cv4;
  myEleChargedIso05v  = cv5;

  myEleChargedIsoNvc03nV = nvcnv3;
  myEleChargedIsoNvc04nV = nvcnv4;
  myEleChargedIsoNvc05nV = nvcnv5;
  myEleChargedIsoNvc03v  = nvcv3;
  myEleChargedIsoNvc04v  = nvcv4;
  myEleChargedIsoNvc05v  = nvcv5;

  myEleNeutralIso03nV = nnv3;
  myEleNeutralIso04nV = nnv4;
  myEleNeutralIso05nV = nnv5;
  myEleNeutralIso03v  = nv3;
  myEleNeutralIso04v  = nv4;
  myEleNeutralIso05v  = nv5;

  myElePhotonsIso03nV = pnv3;
  myElePhotonsIso04nV = pnv4;
  myElePhotonsIso05nV = pnv5;
  myElePhotonsIso03v  = pv3;
  myElePhotonsIso04v  = pv4;
  myElePhotonsIso05v  = pv5;
}

