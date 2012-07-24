#include "../include/RedVecbosTree.h"
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

RedVecbosTree::RedVecbosTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","vecbos analysis tree for Z studies");
  
  myTree->Branch("run", &myRun, "run/I");
  myTree->Branch("lumi", &myLumi, "lumi/I");
  myTree->Branch("event", &myEvent, "event/I");
  myTree->Branch("nJetsHi",             &myNJetsHi,             "nJetsHi/I");  
  myTree->Branch("nPFJetsHi",           &myNPFJetsHi,           "nPFJetsHi/I");  
  myTree->Branch("nTrackJetsHi",        &myNTrackJetsHi,        "nTrackJetsHi/I");
  myTree->Branch("nGenJetsHi",          &myNGenJetsHi,          "nGenJetsHi/I");
  myTree->Branch("nJetsLo",             &myNJetsLo,             "nJetsLo/I");  
  myTree->Branch("nPFJetsLo",           &myNPFJetsLo,           "nPFJetsLo/I");  
  myTree->Branch("nTrackJetsLo",        &myNTrackJetsLo,        "nTrackJetsLo/I");
  myTree->Branch("nGenJetsLo",          &myNGenJetsLo,          "nGenJetsLo/I");
  myTree->Branch("nBTagJets",           myNJetsBTagged,         "nBTagJets[5]/I");
  myTree->Branch("mee",                 &myInvMass,             "mee/F");
  myTree->Branch("mt",                  &myTransvMass,          "mt/F");  
  myTree->Branch("tcmt",                &myTCTransvMass,        "tcmt/F");  
  myTree->Branch("pfmt",                &myPFTransvMass,        "pfmt/F");  
  myTree->Branch("met",                 &myMet,                 "met/F");  
  myTree->Branch("tcmet",               &myTCMet,               "tcmet/F");  
  myTree->Branch("pfmet",               &myPFMet,               "pfmet/F");  
  myTree->Branch("ptW",                 &myPtW,                 "ptW/F");
  myTree->Branch("ptWRecoil",           &myPtWRecoil,           "ptWRecoil/F");
  myTree->Branch("etaW",                &myEtaW,                "etaW/F");
  myTree->Branch("ptZ",                 &myPtZ,                 "ptZ/F");
  myTree->Branch("etaZ",                &myEtaZ,                "etaZ/F");
  myTree->Branch("pfmtZMet",            &myMTZMet,              "pfmtZMet/F");
  myTree->Branch("mZJet",               &myMZJet,               "mZJet/F");
  myTree->Branch("templMetUncorr",      &myTemplMetUncorr,      "templMetUncorr/F");    
  myTree->Branch("templMetCorr",        &myTemplMetCorr,        "templMetCorr/F");    
  myTree->Branch("templTrMassUncorr",   &myTemplTrMassUncorr,   "templTrMassUncorr/F");    
  myTree->Branch("templTrMassCorr",     &myTemplTrMassCorr,     "templTrMassCorr/F");    
  myTree->Branch("caloJetSelectedHi",   &myCaloJetSelectedHi,   "caloJetSelectedHi/I");    
  myTree->Branch("PFJetSelectedHi",     &myPFJetSelectedHi,     "PFJetSelectedHi/I");    
  myTree->Branch("trackJetSelectedHi",  &myTrackJetSelectedHi,  "trackJetSelectedHi/I");    
  myTree->Branch("caloJetSelectedLo",   &myCaloJetSelectedLo,   "caloJetSelectedLo/I");    
  myTree->Branch("PFJetSelectedLo",     &myPFJetSelectedLo,     "PFJetSelectedLo/I");    
  myTree->Branch("trackJetSelectedLo",  &myTrackJetSelectedLo,  "trackJetSelectedLo/I");    
  myTree->Branch("leadingCaloJetPt",    &myLeadingCaloJetPt,    "leadingCaloJetPt/F");
  myTree->Branch("leadingPFJetPt",      &myLeadingPFJetPt,      "leadingPFJetPt/F");
  myTree->Branch("emFracCaloJet",       &myEmFracCaloJet,       "emFracCaloJet/F");
  myTree->Branch("emFracPFJet",         &myEmFracPFJet,         "emFracPFJet/F");
  myTree->Branch("diJetMass",           myDiJetMass,            "diJetMass[2]/F");
  myTree->Branch("diJetPt",             myDiJetPt,              "diJetPt[2]/F");
  myTree->Branch("deltaEtaJet",         myDeltaEtaJet,          "deltaEtaJet[2]/F");
  myTree->Branch("deltaPhiJetMet",      myDeltaPhiJetMet,       "deltaPhiJetMet[2]/F");
  myTree->Branch("ptGen1",              &myPtGen1,              "ptGen1/F");
  myTree->Branch("ptGen2",              &myPtGen2,              "ptGen2/F");
  myTree->Branch("etaGen1",             &myEtaGen1,             "etaGen1/F");
  myTree->Branch("etaGen2",             &myEtaGen2,             "etaGen2/F");
  myTree->Branch("meeGen",              &myMeeGen,              "meeGen/F");
  myTree->Branch("nvtx",                &myNumVtx,              "nvtx/I");
  myTree->Branch("rho",                 &myRho,                 "rho/F");
}

RedVecbosTree::~RedVecbosTree() {delete myFile;}

void RedVecbosTree::addMcTruthInfos() {
  myTree->Branch("promptDecay",    &myPromptDecay, "promptDecay/I");
  myTree->Branch("nPU", &myNumInteractions, "nPU/I");
}

void RedVecbosTree::addPFelectrons() {
  myTree->Branch("PFelePt",        &myPFelePt,       "PFelePt/F");
  myTree->Branch("PFeleEta",       &myPFeleEta,      "PFeleEta/F");
  myTree->Branch("PFelePhi",       &myPFelePhi,      "PFelePhi/F");
  myTree->Branch("PFeleDeltaR",    &myPFeleDeltaR,   "PFeleDeltaR/F");
  myTree->Branch("PFeleMva",       &myPFeleMva,      "PFeleMva/F");
  myTree->Branch("NPFele",         &myNPFele,        "myNPFele/I");
}

void RedVecbosTree::addBTagEVTInfos() {
  myTree->Branch("combinedSecondaryVertexBJetTags",    &myCombinedSecondaryVertexBJetTags,      "combinedSecondaryVertexBJetTags/F");
  myTree->Branch("combinedSecondaryVertexMVABJetTags", &myCombinedSecondaryVertexMVABJetTags,   "combinedSecondaryVertexMVABJetTags/F");
  myTree->Branch("jetBProbabilityBJetTags",            &myJetBProbabilityBJetTags,              "jetBProbabilityBJetTags/F");
  myTree->Branch("jetProbabilityBJetTags",             &myJetProbabilityBJetTags,               "jetProbabilityBJetTags/F");
  myTree->Branch("simpleSecondaryVertexBJetTags",      &mySimpleSecondaryVertexBJetTags,        "simpleSecondaryVertexBJetTags/F");
  myTree->Branch("softMuonBJetTags",                   &mySoftMuonBJetTags,                     "softMuonBJetTags/F");
  myTree->Branch("trackCountingHighPurBJetTags",       &myTrackCountingHighPurBJetTags,         "trackCountingHighPurBJetTags/F");
  myTree->Branch("trackCountingHighEffBJetTags",       &myTrackCountingHighEffBJetTags,         "trackCountingHighEffBJetTags/F");
  myTree->Branch("foundMcB",                           &myFoundMcB,                             "foundMcB/I");
  myTree->Branch("foundMcBmum",                        &myFoundMcBmum,                          "foundMcBmum/I");
  myTree->Branch("nB",                                 &myNumB,                                 "nB/I");
}

void RedVecbosTree::addEventShapeInfos() {
  myTree->Branch("MHTphiJet",           myMHTphiJet,            "MHTphiJet[6]/F");  
  myTree->Branch("MHTphiPFJet",         myMHTphiPFJet,          "MHTphiPFJet[6]/F");  
  myTree->Branch("MHTphiMET",           &myMHTphiMET,           "MHTphiMET/F");  
}

void RedVecbosTree::addElectronInfos() {
  
  myTree->Branch("recoflag", myRecoflag, "recoflag[2]/I");
  myTree->Branch("pt", myPt, "pt[2]/F");
  myTree->Branch("eta", myEta, "eta[2]/F");
  myTree->Branch("phi", myPhi, "phi[2]/F");
  myTree->Branch("classification", myClassification, "classification[2]/I");
  myTree->Branch("nbrem", myNBrems, "nbrem[2]/I");
  myTree->Branch("deta", myDeta, "deta[2]/F");
  myTree->Branch("dphi", myDphi, "dphi[2]/F");
  myTree->Branch("hoe", myHoe, "hoe[2]/F");
  myTree->Branch("see", mySee, "see[2]/F");
  myTree->Branch("spp", mySpp, "spp[2]/F");
  myTree->Branch("eop", myEop, "eop[2]/F");
  myTree->Branch("fbrem", myFbrem, "fbrem[2]/F");
  myTree->Branch("trackerIso", myTrackerIso, "trackerIso[2]/F");
  myTree->Branch("hcalIso", myHcalIso, "hcalIso[2]/F");
  myTree->Branch("ecalJIso", myEcalJIso, "ecalJIso[2]/F");
  myTree->Branch("ecalGTIso", myEcalGTIso, "ecalGTIso[2]/F");
  myTree->Branch("combinedIso", myCombinedIso, "combinedIso[2]/F");
  myTree->Branch("charge", myCharge, "charge[2]/I");
  myTree->Branch("missHits", myMissHits, "missHits[2]/I");
  myTree->Branch("dist", myDist, "dist[2]/F");
  myTree->Branch("dcot", myDcot, "dcot[2]/F");
  myTree->Branch("lh", myLh, "lh[2]/F");
  myTree->Branch("e9esc", myE9ESC, "e9esc[2]/F");
  myTree->Branch("e25esc", myE25ESC, "e25esc[2]/F");
  myTree->Branch("esc", myESC, "esc[2]/F");
  myTree->Branch("eseed", myEseed, "eseed[2]/F");
  myTree->Branch("pin", myPin, "pin[2]/F");
}

void RedVecbosTree::store() { myTree->Fill(); }


void RedVecbosTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedVecbosTree::fillJetMultiplicities(int nJetsHi, int nPFJetsHi, int nTrackJetsHi, int nGenJetsHi, 
                                          int nJetsLo, int nPFJetsLo, int nTrackJetsLo, int nGenJetsLo, 
                                          int caloJetSelectedHi, int PFJetSelectedHi, int trackJetSelectedHi, 
                                          int caloJetSelectedLo, int PFJetSelectedLo, int trackJetSelectedLo, 
                                          float leadingCaloPt, float leadingPFPt, float emFracCaloJet, float emFracPFJet) {
  myNJetsHi            = nJetsHi;
  myNPFJetsHi          = nPFJetsHi;
  myNTrackJetsHi       = nTrackJetsHi;
  myNGenJetsHi         = nGenJetsHi;
  myNJetsLo            = nJetsLo;
  myNPFJetsLo          = nPFJetsLo;
  myNTrackJetsLo       = nTrackJetsLo;
  myNGenJetsLo         = nGenJetsLo;
  myCaloJetSelectedHi  = caloJetSelectedHi;
  myPFJetSelectedHi    = PFJetSelectedHi;
  myTrackJetSelectedHi = trackJetSelectedHi;
  myCaloJetSelectedLo  = caloJetSelectedLo;
  myPFJetSelectedLo    = PFJetSelectedLo;
  myTrackJetSelectedLo = trackJetSelectedLo;
  myLeadingCaloJetPt   = leadingCaloPt;
  myLeadingPFJetPt     = leadingPFPt;
  myEmFracCaloJet      = emFracCaloJet;
  myEmFracPFJet        = emFracPFJet;
}

void RedVecbosTree::fillElectrons(int recoflag[2], float pt[2], float eta[2], float phi[2],
                                  int classification[2], int nbrems[2], float deta[2], float dphi[2], float hoe[2], float see[2], float spp[2], float eop[2], float fbrem[2],
                                  float trackerIso[2], float hcalIso[2], float ecalJIso[2], float ecalGTIso[2], float combinedIso[2], int charge[2],
                                  int missHits[2], float dist[2], float dcot[2], float lh[2], float e9esc[2], float e25esc[2]) {
  
  for(int i=0; i<2; i++) {
    myRecoflag[i] = recoflag[i];
    myPt[i] = pt[i];
    myEta[i] = eta[i];
    myPhi[i] = phi[i];
    myClassification[i] = classification[i];
    myNBrems[i] = nbrems[i];
    myDeta[i] = deta[i];
    myDphi[i] = dphi[i];
    myHoe[i] = hoe[i];
    mySee[i] = see[i];
    mySpp[i] = spp[i];
    myEop[i] = eop[i];
    myFbrem[i] = fbrem[i];
    myTrackerIso[i] = trackerIso[i];
    myHcalIso[i] = hcalIso[i];
    myEcalJIso[i] = ecalJIso[i];
    myEcalGTIso[i] = ecalGTIso[i];
    myCombinedIso[i] = combinedIso[i];
    myCharge[i] = charge[i];
    myMissHits[i] = missHits[i];
    myDist[i] = dist[i];
    myDcot[i] = dcot[i];
    myLh[i] = lh[i];
    myE9ESC[i] = e9esc[i];
    myE25ESC[i] = e25esc[i];
  }
}

void RedVecbosTree::fillMoreElectrons(float esc[2], float eseed[2], float pin[2]) {
  for(int i=0; i<2; i++) {
    myESC[i] = esc[i];
    myEseed[i] = eseed[i];
    myPin[i] = pin[i];
  }
}

void RedVecbosTree::fillPFelectrons(float pfPt, float pfEta, float pfPhi, float pfdR, float pfmva, int npfe) { 
  myPFelePt     = pfPt;
  myPFeleEta    = pfEta;
  myPFelePhi    = pfPhi;
  myPFeleDeltaR = pfdR;
  myPFeleMva    = pfmva;
  myNPFele      = npfe;
}

void RedVecbosTree::fillKinematics(float invmass, float tmass, float tctmass, float pftmass, float met, float tcmet, float pfmet, float ptW, float ptWrecoil, float ptZ, float dijetmass[2], float dijetpt[2], float deltaetajet[2], float deltaphimet[2], float etaW, float etaZ, float mTZMet, float mZJet) {
  
  myInvMass            = invmass;
  myTransvMass         = tmass;
  myTCTransvMass       = tctmass;
  myPFTransvMass       = pftmass;
  myMet                = met;
  myTCMet              = tcmet;
  myPFMet              = pfmet;
  myPtW                = ptW;
  myPtZ                = ptZ;
  myPtWRecoil          = ptWrecoil;
  myDiJetMass[0]       = dijetmass[0];
  myDiJetMass[1]       = dijetmass[1];
  myDiJetPt[0]         = dijetpt[0];
  myDiJetPt[1]         = dijetpt[1];
  myDeltaEtaJet[0]     = deltaetajet[0];
  myDeltaEtaJet[1]     = deltaetajet[1];
  myDeltaPhiJetMet[0]  = deltaphimet[0];
  myDeltaPhiJetMet[1]  = deltaphimet[1];
  myEtaW               = etaW;
  myEtaZ               = etaZ;
  myMTZMet             = mTZMet;
  myMZJet              = mZJet;
}

void RedVecbosTree::fillWTemplates(float metCorr, float metUncorr, float mTCorr, float mTUncorr) {
  myTemplMetCorr      = metCorr;
  myTemplMetUncorr    = metUncorr;
  myTemplTrMassCorr   = mTCorr;
  myTemplTrMassUncorr = mTUncorr;
}

void RedVecbosTree::fillEventShape(std::vector<float> mhtphijet, std::vector<float> mhtphiPFjet, float mhtphimet) {
  assert(mhtphijet.size()==6 && mhtphiPFjet.size()==6);
  for(int i=0; i<6; i++) {
    myMHTphiJet[i]    = mhtphijet[i];    
    myMHTphiPFJet[i]  = mhtphiPFjet[i];
  }
  myMHTphiMET          = mhtphimet;    
}

void RedVecbosTree::fillMcTruth(int promptdecay, int nInteractions) {
  myPromptDecay = promptdecay;
  myNumInteractions = nInteractions;
}

void RedVecbosTree::fillBTagEVT(float combinedSecondaryVertexBJetTags, float combinedSecondaryVertexMVABJetTags,
                                float jetBProbabilityBJetTags, float jetProbabilityBJetTags, float simpleSecondaryVertexBJetTags,
                                float softMuonBJetTags,
                                float trackCountingHighPurBJetTags, float trackCountingHighEffBJetTags, 
				int fB, int fBm, int nB) {

  myCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags;
  myCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags;
  myJetBProbabilityBJetTags = jetBProbabilityBJetTags;
  myJetProbabilityBJetTags = jetProbabilityBJetTags;
  mySimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags;
  mySoftMuonBJetTags = softMuonBJetTags;
  myTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags;
  myTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags;
  myFoundMcB = fB;
  myFoundMcBmum = fBm;
  myNumB = nB;
}

void RedVecbosTree::fillBTagJetMultiplicities(int nJetsBTagged[5]) {
  for(int i=0; i<5; i++) myNJetsBTagged[i] = nJetsBTagged[i];
}

void RedVecbosTree::fillRunInfo(int run, int lumi, int event) {
  myRun = run;
  myLumi = lumi;
  myEvent = event;
}

void RedVecbosTree::fillEventInfo(int nvtx, float rho) {
  myNumVtx = nvtx;
  myRho = rho;
}

void RedVecbosTree::fillGenInfo(float p1, float p2, float e1, float e2, float mm) {
  myPtGen1  = p1;
  myPtGen2  = p2;
  myEtaGen1 = e1;
  myEtaGen2 = e2;
  myMeeGen  = mm;
}
