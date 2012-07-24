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
#include "RazorLeptons.hh"

RazorLeptons::RazorLeptons(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
}

RazorLeptons::RazorLeptons(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

RazorLeptons::~RazorLeptons() {}

void RazorLeptons::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorLeptons::SetWeight(double weight) {
  _weight = weight;
}

void RazorLeptons::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  cout << "ciaoo" << endl;

  // Ele Razor Triggers
  int HLTMuMu;
  int HLTEleEle;
  int HLTMuEle;

  double R;
  double RSQ;
  double MR;
  
  // general event info
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double W;
  int BOX_NUM;
  int nSelPFJets;

  // prepare the output tree
  TTree* outTree = new TTree("EVENTS", "EVENTS");
  
  outTree->Branch("HLTMuMu", &HLTMuMu, "HLTMuMu/I");
  outTree->Branch("HLTMuEle", &HLTMuEle, "HLTMuEle/I");
  outTree->Branch("HLTEleEle", &HLTEleEle, "HLTEleEle/I");
  
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("W", &W, "W/D");
  outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
  outTree->Branch("nSelPFJets", &nSelPFJets, "nSelPFJets/I");
  
  // kinematic block
  outTree->Branch("R", &R, "R/D");
  outTree->Branch("RSQ", &RSQ, "RSQ/D");
  outTree->Branch("MR", &MR, "MR/D");
  
  outTree->Branch("pTLep1", &pTLep1, "pTLep1/D");
  outTree->Branch("etaLep1", &etaLep1, "etaLep1/D");
  outTree->Branch("phiLep1", &phiLep1, "phiLep1/D");
  outTree->Branch("idLep1", &idLep1, "idLep1/I");
  
  outTree->Branch("pTLep2", &pTLep2, "pTLep2/D");
  outTree->Branch("etaLep2", &etaLep2, "etaLep2/D");
  outTree->Branch("phiLep2", &phiLep2, "phiLep2/D");
  outTree->Branch("idLep2", &idLep2, "idLep2/I");
  
  //  double _weight = 1.;
  unsigned int lastLumi = 0;
  unsigned int lastRun = 0;

  std::vector<std::string> maskHLTMuMu; 
  maskHLTMuMu.push_back("HLT_DoubleMu7_v");
  maskHLTMuMu.push_back("HLT_Mu13_Mu8_v");
  maskHLTMuMu.push_back("HLT_Mu17_Mu8_v");
  maskHLTMuMu.push_back("HLT_Mu17_TkMu8_v");

  //maskHLTMuMu.push_back("1-165208:HLT_DoubleMu7_v");
  //maskHLTMuMu.push_back("165364-178419:HLT_Mu13_Mu8_v");
  //maskHLTMuMu.push_back("178420-999999:HLT_Mu17_Mu8_v");
  //maskHLTMuMu.push_back("175832-999999:HLT_Mu17_TkMu8_v");

  std::vector<std::string> maskHLTMuEle; 
  maskHLTMuEle.push_back("HLT_Mu17_Ele8_CaloIdL_v");
  maskHLTMuEle.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v");
  maskHLTMuEle.push_back("HLT_Mu8_Ele17_CaloIdL_v");
  maskHLTMuEle.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v");

  //maskHLTMuEle.push_back("1-175972:HLT_Mu17_Ele8_CaloIdL_v");
  //maskHLTMuEle.push_back("175973-999999:HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v");
  //maskHLTMuEle.push_back("1-167913:HLT_Mu8_Ele17_CaloIdL_v");
  //maskHLTMuEle.push_back("167914-999999:HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v");

  std::vector<std::string> maskHLTEleEle; 
  maskHLTEleEle.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");
  maskHLTEleEle.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

  //maskHLTEleEle.push_back("1-170901:HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");
  //maskHLTEleEle.push_back("171050-999999:HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

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
      setRequiredTriggers(maskHLTMuMu);   reloadTriggerMask(true); HLTMuMu   = hasPassedHLT();
      setRequiredTriggers(maskHLTMuEle);  reloadTriggerMask(true); HLTMuEle  = hasPassedHLT();
      setRequiredTriggers(maskHLTEleEle); reloadTriggerMask(true); HLTEleEle = hasPassedHLT();
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
    
    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
        if(nPV<1) continue;
    for(int i=0; i < nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 


    int goodPFevent = true;
    int goodCaloevent = true;

    // HCAL FLAGS
    if(_isData && !eventPassHcalFilter()) continue;
    //ECALTPFilterFlag
    //    if(_isData && METFlags << 0 == 0) continue;
    //drBoundary      
    //    if(_isData && METFlags << 1 == 0) continue;
    // drDead         
    //    if(_isData && METFlags << 2 == 0) continue;
    // CSCHaloFilterFlag
    //    if(_isData && METFlags << 3 == 0) continue;
    // trackerFailureFilterFlag
    //    if(_isData && METFlags << 4 == 0) continue;
    // BE ECAL flag            
    //    if(_isData && METFlags << 5 == 0) continue;

    // look for the leptons
    TLorentzVector myMuLoose(0.1, 0., 0., 0.1);
    TLorentzVector myMuTight(0.1, 0., 0., 0.1);
    int iMuLoose = -99;
    int iMuTight = -99;

    for(int i=0; i<nMuon; i++) {
      if(muonPassTight(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 20.) {
	  if(thisMu.Pt() > myMuTight.Pt()) {
	    myMuTight = thisMu;
	    iMuTight = i;
	  }
	}
      }
    }

    for(int i=0; i<nMuon; i++) {
      if(muonPassTight(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 20. && iMuTight != i) {
	  if(thisMu.Pt() > myMuLoose.Pt()) {
	    myMuLoose = thisMu;
	    iMuLoose = i;
	  }
	}
      }
    }

    TLorentzVector myEle1WP80(0.1, 0., 0., 0.1);
    TLorentzVector myEle2WP80(0.1, 0., 0., 0.1);
    int iEle1WP80 = -99;
    int iEle2WP80 = -99;

    for(int i=0; i<nEle; i++) {
      if(electronPassWP80(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 20.) {
	  if(thisEle.Pt() > myEle1WP80.Pt()) {
	    myEle1WP80 = thisEle;
	    iEle1WP80 = i;
	  }
	}
      }
    }

    for(int i=0; i<nEle; i++) {
      if(electronPassWP80(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 20. && i != iEle2WP80) {
	  if(thisEle.Pt() > myEle2WP80.Pt()) {
	    myEle2WP80 = thisEle;
	    iEle2WP80 = i;
	  }
	}
      }
    }

    int nGoodLep = 0;
    if(iMuTight>=0) nGoodLep++;
    if(iMuLoose>=0) nGoodLep++;
    if(iEle1WP80>-0) nGoodLep++;
    if(iEle2WP80>-0) nGoodLep++;

    if(nGoodLep<2) continue;

    // Fill the tree per box
    pTLep1 = -999.;
    etaLep1 = -999.;
    phiLep1 = -999.;
    idLep1 = -999.;

    pTLep2 = -999.;
    etaLep2 = -999.;
    phiLep2 = -999.;
    idLep2 = -999.;    

    BOX_NUM = -99;
    // mue  0
    // mumu 1 
    // ee   2 

    // MUELE BOX
    if(iEle1WP80 >= 0 && iMuTight >=0) { 
      FirstLepton(myMuTight, chargeMuon[iMuTight]*13);
      SecondLepton(myEle1WP80, chargeEle[iEle1WP80]*11);
      BOX_NUM = 0;
    } else if(iEle1WP80 >= 0 && iMuLoose >=0) {
      FirstLepton(myMuLoose, chargeMuon[iMuLoose]*13);
      SecondLepton(myEle1WP80, chargeEle[iEle1WP80]*11);
      BOX_NUM = 0;
    } else if(iMuTight >=0 && iMuLoose >=0) {
      FirstLepton(myMuTight, chargeMuon[iMuTight]*13);
      SecondLepton(myMuLoose, chargeMuon[iMuLoose]*13);
      BOX_NUM = 1;
    } else if(iEle1WP80 >= 0 && iEle2WP80 >= 0) {
      FirstLepton(myEle1WP80, chargeEle[iEle1WP80]*11);
      SecondLepton(myEle2WP80, chargeEle[iEle2WP80]*11);
      BOX_NUM = 2;
    }

    if(BOX_NUM < 0) {
      cout << "PROBLEM HERE!!!!!!!!" << endl;
    }
    
    TLorentzVector PFLep1;
    TLorentzVector PFLep2;
    PFLep1.SetPtEtaPhiM(pTLep1, etaLep1, phiLep1, 0.0);
    PFLep2.SetPtEtaPhiM(pTLep2, etaLep2, phiLep2, 0.0);
    
    // dummy values
    R = -99999.;
    RSQ = -99999.;
    MR = -99999.;
    
    // use PFMET
    TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);
    
    // compute boost
    double num = PFLep1.P()-PFLep2.P();
    double den = PFLep1.Pz()-PFLep2.Pz();      
    
    double MT = CalcMTR(PFLep1, PFLep2, MET);
    double variable = -999999.;
    double Rvariable = -999999.;
    variable = CalcGammaMRstar(PFLep1, PFLep2);
    if(variable > 0) Rvariable = MT/variable;
    
    // fill the tree
    R = Rvariable;
    RSQ = Rvariable*Rvariable;
    MR = variable;    

    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;
    W = _weight;
    
    outTree->Fill();
    
  }
  
  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  file->Close();
}
  
void RazorLeptons::FirstLepton(TLorentzVector lep, int id) {
  pTLep1 = lep.Pt();
  etaLep1 = lep.Eta();
  phiLep1 = lep.Phi();
  idLep1 = id;
}

void RazorLeptons::SecondLepton(TLorentzVector lep, int id) {
  pTLep1 = lep.Pt();
  etaLep1 = lep.Eta();
  phiLep1 = lep.Phi();
  idLep1 = id;
}
