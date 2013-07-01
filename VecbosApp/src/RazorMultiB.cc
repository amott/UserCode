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
#include "AnalysisSelector.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "RazorMultiB.hh"

RazorMultiB::RazorMultiB(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

RazorMultiB::RazorMultiB(TTree *tree, string jsonFile, bool goodRunLS, bool isData) : Vecbos(tree) {
  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;
  _isSMS = false;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
    InitEventFlag("/afs/cern.ch/user/w/woodson/public/WEIGHT/AllBadABCDNEWTAUID.txt");
  }

}


RazorMultiB::RazorMultiB(TTree *tree, string jsonFile, bool goodRunLS, bool isData, string smsName) : Vecbos(tree) {
  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;
  _isSMS = false;
  if (smsName!="none") _isSMS = true;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
    InitEventFlag("/afs/cern.ch/user/w/woodson/public/WEIGHT/AllBadABCDNEWTAUID.txt");
  }
}

RazorMultiB::~RazorMultiB() {}

void RazorMultiB::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorMultiB::SetWeight(double weight){
  _weight=weight;
}

void RazorMultiB::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  int HLT_DoubleMu;
  int HLT_DoubleEle;
  int HLT_MuEle;

  bool ECALTPFilterFlag;
  bool drBoundary;
  bool drDead;
  bool CSCHaloFilterFlag;
  bool trackerFailureFilterFlag;
  bool BEECALFlag;

  // PF block
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double PFR;
  double RSQ;
  double MR;
  double MRT;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  int    nBtag;
  int    nBtag_medium;
  int    nBtag_tight;
  double W=_weight;
  float mg, mchi;
  int    BOX_NUM;
  int    ss;
  int    nPV;
  
  // fast-hemispheres
  double FJR;
  double FJRSQ;
  double FJMR;
  double FJMRT;
  double pTFJHem1;
  double etaFJHem1;
  double phiFJHem1;
  double pTFJHem2;
  double etaFJHem2;
  double phiFJHem2;

  // gen level info
  double pT1, pT2, eta1, eta2, phi1, phi2;
  int idMc1, idMothMc1, idGrandMothMc1;
  int idMc2, idMothMc2, idGrandMothMc2;
  // ttbar decay: 0 = nolep, 1 = semilep; 2 = fully lep
  int nLepTopDecay;

  // PFElectron Block
  double pfElectron_pt;
  double pfElectron_eta;
  double pfElectron_phi;
  double pfElectron_energy;
  double pfElectron_mass = 0.511/1000;

  // PFMuon Block
  double pfMuon_pt;
  double pfMuon_eta;
  double pfMuon_phi;
  double pfMuon_energy;
  double pfMuon_mass = 105.7/1000;

  int nIsolatedPFJets;
  int nPFJets;
  int nBtag_lead4jets;
  int nBtag_medium_lead4jets;
  int nBtag_tight_lead4jets;
  int nBtag_TCHPT;
  int nBtag_TCHPT_lead4jets;
  double Mll;

  // New Razor Variables
  double MR_pTcorr;
  double gammaR;
  double shatR_bl;
  double dPhiCM;
  double EB1;
  double EB2;
  double CosThetaB1;
  double CosThetaB2;
  double TopMass1;
  double TopMass2;
  double EL1;
  double EL2;
  double CosThetaL1;
  double CosThetaL2;
  double GluinoMass1;
  double GluinoMass2;
  double VisHemMass1;
  double VisHemMass2;
  double TotalHemMass1;
  double TotalHemMass2;
  double TopHemMass1;
  double TopHemMass2;
  double MR1;
  double MR2;
  //here are the ones including trasverse masses
  double TotalHemMass1Trans;
  double TotalHemMass2Trans;
  double TopHemMass1Trans;
  double TopHemMass2Trans;
  double MR1Trans;
  double MR2Trans;
  
  //	
  double MetMag;  	
  double HT;
  double MHT_x;
  double MHT_y;
  double MHT;  

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
  outTree->Branch("ss", &ss, "ss/I");
  
  // HLT bits
  outTree->Branch("HLT_DoubleMu", &HLT_DoubleMu, "HLT_DoubleMu/I"); 
  outTree->Branch("HLT_DoubleEle", &HLT_DoubleEle, "HLT_DoubleEle/I"); 
  outTree->Branch("HLT_MuEle", &HLT_MuEle, "HLT_MuEle/I"); 

  // PF block
  outTree->Branch("PFR", &PFR, "PFR/D");
  outTree->Branch("RSQ", &RSQ, "RSQ/D");
  outTree->Branch("MR", &MR, "MR/D");
  outTree->Branch("MRT", &MRT, "MRT/D");
  outTree->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
  outTree->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
  outTree->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
  outTree->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
  outTree->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
  outTree->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
  outTree->Branch("nBtag", &nBtag, "nBtag/I");
  outTree->Branch("nBtag_medium", &nBtag_medium, "nBtag_medium/I");
  outTree->Branch("nBtag_tight", &nBtag_tight, "nBtag_tight/I");
  outTree->Branch("nBtag_TCHPT", &nBtag_TCHPT, "nBtag_TCHPT/I");
  outTree->Branch("W", &W, "W/D");
  outTree->Branch("mg", &mg, "mg/F");
  outTree->Branch("mchi", &mchi, "mchi/F");
  outTree->Branch("nPV", &nPV, "nPV/I");

  // fast-hemispheres
  outTree->Branch("FJR", &FJR, "FJR/D");
  outTree->Branch("FJRSQ", &FJRSQ, "FJRSQ/D");
  outTree->Branch("FJMR", &FJMR, "FJMR/D");
  outTree->Branch("FJMRT", &FJMRT, "FJMRT/D");
  outTree->Branch("pTFJHem1", &pTFJHem1, "pTFJHem1/D");
  outTree->Branch("etaFJHem1", &etaFJHem1, "etaFJHem1/D");
  outTree->Branch("phiFJHem1", &phiFJHem1, "phiPFJHem1/D");
  outTree->Branch("pTFJHem2", &pTFJHem2, "pTFJHem2/D");
  outTree->Branch("etaFJHem2", &etaFJHem2, "etaFJHem2/D");
  outTree->Branch("phiFJHem2", &phiFJHem2, "phiFJHem2/D");
  

  // New Razor Variables
  outTree->Branch("MR_pTcorr", &MR_pTcorr, "MR_pTcorr/D");
  outTree->Branch("gammaR", &gammaR, "gammaR/D");
  outTree->Branch("dPhiCM", &dPhiCM, "dPhiCM/D");
  outTree->Branch("shatR_bl", &shatR_bl, "shatR_bl/D");
  outTree->Branch("EB1", &EB1, "EB1/D");
  outTree->Branch("EB2", &EB2, "EB2/D");
  outTree->Branch("EL1", &EL1, "EL1/D");
  outTree->Branch("EL2", &EL2, "EL2/D");
  outTree->Branch("CosThetaB1", &CosThetaB1, "CosThetaB1/D");
  outTree->Branch("CosThetaB2", &CosThetaB2, "CosThetaB2/D");
  outTree->Branch("CosThetaL1", &CosThetaL1, "CosThetaL1/D");
  outTree->Branch("CosThetaL2", &CosThetaL2, "CosThetaL2/D");
  outTree->Branch("TopMass1", &TopMass1, "TopMass1/D");
  outTree->Branch("TopMass2", &TopMass2, "TopMass2/D");
  outTree->Branch("GluinoMass1", &GluinoMass1, "GluinoMass1/D");
  outTree->Branch("GluinoMass2", &GluinoMass2, "GluinoMass2/D");
  outTree->Branch("VisHemMass1", &VisHemMass1, "VisHemMass1/D");
  outTree->Branch("VisHemMass2", &VisHemMass2, "VisHemMass2/D");
  outTree->Branch("TotalHemMass1", &TotalHemMass1, "TotalHemMass1/D");
 outTree->Branch("TotalHemMass2", &TotalHemMass2, "TotalHemMass2/D");
 outTree->Branch("TopHemMass1", &TopHemMass1, "TopHemMass1/D");
 outTree->Branch("TopHemMass2", &TopHemMass2, "TopHemMass2/D");
 outTree->Branch("MR1", &MR1, "MR1/D");
 outTree->Branch("MR2", &MR2, "MR2/D");
	outTree->Branch("TotalHemMass1Trans", &TotalHemMass1Trans, "TotalHemMass1Trans/D");
	outTree->Branch("TotalHemMass2Trans", &TotalHemMass2Trans, "TotalHemMass2Trans/D");
	outTree->Branch("TopHemMass1Trans", &TopHemMass1Trans, "TopHemMass1Trans/D");
	outTree->Branch("TopHemMass2Trans", &TopHemMass2Trans, "TopHemMass2Trans/D");
	outTree->Branch("MR1Trans", &MR1Trans, "MR1Trans/D");
	outTree->Branch("MR2Trans", &MR2Trans, "MR2Trans/D");
 outTree->Branch("MetMag", &MetMag, "MetMag/D");
 outTree->Branch("HT", &HT, "HT/D");
 outTree->Branch("MHT", &MHT, "MHT/D");
 outTree->Branch("MHT_x", &MHT_x, "MHT_x/D");
 outTree->Branch("MHT_y", &MHT_y, "MHT_y/D");
 

  //Gen-Level
  outTree->Branch("idMc1", &idMc1, "idMc1/I");
  outTree->Branch("idMothMc1", &idMothMc1, "idMothMc1/I");
  outTree->Branch("idGrandMothMc1", &idGrandMothMc1, "idGrandMothMc1/I");
  outTree->Branch("pT1", &pT1, "pT1/D");
  outTree->Branch("eta1", &eta1, "eta1/D");
  outTree->Branch("phi1", &phi1, "phi1/D");
  outTree->Branch("idMc2", &idMc2, "idMc2/I");
  outTree->Branch("idMothMc2", &idMothMc2, "idMothMc2/I");
  outTree->Branch("idGrandMothMc2", &idGrandMothMc2, "idGrandMothMc2/I");
  outTree->Branch("pT2", &pT2, "pT2/D");
  outTree->Branch("eta2", &eta2, "eta2/D");
  outTree->Branch("phi2", &phi2, "phi2/D");
  outTree->Branch("nLepTopDecay",&nLepTopDecay,"nLepTopDecay/I");

  //Selection
  outTree->Branch("nIsolatedPFJets", &nIsolatedPFJets, "nIsolatedPFJets/I");
  outTree->Branch("nPFJets", &nPFJets, "nPFJets/I");
  outTree->Branch("nBtag_lead4jets", &nBtag_lead4jets, "nBtag_lead4jets/I");
  outTree->Branch("nBtag_medium_lead4jets", &nBtag_medium_lead4jets, "nBtag_medium_lead4jets/I");
  outTree->Branch("nBtag_tight_lead4jets", &nBtag_tight_lead4jets, "nBtag_tight_lead4jets/I");
  outTree->Branch("nBtag_TCHPT_lead4jets", &nBtag_TCHPT_lead4jets, "nBtag_TCHPT_lead4jets/I");
  outTree->Branch("Mll", &Mll, "Mll/D");
  
  double Npassed_In = 0;
  double Npassed_PV = 0;
  //Jets
  double NpassedPF_4Jet = 0;
  //Leptons
  double Npassed_Lept=0;
  //B-tag
  double Npassed_1b=0;
  //Mll Cut
  double Npassed_Mll = 0;
  double Z_mass = 91.1876;

  double weightII = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  std::vector<std::string> maskHLT_DoubleMu; 
  //maskHLT_DoubleMu.push_back("HLT_Mu13_Mu8_v");
  maskHLT_DoubleMu.push_back("HLT_Mu17_Mu8_v");
  
  std::vector<std::string> maskHLT_DoubleEle; 
  maskHLT_DoubleEle.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  
  std::vector<std::string> maskHLT_MuEle; 
  maskHLT_MuEle.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  maskHLT_MuEle.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries= " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;


    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      setRequiredTriggers(maskHLT_DoubleMu); reloadTriggerMask(true); HLT_DoubleMu = hasPassedHLT();
      setRequiredTriggers(maskHLT_DoubleEle); reloadTriggerMask(true); HLT_DoubleEle = hasPassedHLT();
      setRequiredTriggers(maskHLT_MuEle); reloadTriggerMask(true); HLT_MuEle = hasPassedHLT();
    }

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    // Kill bad events in Data
    if(_isData && FailFilters()) continue;
    if(_isData && isFlagged()) continue;
    
    Npassed_In += weightII;
    
    double m0=9999, m12=9999, mc=9999;
     


    if((!_isData) && _isSMS){
      std::vector<float> parameterPoint = ParseEvent();
      mg = parameterPoint[0];
      mchi = parameterPoint[1];
    }
    
    //HLT
    int passedHLT = HLT_DoubleMu+HLT_DoubleEle+HLT_MuEle;
    if (_isData==true) {
      if (passedHLT==0) continue;
      if ((ECALTPFilterFlag==0) || (drBoundary==0) || (drDead==0) || (CSCHaloFilterFlag==0) || (trackerFailureFilterFlag==0) || (BEECALFlag==0)) continue;
    }
    
   // find highest-pT PV [replace with Sagar's code]
    int iPV = passPV();
    //if(iPV<0) continue;
    Npassed_PV += weightII;
    nPV = N_PV_EVENT;

    nIsolatedPFJets = 0;
    nPFJets = 0;
    vector <TLorentzVector> PFJet;
    vector <int> iPFJet;
    vector <TLorentzVector> IsolatedPFJet;
    vector <int> iIsolatedPFJet;
    bool badjet = false;


    for(int i = 0; i < nAK5PFNoPUJet; i++){
      TLorentzVector myJet;
      double px = pxAK5PFNoPUJet[i];
      double py = pyAK5PFNoPUJet[i];
      double pz = pzAK5PFNoPUJet[i];
      double E = sqrt(px*px+py*py+pz*pz);
      myJet.SetPxPyPzE(px,py,pz,E);
  
      if(myJet.Pt() > 30.0 && fabs(myJet.Eta()) < 3.0){
	if ( goodJetID(i) ) {
	  PFJet.push_back(myJet);
	  iPFJet.push_back(i);
	  nPFJets += 1;
	  // Examine if PFJet is isolated from leptons or not and Mll pass the Z mass threshold
	
	} else {
	  PFJet.clear();
	  iPFJet.clear();
	  nPFJets = 0.;
	  break;
	}
      }
    }

    for (int i = 0; i < PFJet.size(); i++){
      TLorentzVector myJet = PFJet[i];
      int iJet = iPFJet[i];
      // Examine if PFJet is isolated from leptons or not and Mll pass the Z mass threshold
      bool isIsolatedJet = true;	
      for(int j=0; j<nMuon; j++) {
	TLorentzVector thisMu(pxMuon[j], pyMuon[j], pzMuon[j], energyMuon[j]);
	if(isTightMuon(j) && (thisMu.Pt()>20.) && (myJet.DeltaR(thisMu)<=0.3)) isIsolatedJet = false;
      }	  
      for(int k=0; k<nEle; k++) {
	TLorentzVector thisEle(pxEle[k], pyEle[k], pzEle[k], energyEle[k]);
	if(isTightElectron(k) && (thisEle.Pt()>20.) && (myJet.DeltaR(thisEle)<=0.3)) isIsolatedJet = false;
      }
      if(isIsolatedJet) {
	IsolatedPFJet.push_back(myJet);
	iIsolatedPFJet.push_back(iJet);
	nIsolatedPFJets += 1;
      }
    }
  
  
    // jet ID
    if (nPFJets==0) continue;

    //4Jet
    /*
    if (isDeltaRIsolated) {
      //if(iIsolatedPFJet.size()<4) continue;
      NpassedPF_4Jet+=weightII;
    } else {
      //if(iPFJet.size() <4) continue;
      NpassedPF_4Jet+=weightII;
    }
    */

    // at least 2 jets
    if (nPFJets < 2) continue;
    
    int flag = 1;
    //bubble sort the jets by Pt()
    for (int i=0; (i < iIsolatedPFJet.size()) && flag; i++){
      TLorentzVector tempvector;
      int tempi;
      flag = 0;
      for (int j=0; j < (iIsolatedPFJet.size()-1); j++){
	if (IsolatedPFJet.at(j+1).Pt() > IsolatedPFJet.at(j).Pt()){
	  tempvector = IsolatedPFJet.at(j);
	  IsolatedPFJet.at(j) = IsolatedPFJet.at(j+1);
	  IsolatedPFJet.at(j+1) = tempvector;
	  
	  tempi = iIsolatedPFJet.at(j);
	  iIsolatedPFJet.at(j) = iIsolatedPFJet.at(j+1);
	  iIsolatedPFJet.at(j+1) = tempi;
	  flag=1;	    
	}
      }  
    }
    
    // b from leading 4 isolated jets
    nBtag_lead4jets = 0;
    nBtag_medium_lead4jets = 0;
    nBtag_tight_lead4jets = 0;
    nBtag_TCHPT_lead4jets = 0;
    if (nIsolatedPFJets >=4) {
      for(int b=0; b< iIsolatedPFJet.size(); b++){
	int n=iIsolatedPFJet.at(b);
	if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n]) && b < 4) nBtag_lead4jets++;  
	if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n]) && b < 4) nBtag_medium_lead4jets++; 
	if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n]) && b < 4) nBtag_tight_lead4jets++;  
	if(trackCountingHighPurBJetTagsAK5PFNoPUJet[n] > 3.41 && b < 4) nBtag_TCHPT_lead4jets++;  
	}
    }    
        
    //Create arrays with the b-discriminators of the jets
    nBtag = 0;
    nBtag_medium = 0;
    nBtag_tight = 0;
    nBtag_TCHPT = 0;
    for(int b=0; b< iIsolatedPFJet.size(); b++){
      int n=iIsolatedPFJet.at(b);
      if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag++; 
      if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag_medium++;
      if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag_tight++;  
      if(trackCountingHighPurBJetTagsAK5PFNoPUJet[n] > 3.41) nBtag_TCHPT++;  
    }
    // 1b requirement
    //if(nBtag == 0) continue; 
    Npassed_1b+=weightII;
    
    // Boxes
    vector<int> iMuTight;
    vector<TLorentzVector> MuTight;
    for(int i=0; i<nMuon; i++) {
      TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
      if(isTightMuon(i) and thisMu.Pt() > 20.) {
	iMuTight.push_back(i);
	MuTight.push_back(thisMu);
      }	  
    }

    vector<int> iEleTight;
    vector<TLorentzVector> EleTight;
    for(int i=0; i<nEle; i++) {
      TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
      if(isTightElectron(i) && thisEle.Pt() > 20.) {
	iEleTight.push_back(i);
	EleTight.push_back(thisEle);
      }
    }
    
    //bubble sort the muons by Pt()
    flag = 1;
    if (iMuTight.size() > 1) {
      for (int i=0; (i < iMuTight.size()) && flag; i++){
	TLorentzVector tempvector;
	int tempi;
	flag = 0;
	for (int j=0; j < (iMuTight.size()-1); j++){
	  if (MuTight.at(j+1).Pt() > MuTight.at(j).Pt()){
	    tempvector = MuTight.at(j);
	  MuTight.at(j) = MuTight.at(j+1);
	  MuTight.at(j+1) = tempvector;
	  
	  tempi = iMuTight.at(j);
	  iMuTight.at(j) = iMuTight.at(j+1);
	  iMuTight.at(j+1) = tempi;
	  flag=1;	    
	  }
	}
      }    
    }   
 
    //bubble sort the electrons by Pt()
    flag = 1;
    if (iEleTight.size() > 1) {
      for (int i=0; (i < iEleTight.size()) && flag; i++){
	TLorentzVector tempvector;
	int tempi;
	flag = 0;
	for (int j=0; j < (iEleTight.size()-1); j++){
	  if (EleTight.at(j+1).Pt() > EleTight.at(j).Pt()){
	    tempvector = EleTight.at(j);
	    EleTight.at(j) = EleTight.at(j+1);
	    EleTight.at(j+1) = tempvector;
	    
	    tempi = iEleTight.at(j);
	    iEleTight.at(j) = iEleTight.at(j+1);
	    iEleTight.at(j+1) = tempi;
	    flag=1;	    
	  }
	}
      }        
    }

    // Determine SS/OS for each dilepton events and veto on Z_mass window of width 10GeV
    int iLepton1;
    int iLepton2;    

    vector<TLorentzVector> DiLepton;
    BOX_NUM = -99;
    if(iMuTight.size() >0 && iEleTight.size()>0){
      iLepton1 = iMuTight[0];     
      iLepton2 = iEleTight[0];
      DiLepton.push_back(MuTight[0]);
      DiLepton.push_back(EleTight[0]);
      Mll = (MuTight[0]+EleTight[0]).M();
      ss = chargeMuon[iLepton1]*chargeEle[iLepton2];
      BOX_NUM = 0;
    } else if(iMuTight.size()>1) {
      iLepton1 = iMuTight[0];  
      iLepton2 = iMuTight[1];
      DiLepton.push_back(MuTight[0]);
      DiLepton.push_back(MuTight[1]);
      Mll = (MuTight[0]+MuTight[1]).M();
      ss = chargeMuon[iLepton1]*chargeMuon[iLepton2];
      BOX_NUM = 1;
    } else if(iEleTight.size()>1) {
      iLepton1 = iEleTight[0];
      iLepton2 = iEleTight[1];
      DiLepton.push_back(EleTight[0]);
      DiLepton.push_back(EleTight[1]);
      Mll = (EleTight[0]+EleTight[1]).M();
      ss = chargeEle[iLepton1]*chargeEle[iLepton2];
      BOX_NUM = 2;
    }
    Npassed_Mll += weightII;

    // two good leptons
    if(BOX_NUM<0) continue;
    Npassed_Lept += weightII;
    

    // create vector of b-jets
    vector<TLorentzVector> BJets;
    vector<float> discBJets;
    vector<int> iBJets;
    vector<TLorentzVector> otherJets;
    vector<int> iotherJets;

    for(int b=0; b < iIsolatedPFJet.size(); b++){
      int n = iIsolatedPFJet.at(b);
      float disc = combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n];
      TLorentzVector tempvector;
      tempvector.SetPxPyPzE(pxAK5PFNoPUJet[n],pyAK5PFNoPUJet[n],pzAK5PFNoPUJet[n],energyAK5PFNoPUJet[n]);
      
      if (disc > 0.244) {
	iBJets.push_back(n);
	discBJets.push_back(disc);
	BJets.push_back(tempvector);
      }
    }
    

    // bubble sort the b-jets by CSV b-tag discrimator values
    flag = 1;
    for(int i=0; (i < iBJets.size()) && flag; i++){
	TLorentzVector tempvector;
	int tempi;
	float tempdisc;
	flag = 0;
	for (int j=0; j < (iBJets.size()-1); j++){
	  if (discBJets.at(j+1) > discBJets.at(j)){
	    tempvector  = BJets.at(j);
	    BJets.at(j) = BJets.at(j+1);
	    BJets.at(j+1) = tempvector;

	    tempi = iBJets.at(j);
	    iBJets.at(j) = iBJets.at(j+1);
	    iBJets.at(j+1) = tempi;

	    tempdisc = discBJets.at(j);
	    discBJets.at(j) = discBJets.at(j+1);
	    discBJets.at(j+1) = tempdisc;
	    flag = 1;
	  }
	}
    }

    // New Razor Variables
    // dummy values 
    HT = -9999.;
    MHT_x = -9999.;
    MHT_y = -9999.;
    MHT = -9999.;
    MetMag = -9999.;
    EB1 = -9999.;
    EB2 = -9999.;
    EL1 = -9999.;
    EL2 = -9999.;
    dPhiCM = -9999.;
    gammaR = -9999.;
    TopMass1 = -9999.;
    TopMass2 = -9999.;
    CosThetaB1 = -9999.;
    CosThetaB2 = -9999.;
    CosThetaL1 = -9999.;
    CosThetaL2 = -9999.;    
    GluinoMass1 = -9999.;
    GluinoMass2 = -9999.;
    VisHemMass1 = -9999.;
    VisHemMass2 = -9999.;
    TotalHemMass1 = -9999.;
    TotalHemMass2 = -9999.;
    TopHemMass1 = -9999.;
    TopHemMass2 = -9999.;
    MR1 = -9999.;
    MR2 = -9999.;
    

    //dummy var
    double thisjet_x;
    double thisjet_y;
    thisjet_x = 0;
    thisjet_y = 0;


    if (DiLepton.size()>=2 && BJets.size()>=2){
      //Start by using PF MET
      TVector3 MET(pxPFMet[2], pyPFMet[2], 0.);
      MetMag = MET.Mag();      
      
      HT =0.; 
      MHT_x = 0.;
      MHT_y = 0.;
     

     TLorentzVector thisjet;


      for (int i=0; i < IsolatedPFJet.size(); i++) {
	
        thisjet = IsolatedPFJet[i];
	
	HT += sqrt(thisjet.X() * thisjet.X() + thisjet.Y() * thisjet.Y());
        MHT_x += thisjet.X();
        MHT_y += thisjet.Y();
	}	      
      MHT = sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
 	
      //two leptons
      TLorentzVector L1 = DiLepton[0];
      TLorentzVector L2 = DiLepton[1];
	
      //and two b's
      TLorentzVector B1, B2;
      int iB1, iB2;
      //first, we choose the hemisphere pairing of the b's and leptons by minimizing invariant masses

      float smallestM2 = 9999999999.;
      for (int i=0; i < BJets.size(); i++) {
	TLorentzVector testB1 = BJets[i];
	for (int j=i+1; j < BJets.size(); j++) {
	  TLorentzVector testB2 = BJets[j];
	  float M2 =  std::min( (B1+L1).M2() + (B2+L2).M2(), (B1+L2).M2() + (B2+L1).M2() );
	  if (M2 <= smallestM2) {
	    smallestM2 = M2;
	    B1 = testB1;
	    iB1 = iBJets[i];
	    B2 = testB2;
	    iB2 = iBJets[j];
	  }
	}
      }

      
      for (int i=0; i < IsolatedPFJet.size(); i++) {
	if (iIsolatedPFJet[i]!=iB1 && iIsolatedPFJet[i]!=iB2){
	  otherJets.push_back(IsolatedPFJet[i]);
	  iotherJets.push_back(iIsolatedPFJet[i]);
	}
      }

      if( (B1+L1).M2() + (B2+L2).M2() > (B1+L2).M2() + (B2+L1).M2() ){
	TLorentzVector temp = B1;
	B1 = B2;
	B2 = temp;
      }	
      
      // Combine jets keeping the Tops (b+l) together and in separate hemispheres:
       vector<TLorentzVector> Tops;
       TLorentzVector Top1 = (B1+L1);
       TLorentzVector Top2 = (B2+L2);
       Tops.push_back(Top1);
       Tops.push_back(Top2);

       vector<TLorentzVector> H12  = CombineJetsTs(otherJets, Tops);
       TLorentzVector H1 = H12[0];
       TLorentzVector H2 = H12[1];
		
      // Now, we must perform longitudinal boost from lab frame to CMz frame
      // in order to make procedure invariant under longitundinal boosts
      TVector3 BL = H1.Vect()+H2.Vect();
      BL.SetX(0.0);
      BL.SetY(0.0);
      BL = (1./(H1.E()+H2.E()))*BL;
	
      // Boost to CMz frame
      Top1.Boost(-BL);
      Top2.Boost(-BL);
      H1.Boost(-BL);
      H2.Boost(-BL);

      VisHemMass1 = H1.M();
      VisHemMass2 = H2.M();
    
      //Now, we must boost to each of the respective hemisphere rest frames
      TVector3 vBETA = Boost_type1(H1,H2);
    
      H1.Boost(-vBETA);
      H2.Boost(vBETA);
      Top1.Boost(-vBETA);
      Top2.Boost(vBETA);

      TLorentzVector N1, N2;
      N1.SetXYZM(H1.X(),H1.Y(),H1.Z(),0.0);
      N2.SetXYZM(H2.X(),H2.Y(),H2.Z(),0.0);

      TotalHemMass1 = (H1+N1).M();
      TotalHemMass2 = (H2+N2).M();

      TopHemMass1 = (Top1+N1).M();
      TopHemMass2 = (Top2+N1).M();
      
      MR1 = 2*H1.Vect().Mag();
      MR2 = 2*H2.Vect().Mag();

      
		
		
  // REDO ABOVE SECTION, INCLUDE TRANSVERSE BOOST AS WELL THOUGH
	
		// Combine jets keeping the Tops (b+l) together and in separate hemispheres:
		//vector<TLorentzVector> Tops;
		Top1 = (B1+L1);
		Top2 = (B2+L2);
		Tops.push_back(Top1);
		Tops.push_back(Top2);
		
		H12  = CombineJetsTs(otherJets, Tops);
		H1 = H12[0];
		H2 = H12[1];
		
		
		
		// Now, we must perform longitudinal boost from lab frame to CMz frame
		// in order to make procedure invariant under longitundinal boosts
		BL = H1.Vect()+H2.Vect();
		TVector3 BTr = H1.Vect() + H2.Vect();
		BL.SetX(0.0);
		BL.SetY(0.0);
		BL = (1./(H1.E()+H2.E()))*BL;
		
		// Boost to CMz frame
		Top1.Boost(-BL);
		Top2.Boost(-BL);
		H1.Boost(-BL);
		H2.Boost(-BL);
		
		// Transverse Boosts go here
		
		BTr.SetZ(0.0);
		Top1.Boost(-BTr);
		Top2.Boost(-BTr);
		H1.Boost(-BTr);
		H2.Boost(-BTr);
		//
		
		
		//continue as before
		//TVector3 vBETA = Boost_type1(H1,H2);
		
		H1.Boost(-vBETA);
		H2.Boost(vBETA);
		Top1.Boost(-vBETA);
		Top2.Boost(vBETA);
		
		N1, N2;
		N1.SetXYZM(H1.X(),H1.Y(),H1.Z(),0.0);
		N2.SetXYZM(H2.X(),H2.Y(),H2.Z(),0.0);
		
		TotalHemMass1Trans = (H1+N1).M();
		TotalHemMass2Trans = (H2+N2).M();
		
		TopHemMass1Trans = (Top1+N1).M();
		TopHemMass2Trans = (Top2+N1).M();
		
		MR1Trans = 2*H1.Vect().Mag();
		MR2Trans = 2*H2.Vect().Mag();
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	// Rogan's Approach to fully leptonic ttbar

      H1 = (B1+L1);
      H2 = (B2+L2);

      // Now, we must perform longitudinal boost from lab frame to CMz frame
      // in order to make procedure invariant under longitundinal boosts
      BL = H1.Vect()+H2.Vect();
      BL.SetX(0.0);
      BL.SetY(0.0);
      BL = (1./(H1.E()+H2.E()))*BL;
	
      // Boost to CMz frame
      B1.Boost(-BL);
      L1.Boost(-BL);
      B2.Boost(-BL);
      L2.Boost(-BL);

      // Now, we will perform a transverse boost from the CMz frame to our
      // approximation of the CM rest frame
     
    

      shatR_bl = shatR(H1+H2,MET);
      TVector3 betaTR = BetaTR(H1+H2,MET);

      // Boost to ~CM frame
      B1.Boost(-betaTR);
      B2.Boost(-betaTR);
      L1.Boost(-betaTR);
      L2.Boost(-betaTR);
      H1.Boost(-betaTR);
      H2.Boost(-betaTR);

      for (int i=0; i < otherJets.size(); i++) {
	otherJets[i].Boost(-BL);
	otherJets[i].Boost(-betaTR);
      }

    
      //at this stage you can calculate a useful angle, the azimuthal angle between 
      //the boost from the CMz to ~CM frame and the visible system of particles
      dPhiCM = (B1+B2+L1+L2).Vect().DeltaPhi(betaTR);

      //Now, we must boost to each of the respective TR frames (or top rest frames in ttbar)
      vBETA = Boost_type1(H1,H2);
      B1.Boost(-vBETA);
      B2.Boost(vBETA);
      L1.Boost(-vBETA);
      L2.Boost(vBETA);




      //Useful variable is gamma associated with that boost, which corresponds to how far off threshold the ttbar event is
      gammaR = 1./sqrt(1.-vBETA.Mag2());
	
      //In the respective top rest frames you can calculate the energy of the B's (sensitive to first mass splitting) and the helicity angle
      //energies
      EB1 = B1.E();
      EB2 = B2.E();
	
      //calculate top helicity angles
      CosThetaB1 = B1.Vect().Dot(vBETA)/(B1.P()*vBETA.Mag());
      CosThetaB2 = B2.Vect().Dot(-vBETA)/(B2.P()*vBETA.Mag());
		
      //calculate top mass - energy of all objects in T-frame assuming massless neutrinos
      TopMass1 = B1.E() + L1.E() + sqrt((B1+L1).Vect().Mag2());
      TopMass2 = B2.E() + L2.E() + sqrt((B2+L2).Vect().Mag2());

      TLorentzVector TopJet1, TopJet2;
      TopJet1.SetXYZM(0.,0.,0.,TopMass1);
      TopJet2.SetXYZM(0.,0.,0.,TopMass2);

      // Boost back to CM frame
      TopJet1.Boost(vBETA);
      TopJet2.Boost(-vBETA);
      // Now COMBINE JETS on all jets in CM frame
      otherJets.push_back(TopJet1);
      otherJets.push_back(TopJet2);
      vector<TLorentzVector> CMHems  = CombineJets(otherJets);

      if (CMHems.size()>=2){
	TLorentzVector CMHem1 = CMHems[0];
	TLorentzVector CMHem2 = CMHems[1];
	GluinoMass1 = CMHem1.M();
	GluinoMass2 = CMHem2.M();
      }
		
      //Now, we must boost to each of the respective WR frames to evaluate the lepton energy
      TVector3 vBeta1 = (-1./(L1.E()+sqrt((L1+B1).Vect().Mag2())))*B1.Vect();
      TVector3 vBeta2 = (-1./(L2.E()+sqrt((L2+B2).Vect().Mag2())))*B2.Vect();
      L1.Boost(-vBeta1);
      L2.Boost(-vBeta2);
      double mygamma1 = 1./sqrt(1.-vBeta1.Mag2());
      double mygamma2 = 1./sqrt(1.-vBeta2.Mag2());
	
      //second mass splitting
      EL1 = L1.E();
      EL2 = L2.E();

      //calculate W helicity angles
      CosThetaL1 = L1.Vect().Dot(vBeta1)/(L1.P()*vBeta1.Mag());
      CosThetaL2 = L2.Vect().Dot(vBeta2)/(L2.P()*vBeta2.Mag());
    }

    // Classic Razor Variables
    // dummy values                                          
    pTPFHem1 = -9999.;
    etaPFHem1 = -9999.;
    phiPFHem1 = -9999.;
    pTPFHem2 = -9999.;
    etaPFHem2 = -9999.;
    phiPFHem2 = -9999.;
    PFR = -99999.;
    MR = -99999.;
    MR_pTcorr = -99999.;
    
    // hemispheres
     vector<TLorentzVector> tmpJet = CombineJets(PFJet);
       
     if(tmpJet.size() >= 2) {
       TLorentzVector PFHem1 = tmpJet[0];
       TLorentzVector PFHem2 = tmpJet[1];
       // PFMET
       TVector3 MET(pxPFMet[2], pyPFMet[2], 0.);
       MRT = CalcMTR(PFHem1, PFHem2, MET);
       double variable = -999999.;
       double Rvariable = -999999.;
       variable = CalcGammaMRstar(PFHem1, PFHem2);
       if(variable >0) Rvariable = MRT/variable;

       // *NEW* pT-corrected version of MR
       MR_pTcorr = shatR(PFHem1+PFHem2,MET);
	 
       // fill the R and hem part of the output tree
       pTPFHem1 = PFHem1.Pt();
       etaPFHem1 = PFHem1.Eta();
       phiPFHem1 = PFHem1.Phi();
       pTPFHem2 = PFHem2.Pt();
       etaPFHem2 = PFHem2.Eta();
       phiPFHem2 = PFHem2.Phi();
       PFR = Rvariable;
       RSQ = Rvariable*Rvariable;
       MR = variable;
     }



     // Classic Razor Variables with FastJet Hemispheres
    // int fjetsize = 9999;
    // double conesize = 0.5;
    // vector<Jet> fJet = FastJetAlgorithmForceTwo(PFJet,conesize,2.);
    // fjetsize = (int) fJet.size();
    // if (fjetsize!=2) cout << "FAILURE: fJet.size() = "  << fjetsize << endl;
    //  vector<TLorentzVector> fastJet;
    //  fastJet.push_back(fJet[0].Get4Vector());
    //  fastJet.push_back(fJet[1].Get4Vector());

    //  // fast-hemispheres 
    //  if(fastJet.size() >= 2) {
    //    TLorentzVector FJHem1 = fastJet[0];
    //    TLorentzVector FJHem2 = fastJet[1];
    //    // PFMET
    //    TVector3 MET(pxPFMet[2], pyPFMet[2], 0.);
    //    FJMRT = CalcMTR(FJHem1, FJHem2, MET);
    //    double variable = -999999.;
    //    double Rvariable = -999999.;
    //    variable = CalcGammaMRstar(FJHem1, FJHem2);
    //    if(variable >0) Rvariable = MRT/variable;
	 
    //    // fill the R and hem part of the output tree
    //    pTFJHem1 = FJHem1.Pt();
    //    etaFJHem1 = FJHem1.Eta();
    //    phiFJHem1 = FJHem1.Phi();
    //    pTFJHem2 = FJHem2.Pt();
    //    etaFJHem2 = FJHem2.Eta();
    //    phiFJHem2 = FJHem2.Phi();
    //    FJR = Rvariable;
    //    FJRSQ = Rvariable*Rvariable;
    //    FJMR = variable;
    //  }

     //Gen-Level
     pT1 = -999;
     eta1 = -999;
     phi1 = -999;     
     pT2 = -999;
     eta2 = -999;
     phi2 = -999;
     idMc1 = -99;
     idMc2 = -99;
     if(!_isData) {
       nLepTopDecay = 0;
       int iL1 = -99;
       int iL2 = -99;
       
       double deltaRGen1;
       double deltaEGen1 = 999999999;
       double deltaRGen2;
       double deltaEGen2 = 999999999;

       for(int i=0; i<nMc; i++) {
	 // TT final state
	 if(abs(idMc[mothMc[i]]) == 24) {
	   if(abs(idMc[i]) >= 11 &&
	      abs(idMc[i]) <= 18) {
	     nLepTopDecay ++;
	   }
	 }
	 // Delta R 0.2 window of selecton dileptons
	 double tempDeltaEGen;
	 deltaRGen1 = sqrt(pow(DiLepton[0].Eta()-etaMc[i],2) + pow(DiLepton[0].Phi()-phiMc[i],2));
	 tempDeltaEGen = abs(DiLepton[0].E()-energyMc[i]);
	 if(deltaRGen1<0.2 && statusMc[i]==3) {
	   if (tempDeltaEGen < deltaEGen1) { 
	     iL1 = i; 
	     deltaEGen1 = tempDeltaEGen;
	   }
	 }
	 deltaRGen2 = sqrt(pow(DiLepton[1].Eta()-etaMc[i],2) + pow(DiLepton[1].Phi()-phiMc[i],2));
	 tempDeltaEGen = abs(DiLepton[1].E()-energyMc[i]);
	 if(deltaRGen2<0.2 && statusMc[i]==3) {
	   if (tempDeltaEGen < deltaEGen2) { 
	     iL2 = i; 
	     deltaEGen2 = tempDeltaEGen;
	   }
	 }
       }	     
       nLepTopDecay = nLepTopDecay/2;
	
       if(iL1>=0) {
	 pT1 = pMc[iL1]*sin(thetaMc[iL1]);
	 eta1 = etaMc[iL1];
	 phi1 = phiMc[iL1];
	 idMc1 = idMc[iL1];
	 idMothMc1 = idMc[mothMc[iL1]];
	 idGrandMothMc1 = idMc[mothMc[mothMc[iL1]]];
       } 
       if(iL2>=0) {
	 pT2 = pMc[iL2]*sin(thetaMc[iL2]);
	 eta2 = etaMc[iL2];
	 phi2 = phiMc[iL2];
	 idMc2 = idMc[iL2];
	 idMothMc2 = idMc[mothMc[iL2]];
	 idGrandMothMc2 = idMc[mothMc[mothMc[iL2]]];
       } 
     }


     // fill output tree
     run = runNumber;
     evNum = eventNumber;
     bx = eventNumber;
     ls = lumiBlock;
     orbit = orbitNumber;
     outTree->Fill();
  }
  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("NpassedPF_4Jet",  &NpassedPF_4Jet,  "NpassedPF_4Jet/D");
  effTree->Branch("Npassed_1b",      &Npassed_1b,      "Npassed_1b/D");
  effTree->Branch("Npassed_Mll",     &Npassed_Mll,     "Npassed_Mll/D");
  effTree->Branch("Npassed_Lept",    &Npassed_Lept,    "Npassed_Lept/D");

  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  effTree->Write();
  //  if(isSMS)FullSMSTree->Write();
  file->Close();
}

struct RazorMultiB::JetConfig{
  fastjet::JetDefinition theJetDef;
  fastjet::ActiveAreaSpec theAreaSpec;
};

vector<Jet> RazorMultiB::FastJetAlgorithmForceTwo(vector<TLorentzVector> InputCollection, double Rparam, double thePtMin){
  string JetFinder = "antikt_algorithm";
  string Strategy = "Best";
  double theDcut = -1;
  double theNjets = -1;
  string UE_subtraction = "no";
  bool theDoSubtraction = false;
  double theGhost_EtaMax = 6.0;
  int theActive_Area_Repeats = 5;
  double theGhostArea = 0.01;

  theJetConfig = new JetConfig;
  theJetConfig->theAreaSpec=fastjet::ActiveAreaSpec(theGhost_EtaMax, theActive_Area_Repeats, theGhostArea);

  fastjet::JetFinder jet_finder = fastjet::antikt_algorithm;

  fastjet::Strategy strategy = fastjet::Best;

  int theMode = 0;
  
  if((theNjets!=-1)&&(theDcut==-1)){
    theMode = 3;
  } else if((theNjets==-1)&&(theDcut!=-1)){
    theMode = 2;
  } else if((theNjets!=-1)&&(theDcut!=-1)){
    theMode = 1;
  } else {
    theMode = 0;
  }

  theJetConfig->theJetDef = fastjet::JetDefinition(jet_finder, Rparam, strategy);

  std::vector<fastjet::PseudoJet> input_vectors;
  int index_ = 0;
  for(int i = 0; i < InputCollection.size(); i++){
    double px = InputCollection[i].Px();
    double py = InputCollection[i].Py();
    double pz = InputCollection[i].Pz();
    double E = InputCollection[i].E();
    fastjet::PseudoJet PsJet(px,py,pz,E);
    PsJet.set_user_index(index_);
    input_vectors.push_back(PsJet);
    index_++;
  }

  vector<Jet> output;
  if(index_ == 0) return output;

  std::vector<fastjet::PseudoJet> theJets;

  // running without subtraction; need to add code for with subtraction
  
  fastjet::ClusterSequence clust_seq(input_vectors, theJetConfig->theJetDef);

  if((theNjets==-1)&&(theDcut==-1)){
    theJets=clust_seq.inclusive_jets(thePtMin);
  } else if((theNjets!=-1)&&(theDcut==-1)){
    theJets=clust_seq.exclusive_jets(theNjets);
  } else if((theNjets==-1)&&(theDcut!=-1)){
    theJets=clust_seq.exclusive_jets(theDcut);
  } else if((theNjets!=-1)&&(theDcut!=-1)){
    theJets=clust_seq.inclusive_jets(thePtMin);
  } else {
    theJets=clust_seq.inclusive_jets(thePtMin);
  }
  while (theJets.size()>2 and Rparam < 3.14*0.5) {
    theJetConfig->theJetDef = fastjet::JetDefinition(jet_finder, Rparam, strategy);  
    fastjet::ClusterSequence clust_seq_temp(input_vectors, theJetConfig->theJetDef);
    theJets=clust_seq.inclusive_jets(thePtMin);
    Rparam += 0.01;
    }
  //cout << "Rparam = " << Rparam << endl;
  // setting up code shamelessly stolen from 
  // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/RecoParticleFlow/PostProcessing/interface/FastJetAlgo.h?revision=1.2&view=markup

  std::vector<fastjet::PseudoJet> output_ = theJets;
  std::vector< std::vector< fastjet::PseudoJet > > ptrvector;
  for (int ct = 0; ct != (int) theJets.size(); ++ct ) {
    ptrvector.push_back(clust_seq.constituents(theJets[ct]));
  }
  int algorithm_ = 2;
  double cont =0;
  double distance_ = Rparam - 0.5;
  /////////// cycle to really force the production of two final jets, implemented for kt and anti-kt only
  while(output_.size() > 2 &&  cont < 100){
    // if (cont!=0) cout << "distance_ = " << distance_ << endl;
    double dimin = 1e09;
    if(algorithm_ == 2)
      dimin = 0;
    unsigned sel1 = 0;
    unsigned sel2 = 0;
    double px1 = 0.,py1 = 0.,pt1 = 0.;
    double phi1 = 0.,eta1 = 0.;
    double di = 0;
    for(int it = 0; it < output_.size(); it++){
      px1 = output_[it].px();
      py1 = output_[it].py();
      pt1 = sqrt( px1*px1 + py1*py1 );
      if(algorithm_ == 2){
	di = pow(1/pt1,2);
	if( di > dimin){
	  dimin = di;
	  sel1 = it;
	}
      }
    }
    double dijmin = 1e09;
    if(algorithm_ == 2)
      dijmin = 0;
    phi1 = output_[sel1].phi();
    eta1 = output_[sel1].eta();
    pt1 = sqrt( output_[sel1].px()*output_[sel1].px() + output_[sel1].py()*output_[sel1].py() );
    for(int it = 0; it < output_.size(); it++){
      if(it != sel1){
	double px,py,pz,pt;
	px = output_[it].px();
	py = output_[it].py();
	pz = output_[it].pz();
	pt = sqrt(px*px+py*py);
	double phi = output_[it].phi();
	double eta = output_[it].eta();
	double r = sqrt( (phi1-phi)*(phi1-phi) + (eta1 -eta)*(eta1-eta) );
	if(algorithm_ == 2 ){
	  if( 1./(pt1*pt1) < 1./(pt*pt) )
	    di = 1./(pt1*pt1);
	  else
	    di = 1./(pt*pt);
	}
	double dij = di*r*r/pow(distance_,2); 
	if(algorithm_ == 2){
	  if(dijmin < dij){
	    dijmin = dij;
	    sel2 = it;
	  }
	}
      }
    }
    distance_ = sqrt(dijmin*pow(distance_,2)/dimin);
    std::vector< fastjet::PseudoJet > temp_;
    std::vector< std::vector< fastjet::PseudoJet > > tempC_;
    std::vector< fastjet::PseudoJet > tempb_;
    temp_ = output_;
    output_.clear();
    tempC_ = ptrvector;
    ptrvector.clear();
    for(int it = 0; it < temp_.size(); it++){
      if( it != sel1 && it != sel2  ){
	output_.push_back(temp_[it]); 
	ptrvector.push_back(tempC_[it]);
      }
      else if(it == sel1){
	fastjet::PseudoJet psj(temp_[sel1].px() + temp_[sel2].px(),
			       temp_[sel1].py() + temp_[sel2].py(),
			       temp_[sel1].pz() + temp_[sel2].pz(),
			       temp_[sel1].e() + temp_[sel2].e() );
	unsigned it1;
	for(it1 = 0; it1 < tempC_[sel1].size() ; it1++){
	  tempb_.push_back(tempC_[sel1][it1]);
	}
	for(it1 =0; it1 < tempC_[sel2].size() ; it1++){
	  tempb_.push_back(tempC_[sel2][it1]);
	}
	ptrvector.push_back(tempb_);
	output_.push_back(psj);
      }
    }
    cont++;
    if(cont == 100) 
      cout << "If this is printed out, you just tried to create a final state of two jets using an inclusive algorithm a hundred times, please check the logic, something is going on!" << endl;
  }
  
  
  
  //here, for the reco jets, need to loop through constituents to get fractions
  for(std::vector<fastjet::PseudoJet>::const_iterator itJet=output_.begin();
      itJet!=output_.end();itJet++){
    
    double px = (*itJet).px();
    double py = (*itJet).py();
    double pz = (*itJet).pz();
    double E = (*itJet).E();
    TLorentzVector J(px,py,pz,E);
    output.push_back(Jet(J,0.0,0.0));
  }

  delete theJetConfig;
  return output;
}
  

TVector3 RazorMultiB::Boost_type1(TLorentzVector H1, TLorentzVector H2){
	
	TVector3 vBETA = (1./(H1.E()+H2.E()))*(H1.Vect()-H2.Vect());
	
	return vBETA;
}

TVector3 RazorMultiB::BetaTR(TLorentzVector TOT, TVector3 MET){
	double myshatR = shatR(TOT,MET);
	TVector3 BCM = TOT.Vect()+MET;
	BCM.SetZ(0.0);
	BCM = (1./(sqrt(4.*myshatR*myshatR+BCM.Dot(BCM))))*BCM;
	
	return BCM;
}

double RazorMultiB::shat3D(TLorentzVector TOT, TVector3 P){
	double E = TOT.E();
	double Pz = 0.0;
	
	float MR = sqrt(E*E-Pz*Pz);
	
	
	TVector3 vI = P;
	
	TVector3 vpt = TOT.Vect();
	
	float MR2 = 2.*(MR*MR-vpt.Dot(vI)+MR*sqrt(MR*MR+vI.Dot(vI)-2.*vI.Dot(vpt)));
	
	return sqrt(MR2);
	
}	

double RazorMultiB::shatR(TLorentzVector TOT, TVector3 MET){
	double E = TOT.E();
	double Pz = TOT.Pz();
	
	float MR = sqrt(E*E-Pz*Pz);
	
	
	TVector3 vI = MET+TOT.Vect();
	vI.SetZ(0.0);
	
	TVector3 vpt = TOT.Vect();
	vpt.SetZ(0.0);
	
	float MR2 = 0.5*(MR*MR-vpt.Dot(vI)+MR*sqrt(MR*MR+vI.Dot(vI)-2.*vI.Dot(vpt)));
	
	return sqrt(MR2);
	
}



std::vector<float> RazorMultiB::ParseEvent(){
	
  std::vector<std::string>::const_iterator c_begin = commentLHE->begin();
  std::vector<std::string>::const_iterator c_end = commentLHE->end();
    
  float mg, mchi;
  for( std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
    size_t found = (*cit).find("model");
    if( found != std::string::npos)   {    
      size_t foundLength = (*cit).size();
      found = (*cit).find(" ");
      std::string smaller = (*cit).substr(found+1,foundLength);
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      //cout << smaller;
	
      std::istringstream iss(smaller);
      iss >> mg;
      iss.clear();
      //cout << mg << endl;
	
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
            
      iss.str(smaller);
      iss >> mchi;
      iss.clear();
      //cout << mchi << endl;
    }
  }
    
  std::vector<float> parameterPoint;
  parameterPoint.push_back(mg);
  parameterPoint.push_back(mchi);
  return parameterPoint;
}
  
 


bool RazorMultiB::goodJetID(int i){
  ///////////////////////////////////
  //Now, we do PFNoPU jets with PU corrections 
  ///////////////////////////////////
  TLorentzVector jet;
  double px = pxAK5PFNoPUJet[i];
  double py = pyAK5PFNoPUJet[i];
  double pz = pzAK5PFNoPUJet[i];
  double E = sqrt(px*px+py*py+pz*pz);
  jet.SetPxPyPzE(px,py,pz,E);
       
  bool good_jet = false;
  double EU = uncorrEnergyAK5PFNoPUJet[i];

  double fHAD = (neutralHadronEnergyAK5PFNoPUJet[i]+chargedHadronEnergyAK5PFNoPUJet[i])/EU;

  if(fHAD > 0.99){
    good_jet = false;
  }
  else {
    int nConstituents = chargedHadronMultiplicityAK5PFNoPUJet[i]+neutralHadronMultiplicityAK5PFNoPUJet[i]+photonMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i]+HFHadronMultiplicityAK5PFNoPUJet[i]+HFEMMultiplicityAK5PFNoPUJet[i];
    int chargedMult = chargedHadronMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i];

    float photonFrac = photonEnergyAK5PFNoPUJet[i]/EU;
    float electronFrac = electronEnergyAK5PFNoPUJet[i]/EU;
    float muonFrac = muonEnergyAK5PFNoPUJet[i]/EU;
    float neutralHadFrac = neutralHadronEnergyAK5PFNoPUJet[i]/EU;
    float chargedHadFrac = chargedHadronEnergyAK5PFNoPUJet[i]/EU;
    float HFHadFrac = HFHadronEnergyAK5PFNoPUJet[i]/EU;
    float HFEMFrac = HFEMEnergyAK5PFNoPUJet[i]/EU;

    int photonMult = photonMultiplicityAK5PFNoPUJet[i];
    int electronMult = electronMultiplicityAK5PFNoPUJet[i];
    int muonMult = muonMultiplicityAK5PFNoPUJet[i];
    int neutralHadMult = neutralHadronMultiplicityAK5PFNoPUJet[i];
    int chargedHadMult = chargedHadronMultiplicityAK5PFNoPUJet[i];
    int HFHadMult = HFHadronMultiplicityAK5PFNoPUJet[i];
    int HFEMMult = HFEMMultiplicityAK5PFNoPUJet[i];

    if((neutralHadFrac < 0.99) && (photonFrac < 0.99) && (nConstituents > 1)) {
      //outside of tracker acceptance, these are the only requirementspf
      if (fabs(jet.Eta())>=2.4) good_jet = true;
      //inside of the tracker acceptance, there are extra requirements				     
      else {
	if ((chargedHadFrac > 0.0) && (chargedMult > 0) && (electronFrac < 0.99)) good_jet = true;
      }
    }
  }
  return good_jet;
}


bool RazorMultiB::FailFilters(){
  bool FAIL = false;
  /*
    This is how METFlags is defined in Vecbos:
    http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsToWW2e/src/CmsMetFiller.cc?revision=1.13&view=markup

   *(privateData_->filterBits) = 
   (eeBadScFilterFlag << 8) | (hcalLaserEventFilterFlag << 7) | (HBHENoiseFilterResultFlag << 6) | 
   (isNotDeadEcalCluster << 5) | (trackerFailureFilterFlag << 4) | (CSCHaloFilterFlag << 3) | 
   ( drDead << 2 ) | ( drBoundary << 1 ) | ECALTPFilterFlag;
  */
  //cout << METFlags << endl;

  if((METFlags >> 0)%2 == 0) FAIL = true; //ecal dead cell tp
  //if((METFlags >> 1)%2 == 0) FAIL = true; // dr boundary
  //if((METFlags >> 2)%2 == 0) FAIL = true; // dr dead
  if((METFlags >> 3)%2 == 0) FAIL = true; // csc hal
  if((METFlags >> 4)%2 == 0) FAIL = true; // tracker failure
  //if((METFlags >> 5)%2 == 0) FAIL = true; // ecal dead cluster
  if((METFlags >> 6)%2 == 0) FAIL = true; // hbhe noise
  if((METFlags >> 7)%2 == 0) FAIL = true; // hcal laser
  if((METFlags >> 8)%2 == 0) FAIL = true; // bad ee sc
  //if((METFlags >> 9)%2 == 0) FAIL = true; //ecal laser
   
  return FAIL;
}


struct RazorMultiB::EventIndex {
	int RunNumber;
	Long64_t EventNumber;
	
	EventIndex() : RunNumber(0), EventNumber(0) {}
	
	bool operator <(const EventIndex &other) const
	{
		if(RunNumber < other.RunNumber)
			return true;
		if(RunNumber > other.RunNumber)
			return false;
		
		if(EventNumber < other.EventNumber)
			return true;
		if(EventNumber > other.EventNumber)
			return false;
		
		return false;
	}
};


void RazorMultiB::InitEventFlag(char *s_Event){
  ifstream inputfile(s_Event);
  
  int RUN_NUMBER;
  int LS_NUMBER;
  Long64_t EVENT_NUMBER;
  
  EventIndex index;
  
  cout << "Reading bad event list" << endl;
	
  while(!inputfile.eof()){
    inputfile >> RUN_NUMBER >> LS_NUMBER >> EVENT_NUMBER;
    
    index.RunNumber = RUN_NUMBER;
    index.EventNumber = EVENT_NUMBER;
    
    if(index.RunNumber < 0 || index.EventNumber < 0)
      continue;
    
    //Is this event/run-number combo alreading in the map?
    if(EventCounts.find(index) == EventCounts.end()){ //no
      EventCounts.insert(pair<EventIndex, int>(index, 1));
    } else { //yes
      EventCounts[index] = EventCounts[index] + 1;
    }
    index.RunNumber = -1;
    index.EventNumber = -1;
  }
}


bool RazorMultiB::isFlagged(){
  EventIndex index;
  index.EventNumber = eventNumber;
  index.RunNumber = runNumber;
  
  if(EventCounts.find(index) == EventCounts.end()){ //yes
    //cout << "Event not in list" << endl;
    return false;
  } else {
    //cout << "Event IS in list" << endl;
    return true;
    if(EventCounts[index] == 1){
      //cout << "Found the event - all is good in the world" << endl;
      return true;
    }
  } 
}



vector<TLorentzVector> RazorMultiB::CombineJetsTs(vector<TLorentzVector> myjets, vector<TLorentzVector> Ts){
  
  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;
  bool foundGood = false;
  
  int N_comb = 1;
  for(int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }
  
  double M_min = 9999999999.0;
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
	j_temp1 += myjets[count];
      } else {
	j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }    

    

    double M_temp = std::min((j_temp1 + Ts[0]).M2()+(j_temp2 + Ts[1]).M2() , (j_temp1 + Ts[1]).M2()+(j_temp2 + Ts[0]).M2() ) ;
    // smallest mass
    if (M_temp < M_min){
      M_min = M_temp;
      if ((j_temp1 + Ts[0]).M2()+(j_temp2 + Ts[1]).M2() == M_temp ){
	j1 = j_temp1 + Ts[0];
	j2 = j_temp2 + Ts[1];  
      }
      else{ 
      j1 = j_temp1 + Ts[1];
      j2 = j_temp2 + Ts[0];
      }
    }
  }

  // set masses to 0
  // j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  // j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;  
}



