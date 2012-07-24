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
#include "MonoJet.hh"

MonoJet::MonoJet(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

MonoJet::MonoJet(TTree *tree, string jsonFile,bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
  }
}

MonoJet::~MonoJet() {}

void MonoJet::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void MonoJet::SetWeight(double weight){
  _weight=weight;
}

void MonoJet::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  int HLT_CentralJet80_MET80;
  int HLT_CentralJet100_MET80;

  // PF block
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double pTJet1;
  double etaJet1;
  double phiJet1;
  double pTJet2;
  double etaJet2;
  double phiJet2;
  double pTMET;
  double phiMET;
  double pTMETnoMu;
  double phiMETnoMu;
  double MR;
  double RSQ;
  double pTMu1;
  double etaMu1;
  double phiMu1;
  double pTMu2;
  double etaMu2;
  double phiMu2;
  double pTEle1;
  double etaEle1;
  double phiEle1;
  double pTEle2;
  double etaEle2;
  double phiEle2;
  int nBtag;
  int nJet;
  int nElectron;
  int nMu;
  double weight=_weight;
  double mg, mst, mchi;

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");

  // HLT bits
  outTree->Branch("HLT_CentralJet80_MET80", &HLT_CentralJet80_MET80, "HLT_CentralJet80_MET80/I"); 
  outTree->Branch("HLT_CentralJet100_MET80", &HLT_CentralJet100_MET80, "HLT_CentralJet100_MET80/I"); 

  // PF block
  outTree->Branch("pTJet1", &pTJet1, "pTJet1/D");
  outTree->Branch("etaJet1", &etaJet1, "etaJet1/D");
  outTree->Branch("phiJet1", &phiJet1, "phiJet1/D");
  outTree->Branch("pTJet2", &pTJet2, "pTJet2/D");
  outTree->Branch("etaJet2", &etaJet2, "etaJet2/D");
  outTree->Branch("phiJet2", &phiJet2, "phiJet2/D");

  outTree->Branch("pTMu1", &pTMu1, "pTMu1/D");
  outTree->Branch("etaMu1", &etaMu1, "etaMu1/D");
  outTree->Branch("phiMu1", &phiMu1, "phiMu1/D");
  outTree->Branch("pTMu2", &pTMu2, "pTMu2/D");
  outTree->Branch("etaMu2", &etaMu2, "etaMu2/D");
  outTree->Branch("phiMu2", &phiMu2, "phiMu2/D");

  outTree->Branch("pTEle1", &pTEle1, "pTEle1/D");
  outTree->Branch("etaEle1", &etaEle1, "etaEle1/D");
  outTree->Branch("phiEle1", &phiEle1, "phiEle1/D");
  outTree->Branch("pTEle2", &pTEle2, "pTEle2/D");
  outTree->Branch("etaEle2", &etaEle2, "etaEle2/D");
  outTree->Branch("phiEle2", &phiEle2, "phiEle2/D");

  outTree->Branch("pTMET", &pTMET, "pTMET/D");
  outTree->Branch("phiMET", &phiMET, "phiMET/D");
  outTree->Branch("pTMETnoMu", &pTMETnoMu, "pTMETnoMu/D");
  outTree->Branch("phiMETnoMu", &phiMETnoMu, "phiMETnoMu/D");

  outTree->Branch("nBtag", &nBtag, "nBtag/I");
  outTree->Branch("nJet", &nJet, "nJet/I");
  outTree->Branch("nMu",  &nMu, "nMu/I");
  outTree->Branch("nElectron",  &nElectron, "nElectron/I");

  outTree->Branch("RSQ", &RSQ, "RSQ/D");
  outTree->Branch("MR", &MR, "MR/D");

  outTree->Branch("weight", &weight, "weight/D");
  outTree->Branch("mg", &mg, "mg/D");
  outTree->Branch("mst", &mst, "mst/D");
  outTree->Branch("mchi", &mchi, "mchi/D");

  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  std::vector<std::string> maskHLT_CentralJet80_MET80; 
  maskHLT_CentralJet80_MET80.push_back("HLT_CentralJet80_MET80");
  std::vector<std::string> maskHLT_CentralJet100_MET80; 
  maskHLT_CentralJet100_MET80.push_back("HLT_CentralJet100_MET80");
  
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
      setRequiredTriggers(maskHLT_CentralJet80_MET80); reloadTriggerMask(true); HLT_CentralJet80_MET80 = hasPassedHLT();
      setRequiredTriggers(maskHLT_CentralJet100_MET80); reloadTriggerMask(true); HLT_CentralJet100_MET80 = hasPassedHLT();
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
    
    double m0=9999, m12=9999, mc=9999;
     
      // to integrate with others
    if(!_isData && _isSMS){
      //find the simplified model parameters for T1tttt                                                                                              
      std::vector<std::string>::const_iterator c_begin = commentLHE->begin();
      std::vector<std::string>::const_iterator c_end = commentLHE->end();
      for(std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
	size_t found = (*cit).find("T1bbbb");
	if( found != std::string::npos) {
	  size_t foundLength = (*cit).size();
	  found = (*cit).find("=");
	  std::string smaller = (*cit).substr(found+1,foundLength);
	  found = smaller.find("_");
	  smaller = smaller.substr(found+1,smaller.size());
	  
	  std::istringstream iss(smaller);
	  iss >> m0;
	  iss.clear();
	  
	  found = smaller.find("_");
	  smaller = smaller.substr(found+1,smaller.size());
	  iss.str(smaller);
	  iss >> m12;
	  iss.clear();
	  
	  found = smaller.find("_");
	  smaller = smaller.substr(found+1,smaller.size());
	  iss.str(smaller);
	  iss >> mc;
	  iss.clear();
	  
	}
      }
    }
    
    mg=m12;
    mst=m0;
    mchi=mc;
    
    // event filter
    if(_isData && (METFlags >> 0) % 2 == 0) continue;
    //drBoundary
    if(_isData && (METFlags >> 1) % 2 == 0) continue;
    // drDead
    if(_isData && (METFlags >> 2) % 2 == 0) continue;
    // CSCHaloFilterFlag
    if(_isData && (METFlags >> 3) % 2 == 0) continue;
    // trackerFailureFilterFlag
    if(_isData && (METFlags >> 4) % 2 == 0) continue;
    // BE ECAL flag
    if(_isData && (METFlags >> 5) % 2 == 0) continue;

    // find highest-pT PV
    int iHighestPtVTX = -99;
    double HighestPtVTX = -99999.;
    if(nPV<1) continue;
    
    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPtVTX) iHighestPtVTX = i;
    // PV selection
    if(ndofPV[iHighestPtVTX] < 3) continue;
    if(PVzPV[iHighestPtVTX] > 25.) continue; 
    
    vector<TLorentzVector> PFJet;
    vector <int> iPFJet;
    bool badjet = false;
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      double EU = uncorrEnergyAK5PFPUcorrJet[i];
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);
      // Apply jet correction first 
      bool good_jet = false;
      if (myJet.Pt() > 30.0 && fabs(myJet.Eta()) < 3.0) {
        double fHAD = (neutralHadronEnergyAK5PFPUcorrJet[i]+chargedHadronEnergyAK5PFPUcorrJet[i])/EU;
        if(fHAD > 0.99) {
          badjet = true;
          break;
        }
	if (fabs(myJet.Eta()) < 2.4 && chargedEmEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
        if (fabs(myJet.Eta()) < 2.4 && chargedHadronEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
        if (muonEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
        if (fabs(myJet.Eta()) >= 2.4 && neutralHadronEnergyAK5PFPUcorrJet[i]/EU < 0.99 && neutralEmEnergyAK5PFPUcorrJet[i]/EU < 0.99) good_jet = true;
        if (good_jet==false) {
          badjet = true;
          break;
        }
      }

      if (myJet.Pt()>30. && fabs(myJet.Eta())< 2.4) {
	PFJet.push_back(myJet);
	iPFJet.push_back(i);
      }
    }
    // jet ID
    if (badjet == true) continue;
    //1Jet pT>120
    if(int(PFJet.size())<1) continue;

    // dummy values                                          
    pTJet1 = -9999.;
    etaJet1 = -9999.;
    phiJet1 = -9999.;
    pTJet2 = -9999.;
    etaJet2 = -9999.;
    phiJet2 = -9999.;

    int iJet1 = HighestPt(PFJet, -99);
    int iJet2 = HighestPt(PFJet, iJet1);
    if(iJet1 >=0) {
      pTJet1 = PFJet[iJet1].Pt();
      etaJet1 = PFJet[iJet1].Eta();
      phiJet1 = PFJet[iJet1].Phi();
    }
    if(pTJet1<120.) continue;
    if(iJet2 >=0) {
      pTJet2 = PFJet[iJet2].Pt();
      etaJet2 = PFJet[iJet2].Eta();
      phiJet2 = PFJet[iJet2].Phi();
    }
    nJet = PFJet.size();

    // BTAG
    nBtag = 0;
    for(int b=0; b< iPFJet.size(); b++){
      int n=iPFJet.at(b);
      if(combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[n] > 1.74) nBtag++;
    }
    
    // Muons
    vector<int> iMuTight;
    vector<TLorentzVector> MuTight;
    for(int i=0; i<nMuon; i++) {
      TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
      if(muonPassTight(i) && thisMu.Pt() > 5.) {
	iMuTight.push_back(i);
	MuTight.push_back(thisMu);
      }
    }
    // dummy values
    pTMu1 = -9999.;
    etaMu1 = -9999.;
    phiMu1 = -9999.;
    pTMu2 = -9999.;
    etaMu2 = -9999.;
    phiMu2 = -9999.;

    int iMu1 = HighestPt(MuTight, -99);
    int iMu2 = HighestPt(MuTight, iMu1);
    if(iMu1 >=0) {
      pTMu1 = MuTight[iMu1].Pt();
      etaMu1 = MuTight[iMu1].Eta();
      phiMu1 = MuTight[iMu1].Phi();
    }
    if(iMu2 >=0) {
      pTMu2 = MuTight[iMu2].Pt();
      etaMu2 = MuTight[iMu2].Eta();
      phiMu2 = MuTight[iMu2].Phi();
    }
    nMu = iMuTight.size();

    // Electrons
    vector<int> iEleTight;
    vector<TLorentzVector> EleTight;
    for(int i=0; i<nEle; i++) {
      TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
      if(electronPassWP80(i) && thisEle.Pt() > 5.) {
	iEleTight.push_back(i);
	EleTight.push_back(thisEle);
      }
    }
    // dummy values
    pTEle1 = -9999.;
    etaEle1 = -9999.;
    phiEle1 = -9999.;
    pTEle2 = -9999.;
    etaEle2 = -9999.;
    phiEle2 = -9999.;

    int iEle1 = HighestPt(EleTight, -99);
    int iEle2 = HighestPt(EleTight, iEle1);
    if(iEle1 >=0) {
      pTEle1 = EleTight[iEle1].Pt();
      etaEle1 = EleTight[iEle1].Eta();
      phiEle1 = EleTight[iEle1].Phi();
    }
    if(iEle2 >=0) {
      pTEle2 = EleTight[iEle2].Pt();
      etaEle2 = EleTight[iEle2].Eta();
      phiEle2 = EleTight[iEle2].Phi();
    }
    nElectron = iEleTight.size();

    // MET
    TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);
    pTMET = MET.Pt();
    phiMET = MET.Phi();

    // MET without muons
    if(iMu1>=0) MET = MET + MuTight[iMu1].Vect();
    if(iMu2>=0) MET = MET + MuTight[iMu2].Vect();
    pTMETnoMu = MET.Pt();
    phiMETnoMu = MET.Phi();

    // Razor (using METnoMu and the two highest-pT jets)
    RSQ = -99999.;
    MR = -99999.;
    
    if(iJet2>=0) {
       // compute boost
       double num = PFJet[iJet1].P()-PFJet[iJet2].P();
       double den = PFJet[iJet1].Pz()-PFJet[iJet2].Pz();
       // PFMET
       double MRT = CalcMTR(PFJet[iJet1], PFJet[iJet2], MET);
       double variable = -999999.;
       double Rvariable = -999999.;
       variable = CalcGammaMRstar(PFJet[iJet1], PFJet[iJet2]);
       if(variable >0) Rvariable = MRT/variable;
	 
       // fill the R and hem part of the output tree
       RSQ = Rvariable*Rvariable;
       MR = variable;
     }

     // fill output tree
     run = runNumber;
     evNum = eventNumber;
     bx = eventNumber;
     ls = lumiBlock;
     orbit = orbitNumber;
     outTree->Fill();
  }
  
  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  //  if(isSMS)FullSMSTree->Write();
  file->Close();
}
  
int MonoJet::HighestPt(vector<TLorentzVector> v, int ignore) {
  double pT = 0.;
  int iBest = -99;
  for(int i=0; i<int(v.size()); i++) {
    if(v[i].Pt()>pT && i != ignore) {
      pT = v[i].Pt();
      iBest = i;
    }
  }
  return iBest;
}
