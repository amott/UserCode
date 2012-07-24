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
#include "RazorDiMuB.hh"

RazorDiMuB::RazorDiMuB(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

RazorDiMuB::RazorDiMuB(TTree *tree, string jsonFile,bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
  }
}

RazorDiMuB::~RazorDiMuB() {}

void RazorDiMuB::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorDiMuB::SetWeight(double weight){
  _weight=weight;
}

void RazorDiMuB::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  int HLT_DoubleMu;
  int HLT_DoubleEle;
  int HLT_MuEle;

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
  double W=_weight;
  double mst, mchi;
  int    BOX_NUM;
  int    ss;
  int    nPV;

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
  outTree->Branch("W", &W, "W/D");
  outTree->Branch("mst", &mst, "mst/D");
  outTree->Branch("mchi", &mchi, "mchi/D");
  outTree->Branch("nPV", &nPV, "nPV/I");

  double Npassed_In = 0;
  double Npassed_PV = 0;
  //Jets
  double NpassedPF_2Jet = 0;
  //Leptons
  double Npassed_Lept=0;
  //B-tag
  double Npassed_1b=0;

  double weightII = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  std::vector<std::string> maskHLT_DoubleMu; 
  maskHLT_DoubleMu.push_back("HLT_Mu13_Mu8_v");
  maskHLT_DoubleMu.push_back("HLT_Mu17_Mu8_v");
  
  std::vector<std::string> maskHLT_DoubleEle; 
  maskHLT_DoubleEle.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  
  std::vector<std::string> maskHLT_MuEle; 
  maskHLT_MuEle.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  maskHLT_MuEle.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  
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
    
    Npassed_In += weightII;
    
    double m0=9999, m12=9999, mc=9999;
     
    // to integrate with others
    /*
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
    */
    mst=m0;
    mchi=mc;
    
    //HLT
    int passedHLT = HLT_DoubleMu+HLT_DoubleEle+HLT_MuEle;
    if(passedHLT==0 and _isData==true) continue;
    
    // find highest-pT PV [replace with Sagar's code]
    int iPV = passPV();
    if(iPV<0) continue;
    Npassed_PV += weightII;
    nPV = N_PV_EVENT;
    
    vector<TLorentzVector> PFJet;
    vector <int> iPFJet;
    bool badjet = false;

    for(int i=0; i< nAK5PFNoPUJet; i++) {
      double EU = (chargedHadronEnergyAK5PFNoPUJet[i]+neutralHadronEnergyAK5PFNoPUJet[i]+photonEnergyAK5PFNoPUJet[i]+electronEnergyAK5PFNoPUJet[i]+muonEnergyAK5PFNoPUJet[i]);
      TLorentzVector myJet(pxAK5PFNoPUJet[i], pyAK5PFNoPUJet[i], pzAK5PFNoPUJet[i], energyAK5PFNoPUJet[i]);
      // Apply jet correction first 
      bool good_jet = false;
      if (myJet.Pt() > 30.0 && fabs(myJet.Eta()) < 3.0) {
        double fHAD = (neutralHadronEnergyAK5PFNoPUJet[i]+chargedHadronEnergyAK5PFNoPUJet[i])/EU;
        if(fHAD > 0.99) {
          badjet = true;
          break;
        }
	if (fabs(myJet.Eta()) < 2.4 && chargedEmEnergyAK5PFNoPUJet[i] > 0.0 ) good_jet = true;
        if (fabs(myJet.Eta()) < 2.4 && chargedHadronEnergyAK5PFNoPUJet[i] > 0.0 ) good_jet = true;
        if (muonEnergyAK5PFNoPUJet[i] > 0.0 ) good_jet = true;
        if (fabs(myJet.Eta()) >= 2.4 && neutralHadronEnergyAK5PFNoPUJet[i]/EU < 0.99 && neutralEmEnergyAK5PFNoPUJet[i]/EU < 0.99) good_jet = true;
        if (good_jet==false) {
          badjet = true;
          break;
        }
      }

      if (myJet.Pt()>30. && fabs(myJet.Eta())< 3.0) {
	PFJet.push_back(myJet);
	iPFJet.push_back(i);
      }
    } 

    // jet ID
    if (badjet == true) continue;
    //2Jets
    if(int(PFJet.size())<2) continue;
    NpassedPF_2Jet+=weightII;
    
    //Create arrays with the b-discriminators of the jets
    nBtag = 0;
    for(int b=0; b< iPFJet.size(); b++){
      int n=iPFJet.at(b);
      if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag++; 
    }
    
    Npassed_1b+=weightII;
    
    // Boxes
    vector<int> iMuLoose;
    for(int i=0; i<nMuon; i++) {
      TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
      if(isLooseMuon(i) and thisMu.Pt() > 20.) {
	iMuLoose.push_back(i);
      }	  
    }

    vector<int> iEleTight;
    for(int i=0; i<nEle; i++) {
      TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
      if(isTrigElectron(i) && thisEle.Pt() > 20.) {
	iEleTight.push_back(i);
      }
    }

    BOX_NUM = -99;
    if(iMuLoose.size() >0 && iEleTight.size()>0){
      BOX_NUM = 0;
      ss = chargeMuon[iMuLoose[0]]*chargeEle[iEleTight[0]];
    } else if(iMuLoose.size()>1) {
      ss = chargeMuon[iMuLoose[0]]*chargeMuon[iMuLoose[1]];
      BOX_NUM = 1;
    } else if(iEleTight.size()>1) {
      ss = chargeEle[iEleTight[0]]*chargeEle[iEleTight[1]];
      BOX_NUM = 2;
    }
    // two good leptons
    if(BOX_NUM<0) continue;
    Npassed_Lept += weightII;

    // dummy values                                          
    pTPFHem1 = -9999.;
    etaPFHem1 = -9999.;
    phiPFHem1 = -9999.;
    pTPFHem2 = -9999.;
    etaPFHem2 = -9999.;
    phiPFHem2 = -9999.;
    PFR = -99999.;
    RSQ = -99999.;
    MR = -99999.;
    
    // hemispheres
     vector<TLorentzVector> tmpJet = CombineJets(PFJet);
     if(tmpJet.size() >= 2) {
       TLorentzVector PFHem1 = tmpJet[0];
       TLorentzVector PFHem2 = tmpJet[1];

       // PFMET
       TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);
       MRT = CalcMTR(PFHem1, PFHem2, MET);
       double variable = -999999.;
       double Rvariable = -999999.;
       variable = CalcGammaMRstar(PFHem1, PFHem2);
       if(variable >0) Rvariable = MRT/variable;
	 
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
  effTree->Branch("NpassedPF_2Jet",  &NpassedPF_2Jet,  "NpassedPF_2Jet/D");

  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  effTree->Write();
  //  if(isSMS)FullSMSTree->Write();
  file->Close();
}
  

