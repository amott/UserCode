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
#include "RazorBoostedTop.hh"

RazorBoostedTop::RazorBoostedTop(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
}

RazorBoostedTop::RazorBoostedTop(TTree *tree, string json, bool goodRunLS, bool isData,int mod) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;
  MODEL=mod;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

RazorBoostedTop::~RazorBoostedTop() {}

void RazorBoostedTop::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorBoostedTop::SetWeight(double weight) {
  _weight = weight;
}

void RazorBoostedTop::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  // Jet triggers
  int HLT_Jet240;
  int HLT_Jet300;
  int HLT_Jet370;

  // PF  block
  int    passedPF;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double PFR;
  double PFRsq;
  double PFMR;
  double PFHT;
  double PFMET;

  // general event info
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  int nPV;
  double W;

  //PFJets Info for Veto
  double CenJetPt[6];
  //Lepton Info For Veto
  double LeadTightMuonPt;
  double LeadTightElePt;
  double LeadTightTauPt;

  // Top variables are global

  // prepare the output tree
  TTree* outTree[3]; 
  outTree[0] = new TTree("outTree0BT", "outTree0BT");
  outTree[1] = new TTree("outTree1BT", "outTree1BT");
  outTree[2] = new TTree("outTree2BT", "outTree2BT");

  for(int i=0; i<3; i++) {

    outTree[i]->Branch("HLT_Jet240", &HLT_Jet240, "HLT_Jet240/I");
    outTree[i]->Branch("HLT_Jet300", &HLT_Jet300, "HLT_Jet300/I");
    outTree[i]->Branch("HLT_Jet370", &HLT_Jet370, "HLT_Jet370/I");

    outTree[i]->Branch("run", &run, "run/D");
    outTree[i]->Branch("evNum", &evNum, "evNum/D");
    outTree[i]->Branch("bx", &bx, "bx/D");
    outTree[i]->Branch("ls", &ls, "ls/D");
    outTree[i]->Branch("orbit", &orbit, "orbit/D");
    outTree[i]->Branch("nPV", &nPV, "nPV/I");
    outTree[i]->Branch("W", &W, "W/D");
    
    // PF block
    outTree[i]->Branch("passedPF", &passedPF, "passedPF/I");
    outTree[i]->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
    outTree[i]->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
    outTree[i]->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
    outTree[i]->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
    outTree[i]->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
    outTree[i]->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
    outTree[i]->Branch("PFR", &PFR, "PFR/D");
    outTree[i]->Branch("PFRsq", &PFRsq, "PFRsq/D");
    outTree[i]->Branch("PFMR", &PFMR, "PFMR/D");

    outTree[i]->Branch("PFHT", &PFHT, "PFHT/D");
    outTree[i]->Branch("PFMET", &PFMET, "PFMET/D");
    
    
    
    // first Top
    outTree[i]->Branch("PFT1Pt",   &PFT1Pt, "PFT1Pt/D");
    outTree[i]->Branch("PFT1Eta",  &PFT1Eta, "PFT1Eta/D");
    outTree[i]->Branch("PFT1Phi",  &PFT1Phi, "PFT1Phi/D");
    outTree[i]->Branch("PFT1Mass", &PFT1Mass, "PFT1Mass/D");
    outTree[i]->Branch("MergedT1", &MergedT1, "MergedT1/I");

    // second Top
    outTree[i]->Branch("PFT2Pt",   &PFT2Pt, "PFT2Pt/D");
    outTree[i]->Branch("PFT2Eta",  &PFT2Eta, "PFT2Eta/D");
    outTree[i]->Branch("PFT2Phi",  &PFT2Phi, "PFT2Phi/D");
    outTree[i]->Branch("PFT2Mass", &PFT2Mass, "PFT2Mass/D");
    outTree[i]->Branch("MergedT2", &MergedT2, "MergedT2/I");

    outTree[i]->Branch("CenJetPt", CenJetPt, "CenJetPt[6]/D");

    outTree[i]->Branch("LeadTightMuonPt",&LeadTightMuonPt,"LeadTightMuonPt/D");
    outTree[i]->Branch("LeadTightElePt",&LeadTightElePt,"LeadTightElePt/D");
    outTree[i]->Branch("LeadTightTauPt",&LeadTightTauPt,"LeadTightTauPt/D");

    outTree[i]->Branch("m0", &m0, "m0/D");    
    outTree[i]->Branch("m12", &m12, "m12/D");    
    outTree[i]->Branch("tanb", &tanb, "tanb/D");    
    outTree[i]->Branch("A0", &A0, "A0/D");    
    outTree[i]->Branch("mu", &mu, "mu/D");    
  }

  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // hadronic razor triggers
  std::vector<std::string> maskHLT_Jet240; maskHLT_Jet240.push_back("HLT_Jet240_v");
  std::vector<std::string> maskHLT_Jet300; maskHLT_Jet300.push_back("HLT_Jet300_v");
  std::vector<std::string> maskHLT_Jet370; maskHLT_Jet370.push_back("HLT_Jet370_v");

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
      
      // hadronic razor triggers
      setRequiredTriggers(maskHLT_Jet240); reloadTriggerMask(true); HLT_Jet240 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Jet300); reloadTriggerMask(true); HLT_Jet300 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Jet370); reloadTriggerMask(true); HLT_Jet370 = hasPassedHLT();
      
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
    //if(nPV<1) continue;
    //for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    //if(ndofPV[iHighestPt] < 3) continue;
    //if(PVzPV[iHighestPt] > 25.) continue; 

    // HCAL FLAGS
    if(_isData && !eventPassHcalFilter()) continue;
    //ECALTPFilterFlag 
    if(_isData && METFlags << 0 == 0) continue;
    //drBoundary 
    if(_isData && METFlags << 1 == 0) continue;
    // drDead
    if(_isData && METFlags << 2 == 0) continue;
    // CSCHaloFilterFlag
    if(_isData && METFlags << 3 == 0) continue;
    // trackerFailureFilterFlag
    if(_isData && METFlags << 4 == 0) continue;
    // BE ECAL flag 
    if(_isData && METFlags << 5 == 0) continue;

    // Jet selection 
    bool goodPFevent = true;
    vector<TLorentzVector> PFPUcorrJet; 
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      // to avoid messages of 0 pT
      if(sqrt(pow(pxAK5PFPUcorrJet[i],2.)+pow(pyAK5PFPUcorrJet[i],2.))<10.) continue;
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);
      if(myJet.Pt()>40. && fabs(etaAK5PFPUcorrJet[i])< 2.4) {
	PFPUcorrJet.push_back(myJet);
      }
      // check if the event is good (from PFJets ID point of view)                                                              
      double EU = uncorrEnergyAK5PFPUcorrJet[i];
      // Apply jet correction first                                                                                             
      bool good_jet = false;
      PFHT=0;
      if (myJet.Pt() > 40.0 && fabs(etaAK5PFPUcorrJet[i]) < 3.0) {
	double fHAD = (neutralHadronEnergyAK5PFPUcorrJet[i]+chargedHadronEnergyAK5PFPUcorrJet[i])/EU;
	if(fHAD > 0.99) {
	  goodPFevent = false;
	  break;
	}
	if (neutralEmEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (fabs(etaAK5PFPUcorrJet[i])  < 2.4 && chargedEmEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (fabs(etaAK5PFPUcorrJet[i])  < 2.4 && chargedHadronEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (muonEnergyAK5PFPUcorrJet[i] > 0.0 ) good_jet = true;
	if (good_jet==false) {
	  goodPFevent = false;
	  break;
	}
	PFHT+=myJet.Pt();
      }
    }

    if(goodPFevent == false) continue;
    // at least two jets
    if(PFPUcorrJet.size()<2) continue;

    // find the 6 leading pT jets in the event
    for(int i=0;i<6;i++) CenJetPt[i]=0;

    int nJets=0;
    for(int iJet=0;iJet<nAK5PFPUcorrJet;iJet++){
      // check if the event is good (from PFJets ID point of view)                                                              
      double EU = uncorrEnergyAK5PFPUcorrJet[iJet];
      // Apply jet correction first                                                                                             
      bool good_jet = false;
      if (fabs(etaAK5PFPUcorrJet[iJet]) >= 3.0)  continue; //central jets only
      double fHAD = (neutralHadronEnergyAK5PFPUcorrJet[iJet]+chargedHadronEnergyAK5PFPUcorrJet[iJet])/EU;
      if(fHAD > 0.99) continue; //reject HCal noise

      if (neutralEmEnergyAK5PFPUcorrJet[iJet] > 0.0 ) good_jet = true;
      if (fabs(etaAK5PFPUcorrJet[iJet])  < 2.4 && chargedEmEnergyAK5PFPUcorrJet[iJet] > 0.0 ) good_jet = true;
      if (fabs(etaAK5PFPUcorrJet[iJet])  < 2.4 && chargedHadronEnergyAK5PFPUcorrJet[iJet] > 0.0 ) good_jet = true;
      if (muonEnergyAK5PFPUcorrJet[iJet] > 0.0 ) good_jet = true;
      if (good_jet==false) continue; //need one of the ID parameters to be true
      CenJetPt[nJets++] = sqrt(pow(pxAK5PFPUcorrJet[iJet],2.)+pow(pyAK5PFPUcorrJet[iJet],2.)); //fill the jet pT and increment counter
      if(nJets>=6) break;
    }

    LeadTightMuonPt=0;
    for(int i=0;i<nMuon;i++){
      if(muonPassTight(i)){
	double thisPt = sqrt(pxMuon[i]*pxMuon[i]+pyMuon[i]*pyMuon[i]);
	if(thisPt>LeadTightMuonPt) LeadTightMuonPt=thisPt;
      }
    }

    LeadTightElePt=0;
    for(int i=0;i<nEle;i++){
      if(electronPassWP80(i)){
	double thisPt = sqrt(pxEle[i]*pxEle[i]+pyEle[i]*pyEle[i]);	
	if(thisPt>LeadTightElePt) LeadTightElePt=thisPt;
      }
    }

    LeadTightTauPt=0;
    // Find the best HPS Tau
    int ntauHPS = 0;
    TLorentzVector myTauHPS(0.1, 0., 0., 0.1);
    int iTauHPS = -99;
    for (int i=0; i<nPFTau; i++) {
      if (abs(thehpsTauDiscrByTightElectronRejectionPFTau[i]-1.)<0.001 && abs(thehpsTauDiscrByTightMuonRejectionPFTau[i]-1.)<0.0001 &&
     	  abs(thehpsTauDiscrByLooseIsolationPFTau[i]-1.)<0.0001) {
     	double thisPt = sqrt(pxPFTau[i]*pxPFTau[i]+pyPFTau[i]*pyPFTau[i]);
	if(thisPt>LeadTightTauPt) LeadTightTauPt = thisPt;
      }
    }


    // use PFMET
    TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);

    // dummy values
    passedPF = 0;
    pTPFHem1 = -9999;
    etaPFHem1 = -9999;
    phiPFHem1 = -9999;
    pTPFHem2 = -9999;
    etaPFHem2 = -9999;
    phiPFHem2 = -9999;
    PFR = -99999.;
    PFRsq = -99999.;
    PFMR = -99999.;

    // hemispheres and Razor selection
    double RsqMin = 0.0; // to set after first plots
    double mRmin = 600.0;  // to set after first plots
    vector<TLorentzVector> tmpJet = CombineJets(PFPUcorrJet);
    if(tmpJet.size() <2) continue;
    
    TLorentzVector PFHem1 = tmpJet[0];
    TLorentzVector PFHem2 = tmpJet[1];

    // compute boost
    double num = PFHem1.P()-PFHem2.P();
    double den = PFHem1.Pz()-PFHem2.Pz();      
    
    double MT = CalcMTR(PFHem1, PFHem2, MET);
    double variable = -999999.;
    double Rvariable = -999999.;
    variable = CalcGammaMRstar(PFHem1, PFHem2);
    if(variable >0) Rvariable = MT/variable;

    // fill the tree
    passedPF = 1;
    pTPFHem1 = PFHem1.Pt();
    etaPFHem1 = PFHem1.Eta();
    phiPFHem1 = PFHem1.Phi();
    pTPFHem2 = PFHem2.Pt();
    etaPFHem2 = PFHem2.Eta();
    phiPFHem2 = PFHem2.Phi();
    PFR = Rvariable;
    PFRsq = Rvariable*Rvariable;
    PFMR = variable;    

    PFMET = MET.Pt();

    // baseline Razor selection
    if(PFMR<mRmin) continue;
    if(PFRsq<RsqMin) continue;

    ////////////////////////////////////////////////////////////
    // look for the merged tops
    ////////////////////////////////////////////////////////////
    int iTJet[2];

    // look for the first largest-mass jet
    iTJet[0] = -99;
    double massOne = -9999999999.;
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(PFPUcorrJet[i].M() > massOne) {
	massOne = PFPUcorrJet[i].M();
	iTJet[0] = i;
      }
    }

    // look for the second largest-mass jet
    iTJet[1] = -99;
    double massTwo = -9999999999.; 
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(i == iTJet[0]) continue;
      if(PFPUcorrJet[i].M() > massTwo) {
	massTwo = PFPUcorrJet[i].M();
	iTJet[1] = i;
      }
    }

    SetFirstTop(PFPUcorrJet[iTJet[0]], true);
    SetSecondTop(PFPUcorrJet[iTJet[1]], true);

    // fill the boxes
    // 0: No  Heavy Jets  
    // 1: One Heavy Jet  
    // 2: Two Heavy Jets  
    double pTMinTop = 140.;

    /*
    if(PFT1Mass > pTMinTop && PFT2Mass > pTMinTop)  outTree[0]->Fill();
    else if(PFT1Mass > pTMinTop) {
      MergedT2 = false;
      outTree[1]->Fill();
    } else if(PFT2Mass > pTMinTop) {
      MergedT1 = false;
      outTree[1]->Fill();
    } else {
      MergedT2 = false;
      MergedT1 = false;
      outTree[0]->Fill();
    }
    */
    if(MODEL!=-1) PARSE_EVENT();
    outTree[0]->Fill();

    // cleanup
    PFPUcorrJet.clear();
  }

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  for(int i=0; i<3; i++) 
    outTree[i]->Write();
  file->Close();
}
  
void RazorBoostedTop::SetFirstTop(TLorentzVector myTop, int merged) {
  PFT1Pt = myTop.Pt();
  PFT1Eta = myTop.Eta();
  PFT1Phi = myTop.Phi();
  PFT1Mass = myTop.M();
  MergedT1  = merged;
}

void RazorBoostedTop::SetSecondTop(TLorentzVector myTop, int merged) {
  PFT2Pt = myTop.Pt();
  PFT2Eta = myTop.Eta();
  PFT2Phi = myTop.Phi();
  PFT2Mass = myTop.M();
  MergedT2  = merged;
}

vector<TLorentzVector> RazorBoostedTop::CombineJets(vector<TLorentzVector> myjets){
  
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
    double M_temp = j_temp1.M2()+j_temp2.M2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  // set masses to 0
  //j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  //j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);

  return mynewjets;  
}

void RazorBoostedTop::PARSE_EVENT(){
	
  std::vector<std::string>::const_iterator c_begin = commentLHE->begin();
  std::vector<std::string>::const_iterator c_end = commentLHE->end();

  m0=m12=tanb=A0=0;
  mu=1.0;
  double signMu;
  for( std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
    //cout << (*cit) << endl;
    size_t found = (*cit).find("model");
    if( found != std::string::npos)   {    
      //         std::cout << *cit << std::endl;  

      if(MODEL == 0){// mSUGRA
	size_t foundLength = (*cit).size();
	found = (*cit).find("=");
	std::string smaller = (*cit).substr(found+1,foundLength);
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	std::istringstream iss(smaller);
	iss >> m0;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> m12;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> tanb;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> A0;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> signMu;
	iss.clear();
	mu *= signMu;
      }
     
      if(MODEL == 1){ //T1bbbb
	size_t foundLength = (*cit).size();
	found = (*cit).find("T1bbbb");
	std::string smaller = (*cit).substr(found+1,foundLength);
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	std::istringstream iss(smaller);
	iss >> m0;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	iss.str(smaller);
	iss >> m0;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> m12;
	iss.clear();
      }
      if(MODEL == 2){ //T2bb
	size_t foundLength = (*cit).size();
	found = (*cit).find("T2bb");
	std::string smaller = (*cit).substr(found+1,foundLength);
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	std::istringstream iss(smaller);
	iss >> m0;
	iss.clear();
	//
	
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> m12;
	iss.clear();
      }
      if(MODEL == 3){ //T2tt
	size_t foundLength = (*cit).size();
	found = (*cit).find("T2tt");
	std::string smaller = (*cit).substr(found+1,foundLength);
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	std::istringstream iss(smaller);
	iss >> m0;
	iss.clear();
	//
	
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> m12;
	iss.clear();
      }
      if(MODEL == 5){ //T2
	size_t foundLength = (*cit).size();
	found = (*cit).find("T2tt");
	std::string smaller = (*cit).substr(found+1,foundLength);
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	std::istringstream iss(smaller);
	iss >> m0;
	iss.clear();
	//
	
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> m12;
	iss.clear();
      }
      if(MODEL == 4){ //T1tttt
	size_t foundLength = (*cit).size();
	found = (*cit).find("T1tttt");
	std::string smaller = (*cit).substr(found+1,foundLength);
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	std::istringstream iss(smaller);
	iss >> m0;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	//
	iss.str(smaller);
	iss >> m0;
	iss.clear();
	//
	found = smaller.find("_");
	smaller = smaller.substr(found+1,smaller.size());
	iss.str(smaller);
	iss >> m12;
	iss.clear();
      }
    }
  }
  //cout << MODEL << " " << m0 << " " << m12 << endl;	
}
