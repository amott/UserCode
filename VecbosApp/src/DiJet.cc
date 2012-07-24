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
#include "DiJet.hh"

DiJet::DiJet(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  
}

DiJet::DiJet(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;

  //To read good run list!
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "config/vecbos/json/lumiSummary_HTRun2011Av2.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

DiJet::~DiJet() {}

void DiJet::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");

  // PF block
  int    pass2PFJet;
  double nPFjets;
  double pTPFJet1;
  double etaPFJet1;
  double phiPFJet1;
  double pTPFJet2;
  double etaPFJet2;
  double phiPFJet2;
  double massPFJet;
  double etPFMET;
  double phiPFMET;

  double massDiPF1p1;
  double pTPF1p1_1;
  double etaPF1p1_1;
  double phiPF1p1_1;
  double massPF1p1_1;
  double pTPF1p1_2;
  double etaPF1p1_2;
  double phiPF1p1_2;
  double massPF1p1_2;

  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;

  // Prescaled Jet Triggers
  int HLT_DiJetAve30;
  int HLT_DiJetAve60;
  int HLT_Jet30;
  int HLT_Jet60;
  // HT triggers
  int HLT_HT150;
  int HLT_HT160;
  int HLT_HT200;
  int HLT_HT240;
  int HLT_HT250;
  int HLT_HT260;
  int HLT_HT300;
  int HLT_HT400;
  int HLT_HT450;
  int HLT_HT500;
  int HLT_HT550;
  // Fatjet Mass
  int HLT_FatJetMass750_DR1p1_Deta2p0;
  int HLT_FatJetMass850_DR1p1_Deta2p0;

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");

  // PF block
  outTree->Branch("pass2PFJet", &pass2PFJet, "pass2PFJet/I");
  outTree->Branch("nPFjets", &nPFjets, "nPFjets/D");
  outTree->Branch("pTPFJet1", &pTPFJet1, "pTPFJet1/D");
  outTree->Branch("etaPFJet1", &etaPFJet1, "etaPFJet1/D");
  outTree->Branch("phiPFJet1", &phiPFJet1, "phiPFJet1/D");
  outTree->Branch("pTPFJet2", &pTPFJet2, "pTPFJet2/D");
  outTree->Branch("etaPFJet2", &etaPFJet2, "etaPFJet2/D");
  outTree->Branch("phiPFJet2", &phiPFJet2, "phiPFJet2/D");
  outTree->Branch("massPFJet", &massPFJet, "massPFJet/D");
  outTree->Branch("etPFMET", &etPFMET, "etPFMET/D");
  outTree->Branch("phiPFMET", &phiPFMET, "phiPFMET/D");

  outTree->Branch("massDiPF1p1", &massDiPF1p1, "massDiPF1p1/D");
  outTree->Branch("pTPF1p1_1", &pTPF1p1_1, "pTPF1p1_1/D");
  outTree->Branch("etaPF1p1_1", &etaPF1p1_1, "etaPF1p1_1/D");
  outTree->Branch("phiPF1p1_1", &phiPF1p1_1, "phiPF1p1_1/D");
  outTree->Branch("massPF1p1_1", &massPF1p1_1, "massPF1p1_1/D");
  outTree->Branch("pTPF1p1_2", &pTPF1p1_2, "pTPF1p1_2/D");
  outTree->Branch("etaPF1p1_2", &etaPF1p1_2, "etaPF1p1_2/D");
  outTree->Branch("phiPF1p1_2", &phiPF1p1_2, "phiPF1p1_2/D");
  outTree->Branch("massPF1p1_2", &massPF1p1_2, "massPF1p1_2/D");

  // Prescaled Jet Triggers
  outTree->Branch("HLT_DiJetAve30", &HLT_DiJetAve30, "HLT_DiJetAve30/I");
  outTree->Branch("HLT_DiJetAve60", &HLT_DiJetAve60, "HLT_DiJetAve60/I");
  outTree->Branch("HLT_Jet30", &HLT_Jet30, "HLT_Jet30/I");
  outTree->Branch("HLT_Jet60", &HLT_Jet60, "HLT_Jet60/I");
  // HT triggers
  outTree->Branch("HLT_HT150", &HLT_HT150, "HLT_HT150/I");
  outTree->Branch("HLT_HT160", &HLT_HT160, "HLT_HT160/I");
  outTree->Branch("HLT_HT200", &HLT_HT200, "HLT_HT200/I");
  outTree->Branch("HLT_HT240", &HLT_HT240, "HLT_HT240/I");
  outTree->Branch("HLT_HT250", &HLT_HT250, "HLT_HT250/I");
  outTree->Branch("HLT_HT260", &HLT_HT260, "HLT_HT260/I");
  outTree->Branch("HLT_HT300", &HLT_HT300, "HLT_HT300/I");
  outTree->Branch("HLT_HT400", &HLT_HT400, "HLT_HT400/I");
  outTree->Branch("HLT_HT450", &HLT_HT450, "HLT_HT450/I");
  outTree->Branch("HLT_HT500", &HLT_HT500, "HLT_HT500/I");
  outTree->Branch("HLT_HT550", &HLT_HT550, "HLT_HT550/I");
  // Fatjet Mass
  outTree->Branch("HLT_FatJetMass750_DR1p1_Deta2p0", &HLT_FatJetMass750_DR1p1_Deta2p0, "HLT_FatJetMass750_DR1p1_Deta2p0/I");
  outTree->Branch("HLT_FatJetMass850_DR1p1_Deta2p0", &HLT_FatJetMass850_DR1p1_Deta2p0, "HLT_FatJetMass850_DR1p1_Deta2p0/I");

  TH1D* massAll = new TH1D("massAll", "massAll", 500., 0., 3500.);
  TH1D* massMU  = new TH1D("massMU", "massMU", 500., 0., 3500.);

  // prepare vectors for efficiency
  int NpassedPF_In = 0;
  int NpassedPF_PV = 0;
  int NpassedPF_2Jet = 0;
  int NpassedPF_Mass = 0;

  double weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  std::vector<std::string> mask; mask.push_back("_v");
  // Prescaled Jet Triggers
  std::vector<std::string> maskHLT_DiJetAve30; mask.push_back("HLT_DiJetAve30_v");
  std::vector<std::string> maskHLT_DiJetAve60; mask.push_back("HLT_DiJetAve60_v");
  std::vector<std::string> maskHLT_Jet30; mask.push_back("HLT_Jet30_v");
  std::vector<std::string> maskHLT_Jet60; mask.push_back("HLT_Jet60_v");
  // HT triggers
  std::vector<std::string> maskHLT_HT150; mask.push_back("HLT_HT150_v");
  std::vector<std::string> maskHLT_HT160; mask.push_back("HLT_HT160_v");
  std::vector<std::string> maskHLT_HT200; mask.push_back("HLT_HT200_v");
  std::vector<std::string> maskHLT_HT240; mask.push_back("HLT_HT240_v");
  std::vector<std::string> maskHLT_HT250; mask.push_back("HLT_HT250_v");
  std::vector<std::string> maskHLT_HT260; mask.push_back("HLT_HT260_v");
  std::vector<std::string> maskHLT_HT300; mask.push_back("HLT_HT300_v");
  std::vector<std::string> maskHLT_HT400; mask.push_back("HLT_HT400_v");
  std::vector<std::string> maskHLT_HT450; mask.push_back("HLT_HT450_v");
  std::vector<std::string> maskHLT_HT500; mask.push_back("HLT_HT500_v");
  std::vector<std::string> maskHLT_HT550; mask.push_back("HLT_HT550_v");
  // Fatjet Mass
  std::vector<std::string> maskHLT_FatJetMass750_DR1p1_Deta2p0; mask.push_back("HLT_FatJetMass750_DR1p1_Deta2p0_v");
  std::vector<std::string> maskHLT_FatJetMass850_DR1p1_Deta2p0; mask.push_back("HLT_FatJetMass850_DR1p1_Deta2p0_v");

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }  continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    if(_isData) {

      // Prescaled Jet Triggers
      setRequiredTriggers(maskHLT_DiJetAve30);  reloadTriggerMask(true); HLT_DiJetAve30 = hasPassedHLT();
      setRequiredTriggers(maskHLT_DiJetAve60);  reloadTriggerMask(true); HLT_DiJetAve60 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Jet30);  reloadTriggerMask(true); HLT_Jet30 = hasPassedHLT();
      setRequiredTriggers(maskHLT_Jet60);  reloadTriggerMask(true); HLT_Jet60 = hasPassedHLT();
      // HT triggers
      setRequiredTriggers(maskHLT_HT150);  reloadTriggerMask(true); HLT_HT150 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT160);  reloadTriggerMask(true); HLT_HT160 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT200);  reloadTriggerMask(true); HLT_HT200 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT240);  reloadTriggerMask(true); HLT_HT240 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT250);  reloadTriggerMask(true); HLT_HT250 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT260);  reloadTriggerMask(true); HLT_HT260 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT300);  reloadTriggerMask(true); HLT_HT300 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT400);  reloadTriggerMask(true); HLT_HT400 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT450);  reloadTriggerMask(true); HLT_HT450 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT500);  reloadTriggerMask(true); HLT_HT500 = hasPassedHLT();
      setRequiredTriggers(maskHLT_HT550);  reloadTriggerMask(true); HLT_HT550 = hasPassedHLT();
      // Fatjet Mass 
      setRequiredTriggers(maskHLT_FatJetMass750_DR1p1_Deta2p0);  reloadTriggerMask(true); HLT_FatJetMass750_DR1p1_Deta2p0 = hasPassedHLT();
      setRequiredTriggers(maskHLT_FatJetMass850_DR1p1_Deta2p0);  reloadTriggerMask(true); HLT_FatJetMass850_DR1p1_Deta2p0 = hasPassedHLT();

    }

    NpassedPF_In += int(weight);


    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 15.) continue; 

    NpassedPF_PV += int(weight);
    
    vector<TLorentzVector> PFJet;
    for(int i=0; i< nAK5PFNoPUJet; i++) {
      TLorentzVector myJet(pxAK5PFNoPUJet[i], pyAK5PFNoPUJet[i], pzAK5PFNoPUJet[i], energyAK5PFNoPUJet[i]);   
      if(myJet.Pt()>30. && fabs(myJet.Eta())< 2.5) PFJet.push_back(myJet);
    }

    if(PFJet.size()<2)  continue;

    // simple jets
    NpassedPF_2Jet += int(weight);
    int iFirstJet = HighestPtJet(PFJet, -99);
    int iSecondJet = HighestPtJet(PFJet, iFirstJet);
    TLorentzVector PFJet1 = PFJet[iFirstJet];
    TLorentzVector PFJet2 = PFJet[iSecondJet];
    TLorentzVector DiJet = PFJet1+PFJet2;

    // jet recovery
    vector<TLorentzVector> Jet25Rec = JetRecovery(PFJet, iFirstJet, iSecondJet, 1.1, 30.);
    TLorentzVector PF1p1_1 = Jet25Rec[0];
    TLorentzVector PF1p1_2 = Jet25Rec[1];
    TLorentzVector Di1p1_ = PF1p1_1 + PF1p1_2;
    
    massAll->Fill(Di1p1_.M());

    for(int i = 0; i<nMuon; i++) {
      TLorentzVector mu(pxMuon[i], pyMuon[i], pzMuon[i], 0.);
      if(mu.Pt()<30.) continue;
      if(fabs(mu.Eta())>2.4) continue;
      massMU->Fill(Di1p1_.M());
    }

    pass2PFJet = 0;
    nPFjets = PFJet.size();
    pTPFJet1 = PFJet1.Pt();
    etaPFJet1 = PFJet1.Eta();
    phiPFJet1 = PFJet1.Phi();
    pTPFJet2 = PFJet2.Pt();
    etaPFJet2 = PFJet2.Eta();
    phiPFJet2 = PFJet2.Phi();
    massPFJet = DiJet.M();
    etPFMET = energyPFMet[0];
    phiPFMET = phiPFMet[0];
    
    massDiPF1p1 = (PF1p1_1+PF1p1_2).M();
    pTPF1p1_1   = PF1p1_1.Pt();
    etaPF1p1_1  = PF1p1_1.Eta();
    phiPF1p1_1  = PF1p1_1.Phi();
    massPF1p1_1 = PF1p1_1.M();
    pTPF1p1_2   = PF1p1_2.Pt();
    etaPF1p1_2  = PF1p1_2.Eta();
    phiPF1p1_2  = PF1p1_2.Phi();
    massPF1p1_2 = PF1p1_1.M();
    
    // fill output tree
    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;

    // PF block
    if(fabs(PFJet1.Eta()-PFJet2.Eta())<1.3) { 
      pass2PFJet = 1;
      NpassedPF_Mass += int(weight);
    }
    outTree->Fill();

  }

  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
    
  effTree->Branch("NpassedPF_In",      &NpassedPF_In,      "NpassedPF_In/I");
  effTree->Branch("NpassedPF_PV",      &NpassedPF_PV,      "NpassedPF_PV/I");
  effTree->Branch("NpassedPF_2Jet",      &NpassedPF_2Jet,      "NpassedPF_2Jet/I");
  effTree->Branch("NpassedPF_Mass",      &NpassedPF_Mass,      "NpassedPF_Mass/I");
  
  effTree->Fill();
  effTree->Write();
  outTree->Write();
  file->Close();
}
  
int DiJet::HighestPtJet(vector<TLorentzVector> Jet, int firstJet) {

  int index=-99;
  double pT=-999999.;
  for(int i=0; i<Jet.size(); i++) {
    if(i == firstJet) continue;
    if(Jet[i].Pt()>pT) {
      pT = Jet[i].Pt();
      index = i;
    }
  }
  return index;
}

vector<TLorentzVector> DiJet::JetRecovery(vector<TLorentzVector> Jet, int iJ1, int iJ2, double dR, double pTthreshold) {

  TLorentzVector JR1 = Jet[iJ1];
  TLorentzVector JR2 = Jet[iJ2];
  
  //  recover the second jet first
  for(int i=0; i < Jet.size(); i++) {
    if(i == iJ1 || i == iJ2) continue;
    double dR1 = fabs(Jet[iJ1].DeltaR(Jet[i]));
    double dR2 = fabs(Jet[iJ2].DeltaR(Jet[i]));
    if(dR1 < dR2) { // closest to first jet
      if(dR1 < dR && Jet[i].Pt() > pTthreshold) {
	JR1 = JR1 + Jet[i];
      }
    } else { // closest to second jet
      if(dR2 < dR && Jet[i].Pt() > pTthreshold) {
	JR2 = JR2 + Jet[i];
      }
    } 
  }
  
  vector<TLorentzVector> newJets;
  if(JR1.Pt() > JR2.Pt()) {
    newJets.push_back(JR1);
    newJets.push_back(JR2);
  } else {
    newJets.push_back(JR2);
    newJets.push_back(JR1);
  }
  
  return newJets;
}

