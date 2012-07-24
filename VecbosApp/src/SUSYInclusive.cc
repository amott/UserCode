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
#include "SUSYInclusive.hh"

SUSYInclusive::SUSYInclusive(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  
}

SUSYInclusive::SUSYInclusive(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;

  //To read good run list!
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

SUSYInclusive::~SUSYInclusive() {}

void SUSYInclusive::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void SUSYInclusive::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;


  int HLTbit;
  int HLTbitEmulator;

  // PF block
  int    pass2PFJet;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double massPFHem;
  int PFgoodR;
  double PFR;
  double PFMR;

  // CALO block
  int    pass2CALOJet;
  double pTCALOHem1;
  double etaCALOHem1;
  double phiCALOHem1;
  double pTCALOHem2;
  double etaCALOHem2;
  double phiCALOHem2;
  double massCALOHem;
  int CALOgoodR;
  double CALOR;
  double CALOMR;

  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("HLTbit", &HLTbit, "HLTbit/I");
  outTree->Branch("HLTbitEmulator", &HLTbitEmulator, "HLTbitEmulator/I");

  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  

  // PF block
  outTree->Branch("pass2PFJet", &pass2PFJet, "pass2PFJet/I");
  outTree->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
  outTree->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
  outTree->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
  outTree->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
  outTree->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
  outTree->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
  outTree->Branch("massPFHem", &massPFHem, "massPFHem/D");
  outTree->Branch("PFgoodR", &PFgoodR, "PFgoodR/I");
  outTree->Branch("PFR", &PFR, "PFR/D");
  outTree->Branch("PFMR", &PFMR, "PFMR/D");

  // CALO block
  outTree->Branch("pass2CALOJet", &pass2CALOJet, "pass2CALOJet/I");
  outTree->Branch("pTCALOHem1", &pTCALOHem1, "pTCALOHem1/D");
  outTree->Branch("etaCALOHem1", &etaCALOHem1, "etaCALOHem1/D");
  outTree->Branch("phiCALOHem1", &phiCALOHem1, "phiCALOHem1/D");
  outTree->Branch("pTCALOHem2", &pTCALOHem2, "pTCALOHem2/D");
  outTree->Branch("etaCALOHem2", &etaCALOHem2, "etaCALOHem2/D");
  outTree->Branch("phiCALOHem2", &phiCALOHem2, "phiCALOHem2/D");
  outTree->Branch("massCALOHem", &massCALOHem, "massCALOHem/D");
  outTree->Branch("CALOgoodR", &CALOgoodR, "CALOgoodR/I");
  outTree->Branch("CALOR", &CALOR, "CALOR/D");
  outTree->Branch("CALOMR", &CALOMR, "CALOMR/D");

  // prepare vectors for efficiency
  double NpassedPF_In = 0;
  double NpassedPF_HLT = 0; 
  double NpassedPF_PV = 0;
  double NpassedPF_2Jet = 0;
  double NpassedPF_CenJet = 0;
  double NpassedPF_beta = 0;
  double NpassedPF_DeltaPhi = 0;
  double NpassedPF_R = 0;

  double NpassedCALO_In = 0;
  double NpassedCALO_HLT = 0;
  double NpassedCALO_PV = 0;
  double NpassedCALO_2Jet = 0;
  double NpassedCALO_CenJet = 0;
  double NpassedCALO_beta = 0;
  double NpassedCALO_DeltaPhi = 0;
  double NpassedCALO_R = 0;

  double weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask(true);
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
    
    NpassedPF_In += weight;
    NpassedCALO_In += weight;
    
    // HLT 
    vector<int> myHLT = getHLTOutput();
    HLTbit = 0;
    if(myHLT.size() == 1) HLTbit = myHLT[0];
    HLTbitEmulator = HTtrigger(150.);

    NpassedPF_HLT += weight;
    NpassedCALO_HLT += weight;

    // vbtf muon
    bool goodmu =  false;
    for(int jj=0; jj<nMuon;jj++){
      double muonPt=sqrt(pxMuon[jj]*pxMuon[jj]+pyMuon[jj]*pyMuon[jj]);
	if(isTightMuon(jj)  && muonPt>15.){
	  goodmu = true;
      }
    }
    if(!goodmu) continue;

    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    NpassedPF_PV += weight;
    NpassedCALO_PV += weight;

    // Jets selection
    vector<TLorentzVector> PFJet;
    vector<TLorentzVector> CALOJet;
    
    for(int i=0; i< nAK5PFJet; i++) {
      TLorentzVector myJet(pxAK5PFJet[i], pyAK5PFJet[i], pzAK5PFJet[i], energyAK5PFJet[i]);   
      if(myJet.Pt()>30. && fabs(myJet.Eta())< 2.4) PFJet.push_back(myJet);
    }

    for(int i=0; i< nAK5Jet; i++) {
      TLorentzVector myJet(pxAK5Jet[i], pyAK5Jet[i], pzAK5Jet[i], energyAK5Jet[i]);   
      if(myJet.Pt()>30. && fabs(myJet.Eta())< 2.4) CALOJet.push_back(myJet);
    }

    if(CALOJet.size() < 2 && PFJet.size() < 2) continue;

    // dummy values
    pass2PFJet = 0;
    pTPFHem1 = -9999;
    etaPFHem1 = -9999;
    phiPFHem1 = -9999;
    pTPFHem2 = -9999;
    etaPFHem2 = -9999;
    phiPFHem2 = -9999;
    massPFHem = -99;
    PFgoodR = -1;
    PFR = -99999.;
    PFMR = -99999.;
    
    if(PFJet.size()>=2) {

     // simple jets
      NpassedPF_2Jet += weight;

      // first jet centrality
      if(fabs(PFJet[HighestPtJet(PFJet,-999)].Eta())> 1.7) continue;

      // simple jets
      NpassedPF_CenJet += weight;

      // hemispheres
      vector<TLorentzVector> tmpJet = CombineJets(PFJet);
      TLorentzVector PFHem1 = tmpJet[0];
      TLorentzVector PFHem2 = tmpJet[1];
      TLorentzVector DiHem = PFHem1 + PFHem2;

      // compute boost
      double num = PFHem1.P()-PFHem2.P();
      double den = PFHem1.Pz()-PFHem2.Pz();
      
      double beta;
      if(fabs(num) < fabs(den)) {
	// R good, R' bad
	beta = num/den;
	PFgoodR = 1;
      } else if(fabs(num) > fabs(den)) {
	beta = den/num;
	PFgoodR = 0;
      }

      if(fabs(beta)< 0.99) {
	NpassedPF_beta += weight;
	if(fabs(PFHem1.DeltaPhi(PFHem2)) < 2.8) {
	  NpassedPF_DeltaPhi += weight;

	  TVector3 MET(pxMet[0], pyMet[0], 0.);
	  double MT = CalcMT(PFHem1, PFHem2, MET);
	  double variable = -999999.;
	  double Rvariable = -999999.;
	  if(PFgoodR==1) {
	    // variable is R
	    variable = CalcMR(PFHem1, PFHem2);
	    if(variable >0) Rvariable = MT/variable;
	  } else if(PFgoodR==0) {
	    // variable is R'
	    variable = CalcMRP(PFHem1, PFHem2, MET);
	    if(variable >0) Rvariable = MT/variable;
	  }

	  // fill the tree
	  pass2PFJet = 1;
	  pTPFHem1 = PFHem1.Pt();
	  etaPFHem1 = PFHem1.Eta();
	  phiPFHem1 = PFHem1.Phi();
	  pTPFHem2 = PFHem2.Pt();
	  etaPFHem2 = PFHem2.Eta();
	  phiPFHem2 = PFHem2.Phi();
	  massPFHem = DiHem.M();
	  PFR = Rvariable;
	  PFMR = variable;
	}
      }
    }
    
    // dummy values
    pass2CALOJet = 0;
    pTCALOHem1 = -9999;
    etaCALOHem1 = -9999;
    phiCALOHem1 = -9999;
    pTCALOHem2 = -9999;
    etaCALOHem2 = -9999;
    phiCALOHem2 = -9999;
    massCALOHem = -99;
    CALOgoodR = -1;
    CALOR = -99999.;
    CALOMR = -99999.;
    
    if(CALOJet.size()>=2) {

      // simple jets
      NpassedCALO_2Jet += weight;

      // first jet centrality
      if(fabs(CALOJet[HighestPtJet(CALOJet,-999)].Eta())> 1.7) continue;

      // simple jets
      NpassedCALO_CenJet += weight;

      // hemispheres
      vector<TLorentzVector> tmpJet = CombineJets(CALOJet);
      TLorentzVector CALOHem1 = tmpJet[0];
      TLorentzVector CALOHem2 = tmpJet[1];
      TLorentzVector DiHem = CALOHem1 + CALOHem2;

      // compute boost
      double num = CALOHem1.P()-CALOHem2.P();
      double den = CALOHem1.Pz()-CALOHem2.Pz();
      
      double beta;
      if(fabs(num) < fabs(den)) {
	// R good, R' bad
	beta = num/den;
	CALOgoodR = 1;
      } else if(fabs(num) > fabs(den)) {
	beta = den/num;
	CALOgoodR = 0;
      }
      if(fabs(beta)<0.99) {
	NpassedCALO_beta += weight;
	if(fabs(CALOHem1.DeltaPhi(CALOHem2)) <= 2.8) {
	  NpassedCALO_DeltaPhi += weight;
	  
	  TVector3 MET(pxMet[0], pyMet[0], 0.);
	  double MT = CalcMT(CALOHem1, CALOHem2, MET);
	  double variable = -999999.;
	  double Rvariable = -999999.;
	  if(CALOgoodR==1) {
	    // variable is R
	    variable = CalcMR(CALOHem1, CALOHem2);
	    if(variable >0) Rvariable = MT/variable;
	  } else if(CALOgoodR==0) {
	    // variable is R'
	    variable = CalcMRP(CALOHem1, CALOHem2, MET);
	    if(variable >0) Rvariable = MT/variable;
	  }
	  
	  // fill the tree
	  pass2CALOJet = 1;
	  pTCALOHem1 = CALOHem1.Pt();
	  etaCALOHem1 = CALOHem1.Eta();
	  phiCALOHem1 = CALOHem1.Phi();
	  pTCALOHem2 = CALOHem2.Pt();
	  etaCALOHem2 = CALOHem2.Eta();
	  phiCALOHem2 = CALOHem2.Phi();
	  massCALOHem = DiHem.M();
	  CALOR = Rvariable;
	  CALOMR = variable;
	}
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

  /*
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
    
  effTree->Branch("NpassedPF_In",      &NpassedPF_In,      "NpassedPF_In/D");
  effTree->Branch("NpassedPF_HLT",      &NpassedPF_HLT,      "NpassedPF_HLT/D");
  effTree->Branch("NpassedPF_PV",      &NpassedPF_PV,      "NpassedPF_PV/D");
  effTree->Branch("NpassedPF_2Jet",      &NpassedPF_2Jet,      "NpassedPF_2Jet/D");
  effTree->Branch("NpassedPF_CenJet",      &NpassedPF_CenJet,      "NpassedPF_CenJet/D");
  effTree->Branch("NpassedPF_beta",      &NpassedPF_beta,      "NpassedPF_beta/D");
  effTree->Branch("NpassedPF_DeltaPhi",      &NpassedPF_DeltaPhi,      "NpassedPF_DeltaPhi/D");
  effTree->Branch("NpassedPF_R",      &NpassedPF_R,      "NpassedPF_R/D");
  
  effTree->Branch("NpassedCALO_In",      &NpassedCALO_In,      "NpassedCALO_In/D");
  effTree->Branch("NpassedCALO_HLT",      &NpassedCALO_HLT,      "NpassedCALO_HLT/D");
  effTree->Branch("NpassedCALO_PV",      &NpassedCALO_PV,      "NpassedCALO_PV/D");
  effTree->Branch("NpassedCALO_2Jet",      &NpassedCALO_2Jet,      "NpassedCALO_2Jet/D");
  effTree->Branch("NpassedCALO_CenJet",      &NpassedCALO_CenJet,      "NpassedCALO_CenJet/D");
  effTree->Branch("NpassedCALO_beta",      &NpassedCALO_beta,      "NpassedCALO_beta/D");
  effTree->Branch("NpassedCALO_DeltaPhi",      &NpassedCALO_DeltaPhi,      "NpassedCALO_DeltaPhi/D");
  effTree->Branch("NpassedCALO_R",      &NpassedCALO_R,      "NpassedCALO_R/D");
  effTree->Fill();
  */

  cout << outTree->GetEntries() << endl;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  //  effTree->Write();
  outTree->Write();
  file->Close();
}
  
int SUSYInclusive::HighestPtJet(vector<TLorentzVector> Jet, int firstJet) {

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

vector<TLorentzVector> SUSYInclusive::JetRecovery(vector<TLorentzVector> Jet, int iJ1, int iJ2, int dRRatio, double pTthreshold) {

  TLorentzVector JR1 = Jet[iJ1];
  TLorentzVector JR2 = Jet[iJ2];
  
  //  recover the second jet first
  for(int i=0; i < Jet.size(); i++) {
    if(i == iJ1 || i == iJ2) continue;
    double dR1 = Jet[iJ1].DeltaR(Jet[i]);
    double dR2 = Jet[iJ2].DeltaR(Jet[i]);
    if(dR1 < dR2) { // closest to first jet
      if(dR1 < 0.7*dRRatio && Jet[i].Pt() > pTthreshold) {
	JR1 = JR1 + Jet[i];
      }
    } else { // closest to second jet
      if(dR2 < 0.7*dRRatio && Jet[i].Pt() > pTthreshold) {
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

vector<TLorentzVector> SUSYInclusive::CombineJets(vector<TLorentzVector> myjets){

  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;

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
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }
  
  j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }

  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;
  
}

double SUSYInclusive::SumPt(int iMu) {
  double eta0 = etaMuon[iMu];
  double phi0 = phiMuon[iMu];
  double sumPt_tmp = 0;
  for(int i=0; i< nTrack; i++) {
    if(i == trackIndexMuon[iMu]) continue; // take out the muon
    if(trackValidHitsTrack[i] <5) continue;                                     // minimum number of hits  XXX
    if(eleDxyPV(PVxPV[0], PVyPV[0], PVzPV[0], trackVxTrack[i], trackVyTrack[i], trackVzTrack[i], pxTrack[i], pyTrack[i], pzTrack[i]) > 0.04) continue; // track incompatible with the vertex on (x,y)
    if(eleDszPV(PVxPV[0], PVyPV[0], PVzPV[0], trackVxTrack[i], trackVyTrack[i], trackVzTrack[i], pxTrack[i], pyTrack[i], pzTrack[i]) > 0.12) continue; // track incompatible with the vertex on z
    TVector3 v(pxTrack[i], pyTrack[i], pzTrack[i]);
    if(sqrt(pow(v.Eta()-eta0,2.)+pow(v.Phi()-phi0,2.)) > 0.5) continue; // track outside the cone
    if(v.Pt() < 0.500) continue;     // minimum pT             
    if(v.Pt() > 500.) continue;     // maximum pT             
    sumPt_tmp += v.Pt();
  }
  return sumPt_tmp; 
}

bool SUSYInclusive::isTightMuon(int iMu){
  
  bool ret = false;
  Utils anaUtils;
  bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
  bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
  bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);
  
  if(isMuGlobal && isMuGlobalPrompt && isMuTracker){

    double pt = sqrt(pxMuon[iMu]*pxMuon[iMu]+pyMuon[iMu]*pyMuon[iMu]);
    int iTrack = trackIndexMuon[iMu];
    if(numberOfValidStripTIBHitsTrack[iTrack]+
       numberOfValidStripTIDHitsTrack[iTrack]+
       numberOfValidStripTOBHitsTrack[iTrack]+
       numberOfValidStripTECHitsTrack[iTrack] > 10){
      
      if(numberOfValidPixelBarrelHitsTrack[iTrack] > 0 ||
	 numberOfValidPixelEndcapHitsTrack[iTrack] > 0){
	if(fabs(transvImpactParTrack[iTrack]) < 0.2){
	  //if(trackNormalizedChi2GlobalMuonTrack[iTrack] < 10){                                                                                           

	  double IECAL = emEt03Muon[iMu];
	  double IHCAL = hadEt03Muon[iMu];
	  double ITRK  = sumPt03Muon[iMu];
	  
	  if((IECAL+IHCAL+ITRK)/pt < 0.15){

	    ret = true;
	    
	  }
	}
      }
    }
  }

  return ret;
}

bool SUSYInclusive::HTtrigger(double HTmin){
 
  TVector3 TRIG_MET;
  TRIG_MET.SetXYZ(pxMet[0],pyMet[0],0.0); //TVector3 type
  TLorentzVector TRIG_j1;
  TRIG_j1.SetPxPyPzE(0.0,0.0,0.0,0.0); // ditto 
  TLorentzVector TRIG_j2;
  TRIG_j2.SetPxPyPzE(0.0,0.0,0.0,0.0); // ditto
  double modHT20 = 0.0;
  double pt_0 = 0.0;
  double pt_1 = 0.0;
  for(int i = 0; i < nAK5Jet; i++){
    TLorentzVector jet;
    double px = pxAK5Jet[i];
    double py = pyAK5Jet[i];
    double pz = pzAK5Jet[i];
    double E = sqrt(px*px+py*py+pz*pz);
    double scale = 1.;
    scale = uncorrEnergyAK5Jet[i]/energyAK5Jet[i]; // get the uncorrected energy 
    
    jet.SetPxPyPzE(px*scale,py*scale,pz*scale,E*scale); // massless jets
    
    //if(fabs(jet.Eta()) >= 3.0) continue; //for now, unrestricted eta b/c this is what is implemented @ trigger
    //Get the two leading jets (for emulation of triggers based on these - not for the HT trigger
    if(jet.Pt() > pt_0){
      pt_1 = pt_0;
      TRIG_j2 = TRIG_j1;
      pt_0 = jet.Pt();
      TRIG_j1 = jet;
    } else {
      if(jet.Pt() > pt_1){
	pt_1 = jet.Pt();
	TRIG_j2 = jet;
      }
    }
    //Include only uncorredted jets with pt > 20 GeV in HT calculation
    if(jet.Pt() > 20.0){
      if(jet.Pt() > 20.) modHT20 += jet.Pt();
    }
  }
  bool passed = (modHT20> HTmin ? true : false);
  return passed;
}

