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
#include "SUSYMultiB.hh"

SUSYMultiB::SUSYMultiB(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
}

SUSYMultiB::SUSYMultiB(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "config/vecbos/json/Golden_1Oct.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

SUSYMultiB::~SUSYMultiB() {}

void SUSYMultiB::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void SUSYMultiB::SetWeight(double weight){
  _weight=weight;
}

int SUSYMultiB::HighestPtJet(vector<TLorentzVector> Jet, int firstJet) {

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

bool SUSYMultiB::isLooseMuon(int iMu){

  bool ret = false;
  Utils anaUtils;
  bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);

  if(isMuGlobal){

    int iTrack = trackIndexMuon[iMu];
    if(numberOfValidStripTIBHitsTrack[iTrack]+
       numberOfValidStripTIDHitsTrack[iTrack]+
       numberOfValidStripTOBHitsTrack[iTrack]+
       numberOfValidStripTECHitsTrack[iTrack] > 10){

      ret = true;

    }
  }
  return ret;
}

bool SUSYMultiB::isGlobalMuon(int iMu){
  bool ret = false;
  Utils anaUtils;
  bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
  bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
  bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);

  if(isMuGlobal && isMuGlobalPrompt && isMuTracker)ret=true;

  return ret;
}


bool SUSYMultiB::isTightMuon(int iMu){
  
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

	  double CombinedIso= (IECAL+IHCAL+ITRK) - rhoFastjet * TMath::Pi() * 0.3 * 0.3;

	  if(CombinedIso/pt < 0.3){
	    ret = true;
	  }
	}
      }
    }
  }

  return ret;
}


bool SUSYMultiB::is80Electron(int iEle){

  Utils anaUtils;

  //is an ECAL driven electron                                                                                                                               
  bool isECALdriven = anaUtils.electronRecoType(recoFlagsEle[iEle], bits::isEcalDriven);

  if(isECALdriven == false) return false;

  bool isBarrel = true;
  bool isEndcap = true;
  int iSC = superClusterIndexEle[iEle];
  double ETA = etaSC[iSC];
  double ET = energySC[iSC]/cosh(etaSC[iSC]);

  if(ET <= 20.0) return false;

  //is the electron in the ECAL fiducial region?                                                                                                             
  if(fabs(ETA) > 1.4442) isBarrel = false;
  if(fabs(ETA) > 2.5 || fabs(ETA) < 1.566) isEndcap = false;

  if(isBarrel == false && isEndcap == false) return false;

  double trackISO = dr03TkSumPtEle[iEle];
  double ECALISO = dr03EcalRecHitSumEtEle[iEle];
  double HCALISO = dr03HcalTowerSumEtEle[iEle];

  double pt = sqrt(pxEle[iEle]*pxEle[iEle]+pyEle[iEle]*pyEle[iEle]);

  double sigietaieta = covIEtaIEtaSC[iSC];
  double dphi = deltaPhiAtVtxEle[iEle];
  double deta = deltaEtaAtVtxEle[iEle];
  double HoE = hOverEEle[iEle];
  double convDist = convDistEle[iEle];
  double convDcot = convDcotEle[iEle];

  if(fabs(convDist) <= 0.02 && fabs(convDcot) <= 0.02) return false;

  int iTrack = gsfTrackIndexEle[iEle];

  if(expInnerLayersGsfTrack[iTrack] > 0) return false;

  trackISO /= ET;
  ECALISO /= ET;
  HCALISO /= ET;


  if(isBarrel){
    //ISO WP80                                                                                                                                               
    //    if(trackISO > 0.09) return false;
    //if(ECALISO > 0.07) return false;
    //if(HCALISO > 0.10) return false;

    //ID WP80                                                                                                                                                
    if(fabs(dphi) > 0.06) return false;
    if(fabs(deta) > 0.004) return false;
    if(HoE > 0.04) return false;
    if(sigietaieta > 0.01) return false;
  } else {
    //ISO WP80                                                                                                                                               
    //if(trackISO > 0.04) return false;
    //if(ECALISO > 0.05) return false;
    //if(HCALISO > 0.025) return false;

    //ID WP80                                                                                                                                                
    if(fabs(dphi) > 0.03) return false;
    //if(fabs(deta) > 0.007) return false;                                                                                                                   
    if(HoE > 0.025) return false;
    if(sigietaieta > 0.03) return false;
  }

  return true;
}


bool SUSYMultiB::is95Electron(int iEle){

  Utils anaUtils;

  //is an ECAL driven electron                                                                                                                               
  bool isECALdriven = anaUtils.electronRecoType(recoFlagsEle[iEle], bits::isEcalDriven);
  if(isECALdriven == false) return false;

  bool isBarrel = true;
  bool isEndcap = true;
  int iSC = superClusterIndexEle[iEle];
  double ETA = etaSC[iSC];
  double ET = energySC[iSC]/cosh(etaSC[iSC]);

  if(ET <= 20.0) return false;

  //is the electron in the ECAL fiducial region?                                                                                                             
  if(fabs(ETA) > 1.4442) isBarrel = false;
  if(fabs(ETA) > 2.5 || fabs(ETA) < 1.566) isEndcap = false;

  if(isBarrel == false && isEndcap == false) return false;

  double trackISO = dr03TkSumPtEle[iEle];
  double ECALISO = dr03EcalRecHitSumEtEle[iEle];
  double HCALISO = dr03HcalTowerSumEtEle[iEle];

  double pt = sqrt(pxEle[iEle]*pxEle[iEle]+pyEle[iEle]*pyEle[iEle]);

  double sigietaieta = covIEtaIEtaSC[iSC];
  double deta = deltaEtaAtVtxEle[iEle];
  double HoE = hOverEEle[iEle];
  int iTrack = gsfTrackIndexEle[iEle];

  if(expInnerLayersGsfTrack[iTrack] > 1) return false;

  trackISO /= ET;
  ECALISO /= ET;
  HCALISO /= ET;


  if(isBarrel){
    //ISO WP95                                                                                                                                               
    if(trackISO > 0.15) return false;
    if(ECALISO > 2.0) return false;
    if(HCALISO > 0.12) return false;

    //ID WP95                                                                                                                                                

    if(fabs(deta) > 0.007) return false;
    if(HoE > 0.15) return false;
    if(sigietaieta > 0.01) return false;
  } else {
    //ISO WP95                                                                                                                                               
    if(trackISO > 0.08) return false;
    if(ECALISO > 0.06) return false;
    if(HCALISO > 0.05) return false;

    //ID WP95                                                                                                                                                

    if(HoE > 0.07) return false;
    if(sigietaieta > 0.03) return false;
  }

  return true;
}


void SUSYMultiB::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  int HLT_Mu9;
  int HLT_Mu11;
  int HLT_Mu15;

  int HLT_Ele10_LW_L1R;
  int HLT_Ele15_SW_L1R;
  int HLT_Ele15_LW_L1R;
  int HLT_Ele15_SW_CaloEleId_L1R; 
  int HLT_Ele17_SW_CaloEleId_L1R; 
  int HLT_Ele17_SW_TightEleId_L1R;
  int HLT_Ele17_SW_TighterEleIdIsol_L1R_v2;
  int HLT_Ele17_SW_TighterEleIdIsol_L1R_v3;

  // PF block
  double MET;
  double ST;
  double SLT;
  int nJets;
  int JetMultiplicity;
  int EleFakeJetMultiplicity;
  int MuFakeJetMultiplicity;
  int NLooseMuons, NTightMuons; 
  int N80Eles, N95Eles;
  int nEleJets;
  int nBTagJets;
  double ptJets[100];
  double HT;
  double EleJetDR[1000];
  int PFgoodR;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double massPFHem;
  double PFR;
  double PFMR;
  double PFMT;
  double phiMETLept;
  double etaMETLept;
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double CSV_TopFour;
  double TCHE_TopFour;
  double SSV_TopFour;
  double JP_TopFour;
  double JBP_TopFour;
  double MuIso[2];
  double EleTrackIso[2], EleECALIso[2], EleHCALIso[2];
  int n_PV;
  double weight=_weight;
  double mg, mst, mchi;

  //prepare tree for SMS eff
  TTree* FullSMSTree = new TTree("FullSMSTree", "FullSMSTree");
  FullSMSTree->Branch("mg", &mg, "mg/D");
  FullSMSTree->Branch("mchi", &mchi, "mchi/D");

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");

  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("nPV", &nPV, "nPV/I");

  // PF block
  outTree->Branch("MET", &MET, "MET/D");
  outTree->Branch("phiMETLept", &phiMETLept, "phiMETLept/D");
  outTree->Branch("etaMETLept", &etaMETLept, "etaMETLept/D");
  outTree->Branch("ST", &ST, "ST/D");
  outTree->Branch("SLT", &SLT, "SLT/D");
  outTree->Branch("HT", &HT, "HT/D");
  outTree->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
  outTree->Branch("nBTagJets", &nBTagJets, "nBTagJets/I");
  outTree->Branch("N80Eles", &N80Eles, "N80Eles/I");
  outTree->Branch("N95Eles", &N95Eles, "N95Eles/I");
  outTree->Branch("NLooseMuons", &NLooseMuons, "NLooseMuons/I");
  outTree->Branch("NTightMuons", &NTightMuons, "NTightMuons/I");
  outTree->Branch("EleFakeJetMultiplicity", &EleFakeJetMultiplicity, "EleFakeJetMultiplicity/I");
  outTree->Branch("MuFakeJetMultiplicity", &MuFakeJetMultiplicity, "MuFakeJetMultiplicity/I");
  outTree->Branch("PFR", &PFR, "PFR/D");
  outTree->Branch("PFMR", &PFMR, "PFMR/D");
  outTree->Branch("PFMT", &PFMT, "PFMT/D");
  outTree->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
  outTree->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
  outTree->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
  outTree->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
  outTree->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
  outTree->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
  outTree->Branch("massPFHem", &massPFHem, "massPFHem/D");
  outTree->Branch("PFgoodR", &PFgoodR, "PFgoodR/I");
  outTree->Branch("CSV_TopFour", &CSV_TopFour, "CSV_TopFour/D");
  outTree->Branch("TCHE_TopFour", &TCHE_TopFour, "TCHE_TopFour/D");
  outTree->Branch("SSV_TopFour", &SSV_TopFour, "SSV_TopFour/D");
  outTree->Branch("JP_TopFour", &JP_TopFour, "JP_TopFour/D");
  outTree->Branch("JBP_TopFour", &JBP_TopFour, "JBP_TopFour/D");
  outTree->Branch("MuIso", MuIso, "MuIso[NLooseMuons]/D");
  outTree->Branch("EleTrackIso", EleTrackIso, "EleTrackIso[N95Eles]/D");
  outTree->Branch("EleECALIso", EleECALIso, "EleECALIso[N95Eles]/D");
  outTree->Branch("EleHCALIso", EleHCALIso, "EleHCALIso[N95Eles]/D");
  outTree->Branch("weight", &weight, "weight/D");
  outTree->Branch("mg", &mg, "mg/D");
  outTree->Branch("mst", &mst, "mst/D");
  outTree->Branch("mchi", &mchi, "mchi/D");

  // prepare vectors for efficiency
  double Npassed_In = 0;
  //HLT
  double Npassed_HLT = 0;
  
  double Npassed_PV = 0;
  double Npassed_MET= 0;
  double Npassed_ST= 0;

  //Jets
  double NpassedPF_4Jet = 0;
  
  //Hemispheres
  double NpassedPF_beta=0;
  double NpassedPF_DeltaPhi=0;
  //Leptons
  double Npassed_NoLept=0;
  //B-tag
  double Npassed_3b=0;

  double weightII = 1.;
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
    
    if(_isData)reloadTriggerMask(runNumber);
    
    Npassed_In += weightII;
    
    double m0=9999, m12=9999, mc=9999;
     
    if(!_isData && isSMS){

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

    //HLT
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();

    if(!_isData)passedHLT = 1;
    if(!passedHLT)continue;

    Npassed_HLT += weightII;

    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    Npassed_PV += weightII;

    vector<TLorentzVector> PFJet;

    int iTightMuons[10], tightMuonCounter=0;
    int iLooseMuons[10], looseMuonCounter=0;
    int i95Eles[10], Ele95Counter=0;
    int i80Eles[10], Ele80Counter=0;

    double etaLept;
    double phiLept;
    double ptLept;

    //Muons
    SLT=0.;
    for(int j=0; j<nMuon;j++){
      double muonPt=sqrt(pxMuon[j]*pxMuon[j]+pyMuon[j]*pyMuon[j]);      
      if(isTightMuon(j) && muonPt>25.){
	iTightMuons[tightMuonCounter]=j;
	etaLept=etaMuon[j];
	phiLept=phiMuon[j];
	ptLept=muonPt;
	tightMuonCounter++;
      }
    }

    double IECAL; 
    double IHCAL;
    double ITRK;

    for(int j=0; j<nMuon;j++){
      double muonPt=sqrt(pxMuon[j]*pxMuon[j]+pyMuon[j]*pyMuon[j]);      
      if(isLooseMuon(j) && muonPt>20.){
	iLooseMuons[looseMuonCounter]=j;
	IECAL = emEt03Muon[j];
	IHCAL = hadEt03Muon[j];
	ITRK  = sumPt03Muon[j];
	MuIso[looseMuonCounter]=((IECAL+IHCAL+ITRK)/muonPt);
	looseMuonCounter++;
	SLT+=muonPt;
      }
    }

    //Electrons
     for(int k=0; k< nEle; k++){
        double ElePt=sqrt(pxEle[k]*pxEle[k]+pyEle[k]*pyEle[k]);
	if(is80Electron(k)&& ElePt>20.){
	  i80Eles[Ele80Counter]=k;
	  etaLept=etaEle[k];
	  phiLept=phiEle[k];
	  ptLept=ElePt;
	  Ele80Counter++;
	}
     }

     int iSC;
     double ET;
     double trackISO;
     double ECALISO;
     double HCALISO;

     for(int k=0; k< nEle; k++){
        double ElePt=sqrt(pxEle[k]*pxEle[k]+pyEle[k]*pyEle[k]);
	if(is95Electron(k) && ElePt>20.){
	  i95Eles[Ele95Counter]=k;
	  iSC=superClusterIndexEle[k];
	  ET=energySC[iSC]/cosh(etaSC[iSC]);
	  trackISO = (dr03TkSumPtEle[k])/ET;
	  ECALISO = (dr03EcalRecHitSumEtEle[k])/ET;
	  HCALISO = (dr03HcalTowerSumEtEle[k])/ET;
	  EleTrackIso[Ele95Counter]=trackISO;
	  EleECALIso[Ele95Counter]=ECALISO;
	  EleHCALIso[Ele95Counter]=HCALISO;
	  Ele95Counter++;
	  SLT+=ElePt;
	}
     }

     //Jets                                                                                                                             
     HT=0.;
     vector <int> iPFJet;
     for(int i=0; i< nAK5PFPUcorrJet; i++) {
       TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);
       if(myJet.Pt()>30. && fabs(myJet.Eta())< 2.4) {
	 HT+=myJet.Pt();
	 PFJet.push_back(myJet);
	 iPFJet.push_back(i);
       }
     }

     double DREleJet;
     vector<int> BadPFEleJet;
     int counter, double_check;
     for(int l=0; l<Ele80Counter; l++){
       for(int n=0; n < PFJet.size(); n++){
	 TLorentzVector myJet = PFJet.at(n);
	 DREleJet=sqrt((etaEle[i80Eles[l]]-myJet.Eta())*(etaEle[i80Eles[l]]-myJet.Eta())+(phiEle[i80Eles[l]]-myJet.Phi())*(phiEle[i80Eles[l]]-myJet.Phi()));
	 double_check=0;
	 counter=0;
	 if(DREleJet < 0.5){
	   if(BadPFEleJet.size()!=0){
	     for(int f=0; f<BadPFEleJet.size(); f++){
	       if(n!=BadPFEleJet.at(f))counter++;
	     }
	     if(counter == BadPFEleJet.size())double_check=1;
	     if(double_check!=0)BadPFEleJet.push_back(n);
	   }else{
	     BadPFEleJet.push_back(n);
	   }
	 }
       }
     }

     double DRMuJet;
     vector<int> BadPFMuJet;
     int counterII, double_checkII, counter_Ele, double_checkEle;
     for(int l=0; l<looseMuonCounter; l++){
       for(int n=0; n < PFJet.size(); n++){
	 TLorentzVector myJet = PFJet.at(n);
	 DRMuJet=sqrt((etaMuon[iLooseMuons[l]]-myJet.Eta())*(etaMuon[iLooseMuons[l]]-myJet.Eta())+(phiMuon[iLooseMuons[l]]-myJet.Phi())*(phiMuon[iLooseMuons[l]]-myJet.Phi()));
	 double_checkII=0;
	 counterII=0;
	 counter_Ele=0;
	 double_checkEle=0;
	 if(DRMuJet < 0.5){
	   for(int z=0; z<BadPFEleJet.size(); z++){
             if(n!=BadPFEleJet.at(z))counter_Ele++;
           }
	   if(counter_Ele==BadPFEleJet.size())double_checkEle=1;
	   if(BadPFMuJet.size()!=0){
	   for(int f=0; f<BadPFMuJet.size(); f++){
	     if(n!=BadPFMuJet.at(f))counterII++;
	   }
	   if(counterII == BadPFMuJet.size())double_checkII=1;
	   if(double_checkII!=0 && double_checkEle!=0)BadPFMuJet.push_back(n);
	   }else{
	     if(double_checkEle!=0)BadPFMuJet.push_back(n);
	   }
	 }
       }
     }

     //Subtract fake jets pT from HT
     int c;
     for(int z=0; z< BadPFEleJet.size(); z++){
       c=BadPFEleJet.at(z);
       TLorentzVector myJet(pxAK5PFPUcorrJet[c], pyAK5PFPUcorrJet[c], pzAK5PFPUcorrJet[c], energyAK5PFPUcorrJet[c]);
       HT-=myJet.Pt();
     }
     for(int a=0; a < BadPFMuJet.size(); a++){
       c=BadPFMuJet.at(a);
       TLorentzVector myJet(pxAK5PFPUcorrJet[c], pyAK5PFPUcorrJet[c], pzAK5PFPUcorrJet[c], energyAK5PFPUcorrJet[c]);
       HT-=myJet.Pt();
     }

     //Create arrays with the b-discriminators of the jets, in descending order
     int n;
     vector <float> CSV;
     vector <float> TCHE;
     vector <float> SSV;
     vector <float> JP;
     vector <float> JBP;
     nBTagJets=0;
     for(int b=0; b< iPFJet.size(); b++){
       n=iPFJet.at(b);
       int b_counter=0;
       for(int f=0; f<BadPFMuJet.size(); f++){
	 if(n!=BadPFMuJet.at(f))b_counter++;
       }
       for(int d=0; d<BadPFEleJet.size(); d++){
	 if(n!=BadPFEleJet.at(d))b_counter++;
       }
       if(b_counter==(BadPFEleJet.size()+BadPFMuJet.size())){
	 if(trackCountingHighEffBJetTagsAK5PFPUcorrJet[n] > 3.3)nBTagJets++;
	 CSV.push_back(combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet[n]);
	 TCHE.push_back(trackCountingHighEffBJetTagsAK5PFPUcorrJet[n]);
	 SSV.push_back(simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[n]);
	 JP.push_back(jetProbabilityBJetTagsAK5PFPUcorrJet[n]);
	 JBP.push_back(jetBProbabilityBJetTagsAK5PFPUcorrJet[n]);
       }
     }

     BubbleSort(CSV);
     BubbleSort(TCHE);
     BubbleSort(SSV);
     BubbleSort(JP);
     BubbleSort(JBP);
     //Put the sum of the four highest discr. values into the tree
     CSV_TopFour=-9999;
     TCHE_TopFour=-9999;
     SSV_TopFour=-9999;
     JP_TopFour=-9999;
     JBP_TopFour=-9999;
     if(CSV.size() >= 4)CSV_TopFour=(CSV.at(0)+CSV.at(1)+CSV.at(2)+CSV.at(3));
     if(TCHE.size()>= 4)TCHE_TopFour=(TCHE.at(0)+TCHE.at(1)+TCHE.at(2)+TCHE.at(3));
     if(SSV.size() >= 4)SSV_TopFour=(SSV.at(0)+SSV.at(1)+SSV.at(2)+SSV.at(3));
     if(JP.size()>= 4)JP_TopFour=(JP.at(0)+JP.at(1)+JP.at(2)+JP.at(3));
     if(JBP.size()>= 4)JBP_TopFour=(JBP.at(0)+JBP.at(1)+JBP.at(2)+JBP.at(3));


     //Fill tree for SMS efficiencies
     if(isSMS)FullSMSTree->Fill();

     //Selection

     //MET
     if(sqrt(pxPFMet[0]*pxPFMet[0]+pyPFMet[0]*pyPFMet[0])<150.)continue;
     Npassed_MET+=weightII;
     //ST
     ST=MET+SLT+HT;
     if(ST<1000.)continue;
     Npassed_ST+=weightII;
     //4Jets
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size())<4)continue;
     NpassedPF_4Jet+=weightII;
     //3b
     if(nBTagJets < 3)continue;
     Npassed_3b+=weightII;

     // dummy values                                          
     pTPFHem1 = -9999.;
     etaPFHem1 = -9999.;
     phiPFHem1 = -9999.;
     pTPFHem2 = -9999.;
     etaPFHem2 = -9999.;
     phiPFHem2 = -9999.;
     massPFHem = -99.;
     PFgoodR = -1;
     PFR = -99999.;
     PFMR = -99999.;

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
     
     if(fabs(beta)< 0.99){
       NpassedPF_beta += weightII;
     
       if(fabs(PFHem1.DeltaPhi(PFHem2)) < 2.8){
	 NpassedPF_DeltaPhi += weightII;

	 TVector3 MET(pxMet[0], pyMet[0], 0.);
	 double MT = CalcMTR(PFHem1, PFHem2, MET);
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
	 
	 // fill the R and hem part of the output tree
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

     // fill output tree
     n_PV=nPV;
     MET=sqrt(pxPFMet[0]*pxPFMet[0]+pyPFMet[0]*pyPFMet[0]);
     etaMETLept=etaPFMet[0]-etaLept;
     phiMETLept=phiPFMet[0]-phiLept;
     PFMT=sqrt(2*MET*ptLept*(1-cos(phiMETLept)));
     nJets=nAK5PFJet;
     nEleJets=0.;
     N80Eles = Ele80Counter;
     N95Eles = Ele95Counter;
     NLooseMuons = looseMuonCounter;
     NTightMuons = tightMuonCounter;
  
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
  effTree->Branch("Npassed_HLT",      &Npassed_HLT,      "Npassed_HLT/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("NpassedPF_4Jet",      &NpassedPF_4Jet,      "NpassedPF_4Jet/D");
  effTree->Branch("Npassed_NoLept",      &Npassed_NoLept,      "Npassed_NoLept/D");
  effTree->Branch("NpassedPF_beta",      &NpassedPF_beta,      "NpassedPF_beta/D");
  effTree->Branch("NpassedPF_DeltaPhi",      &NpassedPF_DeltaPhi,      "NpassedPF_DeltaPhi/D");
  effTree->Branch("Npassed_MET",      &Npassed_MET,      "Npassed_MET/D");
  effTree->Branch("Npassed_ST",      &Npassed_ST,      "Npassed_ST/D");
  effTree->Branch("Npassed_3b",      &Npassed_3b,      "Npassed_3b/D");

  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  effTree->Write();
  if(isSMS)FullSMSTree->Write();
  file->Close();
}
  
vector<TLorentzVector> SUSYMultiB::CombineJets(vector<TLorentzVector> myjets){

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

void SUSYMultiB::BubbleSort(vector <float> &num)
{
  int i, j, flag = 1;    // set flag to 1 to start first pass
  float temp;             // holding variable
  int numLength = num.size(); 
  for(i = 1; (i <= numLength) && flag; i++)
    {
      flag = 0;
      for (j=0; j < (numLength -1); j++)
	{
	  if (num[j+1] > num[j])      // ascending order simply changes to <
	    { 
	      temp = num[j];             // swap elements
	      num[j] = num[j+1];
	      num[j+1] = temp;
	      flag = 1;               // indicates that a swap occurred.
	    }
	}
    }
  return;   //arrays are passed to functions by address; nothing is returned
}
