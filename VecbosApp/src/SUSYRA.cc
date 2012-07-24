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
#include "RA4Selector.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "SUSYRA.hh"

SUSYRA::SUSYRA(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;

}

SUSYRA::SUSYRA(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;

  //To read good run list!                                                                                                                                                                                                                   
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

SUSYRA::~SUSYRA() {}

void SUSYRA::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

bool SUSYRA::isGlobalMuon(int iMu){
  bool ret = false;
  Utils anaUtils;
  bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
  bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
  bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);

  if(isMuGlobal && isMuGlobalPrompt && isMuTracker)ret=true;

  return ret;
}

bool SUSYRA::isVetoMuon(int iMu){

  bool ret=false;
  bool isGlobal = isGlobalMuon(iMu);
  if(isGlobal){
     double pt = sqrt(pxMuon[iMu]*pxMuon[iMu]+pyMuon[iMu]*pyMuon[iMu]);
    double IECAL = emEt03Muon[iMu];
    double IHCAL = hadEt03Muon[iMu];
    double ITRK  = sumPt03Muon[iMu];
    
    if((IECAL+IHCAL+ITRK)/pt < 0.2){                                                                                                                                                                                          
      ret = true;
    } 
  }

  return ret;
  
}

bool SUSYRA::isRA4Muon(int iMu){

  bool ret = false;
  Utils anaUtils;
  bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
  bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);

  if(isMuGlobalPrompt && isMuTracker){

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

	  if((IECAL+IHCAL+ITRK)/pt < 0.10){                                                                                                                                                                                          
	    ret = true;
	  }                                                                                                                                                                                                                              
        }
      }
    }
  }

  return ret;
}

bool SUSYRA::is80Electron(int iEle){

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
    if(trackISO > 0.09) return false;
    if(ECALISO > 0.07) return false;
    if(HCALISO > 0.10) return false;

    //ID WP80                                                                                                                                                
    if(fabs(dphi) > 0.06) return false;
    if(fabs(deta) > 0.004) return false;
    if(HoE > 0.04) return false;
    if(sigietaieta > 0.01) return false;
  } else {
    //ISO WP80                                                                                                                                               
    if(trackISO > 0.04) return false;
    if(ECALISO > 0.05) return false;
    if(HCALISO > 0.025) return false;

    //ID WP80                                                                                                                                                
    if(fabs(dphi) > 0.03) return false;
    //if(fabs(deta) > 0.007) return false;                                                                                                                   
    if(HoE > 0.025) return false;
    if(sigietaieta > 0.03) return false;
  }

  return true;
}


bool SUSYRA::is95Electron(int iEle){

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


void SUSYRA::Loop(string outFileName, int start, int stop) {
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

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");

  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");

  // PF block
  outTree->Branch("MET", &MET, "MET/D");
  outTree->Branch("phiMETLept", &phiMETLept, "phiMETLept/D");
  outTree->Branch("etaMETLept", &etaMETLept, "etaMETLept/D");
  outTree->Branch("ST", &ST, "ST/D");
  outTree->Branch("SLT", &SLT, "SLT/D");
  outTree->Branch("HT", &HT, "HT/D");
  outTree->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
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

  // prepare vectors for efficiency
  double Npassed_In = 0;
  //HLT
  double Npassed_HLT = 0;
  
  double Npassed_PV = 0;
  
  //Jets
  double NpassedPF_4Jet = 0;
  //Leptons
  double Npassed_TightMu=0;
  double Npassed_Ele=0;
 
  double Npassed_MET=0;

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
    
    if(_isData)reloadTriggerMask(true);
    //    else if (jentry == 0) reloadTriggerMask();

    Npassed_In += weight;

    //HLT
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();

    if(!_isData)passedHLT = 1;
    if(!passedHLT)continue;

    Npassed_HLT += weight;

    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    Npassed_PV += weight;

    vector<TLorentzVector> PFJet;

    int iTightMuons[10], tightMuonCounter=0;
    int iLooseMuons[10], looseMuonCounter=0;
    int i95Eles[10], Ele95Counter=0;
    int i80Eles[10], Ele80Counter=0;

    double etaLept;
    double phiLept;

    //Muons
    SLT=0.;
    for(int j=0; j<nMuon;j++){
      double muonPt=sqrt(pxMuon[j]*pxMuon[j]+pyMuon[j]*pyMuon[j]);      
      if(isRA4Muon(j) && muonPt>20. && fabs(etaMuon[j])< 2.1){
	iTightMuons[tightMuonCounter]=j;
	etaLept=etaMuon[j];
	phiLept=phiMuon[j];
	tightMuonCounter++;
      }
    }

    for(int j=0; j<nMuon;j++){
      double muonPt=sqrt(pxMuon[j]*pxMuon[j]+pyMuon[j]*pyMuon[j]);      
      if(isVetoMuon(j) && muonPt>10. && fabs(etaMuon[j])<2.5){
	iLooseMuons[looseMuonCounter]=j;
	looseMuonCounter++;
	SLT+=muonPt;
      }
    }

    //Electrons
     for(int k=0; k< nEle; k++){
        double ElePt=sqrt(pxEle[k]*pxEle[k]+pyEle[k]*pyEle[k]);
	if(is80Electron(k)&& ElePt>25.){
	  i80Eles[Ele80Counter]=k;
	  etaLept=etaEle[k];
	  phiLept=phiEle[k];
	  Ele80Counter++;
	}
     }

     for(int k=0; k< nEle; k++){
        double ElePt=sqrt(pxEle[k]*pxEle[k]+pyEle[k]*pyEle[k]);
	if(is95Electron(k) && ElePt>15.){
	  i95Eles[Ele95Counter]=k;
	  Ele95Counter++;
	  SLT+=ElePt;
	}
     }

     //Jets                                                                                                                                                          
     HT=0.;
     for(int i=0; i< nAK5PFJet; i++) {
       TLorentzVector myJet(pxAK5PFJet[i], pyAK5PFJet[i], pzAK5PFJet[i], energyAK5PFJet[i]);
       HT+=myJet.Pt();
       if(myJet.Pt()>30. && fabs(myJet.Eta())< 2.4) {
	 PFJet.push_back(myJet);
       }
     }

     double DREleJet;
     vector<int> BadPFEleJet;
     int counter, double_check;
     for(int l=0; l<Ele95Counter; l++){
       for(int n=0; n < PFJet.size(); n++){
	 TLorentzVector myJet = PFJet.at(n);
	 DREleJet=sqrt((etaEle[i95Eles[l]]-myJet.Eta())*(etaEle[i95Eles[l]]-myJet.Eta())+(phiEle[i95Eles[l]]-myJet.Phi())*(phiEle[i95Eles[l]]-myJet.Phi()));
	 double_check=0;
	 counter=0;
	 if(DREleJet < 0.3){
	   for(int f=0; f<BadPFEleJet.size(); f++){
	     if(n!=BadPFEleJet.at(f))counter++;
	   }
	   if(counter == BadPFEleJet.size())double_check=1;
	 }
	 if(double_check!=0 || BadPFEleJet.size()==0)BadPFEleJet.push_back(n);
       }
     }

     double DRMuJet;
     vector<int> BadPFMuJet;
     int counterII, double_checkII;
     for(int l=0; l<looseMuonCounter; l++){
       for(int n=0; n < PFJet.size(); n++){
	 TLorentzVector myJet = PFJet.at(n);
	 DRMuJet=sqrt((etaMuon[iLooseMuons[l]]-myJet.Eta())*(etaMuon[iLooseMuons[l]]-myJet.Eta())+(phiMuon[iLooseMuons[l]]-myJet.Phi())*(phiMuon[iLooseMuons[l]]-myJet.Phi())) \
	   ;
	 double_checkII=0;
	 counterII=0;
	 if(DRMuJet < 0.1){
	   for(int f=0; f<BadPFMuJet.size(); f++){
	     if(n!=BadPFMuJet.at(f))counterII++;
	   }
	   if(counterII == BadPFMuJet.size())double_checkII=1;
	 }
	 if(double_checkII!=0 || BadPFMuJet.size()==0)BadPFMuJet.push_back(n);
       }
     }


     //Selection

#if RA4Selector == 1
     if(tightMuonCounter!=1)continue;
     if(looseMuonCounter!=1)continue;
     if(Ele80Counter!=0)continue;
     if(Ele95Counter!=0)continue;
     Npassed_TightMu+=weight;
#endif     

#if RA4Selector == 2
     if(tightMuonCounter!=0)continue;
     if(looseMuonCounter!=0)continue;
     if(Ele80Counter!=1)continue;
     if(Ele95Counter!=1)continue;
     Npassed_Ele+=weight;
#endif

#if isJustScaling == false
     if(PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size() < 4) continue;
     NpassedPF_4Jet += weight;
#endif

     //MET
     if(sqrt(pxPFMet[0]*pxPFMet[0]+pyPFMet[0]*pyPFMet[0])< 100.)continue;
     Npassed_MET+=weight;

     // fill output tree
     MET=sqrt(pxPFMet[0]*pxPFMet[0]+pyPFMet[0]*pyPFMet[0]);
     etaMETLept=etaPFMet[0]-etaLept;
     phiMETLept=phiPFMet[0]-phiLept;
     ST=MET+SLT+HT;
     nJets=nAK5PFJet;
     nEleJets=0.;
     JetMultiplicity=PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size();
     EleFakeJetMultiplicity=BadPFEleJet.size();
     MuFakeJetMultiplicity=BadPFMuJet.size();
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
  effTree->Branch("Npassed_TightMu",      &Npassed_TightMu,      "Npassed_TightMu/D");
  effTree->Branch("Npassed_Ele",      &Npassed_Ele,      "Npassed_Ele/D");
  effTree->Branch("Npassed_MET",      &Npassed_MET,      "Npassed_MET/D");

  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  effTree->Write();
  outTree->Write();
  file->Close();
}
  
vector<TLorentzVector> SUSYRA::CombineJets(vector<TLorentzVector> myjets){

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
