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
#include "RazorDiPhoton.hh"

RazorDiPhoton::RazorDiPhoton(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
}

RazorDiPhoton::RazorDiPhoton(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

RazorDiPhoton::~RazorDiPhoton() {}

void RazorDiPhoton::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorDiPhoton::SetWeight(double weight) {
  _weight = weight;
}

void RazorDiPhoton::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  // prescaled single-photon triggers
  int HLT_Photon50_CaloIdVL;
  int HLT_Photon50_CaloIdVL_IsoL;
  // double-photon triggers
  int HLT_DoublePhoton33;
  int HLT_DoubleEle33_CaloIdL;
  int HLT_DoublePhoton33_HEVT;
  int HLT_DoublePhoton50;
  // photon+Had trigger
  int HLT_Photon70_CaloIdL_MHT30;
  int HLT_Photon70_CaloIdL_HT200;
  // razor diphoton
  int HLT_DoublePhoton40_MR150;
  int HLT_DoublePhoton40_R014_MR150;
  // razor photon+jets
  int HLT_Photon40_R005_MR150;
  int HLT_Photon40_R014_MR450;
  int HLT_Photon40_R014_MR500;
  int HLT_Photon40_R017_MR500;
  int HLT_Photon40_R020_MR300;
  int HLT_Photon40_R020_MR350;
  int HLT_Photon40_R023_MR350;
  int HLT_Photon40_R025_MR200;
  int HLT_Photon40_R025_MR250;
  int HLT_Photon40_R029_MR250;
  int HLT_Photon40_R038_MR150;
  int HLT_Photon40_R038_MR200;
  int HLT_Photon40_R042_MR200;

  // photon block
  double pTPh1;
  double etaPh1;
  double phiPh1;
  int    noPhoton1;
  int    isEMPh1;
  int    isLoosePh1;
  int    isTightPh1;
  double pTPh2;
  double etaPh2;
  double phiPh2;
  int    noPhoton2;
  int    isEMPh2;
  int    isLoosePh2;
  int    isTightPh2;

  // PF  block
  int    passedPF;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double PFR;
  double PFMR;

  // calo block
  int    passedCalo;
  double pTCaloHem1;
  double etaCaloHem1;
  double phiCaloHem1;
  double pTCaloHem2;
  double etaCaloHem2;
  double phiCaloHem2;
  double CaloR;
  double CaloMR;

  // general event info
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double W;

  // Muon Block
  double pTMu1;
  double etaMu1;
  double phiMu1;
  double chargeMu1;
  double pTMu2;
  double etaMu2;
  double phiMu2;
  double chargeMu2;

  // Electron Block
  double pTEle1;
  double etaEle1;
  double phiEle1;
  double chargeEle1;
  double pTEle2;
  double etaEle2;
  double phiEle2;
  double chargeEle2;


  // prepare the output tree
  TTree* outTree[6]; 
  outTree[0] = new TTree("outTreeHad", "outTreeHad");
  outTree[1] = new TTree("outTreeEle", "outTreeEle");
  outTree[2] = new TTree("outTreeMu", "outTreeMu");
  outTree[3] = new TTree("outTreeMuMu", "outTreeMuMu");
  outTree[4] = new TTree("outTreeMuEle", "outTreeMuEle");
  outTree[5] = new TTree("outTreeEleEle", "outTreeEleEle");

  for(int i=0; i<6; i++) {
    outTree[i]->Branch("HLT_Photon50_CaloIdVL", &HLT_Photon50_CaloIdVL, "HLT_Photon50_CaloIdVL/I");
    outTree[i]->Branch("HLT_Photon50_CaloIdVL_IsoL", &HLT_Photon50_CaloIdVL_IsoL, "HLT_Photon50_CaloIdVL_IsoL/I");
    outTree[i]->Branch("HLT_DoublePhoton33", &HLT_DoublePhoton33, "HLT_DoublePhoton33/I");
    outTree[i]->Branch("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL, "HLT_DoubleEle33_CaloIdL/I");
    outTree[i]->Branch("HLT_DoublePhoton33_HEVT", &HLT_DoublePhoton33_HEVT, "HLT_DoublePhoton33_HEVT/I");
    outTree[i]->Branch("HLT_DoublePhoton50", &HLT_DoublePhoton50, "HLT_DoublePhoton50/I");
    outTree[i]->Branch("HLT_Photon70_CaloIdL_MHT30", &HLT_Photon70_CaloIdL_MHT30, "HLT_Photon70_CaloIdL_MHT30/I");
    outTree[i]->Branch("HLT_Photon70_CaloIdL_HT200", &HLT_Photon70_CaloIdL_HT200, "HLT_Photon70_CaloIdL_HT200/I");
    outTree[i]->Branch("HLT_DoublePhoton40_MR150", &HLT_DoublePhoton40_MR150, "HLT_DoublePhoton40_MR150/I");
    outTree[i]->Branch("HLT_DoublePhoton40_R014_MR150", &HLT_DoublePhoton40_R014_MR150, "HLT_DoublePhoton40_R014_MR150/I");
    outTree[i]->Branch("HLT_Photon40_R005_MR150", &HLT_Photon40_R005_MR150, "HLT_Photon40_R005_MR150/I");
    outTree[i]->Branch("HLT_Photon40_R014_MR450", &HLT_Photon40_R014_MR450, "HLT_Photon40_R014_MR450/I");
    outTree[i]->Branch("HLT_Photon40_R014_MR500", &HLT_Photon40_R014_MR500, "HLT_Photon40_R014_MR500/I");
    outTree[i]->Branch("HLT_Photon40_R017_MR500", &HLT_Photon40_R017_MR500, "HLT_Photon40_R017_MR500/I");
    outTree[i]->Branch("HLT_Photon40_R020_MR300", &HLT_Photon40_R020_MR300, "HLT_Photon40_R020_MR300/I");
    outTree[i]->Branch("HLT_Photon40_R020_MR350", &HLT_Photon40_R020_MR350, "HLT_Photon40_R020_MR350/I");
    outTree[i]->Branch("HLT_Photon40_R023_MR350", &HLT_Photon40_R023_MR350, "HLT_Photon40_R023_MR350/I");
    outTree[i]->Branch("HLT_Photon40_R025_MR200", &HLT_Photon40_R025_MR200, "HLT_Photon40_R025_MR200/I");
    outTree[i]->Branch("HLT_Photon40_R025_MR250", &HLT_Photon40_R025_MR250, "HLT_Photon40_R025_MR250/I");
    outTree[i]->Branch("HLT_Photon40_R029_MR250", &HLT_Photon40_R029_MR250, "HLT_Photon40_R029_MR250/I");
    outTree[i]->Branch("HLT_Photon40_R038_MR150", &HLT_Photon40_R038_MR150, "HLT_Photon40_R038_MR150/I");
    outTree[i]->Branch("HLT_Photon40_R038_MR200", &HLT_Photon40_R038_MR200, "HLT_Photon40_R038_MR200/I");
    outTree[i]->Branch("HLT_Photon40_R042_MR200", &HLT_Photon40_R042_MR200, "HLT_Photon40_R042_MR200/I");

    outTree[i]->Branch("run", &run, "run/D");
    outTree[i]->Branch("evNum", &evNum, "evNum/D");
    outTree[i]->Branch("bx", &bx, "bx/D");
    outTree[i]->Branch("ls", &ls, "ls/D");
    outTree[i]->Branch("orbit", &orbit, "orbit/D");
    outTree[i]->Branch("W", &W, "W/D");
    
    outTree[i]->Branch("pTPh1", &pTPh1, "pTPh1/D");
    outTree[i]->Branch("etaPh1", &etaPh1, "etaPh1/D");
    outTree[i]->Branch("phiPh1", &phiPh1, "phiPh1/D");
    outTree[i]->Branch("noPhoton1", &noPhoton1, "noPhoton1/I");
    outTree[i]->Branch("isEMPh1", &isEMPh1, "isEMPh1/I");
    outTree[i]->Branch("isLoosePh1", &isLoosePh1, "isLoosePh1/I");
    outTree[i]->Branch("isTightPh1", &isTightPh1, "isTightPh1/I");
    outTree[i]->Branch("pTPh2", &pTPh2, "pTPh2/D");
    outTree[i]->Branch("etaPh2", &etaPh2, "etaPh2/D");
    outTree[i]->Branch("phiPh2", &phiPh2, "phiPh2/D");
    outTree[i]->Branch("noPhoton2", &noPhoton2, "noPhoton2/I");
    outTree[i]->Branch("isEMPh2", &isEMPh2, "isEMPh2/I");
    outTree[i]->Branch("isLoosePh2", &isLoosePh2, "isLoosePh2/I");
    outTree[i]->Branch("isTightPh2", &isTightPh2, "isTightPh2/I");
    
    // PF block
    outTree[i]->Branch("passedPF", &passedPF, "passedPF/I");
    outTree[i]->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
    outTree[i]->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
    outTree[i]->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
    outTree[i]->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
    outTree[i]->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
    outTree[i]->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
    outTree[i]->Branch("PFR", &PFR, "PFR/D");
    outTree[i]->Branch("PFMR", &PFMR, "PFMR/D");
    
    // Calo block
    outTree[i]->Branch("passedCalo", &passedCalo, "passedCalo/I");
    outTree[i]->Branch("pTCaloHem1", &pTCaloHem1, "pTCaloHem1/D");
    outTree[i]->Branch("etaCaloHem1", &etaCaloHem1, "etaCaloHem1/D");
    outTree[i]->Branch("phiCaloHem1", &phiCaloHem1, "phiCaloHem1/D");
    outTree[i]->Branch("pTCaloHem2", &pTCaloHem2, "pTCaloHem2/D");
    outTree[i]->Branch("etaCaloHem2", &etaCaloHem2, "etaCaloHem2/D");
    outTree[i]->Branch("phiCaloHem2", &phiCaloHem2, "phiCaloHem2/D");
    outTree[i]->Branch("CaloR", &CaloR, "CaloR/D");
    outTree[i]->Branch("CaloMR", &CaloMR, "CaloMR/D");
    
    // First muon for Mu,  DiMu and MuEle
    if(i ==2 || i == 3 || i == 4) {
      outTree[i]->Branch("pTMu1", &pTMu1, "pTMu1/D");
      outTree[i]->Branch("etaMu1", &etaMu1, "etaMu1/D");
      outTree[i]->Branch("phiMu1", &phiMu1, "phiMu1/D");
      outTree[i]->Branch("chargeMu1", &chargeMu1, "chargeMu11/D");
    }
    // Second MU
    if(i == 3) {
      outTree[i]->Branch("pTMu2", &pTMu2, "pTMu2/D");
      outTree[i]->Branch("etaMu2", &etaMu2, "etaMu2/D");
      outTree[i]->Branch("phiMu2", &phiMu2, "phiMu2/D");
      outTree[i]->Branch("chargeMu2", &chargeMu2, "chargeMu21/D");
    }
    // First Ele for Ele, DiEle and MuEle
    if(i ==1 || i == 4 || i == 5) {
      outTree[i]->Branch("pTEle1", &pTEle1, "pTEle1/D");
      outTree[i]->Branch("etaEle1", &etaEle1, "etaEle1/D");
      outTree[i]->Branch("phiEle1", &phiEle1, "phiEle1/D");
      outTree[i]->Branch("chargeEle1", &chargeEle1, "chargeEle11/D");
    }
    // Second Ele
    if(i == 5) {
      outTree[i]->Branch("pTEle2", &pTEle2, "pTEle2/D");
      outTree[i]->Branch("etaEle2", &etaEle2, "etaEle2/D");
      outTree[i]->Branch("phiEle2", &phiEle2, "phiEle2/D");
      outTree[i]->Branch("chargeEle2", &chargeEle2, "chargeEle21/D");
    }
  }


  // prepare vectors for efficiency (variables for effTree)
  double Npassed_In = 0;
  double Npassed_HLT = 0; 
  double Npassed_PV = 0;
  double Npassed_2Ph = 0;
  double Npassed_PFHem = 0;
  double Npassed_CaloHem = 0;

  //  double _weight = 1.;
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
    
    Npassed_In += _weight;

    // HLT 
    //   vector<int> myHLT = getHLTOutput();
    //    HLTbit = 0;
    //    if(myHLT.size() == 1) HLTbit = myHLT[0];
    //    HLTbitEmulator = HTtrigger(150.);
    Npassed_HLT += _weight;

    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    Npassed_PV += _weight;

    int goodPFevent = true;
    int goodCaloevent = true;

    // HCAL FLAGS
    if(!eventPassHcalFilter()) continue;

    // Jet selection 
    vector<TLorentzVector> PFPUcorrJet;    
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);   
      if(myJet.Pt()>40. && fabs(myJet.Eta())< 2.4) 
	//	if(myJet.DeltaR(PFPhoton[iPh1]) > 0.35 && myJet.DeltaR(PFPhoton[iPh2]) > 0.35)
	PFPUcorrJet.push_back(myJet);
    }
    
    vector<TLorentzVector> CaloJet;    
    for(int i = 0; i < nAK5Jet; i++){
      TLorentzVector jet;
      double px = pxAK5Jet[i];
      double py = pyAK5Jet[i];
      double pz = pzAK5Jet[i];
      double E = sqrt(px*px+py*py+pz*pz);
      jet.SetPxPyPzE(px,py,pz,E);
      
      if(jet.Pt() < 40. || fabs(jet.Eta()) >= 3.0) continue; //jet kinematic cuts
      CaloJet.push_back(jet);
      
      // if(nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98){ //fails jet id
      // 	goodCaloevent = false;
      // }
      // if(fabs(jet.Eta()) < 2.55 && emFracAK5Jet[i] <= 0.01){ //fails jet id
      // 	goodCaloevent = false;
      // }
      // if(fabs(jet.Eta()) >= 2.55){ //fails jet ID
      // 	if(jet.Pt() > 80. && emFracAK5Jet[i] >= 1.){
      // 	  goodCaloevent = false;	
      // 	}
      // 	if(emFracAK5Jet[i] <= -0.9){
      // 	  goodCaloevent = false;
      // 	}
      //      }
    }

    vector<TLorentzVector> PFPhoton;
    vector<bool> isEM;
    vector<bool> isLoose;
    vector<bool> isTight;
    vector<bool> notPhoton;
    for(int i=0; i< nPho; i++) {
      // ET cut
      double pT = energyPho[i]*sin(thetaPho[i]);
      if(pT  >= 40.) {

	/*	
	  bool isEM_i = true;
	  bool isLoose_i = true;
	  bool isTight_i = true;
	  PhotonIdBarrel(i,isEM_i,isLoose_i,isTight_i);  ////// rearrange
	*/

	bool isEM_i = false;
	bool isLoose_i = false;
	bool isTight_i = false;
	
	if ( PhotonIdBarrelisEM(i)){ isEM_i = true;} // == true
	if ( PhotonIdBarrelisLoose(i)){ isLoose_i = true;} // == true
	if ( PhotonIdBarrelisTight(i)){ isTight_i = true;} // == true

	//if isTight is OK isEM and isLoose too//

	//cout<<" isEM = "<<isEM_i<<"//// isLoose= "<<isLoose_i<<"//// isTight= "<<isTight_i<<endl; 

	bool notPhoton_i = AntiPhotonIdBarrel(i);
	// e2/e9
	int iSC = superClusterIndexPho[i];
	double e2OVERe9 = (eMaxSC[iSC]+e2ndSC[iSC])/e3x3SC[iSC];
	//	if(e2OVERe9 <= 0.95) {
	// one jet OR timing [for no-jet events]
	//	  double time = fabs(timeSC[superClusterIndexPho[i]]);
	//	  if(PFPUcorrJet.size()+CaloJet.size() > 0 || time < 3.) {
	TLorentzVector myPhoton;
	myPhoton.SetPtEtaPhiE(energyPho[i]*sin(thetaPho[i]), etaPho[i], phiPho[i], energyPho[i]);
	PFPhoton.push_back(myPhoton);
	isLoose.push_back(isLoose_i);
	isTight.push_back(isTight_i);
	isEM.push_back(isEM_i);
	notPhoton.push_back(notPhoton_i);
	//	  }
	//	}
      }
    }

    if(PFPhoton.size() < 2 ) continue; // keep at least 2 objects

    // select two highest-pT photons
    // and select what kind of photons we are : isEM[i] or isLoose[i] or isTight[i] 
    int iPh1 = -99;
    int iPh2 = -99;
    double phMax1 = 0.;
    double phMax2 = 0.;
    for(int i=0; i<PFPhoton.size(); i++) {
      if(PFPhoton[i].Pt() > phMax1 && isEM[i]) {
	phMax1 = PFPhoton[i].Pt();
	iPh1 = i;
      }
    }

    for(int i=0; i<PFPhoton.size(); i++) {
      if(PFPhoton[i].Pt() > phMax2 && i != iPh1 && isEM[i]) {
	phMax2 = PFPhoton[i].Pt();
	iPh2 = i;
      }
    }

    //    if(iPh1 == -99 && iPh2 == -99) continue;
    Npassed_2Ph += _weight;

    // if only one found, select the highest-pT fake
    // the first is always the good one
    if(iPh1 == -99) {
      iPh1 = iPh2;
      iPh2 = -99;
    }

    // if no good photon was found 
    // use the highest-pT fake
    if(iPh1 == -99) {
      phMax1 = 0.;
      for(int i=0; i<PFPhoton.size(); i++) {
	if(PFPhoton[i].Pt() > phMax1) {
	  phMax2 = PFPhoton[i].Pt();
	  iPh2 = i;
	}
      }
    }

    if(iPh2 == -99) {
      phMax2 = 0.;
      for(int i=0; i<PFPhoton.size(); i++) {
	if(PFPhoton[i].Pt() > phMax2 && i != iPh1 ) {
	  phMax2 = PFPhoton[i].Pt();
	  iPh2 = i;
	}
      }
    }
    
    if(iPh1 == -99) {
      pTPh1  = -99.;
      etaPh1 = -99.;
      phiPh1 = -99.;
      noPhoton1 = -99;
      isEMPh1 = -99;
      isLoosePh1 = -99;
      isTightPh1 = -99;
    } else {
      pTPh1  = double(PFPhoton[iPh1].Pt());
      etaPh1 = double(PFPhoton[iPh1].Eta());
      phiPh1 = double(PFPhoton[iPh1].Phi());
      noPhoton1 = int(notPhoton[iPh1]);
      isEMPh1 = int(isEM[iPh1]);
      isLoosePh1 = int(isLoose[iPh1]);
      isTightPh1 = int(isTight[iPh1]);
    }

    if(iPh2 == -99) {
      // this is gamma+jets
      pTPh2  = -99.;
      etaPh2 = -99.;
      phiPh2 = -99.;
      noPhoton2 = -99;
      isEMPh2 = -99;
      isLoosePh2 = -99;
      isTightPh2 = -99;
    } else {
      pTPh2  = double(PFPhoton[iPh2].Pt());
      etaPh2 = double(PFPhoton[iPh2].Eta());
      phiPh2 = double(PFPhoton[iPh2].Phi());
      noPhoton2 = int(notPhoton[iPh2]);
      isEMPh2 = int(isEM[iPh2]);
      isLoosePh2 = int(isLoose[iPh2]);
      isTightPh2 = int(isTight[iPh2]);
    }

    //cout<<" is Tight = "<<isTightPh1<<endl;

    // photon separation
    //    if(PFPhoton[iPh1].DeltaR(PFPhoton[iPh2])<0.8) continue;
    //    if(DeltaPhi(PFPhoton[iPh1].Phi(), PFPhoton[iPh2].Phi())<0.05) continue;

    // add the photons to the jet list;
    // just use the jets instead???? 
    //    PFPUcorrJet.push_back(PFPhoton[iPh1]);
    //    PFPUcorrJet.push_back(PFPhoton[iPh2]);
    
    // use PFMET (missing transverse energy)
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
    PFMR = -99999.;

    // hemispheres
    vector<TLorentzVector> tmpJet = CombineJets_R(PFPUcorrJet);
    if(tmpJet.size() >= 2) {
      Npassed_PFHem += _weight;
      
      TLorentzVector PFHem1 = tmpJet[0];
      TLorentzVector PFHem2 = tmpJet[1];
      
      // compute boost
      double num = PFHem1.P()-PFHem2.P();
      double den = PFHem1.Pz()-PFHem2.Pz();      
      double beta = num/den;
      
      double MT = CalcMTR(PFHem1, PFHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcMR(PFHem1, PFHem2);
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
      PFMR = variable;   
    }

    // dummy values
    passedCalo = 0;
    pTCaloHem1 = -9999;
    etaCaloHem1 = -9999;
    phiCaloHem1 = -9999;
    pTCaloHem2 = -9999;
    etaCaloHem2 = -9999;
    phiCaloHem2 = -9999;
    CaloR = -99999.;
    CaloMR = -99999.;

    // hemispheres
    tmpJet = CombineJets_R(CaloJet);
    if(tmpJet.size() >= 2) {
      Npassed_CaloHem += _weight;
    
      TLorentzVector CaloHem1 = tmpJet[0];
      TLorentzVector CaloHem2 = tmpJet[1];
      
      // compute boost
      double num = CaloHem1.P()-CaloHem2.P();
      double den = CaloHem1.Pz()-CaloHem2.Pz();      
      double beta = num/den;
      
      double MT = CalcMTR(CaloHem1, CaloHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcMR(CaloHem1, CaloHem2);
      if(variable >0) Rvariable = MT/variable;
      
      // fill the tree
      passedCalo = 1;
      pTCaloHem1 = CaloHem1.Pt();
      etaCaloHem1 = CaloHem1.Eta();
      phiCaloHem1 = CaloHem1.Phi();
      pTCaloHem2 = CaloHem2.Pt();
      etaCaloHem2 = CaloHem2.Eta();
      phiCaloHem2 = CaloHem2.Phi();
      CaloR = Rvariable;
      CaloMR = variable;    
    }

    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;
    W = _weight;

    // Muons case

    int nmuLoose = 0;
    int nmuTight = 0;
    TLorentzVector myMuLoose(0.1, 0., 0., 0.1);
    TLorentzVector myMuTight(0.1, 0., 0., 0.1);
    int iMuLoose = -99;
    int iMuTight = -99;

    for(int i=0; i<nMuon; i++) {
      if(muonPassTight(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 10.) {
	  nmuTight++;
	  if(thisMu.Pt() > myMuTight.Pt()) {
	    myMuTight = thisMu;
	    iMuTight = i;
	  }
	}
      }
    }

    for(int i=0; i<nMuon; i++) {
      if(muonPassLoose(i)) {
	TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
	if(thisMu.Pt() > 10. && iMuTight != i) {
	  // we don't count the best Tight mu among the Loose muons
	  nmuLoose++;
	  if(thisMu.Pt() > myMuLoose.Pt()) {
	    myMuLoose = thisMu;
	    iMuLoose = i;
	  }
	}
      }
    }

    // Electrons case

    int neleWP95 = 0;
    int neleWP80 = 0;
    TLorentzVector myEleWP95(0.1, 0., 0., 0.1);
    TLorentzVector myEleWP80(0.1, 0., 0., 0.1);
    int iEleWP95 = -99;
    int iEleWP80 = -99;

    for(int i=0; i<nEle; i++) {
      if(electronPassWP80(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 12.) {
	  neleWP80++;
	  if(thisEle.Pt() > myEleWP80.Pt()) {
	    myEleWP80 = thisEle;
	    iEleWP80 = i;
	  }
	}
      }
    }

    for(int i=0; i<nEle; i++) {
      if(electronPassWP95(i)) {
	TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
	if(thisEle.Pt() > 12. && i != iEleWP80) {
	  // we don't count the best WP80 Ele among the WP95 muons
	  neleWP95++;
	  if(thisEle.Pt() > myEleWP95.Pt()) {
	    myEleWP95 = thisEle;
	    iEleWP95 = i;
	  }
	}
      }
    }

    // Fill the tree per box
    pTMu1 = -999.;
    etaMu1 = -999.;
    phiMu1 = -999.;
    chargeMu1 = -999.;
    pTMu2 = -999.;
    etaMu2 = -999.;
    phiMu2 = -999.;
    chargeMu2 = -999.;    
    pTEle1 = -999.;
    etaEle1 = -999.;
    phiEle1 = -999.;
    chargeEle1 = -999.;
    pTEle2 = -999.;
    etaEle2 = -999.;
    phiEle2 = -999.;
    chargeEle2 = -999.;

    if(nmuTight>=1) { 
      // tight Mu leg
      pTMu1 = myMuTight.Pt();
      etaMu1 = myMuTight.Eta();
      phiMu1 = myMuTight.Phi();
      chargeMu1 = chargeMuon[iMuTight];
      if(nmuLoose>=1) { // DiMuon Box
	// loose Mu leg
	pTMu2 = myMuLoose.Pt();
	etaMu2 = myMuLoose.Eta();
	phiMu2 = myMuLoose.Phi();
	chargeMu2 = chargeMuon[iMuLoose];	  
	outTree[3]->Fill();
      } else if(neleWP80>=1) { // MuEle Box
	// tight Ele leg
	pTEle1 = myEleWP80.Pt();
	etaEle1 = myEleWP80.Eta();
	phiEle1 = myEleWP80.Phi();
	chargeEle1 = chargeEle[iEleWP80];
	outTree[4]->Fill();
      } else { // Mu Box
	outTree[2]->Fill();
      }
    } else if(neleWP80>=1) {
      // tight Ele leg
      pTEle1 = myEleWP80.Pt();
      etaEle1 = myEleWP80.Eta();
      phiEle1 = myEleWP80.Phi();
      chargeEle1 = chargeEle[iEleWP80];
      if(neleWP95>=1) { // DiEleBox
	// loose Ele leg
	pTEle2 = myEleWP95.Pt();
	etaEle2 = myEleWP95.Eta();
	phiEle2 = myEleWP95.Phi();
	chargeEle2 = chargeEle[iEleWP95];
	outTree[5]->Fill();
      } else { // Ele Box
	outTree[1]->Fill();
      }
    } else { // Had box
      outTree[0]->Fill();
    }

  }

  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
    
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_HLT",      &Npassed_HLT,      "Npassed_HLT/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("Npassed_2Ph",      &Npassed_2Ph,      "Npassed_2Ph/D");
  effTree->Branch("Npassed_PFHem",      &Npassed_PFHem,      "Npassed_PFHem/D");
  effTree->Branch("Npassed_CaloHem",      &Npassed_CaloHem,      "Npassed_CaloHem/D");
  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  effTree->Write();
  for(int i=0; i<6; i++) 
    outTree[i]->Write();
  file->Close();
}
  
vector<TLorentzVector> RazorDiPhoton::CombineJets_R(vector<TLorentzVector> myjets){
  
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
      // R selection
      double beta = fabs(j_temp1.P()-j_temp2.P())/fabs(j_temp1.Pz()-j_temp2.Pz());
      if(beta < 0.99) {
	// DeltaPhi selection
	if(fabs(j_temp1.DeltaPhi(j_temp2)) < 2.8) {
	  foundGood = true;
	  M_min = M_temp;
	  j1 = j_temp1;
	  j2 = j_temp2;
	}
      }
    }
  }

  if(!foundGood) return mynewjets;
  // set masses to 0
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

bool RazorDiPhoton::IsPhotonBarrel(int iPh) {
  bool isBarrel = true;
  // to fix
  double etaSeed = etaPho[iPh];
  //  if(1.479 - fabs(etaSeed) < 0.1) isBarrel = false;
  if(fabs(etaSeed) > 1.5) isBarrel = false;
  return isBarrel;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

bool RazorDiPhoton::PhotonIdBarrelisEM(int i) {
  return PhotonIdBarrel(i, "isEM");
}

bool RazorDiPhoton::PhotonIdBarrelisLoose(int i) {
  return PhotonIdBarrel(i, "isLoose");   
}

bool RazorDiPhoton::PhotonIdBarrelisTight(int i) {
  return PhotonIdBarrel(i, "isTight");
}

bool RazorDiPhoton::PhotonIdBarrel(int i, TString Selector) {
  int iSC = superClusterIndexPho[i];
  double pT = energyPho[i]*sin(thetaPho[i]);  
  bool passed = false;

  if (Selector == TString("isTight")){

    /*
    // Photon ID criteria
    if(IsPhotonBarrel(i));
    // Jurassic Isolation
    if(ecalRecHitSumEtConeDR04SC[iSC] <= 4.2+0.006*pT);
    // Tower-based HCAL isolation
    if(hcalTowerSumEtConeDR04SC[iSC] <= 2.2+0.0025*pT);
    //  H/E
    if(hOverESC[iSC] >= 0.05);
    // hollow cone track isolation
    //if(dr04HollowTkSumPtPho[i]<= 3.5+0.001*pT);
    if(dr04HollowTkSumPtPho[i]<= 2.0+0.001*pT);
    // eta width
    float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
    if(sigmaIetaIeta <= 0.013);     
    // track veto is optional... 
    if(!hasPixelSeedPho[i]);
    */

    float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
    
    if ( IsPhotonBarrel(i) && ecalRecHitSumEtConeDR04SC[iSC] < 4.2+0.006*pT && hcalTowerSumEtConeDR04SC[iSC] < 2.2+0.0025*pT && hOverESC[iSC] < 0.05 && dr04HollowTkSumPtPho[i]< 2.0+0.001*pT  && sigmaIetaIeta < 0.013 && !hasPixelSeedPho[i])


      passed  = true;

  }else if (Selector == "isLoose"){

    /*
    // Photon ID criteria
    if(IsPhotonBarrel(i));
    // Jurassic Isolation
    if(ecalRecHitSumEtConeDR04SC[iSC] <= 4.2+0.006*pT);
    // Tower-based HCAL isolation
    if(hcalTowerSumEtConeDR04SC[iSC] <= 2.2+0.0025*pT);
    //  H/E
    if(hOverESC[iSC] >= 0.05);
    // hollow cone track isolation
    if(dr04HollowTkSumPtPho[i]<= 3.5+0.001*pT);
    // if(dr04HollowTkSumPtPho[i]<= 2.0+0.001*pT);
    // eta width
    // float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
    // if(sigmaIetaIeta <= 0.013);     
    // track veto is optional... 
    if(!hasPixelSeedPho[i]);
    */

    if ( IsPhotonBarrel(i) && ecalRecHitSumEtConeDR04SC[iSC] < 4.2+0.006*pT && hcalTowerSumEtConeDR04SC[iSC] < 2.2+0.0025*pT && hOverESC[iSC] < 0.05 && dr04HollowTkSumPtPho[i]< 3.5+0.001*pT  && !hasPixelSeedPho[i])
 
      passed  = true;

  }else if (Selector == "isEM"){


    /*
  //// Photon ID criteria
  if(IsPhotonBarrel(i));
  //// Jurassic Isolation
  if(ecalRecHitSumEtConeDR04SC[iSC] <= 4.2+0.006*pT);
  //// Tower-based HCAL isolation
  if(hcalTowerSumEtConeDR04SC[iSC] <= 2.2+0.0025*pT);
  ////  H/E
  if(hOverESC[iSC] >= 0.05);
  //// hollow cone track isolation
  //if(dr04HollowTkSumPtPho[i]<= 3.5+0.001*pT);
  //if(dr04HollowTkSumPtPho[i]<= 2.0+0.001*pT);
  //// eta width
  // float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
  // if(sigmaIetaIeta <= 0.013);     
  //// track veto is optional... 
  if(!hasPixelSeedPho[i]);
    */

    if ( IsPhotonBarrel(i) && ecalRecHitSumEtConeDR04SC[iSC] < 4.2+0.006*pT && hcalTowerSumEtConeDR04SC[iSC] < 2.2+0.0025*pT && hOverESC[iSC] < 0.05 && !hasPixelSeedPho[i])

      passed  = true;

  }

  return passed;

}


////////////////////////Maurozio////////////////////////
/*
  void RazorDiPhoton::PhotonIdBarrel(int i, TString Selector) {
  int iSC = superClusterIndexPho[i];
  double pT = energyPho[i]*sin(thetaPho[i]);  
  bool passed = false;
  // Photon ID criteria
  isEM = true;
  isLoose = true;
  isTight = true;
  if(!IsPhotonBarrel(i)) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  // Jurassic Isolation
  if(ecalRecHitSumEtConeDR04SC[iSC] >= 4.2+0.006*pT) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  // Tower-based HCAL isolation
  if(hcalTowerSumEtConeDR04SC[iSC] >= 2.2+0.0025*pT) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  //  H/E
  if(hOverESC[iSC] >= 0.05) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  // hollow cone track isolation
  if(dr04HollowTkSumPtPho[i]>= 3.5+0.001*pT) {
  isLoose = false;
  }
  if(dr04HollowTkSumPtPho[i]>= 2.0+0.001*pT) {
  isTight = false;
  }
  // eta width
  float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
  if(sigmaIetaIeta >= 0.013) isTight = false;  
  
  // track veto is optional... 
  if(hasPixelSeedPho[i]) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }
  }
*/

int RazorDiPhoton::AntiPhotonIdBarrel(int i) {
  int iSC = superClusterIndexPho[i];
  double pT = energyPho[i]*sin(thetaPho[i]);
  // Photon ID criteria
  bool isEM = true;
  if(!IsPhotonBarrel(i)) {
    isEM = false;
  }
  // Jurassic Isolation
  if(ecalRecHitSumEtConeDR04SC[iSC] >= 4.2+0.006*pT) {
    isEM = false;
  }
  // Tower-based HCAL isolation
  if(hcalTowerSumEtConeDR04SC[iSC] >= 2.2+0.0025*pT) {
    isEM = false;
  }
  //  H/E
  if(hOverESC[iSC] >= 0.05) {
    isEM = false;
  }
  // hollow cone track isolation AND sigmaIetaIeta
  float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
  if(dr04HollowTkSumPtPho[i]< 2.0+0.001*pT && sigmaIetaIeta < 0.013) {
    isEM = false;
  }
  // track veto [we want photon-like jets, not electrons]
  if(hasPixelSeedPho[i]) isEM = false;
  return isEM;
}
