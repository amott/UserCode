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
#include "VBTFLeptEff.hh"


VBTFLeptEff::VBTFLeptEff(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  
}

VBTFLeptEff::VBTFLeptEff(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;

  //To read good run list!
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

VBTFLeptEff::~VBTFLeptEff() {}

bool VBTFLeptEff::MatchEle(int eEle) {
  bool matched = false;
  TLorentzVector myRecoEle(pxEle[eEle], pyEle[eEle], pzEle[eEle], energyEle[eEle]);
  for(int i=0; i<nMc; i++) {
    if(statusMc[i] == 1 && abs(idMc[i]) == 11) {
      TLorentzVector myGenEle;
      myGenEle.SetPtEtaPhiE(pMc[i]*sin(thetaMc[i]), etaMc[i], phiMc[i], energyMc[i]);
      if( myRecoEle.DeltaR(myGenEle)<0.1) matched = true;
    }
  }
  return matched;
}

void VBTFLeptEff::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

bool VBTFLeptEff::isGlobalMuon(int iMu){
  bool ret = false;
  Utils anaUtils;
  bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
  bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
  bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);

  if(isMuGlobal && isMuGlobalPrompt && isMuTracker)ret=true;

  return ret;
}

void VBTFLeptEff::Loop(string outFileName, int start, int stop) {
 
  if(fChain == 0) return;

  double Npassed_In=0;
  double Npassed_HLT=0;
  double Npassed_PV=0;
  double Npassed_MuID=0;
  double Npassed_MuIso=0;
  double Npassed_MuStrip=0;
  double Npassed_MuPixel=0;
  double Npassed_MupT=0;
  double Npassed_MuIP=0;
  double Npassed_MuEta=0;
  double Npassed_TrueEle=0;
  double Npassed_EleEcalDriven=0;
  double Npassed_EleEt=0;
  double Npassed_EleEta=0;
  double Npassed_EleBarrel=0;
  double Npassed_EleEndcap=0;
  double Npassed_Ele80BarrelIso=0;
  double Npassed_Ele80BarrelID=0;
  double Npassed_Ele80BarrelSig=0;
  double Npassed_Ele80BarrelTrackCut=0;
  double Npassed_Ele80EndcapIso=0;
  double Npassed_Ele80EndcapID=0;
  double Npassed_Ele80EndcapSig=0;
  double Npassed_Ele80EndcapTrackCut=0;
  double Npassed_Ele95BarrelIso=0;
  double Npassed_Ele95BarrelID=0;
  double Npassed_Ele95BarrelSig=0;
  double Npassed_Ele95BarrelTrackCut=0;
  double Npassed_Ele95EndcapIso=0;
  double Npassed_Ele95EndcapID=0;
  double Npassed_Ele95EndcapSig=0;
  double Npassed_Ele95EndcapTrackCut=0;

  double weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  double MupT[100], MuIP[100], EleET[100];
  int nMu, nElectron;

  // out tree
  TTree* outTree = new TTree("outTree","outTree");
  outTree->Branch("nMu", &nMu, "nMu/I");
  outTree->Branch("nElectron", &nElectron, "nElectron/I");
  outTree->Branch("MuIP", MuIP, "MuIP[nMu]/D");
  outTree->Branch("MupT", MupT, "MupT[nMu]/D");
  outTree->Branch("EleET", EleET, "EleET[nElectron]/D");

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


#if AnalysisSelector==1
    
    //nMu=nMuon;
    int iMupT=0;
    for(int iMu=0; iMu<nMuon; iMu++){

      Utils anaUtils;
      bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
      bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
      bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);
      
      if(isMuGlobal && isMuGlobalPrompt && isMuTracker){
	
	Npassed_MuID+=weight;
	
	double pt = sqrt(pxMuon[iMu]*pxMuon[iMu]+pyMuon[iMu]*pyMuon[iMu]);

	int iTrack = trackIndexMuon[iMu];

	//	MuIP[iMu]=fabs(transvImpactParTrack[iTrack]);
        //MupT[iMu]=pt;

	if(numberOfValidStripTIBHitsTrack[iTrack]+
	   numberOfValidStripTIDHitsTrack[iTrack]+
	   numberOfValidStripTOBHitsTrack[iTrack]+
	   numberOfValidStripTECHitsTrack[iTrack] > 10){
	  
	  Npassed_MuStrip+=weight;
	  
	  if(numberOfValidPixelBarrelHitsTrack[iTrack] > 0 ||
	     numberOfValidPixelEndcapHitsTrack[iTrack] > 0){
	    
	    Npassed_MuPixel+=weight;
	    
	    if(fabs(transvImpactParTrack[iTrack]) < 0.2){
	      //if(trackNormalizedChi2GlobalMuonTrack[iTrack] < 10){                                                                                           
	      Npassed_MuIP+=weight;
	      
	      double IECAL = emEt03Muon[iMu];
	      double IHCAL = hadEt03Muon[iMu];
	      double ITRK  = sumPt03Muon[iMu];
	      if(pt > 20.){
		Npassed_MupT+=weight;

	      if(fabs(etaMuon[iMu])< 2.1){
		Npassed_MuEta+=weight;
		pt = sqrt(pxMuon[iMu]*pxMuon[iMu]+pyMuon[iMu]*pyMuon[iMu]);
		iTrack = trackIndexMuon[iMu];
		MuIP[iMupT]=fabs(transvImpactParTrack[iTrack]);
		MupT[iMupT]=pt;
		iMupT++;
		//		if(pt > 20.){
		//		  Npassed_MupT+=weight;
		  if((IECAL+IHCAL+ITRK)/pt < 0.15){
		    Npassed_MuIso+=weight;  
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    nMu=iMupT;
#endif

#if AnalysisSelector==3
    nElectron=0;

    for(int iEle=0; iEle<nEle; iEle++){
    Utils anaUtils;

    // consider only real electrons
    // and remove fakes
    if(!MatchEle(iEle)) continue;
    Npassed_TrueEle+= weight;

    //is an ECAL driven electron                                                                                                                               
    bool isECALdriven = anaUtils.electronRecoType(recoFlagsEle[iEle], bits::isEcalDriven);
    
    if(isECALdriven == true){

      Npassed_EleEcalDriven+=weight;
      
      bool isBarrel = true;
      bool isEndcap = true;
      int iSC = superClusterIndexEle[iEle];
      double ETA = etaSC[iSC];
      double ET = energySC[iSC]/cosh(etaSC[iSC]);
      EleET[nElectron]=ET;
      nElectron++;
      if(ET > 20.){
	
	Npassed_EleEt+=weight;
	
	//is the electron in the ECAL fiducial region?                                                                                                             
	if(fabs(ETA) > 1.4442) isBarrel = false;
	if(fabs(ETA) > 2.5 || fabs(ETA) < 1.566) isEndcap = false;
	
	if(isBarrel == true || isEndcap == true){
	  
	  Npassed_EleEta+=weight;
	  
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
	  
	  trackISO /= ET;
	  ECALISO /= ET;
	  HCALISO /= ET;
	  
	  if(isBarrel){
	    Npassed_EleBarrel+=weight;
	    //ISO WP80
	    if(trackISO < 0.09 && ECALISO < 0.07 && HCALISO < 0.10){
	      Npassed_Ele80BarrelIso+=weight;
	      if(fabs(dphi) < 0.06){
		if(fabs(deta) < 0.004){
		  if(HoE < 0.04){
		    Npassed_Ele80BarrelID+=weight;
		    if(sigietaieta < 0.01){
		      Npassed_Ele80BarrelSig+=weight;
		    if(fabs(convDist) <= 0.02 && fabs(convDcot) <= 0.02){
		      int iTrack = gsfTrackIndexEle[iEle];
		      if(expInnerLayersGsfTrack[iTrack] > 0)Npassed_Ele80BarrelTrackCut+=weight;
		    }
		    }
		  }
		}
	      }
	    }
	  } else {
	    Npassed_EleEndcap+=weight;
	    //ISO WP80
	    if(trackISO < 0.04 && ECALISO < 0.05 && HCALISO < 0.025){
	      Npassed_Ele80EndcapIso+=weight;
	      //ID WP80                      
	      if(fabs(dphi) < 0.03){
		if(HoE < 0.025){
		  Npassed_Ele80EndcapID+=weight;
		  if(sigietaieta < 0.03){
		    Npassed_Ele80EndcapSig+=weight;
		    
		    if(fabs(convDist) <= 0.02 && fabs(convDcot) <= 0.02){
		      int iTrack = gsfTrackIndexEle[iEle];
		      if(expInnerLayersGsfTrack[iTrack] > 0)Npassed_Ele80EndcapTrackCut+=weight;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
#endif

#if AnalysisSelector==4

for(int iEle=0; iEle<nEle; iEle++){

  Utils anaUtils;

  //is an ECAL driven electron                                                                                                                               
  bool isECALdriven = anaUtils.electronRecoType(recoFlagsEle[iEle], bits::isEcalDriven);

  // and remove fakes
  if(!MatchEle(iEle)) continue;
  Npassed_TrueEle+= weight;

  if(isECALdriven == true){

    Npassed_EleEcalDriven+=weight;

    bool isBarrel = true;
    bool isEndcap = true;
    int iSC = superClusterIndexEle[iEle];
    double ETA = etaSC[iSC];
    double ET = energySC[iSC]/cosh(etaSC[iSC]);
    
    if(ET > 20.0){

      Npassed_EleEt+=weight;

      //is the electron in the ECAL fiducial region?                                                                                                             
      if(fabs(ETA) > 1.4442) isBarrel = false;
      if(fabs(ETA) > 2.5 || fabs(ETA) < 1.566) isEndcap = false;

      if(isBarrel == true || isEndcap == true){

	Npassed_EleEta+=weight;

	double trackISO = dr03TkSumPtEle[iEle];
	double ECALISO = dr03EcalRecHitSumEtEle[iEle];
	double HCALISO = dr03HcalTowerSumEtEle[iEle];
	
	double pt = sqrt(pxEle[iEle]*pxEle[iEle]+pyEle[iEle]*pyEle[iEle]);
	
	double sigietaieta = covIEtaIEtaSC[iSC];
	double dphi = deltaPhiAtVtxEle[iEle];
	double deta = deltaEtaAtVtxEle[iEle];
	double HoE = hOverEEle[iEle];

	trackISO /= ET;
	ECALISO /= ET;
	HCALISO /= ET;
	
	if(isBarrel){
	  Npassed_EleBarrel+=weight;
	  //ISO WP95
	  if(trackISO < 0.15 && ECALISO < 2.0 && HCALISO < 0.12){
	    Npassed_Ele95BarrelIso+=weight;
	    //ID WP95 
	    if(fabs(dphi) < 0.007){
	      if(fabs(deta) < 0.15){
		if(HoE < 0.01){
		  Npassed_Ele95BarrelID+=weight;
		  if(sigietaieta < 0.01){
		    Npassed_Ele95BarrelSig+=weight;
		    int iTrack = gsfTrackIndexEle[iEle];
		    if(expInnerLayersGsfTrack[iTrack] > 1)Npassed_Ele95BarrelTrackCut+=weight;
		  }
		}
	      }
	    }
	  }
	} else {
	  Npassed_EleEndcap+=weight;
	  //ISO WP95
	  if(trackISO < 0.08 && ECALISO < 0.06 && HCALISO < 0.05){
	    Npassed_Ele95EndcapIso+=weight;
	    //ID WP95                      
	    if(HoE < 0.07){
	      Npassed_Ele95EndcapID+=weight;
	      if(sigietaieta < 0.03){
		Npassed_Ele95EndcapSig+=weight;
		int iTrack = gsfTrackIndexEle[iEle];
		if(expInnerLayersGsfTrack[iTrack] > 1)Npassed_Ele95EndcapTrackCut+=weight;
	      }
	    }
	  }
	}
      }
    }
  }
 }
#endif

 outTree->Fill();

  }

  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
    
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_HLT",      &Npassed_HLT,      "Npassed_HLT/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("Npassed_MuIP",      &Npassed_MuIP,      "Npassed_MuIP/D");
  effTree->Branch("Npassed_MuID",      &Npassed_MuID,      "Npassed_MuID/D");
  effTree->Branch("Npassed_MuIso",      &Npassed_MuIso,      "Npassed_MuIso/D");
  effTree->Branch("Npassed_MuPixel",      &Npassed_MuPixel,      "Npassed_MuPixel/D");
  effTree->Branch("Npassed_MuStrip",      &Npassed_MuStrip,"Npassed_MuStrip/D");
  effTree->Branch("Npassed_MupT",      &Npassed_MupT, "Npassed_MupT/D");
  effTree->Branch("Npassed_MuEta",      &Npassed_MuEta, "Npassed_MuEta/D");
  effTree->Branch("Npassed_TrueEle", &Npassed_TrueEle, "Npassed_TrueEle/D");
  effTree->Branch("Npassed_EleEcalDriven", &Npassed_EleEcalDriven,"Npassed_EleEcalDriven/D");
  effTree->Branch("Npassed_EleEt", &Npassed_EleEt, "Npassed_EleEt/D");
  effTree->Branch("Npassed_EleEta", &Npassed_EleEta, "Npassed_EleEta/D"); 
  effTree->Branch("Npassed_EleBarrel", &Npassed_EleBarrel ,"Npassed_EleBarrel/D" ); 
  effTree->Branch("Npassed_EleEndcap", &Npassed_EleEndcap ,"Npassed_EleEndcap/D" );
  effTree->Branch("Npassed_Ele80BarrelIso", &Npassed_Ele80BarrelIso ,"Npassed_Ele80BarrelIso/D" ); 
  effTree->Branch("Npassed_Ele80BarrelID", &Npassed_Ele80BarrelID ,"Npassed_Ele80BarrelID/D" ); 
  effTree->Branch("Npassed_Ele80BarrelSig", &Npassed_Ele80BarrelSig ,"Npassed_Ele80BarrelSig/D" ); 
  effTree->Branch("Npassed_Ele80BarrelTrackCut", &Npassed_Ele80BarrelTrackCut ,"Npassed_Ele80BarrelTrackCut/D" ); 
  effTree->Branch("Npassed_Ele80EndcapIso", &Npassed_Ele80EndcapIso ,"Npassed_Ele80EndcapIso/D" ); 
  effTree->Branch("Npassed_Ele80EndcapID", &Npassed_Ele80EndcapID ,"Npassed_Ele80EndcapID/D" ); 
  effTree->Branch("Npassed_Ele80EndcapSig", &Npassed_Ele80EndcapSig ,"Npassed_Ele80EndcapSig/D" ); 
  effTree->Branch("Npassed_Ele80EndcapTrackCut", &Npassed_Ele80EndcapTrackCut ,"Npassed_Ele80EndcapTrackCut/D" ); 
  effTree->Branch("Npassed_Ele95BarrelIso", &Npassed_Ele95BarrelIso ,"Npassed_Ele95BarrelIso/D" ); 
  effTree->Branch("Npassed_Ele95BarrelID", &Npassed_Ele95BarrelID ,"Npassed_Ele95BarrelID/D" ); 
  effTree->Branch("Npassed_Ele95BarrelSig", &Npassed_Ele95BarrelSig ,"Npassed_Ele95BarrelSi/D" ); 
  effTree->Branch("Npassed_Ele95BarrelTrackCut", &Npassed_Ele95BarrelTrackCut ,"Npassed_Ele95BarrelTrackCut/D" ); 
  effTree->Branch("Npassed_Ele95EndcapIso", &Npassed_Ele95EndcapIso ,"Npassed_Ele95EndcapIso/D" ); 
  effTree->Branch("Npassed_Ele95EndcapID", &Npassed_Ele95EndcapID ,"Npassed_Ele95EndcapID/D" ); 
  effTree->Branch("Npassed_Ele95EndcapSig", &Npassed_Ele95EndcapSig ,"Npassed_Ele95EndcapSig/D" ); 
  effTree->Branch("Npassed_Ele95EndcapTrackCut", &Npassed_Ele95EndcapTrackCut, "Npassed_Ele95EndcapTrackCut/D"); 

  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  effTree->Write();
  outTree->Write();
  file->Close();

}
