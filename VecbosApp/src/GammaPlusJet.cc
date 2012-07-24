// std includes
#include <iostream>
#include <string>
#include <vector>
#include <valarray>

//for absolute value function
#include <stdio.h>
#include <math.h>

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
#include "GammaPlusJet.hh"

GammaPlusJet::GammaPlusJet(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;

  _xsec = 1.;
  _Lumi = tree->GetEntries();

}

GammaPlusJet::GammaPlusJet(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;

  _xsec = 1.;
  _Lumi = tree->GetEntries();

  //To read good run list!
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

void GammaPlusJet::SetLuminosity(double Lumi){
  _Lumi = Lumi;
}

void GammaPlusJet::SetXsection(double xsec){
  _xsec = xsec;
}

GammaPlusJet::~GammaPlusJet(){}

void GammaPlusJet::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void GammaPlusJet::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");

  //////////////////////////////////////
  //       read trigger information   //
  //////////////////////////////////////

  TriggerMask mask(_treeCond);
  mask.requireTrigger("HLT_Photon10_L1R");
  mask.requireTrigger("HLT_Photon15_L1R");

  vector<int> requiredTrigger = mask.getBits();
  vector<int> requiredTriggerPh10;
  vector<int> requiredTriggerPh15;
  requiredTriggerPh10.push_back(requiredTrigger[0]);
  requiredTriggerPh15.push_back(requiredTrigger[1]);

  // event information
  double weight; 

  // photon
  double etaPhoton;
  double phiPhoton;
  double ptPhoton;
  double matchPh;
  // MET
  double MET;
  double phiMET;
  // Jets
  double njets;
  // HLT bits
  double HLT_Photon10;
  double HLT_Photon15;

  // photon info
  double nPhotons;

 // event ID
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;

  //Variables Natasha added for Photon Analysis (June 22 2010)
  double EcalSum03SC;
  double EcalSum04SC;
  double ecalRecHit03SC;
  double hcalTow03SC;
  double trkSumPt03SC;
  double ecalRecSum04SC;
  double hcalTowSum04SC;
  double trkSumPt04SC;
  double trkIndEle;
  double covEtaEtaSC;
  double covEtaPhiSC;
  double covPhiPhiSC;
  double e4SwiCrossSC;
  double eMxSC;
  double tSC;
  double HE;  //variable Emanual suggested HE = HOverESC/eSC
  double HvrESC;
  double eSC;
  
// prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("weight",     &weight,     "weight/D");
  outTree->Branch("etaPhoton",     &etaPhoton,     "etaPhoton/D");
  outTree->Branch("ptPhoton",      &ptPhoton,      "ptPhoton/D");
  outTree->Branch("phiPhoton",     &phiPhoton,     "phiPhoton/D");
  outTree->Branch("matchPh",     &matchPh,     "matchPh/D");
  outTree->Branch("MET",        &MET,        "MET/D");
  outTree->Branch("phiMET",     &phiMET,     "phiMET/D");
  outTree->Branch("njets",      &njets,      "njets/D");
  outTree->Branch("nPhotons",      &nPhotons,      "nPhotons/D");
  outTree->Branch("HLT_Photon10", &HLT_Photon10, "HLT_Photon10/D");
  outTree->Branch("HLT_Photon15", &HLT_Photon15, "HLT_Photon15/D");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");

//Variables Natasha added to output tree to determine variable cutoffs for Photon Analysis (June 22, 2010)
  
  outTree->Branch("HvrESC", &HvrESC, "HvrESC/D");
  outTree->Branch("eSC", &eSC, "eSC/D");
  outTree->Branch("HE", &HE, "HE/D");
  outTree->Branch("EcalSum03SC", &EcalSum03SC, "EcalSum03SC/D");
  outTree->Branch("EcalSum04SC", &EcalSum04SC, "EcalSum04SC/D");
  outTree->Branch("ecalRecHit03SC", &ecalRecHit03SC, "ecalRecHit03SC/D");
  outTree->Branch("hcalTow03SC", &hcalTow03SC, "hcalTow03SC/D");
  outTree->Branch("trkSumPt03SC", &trkSumPt03SC, "trkSumPt03SC/D");
  outTree->Branch("ecalRecSum04SC", &ecalRecSum04SC, "ecalRecSum04SC/D");
  outTree->Branch("hcalTowSum04SC", &hcalTowSum04SC, "hcalTowSum04SC/D");
  outTree->Branch("trkSumPt04SC", &trkSumPt04SC, "trkSumPt04SC/D");
  outTree->Branch("trkIndEle", &trkIndEle, "trkIndEle/D");
  outTree->Branch("covEtaEtaSC", &covEtaEtaSC, "covEtaEtaSC/D");
  outTree->Branch("covEtaPhiSC", &covEtaPhiSC, "covEtaPhiSC/D");
  outTree->Branch("covPhiPhiSC", &covPhiPhiSC, "covPhiPhiSC/D");
  outTree->Branch("e4SwiCrossSC", &e4SwiCrossSC, "e4SwiCrossSC/D");
  outTree->Branch("eMxSC", &eMxSC, "eMxSC/D");
  outTree->Branch("tSC", &tSC, "tSC/D");







  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    weight = 1.;
    if(!_isData) weight = _Lumi*_xsec/double(nentries);

    //////////////////////
    // photon selection //
    //////////////////////

    vector<int> goodPhoton;

    for(int i=0; i<nSC; i++) 
    {

      cout << covIEtaIEtaSC[i] << endl;

      if(!isElectron(i)) //not an electron
      {
	if( fabs(etaSC[i])<3.0 ) //|eta| cut < 3
	{
	  if( energySC[i]*sin(thetaSC[i])>50.0 ) 
	  { 
	    if(1)  // cluster shape : eMaxSC/e3x3SC is ? use covlEtalEta
	    {
	      goodPhoton.push_back(i);
	    }
	  }
	}
      }
    }

    if(int(goodPhoton.size()) < 1 ) continue; // no isolated photon found

    // best candidate selection
    
    // by highest energy
    int iGoodPhoton = HighestPtSC(goodPhoton);
    // by more isolation


    // MC matching
    matchPh = 0.;
    if(!_isData) {
      for(i =0; i<nMc; i++) {
	// status
	// id
	//dR
	matchPh = 1.;
      }
    }

    //////////////////////
    //   jet selection  //
    //////////////////////
    
    njets = 0;
    for(int j=0; j<nAK5Jet; j++) {
      if(sqrt(pow(double(pxAK5Jet[j]),2.) + pow(double(pyAK5Jet[j]),2.)) > 30.) { // pT requirement > 30 GeV
	if(fabs(etaAK5Jet[j])<3.) { // |eta requirement|< 3  (fabs(double) )
	  if((DeltaR(etaSC[iGoodPhoton], phiSC[iGoodPhoton],      //isolation
		     etaAK5Jet[j], phiAK5Jet[j]) > 0.35)) 
	    {
	      njets++;
	    }
	}
      }
    }

    // write the tree
    etaPhoton = etaSC[iGoodPhoton];
    phiPhoton = phiSC[iGoodPhoton];
    ptPhoton = energySC[iGoodPhoton]*sin(thetaSC[iGoodPhoton]);
    MET = sqrt(pow(double(pxMet[0]),2.)+pow(double(pyMet[0]),2.));
    phiMET = phiMet[0];
    

    //where Natasha wrote more to tree for Photon Analysis (June 22 2010)
    
   
    HvrESC = hOverESC[iGoodPhoton];
    eSC = energySC[iGoodPhoton];
    HE = HvrESC/eSC;
    EcalSum03SC = scBasedEcalSum03SC[iGoodPhoton];
    EcalSum04SC = scBasedEcalSum04SC[iGoodPhoton];
    ecalRecHit03SC = ecalRecHitSumEtConeDR03SC[iGoodPhoton];
    hcalTow03SC = hcalTowerSumEtConeDR03SC[iGoodPhoton];
    trkSumPt03SC = trkSumPtSolidConeDR03SC[iGoodPhoton];
    ecalRecSum04SC = ecalRecHitSumEtConeDR04SC[iGoodPhoton];
    hcalTowSum04SC = hcalTowerSumEtConeDR04SC[iGoodPhoton];
    trkSumPt04SC = trkSumPtSolidConeDR04SC[iGoodPhoton];
    trkIndEle = trackIndexEle[iGoodPhoton];
    covEtaEtaSC = covIEtaIEtaSC[iGoodPhoton];
    covEtaPhiSC = covIEtaIPhiSC[iGoodPhoton];

    covPhiPhiSC = covIPhiIPhiSC[iGoodPhoton];
    e4SwiCrossSC = e4SwissCrossSC[iGoodPhoton];
    eMxSC = eMaxSC[iGoodPhoton];
    tSC = timeSC[iGoodPhoton];
    

    // HLT 
    Utils anaUtils;
    HLT_Photon10 = int(anaUtils.getTriggersOR(requiredTriggerPh10, firedTrg));
    HLT_Photon15 = int(anaUtils.getTriggersOR(requiredTriggerPh15, firedTrg));

    // photon info
    nPhotons = double(goodPhoton.size());

    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;

    outTree->Fill();

    // clean the vectors
    goodPhoton.clear();
  }
  
  outTree->Write();
  file->Close();
  
}

int GammaPlusJet::HighestPtSC(vector<int> indPhotons) {
  double maxPt = 0.;
  int iHpt = -99;
  for(int i=0; i< indPhotons.size(); i++) {
    int iPh = indPhotons[i];
    if(energySC[iPh]*sin(thetaSC[iPh]) > maxPt) {
      maxPt = energySC[iPh]*sin(thetaSC[iPh]);
      iHpt = iPh;
    }
  }

  return iHpt;
}

bool GammaPlusJet::isElectron(int iPh) {
  bool isEle = false;
  for(int iEle = 0; iEle<nEle; iEle++) {
    if(superClusterIndexEle[iEle] == iPh) {
      isEle = true;
    }
  }

  return isEle;

}
