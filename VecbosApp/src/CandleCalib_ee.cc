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
#include "CandleCalib_ee.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

CandleCalibee::CandleCalibee(TTree *tree) : Vecbos(tree) {}

CandleCalibee::~CandleCalibee(){}

void CandleCalibee::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  TFile* triggemask = new TFile("triggermask.root");
  TTree* treeCond = (TTree*) triggemask->Get("Conditions");
  TriggerMask mask(treeCond);

  mask.requireTrigger("HLT_Ele10_SW_L1R");
  mask.requireTrigger("HLT_Ele15_SW_L1R");
  mask.requireTrigger("HLT_Ele15_LW_L1R");
  mask.requireTrigger("HLT_IsoEle15_L1I");
  mask.requireTrigger("HLT_IsoEle18_L1R");
  mask.requireTrigger("HLT_IsoEle15_LW_L1I");
  mask.requireTrigger("HLT_LooseIsoEle15_LW_L1R");

  std::vector<int> requiredTriggers = mask.getBits();

  double njets;
  double ntrackjets;
  double npartons;
  double processID;
  double mll;
  double mllGEN;
  double Zpt;
  double ZptGEN;
  double etaEle1;
  double etaEle2;
  double ptEle1;
  double ptEle2;
  double phiEle1;
  double phiEle2;
  double genMET;
  double genphiMET;
  double MET;
  double phiMET;
  double ptJet1;
  double etaJet1;
  double phiJet1;
  double ptJet2;
  double etaJet2;
  double phiJet2;
  double ptJet3;
  double etaJet3;
  double phiJet3;
  double emfJet1;
  double emfJet2;
  double emfJet3;
  double sumPtJet;
  double sumPtOverPtEle1;
  double sumPtOverPtEle2;
  double MHTphiMET;
  double MHTphiJet;
  double TSphiMET;
  double sinMHTphiMET;
  double sinTSphiMET;
  double HLT_Ele10_SW_L1R;
  double HLT_Ele15_SW_L1R;
  double HLT_Ele15_LW_L1R;
  double HLT_IsoEle15_L1I;
  double HLT_IsoEle18_L1R;
  double HLT_IsoEle15_LW_L1I;
  double HLT_LooseIsoEle15_LW_L1R;
  double mJet;
  double sinJMET1;
  double sinJMET2;
  double sinJMET3;
  double sinJMET4;

  //   double weight; THIS IS GLOBAL
  TH1D* HNcaloJets = new TH1D("NCaloJets", "NCaloJets", 6, 0., 6.);
  TH1D* HNtrackJets = new TH1D("NTrackJets", "NTrackJets", 6, 0., 6.);
  TH1D* HeleIDHardestEle = new TH1D("HeleIDHardestEle", "HeleIDHardestEle", 100, 0., 1.);

  TH1D* aveNtrack = new TH1D("aveNtrack", "aveNtrack", 100, 0., 100.);
  TH1D* fNtrack = new TH1D("fNtrack", "fNtrack", 50, 0., 1.);
  TH1D* dxytrack = new TH1D("dxytrack", "dxytrack", 100, -0.1, 0.1);
  TH1D* dztrack = new TH1D("dztrack", "dztrack", 100, -0.3, 0.3);

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("mllGEN",        &mllGEN,        "mllGEN/D");
  outTree->Branch("ZptGEN",        &ZptGEN,        "ZptGEN/D");
  outTree->Branch("Zpt",        &Zpt,        "Zpt/D");
  outTree->Branch("njets",      &njets,      "njets/D");
  outTree->Branch("ntrackjets", &ntrackjets, "ntrackjets/D");
  outTree->Branch("npartons",   &npartons,   "npartons/D");
  outTree->Branch("processID",  &processID,   "processID/D");
  outTree->Branch("mll",        &mll,        "mll/D");
  outTree->Branch("etaEle1",     &etaEle1,     "etaEle1/D");
  outTree->Branch("etaEle2",     &etaEle2,     "etaEle2/D");
  outTree->Branch("ptEle1",      &ptEle1,      "ptEle1/D");
  outTree->Branch("ptEle2",      &ptEle2,      "ptEle2/D");
  outTree->Branch("phiEle1",     &phiEle1,     "phiEle1/D");
  outTree->Branch("phiEle2",     &phiEle2,     "phiEle2/D");
  outTree->Branch("genMET",     &genMET,     "genMET/D");
  outTree->Branch("genphiMET",  &genphiMET,  "genphiMET/D");
  outTree->Branch("MET",        &MET,        "MET/D");
  outTree->Branch("phiMET",     &phiMET,     "phiMET/D");
  outTree->Branch("ptJet1",     &ptJet1,     "ptJet1/D");
  outTree->Branch("etaJet1",    &etaJet1,    "etaJet1/D");
  outTree->Branch("phiJet1",    &phiJet1,    "phiJet1/D");
  outTree->Branch("ptJet2",     &ptJet2,     "ptJet2/D");
  outTree->Branch("etaJet2",    &etaJet2,    "etaJet2/D");
  outTree->Branch("phiJet2",    &phiJet2,    "phiJet2/D");
  outTree->Branch("ptJet3",     &ptJet3,     "ptJet3/D");
  outTree->Branch("etaJet3",    &etaJet3,    "etaJet3/D");
  outTree->Branch("phiJet3",    &phiJet3,    "phiJet3/D");
  outTree->Branch("emfJet1",    &emfJet1,    "emfJet1/D");
  outTree->Branch("emfJet2",    &emfJet2,    "emfJet2/D");
  outTree->Branch("emfJet3",    &emfJet3,    "emfJet3/D");
  outTree->Branch("sumPtJet",   &sumPtJet,   "sumPtJet/D");
  outTree->Branch("sumPtOverPtEle1",   &sumPtOverPtEle1,   "sumPtOverPtEle1/D");
  outTree->Branch("sumPtOverPtEle2",   &sumPtOverPtEle2,   "sumPtOverPtEle2/D");
  outTree->Branch("weight",     &weight,     "weight/D");
  outTree->Branch("MHTphiMET",  &MHTphiMET,  "MHTphiMET/D");
  outTree->Branch("MHTphiJet",  &MHTphiJet,  "MHTphiJet/D");
  outTree->Branch("TSphiMET",   &TSphiMET,   "TSphiMET/D");
  outTree->Branch("sinMHTphiMET",   &sinMHTphiMET,   "sinMHTphiMET/D");
  outTree->Branch("sinTSphiMET",    &sinTSphiMET,    "sinTSphiMET/D");
  outTree->Branch("HLT_Ele10_SW_L1R",         &HLT_Ele10_SW_L1R,         "HLT_Ele10_SW_L1R/D");
  outTree->Branch("HLT_Ele15_SW_L1R",         &HLT_Ele15_SW_L1R,         "HLT_Ele15_SW_L1R/D");
  outTree->Branch("HLT_Ele15_LW_L1R",         &HLT_Ele15_LW_L1R,         "HLT_Ele15_LW_L1R/D");
  outTree->Branch("HLT_IsoEle15_L1I",         &HLT_IsoEle15_L1I,         "HLT_IsoEle15_L1I/D");
  outTree->Branch("HLT_IsoEle18_L1R",         &HLT_IsoEle18_L1R,         "HLT_IsoEle18_L1R/D");
  outTree->Branch("HLT_IsoEle15_LW_L1I",      &HLT_IsoEle15_LW_L1I,      "HLT_IsoEle15_LW_L1I/D");
  outTree->Branch("HLT_LooseIsoEle15_LW_L1R", &HLT_LooseIsoEle15_LW_L1R, "HLT_LooseIsoEle15_LW_L1R/D");
  outTree->Branch("mJet", &mJet, "mJet/D");
  outTree->Branch("sinJMET1", &sinJMET1, "sinJMET1/D");
  outTree->Branch("sinJMET2", &sinJMET2, "sinJMET2/D");
  outTree->Branch("sinJMET3", &sinJMET3, "sinJMET3/D");
  outTree->Branch("sinJMET4", &sinJMET4, "sinJMET4/D");

  // prepare vectors for efficiency
  vector<string> cuts;
  vector<double> Npassed;
  vector<double> Npassed_2j;
  vector<double> Npassed_3j;
  vector<double> Npassed_4j;
  
  cuts.push_back("One jet, no other cuts");
  cuts.push_back("Electron Selection [PV independent]");
  cuts.push_back("Z selection");
  cuts.push_back("Electron |dz| cut [PV dependent]");
  cuts.push_back("Electron isolation [PV dependent]");
  cuts.push_back("|sin(MHTphiMET)|<0.85");
  cuts.push_back("At least 1 calo jet");
  cuts.push_back("At least 1 track jet");

  for(int i=0; i< int(cuts.size())+1; i++) {
    Npassed.push_back(0.);
    Npassed_2j.push_back(0.);
    Npassed_3j.push_back(0.);
    Npassed_4j.push_back(0.);
  }

  int nCount = 0;

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
    
    // W/Z+jet        weight = 1./3.;
    // ttbar    weight = 1./30.;

    // QCD  EM enriched Pt20to30          weight = 100.*0.0080*0.40*pow(10.,9.)/(nentries*1.);
    // QCD  EM enriched Pt30to80  weight = 100.*0.0470*0.01*pow(10.,9.)/(nentries*1.);
    // QCD  EM enriched Pt80to170 weight = 100.*0.15*0.0019*pow(10.,9.)/(nentries*1.);

    // QCD  bce Pt20to30  weight = 100.*0.00048*0.40*pow(10.,9.)/(nentries*1.);
    // QCD  bce Pt30to80  weight = 100.*0.0024*0.01*pow(10.,9.)/(nentries*1.);
    // QCD  bce Pt80to170 weight = 100.*0.012*0.0019*pow(10.,9.)/(nentries*1.);

    weight = 1.;
    npartons = 0;
    processID = 0.;

    // compute how many events have at least 1 jet
    // Did the Z decay to ee?
    // Look for the gen-level electron from Z
    int iGene1 = -99;
    for(int i=0; i<nMc; i++) 
      if(statusMc[i] == 3) // stable particle
   	if(abs(idMc[i]) == 11) // it is an electron
	  //	  if(abs(idMc[mothMc[i]]) == 23) // its mother is a Z
  	    iGene1 = i;
 
    int iGene2 = -99;
    for(int i=0; i<nMc; i++) 
      if(statusMc[i] == 3) // stable particle
   	if(abs(idMc[i]) == 11) // it is an electron
	  //	  if(abs(idMc[mothMc[i]]) == 23) // its mother is a Z
	    if(i != iGene1)
	      iGene2 = i;
    // RecHit thresholds for CaloTowers                                                                                                                                                                                                    
    vector<float> thresh;
    thresh.push_back(.9);       //HB 
    thresh.push_back(99999.0);      //etc...
    thresh.push_back(1.4);
    thresh.push_back(1.4);
    thresh.push_back(1.2);
    thresh.push_back(1.8);
    thresh.push_back(.09);
    thresh.push_back(.45);
    thresh.push_back(.2);
    thresh.push_back(.45);
    CaloF = 0.06;
    double jet_et_cut = 30.;
    double jet_pt_cut = 15.;

    int nj = 0;
    if(!usePF) {
      // determine the number of CALO jets in the event
      int nJetCalo = 0;
      for(int i=0; nAK5Jet; i++) {
	TVector3 jet(pxAK5Jet[i], pyAK5Jet[i], pzAK5Jet[i]);
	bool matchesEle = false;
	for(int j=0; nEle; j++) { 
	  TVector3 ele(pxEle[j], pyEle[j], pzEle[j]);
	  // angular matching
	  if(ele.DeltaR(jet) < 0.3) matchesEle = true;
	}
	if(matchesEle) nJetCalo++;
      }
      nj = nJetCalo;
    } else {
      // determine the number of PF jets
      int nPFJet = 0;
      for(int i=0; nAK5PFJet; i++) {
	TVector3 jet(pxAK5PFJet[i], pyAK5PFJet[i], pzAK5PFJet[i]);
	bool matchesEle = false;
	for(int j=0; nPFEle; j++) {
	  TVector3 ele(pxPFEle[j], pyPFEle[j], pzPFEle[j]);
	  // angular matching                                                                                                                                                                                                                
	  if(ele.DeltaR(jet) < 0.3) matchesEle = true;
	}
	if(matchesEle) nPFJet++;
      }
      nj = nPFJet;
    }

    // TO BE FIXED IN A REAL WAY!!!!
    if(nj>0) Npassed[0]+= weight;
    if(nj>1) Npassed_2j[0] += weight;
    if(nj>2) Npassed_3j[0] += weight;
    if(nj>3) Npassed_4j[0] += weight;
    
    ///////////////////////
    // electron selection
    ///////////////////////

    vector<int> goodEle;

    // number of valid hits
    Utils anaUtils;
    for(int i=0; i<nEle; i++) {
      // pT cut
      if(pTEle(i) < 10.) continue;
      // electronID
      bool eleId = anaUtils.electronIdVal(eleIdCutsEle[i], bits::eleIdLoose);
      if(eleId == false) continue;
      goodEle.push_back(i);
    }

    if(int(goodEle.size()) < 2 ) continue; // no Ele pair found
    if(nj > 0) Npassed[1]+= weight;;
    if(nj > 1) Npassed_2j[1] += weight;
    if(nj > 2) Npassed_3j[1] += weight;
    if(nj > 3) Npassed_4j[1] += weight;

    ///////////////////
    // Z selection
    ///////////////////

    // make Z candidates;
    for(int i=0; i<int(goodEle.size()); i++) {
      TLorentzVector Ele1(pxEle[goodEle[i]], pyEle[goodEle[i]], pzEle[goodEle[i]], energyEle[goodEle[i]]);
      for(int j=i; j<int(goodEle.size()); j++) {
        TLorentzVector Ele2(pxEle[goodEle[j]], pyEle[goodEle[j]], pzEle[goodEle[j]], energyEle[goodEle[j]]);
        TLorentzVector Z = Ele1+Ele2;
        if(Z.M()> 0. && Z.M()<10000. && chargeEle[i]*chargeEle[j]<0 ) {
          Zcand.push_back(Z);
          if(Ele1.Pt() > Ele2.Pt()) {
            iZDaugh1.push_back(goodEle[i]);
            iZDaugh2.push_back(goodEle[j]);
          } else {
	    iZDaugh1.push_back(goodEle[j]);
            iZDaugh2.push_back(goodEle[i]);
          }
        }
      }
    }

    for(int i=0; i<int(Zcand.size()); i++) {
      int iEle1 = iZDaugh1[i];
      int iEle2 = iZDaugh2[i];
      if(Zcand[i].M() < 60. || Zcand[i].M() > 120. || // Z mass window
	 pTEle(iEle1) < 20.)  {  // Pt 1st leg
	EraseZ(i);
      }
    }
    
    if(int(Zcand.size()) == 0) continue; // no Z found
    if(nj > 0) Npassed[2]+= weight;
    if(nj > 1) Npassed_2j[2] += weight;
    if(nj > 2) Npassed_3j[2] += weight;
    if(nj > 3) Npassed_4j[2] += weight;
    
    //////////////////////////////////////
    // Electron |dz| cut [PV dependent]
    //////////////////////////////////////
    for(int i=0; i<Zcand.size(); i++) {
      int iEle1 = iZDaugh1[i];
      int iEle2 = iZDaugh2[i];
      bool passed = true;
      if(CalcDzPV(iEle1, BestPV(i)) > 0.12) passed = false; // |dz| cut for 1st leg
      if(CalcDzPV(iEle2, BestPV(i)) > 0.12) passed = false; // |dz| cut for 2nd leg 
      if(!passed) EraseZ(i);
    }

    if(int(Zcand.size()) == 0) continue;
    if(nj > 0) Npassed[3]+= weight;
    if(nj > 1) Npassed_2j[3] += weight;
    if(nj > 2) Npassed_3j[3] += weight;
    if(nj > 3) Npassed_4j[3] += weight;

    //////////////////////////////////////
    // Electron Dxy cut
    //////////////////////////////////////
    for(int i=0; int(i<Zcand.size()); i++) {
      int iEle1 = iZDaugh1[i];
      int iEle2 = iZDaugh2[i];
      // Dxy significance cut
      bool passed = true;
      if(CalcDxyPV(iEle1, BestPV(i)) > 0.04) passed = false; // |dz| cut for 1st leg        
      if(CalcDxyPV(iEle2, BestPV(i)) > 0.04) passed = false; // |dz| cut for 2nd leg        
      if(!passed) EraseZ(i);
    }

    if(int(Zcand.size()) == 0) continue;
    if(nj > 0) Npassed[4]+= weight;
    if(nj > 1) Npassed_2j[4] += weight;
    if(nj > 2) Npassed_3j[4] += weight;
    if(nj > 3) Npassed_4j[4] += weight;

    //////////////////////////////////////
    // Muon isolation [PV dependent]
    //////////////////////////////////////
    for(int i=0; i<int(Zcand.size()); i++) {
      int iEle1 = iZDaugh1[i];
      int iEle2 = iZDaugh2[i];
      bool passed = true;
      if(SumPt(iEle1, i)) passed = false; // Isolation on 1st leg
      if(SumPt(iEle2, i)) passed = false; // Isolation on 2nd leg
      if(!passed) EraseZ(i);
    }

    if(int(Zcand.size()) == 0) continue;
    if(nj > 0) Npassed[5]+= weight;
    if(nj > 1) Npassed_2j[5] += weight;
    if(nj > 2) Npassed_3j[5] += weight;
    if(nj > 3) Npassed_4j[5] += weight;

    // best candidate selection
    double highestpt1 = 0.;
    double highestpt2 = 0.;
    int iEhighestpt1 = -99;
    int iEhighestpt2 = -99;

    int ibestZ = -99;
    for(int i=0; i< int(Zcand.size()); i++) {
      double pt1 = pTEle(iZDaugh1[i]);
      double pt2 = pTEle(iZDaugh2[i]);
      if(pt1 > highestpt1 || // the highest pT is larger
	 ( iEhighestpt1 ==  iZDaugh1[i] &&  pt2 >= highestpt2) ) { // the second highest pT is larger
	highestpt1 = pt1;
	highestpt2 = pt2;
	iEhighestpt1 = iZDaugh1[i];
	iEhighestpt2 = iZDaugh2[i];
	ibestZ = i;
      }
    }

    // erase the other candidates
    for(int i=0; i<int(Zcand.size()); i++) if(i != ibestZ) EraseZ(i);
    int iBestPV = BestPV(0);

    // make jets    
    double zPV = PVzPV[iBestPV];
    vector<CaloTower> calotowers = CreateCaloTowers(thresh, zPV, 3);
    
    int nfound = 0;
    int iCT = 0;
    while(nfound == 0 && iCT< int(calotowers.size())) {
      if(calotowers[iCT].et() > 0.5) nfound = 1;
      iCT++;
    }
    
    vector<Jet> jet_SIS_calo;
    
    if(nfound == 1) {
      //SISCone = jet_SIS_calo
      jet_SIS_calo = SortJet(SISCone(calotowers, 0.5, 0.0));
    }
    
    // Make track jets
    vector<CaloTower> track_collection;
    TVector3 vPV(PVxPV[iBestPV],PVyPV[iBestPV],PVzPV[iBestPV]);
    for(int i = 0; i < nTrack; i++){
      TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
      vT = vT-vPV;
      //remove reco Z muons from track collection
      int iEle1 = iZDaugh1[i];
      int iEle2 = iZDaugh2[i];
      if(i == trackIndexEle[iEle1]) continue;
      if(i == trackIndexEle[iEle2]) continue;
      if(fabs(trackDxyTrack[i]) > 0.04) continue;
      if(fabs(vT.z()) > 0.1) continue;
      if(fabs(vT.Mag()) > 0.1) continue;
      if(trackNormalizedChi2Track[i] > 20.0) continue;
      if(trackValidHitsTrack[i] < 5) continue;
      TVector3 v(pxTrack[i], pyTrack[i], pzTrack[i]);
      double pt = v.Pt();
      if(pt < 0.5) continue;
      if(pt > 500.) continue;
      if(fabs(v.Eta()) > 2.4) continue;
      track_collection.push_back(CaloTower(pt*cosh(v.Eta()), 0.0, v, v, v));
    }

    nfound = 0;
    iCT = 0;
    while(nfound == 0 && iCT < int(track_collection.size())) {
      if(track_collection[iCT].et() > 0.5) nfound = 1;
      iCT++;
    }
    
    vector<Jet> jet_SIS_track;
    
    if(nfound == 1) {
      //SISCone = jet_SIS_track
      jet_SIS_track = SortJet(SISCone(track_collection, 0.5, 0.0));
    }

    // compute the event shape variables

    CoolTools* CT = new CoolTools();
    double Ztruemass = 91.1876;
    TLorentzVector ZfourVect(Zcand[0].Px(), Zcand[0].Py(), Zcand[0].Pz(), sqrt(pow(Zcand[0].P(),2.)+pow(Ztruemass,2.)));
    TVector3 b = ZfourVect.BoostVector();
    b.SetZ(0.0);

    vector<TLorentzVector> calovectors = CT->Get4Vectors(CT->BoostJets(CT->CaloTowers2Jets(calotowers, 0),b));
    TLorentzVector e1;
    TLorentzVector e2;
    e1.SetPxPyPzE(-pxEle[iZDaugh1[0]],
		  -pyEle[iZDaugh1[0]],
		  -pzEle[iZDaugh1[0]],
		  energyEle[iZDaugh1[0]]);
    e2.SetPxPyPzE(-pxEle[iZDaugh2[0]],
		  -pyEle[iZDaugh2[0]],
		  -pzEle[iZDaugh2[0]],
		  energyEle[iZDaugh2[0]]);
    e1.Boost(-b); 
    e2.Boost(-b); 
    calovectors.push_back(e1);
    calovectors.push_back(e2);

    CT->CalcTSphericity(calovectors);

    TLorentzVector MHT = CreateMET(calovectors).metVector();

    TLorentzVector Tmet    = CreateMET(calotowers).metVector();
    TLorentzVector TmetDEF = CreateMET(calotowers).metVector();

    double phiMET2 = atan2(Tmet.Py()+Zcand[0].Py(), Tmet.Px()+Zcand[0].Px());

    MHTphiMET = DeltaPhi(phiMET2,MHT.Phi());
    TSphiMET  = DeltaPhi(phiMET2,CT->TSphericityAxis().Phi());
    // we need the calojet to be there. This is a problem
    // for nTrackJet >0  and nCaloJet = 0 
    if(jet_SIS_calo.size()> 0) {
      MHTphiJet = DeltaPhi(jet_SIS_calo[0].Phi(),MHT.Phi());
    } else {
      MHTphiJet = 1.;
    }
    calotowers.clear();
    
    sinMHTphiMET = sin(MHTphiMET);
    sinTSphiMET = sin(TSphiMET);
    if(fabs(sinTSphiMET) > 0.85) continue;
    if(nj > 0) Npassed[6]+= weight;
    if(nj > 1) Npassed_2j[6] += weight;
    if(nj > 2) Npassed_3j[6] += weight;
    if(nj > 3) Npassed_4j[6] += weight;
    
    // Calo Jet counting
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double e = 0.;
    int nJetsCalo = 0;
    double eta_cut = 3.0;
    for(int i=0; i<int(jet_SIS_calo.size()); i++) 
      if(IsItTheElectron(jet_SIS_calo[i], iZDaugh1[0]) == false &&
	 IsItTheElectron(jet_SIS_calo[i], iZDaugh2[0]) == false) 
	if(fabs(jet_SIS_calo[i].eta()) <eta_cut && jet_SIS_calo[i].pt() > jet_et_cut) {
	  HeleIDHardestEle->Fill(jet_SIS_calo[i].EmFrac());
	  nJetsCalo++;
	  x += jet_SIS_calo[i].px();
	  y += jet_SIS_calo[i].py();
	  z += jet_SIS_calo[i].pz();
	  e += jet_SIS_calo[i].e();
	}
  
    sumPtJet = sqrt(x*x+y*y);

    mJet = TLorentzVector(x,y,z,e).M();
    
    TLorentzVector zpt(ZfourVect.X(), ZfourVect.Y(), 0., ZfourVect.E());

    TLorentzVector j(0,0,0,0);
    if(nJetsCalo>=1) {
      j+= TLorentzVector(jet_SIS_calo[0].px(), jet_SIS_calo[0].py(), jet_SIS_calo[0].pz(), jet_SIS_calo[0].e());
      sinJMET1 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), phiMET2);
      //      sinJMET1 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), zpt.Vect().Phi());
    } else sinJMET2 = -90.;

    if(nJetsCalo>=2) {
      j+= TLorentzVector(jet_SIS_calo[1].px(), jet_SIS_calo[1].py(), jet_SIS_calo[1].pz(), jet_SIS_calo[1].e());
      sinJMET2 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), phiMET2);
      //      sinJMET1 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), zpt.Vect().Phi());
    } else sinJMET2 = -90.;

    if(nJetsCalo>=3) {
      j+= TLorentzVector(jet_SIS_calo[2].px(), jet_SIS_calo[2].py(), jet_SIS_calo[2].pz(), jet_SIS_calo[2].e());
      sinJMET3 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), phiMET2);
      //      sinJMET1 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), zpt.Vect().Phi());
    } else sinJMET3 = -90.;

    if(nJetsCalo>=4) {
      j+= TLorentzVector(jet_SIS_calo[3].px(), jet_SIS_calo[3].py(), jet_SIS_calo[3].pz(), jet_SIS_calo[3].e());
      sinJMET4 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), phiMET2);
      //      sinJMET1 = DeltaPhi(TLorentzVector(j.X(),j.Y(),0.,j.E()).Vect().Phi(), zpt.Vect().Phi());
    } else sinJMET4 = -90.;


    // Track jet counting
    int nJetsTrack = 0;
    for(int i=0; i<int(jet_SIS_track.size()); i++)
      if(fabs(jet_SIS_track[i].eta()) <2.4 &&
      //      if(fabs(jet_SIS_track[i].eta()) <eta_cut &&
	 jet_SIS_track[i].pt() > jet_pt_cut) 
	nJetsTrack++;
    
    if(nJetsCalo > 0) Npassed[7]+= weight;
    if(nJetsCalo > 1) Npassed_2j[7] += weight;
    if(nJetsCalo > 2) Npassed_3j[7] += weight;
    if(nJetsCalo > 3) Npassed_4j[7] += weight;

    if(nJetsTrack > 0) Npassed[8]+= weight;
    if(nJetsTrack > 1) Npassed_2j[8] += weight;
    if(nJetsTrack > 2) Npassed_3j[8] += weight;
    if(nJetsTrack > 3) Npassed_4j[8] += weight;

    // write the tree
    njets = double(nJetsCalo);
    ntrackjets = double(nJetsTrack);
    mll = Zcand[0].M();
    int iEle1 = iZDaugh1[0];
    int iEle2 = iZDaugh2[0];
    etaEle1 = etaEle[iEle1];
    etaEle2 = etaEle[iEle2];
    phiEle1 = phiEle[iEle1];
    phiEle2 = phiEle[iEle2];
    ptEle1 = pTEle(iEle1);
    ptEle2 = pTEle(iEle2);
    Zpt = sqrt(pow(ptEle1*cos(phiEle1)+ptEle2*cos(phiEle2),2.)+pow(ptEle1*sin(phiEle1)+ptEle2*sin(phiEle2),2.));
    MET = TmetDEF.Et();
    phiMET = TmetDEF.Phi();
    if(njets > 0) {
      ptJet1 = jet_SIS_calo[0].et();
      etaJet1 = jet_SIS_calo[0].eta();
      phiJet1 = jet_SIS_calo[0].phi();
      emfJet1 = jet_SIS_calo[0].EmFrac();
    } else {
      ptJet1 = 31.; // default dummy variable for sure in the range
      etaJet1 = 1.; // default dummy variable for sure in the range
      phiJet1 = 1.; // default dummy variable for sure in the range
      emfJet1 = 1.;;
    }
    if(njets > 1) {
      ptJet2 = jet_SIS_calo[1].et();
      etaJet2 = jet_SIS_calo[1].eta();
      phiJet2 = jet_SIS_calo[1].phi();
      emfJet2 = jet_SIS_calo[1].EmFrac();
    } else {
      ptJet2 = 31.; // default dummy variable for sure in the range
      etaJet2 = 1.; // default dummy variable for sure in the range
      phiJet2 = 1.; // default dummy variable for sure in the range
      emfJet2 = 1.;
    }      
    if(njets > 2) {
      ptJet3 = jet_SIS_calo[2].et();
      etaJet3 = jet_SIS_calo[2].eta();
      phiJet3 = jet_SIS_calo[2].phi();
      emfJet3 = jet_SIS_calo[2].EmFrac();
    } else {
      ptJet3 = 31.; // default dummy variable for sure in the range
      etaJet3 = 1.; // default dummy variable for sure in the range
      phiJet3 = 1.; // default dummy variable for sure in the range
      emfJet3 = 1.;
    }
    sumPtOverPtEle1 = SumPt(iEle1, 0)/pTEle(iEle1);
    sumPtOverPtEle2 = SumPt(iEle2, 0)/pTEle(iEle2);

    genMET    = etGenMet[0];
    genphiMET = phiGenMet[0];

    // weight  is global

    // fill trigger bits

    HLT_Ele10_SW_L1R =  double(firedTrg[requiredTriggers[0]]);
    HLT_Ele15_SW_L1R =  double(firedTrg[requiredTriggers[1]]);
    HLT_Ele15_LW_L1R =  double(firedTrg[requiredTriggers[2]]);
    HLT_IsoEle15_L1I =  double(firedTrg[requiredTriggers[3]]);
    HLT_IsoEle18_L1R =  double(firedTrg[requiredTriggers[4]]);
    HLT_IsoEle15_LW_L1I =  double(firedTrg[requiredTriggers[5]]);
    HLT_LooseIsoEle15_LW_L1R =  double(firedTrg[requiredTriggers[6]]);

    outTree->Fill();

    //make the plots
    for(int i=0; i<jet_SIS_calo.size(); i++) {
      int nTkgood =0;
      int nTk =0;
      for(int j=0; j<nTrack; j++) {
	TVector3 v(pxTrack[j], pyTrack[j], pzTrack[j]);
	if(DeltaR(jet_SIS_calo[i].eta(), jet_SIS_calo[i].phi(), v.Eta(), v.Phi()) < 0.5) {
	  TVector3 vT(trackVxTrack[j],trackVyTrack[j],trackVzTrack[j]);
	  vT = vT-vPV;
	  dztrack->Fill(vT.z());
	  dxytrack->Fill(trackDxyPVTrack[j]);
	  nTk++;
	  double pt = sqrt(pxTrack[j]*pxTrack[j]+pyTrack[j]*pyTrack[j]);
	  if(fabs(vT.z()) < 0.03&&
	     trackDxyPVTrack[j] < 0.02 &&
	     pt > 0.5 &&
	     pt < 500.)
	    nTkgood++;
	}
      }
      aveNtrack->Fill(nTkgood*1.0001);
      fNtrack->Fill(nTkgood*1./nTk);
    }
  }

  TTree* nevTree = new TTree("nevTree", "nevTree");
  double nEv_nocuts    = Npassed[0];
  double nEv_MuonPVind = Npassed[1];
  double nEv_Zsel      = Npassed[2];
  double nEv_MuonDz    = Npassed[3];
  double nEv_Electron  = Npassed[4];
  double nEv_OnlyOneZ  = Npassed[5];
  double nEv_sinphi    = Npassed[6];
  double nEv_1calojet  = Npassed[7];
  double nEv_1trackjet = Npassed[8];

  double nEv_nocuts_2j    = Npassed_2j[0];
  double nEv_MuonPVind_2j = Npassed_2j[1];
  double nEv_Zsel_2j      = Npassed_2j[2];
  double nEv_MuonDz_2j    = Npassed_2j[3];
  double nEv_Electron_2j  = Npassed_2j[4];
  double nEv_OnlyOneZ_2j  = Npassed_2j[5];
  double nEv_sinphi_2j    = Npassed_2j[6];
  double nEv_1calojet_2j  = Npassed_2j[7];
  double nEv_1trackjet_2j = Npassed_2j[8];

  double nEv_nocuts_3j    = Npassed_3j[0];
  double nEv_MuonPVind_3j = Npassed_3j[1];
  double nEv_Zsel_3j      = Npassed_3j[2];
  double nEv_MuonDz_3j    = Npassed_3j[3];
  double nEv_Electron_3j  = Npassed_3j[4];
  double nEv_OnlyOneZ_3j  = Npassed_3j[5];
  double nEv_sinphi_3j    = Npassed_3j[6];
  double nEv_1calojet_3j  = Npassed_3j[7];
  double nEv_1trackjet_3j = Npassed_3j[8];

  double nEv_nocuts_4j    = Npassed_4j[0];
  double nEv_MuonPVind_4j = Npassed_4j[1];
  double nEv_Zsel_4j      = Npassed_4j[2];
  double nEv_MuonDz_4j    = Npassed_4j[3];
  double nEv_Electron_4j  = Npassed_4j[4];
  double nEv_OnlyOneZ_4j  = Npassed_4j[5];
  double nEv_sinphi_4j    = Npassed_4j[6];
  double nEv_1calojet_4j  = Npassed_4j[7];
  double nEv_1trackjet_4j = Npassed_4j[8];


  nevTree->Branch("nEv_nocuts",       &nEv_nocuts,    "nEv_nocuts/D");
  nevTree->Branch("nEv_MuonPVind",    &nEv_MuonPVind,    "nEv_MuonPVind/D");
  nevTree->Branch("nEv_Zsel",         &nEv_Zsel,    "nEv_Zsel/D");
  nevTree->Branch("nEv_MuonDz",       &nEv_MuonDz,    "nEv_MuonDz/D");
  nevTree->Branch("nEv_Electron",     &nEv_Electron,    "nEv_Electron/D");
  nevTree->Branch("nEv_OnlyOneZ",     &nEv_OnlyOneZ,    "nEv_OnlyOneZ/D");
  nevTree->Branch("nEv_sinphi",       &nEv_sinphi,    "nEv_sinphi/D");
  nevTree->Branch("nEv_1calojet",     &nEv_1calojet,    "nEv_1calojet/D");
  nevTree->Branch("nEv_1trackjet",    &nEv_1trackjet,    "nEv_1trackjet/D");

  nevTree->Branch("nEv_nocuts_2j",       &nEv_nocuts_2j,    "nEv_nocuts_2j/D");
  nevTree->Branch("nEv_MuonPVind_2j",    &nEv_MuonPVind_2j, "nEv_MuonPVind_2j/D");
  nevTree->Branch("nEv_Zsel_2j",         &nEv_Zsel_2j,      "nEv_Zsel_2j/D");
  nevTree->Branch("nEv_MuonDz_2j",       &nEv_MuonDz_2j,    "nEv_MuonDz_2j/D");
  nevTree->Branch("nEv_Electron_2j",     &nEv_Electron_2j,  "nEv_Electron_2j/D");
  nevTree->Branch("nEv_OnlyOneZ_2j",     &nEv_OnlyOneZ_2j,  "nEv_OnlyOneZ_2j/D");
  nevTree->Branch("nEv_sinphi_2j",       &nEv_sinphi_2j,    "nEv_sinphi_2j/D");
  nevTree->Branch("nEv_1calojet_2j",     &nEv_1calojet_2j,  "nEv_1calojet_2j/D");
  nevTree->Branch("nEv_1trackjet_2j",    &nEv_1trackjet_2j, "nEv_1trackjet_2j/D");

  nevTree->Branch("nEv_nocuts_3j",       &nEv_nocuts_3j,    "nEv_nocuts_3j/D");
  nevTree->Branch("nEv_MuonPVind_3j",    &nEv_MuonPVind_3j, "nEv_MuonPVind_3j/D");
  nevTree->Branch("nEv_Zsel_3j",         &nEv_Zsel_3j,      "nEv_Zsel_3j/D");
  nevTree->Branch("nEv_MuonDz_3j",       &nEv_MuonDz_3j,    "nEv_MuonDz_3j/D");
  nevTree->Branch("nEv_Electron_3j",     &nEv_Electron_3j,  "nEv_Electron_3j/D");
  nevTree->Branch("nEv_OnlyOneZ_3j",     &nEv_OnlyOneZ_3j,  "nEv_OnlyOneZ_3j/D");
  nevTree->Branch("nEv_sinphi_3j",       &nEv_sinphi_3j,    "nEv_sinphi_3j/D");
  nevTree->Branch("nEv_1calojet_3j",     &nEv_1calojet_3j,  "nEv_1calojet_3j/D");
  nevTree->Branch("nEv_1trackjet_3j",    &nEv_1trackjet_3j, "nEv_1trackjet_3j/D");

  nevTree->Branch("nEv_nocuts_4j",       &nEv_nocuts_4j,    "nEv_nocuts_4j/D");
  nevTree->Branch("nEv_MuonPVind_4j",    &nEv_MuonPVind_4j, "nEv_MuonPVind_4j/D");
  nevTree->Branch("nEv_Zsel_4j",         &nEv_Zsel_4j,      "nEv_Zsel_4j/D");
  nevTree->Branch("nEv_MuonDz_4j",       &nEv_MuonDz_4j,    "nEv_MuonDz_4j/D");
  nevTree->Branch("nEv_Electron_4j",     &nEv_Electron_4j,  "nEv_Electron_4j/D");
  nevTree->Branch("nEv_OnlyOneZ_4j",     &nEv_OnlyOneZ_4j,  "nEv_OnlyOneZ_4j/D");
  nevTree->Branch("nEv_sinphi_4j",       &nEv_sinphi_4j,    "nEv_sinphi_4j/D");
  nevTree->Branch("nEv_1calojet_4j",     &nEv_1calojet_4j,  "nEv_1calojet_4j/D");
  nevTree->Branch("nEv_1trackjet_4j",    &nEv_1trackjet_4j, "nEv_1trackjet_4j/D");

  nevTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  HeleIDHardestEle->Write();
  aveNtrack->Write();
  fNtrack->Write();
  dztrack->Write();
  dxytrack->Write();  
  outTree->Write();
  nevTree->Write();
  file->Close();
  
}

int CandleCalibee::BestPV(int bestzIdx) {
  // the best primary vertex is the vtx associated to the hardest muon 
  return vtxIndexTrack[trackIndexEle[iZDaugh1[bestzIdx]]];
}


 double CandleCalibee::pTEle(int i) {
   return sqrt(pxEle[i]*pxEle[i]+pyEle[i]*pyEle[i]);
 }

double CandleCalibee::SumPt(int iEle, int iZ) {
   double eta0 = etaEle[iEle];
  double phi0 = phiEle[iEle];
  double sumPt_tmp = 0;
  for(int i=0; i< nTrack; i++) {
    if(i == trackIndexEle[iEle]) continue; // take out the electron
    if(trackValidHitsTrack[i] <5) continue;                                     // minimum number of hits  XXX    
    if(fabs(trackDxyTrack[i]/trackDxyErrorTrack[i]) > 5.) continue;    // track incompatible with the vertex on (x,y) 
    if(fabs(PVzPV[BestPV(iZ)]-trackVzTrack[i]) > 0.1) continue;              // track incompatible with the vertex on z 
    TVector3 v(pxTrack[i], pyTrack[i], pzTrack[i]);
    if(sqrt(pow(v.Eta()-eta0,2.)+pow(v.Phi()-phi0,2.)) > 0.5) continue; // track outside the cone 
    if(v.Pt() < 0.500) continue;     // minimum pT 
    if(v.Pt() > 500.) continue;     // maximum pT 
    sumPt_tmp += v.Pt();
  }
  return sumPt_tmp; // -1 because the muon pT was not explicitely taken out
}

double CandleCalibee::DeltaPhi_PiHalf(double phi1, double phi2) {
  double dp = fabs(DeltaPhi(phi1, phi2));
  if(dp > asin(1.)) 
    dp = asin(1.)*2. - dp;
  return dp;
}

double CandleCalibee::CalcDxyPV(int iEle, int iPV) {
  int iTrack = trackIndexEle[iEle];
  return eleDxyPV(PVxPV[iPV], PVyPV[iPV], PVzPV[iPV], 
		  trackVxTrack[iTrack], trackVyTrack[iTrack], trackVzTrack[iTrack],
		  pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack]);
}

double CandleCalibee::CalcErrDxyPV(int iEle, int iPV) {
  // to fix
  return trackDxyErrorTrack[trackIndexEle[iEle]];
}

double CandleCalibee::CalcDzPV(int iEle, int iPV) {
  int iTrack = trackIndexEle[iEle];
  int iTrackPV = vtxIndexTrack[iTrack];
  return eleDszPV(PVxPV[iPV], PVyPV[iPV], PVzPV[iPV], 
		  trackVxTrack[iTrack], trackVyTrack[iTrack], trackVzTrack[iTrack],
		  pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack]);
}

double CandleCalibee::CalcErrDzPV(int iEle, int iPV) {
  return trackDzErrorTrack[trackIndexEle[iEle]];
}

void CandleCalibee::EraseZ(int iZ) {
  Zcand.erase(Zcand.begin()+iZ);
  iZDaugh1.erase(iZDaugh1.begin()+iZ);
  iZDaugh2.erase(iZDaugh2.begin()+iZ);
}

bool CandleCalibee::IsItTheElectron(Jet jet, int iEle) {
  bool itis = false;
  //  if((DeltaR(jet.eta(), jet.phi(), double(etaEle[iEle]), double(phiEle[iEle])) < 0.3) ||                                                                                                                                               
  //     (jet.EmFrac() > 0.95)) itis = true;                                                                                                                                                                                               
  if(DeltaR(jet.eta(), jet.phi(), double(etaEle[iEle]), double(phiEle[iEle])) < 0.3) itis = true;
  return itis;
}


