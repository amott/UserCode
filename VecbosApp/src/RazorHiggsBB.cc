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
#include "RazorHiggsBB.hh"

RazorHiggsBB::RazorHiggsBB(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
}

RazorHiggsBB::RazorHiggsBB(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

RazorHiggsBB::~RazorHiggsBB() {}

void RazorHiggsBB::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorHiggsBB::SetWeight(double weight) {
  _weight = weight;
}

void RazorHiggsBB::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  // prescaled RazorHiggsBB Triggers
  int HLT_R014_MR150;
  int HLT_R020_MR150;
  int HLT_R025_MR150;

  // hadronic razor triggers
  int HLT_R020_MR550;
  int HLT_R025_MR450;
  int HLT_R033_MR350;
  int HLT_R038_MR250;

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

  // general event info
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  int nPV;
  double W;

  // Higgs variables are global

  // prepare the output tree
  TTree* outTree[6]; 
  outTree[0] = new TTree("outTree2H2M", "outTree2H2M");
  outTree[1] = new TTree("outTree2H1M", "outTree2H1M");
  outTree[2] = new TTree("outTree2H0M", "outTree2H0M");
  outTree[3] = new TTree("outTree1H1M", "outTree1H1M");
  outTree[4] = new TTree("outTree1H0M", "outTree1H0M");
  outTree[5] = new TTree("outTree0H", "outTree0H");

  for(int i=0; i<6; i++) {

    outTree[i]->Branch("HLT_R020_MR550", &HLT_R020_MR550, "HLT_R020_MR550/I");
    outTree[i]->Branch("HLT_R025_MR450", &HLT_R025_MR450, "HLT_R025_MR450/I");
    outTree[i]->Branch("HLT_R033_MR350", &HLT_R033_MR350, "HLT_R033_MR350/I");
    outTree[i]->Branch("HLT_R038_MR250", &HLT_R038_MR250, "HLT_R038_MR250/I");

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
    
    if(i<5) {
      // first Higgs
      outTree[i]->Branch("PFH1Pt",   &PFH1Pt, "PFH1Pt/D");
      outTree[i]->Branch("PFH1Eta",  &PFH1Eta, "PFH1Eta/D");
      outTree[i]->Branch("PFH1Phi",  &PFH1Phi, "PFH1Phi/D");
      outTree[i]->Branch("PFH1Mass", &PFH1Mass, "PFH1Mass/D");
      outTree[i]->Branch("MergedH1", &MergedH1, "MergedH1/I");
    }

    if(i<3) {
      // second Higgs
      outTree[i]->Branch("PFH2Pt",   &PFH2Pt, "PFH2Pt/D");
      outTree[i]->Branch("PFH2Eta",  &PFH2Eta, "PFH2Eta/D");
      outTree[i]->Branch("PFH2Phi",  &PFH2Phi, "PFH2Phi/D");
      outTree[i]->Branch("PFH2Mass", &PFH2Mass, "PFH2Mass/D");
      outTree[i]->Branch("MergedH2", &MergedH2, "MergedH2/I");
    }
  }

  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // hadronic razor triggers
  std::vector<std::string> maskHLT_R020_MR550; maskHLT_R020_MR550.push_back("HLT_R020_MR550_v");
  std::vector<std::string> maskHLT_R025_MR450; maskHLT_R025_MR450.push_back("HLT_R025_MR450_v");
  std::vector<std::string> maskHLT_R033_MR350; maskHLT_R033_MR350.push_back("HLT_R033_MR350_v");
  std::vector<std::string> maskHLT_R038_MR250; maskHLT_R038_MR250.push_back("HLT_R038_MR250_v");

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
      setRequiredTriggers(maskHLT_R020_MR550); reloadTriggerMask(true); HLT_R020_MR550 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R025_MR450); reloadTriggerMask(true); HLT_R025_MR450 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R033_MR350); reloadTriggerMask(true); HLT_R033_MR350 = hasPassedHLT();
      setRequiredTriggers(maskHLT_R038_MR250); reloadTriggerMask(true); HLT_R038_MR250 = hasPassedHLT();
      
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
    vector<double> PFPUcorrJetBTAG; 
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      // to avoid messages of 0 pT
      if(sqrt(pow(pxAK5PFPUcorrJet[i],2.)+pow(pyAK5PFPUcorrJet[i],2.))<10.) continue;
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);
      if(myJet.Pt()>40. && fabs(etaAK5PFPUcorrJet[i])< 2.4) {
	PFPUcorrJet.push_back(myJet);
	PFPUcorrJetBTAG.push_back(trackCountingHighEffBJetTagsAK5Jet[i]);
      }
      // check if the event is good (from PFJets ID point of view)                                                              
      double EU = uncorrEnergyAK5PFPUcorrJet[i];
      // Apply jet correction first                                                                                             
      bool good_jet = false;
      if (myJet.Pt() > 40.0 && fabs(etaAK5PFPUcorrJet[i]) < 3.0) {
	double fHAD = (neutralHadronEnergyAK5PFPUcorrJet[i]+chargedHadronEnergyAK5PFPUcorrJet[i])/EU;
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
	}
      }
    }
    if(goodPFevent == false) continue;

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
    double mRmin = 0.0;  // to set after first plots
    if(PFPUcorrJet.size()<2) continue;
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

    // baseline Razor selection
    if(PFMR<RsqMin) continue;
    if(PFRsq<RsqMin) continue;

    ////////////////////////////////////////////////////////////
    // look for the highest four btag jet
    ////////////////////////////////////////////////////////////
    int iBJet[4];

    iBJet[0] = -99;
    double btagOne = 0.;
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(trackCountingHighEffBJetTagsAK5Jet[i] > btagOne) {
	btagOne = PFPUcorrJetBTAG[i];
	iBJet[0] = i;
      }
    }

    // look for the second highest btag
    iBJet[1] = -99;
    double btagTwo = 0;
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(i == iBJet[0]) continue;
      if(trackCountingHighEffBJetTagsAK5Jet[i] > btagTwo) {
	btagTwo = PFPUcorrJetBTAG[i];
	iBJet[1] = i;
      }
    }

    // look for the third highest btag
    iBJet[2] = -99;
    btagTwo = 0;
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(i == iBJet[0]) continue;
      if(i == iBJet[1]) continue;
      if(trackCountingHighEffBJetTagsAK5Jet[i] > btagTwo) {
	btagTwo = PFPUcorrJetBTAG[i];
	iBJet[2] = i;
      }
    }

    // look for the fourth highest btag
    iBJet[3] = -99;
    btagTwo = 0;
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(i == iBJet[0]) continue;
      if(i == iBJet[1]) continue;
      if(i == iBJet[2]) continue;
      if(trackCountingHighEffBJetTagsAK5Jet[i] > btagTwo) {
	btagTwo = PFPUcorrJetBTAG[i];
	iBJet[3] = i;
      }
    }

    ////////////////////////////////////////////////////////////
    // at least one BTag medium WP
    if(btagOne < 3.3) continue; 
    ////////////////////////////////////////////////////////////

    // Are any of the best four btagged jets merged Higgses?
    // by following the list order, we pick by default the most btag candidate
    double mergedMassMin = 70.; 
    vector<TLorentzVector> mergedHiggs;
    vector<double> mergedHiggsBTAG;
    vector<TLorentzVector> unmergedBJets;
    vector<double> unmergedBJetsBTAG;
    if(PFPUcorrJet[iBJet[0]].M() > mergedMassMin) {
      mergedHiggs.push_back(PFPUcorrJet[iBJet[0]]);
      mergedHiggsBTAG.push_back(PFPUcorrJetBTAG[iBJet[0]]);
    } else {
      unmergedBJets.push_back(PFPUcorrJet[iBJet[0]]);
      unmergedBJetsBTAG.push_back(PFPUcorrJetBTAG[iBJet[0]]);
    }
    if(PFPUcorrJet[iBJet[1]].M() > mergedMassMin) {
      mergedHiggs.push_back(PFPUcorrJet[iBJet[1]]);
      mergedHiggsBTAG.push_back(PFPUcorrJetBTAG[iBJet[1]]);      
    } else {
      unmergedBJets.push_back(PFPUcorrJet[iBJet[1]]);
      unmergedBJetsBTAG.push_back(PFPUcorrJetBTAG[iBJet[1]]);
    }
    if(iBJet[2] >= 0) {
      if(PFPUcorrJet[iBJet[2]].M() > mergedMassMin) {
	mergedHiggs.push_back(PFPUcorrJet[iBJet[2]]);
	mergedHiggsBTAG.push_back(PFPUcorrJetBTAG[iBJet[2]]);      
      } else {
	unmergedBJets.push_back(PFPUcorrJet[iBJet[2]]);
	unmergedBJetsBTAG.push_back(PFPUcorrJetBTAG[iBJet[2]]);
      }
    }
    if(iBJet[3] >= 0) {
      if(PFPUcorrJet[iBJet[3]].M() > mergedMassMin) {
	mergedHiggs.push_back(PFPUcorrJet[iBJet[3]]);
	mergedHiggsBTAG.push_back(PFPUcorrJetBTAG[iBJet[3]]);      
      } else {
	unmergedBJets.push_back(PFPUcorrJet[iBJet[3]]);
	unmergedBJetsBTAG.push_back(PFPUcorrJetBTAG[iBJet[3]]);
      }
    }

    // all the other jets
    vector<TLorentzVector> unmergedJets;
    vector<double> unmergedJetsBTAG;
    for(int i=0; i<PFPUcorrJet.size();i++) {
      if(i != iBJet[0] && i != iBJet[1] &&
	 i != iBJet[2] && i != iBJet[3]) {
	unmergedJets.push_back(PFPUcorrJet[i]);
	unmergedJetsBTAG.push_back(PFPUcorrJetBTAG[i]);
      }
    }

    ////////////////////////////////////////////
    // the dijet Higgs candidates
    ////////////////////////////////////////////

    vector<TLorentzVector> unmergedHiggs;
    double unmergedMassMin = 70.;

    // combine unmerged bjets into dijet Higgs candidates
    // - take the first entry in the list
    // - loop over the others
    // find the best-b pair with mass>threshold
    // REMOVE the two jets if found
    // REMOVE the first jet if not found
    // THIS NEEDS at least two jets left in the bjet list
    while(unmergedBJets.size() >= 2) {
      int i2 = -99;
      for(int i=1; i<unmergedBJets.size(); i++) {
	TLorentzVector Higgs = unmergedBJets[0]+unmergedBJets[i];
	// by following the list order we form the highest-btag combinations first
	if(Higgs.M() > unmergedMassMin and i2<0) {
	  i2 = i;
	}
      }
      if(i2>0) {
	// one combination found
	unmergedHiggs.push_back(unmergedBJets[0]+unmergedBJets[i2]);
	// remove the second leg 
	unmergedBJets.erase(unmergedBJets.begin()+i2);
	unmergedBJetsBTAG.erase(unmergedBJetsBTAG.begin()+i2);
      }
      // remove the first leg
      unmergedBJets.erase(unmergedBJets.begin());
      unmergedBJetsBTAG.erase(unmergedBJetsBTAG.begin());
    }

    // for each bjet left in the list
    // combine it with one non-bjet such that 
    // the mass is above threshold
    // pick the candidate with the highest sum of btags
    // in case more candidates are present
    for(int i1 = 0; i1< unmergedBJets.size(); i1++) {
      int iBest2 = -99;
      double bestBTAG = 0.;
      // take a second leg from non-bjet
      for(int i2=0; i2<unmergedJets.size();i2++) {
	// form the Higgs candidate
	TLorentzVector Higgs = unmergedJets[i2]+unmergedBJets[i1];
	if(Higgs.M() > unmergedMassMin && unmergedJetsBTAG[i2] > bestBTAG) {
	  bestBTAG = unmergedJetsBTAG[i2];
	  iBest2 = i2;
	}
      }
      // if a candidate was found
      if(iBest2>=0) {
	// - append it to the merged Higgs 
	unmergedHiggs.push_back(unmergedBJets[i1]+unmergedJets[iBest2]);
	//-  remove the second leg from the nobjet list
	unmergedBJets.erase(unmergedBJets.begin()+iBest2);
	unmergedBJetsBTAG.erase(unmergedBJetsBTAG.begin()+iBest2);
      }
    }

    // initialization
    PFH1Pt = -99.;
    PFH1Eta = -99.;
    PFH1Phi = -99.;
    PFH1Mass = -99.;
    MergedH1 = -99;
    PFH2Pt = -99.;
    PFH2Eta = -99.;
    PFH2Phi = -99.;
    PFH2Mass = -99.;
    MergedH2 = -99;

    if(mergedHiggs.size()>=2) {
      // 2H2Me box
      SetFirstHiggs(mergedHiggs[0], true);
      SetSecondHiggs(mergedHiggs[1], true);
      outTree[0]->Fill();
    } else if(mergedHiggs.size()>=1 and unmergedHiggs.size()>=1) {
      // 2H1Me box
      SetFirstHiggs(mergedHiggs[0], true);
      SetSecondHiggs(unmergedHiggs[0], false);
      outTree[1]->Fill();
    } else if(unmergedHiggs.size()>=2) {
      // 2H0Me box
      SetFirstHiggs(unmergedHiggs[0], false);
      SetSecondHiggs(unmergedHiggs[1], false);
      outTree[2]->Fill();
    } else if(mergedHiggs.size()>0) {
      // 1H1Me
      SetFirstHiggs(mergedHiggs[0], true);
      outTree[3]->Fill();
    } else if(unmergedHiggs.size()>0) {
      SetFirstHiggs(unmergedHiggs[0], false);
      outTree[4]->Fill();
    } else {
      // 0H
      outTree[5]->Fill();
    }
  }

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  for(int i=0; i<6; i++) 
    outTree[i]->Write();
  file->Close();
}
  
void RazorHiggsBB::SetFirstHiggs(TLorentzVector myHiggs, int merged) {
  PFH1Pt = myHiggs.Pt();
  PFH1Eta = myHiggs.Eta();
  PFH1Phi = myHiggs.Phi();
  PFH1Mass = myHiggs.M();
  MergedH1  = merged;
}

void RazorHiggsBB::SetSecondHiggs(TLorentzVector myHiggs, int merged) {
  PFH2Pt = myHiggs.Pt();
  PFH2Eta = myHiggs.Eta();
  PFH2Phi = myHiggs.Phi();
  PFH2Mass = myHiggs.M();
  MergedH2  = merged;
}

vector<TLorentzVector> RazorHiggsBB::CombineJets(vector<TLorentzVector> myjets){
  
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

