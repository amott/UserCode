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
#include "RazorDMAnalysis.hh"

Float_t Jet_Min_Pt = 70;

RazorDMAnalysis::RazorDMAnalysis(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

RazorDMAnalysis::RazorDMAnalysis(TTree *tree, string jsonFile,bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
  }
}

RazorDMAnalysis::~RazorDMAnalysis() {}

void RazorDMAnalysis::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorDMAnalysis::SetWeight(double weight){
  _weight=weight;
}

void RazorDMAnalysis::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  int HLT_Razor;
  int HLT_Razor_prescaled;

  bool ECALTPFilterFlag;
  bool drBoundary;
  bool drDead;
  bool CSCHaloFilterFlag;
  bool trackerFailureFilterFlag;
  bool BEECALFlag; 

  // PF block
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double R;
  double RSQ;
  double MR;
  double MRT;
  double pTHem1;
  double etaHem1;
  double phiHem1;
  double pTHem2;
  double etaHem2;
  double phiHem2;
  int    nBtag;
  int    nBtagTight;
  double W=_weight;
  double mst, mchi;
  int    BOX_NUM;
  Double_t Mu_Px_[2], Mu_Py_[2], Mu_Pz_[2], Mu_E_[2];
  //int    ss;
  int    nPV;
  int Jet_Multiplicity;

  // gen level info
  double pT1, pT2, eta1, eta2, phi1, phi2;
  int i1, i2;
  // ttbar decay: 0 = nolep, 1 = semilep; 2 = fully lep
  int iTopDecay;

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
  //outTree->Branch("ss", &ss, "ss/I");

  // HLT bits
  outTree->Branch("HLT_Razor", &HLT_Razor, "HLT_Razor/I"); 
  outTree->Branch("HLT_Razor_prescaled", &HLT_Razor_prescaled, "HLT_Razor_prescaled/I"); 

  //  block
  outTree->Branch("R", &R, "R/D");
  outTree->Branch("RSQ", &RSQ, "RSQ/D");
  outTree->Branch("MR", &MR, "MR/D");
  outTree->Branch("MRT", &MRT, "MRT/D");
  outTree->Branch("pTHem1", &pTHem1, "pTHem1/D");
  outTree->Branch("etaHem1", &etaHem1, "etaHem1/D");
  outTree->Branch("phiHem1", &phiHem1, "phiHem1/D");
  outTree->Branch("pTHem2", &pTHem2, "pTHem2/D");
  outTree->Branch("etaHem2", &etaHem2, "etaHem2/D");
  outTree->Branch("phiHem2", &phiHem2, "phiHem2/D");
  outTree->Branch("nBtag", &nBtag, "nBtag/I");
  outTree->Branch("nBtagTight", &nBtagTight, "nBtagTight/I");
  outTree->Branch("W", &W, "W/D");
  outTree->Branch("mst", &mst, "mst/D");
  outTree->Branch("mchi", &mchi, "mchi/D");
  outTree->Branch("nPV", &nPV, "nPV/I");

  outTree->Branch("i1", &i1, "i1/I");
  outTree->Branch("pT1", &pT1, "pT1/D");
  outTree->Branch("eta1", &eta1, "eta1/D");
  outTree->Branch("phi1", &phi1, "phi1/D");
  outTree->Branch("i2", &i2, "i2/I");
  outTree->Branch("pT2", &pT2, "pT2/D");
  outTree->Branch("eta2", &eta2, "eta2/D");
  outTree->Branch("phi2", &phi2, "phi2/D");

  outTree->Branch("N_Jets", &Jet_Multiplicity, "Jet_Multiplicity/I");

  outTree->Branch("Mu_Px",&Mu_Px_,"Mu_Px_[2]/D");
  outTree->Branch("Mu_Py",&Mu_Py_,"Mu_Py_[2]/D");
  outTree->Branch("Mu_Pz",&Mu_Pz_,"Mu_Pz_[2]/D");
  outTree->Branch("Mu_E",&Mu_E_,"Mu_E_[2]/D");

  outTree->Branch("iTopDecay", &iTopDecay, "iTopDecay/I");

  double Npassed_In = 0;
  double Npassed_PV = 0;
  //Jets
  double Npassed_2Jet = 0;
  //Leptons
  double Npassed_LepVeto=0;
  //B-tag
  double Npassed_0btag=0;

  double weightII = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  std::vector<std::string> maskHLT_Razor; 
  maskHLT_Razor.push_back("HLT_RsqMR55_Rsq0p09_MR150");
  maskHLT_Razor.push_back("HLT_RsqMR60_Rsq0p09_MR150");
  maskHLT_Razor.push_back("HLT_RsqMR65_Rsq0p09_MR150");
  
  std::vector<std::string> maskHLT_Razor_prescaled; 
  maskHLT_Razor_prescaled.push_back("HLT_RsqMR40_Rsq0p04");
  //  maskHLT_Razor_prescaled.push_back("HLT_RsqMR45_Rsq0p09");
  
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
      setRequiredTriggers(maskHLT_Razor); reloadTriggerMask(true); HLT_Razor = hasPassedHLT();
      setRequiredTriggers(maskHLT_Razor_prescaled); reloadTriggerMask(true); HLT_Razor_prescaled = hasPassedHLT();

      ECALTPFilterFlag = (METFlags >> 0)%2;
      drBoundary = (METFlags >> 1)%2;
      drDead = (METFlags >> 2)%2;
      CSCHaloFilterFlag = (METFlags >> 3)%2;
      trackerFailureFilterFlag = (METFlags >> 4)%2;
      BEECALFlag = (METFlags >> 5)%2; 
    }

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
    
    Npassed_In += weightII;
    
    double m0=9999, m12=9999, mc=9999;
     
    // to integrate with others
    /*
    if(!_isData && _isSMS){
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
    */
    mst=m0;
    mchi=mc;
    
    //HLT and Data Filter
    int passedHLT = HLT_Razor + HLT_Razor_prescaled;

    if (_isData==true) {
      if (passedHLT==0) continue;
      if ((ECALTPFilterFlag==0) || (drBoundary==0) || (drDead==0) || (CSCHaloFilterFlag==0) || (trackerFailureFilterFlag==0) || (BEECALFlag==0)) continue;
    }
    
    // find highest-pT PV [replace with Sagar's code]
    int iPV = passPV();
    if( iPV < 0 ) continue;//If negative no PV found to pass the cuts
    Npassed_PV += weightII;
    nPV = N_PV_EVENT;
    
    vector<TLorentzVector> Jet;
    vector <int> iJet;
    bool badjet = false;    
    for( int i = 0; i < nAK5Jet; i++ ) {
      TLorentzVector jet;
      double px = pxAK5Jet[i];
      double py = pyAK5Jet[i];
      double pz = pzAK5Jet[i];
      double E = sqrt( px*px + py*py + pz*pz );//massless Jet
      jet.SetPxPyPzE( px, py, pz, E );
      double pt = sqrt( px*px + py*py );
      if( fabs( jet.Eta() ) >= 3.0 ) continue;
      if( nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98 ){//What are this variables for???
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) < 2.55 && emFracAK5Jet[i] <= 0.01 ){
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) >= 2.55 && jet.Pt() > 80. && emFracAK5Jet[i] >= 1. ){
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) >= 2.55 && emFracAK5Jet[i] <= -0.9 ){
	badjet = true;
	break;
      }
      
      if ( jet.Pt() > 40. && fabs( jet.Eta() ) < 3.0 ) {
	Jet.push_back(jet);
	iJet.push_back(i);
      }
    }

    // Number of Jets
    Jet_Multiplicity = Jet.size();

    // jet ID
    if (badjet == true) continue;// If any Jet is bad (see loop before) event is rejected

    // >= 2Jets with pT> 70 GeV
    if( int( Jet.size() ) < 2 ) continue;//At least 2 Jets
    int iFirstJet = HighestPt(Jet, -99);
    int iSecondJet = HighestPt(Jet, iFirstJet);    
    if( Jet[iSecondJet].Pt() < 70. ) continue;//First and second most energetic Jets Must have Pt > 70 GeV
    Npassed_2Jet+=weightII;
    
    //count btagged jets
    nBtag = 0;
    nBtagTight = 0;
    
    for( int b = 0; b < iJet.size(); b++ ){
      if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5Jet[iJet[b]] ) ) nBtag++; //Loose	
      if( pfJetPassCSVT( combinedSecondaryVertexBJetTagsAK5Jet[ iJet[b] ] ) ) nBtagTight++; //Tight
    }
    //if(nBtag>0) continue;
    Npassed_0btag+=weightII;//What is the meaning of this variable
    
    TVector3 Muon_MET_Correction(0,0,0);

    // Boxes
    vector<int> iMuLoose;
    vector<TLorentzVector> MuLoose;

    for( int i = 0; i < nMuon; i++ ) {
      TLorentzVector thisMu( pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i] );
      if( ( isLooseMuon(i) ) && ( thisMu.Pt() > 15. ) ) {
	iMuLoose.push_back(i);
	MuLoose.push_back(thisMu);
      }
    }

    Mu_Px_[0] = Mu_Py_[0] = Mu_Pz_[0] = Mu_E_[0] = Mu_Px_[1] = Mu_Py_[1] = Mu_Pz_[1] = Mu_E_[1] =  -9999;

    // Correct MET for Presence of Loose Muons with pT>15 GeV
    int iFirstMuon = HighestPt(MuLoose, -99);
    if(iFirstMuon >= 0) {
      Muon_MET_Correction.SetX(Muon_MET_Correction.X() + pxMuon[iFirstMuon]);
      Muon_MET_Correction.SetY(Muon_MET_Correction.Y() + pyMuon[iFirstMuon]);

      Mu_Px_[0] = pxMuon[iFirstMuon];
      Mu_Py_[0] = pyMuon[iFirstMuon];
      Mu_Pz_[0] = pzMuon[iFirstMuon];
      Mu_E_[0] = energyMuon[iFirstMuon];
    }

    int iSecondMuon = HighestPt(MuLoose, iFirstMuon);
    if(iSecondMuon >= 0) {
      Muon_MET_Correction.SetX(Muon_MET_Correction.X() + pxMuon[iSecondMuon]); 
      Muon_MET_Correction.SetY(Muon_MET_Correction.Y() + pyMuon[iSecondMuon]);

      Mu_Px_[1] = pxMuon[iSecondMuon];
      Mu_Py_[1] = pyMuon[iSecondMuon];
      Mu_Pz_[1] = pzMuon[iSecondMuon];
      Mu_E_[1] = energyMuon[iSecondMuon];
    }

    // Tight Electron ID
    vector<int> iEleTight;
    for(int i=0; i<nEle; i++) {
      TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
      if(isTrigElectron(i) && thisEle.Pt() > 15.) {
	iEleTight.push_back(i);//Fill if the electron is tight and Pt > 15 GeV
      }
    }
    
    // Tight Tau ID
    vector<int> iTauTight;
    for ( int i = 0; i < nPFTau; i++) {
      TLorentzVector thisTau( pxPFTau[i], pyPFTau[i], pzPFTau[i], energyPFTau[i] );
      if ( isTightTau(i) && thisTau.Pt() > 20. ) {
	iTauTight.push_back(i);//Fill if the tau is tight and Pt > 20. GeV
      }
    }

    // Electron VETO
    if(iEleTight.size()>0) continue;

    // TAU VETO
    if (iTauTight.size()>0) continue;

    //This will give the number of events that pass until this point but not necesarily was recorded
    //This is important!! so "Npassed_LepVeto" and  "Npassed_2Jet" and the total number of recorded events
    //can be different. Ok (this is kind of arbitrary because it depends on the order it was done (probably problem is here))
    Npassed_LepVeto += weightII;//Record how many of th

    //cout << "NBTAG IS " << nBtag << endl;

    // BOX NUMBER
    BOX_NUM = -99;
    if(iMuLoose.size() == 0){
      BOX_NUM = 0;
    } else if(iMuLoose.size() == 1) {
      BOX_NUM = 1;
    } else if(iMuLoose.size() >= 2) {
      BOX_NUM = 2;
    }

    // dummy values                                          
    pTHem1 = -9999.;
    etaHem1 = -9999.;
    phiHem1 = -9999.;
    pTHem2 = -9999.;
    etaHem2 = -9999.;
    phiHem2 = -9999.;
    R = -99999.;
    RSQ = -99999.;
    MR = -99999.;
    
    // hemispheres
     vector<TLorentzVector> tmpJet = CombineJets(Jet);
     if( tmpJet.size() >= 2 ) {
       TLorentzVector Hem1 = tmpJet[0];
       TLorentzVector Hem2 = tmpJet[1];

       // PFMET
       TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);
       /*TVector3 MET(0,0, 0.);
       for(int i=0; i<Jet.size(); i++) {
	 MET = MET - TVector3(Jet[i].Px(), Jet[i].Py(), 0.);
	 }*/

       MRT = CalcMTR(Hem1, Hem2, MET + Muon_MET_Correction);
       //MRT = CalcMTR(Hem1, Hem2, MET);
       double variable = -999999.;
       double Rvariable = -999999.;
       variable = CalcGammaMRstar(Hem1, Hem2);
       if(variable >0) Rvariable = MRT/variable;
	 
       // fill the R and hem part of the output tree
       pTHem1 = Hem1.Pt();
       etaHem1 = Hem1.Eta();
       phiHem1 = Hem1.Phi();
       pTHem2 = Hem2.Pt();
       etaHem2 = Hem2.Eta();
       phiHem2 = Hem2.Phi();
       R = Rvariable;
       RSQ = Rvariable*Rvariable;
       MR = variable;
     }

     //gen-level info
     pT1 = -999;
     eta1 = -999;
     phi1 = -999;     
     pT2 = -999;
     eta2 = -999;
     phi2 = -999;
     i1 = -99;
     i2 = -99;
     if(!_isData) {
       iTopDecay = 0;
       int iL1 = -99;
       int iL2 = -99;
       for(int i=0; i<nMc; i++) {
	 // TT final state
	 if(abs(idMc[mothMc[i]]) == 24) {
	   if(idMc[i] >= 11 &&
	      idMc[i] <= 18) {
	     iTopDecay ++;
	   }
	 }
	 // Z daughters
	 if(idMc[mothMc[i]] == 23) {
	   if(iL1 <0) iL1 = i;
	   else if(iL2<0) iL2 = i;
	 }
       }
       iTopDecay = iTopDecay/2;
	
       if(iL1>=0) {
	 pT1 = pMc[iL1]*sin(thetaMc[iL1]);
	 eta1 = etaMc[iL1];
	 phi1 = phiMc[iL1];
	 i1 = idMc[iL1];
       } 
       if(iL2>=0) {
	 pT2 = pMc[iL2]*sin(thetaMc[iL2]);
	 eta2 = etaMc[iL2];
	 phi2 = phiMc[iL2];
	 i2 = idMc[iL2];
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
  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("Npassed_2Jet",  &Npassed_2Jet,  "Npassed_2Jet/D");
  effTree->Branch("Npassed_0btag",  &Npassed_0btag,  "Npassed_0btag/D");
  effTree->Branch("Npassed_LepVeto",  &Npassed_LepVeto,  "Npassed_LepVeto/D");

  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  effTree->Write();
  //  if(isSMS)FullSMSTree->Write();
  file->Close();
}
  
int RazorDMAnalysis::HighestPt(vector<TLorentzVector> p, int iHIGHEST) {

  int iH = -99;
  double highestPT = 0.;
  for(int i=0; i<p.size();i++) {
    if((p[i].Pt()>= highestPT) && (i != iHIGHEST)) {
      iH = i;
      highestPT = p[i].Pt();
    }
  }
  return iH;
}
