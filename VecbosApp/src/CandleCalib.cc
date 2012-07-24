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
#include "CandleCalib.hh"

CandleCalib::CandleCalib(TTree *tree) : Vecbos(tree) {}

CandleCalib::~CandleCalib(){}

vector<TH1D*> CandleCalib::CreateHistos1D(string dirname){
  vector<TH1D*> histos;
  string name;

  name = dirname+"_Zmass";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 140, 40., 110.));
  
  name = dirname+"_pTMu1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 1000.));
  
  name = dirname+"_pTMu2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 1000.));
  
  name = dirname+"_dxySigMu1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 40, 0., 10.));
  
  name = dirname+"_dxySigMu2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 40, 0., 10.));
  
  name = dirname+"_dzPVMu1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., .1));
  
  name = dirname+"_dzPVMu2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., .1));
  
  name = dirname+"_SumPtOPtMu1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 1.));
  
  name = dirname+"_SumPtOPtMu2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, 0., 1.));
  
  return histos;
}

void CandleCalib::FillHistos(vector<TH1D*> histos, int iZ) {
  int iMu1 = iZDaugh1[iZ];
  int iMu2 = iZDaugh2[iZ];
  int iTrack1 = trackIndexMuon[iMu1];
  int iTrack2 = trackIndexMuon[iMu2];
  int iPV = BestPV(iZ);
  histos[1]->Fill(pTMuon(iMu1),weight);
  histos[2]->Fill(pTMuon(iMu2),weight);
  histos[3]->Fill(transvImpactParTrack[iTrack1]/transvImpactParErrorTrack[iTrack1], weight);
  histos[4]->Fill(transvImpactParTrack[iTrack2]/transvImpactParErrorTrack[iTrack2], weight);
  histos[5]->Fill(impactPar3DTrack[iTrack1]/impactPar3DErrorTrack[iTrack1],weight);
  histos[6]->Fill(impactPar3DTrack[iTrack2]/impactPar3DErrorTrack[iTrack2],weight);
  histos[7]->Fill(SumPt(iMu1, iZ)/pTMuon(iMu1),weight);
  histos[8]->Fill(SumPt(iMu2, iZ)/pTMuon(iMu2),weight);
}

void CandleCalib::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");

  double njets;
  double ntrackjets;
  double npartons;
  double processID;
  double mll;
  double etaMu1;
  double etaMu2;
  double flagMu1;
  double flagMu2;
  double ptMu1;
  double ptMu2;
  double phiMu1;
  double phiMu2;
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
  double sumPtJet;
  double sumPtOverPtMu1;
  double sumPtOverPtMu2;
  double MHTphiMET;
  double TTphiMET;
  double TMphiMET;
  double MHTphiJet;
  double TTphiJet;
  double TMphiJet;
  double sinMHTphiMET;
  double HLT_Mu9;
  double HLT_Mu11;
  double HLT_Mu15;
  double HLT_DoubleMu3;

  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;

  //   double weight; THIS IS GLOBAL

  TH1D* HNcaloJets = new TH1D("NCaloJets", "NCaloJets", 6, 0., 6.);
  TH1D* HNtrackJets = new TH1D("NTrackJets", "NTrackJets", 6, 0., 6.);

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("njets",      &njets,      "njets/D");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("ntrackjets", &ntrackjets, "ntrackjets/D");
  outTree->Branch("npartons",   &npartons,   "npartons/D");
  outTree->Branch("processID",  &processID,   "processID/D");
  outTree->Branch("mll",        &mll,        "mll/D");
  outTree->Branch("etaMu1",     &etaMu1,     "etaMu1/D");
  outTree->Branch("etaMu2",     &etaMu2,     "etaMu2/D");
  outTree->Branch("flagMu1",     &flagMu1,     "flagMu1/D");
  outTree->Branch("flagMu2",     &flagMu2,     "flagMu2/D");
  outTree->Branch("ptMu1",      &ptMu1,      "ptMu1/D");
  outTree->Branch("ptMu2",      &ptMu2,      "ptMu2/D");
  outTree->Branch("phiMu1",     &phiMu1,     "phiMu1/D");
  outTree->Branch("phiMu2",     &phiMu2,     "phiMu2/D");
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
  outTree->Branch("sumPtJet",   &sumPtJet,   "sumPtJet/D");
  outTree->Branch("sumPtOverPtMu1",   &sumPtOverPtMu1,   "sumPtOverPtMu1/D");
  outTree->Branch("sumPtOverPtMu2",   &sumPtOverPtMu2,   "sumPtOverPtMu2/D");
  outTree->Branch("weight",     &weight,     "weight/D");
  outTree->Branch("MHTphiMET",  &MHTphiMET,  "MHTphiMET/D");
  outTree->Branch("TTphiMET",   &TTphiMET,   "TTphiMET/D");
  outTree->Branch("TMphiMET",   &TMphiMET,   "TTMphiMET/D");
  outTree->Branch("MHTphiJet",  &MHTphiJet,  "MHTphiJet/D");
  outTree->Branch("TTphiJet",   &TTphiJet,   "TTphiJet/D");
  outTree->Branch("TMphiJet",   &TMphiJet,   "TTMphiJet/D");
  outTree->Branch("sinMHTphiMET",   &sinMHTphiMET,   "sinMHTphiMET/D");
  outTree->Branch("HLT_Mu9", &HLT_Mu9, "HLT_Mu9/D");
  outTree->Branch("HLT_Mu11", &HLT_Mu11, "HLT_Mu11/D");
  outTree->Branch("HLT_Mu15", &HLT_Mu15, "HLT_Mu15/D");
  outTree->Branch("HLT_DoubleMu3", &HLT_DoubleMu3, "HLT_DoubleMu3/D");

  // prepare vectors for efficiency
  vector<string> cuts;
  vector<double> Npassed;
  
  cuts.push_back("Two Muons");
  cuts.push_back("Muon Selection [PV independent]");
  cuts.push_back("Z selection");
  cuts.push_back("Muon |dz| cut [PV dependent]");
  cuts.push_back("Muon isolation [PV dependent]");
  cuts.push_back("Only one Z");
  cuts.push_back("|sin(MHTphiMET)|<0.85");
  cuts.push_back("At least 1 calo jet");
  cuts.push_back("At least 1 track jet");

  for(int i=0; i< int(cuts.size())+1; i++) 
    Npassed.push_back(0.);

  int nCount =0;

  vector< vector<TH1D*> > Histos;
  Histos.push_back(CreateHistos1D("noCut"));
  for(int i=0; i< int(cuts.size()); i++) Histos.push_back(CreateHistos1D(cuts[i]));

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

    // check the alpgen ID
    //    if(AlpgenIdSelection(genAlpgenID, "Zjets") == false) continue;
    //    npartons = int(genAlpgenID)%1000;
    //    weight = genWeight; // times your k factor
    //    processID = genAlpgenID;
    
    //    if(nPV == 0) continue; // THIS IS A PATCH: we need a primary vertex

    // for ppMuX
    //    weight = 1.;
    //    npartons = 0;
    //    processID = 0.;
    
    // Did the Z decay to mu mu?
    // Look for the gen-level muon from Z
    //     int iGenmu = -99;
    //     for(int i=0; i<nMc; i++) 
    //       if(statusMc[i] == 3) // stable particle
    //   	if(abs(idMc[i]) == 13) // it is a muon
    //   	  if(abs(idMc[mothMc[i]]) == 23) // its mother is a Z
    //   	    iGenmu = i;
    //     if(iGenmu == -99) continue;
    
    //////////////////////////////
    // to cut at N Z events
    //////////////////////////////
    //     if(nCount == 1000) {
    //       jentry = stop;
    //       continue;
    //     }
    //////////////////////////////

    Npassed[0] += weight;
    nCount++;
    //for(int i=0; i<int(Zcand.size()); i++) FillHistos(Histos[0], i);
    
    // At least two reconstructed muons
    if(nMuon < 2) continue;
    Npassed[1]+= weight;
    // for(int i=0; i<int(Zcand.size()); i++) FillHistos(Histos[1], i);

    ///////////////////
    // muon selection
    ///////////////////

    vector<int> goodMuons;
    vector<int> flagMuons;
    Utils anaUtils;
    for(int i=0; i<nMuon; i++) {
      // is it global?
      bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[i], bits::AllGlobalMuons);
      flagMuons.push_back(int(isMuGlobal));
      //      if(!isMuGlobal) continue;      
      // number of valid hits
      int iTrack = trackIndexMuon[i];
      //      if(trackValidHitsTrack[iTrack] < 3) continue;    
      // pT cut
      //      if(pTMuon(i) < 3.) continue;
      // Dxy significance cut
      //      if(fabs(trackDxyPVTrack[iTrack]/trackDxyErrorTrack[iTrack]) > 10.) continue;
      goodMuons.push_back(i);
    }

    if(int(goodMuons.size()) < 2 ) continue; // no muon pair found
    Npassed[2]+= weight;;
    for(int i=0; i<int(Zcand.size()); i++) FillHistos(Histos[2],i);

    ///////////////////
    // Z selection
    ///////////////////

    // make Z candidates;
    for(int i=0; i<int(goodMuons.size()); i++) {
      TLorentzVector Muon1(pxMuon[goodMuons[i]], pyMuon[goodMuons[i]], pzMuon[goodMuons[i]], energyMuon[goodMuons[i]]);
      for(int j=i+1; j<int(goodMuons.size()); j++) {
	TLorentzVector Muon2(pxMuon[goodMuons[j]], pyMuon[goodMuons[j]], pzMuon[goodMuons[j]], energyMuon[goodMuons[j]]);
	TLorentzVector Z = Muon1+Muon2;
	if(Z.M()> 0. && Z.M()<10000. && chargeMuon[goodMuons[i]]*chargeMuon[goodMuons[j]]<0 ) {
	  Zcand.push_back(Z);
	  //	  if(Muon1.Pt() > Muon2.Pt()) {
	    iZDaugh1.push_back(goodMuons[i]);
	    iZDaugh2.push_back(goodMuons[j]);
	    //	  } else {
	    //	    iZDaugh1.push_back(goodMuons[j]);
	    //	    iZDaugh2.push_back(goodMuons[i]);
	    //	  }
	}
      }
    }

    for( int i = 0; i< Zcand.size(); i++) {
      // write the tree
      mll = Zcand[i].M();
      int iMu1 = iZDaugh1[i];
      int iMu2 = iZDaugh2[i];
      flagMu1 = flagMuons[iMu1];
      flagMu2 = flagMuons[iMu2];
      etaMu1 = etaMuon[iMu1];
      etaMu2 = etaMuon[iMu2];
      phiMu1 = phiMuon[iMu1];
      phiMu2 = phiMuon[iMu2];
      ptMu1 = pTMuon(iMu1);
      ptMu2 = pTMuon(iMu2);
      MET = sqrt(pow(double(pxMet[0]),2.)+pow(double(pyMet[0]),2.));   //<-- you might want to use PFmet here
      phiMET = phiMet[0];
      run = runNumber;
      evNum = eventNumber;
      bx = eventNumber;
      ls = lumiBlock;
      orbit = orbitNumber;
      outTree->Fill();
    }

    // clean the vectors
    Zcand.clear();
    iZDaugh1.clear();
    iZDaugh1.clear();
  }
  
  WriteHistos(Histos[0], file, "noCut");
  for(int i=0; i< int(cuts.size()); i++) 
    WriteHistos(Histos[i+1], file, cuts[i]);
  outTree->Write();
  file->Close();
  
}

int CandleCalib::BestPV(int bestzIdx) {
  // find the highestpT PV
  double maxpT = -9999.;
  for(int i = 0; i < nPV; i++){
    if(SumPtPV[i] > maxpT) {
      iPV = i;
      maxpT = SumPtPV[i];
    }
  }
  return iPV;
}

double CandleCalib::pTMuon(int i) {
  return sqrt(pxMuon[i]*pxMuon[i]+pyMuon[i]*pyMuon[i]);
}

double CandleCalib::SumPt(int iMu, int iZ) {
  double eta0 = etaMuon[iMu];
  double phi0 = phiMuon[iMu];
  double sumPt_tmp = 0;
  for(int i=0; i< nTrack; i++) {
    if(i == trackIndexMuon[iMu]) continue; // take out the muon
    if(trackValidHitsTrack[i] <5) continue;                                     // minimum number of hits  XXX
    if(fabs(transvImpactParTrack[i]/transvImpactParErrorTrack[i]) > 5.) continue;    // track incompatible with the vertex on (x,y)
    if(fabs(PVzPV[BestPV(iZ)]-trackVzTrack[i]) > 0.1) continue;              // track incompatible with the vertex on z
    TVector3 v(pxTrack[i], pyTrack[i], pzTrack[i]);
    if(sqrt(pow(v.Eta()-eta0,2.)+pow(v.Phi()-phi0,2.)) > 0.5) continue; // track outside the cone
    if(v.Pt() < 0.500) continue;     // minimum pT             
    if(v.Pt() > 500.) continue;     // maximum pT             
    sumPt_tmp += v.Pt();
  }
  return sumPt_tmp; 
}

void CandleCalib::EraseZ(int iZ) {
  Zcand.erase(Zcand.begin()+iZ);
  iZDaugh1.erase(iZDaugh1.begin()+iZ);
  iZDaugh2.erase(iZDaugh2.begin()+iZ);
}

double CandleCalib::DeltaPhi_PiHalf(double phi1, double phi2) {
  double dp = fabs(DeltaPhi(phi1, phi2));
  if(dp > asin(1.)) 
    dp = asin(1.)*2. - dp;
  return dp;
}


