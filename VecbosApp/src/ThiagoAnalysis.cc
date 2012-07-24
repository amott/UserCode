// std includes
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TLorentzVector.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "ThiagoAnalysis.hh"

/// Z Mass
const double trueZmass = 91.1876;
/// W mass
const double trueWmass = 80.403;
/// Z pdgId
const int Zid = 23;
/// W pdgId
const int Wid = 24;
/// top pdgId
const int topid = 6;
/// Status 3 particles
const int status3 = 3;

ThiagoAnalysis::ThiagoAnalysis(TTree *tree) 
  : Vecbos(tree), globalLeadingJets(10), globalLeadingMuons(5)
{
  globalLeadingJets.clear();
  globalLeadingMuons.clear();
}

ThiagoAnalysis::~ThiagoAnalysis(){}

vector<TH1D*> ThiagoAnalysis::CreateHistos(string dirname){
  
  vector<TH1D*> histos;
  string name;

  name=dirname+"_numjets";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 10, -0.5, 9.5));
  name=dirname+"_jet_et_1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_jet_et_2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_jet_et_3";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_jet_et_4";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_jet_et_5";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_jet_eta_1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_jet_eta_2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_jet_eta_3";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_jet_eta_4";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_jet_eta_5";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_jet_phi_1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_jet_phi_2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_jet_phi_3";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_jet_phi_4";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_jet_phi_5";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  // 16
  name=dirname+"_nummuons";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 10, -0.5, 9.5));
  name=dirname+"_muon_et_1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_muon_et_2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_muon_eta_1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_muon_eta_2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 100, -5., 5.));
  name=dirname+"_muon_phi_1";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_muon_phi_2";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  // 7
  name=dirname+"_dimuon_mass";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 1000, 0., 1000.));
  name=dirname+"_dimuon_pt";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_multijets_mass";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 1000, 0., 1000.));
  name=dirname+"_multijets_pt";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_boson_mass";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 1000, 0., 1000.));
  name=dirname+"_boson_pt";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  // 6
  name=dirname+"_MET";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 5000, 0., 5000.));
  name=dirname+"_deltaphi_muon_jet";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_deltaphi_muon_MET";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  name=dirname+"_deltaphi_jet_MET";
  histos.push_back(new TH1D(name.c_str(),name.c_str(), 72, -TMath::Pi(), TMath::Pi()));
  // 4
  return histos;
  // Total 33 histograms.
}

void ThiagoAnalysis::FillHistos(vector<TH1D*> histos){
  
  int ih = 0;
  int njets = globalLeadingJets.size();
  histos.at(ih)->Fill(njets);
  ++ih;

  // Fifteen histograms for the jets.
  int posEt = 0; 
  int posEta = 1; // 0 = et, 1 = eta, 2 = phi 
  int posPhi = 2;
  int maxNumberJets = 5;
  for(int ijet = 0; (ijet != njets && ijet != maxNumberJets); ++ijet) {
    histos.at(ih + posEt*maxNumberJets + ijet)->Fill(etSisConeJet[globalLeadingJets[ijet]]);
    histos.at(ih + posEta*maxNumberJets + ijet)->Fill(etaSisConeJet[globalLeadingJets[ijet]]);
    histos.at(ih + posPhi*maxNumberJets + ijet)->Fill(phiSisConeJet[globalLeadingJets[ijet]]);
  }
  ih = 16;

  int nmuons = globalLeadingMuons.size();
  histos.at(ih)->Fill(nmuons);
  ++ih;

  // Six histograms for the muons.
  int maxNumberMuons = 2;
  for(int imuon = 0; (imuon != nmuons && imuon != maxNumberMuons); ++imuon) {
    histos.at(ih + posEt*maxNumberMuons + imuon)->Fill(etMuon[globalLeadingMuons[imuon]]);
    histos.at(ih + posEta*maxNumberMuons + imuon)->Fill(etaMuon[globalLeadingMuons[imuon]]);
    histos.at(ih + posPhi*maxNumberMuons + imuon)->Fill(phiMuon[globalLeadingMuons[imuon]]);
  }
  ih = 23;

  // Two histograms for the dimuon.
  TLorentzVector dimuon;
  if(nmuons > 1) {
    dimuon.SetPxPyPzE(pxMuon[globalLeadingMuons[0]] + pxMuon[globalLeadingMuons[1]],
		      pyMuon[globalLeadingMuons[0]] + pyMuon[globalLeadingMuons[1]],
		      pzMuon[globalLeadingMuons[0]] + pzMuon[globalLeadingMuons[1]],
		      energyMuon[globalLeadingMuons[0]] + energyMuon[globalLeadingMuons[1]]);
    histos.at(ih)->Fill(dimuon.M());
    histos.at(ih+1)->Fill(dimuon.Perp());
  }
  ih = 25;

  // Two histograms for the jet system.
  TLorentzVector sumjets = getBosonFromJets();
  if(njets>0)
    histos.at(ih)->Fill(sumjets.M());
  ++ih;
  if(njets>0)
    histos.at(ih)->Fill(sumjets.Perp());
  ++ih;

  int trueZ = getTrueZ();
  int trueW = getTrueW();

  // Guard against events with no vector bosons.
  if (trueZ == 0 && trueW == 0) {
    cerr << "error: no MC W or Z boson found!!! exiting..." << endl;
    exit(1);
  }
  int theBoson = trueZ + trueW; // What a hack... but one of the two is going to be zero, always...

  // Two histograms for the true boson.
  double massBoson = massMc[theBoson];
  double ptBoson = pMc[theBoson]*sin(thetaMc[theBoson]);
  histos.at(ih)->Fill(massBoson);
  ++ih;
  histos.at(ih)->Fill(ptBoson);
  ++ih;

  // One histogram for the missing energy in its dark throne.
  histos.at(ih)->Fill(etMet[0]);
  ++ih;

  // Enough with the LotR references... three deltaphi histos.
  if(njets > 0 && nmuons > 0) {
    double deltaphiMuonJet = DeltaPhi(double(globalLeadingMuons[0]), 
				      double(globalLeadingJets[0])); 
  
    double deltaphiMuonMET = DeltaPhi(double(globalLeadingMuons[0]),
				      double(phiMet[0]));
    
    double deltaphiJetMET  = DeltaPhi(double(globalLeadingJets[0]),
				      double(phiMet[0]));
    histos.at(ih)->Fill(deltaphiMuonJet);
    histos.at(ih+1)->Fill(deltaphiMuonMET);
    histos.at(ih+2)->Fill(deltaphiJetMET);
  }
}

void ThiagoAnalysis::AssignParameters(map<string,double> ndata) {
  jetPtCut = ndata["jetPtCut"];
  muonPtCut = ndata["muonPtCut"];
  muonAcceptance = ndata["muonAcceptance"];
  barrellimit = ndata["barrellimit"];
}

void ThiagoAnalysis::InitParameters() {
  jetPtCut = 30.;
  muonPtCut = 7.;
  muonAcceptance = 2.1;
  barrellimit = 1.3;
}

/// Functions for cuts
int ThiagoAnalysis::NumberOfJets() {
  
  int nPassingJets = 0;
  for(int i = 0; i != nSisConeJet; ++i) {
      double ptJet = momentumSisConeJet[i]*sin(thetaSisConeJet[i]);
      if(ptJet > jetPtCut)
	++nPassingJets;
  }
  return nPassingJets;
}

bool ThiagoAnalysis::BaselineCut(int leadingJet) {
  return(fabs(etaSisConeJet[leadingJet]) < barrellimit);
}

bool ThiagoAnalysis::SelectionCut(int leadingMuon) {
  bool acceptance = fabs(etaMuon[leadingMuon]) < muonAcceptance;
  bool pt = momentumMuon[leadingMuon]*sin(thetaMuon[leadingMuon]) > muonPtCut;
  return(acceptance && pt);
}

/// Functions for getting info from the tree
int ThiagoAnalysis::getTrueZ() {
  int theZ = 0;
  for(int i = 0; i != nMc; ++i) {
    if(idMc[i] == Zid && statusMc[i] == status3) {
      theZ = i;
      break;
    }
  }
  return theZ;
}

int ThiagoAnalysis::getTrueW() {
  int theW = 0;
  for(int i = 0; i != nMc; ++i) {
    if(abs(idMc[i]) == Wid && statusMc[i] == status3) {
      theW = i;
      break;
    }
  }
  return theW;
}

int ThiagoAnalysis::getTrueTop() {
  int theTop = 0;
  for(int i = 0; i != nMc; ++i) {
    if(abs(idMc[i]) == topid && statusMc[i] == status3) {
      theTop = i;
      break;
    }
  }
  return theTop;
}

int ThiagoAnalysis::getBestZ() {
  int theZ = 0;
  double minDifference = 9999.;
  for(int i = 0; i != nZ0ToMuMu; ++i) {
    double thisDifference = fabs(massZ0ToMuMu[i] - trueZmass);
    if(thisDifference < minDifference) {
      theZ = i;
      minDifference = thisDifference;
    }
  }
  return theZ;
}

void ThiagoAnalysis::theSortedJets(vector<int>& sorted){
  sorted.clear();
  vector<pair<double,int> > pT;
  
  for(int i = 0; i != nSisConeJet; ++i) {
    double ptJet = momentumSisConeJet[i]*sin(thetaSisConeJet[i]);
    if(ptJet > jetPtCut)
      pT.push_back(std::make_pair(ptJet,i));
  }
  
  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
   
  for(int i=0; i<pT.size(); i++) 
    sorted.push_back(pT.at(i).second);
}

void ThiagoAnalysis::theSortedMuons(vector<int>& sorted){
  sorted.clear();
  vector<pair<double,int> > pT;
  
  for(int i = 0; i != nMuon; ++i) {
    double ptMuon = momentumMuon[i]*sin(thetaMuon[i]);
    pT.push_back(std::make_pair(ptMuon,i));
  }
  
  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
   
  for(int i=0; i<pT.size(); i++) 
    sorted.push_back(pT.at(i).second);
}

TLorentzVector ThiagoAnalysis::getBosonFromJets() {
  double E = 0.;
  double px = 0.;
  double py = 0.;
  double pz = 0.;
  
  for(int i = 0; i != nSisConeJet; ++i) {
     double ptJet = momentumSisConeJet[i]*sin(thetaSisConeJet[i]);
     if(ptJet > jetPtCut) {
       E += double(energySisConeJet[i]);
       px += double(pxSisConeJet[i]);
       py += double(pySisConeJet[i]);
       pz += double(pzSisConeJet[i]);
     }
  }
  return TLorentzVector(px,py,pz,E);
}

void ThiagoAnalysis::Loop(char* filename) {
  if(fChain == 0) return;

   // Create the four tiers of histogram sets
  vector<TH1D*> original;

  vector<TH1D*> rawdata_1j;
  vector<TH1D*> baseline_1j;
  vector<TH1D*> selection_1j;
  vector<TH1D*> rawdata_2j;
  vector<TH1D*> baseline_2j;
  vector<TH1D*> selection_2j;
  vector<TH1D*> rawdata_3j;
  vector<TH1D*> baseline_3j;
  vector<TH1D*> selection_3j;
  vector<TH1D*> rawdata_4j;
  vector<TH1D*> baseline_4j;
  vector<TH1D*> selection_4j;

  original  = CreateHistos("original");

  rawdata_1j   = CreateHistos("rawdata_1j");
  baseline_1j  = CreateHistos("baseline_1j");
  selection_1j = CreateHistos("selection_1j");
  rawdata_2j   = CreateHistos("rawdata_2j");
  baseline_2j  = CreateHistos("baseline_2j");
  selection_2j = CreateHistos("selection_2j");
  rawdata_3j   = CreateHistos("rawdata_3j");
  baseline_3j  = CreateHistos("baseline_3j");
  selection_3j = CreateHistos("selection_3j");
  rawdata_4j   = CreateHistos("rawdata_4j");
  baseline_4j  = CreateHistos("baseline_4j");
  selection_4j = CreateHistos("selection_4j");
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  //nentries = 10;
  //for (Long64_t jentry=0; jentry<1000;jentry++) {
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    
    // Let us check that there is no top in here.
    if(getTrueTop() != 0) continue;

    // Let us check that there is at least W or Z.
    int trueW = getTrueW();
    int trueZ = getTrueZ();
    if(trueW == 0 && trueZ == 0) continue;

    // choose if you are dealing with W or Z.
    // trueW != 0 --> Want Z
    // trueZ != 0 --> Want W
    if(trueZ != 0) {
      cout << "The world came to a catastrophic end!!!" << endl;
      abort();
    }

    //cout << "trueZ " << trueZ << endl;
    //cout << "trueW " << trueW << endl;

    // clear the information about jets and muons from last entry
    globalLeadingJets.clear();
    globalLeadingMuons.clear();

    // check for muon decay
    bool muonchannel = false;
    int muon_id = 13;
    for(int i = 0; i != 20; ++i)
      if(abs(idMc[i]) == muon_id)
	muonchannel = true;
    if(muonchannel == false) {
      // cout << "nomuon" << endl;
      continue;
    }

    theSortedJets(globalLeadingJets);
    theSortedMuons(globalLeadingMuons);
    // At this point, theSortedJets and theSortedMuons 
    // sizes ARE nJet and nMuon, respectively.
    
    FillHistos(original);
    
    // the RAW data tier.
    int njets = NumberOfJets();
    if(njets > 0)
      FillHistos(rawdata_1j);
    if(njets > 1)
      FillHistos(rawdata_2j);
    if(njets > 2)
      FillHistos(rawdata_3j);
    if(njets > 3)
      FillHistos(rawdata_4j);
    
    
    // the BASE data tier
    bool passBaseline = BaselineCut(globalLeadingJets[0]);
    if(passBaseline) {
      if(njets > 0)
	FillHistos(baseline_1j);
      if(njets > 1)
	FillHistos(baseline_2j);
      if(njets > 2)
	FillHistos(baseline_3j);
      if(njets > 3)
	FillHistos(baseline_4j);
    }

    // the SEL data tier
    bool passSelection = SelectionCut(globalLeadingMuons[0]);
    if(passSelection) {
      if(njets > 0)
	FillHistos(selection_1j);
      if(njets > 1)
	FillHistos(selection_2j);
      if(njets > 2)
	FillHistos(selection_3j);
      if(njets > 3)
	FillHistos(selection_4j);
    }
  }
  
  TFile *file = new TFile(filename,"RECREATE");
  WriteHistos(original,  file, "original");
  WriteHistos(rawdata_1j,   file, "rawdata_1j");
  WriteHistos(baseline_1j,  file, "baseline_1j");
  WriteHistos(selection_1j, file, "selection_1j");
  WriteHistos(rawdata_2j,   file, "rawdata_2j");
  WriteHistos(baseline_2j,  file, "baseline_2j");
  WriteHistos(selection_2j, file, "selection_2j");
  WriteHistos(rawdata_3j,   file, "rawdata_3j");
  WriteHistos(baseline_3j,  file, "baseline_3j");
  WriteHistos(selection_3j, file, "selection_3j");
  WriteHistos(rawdata_4j,   file, "rawdata_4j");
  WriteHistos(baseline_4j,  file, "baseline_4j");
  WriteHistos(selection_4j, file, "selection_4j");
}
