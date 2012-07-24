// std includes
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>

using namespace std;

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

// VecbosApp includes
#include "include/VecbosBase.hh"
#include "include/Vecbos.hh"
#include "include/VtxStudy.hh"

VtxStudy::VtxStudy(TTree *tree) : Vecbos(tree) {
  ptmin = 0.;
  ptmax = 300.;
  etamin = -5.;
  etamax = 5.;
  dRmin = 0.;
  dRmax = 0.5;
  dzmax = 0.2;
}

VtxStudy::~VtxStudy(){}

vector<TH2D*> VtxStudy::CreateHistos(string dirname, string part){
  vector<TH2D*> histos;
  string name;

  // dR plots
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "pt"+part+"GEN_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, ptmin, ptmax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "eta"+part+"GEN_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, etamin, etamax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "pt"+part+"GEN_RECOmu_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, ptmin, ptmax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "eta"+part+"GEN_RECOmu_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, etamin, etamax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "pt"+part+"RECO_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, ptmin, ptmax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "eta"+part+"RECO_"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, etamin, etamax));
  }

  // z distance
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "dz"+part+"VTX"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, -0.2, 0.2));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "dz"+part+"GEN"+part+""+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, -0.2, 0.2));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "dR"+part+"GEN"+part+""+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 500, dRmin, dRmax, 200, 0., 1.));
  }

  // vtx-mu matching
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "pt"+part+"RECO_PVmatched_dz"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., dzmax, 200, ptmin, ptmax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "eta"+part+"RECO_PVmatched_dz"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., dzmax, 200, etamin, etamax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "dz"+part+"RECO_PVmatched_dz"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., dzmax, 200, 0., dzmax));
  }  
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "pt"+part+"RECO_PVmatched_dzPull"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 10., 200, ptmin, ptmax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "eta"+part+"RECO_PVmatched_dzPull"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 10., 200, etamin, etamax));
  }
  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "dzpull"+part+"RECO_PVmatched_dzPull"+string(numjets);
    histos.push_back(new TH2D(name.c_str(),name.c_str(), 100, 0., 10., 200, 0., 10.));
  }
  
  return histos;
}

vector<TH1D*> VtxStudy::CreateHistos1D(string dirname, string part){
  vector<TH1D*> histos;
  string name;

  for(int ij=1; ij<6; ij++) {
    char numjets[5];
    sprintf(numjets,"%ijets",ij);
    name = "dzVTXGEN"+part+""+string(numjets);
    histos.push_back(new TH1D(name.c_str(),name.c_str(), 200, -1.*dzmax, dzmax));
  }

  return histos;
}

void VtxStudy::FillHistos(vector<TH2D*> histos){

  for(int i =0; i<500; i++) {
    double dRi = (i+1.)/500.*0.5;
    for(int ij=1; ij<6; ij++) {
      if(njets>= ij) {
	// Fill Histograms Gen level
	histos[0+(ij-1)]->Fill(dRi*0.999,pMc[iGen]*sin(thetaMc[iGen]));
	histos[5+(ij-1)]->Fill(dRi*0.999,etaMc[iGen]);
	if(dRM < dRi) {
	  // Fill Histograms Gen level if muon is reco
	  histos[10+(ij-1)]->Fill(dRi*0.999,pMc[iGen]*sin(thetaMc[iGen]));
	  histos[15+(ij-1)]->Fill(dRi*0.999,etaMc[iGen]);	
	  // plots of reco quantities
	  histos[20+(ij-1)]->Fill(dRi*0.999,sqrt(pow(double(pxMuon[iRecomu]),2)+pow(double(pyMuon[iRecomu]),2.)));
	  histos[25+(ij-1)]->Fill(dRi*0.999,etaMuon[iRecomu]);	
	  // plots of reco quantities for matched Vtx]
	  if(nPV >0) {
	    double dZ = PVzPV[0] - muTrackVzMuon[iRecomu];
	    histos[30+(ij-1)]->Fill(dRi*0.999,dZ);
	  }
	  double dZ = zMc[iGen] - muTrackVzMuon[iRecomu];
	  histos[35+(ij-1)]->Fill(dRi*0.999,dZ);
	  histos[40+(ij-1)]->Fill(dRi*0.999,dRM);
	}
      }
    }
  }
  
  double dRi = 0.1;
  for(int ij=1; ij<6; ij++) {
    if(njets>= ij) {
      // plots of reco quantities for matched Vtx]
      if(nPV >0) {
	double dZ = fabs(double(PVzPV[0]) - double(muTrackVzMuon[iRecomu]));
	// dz matched
	for(int i=1; i<= 100; i++) {
	  double dZi = i*1./100.*dzmax;
	  if(dZ <= dZi) {
	    histos[45+(ij-1)]->Fill(dZi*0.999,sqrt(pow(double(pxMuon[iRecomu]),2)+pow(double(pyMuon[iRecomu]),2)));
	    histos[50+(ij-1)]->Fill(dZi*0.999,etaMuon[iRecomu]);
	    histos[55+(ij-1)]->Fill(dZi*0.999,dZ);
	  }
	}
	// dz-pull matched
	double dZpull = fabs(double(PVzPV[0]) - double(muTrackVzMuon[iRecomu]))/sqrt(pow(double(PVErrzPV[0]),2.) + pow(double(muTrackDzErrorMuon[iRecomu]),2.));
	for(int i=1; i<= 100; i++) {
	  double dZpulli = i*1./100.*10.;
	  if(dZpull <= dZpulli) {
	    histos[60+(ij-1)]->Fill(dZpulli*0.999,sqrt(pow(double(pxMuon[iRecomu]),2)+pow(double(pyMuon[iRecomu]),2)));
	    histos[65+(ij-1)]->Fill(dZpulli*0.999,etaMuon[iRecomu]);	
	    histos[70+(ij-1)]->Fill(dZpulli*0.999,dZpull);	
	  }
	}
      }
    }
  }
}

void VtxStudy::FillHistos1D(vector<TH1D*> histos){
  for(int ij=1; ij<6; ij++) {
    if(njets>= ij) {
      double dZ = PVzPV[0] - zMc[iGen];
      histos[0+(ij-1)]->Fill(dZ);
    }
  }
}

////////////////////////////////////////////////
// Electrons

void VtxStudy::FillHistosEle(vector<TH2D*> histos){

  for(int i =0; i<500; i++) {
    double dRi = (i+1.)/500.*0.5;
    for(int ij=1; ij<6; ij++) {
      if(njets>= ij) {
	// Fill Histograms Gen level
	histos[0+(ij-1)]->Fill(dRi*0.999,pMc[iGen]*sin(thetaMc[iGen]));
	histos[5+(ij-1)]->Fill(dRi*0.999,etaMc[iGen]);
	if(dRM < dRi) {
	  // Fill Histograms Gen level if ele is reco
	  histos[10+(ij-1)]->Fill(dRi*0.999,pMc[iGen]*sin(thetaMc[iGen]));
	  histos[15+(ij-1)]->Fill(dRi*0.999,etaMc[iGen]);	
	  // plots of reco quantities
	  histos[20+(ij-1)]->Fill(dRi*0.999,sqrt(pow(double(pxEle[iRecoe]),2)+pow(double(pyEle[iRecoe]),2.)));
	  histos[25+(ij-1)]->Fill(dRi*0.999,etaEle[iRecoe]);	
	  // plots of reco quantities for matched Vtx]
	  if(nPV >0) {
	    double dZ = PVzPV[0] - eleTrackVzEle[iRecoe];
	    histos[30+(ij-1)]->Fill(dRi*0.999,dZ);
	  }
	  double dZ = zMc[iGen] - eleTrackVzEle[iRecoe];
	  histos[35+(ij-1)]->Fill(dRi*0.999,dZ);
	  histos[40+(ij-1)]->Fill(dRi*0.999,dRM);
	}
      }
    }
  }

  double dRi = 0.1;
  for(int ij=1; ij<6; ij++) {
    if(njets>= ij) {
      // plots of reco quantities for matched Vtx
      if(nPV >0) {
        double dZ = fabs(double(PVzPV[0]) - double(eleTrackVzEle[iRecoe]));
        // dz matched
        for(int i=1; i<= 100; i++) {
          double dZi = i*1./100.*dzmax;
          if(dZ <= dZi) {
            histos[45+(ij-1)]->Fill(dZi*0.999,sqrt(pow(double(pxEle[iRecoe]),2)+pow(double(pyEle[iRecoe]),2)));
            histos[50+(ij-1)]->Fill(dZi*0.999,etaEle[iRecoe]);
            histos[55+(ij-1)]->Fill(dZi*0.999,dZ);
          }
        }
        // dz-pull matched                                                                                                                        
        double dZpull = fabs(double(PVzPV[0]) - double(eleTrackVzEle[iRecoe]))/sqrt(pow(double(PVErrzPV[0]),2.) + 
										      pow(double(eleTrackDzErrorEle[iRecoe]),2.));
        for(int i=1; i<= 100; i++) {
          double dZpulli = i*1./100.*10.;
          if(dZpull <= dZpulli) {
            histos[60+(ij-1)]->Fill(dZpulli*0.999,sqrt(pow(double(pxEle[iRecoe]),2)+pow(double(pyEle[iRecoe]),2)));
            histos[65+(ij-1)]->Fill(dZpulli*0.999,etaEle[iRecoe]);
	    histos[70+(ij-1)]->Fill(dZpulli*0.999,dZpull);
          }
        }
      }
    }
  }
}




////////////////////////////////////////////////

void VtxStudy::Loop() {
  if(fChain == 0) return;

  vector< vector<TH1D*> > Histos;

  vector<TH2D*> HistoMu = CreateHistos("histosMu", "Mu");
  vector<TH2D*> HistoEle = CreateHistos("histosEle", "Ele");
  vector<TH1D*> Histo1DMu = CreateHistos1D("histos1DMu", "Mu");
  vector<TH1D*> Histo1DEle = CreateHistos1D("histos1DEle", "Ele");
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;

    // check the alpgen ID
    //    if(AlpgenIdSelection(genAlpgenID, "Wjets") == false) continue;

    // determine number of jets pT > 30 
    njets = 0;
    for(int i=0; i< nSisConeJet; i++ )
      if(sqrt(pow(double(pxSisConeJet[i]),2.)+pow(double(pySisConeJet[i]),2.))> 30.) njets ++;
    
    // MUON

    // Look for the gen-level muon from W
    iGen = -99;
    for(int i=0; i<nMc; i++) 
      if(statusMc[i] == 3) // stable particle
	if(abs(idMc[i]) == 13) // it is a muon
	  if(abs(idMc[mothMc[i]]) == 24) // its mother is a W
	    if(etaMc[i]<2.4 && etaMc[i]>-2.4) // eta in the acceptance
	      if(pMc[i]*sin(thetaMc[i]) > 3.) // enough pT to arrive at the Muon stations
		iGen = i;

    if(iGen != -99) {

      // Look for the corresponding reco mu
      dRM = 99999999.;
      iRecomu = -99;
      for(int i=0; i<nMuon; i++) {
	double dR = sqrt(pow(double(etaMuon[i]-etaMc[iGen]),2.) + pow(double(phiMuon[i]-phiMc[iGen]),2.));
	if(dR < dRM) {
	  dRM = dR;
	  iRecomu = i;
	}
      }
      
      // fill plots as a function of jet multiplicity
      FillHistos(HistoMu);
      FillHistos1D(Histo1DMu);
      
    }

    // PM ELECTRONS
    // Look for the gen-level muon from W
    iGen = -99;
    for(int i=0; i<nMc; i++) 
      if(statusMc[i] == 3) // stable particle
	if(abs(idMc[i]) == 11) // it is a muon
	  if(abs(idMc[mothMc[i]]) == 24) // its mother is a W
	    //	    if(etaMc[i]<2.4 && etaMc[i]>-2.4) // eta in the acceptance
	    //	      if(pMc[i]*sin(thetaMc[i]) > 3.) // enough pT to arrive at the Muon stations
	    iGen = i;
   
    if(iGen != -99)
      FillHistos1D(Histo1DEle);

    if(iGen != -99) {
      // Look for the corresponding reco e
      dRM = 99999999.;
      iRecoe = -99;
      for(int i=0; i<nEle; i++) {
	double dR = sqrt(pow(double(etaEle[i]-etaMc[iGen]),2.) + pow(double(phiEle[i]-phiMc[iGen]),2.));
	if(dR < dRM) {
	  dRM = dR;
	  iRecoe = i;
	}
      }
    // fill plots as a function of jet multiplicity
    FillHistosEle(HistoEle);
    }

  }
 
    TFile *file = new TFile("Histograms_mumu.root","RECREATE");
    WriteHistos(HistoMu, file, "histosMu");
    WriteHistos(Histo1DMu, file, "histos1DMu");
    WriteHistos(HistoEle, file, "histosEle");
    WriteHistos(Histo1DEle, file, "histos1DEle");

}


