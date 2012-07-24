#define GenZjets_cxx
#include "GenZjets.hh"
#include "GenVecbos.h"
#include "Vecbos.hh"
#include "Jet.hh"

// Random auxiliary function thrown here.
int GenZjets::numJetsAbovePtAndEta(std::vector<Jet> & theJets, double thePtCut, double theEtaCut) {
  int passingNumJets = 0;
  for (size_t i = 0; i < theJets.size(); ++i) {
    if(theJets.at(i).pt() > thePtCut &&
       fabs(theJets.at(i).eta()) < theEtaCut) passingNumJets++;
  }
  return passingNumJets;
}

GenZjets::GenZjets(TTree *tree, double xsec, double lumi) : GenVecbos(tree) {
  //  GenInit(tree);
  Long64_t nentries = tree->GetEntries();
  _weight = (nentries > 0 ? xsec*lumi/nentries : 0.);
  _EtaMax = 3.;
  _PtZMin =  0.;
}

GenZjets::~GenZjets() {}

void GenZjets::SetEtaMax(double max) { _EtaMax = max;}

void GenZjets::SetPtZMin(double max) { _PtZMin = max;}

void GenZjets::Loop(string outFileName) {

  double ptmin = 5.0;
  double ptmax = 100.;
  int npt = 10;

  double dRmin = 0.3;
  double dRmax = 1.3;
  int ndR = 10;
  
  TH2D* Nevj[6];
  Nevj[0] = new TH2D("Nev0j", "Nev0j", npt, ptmin, ptmax+ptmin, ndR, dRmin, dRmax);
  Nevj[1] = new TH2D("Nev1j", "Nev1j", npt, ptmin, ptmax+ptmin, ndR, dRmin, dRmax);
  Nevj[2] = new TH2D("Nev2j", "Nev2j", npt, ptmin, ptmax+ptmin, ndR, dRmin, dRmax);
  Nevj[3] = new TH2D("Nev3j", "Nev3j", npt, ptmin, ptmax+ptmin, ndR, dRmin, dRmax);
  Nevj[4] = new TH2D("Nev4j", "Nev4j", npt, ptmin, ptmax+ptmin, ndR, dRmin, dRmax);
  Nevj[5] = new TH2D("Nev5j", "Nev5j", npt, ptmin, ptmax+ptmin, ndR, dRmin, dRmax);
  for(int i=0; i!=6; ++i) {
    Nevj[i]->Sumw2();
  }
   if (fChain == 0) return;

   Long64_t nbytes = 0, nb = 0;
   Long64_t nentries = fChain->GetEntries();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue; 
      
      if(jentry%1000==0) std::cout << jentry << " entries done." << std::endl;
      
      // remove Z->tautau decays (we don't like taus)
      bool isThisTau = false;
      for(int i = 0; i< npar; i++) {
	if(status[i] != 3) continue; // stable particle 
	if(fabs(id[i]) == 15) isThisTau =  true;
      }
      if(isThisTau) continue;

      //  Look for the two mu from Z
      bool foundMu1 = false;
      bool foundMu2 = false;
      TLorentzVector Mu1;
      TLorentzVector Mu2;
      vector<TLorentzVector> JetConstituents;
      for(int i = 0; i< npar; i++) {
	if(status[i] != 1) continue; // stable particle 
	if(!foundMu1 && 
	   (fabs(id[i]) == 11 || fabs(id[i]) == 13) ) {
	  foundMu1 = true;
	  Mu1 = TLorentzVector(px[i], py[i], pz[i], E[i]);
	} 
	else if(!foundMu2 && 
		(fabs(id[i]) == 11 || fabs(id[i]) == 13) ) {
	  foundMu2 = true;
	  Mu2 = TLorentzVector(px[i], py[i], pz[i], E[i]);
	} 
	else {
	  //std::cout << "Found stable!" << std::endl;  
	  JetConstituents.push_back(TLorentzVector(px[i], py[i], pz[i], E[i]));
	}
      }

      // something was found
      if(foundMu1 == false || foundMu2 == false) continue;
      if(JetConstituents.size() == 0) continue; 

      // apply Z pT selection
      TLorentzVector Z = Mu1+Mu2;
      if(Z.Pt()<_PtZMin) continue;
      
      // cluster with different dR
      for(int i=0; i<ndR; i++) {
	double dR = dRmin + double(i+0.5)/double(ndR)*(dRmax-dRmin);
	vector<Jet> SisConeGenJets = SISCone(JetConstituents, dR, ptmin);
	// Cut at different pTs
	for(int j=0; j<npt; j++) {
	  double pT = ptmin + double(j+0.5)/double(npt)*(ptmax);
	  int njets = numJetsAbovePtAndEta(SisConeGenJets, pT, _EtaMax);
	  for(int k=0; k<6; k++) {
	    if(njets>=k) Nevj[k]->Fill(pT,dR,_weight);
	  }
	}
      }
   }
   
   TFile* file = new TFile(outFileName.c_str(),"recreate");
   for(int k=0; k<6; k++) Nevj[k]->Write();
   file->Close();
   
}

