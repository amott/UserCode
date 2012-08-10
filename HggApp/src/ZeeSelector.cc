#include <ZeeSelector.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
#define PI atan(1.0)*4
ZeeSelector::ZeeSelector():
  fChain(0),
  valid(false),
  isData_(true),
  Electrons(0)
{
}

ZeeSelector::ZeeSelector(vector<string> fNames, string treeName,string outFName):
  isData_(true),
  Electrons(0)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

ZeeSelector::~ZeeSelector(){
  delete fChain;
}


void ZeeSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int ZeeSelector::init(){
  if(!valid) return -1;
  
  this->setBranchAddresses();
  this->setupOutputTree();
}

void ZeeSelector::Loop(){
  if(!valid) return;

  this->init();
  cout << "Getting Entries ... "  << endl;
  int nEntries = fChain->GetEntries();
  int nbytes = 0;

  for (int eventi=0; eventi<nEntries; eventi++) {
    if(eventi%500==0) {
      cout << ">> Processing Entry " << eventi << "/" << nEntries << endl;
    }    
    nbytes += fChain->GetEvent(eventi);

    // Reset mass reference
    DZmassref = 100;
    isZmass = false;
    // Require the event to have at least two electrons
    if (nEle > 1) {   

      // Choose First Electron
      for (int k=0; k<nEle-1;k++) { 

	// Choose Second Electron
	for (int j=k+1; j<nEle;j++) {   
	  lpass   = 0;
	  tpass   = 0;
	  mvapass = 0;
	  Zeemass = 0;
	  // Verify Opposite Charges
	  if (*Electrons[k].charge != *Electrons[j].charge) {

	    // Calculate Invariant Mass for Z Boson from Two Electrons
	    TLorentzVector Ele1;
	    Ele1.SetPtEtaPhiM(*Electrons[k].correctedEnergy/cosh(*Electrons[k].eta),*Electrons[k].eta,*Electrons[k].phi,0);
	    TLorentzVector Ele2;
	    Ele2.SetPtEtaPhiM(*Electrons[j].correctedEnergy/cosh(*Electrons[j].eta),*Electrons[j].eta,*Electrons[j].phi,0);
	    TLorentzVector Zee = Ele1 + Ele2;
	    Zeemass = Zee.M();

	    PFIsoOverPT1 = (*Electrons[k].dr03ChargedHadronPFIso + max(0.d,(*Electrons[k].dr03NeutralHadronPFIso+*Electrons[k].dr03PhotonPFIso) - pow(0.3, 2)*PI*max(0.f,rho)))/(*Electrons[k].pt);
	    PFIsoOverPT2 = (*Electrons[j].dr03ChargedHadronPFIso + max(0.d,(*Electrons[j].dr03NeutralHadronPFIso+*Electrons[j].dr03PhotonPFIso) - pow(0.3, 2)*PI*max(0.f,rho)))/(*Electrons[j].pt);
	    
	    // Loose Cuts - WP 90
	    if (*Electrons[k].correctedEnergy>25 && *Electrons[j].correctedEnergy>25) {
	      if (((fabs(*Electrons[k].SC.eta)<1.44 && fabs(*Electrons[k].dEtaSCTrackAtVtx)<0.007 && fabs(*Electrons[k].dPhiSCTrackAtVtx)<0.15 && *Electrons[k].SC.sigmaIEtaIEta<0.01 && *Electrons[k].SC.HoverE<0.12) || (fabs(*Electrons[k].SC.eta)>1.52 && fabs(*Electrons[k].dEtaSCTrackAtVtx)<0.009 && fabs(*Electrons[k].dPhiSCTrackAtVtx)<0.1 && *Electrons[k].SC.sigmaIEtaIEta<0.03 && *Electrons[k].SC.HoverE<0.10)) && ((fabs(*Electrons[j].SC.eta)<1.44 && fabs(*Electrons[j].dEtaSCTrackAtVtx)<0.007 && fabs(*Electrons[j].dPhiSCTrackAtVtx)<0.15 && *Electrons[j].SC.sigmaIEtaIEta<0.01 && *Electrons[j].SC.HoverE<0.12) || (fabs(*Electrons[j].SC.eta)>1.52 && fabs(*Electrons[j].dEtaSCTrackAtVtx)<0.009 && fabs(*Electrons[j].dPhiSCTrackAtVtx)<0.1 && *Electrons[j].SC.sigmaIEtaIEta<0.03 && *Electrons[j].SC.HoverE<0.10))) {
		if (PFIsoOverPT1<0.15 && PFIsoOverPT2<0.15 && *Electrons[k].hasMatchedConversion == false && *Electrons[j].hasMatchedConversion == false && *Electrons[k].expInnerLayersHits!=-999 && *Electrons[j].expInnerLayersHits!=-999) {
		  if (*Electrons[k].d0Track != 999 && *Electrons[j].d0Track != 999 && *Electrons[k].dzTrack != 999 && *Electrons[j].dzTrack != 999) {
		    lpass = 1;
		    
		    // Tight Cuts - WP 70
		    if (((fabs(*Electrons[k].SC.eta)<1.44 && fabs(*Electrons[k].dEtaSCTrackAtVtx)<.004 && fabs(*Electrons[k].dPhiSCTrackAtVtx)<0.03) || (fabs(*Electrons[k].SC.eta)>1.52 && fabs(*Electrons[k].dEtaSCTrackAtVtx)<0.005 && fabs(*Electrons[k].dPhiSCTrackAtVtx)<0.02)) && ((fabs(*Electrons[j].SC.eta)<1.44 && fabs(*Electrons[k].dEtaSCTrackAtVtx)<0.004 && fabs(*Electrons[j].dPhiSCTrackAtVtx)<0.03) || (fabs(*Electrons[j].SC.eta)>1.52 && fabs(*Electrons[j].dEtaSCTrackAtVtx)<0.005 && fabs(*Electrons[j].dPhiSCTrackAtVtx)<0.02))) {
		      if (PFIsoOverPT1<0.10 && PFIsoOverPT2<0.10) {
			tpass = 1;
		      };
		    };
		  };
		};
	      };
	    
		    
	      // MVA ID Cuts
	      if (((fabs(*Electrons[k].SC.eta)<0.8 && *Electrons[k].idMVA>0.5) || (fabs(*Electrons[k].SC.eta)>0.8 && fabs(*Electrons[k].SC.eta)<1.479 && *Electrons[k].idMVA>0.120) || (fabs(*Electrons[k].SC.eta)>1.479 && *Electrons[k].idMVA>0.6)) && ((fabs(*Electrons[j].SC.eta)<0.8 && *Electrons[j].idMVA>0.5) || (fabs(*Electrons[j].SC.eta)>0.8 && fabs(*Electrons[j].SC.eta)<1.479 && *Electrons[j].idMVA>0.120) || (fabs(*Electrons[j].SC.eta)>1.479 && *Electrons[j].idMVA>0.6))) {
		mvapass = 1;
	      };
	    };
	  	    
	    // Calculate Difference from True Z Mass 
	    DZmass = fabs(Zeemass - 91.2);
	    // Compare the proximity of uncut Z mass to real Z mass with other electron pairs in event
	    if (DZmass < DZmassref) {
		// Reset the selected Z mass and reference point to this pair
		mass = Zeemass;
		DZmassref = DZmass;          
		Ele1eta = *Electrons[k].SC.eta;
		Ele2eta = *Electrons[j].SC.eta;
		Ele1r9  = *Electrons[k].SC.r9;
		Ele2r9  = *Electrons[j].SC.r9;
		passloose = lpass;
		passtight = tpass;
		passmva   = mvapass;
		Ele1mva = *Electrons[k].idMVA;
		Ele2mva = *Electrons[j].idMVA;
		cout << "made it to selection" << *Electrons[k].SC.eta << *Electrons[j].SC.eta <<*Electrons[k].SC.r9 << *Electrons[j].SC.r9 <<endl;
		isZmass = true;
	    };
	  };
	};
      };
    };
    if (isZmass ==  true) {
      nEleOut = nEle;
      outTree->Fill();
    };
  };

  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  f->Close();
  std::map<std::string,TH1F*>::iterator it;
}



void ZeeSelector::setBranchAddresses(){
  if(!valid) return;
  //Event info
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("evtNumber",&evtNumber);
  //fChain->SetBranchAddress("isRealData",&_isData);
  
  //information for the Electrons
  fChain->SetBranchAddress("nEle",&nEle);
  fChain->SetBranchAddress("Electrons",&Electrons);
 
  // information for the vertex
  fChain->SetBranchAddress("nVtx", &nVtx);
  fChain->SetBranchAddress("rho", &rho);

}

void ZeeSelector::setupOutputTree(){
  outTree = new TTree("ZeeOutput","");
  outTree->Branch("mass",&mass,"mass");
  outTree->Branch("Ele1eta",&Ele1eta,"Ele1eta");
  outTree->Branch("Ele1r9",&Ele1r9,"Ele1r9");
  outTree->Branch("Ele1mva",&Ele1mva,"Ele1mva");
  outTree->Branch("Ele2eta",&Ele2eta,"Ele2eta");
  outTree->Branch("Ele2r9",&Ele2r9,"Ele2r9");
  outTree->Branch("Ele2mva",&Ele2mva,"Ele2mva");
  outTree->Branch("passloose",&passloose,"passloose");  
  outTree->Branch("passtight",&passtight,"passtight");
  outTree->Branch("passmva",&passmva,"passmva");
  outTree->Branch("nEleOut",&nEleOut,"nEleOut");
}








