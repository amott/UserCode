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
	VecbosEle ele1 = Electrons->at(k);

	// Choose Second Electron
	for (int j=k+1; j<nEle;j++) {   
	  lpass   = 0;
	  tpass   = 0;
	  mvapass = 0;
	  Zeemass = 0;
	  VecbosEle ele2 = Electrons->at(j);
	  // Verify Opposite Charges
	  if (ele1.charge != ele2.charge) {

	    // Calculate Invariant Mass for Z Boson from Two Electrons
	    TLorentzVector Ele1;
	    Ele1.SetPtEtaPhiM(ele1.correctedEnergy/cosh(ele1.eta),ele1.eta,ele1.phi,0);
	    TLorentzVector Ele2;
	    Ele2.SetPtEtaPhiM(ele2.correctedEnergy/cosh(ele2.eta),ele2.eta,ele2.phi,0);
	    TLorentzVector Zee = Ele1 + Ele2;
	    Zeemass = Zee.M();

	    PFIsoOverPT1 = (ele1.dr03ChargedHadronPFIso + max(0.d,(ele1.dr03NeutralHadronPFIso+ele1.dr03PhotonPFIso) - pow(0.3, 2)*PI*max(0.f,rho)))/(ele1.pt);
	    PFIsoOverPT2 = (ele2.dr03ChargedHadronPFIso + max(0.d,(ele2.dr03NeutralHadronPFIso+ele2.dr03PhotonPFIso) - pow(0.3, 2)*PI*max(0.f,rho)))/(ele2.pt);
	    
	    // Loose Cuts - WP 90
	    if (ele1.correctedEnergy>25 && ele2.correctedEnergy>25) {
	      if (((fabs(ele1.SC.eta)<1.44 && fabs(ele1.dEtaSCTrackAtVtx)<0.007 && fabs(ele1.dPhiSCTrackAtVtx)<0.15 && ele1.SC.sigmaIEtaIEta<0.01 && ele1.SC.HoverE<0.12) || (fabs(ele1.SC.eta)>1.52 && fabs(ele1.dEtaSCTrackAtVtx)<0.009 && fabs(ele1.dPhiSCTrackAtVtx)<0.1 && ele1.SC.sigmaIEtaIEta<0.03 && ele1.SC.HoverE<0.10)) && ((fabs(ele2.SC.eta)<1.44 && fabs(ele2.dEtaSCTrackAtVtx)<0.007 && fabs(ele2.dPhiSCTrackAtVtx)<0.15 && ele2.SC.sigmaIEtaIEta<0.01 && ele2.SC.HoverE<0.12) || (fabs(ele2.SC.eta)>1.52 && fabs(ele2.dEtaSCTrackAtVtx)<0.009 && fabs(ele2.dPhiSCTrackAtVtx)<0.1 && ele2.SC.sigmaIEtaIEta<0.03 && ele2.SC.HoverE<0.10))) {
		if (PFIsoOverPT1<0.15 && PFIsoOverPT2<0.15 && ele1.hasMatchedConversion == false && ele2.hasMatchedConversion == false && ele1.expInnerLayersHits!=-999 && ele2.expInnerLayersHits!=-999) {
		  if (ele1.d0Track != 999 && ele2.d0Track != 999 && ele1.dzTrack != 999 && ele2.dzTrack != 999) {
		    lpass = 1;
		    
		    // Tight Cuts - WP 70
		    if (((fabs(ele1.SC.eta)<1.44 && fabs(ele1.dEtaSCTrackAtVtx)<.004 && fabs(ele1.dPhiSCTrackAtVtx)<0.03) || (fabs(ele1.SC.eta)>1.52 && fabs(ele1.dEtaSCTrackAtVtx)<0.005 && fabs(ele1.dPhiSCTrackAtVtx)<0.02)) && ((fabs(ele2.SC.eta)<1.44 && fabs(ele1.dEtaSCTrackAtVtx)<0.004 && fabs(ele2.dPhiSCTrackAtVtx)<0.03) || (fabs(ele2.SC.eta)>1.52 && fabs(ele2.dEtaSCTrackAtVtx)<0.005 && fabs(ele2.dPhiSCTrackAtVtx)<0.02))) {
		      if (PFIsoOverPT1<0.10 && PFIsoOverPT2<0.10) {
			tpass = 1;
		      };
		    };
		  };
		};
	      };
	    
		    
	      // MVA ID Cuts
	      if (((fabs(ele1.SC.eta)<0.8 && ele1.idMVA>0.5) || (fabs(ele1.SC.eta)>0.8 && fabs(ele1.SC.eta)<1.479 && ele1.idMVA>0.120) || (fabs(ele1.SC.eta)>1.479 && ele1.idMVA>0.6)) && ((fabs(ele2.SC.eta)<0.8 && ele2.idMVA>0.5) || (fabs(ele2.SC.eta)>0.8 && fabs(ele2.SC.eta)<1.479 && ele2.idMVA>0.120) || (fabs(ele2.SC.eta)>1.479 && ele2.idMVA>0.6))) {
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
		Ele1eta = ele1.SC.eta;
		Ele2eta = ele2.SC.eta;
		Ele1r9  = ele1.SC.r9;
		Ele2r9  = ele2.SC.r9;
		passloose = lpass;
		passtight = tpass;
		passmva   = mvapass;
		Ele1mva = ele1.idMVA;
		Ele2mva = ele2.idMVA;
		cout << "made it to selection" << ele1.SC.eta << ele2.SC.eta <<ele1.SC.r9 << ele2.SC.r9 <<endl;
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








