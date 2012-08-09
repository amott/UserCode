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
  isData_(true)
{
}

ZeeSelector::ZeeSelector(vector<string> fNames, string treeName,string outFName):
  isData_(true)
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
      nbytes += fChain->GetEvent(eventi);
    }
    // Reset mass reference
    DZmassref = 100;
    
    // Require the event to have at least two electrons
    if (nEle > 1) {   

      // Choose First Electron
      for (int k=0; k<nEle-1;k++) { 

	// Choose Second Electron
	for (int j=k+1; j<nEle;j++) {   

	  // Verify Opposite Charges
	  if (charge[k] != charge[j]) {

	    // Calculate Invariant Mass for Z Boson from Two Electrons
	    TLorentzVector Ele1;
	    Ele1.SetPtEtaPhiM(energySC[k]/cosh(eta[k]),eta[k],phi[k],0);
	    TLorentzVector Ele2;
	    Ele2.SetPtEtaPhiM(energySC[j]/cosh(eta[j]),eta[j],phi[j],0);
	    TLorentzVector Zee = Ele1 + Ele2;
	    Zeemass = Zee.M();

	    PFIsoOverPT1 = (chargedPFiso[k] + max(0.d,(neutralPFiso[k]+photonPFiso[k]) - pow(0.3, 2)*PI*max(0.f,rho)))/(PT[k]);
	    PFIsoOverPT2 = (chargedPFiso[j] + max(0.d,(neutralPFiso[j]+photonPFiso[j]) - pow(0.3, 2)*PI*max(0.f,rho)))/(PT[j]);
	    
	    // Loose Cuts - WP 90
	    if (energySC[k]>25 && energySC[j]>25) {
	      if (((fabs(etaSC[k])<1.44 && fabs(dEta[k])<0.007 && fabs(dPhi[k])<0.15 && sigmaIEtaIEta[k]<0.01 && HoverE[k]<0.12) || (fabs(etaSC[k])>1.52 && fabs(dEta[k])<0.009 && fabs(dPhi[k])<0.1 && sigmaIEtaIEta[k]<0.03 && HoverE[k]<0.10)) && ((fabs(etaSC[j])<1.44 && fabs(dEta[j])<0.007 && fabs(dPhi[j])<0.15 && sigmaIEtaIEta[j]<0.01 && HoverE[j]<0.12) || (fabs(etaSC[j])>1.52 && fabs(dEta[j])<0.009 && fabs(dPhi[j])<0.1 && sigmaIEtaIEta[j]<0.03 && HoverE[j]<0.10))) {
		if (PFIsoOverPT1<0.15 && PFIsoOverPT2<0.15 && hasMatchedConversion[k] == false && hasMatchedConversion[j] == false && expInnerLayersHits[k]!=-999 && expInnerLayersHits[j]!=-999) {
		  if (d0[k] != 999 && d0[j] != 999 && dz[k] != 999 && dz[j] != 999) {
		    lpass = 1;
		    
		    // Tight Cuts - WP 70
		    if (((fabs(etaSC[k])<1.44 && fabs(dEta[k])<.004 && fabs(dPhi[k])<0.03) || (fabs(etaSC[k])>1.52 && fabs(dEta[k])<0.005 && fabs(dPhi[k])<0.02)) && ((fabs(etaSC[j])<1.44 && fabs(dEta[k])<0.004 && fabs(dPhi[j])<0.03) || (fabs(etaSC[j])>1.52 && fabs(dEta[j])<0.005 && fabs(dPhi[j])<0.02))) {
		      if (PFIsoOverPT1<0.10 && PFIsoOverPT2<0.10) {
			tpass = 1;
		      };
		    };
		  };
		};
	      };
	    
		    
	      // MVA ID Cuts
	      if (((fabs(etaSC[k])<0.8 && idMVA[k]>0.5) || (fabs(etaSC[k])>0.8 && fabs(etaSC[k])<1.479 && idMVA[k]>0.120) || (fabs(etaSC[k])>1.479 && idMVA[k]>0.6)) && ((fabs(etaSC[j])<0.8 && idMVA[j]>0.5) || (fabs(etaSC[j])>0.8 && fabs(etaSC[j])<1.479 && idMVA[j]>0.120) || (fabs(etaSC[j])>1.479 && idMVA[j]>0.6))) {
		mvapass = 1;
	      };
	    };
	  	    
	    // Calculate Difference from True Z Mass 
	    DZmass = fabs(Zeemass - 91.2);
	    // Compare the proximity of uncut Z mass to real Z mass with other electron pairs in event
	    if (DZmass < DZmassref && (mvapass + lpass + tpass > 0)) {
		// Reset the selected Z mass and reference point to this pair
		mass = Zeemass;
		DZmassref = DZmass;          
		Ele1eta = etaSC[k];
		Ele2eta = etaSC[j];
		Ele1r9  = r9[k];
		Ele2r9  = r9[j];
		passloose = lpass;
		passtight = tpass;
		passmva   = mvapass;
		Ele1mva = idMVA[k];
		Ele2mva = idMVA[j];
	    };
	  };
	  lpass   = 0;
	  tpass   = 0;
	  mvapass = 0;
	  Zeemass = 0;
	};
      };
    };
    nEleOut = nEle;
    outTree->Fill();
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
  fChain->SetBranchAddress("Electrons.charge",charge);
  fChain->SetBranchAddress("Electrons.correctEnergy",energySC);
  fChain->SetBranchAddress("Electrons.SC.eta",etaSC);
  fChain->SetBranchAddress("Electrons.SC.phi",phiSC);
  fChain->SetBranchAddress("Electrons.eta",eta);
  fChain->SetBranchAddress("Electrons.phi",phi);
  fChain->SetBranchAddress("Electrons.SC.HoverE",HoverE);
  fChain->SetBranchAddress("Electrons.SC.dEtaSCTrackAtVtx",dEta);
  fChain->SetBranchAddress("Electrons.dPhiSCTrackAtVtx",dPhi);
  fChain->SetBranchAddress("Electrons.SC.sigmaIEtaIEta",sigmaIEtaIEta);
  fChain->SetBranchAddress("Electrons.d0Track",d0);
  fChain->SetBranchAddress("Electrons.dzTrack", dz);
  fChain->SetBranchAddress("Electrons.SC.idMVA", idMVA);
  fChain->SetBranchAddress("Electrons.SC.r9", r9);
  fChain->SetBranchAddress("Electrons.pt", PT);
  fChain->SetBranchAddress("Electrons.SC.dr03ChargedHadronPFIso", chargedPFiso);
  fChain->SetBranchAddress("Electrons.SC.dr03NeutralHadronPFIso", neutralPFiso);
  fChain->SetBranchAddress("Electrons.SC.dr03PhotonPFIso", photonPFiso);
  fChain->SetBranchAddress("Electrons.hasMatchedConversion", hasMatchedConversion);
  fChain->SetBranchAddress("Electrons.expInnerLayersHits", expInnerLayersHits);
  
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
  outTree->Branch("passtight",&passtight,"passtight");
  outTree->Branch("passmva",&passmva,"passmva");
  outTree->Branch("nEleOut",&nEleOut,"nEleOut");
}








