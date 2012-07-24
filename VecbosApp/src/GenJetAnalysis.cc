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

#include "include/Vecbos.hh"
#include "include/GenJetAnalysis.hh"

GenJetAnalysis::GenJetAnalysis(TTree *tree) : Vecbos(tree) {}

GenJetAnalysis::~GenJetAnalysis(){}

vector<TH1D*> GenJetAnalysis::CreateHistos(string dirname){
  vector<TH1D*> histos;
  string name;

 
  return histos;
}

void GenJetAnalysis::FillHistos(vector<TH1D*> histos){

  int ih = -1;
  
 
  
}

void GenJetAnalysis::Loop() {
  if(fChain == 0) return;

  vector< vector<TH1D*> > Histos;

  //  for(int i=0; i<41; i++) {
  char name[32];
  //  sprintf(name,"0_0_%i",-20+i);
  sprintf(name,"0_0_%i",0);
  Histos.push_back(CreateHistos(name));
  //}
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = 10000;//fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    // check the alpgen ID
    
    vector<TLorentzVector> GenCand;
    for(int i = 0; i < nMc; i++){
      if(statusMc[i] == 1){
	TLorentzVector J;
	J.SetPtEtaPhiE(pMc[i]/cosh(etaMc[i]), etaMc[i], phiMc[i], energyMc[i]);
	GenCand.push_back(J);
      }
    }

    vector<Jet> GenJets = SortJet(SISCone(GenCand, 0.7, 0.0));

    for(int i = 0; i < GenJets.size(); i++){
      cout << GenJets[i].et() << " " << GenJets[i].eta() << endl;
    }
    cout << endl << endl;
  }
  
  TFile *file = new TFile("Histograms.root","RECREATE");
  //for(int i=0; i<41; i++) {
  //  char name[32];
  sprintf(name,"0_0_%i",0);//-20+i);
  //  WriteHistos(Histos[i], file,name);
  WriteHistos(Histos[0], file,name);
  //} 
}
