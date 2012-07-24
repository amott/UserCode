// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// local includes
#include "Vecbos.hh"
#include "TestAnalysis.hh"

TestAnalysis::TestAnalysis(TTree *tree) : Vecbos(tree) {}

TestAnalysis::~TestAnalysis(){}

vector<TH1D*> TestAnalysis::CreateHistos(string dirname){
  vector<TH1D*> histos;
  return histos;
}

void TestAnalysis::FillHistos(vector<TH1D*> histos){
}

void TestAnalysis::Loop(string outFileName) {
  if(fChain == 0) return;
  
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
   
  }
}
