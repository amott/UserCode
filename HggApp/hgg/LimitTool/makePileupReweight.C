#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TString.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
//using namespace std;

void makePileupReweight(TString MCFile, TH1F *dataHist,TString outputName, TString chainName){
  TChain* chain = new TChain(chainName);
  chain->AddFile(MCFile);
  
  TH1F* mcPU = new TH1F("mcPU","",dataHist->GetNbinsX(),0,dataHist->GetNbinsX());
  chain->Project("mcPU","nPU");
  mcPU->Scale(1./mcPU->Integral());
  
  TH1F* pu = (TH1F*)dataHist->Clone("pileupReWeight");
  pu->Divide(mcPU);
  
  TFile *f = new TFile(outputName,"RECREATE");
  pu->Write();
  f->Close();
  delete chain;
  delete mcPU;
  delete pu;
}

void makeAllPU(string mcFileList,TString dataFileName){
  TFile * dataFile = new TFile(dataFileName);
  TH1F* dataHist = (TH1F*)dataFile->Get("pileup");
  dataHist->Scale(1./dataHist->Integral());

  string mcFile;
  ifstream *s = new ifstream(mcFileList.c_str());
  while(s->good()){
    getline(*s,mcFile);
    string outputFile = mcFile.substr(mcFile.find_last_of('/')+1);
    cout << mcFile << " >> " << outputFile << endl;
    outputFile.insert(outputFile.find(".root"),"_puReWeight");
    cout << mcFile << " >> " << outputFile << endl;
    makePileupReweight(mcFile,dataHist,outputFile,"HggReduce");
  }
  dataFile->Close();
  s->close();

}

