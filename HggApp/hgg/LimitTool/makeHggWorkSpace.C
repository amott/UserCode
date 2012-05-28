// Original Author - Doug Berry
#include "TH1F.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "th1fmorph.C"
#include "branchingRatio.cc"

#include <iostream>
#include <sstream>
#include <string>

//#include "th1fmorph.C"
//#include "Normalization.C"

#include <map>

using namespace std;
using namespace RooFit;



//map<TString, TH1F*> map_th1f; 

double scaleLum ; 

string dtoa(double value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

void dofit(TFile* inFile, double fitmass, vector <TString> InterpolationList, TFile* OutputFile, RooWorkspace* WorkSpace, RooWorkspace* WSoutput,int debug=1) {
  if (fitmass>150 || fitmass<110) {
    cout << "Warning!!!!!!!!!!! You must have an input mass between 110 and 150 GeV!" << endl << "Exiting Program!!!!" << endl;
    return;
  }

  if (floor(fitmass)-fitmass<0.00001 && floor(fitmass)-fitmass>0) fitmass=floor(fitmass);
  if (fitmass-ceil(fitmass)>-0.00001 && fitmass-ceil(fitmass)<0) fitmass=ceil(fitmass);
  
  double Masses[20] = {110,120,123,125,130,135,145}; ///right now145 are wrong. ggF is 150 actually
  unsigned int nMassPoint = 7;
  
  double lowerbound = 0;
  double upperbound = 0;
  for (unsigned int i=0; i<= nMassPoint - 1; i++) {
    if (fitmass>Masses[i] && fitmass<Masses[i+1]) {
      lowerbound = Masses[i];
      upperbound = Masses[i+1];
    } else if (fitmass==Masses[i]) {
      cout<<"fitmass " << fitmass <<" "<<i<<endl; 
      
      if(i< nMassPoint-1){
	lowerbound = Masses[i];
	upperbound = Masses[i+1];
      }else{
	lowerbound = Masses[i-1];
	upperbound = Masses[i];
      }
      
    }
  }
  
  
  cout<<"lower/upper bound " << lowerbound <<" "<< upperbound <<endl;
  
  
  TString MassString = dtoa(fitmass);
  TString LowerBoundString = dtoa(lowerbound);
  LowerBoundString.ReplaceAll(".0","");
  TString UpperBoundString = dtoa(upperbound);
  UpperBoundString.ReplaceAll(".0","");
  RooRealVar RooRealMass = *(WorkSpace->var("mass"));
  
  cout<<" LowerBoundString " << LowerBoundString <<" "<< UpperBoundString <<endl;
  
  
  for (unsigned int k=0; k < InterpolationList.size(); k++) {
    
    TString LowerHistName = InterpolationList[k];
    LowerHistName.ReplaceAll("110",LowerBoundString);
    TString UpperHistName = InterpolationList[k];
    UpperHistName.ReplaceAll("110",UpperBoundString);
    TString HistName = InterpolationList[k];
    HistName.ReplaceAll("110",MassString);
    
    TString HistName0 = HistName; 
    
    TString HistTitle = "Interpolated Mass at ";
    HistTitle += dtoa(fitmass);
    HistTitle += "GeV";

    TH1F* LowerHist =  ((TH1F*)inFile->Get(LowerHistName));
    TH1F* UpperHist =  ((TH1F*)inFile->Get(UpperHistName));
    if( LowerHist == NULL || UpperHist == NULL){
      cout<<"non histogram? " << LowerHistName<<" "<< UpperHistName <<endl; 
      exit(1);
    }
    
    //double Normalization =  GetNorm(lowerbound, LowerHist, upperbound, UpperHist, fitmass);
    TH1F* MCHist = ((TH1F*)inFile->Get(HistName));
    

    //    TH1F* InterpolatedHist = (TH1F*) th1fmorph((Char_t*) HistName.Data(),(Char_t*) HistTitle.Data(),LowerHist,UpperHist,lowerbound,upperbound,fitmass,Normalization,0);
    TH1F* InterpolatedHist = (TH1F*) th1fmorph((Char_t*) HistName.Data(),(Char_t*) HistTitle.Data(),LowerHist,UpperHist,lowerbound,upperbound,fitmass,-1,0);
    
  
    if (MCHist!=NULL) { /// different mass point 
      TString ResidualHistName = HistName;
      ResidualHistName += "_Residual";
      TH1F* ResidualHist = (TH1F*) InterpolatedHist->Clone(ResidualHistName.Data());
      ResidualHist->Add(MCHist,-1);
      OutputFile->WriteTObject(ResidualHist);
      ResidualHistName.ReplaceAll("th1f_","");
      RooDataHist RooDataResidual(Form("roohist_%s",ResidualHistName.Data()),ResidualHistName.Data(),RooRealMass,ResidualHist);
      WSoutput->import(RooDataResidual);
    }
    
    
    OutputFile->WriteTObject(InterpolatedHist,InterpolatedHist->GetName());
    HistName.ReplaceAll("th1f_","");
    RooDataHist RooDataInterpolated(Form("roohist_%s",HistName.Data()),HistName.Data(),RooRealMass,InterpolatedHist);
    WSoutput->import(RooDataInterpolated);

    if( MCHist!=NULL) { /// save original 
      HistName = HistName0; 
      HistName += "_Orig";
      MCHist->SetName(HistName.Data());
      OutputFile->WriteTObject(MCHist);
      HistName.ReplaceAll("th1f_","");
      RooDataHist RooData(Form("roohist_%s",HistName.Data()),HistName.Data(),RooRealMass,MCHist);
      WSoutput->import(RooData);
    }
  }
  
}

// /afs/cern.ch/user/c/cmshgg/public/Cert_160404-180252_7TeV_All2011_Nov30ReReco_v1.lumi

//  ==  =  Total :
// | Delivered LS | Delivered(/fb) | Selected LS | Recorded(/fb) |
// ---------------------------------------------------------------
// |       189407 |          5.417 |      162453 |         4.781 |



int phtcorr = 96; 
//int phtcorr = 219;


 void makeHggWorkSpace(double fitmass,double lum,TString InputFileName, TString OutputFileName) {
  
  
  double lumOLD =   5000; /// pb used in testSelection.C
    
  scaleLum = 1; //lum / lumOLD; 
  

  /*   
  TString FileName = "workSpace/CMS-HGG_interpolated_";
  FileName += dtoa(fitmass);
  FileName += TString("_lum") + dtoa(lum);
  FileName.ReplaceAll(".","_");
  FileName += ".root";
  */

  if (fitmass>150 || fitmass<110) {
    cout << "Warning!!!!!!!!!!! You must have an input mass between 110 and 150 GeV!" << endl << "Exiting Program!!!!" << endl;
    exit(1);
  }
  
  TFile* InputFile = new TFile(InputFileName,"read");
  
  TList* HistList = InputFile->GetListOfKeys();
    
  RooWorkspace * WorkSpace = (RooWorkspace*) InputFile->Get("cms_hgg_workspace");
  RooWorkspace * WorkSpaceNew = new RooWorkspace("cms_hgg_workspace");
  
  TFile* OutputFile = new TFile( OutputFileName,"RECREATE");
  OutputFile->cd();
  vector<TString> InterpolationList;
  for(Int_t j=0;j<HistList->GetSize();j++){
    TString HistName = HistList->At(j)->GetName();
    if(HistName.Contains("th1f") && HistName.Contains("110") ){ // just put the histograms for 110, we then change this string to the correct one
      InterpolationList.push_back(HistName);
    }
  }
  
  dofit(InputFile,fitmass, InterpolationList, OutputFile, WorkSpace,WorkSpaceNew);
  WorkSpaceNew->Write();
  
  OutputFile->Close();
  delete WorkSpace;
  delete InputFile;
  delete OutputFile;
  cout << "Done!" << endl;
}

//makeHggWorkSpace(double fitmass,double lum,TString InputFileName, TString OutputFileName)
void MakeAll(double startMass, double stopMass, double stepSize, double lumi, TString InputFileName, TString OutputFolder){

  const double minMass = 110;
  const double maxMass = 150;

  for(double mass = startMass; mass < stopMass; mass+=stepSize){
    if(mass < minMass || mass > maxMass) continue;
    TString OutputFileName = Form("%s/CMS-HGG_interpolated_%0.1f_lum%04.0f.root",OutputFolder.Data(),float(mass),lumi*1000);
    makeHggWorkSpace(mass,lumi,InputFileName, OutputFileName);
  }
}
