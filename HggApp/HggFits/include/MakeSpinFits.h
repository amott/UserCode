#ifndef MakeSpinFits_h
#define MakeSpinFits_h

#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooGlobalFunc.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include "RooStats/SPlot.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooBernstein.h"
#include "RooLinkedListIter.h"
#include "RooCBShape.h"

#include <HggOutputReader2.h>

#include <iostream>
#include <map>
#include <vector>

#include "MakeSpinWorkspace.h"

class MakeSpinFits{
public:
  MakeSpinFits(TString inputFileName, TString outputFileName);
  ~MakeSpinFits();

  RooWorkspace *getWorkspace(){return ws;}

  void setWorkspace(RooWorkspace *inputWs){ws = inputWs;}

  void addMCLabel(TString l){mcLabel.push_back(l);}

  void MakeSignalFit(TString tag,TString mcName);

  void MakeBackgroundFit(TString mcName,TString catTag,float initMass,float range,bool gausPen,TString inputTag="");

  void MakeBackgroundFitCosTBin(TString mcName,TString catTag,float minCosT,float maxCosT);

  void MakeBackgroundOnlyFit(TString mcName,TString catTag);

  void MakeAllBackgroundFits(TString cat, TString mcTag);

  void setAddSWeight(bool b){addSWeight=b;}

  void run();

  void save();

  bool setUseR9(bool b){useR9 = b;}

  const int nCat;
private:
  RooWorkspace *ws;

  std::vector<TString> mcLabel;
  
  bool addSWeight;
			 
  TFile *inputFile;
  TFile *outputFile;

  bool useR9;

  float mean0,meanE0;

  void AddSWeight(TString mcName, TString catTag,TString inputTag);

};

#endif