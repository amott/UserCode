#ifndef MakeSpinWorkspace_h
#define MakeSpinWorkspace_h

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
#include "RooCategory.h"

#include <HggOutputReader2.h>

#include <iostream>
#include <map>
#include <vector>


class MakeSpinWorkspace{
public:
  MakeSpinWorkspace(TString outputFileName);
  ~MakeSpinWorkspace();
  static TH2F* getSelectionMap(int map, bool isData);
  static int   passSelection(TH2F* map,float sigEoE, float etaSC, float pt);
  static int   passSelection(float r9);
  static bool getBaselineSelection(HggOutputReader2* h,int maxI,int minI);
  
  const static int nCat=2;

  void addFile(TString fName,TString l,bool is);

  void setRequireCiC(bool b){requireCiC=b;}
  bool getRequireCiC(){return requireCiC;}
  
  void setSelectionMap(int i){selectionMap=i;}
  int getSelectionMap(){return selectionMap;}
  
  void setRunRange(int rl, int rh){runLow=rl; runHigh=rh;}

  RooWorkspace *getWorkspace(){return ws;}

  void MakeWorkspace();

  static void runOnAll( void (*f)(TString,TString), TString mcName);
  static void runOnAllR9( void (*f)(TString,TString), TString mcName);

  void setUseR9(bool b){useR9 = b;}
  void setTightPt(bool b){tightPt = b;}

private:
  std::vector<TString> fileName,label;
  std::vector<bool> isData;
  RooWorkspace *ws;
  bool requireCiC;
  bool tightPt;
  int selectionMap;

  int runLow,runHigh;

  TFile *outputFile;

  bool useR9;

  void AddToWorkspace(TString inputFile,TString tag, bool isData);

};

#endif