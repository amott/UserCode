#ifndef MakeSpinPlots_h
#define MakeSpinPlots_h
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

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"

#include "MakeSpinWorkspace.h"

#include <map>

class MakeSpinPlots{
public:
  MakeSpinPlots(TString inputFileName,TString outTag);
  ~MakeSpinPlots();
  
  typedef std::pair<TString,TString> tPair;
  typedef std::pair<double,double> dPair;
  typedef std::map<tPair,dPair> paramMap;
  paramMap nSignal, nBackground,fitMean,fitSigEff,fitBkg1Sigma;
  
  void addMCType(TString s);

  RooDataHist* getSimpleBkgSubtraction(TString tag,TString mcName);
  void DrawBlindFit(TString tag,TString mcName);
  void DrawFit(TString tag,TString mcName);
  void DrawSpinBackground(TString tag,TString mcName,bool signal);
  void PlotSignalFits(TString tag,TString mcName);


  void runAll(TString tag,TString mcName);
  void runAll(TString mcName);

  const int nCat;

  void setLumi(float l){lumi = l;}
  void setBasePath(TString s){basePath = s;}

  void setUseR9(bool b){useR9=b;}

private:
  TFile *inputFile;
  RooWorkspace *ws;

  float lumi;
  
  TString basePath;
  TString outputTag;
  std::vector<TString> mcNames;

  bool useR9;
  bool combinedFit;

  void getFitValues(TString tag,TString mcName);

  TStyle *vecbosStyle;
  void setStyle();
};

#endif
