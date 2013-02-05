#ifndef MakeSpinPlots_h
#define MakeSpinPlots_h
/*!
  Used to make the output plots of data, MC and fits for the Hgg analysis.
  Requires a workspace file output from MakeSpinFits

  Author: Alex Mott (Caltech)
  Date: Jan 2013
*/

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
#include "RooHist.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TBox.h"
#include "TLine.h"
#include "TF1.h"
#include "TText.h"

#include "MakeSpinFits.h"

#include <map>

class MakeSpinPlots{
public:
  MakeSpinPlots(TString inputFileName,TString outTag); //!< Constructor requires path to input workspace file and a tag which will be inserted into every output filename
  ~MakeSpinPlots();

  typedef std::pair<TString,TString> tPair;
  typedef std::pair<double,double> dPair;
  typedef std::map<tPair,dPair> paramMap;
  
  void DrawBlindFit(TString tag,TString mcName); //!< draw blinded mass distribution of the given category
  void DrawFit(TString tag,TString mcName); //!< draw non-blinded mass distributions
  void DrawIndFit(TString tag,TString mcName); //!< draw non-blinded mass distributions using the floated category fits
  void DrawSpinBackground(TString tag,TString mcName,bool signal); //!< Draw cos(theta) SPlots
  void DrawSpinSubBackground(TString tag,TString mcName,bool signal); //!< Draw cos(theta) background-subtracted plots
  void DrawSpinSubTotBackground(TString mcName,bool signal); //!< Draw inclusive cos(theta) background-subtracted plots
  void PlotSignalFits(TString tag,TString mcName); //!< Draw the signal fits


  void runAll(TString tag,TString mcName); //!< make all plots for given category and MC type (for the signal fits)
  void runAll(TString mcName); //!< Draw all plots for all categories

  void setLumi(float l){lumi = l;} //!< set the lumi for the plots to display
  void setBasePath(TString s){basePath = s;} //!< set the base path to which to save the figures

  void printYields(const char* mcType);
  void MakeChannelComp(const char* mcType);
protected:
  TFile *inputFile;
  RooWorkspace *ws;

  float lumi;

  paramMap nSignal, nBackground,fitMean,fitSigEff,fitBkg1Sigma;
  
  TString basePath;
  TString outputTag;
  std::vector<TString> mcNames;
  std::vector<TString> catNames;

  void getFitValues(TString tag,TString mcName); //!< fill the maps with fitted yields from the signal

  TStyle *vecbosStyle;
  void setStyle(); //!< setup the plot style 
};

#endif
