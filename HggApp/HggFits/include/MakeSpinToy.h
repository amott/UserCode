#ifndef MakeSpinToy_h
#define MakeSpinToy_h
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


#include <vector>
#include <TRandom3.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <math.h>

#include "MakeSpinFits.h"
#include "MakeSpinPlots.h"
#include "MakeSpinWorkspace.h"

class MakeSpinToy{
public:
  MakeSpinToy(TString fileName,TString wsName="cms_hgg_spin_workspace");
  void runMCStudy(int Ntoys,float lumi,TString cat);

  void setTargetLumi(float l){targetLumi = l;}
  void setNominalLumi(float l){nominalLumi = l;}

  enum genType{Data,Hgg125,ALT125};
  
  float getExpEvents(float lumi, TString cat,TString mcType,RooWorkspace* toyws=0);

  double computeLL(RooAbsPdf* pdf, RooAbsData* data,RooRealVar* var, int rebin=1);

  double* run1(genType gen, int& N);

  void runN(int N);

  void save(TString outputFile);

  TH1F* getHistogram(RooAbsData* data, TString title,int rebin=1);

  RooWorkspace *ws;
  RooRealVar* cosT, *mass;
  RooRealVar *S,*GenMinusFit;
  RooDataSet *S_TruthHgg, *S_TruthALT, *S_TruthData;
  RooDataSet *S_tot_TruthHgg, *S_tot_TruthALT, *S_tot_TruthData;
  RooDataSet *S_splot_TruthHgg, *S_splot_TruthALT, *S_splot_TruthData;
  RooDataSet *S_2D_TruthHgg, *S_2D_TruthALT, *S_2D_TruthData;
  RooDataSet *S_2DFIT_TruthHgg, *S_2DFIT_TruthALT, *S_2DFIT_TruthData;

  RooHistPdf *altPdf,*hggPdf,*bkgPdf;

  TTree* makeForCombineTool(TString treeName, RooAbsData* hggData, RooAbsData* altData,RooAbsData* dataData=0);

  void setSaveWorkspaces(bool b){saveWorkspaces = b;}

  void setDoData(bool b);

  void setMCComparison(TString s){mcLabels[2]=s;}

  const static int nBins = 9;

  const int nCat;
  TObjArray toyWSs;
protected:
  float targetLumi;
  float nominalLumi;
  void generateToyWorkspace(RooWorkspace& toyws, const char* cat,genType gen,float nSigTot);
  void generateToyWorkspace(RooWorkspace& toyws,genType gen);
  
  bool useR9;
  bool saveWorkspaces;
  bool doData;
  std::vector<TString> catLabels;

  TString mcLabels[3];
};

#endif
