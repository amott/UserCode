#ifndef MakeSpinSPlot_h
#define MakeSpinSPlot_h
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
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"

#include "TMatrixT.h"

class MakeSpinSPlot{
public:

  MakeSpinSPlot(RooAbsData *data);
  
  void addSpecies(TString name, RooAbsPdf* pdf, double expYield);

  void calculate();

  RooDataSet* getSWeightDataSet(){ return __sWeightDataSet; }

  RooArgSet* getSWeightVars(){return __sWeightVars;}
protected:
  void computeCovMatrix();

  double computeDenom();

  void computeSWeight();
  RooDataSet* __sWeightDataSet;
  RooArgSet* __sWeightVars;

  int __nSpec;
  std::vector<RooAbsPdf*> __pdfs;
  std::vector<double> __expectedYields;
  std::vector<TString> __speciesNames;

  RooAbsData *__dataSet;
  
  
  TMatrixD *__covMatrix;
};

#endif
