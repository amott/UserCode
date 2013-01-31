#ifndef MakeSpinFits_h
#define MakeSpinFits_h

//! class to make all the Hgg Fits

/*!

This class takes in a properly formatted workspace output from MakeSpinWorkspace and runs all the fits necessary
to extract the signal yield and cos(theta) distributions from the data.

Author: Alex Mott (Caltech)
Date: Jan 2013
*/

#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>

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
#include "RooProdPdf.h"

#include <HggOutputReader2.h>

#include <iostream>
#include <map>
#include <vector>
#include <assert.h>


#include "MakeSpinWorkspace.h"
#include "MakeSpinSPlot.h"

class MakeSpinFits{
public:
  MakeSpinFits(); //!< default constructor
  /*!
    If the default constructor is used, the input workspace must be input manually
    with the setWorkspace() method.
    Saving the workspace must then be done manually (the save() method will perform no action)
    This is used most often for fitting toy data
  */
  MakeSpinFits(TString inputFileName, TString outputFileName); //!< constructor requires the paths to the input and output files for the workspaces
  /*!
    \param inputFileName should point to a file containing a workspace named "cms_hgg_spin_workspace" that is 
    the output from MakeSpinWorkspace.
    \param outputFileName points to the location where a new file will be created containing all the fits as well as the original data for making plots and further testing
   */
  ~MakeSpinFits();

  RooWorkspace *getWorkspace(){return ws;} //<! returns the working RooWorkspace

  void setWorkspace(RooWorkspace *inputWs){ //<! set the input workspace 
    ws = inputWs;
    getLabels("labels",&mcLabel,ws);
    getLabels("evtcat",&catLabels,ws);
  }

  void MakeSignalFit(TString tag,TString mcName); //<! used to make the signal fits
  void MakeSignalFitForFit(TString tag, TString mcName); //<! Copies the signal fits to knew pdfs whose parameters can be floated in a fit to data without distrubing the original fits

  void MakeCombinedSignalSpin(TString mcName); //<! Make RooHistPdfs of the cos(theta) distribution for inclusive signal samples

  void MakeBackgroundOnlyFit(TString catTag); //<! Make a background only fit to data in single category the type of fit is controlled by the fitType member

  void MakeCombinedSignalTest(TString mcName); //<! Make a combined simultaneous fit to S+B
  /*!
    Performs a simultaneous signal hypothesis test and extracts a total number of signal events.
    The fraction of signal in the different categories is fixed in the fit (by the signal MC sample) and only the total yield is allowed to float
  */
  void Make2DCombinedSignalTest(TString massMcName,TString costMcName); //<! Make a 2D S+B fit in the mass X cos(theta) plane
  /*!
    Prforms a simultaneous signal hypothesis test in using 2D pdfs of mass X cos(theta).  The cos(theta) pdfs are RooHistPdfs taken from the signal MC and data sidebands.
    As with MakeCombinedSignalTest() the fraction of events in each category is fixed.

    The signal model for mass and cos(theta) can be specified indepedantly to allow for modelling (for instance) a spin-2 signal with a spin-0-like lineshape
    \param massMcName the name of the signal model to use for the Mass dimension of the fitting model.
    \param costMcName the name of the signal model to use for the cos(theta) dimension of the fitting model.
  */

  void MakeFloatingSignalTest(TString mcName); //<! Make a S+B signal fit with the category yields floated independently
  /*!
    Prforms a S+B fit with the individual category signal yields unconstrained.  Useful for channel compatibility measurement
    \sa MakeCombinedSignalTest()
  */
  void Make2DFloatingSignalTest(TString massMcName,TString costMcName); //<! Make a 2D S+B test with category yields floated independently
  /*!
    Peorms a 2D S+B fit with the individual category signal yields unconstrained.  Useful for channel compatibility measurement
    \sa Make2DCombinedSignalTest()
  */

  void getSimpleBkgSubtraction(TString mcName,TString tag); //<! makes the background-subtracted cos(theta) distribution from data in the specified category
  /*! 
    Builds a background-sbutracted RooDataHist for data in the category specified by tag and using the yields and line-shape specified by mcName
  */
  void getSimpleTotalBkgSubtraction(TString mcName); //<! builds the background-subtracted cos(theta) distribution for the inclusive dataset
  void setAddSWeight(bool b){addSWeight=b;} //<! specify whether to add the SWeighted dataset to the output workspace

  void run(); //<! run all fits in the correct order

  void save(); //<! save the output workspace

  static float computeFWHM(RooAbsPdf* pdf, float mean, RooRealVar* var); //<! compute the Full Width at Half Maximum for a pdf

  enum BkgFitType{kExp,kPoly}; //<! allowed types for background fit

  void setBkgFit(BkgFitType t){fitType=t;} //<! specify which type of background fit to use

  void AddCombinedBkgOnlySWeight(TString mcName); //<! add the SWeighted datasets from the combined fit

  static void getLabels(const char *varName, std::vector<TString> *lblVec,RooWorkspace *w);
protected:
  RooWorkspace *ws;

  std::vector<TString> mcLabel;

  std::vector<TString> catLabels;

  bool addSWeight;
			 
  TFile *inputFile;
  TFile *outputFile;

  float mean0,meanE0;

  BkgFitType fitType;
};

#endif
