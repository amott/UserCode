#ifndef MakeSpinWorkspace_h
#define MakeSpinWorkspace_h
//! class to make the RooWorkspaces for the Hgg fits

/*! 

This class takes in an output tree from the HggSelector and makes a RooWorkspace for doing the HggFits.
We allow either CiC (R9) or sigma_E/E selection with different maps to divide the data into 
categories.  A baseline selection is applied (though, this should be redundant, as the selection is also
applied by the selector).

The class takes one file as the data tree and an arbitrary number of signal monte-calo samples.

Author: Alex Mott (Caltech)
Date: Jan 2013

*/

#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TList.h>
#include <TObjArray.h>
#include <TH3F.h>

//All RooFit includes
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include "RooDataHist.h"
#include "RooLinkedListIter.h"
#include "RooCategory.h"

#include <HggOutputReader2.h>
#include <MixSpinDatasets.h>

#include <iostream>
#include <map>
#include <vector>


class MakeSpinWorkspace{
public:
  MakeSpinWorkspace(TString outputFileName); //!< Constructor requires the name of the file to which the output fill be written
  ~MakeSpinWorkspace();
  static TH2F* getSelectionMap(int map, bool isData); //!< takes the index of the map required and returns a TH2F* map of the selection cut as a function of pT and eta
  static int   passSelection(TH2F* map,float sigEoE, float etaSC, float pt);  //!< Takes a selection map and photon information and returns the category to which the photon is assigned
  static int   passSelection(float r9); //!< For the CiC (R9) selection: takes the r9 of the photon and returns the photon's category
  static bool getBaselineSelection(HggOutputReader2* h,int maxI,int minI); //!< Takes the current state of the input tree and which indices correspond to the leading and trailing photons and returns whether the event passes the preselection
  
  const static int nCat=2; //!< currently define only two categories (for both CiC and R9 categorization)

  void addFile(TString fName,TString l,bool is); //!< takes a file name, label, and a bool specifying whether this corresponds to data and adds this to the list of files to process
  /*!
    \param fName the path to the input file
    \param l the label associated with this file
    \param is true = data, false = MC
  */
  void setRequireCiC(bool b){requireCiC=b;} //!< specify whether to require the photons to have passed the CiC selection
  bool getRequireCiC(){return requireCiC;}  //!< returns whether CiC will be required 
  
  void setSelectionMap(int i){selectionMap=i;} //!< specify the selection map to use for the sigma_E/E analysis
  int getSelectionMap(){return selectionMap;} //!< returns the selection map
  
  void setRunRange(int rl, int rh){runLow=rl; runHigh=rh;} //!< set the run range to consider for the data

  RooWorkspace *getWorkspace(){return ws;} //!< returns the workspace in its current state

  void MakeWorkspace(); //!< Runs the workspace maker on all input files and saves the resulting workspace

  //static void runOnAll( void (*f)(TString,TString), TString mcName);
  //static void runOnAllR9( void (*f)(TString,TString), TString mcName);

  void setUseR9(bool b){useR9 = b;} //!< Specify whether to use the CiC (R9) categorization
  void setTightPt(bool b){tightPt = b;} //!< Specify whether to use the tight pt/m cuts

  void setEfficiencyCorrectionFile(TString f){EfficiencyCorrectionFile = f;}

  void setMixDatasets(){mixer = new MixSpinDatasets(ws);}
  MixSpinDatasets* getMixer(){return mixer;}

protected:
  std::vector<TString> fileName,label; //!< lists of input file names and corresponding labels
  std::vector<bool> isData;            //!< list of bools specifying whether the files correspond to data
  RooWorkspace *ws;                    //!< workspace for output
  RooCategory* labels;                 //!< list of labels to store inside the RooWorkspace
  bool requireCiC;                     //!< whether to require the photons to pass CiC (default: true)
  bool tightPt;                        //!< whether to require tight pt/m cuts (default: false)
  int selectionMap;                    //!< selection map number to use

  int runLow,runHigh;                  //!< run range to use (default 0-999999)

  TFile *outputFile;                   //!< pointer to output file

  bool useR9;                          //!< whether to use CiC (R9) or sigma_E/E cateogries (default: false)

  void AddToWorkspace(TString inputFile,TString tag, bool isData); //!< takes a file and its labels and adds to the workspace

  TString EfficiencyCorrectionFile;

  float getEffWeight(TFile *effFile, float eta, float pt, float phi, float r9); //!< returns the weight for this photon from the efficiency correction file

  MixSpinDatasets *mixer;

};

#endif
