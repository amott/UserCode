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
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>

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
#include "assert.h"

class MakeSpinWorkspace{
public:
  MakeSpinWorkspace(TString outputFileName); //!< Constructor requires the name of the file to which the output fill be written
  ~MakeSpinWorkspace();
  static TH2F* getSelectionMap(int map, bool isData); //!< takes the index of the map required and returns a TH2F* map of the selection cut as a function of pT and eta
  int   passSelection(TH2F* map,float sigEoE, float etaSC, float pt);  //!< Takes a selection map and photon information and returns the category to which the photon is assigned
  int   passSelection(float r9); //!< For the CiC (R9) selection: takes the r9 of the photon and returns the photon's category
  bool getBaselineSelection(HggOutputReader2* h,int maxI,int minI,float mass); //!< Takes the current state of the input tree and which indices correspond to the leading and trailing photons and returns whether the event passes the preselection
  
  const static int nCat=2; //!< currently define only two categories (for both CiC and R9 categorization)

  /*!
    \param fName the path to the input file
    \param l the label associated with this file
    \param is true = data, false = MC
  */
  void addFile(TString fName,TString l,bool is); //!< takes a file name, label, and a bool specifying whether this corresponds to data and adds this to the list of files to process
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

  void setUseUncorrMass(bool b=true){useUncorrMass=b;} //!< selected whether to use the default photon energies (no cluster corrections or scaling/smearing) for the mass
  void setEfficiencyCorrectionFile(TString f_data,TString f_MC) //!< specify the file to do MC-based data efficiency corrections
  {
    EfficiencyCorrectionFile_Data = f_data;
    EfficiencyCorrectionFile_MC = f_MC;    
  } 

  void setMixDatasets(){mixer = new MixSpinDatasets(ws);} //!< specify that we will be mixing MC datasets
  MixSpinDatasets* getMixer(){return mixer;} //!< return the mixing module for customization

  void setKFactorFile(TString s){filenameKFactor = s;} //!< specify the file containing the KFactor weights for the higgs signal
  void setRescaleFile(TString s){filenameRescaleFactor = s;} //!< specify the file containing the data/MC weights for the signal MC

  void setMassRange(float min, float max){mMin=min; mMax=max;}
  void setPtCuts(float c1, float c2){pt1Min=c1; pt2Min=c2;}
  void setPt1Cut(float c){pt1Min=c;}
  void setPt2Cut(float c){pt2Min=c;}
protected:
  std::vector<TString> fileName,label; // lists of input file names and corresponding labels
  std::vector<bool> isData;            // list of bools specifying whether the files correspond to data
  RooWorkspace *ws;                    // workspace for output
  RooCategory* labels;                 // list of labels to store inside the RooWorkspace
  TFile *outputFile;                   //!< pointer to output file

  //selections
  bool requireCiC;                     // whether to require the photons to pass CiC (default: true)
  bool tightPt;                        // whether to require tight pt/m cuts (default: false)
  int selectionMap;                    // selection map number to use

  int runLow,runHigh;                  // run range to use (default 0-999999)

  float mMin, mMax;
  float pt1Min, pt2Min;

  bool useR9;                          // whether to use CiC (R9) or sigma_E/E cateogries (default: false)
  bool useUncorrMass;

  void AddToWorkspace(TString inputFile,TString tag, bool isData); // takes a file and its labels and adds to the workspace

  TString EfficiencyCorrectionFile_Data;
  TString EfficiencyCorrectionFile_MC;

  std::vector<TH3F*> effMaps;
  float getEffWeight(float eta, float pt, float phi, float r9); // returns the weight for this photon from the efficiency correction file

  MixSpinDatasets *mixer;

  float getEfficiency(HggOutputReader2 &h, int massPoint);

  TString filenameKFactor;
  TFile *fileKFactor;
  float getKFactor(HggOutputReader2 &h, int massPoint);

  TString filenameRescaleFactor;
  TFile *fileRescaleFactor;
  float getRescaleFactor(HggOutputReader2 &h);

  float getEffFromTGraph(TGraphAsymmErrors* e,float pt);


};

#endif
