//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

/// The VtxStudy class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef VtxStudy_h
#define VtxStudy_h

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <vector>

#include "VecbosBase.hh"
#include "Vecbos.hh"

using namespace std;

class VtxStudy : public Vecbos{
public:

  VtxStudy(TTree *tree=0); /// Class Constructor
  virtual ~VtxStudy();     /// Class Destructor
  /// The function to run on each events
  void Loop();
    
private:
  
  /// Create a set of histograms. Takes as input the name of the 
  /// directory were to write the histograms, also used to give 
  /// names to the histograms (to avoid memory problems if used
  /// for more than a set of histograms)
  vector<TH2D*> CreateHistos(string dirname, string part); 
  vector<TH1D*> CreateHistos1D(string dirname, string part); 

  /// Fill the histograms passes as input
  void FillHistos(vector<TH2D*> histos);
  void FillHistosEle(vector<TH2D*> histos);
  void FillHistos1D(vector<TH1D*> histos);

  double dRM;
  int iRecomu, iRecoe, iGen, njets;

  double ptmin;
  double ptmax;
  double etamin;
  double etamax;
  double dRmin;
  double dRmax;
  double dzmax;


};
#endif
