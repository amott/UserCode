//-------------------------------------------------------
// Description:
//    Search of RS gravitons
// Authors:
//
//-------------------------------------------------------

/// The RSAnalysis class is dedicated to perform searches or RS gravitons.
/// It is derived from VecbosBase and it is based on the same
/// structure of input ntuples. It includes all the tools that
/// are common to RS analyses. In particular, ot os focused to
/// RS gravitosn -> ZZ with at least one of the two Z's decaying to
/// high energetic jets

#ifndef RSAnalysis_h
#define RSAnalysis_h

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

class RSAnalysis : public Vecbos{
public:

  RSAnalysis(TTree *tree=0); /// Class Constructor
  virtual ~RSAnalysis();     /// Class Destructor

  void Loop(string outfile);
 
private:

  bool _verbose; /// verbosity of the printouts

  vector<CaloTower> c_uncorr;
  vector<CaloTower> c_uncorr_v;
  vector<CaloTower> c_fixed;
  vector<CaloTower> c_var;
  
  
  vector<Jet> j_uncorr;
  vector<Jet> j_corr_0;
  vector<Jet> j_corr_1;
  vector<Jet> j_corr_2;
  vector<Jet> j_corr_3;
  
  vector<Jet> j_uncorr_all;
  vector<Jet> j_corr_0_all;
  vector<Jet> j_corr_1_all;
  vector<Jet> j_corr_2_all;
  vector<Jet> j_corr_3_all;
  
};
#endif
