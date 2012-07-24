//-------------------------------------------------------
// Description:
//    Class for GenJet search analyses
// Authors:
//
//-------------------------------------------------------

/// The EventShapeAnalysis class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef EventShapeAnalysis_h
#define EventShapeAnalysis_h

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


class EventShapeAnalysis : public Vecbos{
public:

  EventShapeAnalysis(TTree *tree=0); /// Class Constructor
  virtual ~EventShapeAnalysis();     /// Class Destructor
  /// The function to run on each events
  void EventShapeAnalysis::Loop(string outname="test.root");
    
private:
  string outfilename;
  /// Create a set of histograms. Takes as input the name of the 
  /// directory were to write the histograms, also used to give 
  /// names to the histograms (to avoid memory problems if used
  /// for more than a set of histograms)
  vector<TProfile*> CreateHistosT(string dirname); 
  vector<TH1D*>     CreateHistosD(string dirname); 
  vector<TH2D*>     CreateHistos2D(string dirname);

  /// Fill the histograms passes as input
  void FillHistos(vector<TH1D*> h_d, vector<TProfile*> h_p, vector<TH2D*> h_2D);

  double x0; ///< the X coordinate of the new vertex
  double y0; ///< the Y coordinate of the new vertex
  double z0; ///< the Z coordinate of the new vertex

  bool W_pass;
  bool Z_pass;
  int lead_mu_index;
  double lead_mu_pt;

  int z_index;

  CoolTools *CT[4];

  static const double Ztruemass = 91.1876;
  static const double Wtruemass = 80.403;

  vector<CaloTower> c_calo;
  vector<CaloTower> track_collection;
  vector<Jet> c_all;
  vector<Jet> c_all_boost;

  vector<Jet> calo_jet;
  vector<Jet> track_jet;
  vector<Jet> boost_calo_jet;
  vector<Jet> boost_track_jet;
  vector<Jet> j_dum;

  TLorentzVector MHT[4];
  TLorentzVector lead_mu;
  TLorentzVector lead_j;
  TLorentzVector lead_track;
  TLorentzVector Z;

  int isample; //0 ttbar 1 W 2 Z
  int N_calo_jet;
  int N_track_jet;

  vector<vector<TLorentzVector> > JETS;

  
};
#endif
