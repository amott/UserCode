//-------------------------------------------------------
// Description:
//    Search for RS gravitons
// Authors:
//    Thiago Tomei, SPRACE
//-------------------------------------------------------

/// The RSZZAnalysis class is dedicated to searches for RS gravitons.
/// It is derived from VecbosBase and it is based on the same
/// structure of input ntuples. It includes all the tools that
/// are common to RS analyses. In particular, it is focused to
/// RS gravitons -> ZZ with at least one of the two Z's decaying to
/// high energetic jets

#ifndef RSZZAnalysis_h
#define RSZZAnalysis_h

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <vector>

#include "Vecbos.hh"

using namespace std;

class RSZZAnalysis : public Vecbos{
public:

  RSZZAnalysis(TTree *tree=0); /// Class Constructor
  virtual ~RSZZAnalysis();     /// Class Destructor

  void Loop(string outfile);
 
private:

  vector<float> caloThresholds; ///< Calorimetric thresholds.

  /// Minimum number of jets.
  int minNumJets;
  /// Jet eT cuts.
  double jetEtCut, j1EtCut, j2EtCut;
  /// Jet mass cuts.
  double j1MinMass, j1MaxMass, j2MinMass, j2MaxMass;
  /// Track pT cut.
  double trackPtCut;
  /// Jet number of tracks cut.
  int jetTracksCut;
  /// Jet radius
  double jetRadius;
 
  /// Overloading from VecBos.
  virtual void AssignParameters(map<string, double>);
  virtual void AssignParameters(map<string, string>);
  virtual void InitParameters();

  /// Pointer to members.
  TLorentzVector* graviton;
  TLorentzVector* Z1;
  TLorentzVector* Z2;
  double* theMET;

  /// Useful typedef - associates a Jet with a vector of TLorentzVectors - tracks, calotowers...
  typedef pair<Jet,vector<TLorentzVector> > JetWithAssociation;
  
  /// Pointer to particles and jet collections.
  vector<Jet>* pGenJets;
  vector<Jet>* pGoodJets;
  vector<TLorentzVector>* pCharged;
  vector<TLorentzVector>* pGoodTracks;

  /// Files and histograms.
  TFile* theFile;
//   TH1D* flowJet1;
//   TH1D* antiFlowJet1;
//   TH1D* aFOverFJet1;
//   TH1D* flowJet1_AC;
//   TH1D* antiFlowJet1_AC;
//   TH1D* aFOverFJet1_AC;
//   TH1D* flowJet2;
//   TH1D* antiFlowJet2;
//   TH1D* aFOverFJet2;
//   TH1D* flowJet2_AC;
//   TH1D* antiFlowJet2_AC;
//   TH1D* aFOverFJet2_AC;
  
//   vector<TH1D*> histosJet1;
//   vector<TH1D*> histosJet2;
//   vector<TH2D*> histos2DJet1;
//   vector<TH2D*> histos2DJet2;
//   vector<TH1D*> histosJet1_AC;
//   vector<TH1D*> histosJet2_AC;
//   vector<TH2D*> histos2DJet1_AC;
//   vector<TH2D*> histos2DJet2_AC;
  TH2D* H_jm1Xjm2;
  TH1D* H_j1_manytracks_et;
  TH1D* H_j1_manytracks_eta;
  TH1D* H_j2_manytracks_et;
  TH1D* H_j2_manytracks_eta;

  /// Function to create histograms specific for my fat jets.
  vector<TH1D*> CreateJetHistos(string dirname);
  /// Function to fill histograms for my fat jets.
  void FillJetHistos(vector<TH1D*>& theHistos, JetWithAssociation& theJet);
  
  /// Function to create histograms for a TLorentzVector*.
  vector<TH1D*> CreateBasicHistos(string dirname);
  /// Function to fill the histograms.
  void FillBasicHistos(vector<TH1D*>& theHistos, TLorentzVector* theCandidate);
  
  /// Function to create 2D histograms specific for my fat jets.
  vector<TH2D*> Create2DHistos(string dirname);
  /// Function to fill 2D histograms for my fat jets. 
  void Fill2DHistos(vector<TH2D*>& theHistos, JetWithAssociation& theJet);

  /// Function to get the hardest PV in the event.
  int GetHardestPV(double threshold = 0);
  /// Function to filter jets given a et cut.
  vector<Jet> FilterJetsByEt(vector<Jet> theJets, double jetEtCut);
  /// Function to count tracks around a TLorentzVector.
  int RSZZAnalysis::NumTracks(TLorentzVector v, double r);
  /// Function to count charged particles around a TLorentzVector.
  int RSZZAnalysis::NumCharged(TLorentzVector v, double r);
};

#endif
