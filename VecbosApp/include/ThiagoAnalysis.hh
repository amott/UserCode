//-------------------------------------------------------
// Description:
//    Class for Alpgen Validation
// Authors:
//    Thiago Tomei & VecBos Team
//-------------------------------------------------------

/// The ThiagoAnalysis class is compatible with Thiago's
/// analysis in Chowder soup back in CSA07. It defines 
/// three tiers of requirements, and plots some basic 
/// quantities.

#ifndef ThiagoAnalysis_h
#define ThiagoAnalysis_h

#include "Vecbos.hh"

class ThiagoAnalysis : public Vecbos{
public:

  ThiagoAnalysis(TTree *tree=0); /// Class Constructor
  virtual ~ThiagoAnalysis();     /// Class Destructor
  /// The function to run on each events
  void ThiagoAnalysis::Loop(char* filename);
    
private:

  /// Global data members that change from event to event;
  vector<int> globalLeadingJets;
  vector<int> globalLeadingMuons;
  
  /// Create a set of histograms. Takes as input the name of the 
  /// directory were to write the histograms, also used to give 
  /// names to the histograms (to avoid memory problems if used
  /// for more than a set of histograms)
  vector<TH1D*> CreateHistos(string dirname); 

  /// Fill the histograms passed as input
  void FillHistos(vector<TH1D*> histos);

  /// Overloading from Vecbos
  virtual void AssignParameters(map<string, double>);
  virtual void InitParameters();

  /// Jet minimum pT
  double jetPtCut;
  /// Muon minimum eT
  double muonPtCut;
  /// Muon acceptance
  double muonAcceptance;

  /// Gives the number of jets in this event.
  int NumberOfJets();
  /// Makes baseline cut in leading jet.
  bool BaselineCut(int leadingJet);
  /// Makes selection cut in leading muon.
  bool SelectionCut(int leadingMuon);

  /// Gets true Z.
  int getTrueZ();
  /// Gets true W.
  int getTrueW();
  /// Gets true top.
  int getTrueTop();
  /// Gets best reconstructed Z.
  int getBestZ();
  /// Gets collection of jets from the event;
  void theSortedJets(vector<int>& jets);
  /// Gets collection of muons from the event;
  void theSortedMuons(vector<int>& muons);
  /// Gets reconstructed boson from the jets;
  TLorentzVector getBosonFromJets();
  
};
#endif
