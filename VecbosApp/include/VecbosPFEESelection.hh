//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed W-Z + jets
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef VecbosPFEESelection_h
#define VecbosPFEESelection_h

#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TMVA/Reader.h"

#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "include/Vecbos.hh"
#include "VecbosBase.hh"
#include "include/RedVecbosPFTree.h"
#include "include/RedEleIdIsolPFTree.hh"
#include "include/McTruthEvent.hh"
#include "include/CutBasedSelectorEE.hh"

class VecbosPFEESelection : public Vecbos {
public:
  
  //! constructor
  VecbosPFEESelection(TTree *tree=0);
  //! destructor
  virtual ~VecbosPFEESelection();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies();
  //! set the list of the required triggers
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  //! use a prefix for the output 
  void setPrefix(const char *prefix) { m_prefix = prefix; }
  //! set the signal type to select based on MC-truth
  void setSignal(int signal) { m_signal = signal; }

private:

  //! return the two highest pt electrons  
  std::pair<int, int> getBestGoodElePair(std::vector<int> goodElectrons);
  //! set the event kinematics
  void setKinematics(int theEle1, int theEle2);
  //! set the best (di)-electron variables
  void setElectrons();
  //! compute invariant mass
  float getMee(int theEle1, int theEle2);

  //! make a list of SIS cone jets from calotowers sorted in Pt
  std::vector<Jet> buildSortedSISConeCaloJets();
  
 //! returns the output of the custom cut electron ID
  bool isEleID(int eleIndex, bool wantTight);
  //! returns the list of the calotowers 
  std::vector<CaloTower> getCaloTowers(float zpv, int type);
  //! returns true if the W electron combined with a reco electron makes Z mass
  bool foundZCandidate(int eleIndex, std::vector<int> accEles);

  //! fill the tree for electron ID / isolation optimization
  void FillEleIdIsolPFTree(int eleIndex, int eleid, int isol, int conv);

  //! selections common to all analyses
  void ConfigCommonSelections(Selection* _selection);
  //! selections for Z->ee analysis
  void ConfigZeeSelections(Selection* _selection, Counters *_counter);
  void displayZeeEfficiencies(Counters _counter);
  //! selections for W->enu analysis
  void ConfigWenuSelections(Selection* _selection, Counters *_counter);
  void displayWenuEfficiencies(Counters _counter);

  //! methods for B veto
  bool isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits);
  vector<float> jetBTagVariables(Jet jet);
  float jetDxy(Jet jet);
  float jetDsz(Jet jet);
  float invMassTracks(Jet jet);
  int tracksMultiplicity(Jet jet);
  float fractionalEnergy(Jet jet);
  void calcEventBVetoVariables(std::vector<Jet> jets, std::vector<BTagJet> btags);
  
  static const int NMET = 50;
  static const int NMu = 50;

  //! the members
  int theEle1_, theEle2_;
  TLorentzVector p4Ele1_, p4Ele2_;
  TVector3 pT3Ele1_, pT3Ele2_, p3Met_, U_, pt3W_, pt3Z_, sumCalo_;
  TVector3 p3CaloMet_, p3TCMet_, p3PFMet_;
  float WCalomT_, WmT_, WTCmT_, WPFmT_, mee_;
  int closestPV;

  //! count muons in the event
  int CountMuons(float etaMax, float ptMin);

  //! be verbose during runtime
  bool _verbose;
  //! the list of required triggers
  vector<int> m_requiredTriggers;
  //! to compute final quantities
  int theBestZee; 
  //! vectors to store indices of best candidates
  std::vector<int> goodElectrons;

  //! vector of isolation variables for reconstructed electrons
  std::vector<float> chargedIsolation, photonIsolation, neutralIsolation;
  //! input variables to isolation Fisher 
  float chargedIsolRel, photonIsolRel, neutralIsolRel;
  //! vector of combined isolation for reconstructed electrons
  std::vector<float> combinedIsolation;

  //! vectors of cleaned jets
  std::vector<Jet> goodWJets_, goodZJets_, goodWPFJets_, goodZPFJets_;
  //! vectors of cleaned jets for b-veto
  std::vector<Jet> m_cleanedBvetoWJets, m_cleanedBvetoZJets, m_cleanedBvetoWPFJets, m_cleanedBvetoZPFJets;
  std::vector<BTagJet> m_cleanedBvetoWBTagJets, m_cleanedBvetoZBTagJets, m_cleanedBvetoWBTagPFJets, m_cleanedBvetoZBTagPFJets;
  //! EVT b-tagging
  BTagJet m_maxBTagEvt;

  //! the CMS calorimeters
  enum { ecal=0, hcal=1 };

  //! the processes
  enum { wjets=0, zjets=1, wother=2, zother=3, all=4 };
  int m_signal;

  //! the calotowers of the event
  std::vector<CaloTower> _theCaloTowers;
  std::vector<CaloTower> _theCaloTowersForShape;
  std::vector<CaloTower> _theCaloTowersForPFJets;
  
  //! the SISCone jets of the event
  std::vector<Jet> _theSISConeCaloJets;
  std::vector<Jet> _theSISConePFJets;

  //! the antiKT jets of the event
  std::vector<Jet> _theAK5CaloJets;
  std::vector<Jet> _theAK5PFJets;
  std::vector<Jet> _theAK5GenJets;

  //! the btag collection (only for calojets, so far)
  std::vector<BTagJet> _theSISConeBTagCaloJets;
  std::vector<BTagJet> _theAK5BTagCaloJets;
 
  //!counters and selections
  Selection *_commonSel,       *_zeeSel,         *_wenuSel; 
  Counters _zeeCounter,      _wenuCounter;    
  Counters _zeeCounterZjets, _wenuCounterZjets; 
  Counters _zeeCounterWjets, _wenuCounterWjets; 
  Counters _zeeCounterTTbar, _wenuCounterTTbar; 
  Counters *_zeePCounter,      *_wenuPCounter;    
  Counters *_zeePCounterZjets, *_wenuPCounterZjets; 
  Counters *_zeePCounterWjets, *_wenuPCounterWjets; 
  Counters *_zeePCounterTTbar, *_wenuPCounterTTbar; 

  //! number of reco jets (with reco electrons removal)
  int howManyWJets, howManyZJets;
  int howManyForEffWJets, howManyForEffZJets;
  int howManyWPFJets, howManyZPFJets;
  int howManyForEffWPFJets, howManyForEffZPFJets;

  //! counters by jet multiplicity for signal efficiency (with MC truth electrons removal)
  Counters m_wenuJetbinCounter[6], m_zeeJetbinCounter[6];
  Counters m_wenuRecoJetbinCounter[6], m_zeeRecoJetbinCounter[6];
  Counters m_wenuPFJetbinCounter[6], m_zeePFJetbinCounter[6];
  Counters m_wenuRecoPFJetbinCounter[6], m_zeeRecoPFJetbinCounter[6];
  
  //! reduced trees
  RedVecbosPFTree    *myOutTree_Zee;
  RedVecbosPFTree    *myOutTree_Wenu;
  RedEleIdIsolPFTree *myOutEleIdIsolPFTree_Wenu;

  //! name of rootfile with dataset
  std::string _datasetName;

  //! WToENuDecay = 1 if W->e nu promptly (exclude other W decays, as W->tau nu)
  int WToENuDecay;
  int ZToEEDecay;
  //! pthat for the case of photon+jet
  float photonj_pthat;
  //! the prefix to add to the output files
  const char *m_prefix;
  char namefile[200];
  //! shower shape variables 
  vector<TLorentzVector> CTinput;

  //! counters to compute number of jets at the beginning of the selection
  double m_wInitJets[6];
  double m_zInitJets[6];
  double m_wEndJets[6];
  double m_zEndJets[6];
  int nwjets_mc, nzjets_mc, nwpfjets_mc, nzpfjets_mc, nwgenjets_mc, nzgenjets_mc;

  //! best (di)-electrons variables
  int recoflag_[2];
  float pt_[2], eta_[2];
  float mva_[2];
  float chargedIso_[2], neutralIso_[2], photonIso_[2], combinedIso_[2];
  int charge_[2];
  
  int isData_;
  McTruthEvent mcevent;

};
#endif
