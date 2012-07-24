//-------------------------------------------------------
// Description:
//    Class for selection of top events
//-------------------------------------------------------

#ifndef TopControlSample_h
#define TopControlSample_h

#include <vector>
#include <TVector3.h>
#include <TMath.h>
#include <TF2.h>
#include <TLorentzVector.h>
#include "TMVA/Reader.h"

#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "include/Vecbos.hh"
#include "VecbosBase.hh"
#include "include/McTruthEvent.hh"
#include "include/CutBasedSelectorEE.hh"
#include "include/RedTopTree.h"

class TopControlSample : public Vecbos {
public:
  
  //! constructor
  TopControlSample(TTree *tree=0);
  //! destructor
  virtual ~TopControlSample();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies();
  //! use a prefix for the output 
  void setPrefix(const char *prefix) { m_prefix = prefix; }
  //! set the signal type to select based on MC-truth
  void setSignal(int signal) { m_signal = signal; }

private:

  //! set the MC truth tree
  void setMcTruthTree();
  //! return the two highest pt electrons or muons 
  std::pair<int,int> getBestGoodElePair(std::vector<int> goodElectrons);
  std::pair<int,int> getBestGoodMuonPair(std::vector<int> goodMuons);
  //! set the event kinematics
  void setKinematics2Ele(int theEle1, int theEle2);
  void setKinematics2Mu(int theMu1, int theMu2);
  void setKinematics1Ele1Mu(int theEle1, int theMu1);

  //! kine selection and counter
  void ConfigureSelection(Selection* _selection, Counters *_counter);
  Selection *_selection;  
  Counters  _counterEE, _counterMM, _counterEM;
  Counters *_pCounterEE, *_pCounterMM, *_pCounterEM;
  
  //! efficiencies
  void displayRecoEfficiencies(Counters _counter);

  //! returns the output of the custom cut electron or muon ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);
  void isMuonID(int muonIndex, bool *muonIdOutput);
  double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);

  //! projected met
  float GetProjectedMet(TVector3 p1, TVector3 p2);

  //! to check bjets
  bool countBTagJets(Jet myJet);
  int countAllBTagJets(std::vector<Jet> jets);
  int foundBquarks();
  int foundBhadrons();

  //! reduced tree
  RedTopTree *myOutTree;

  // max
  static const int NMET = 50;
  static const int NMu = 50;

  //! the members
  TLorentzVector p4Ele1_, p4Ele2_;
  TLorentzVector p4Mu1_,  p4Mu2_;
  TVector3 pT3Ele1_, pT3Ele2_;
  TVector3 pT3Mu1_,  pT3Mu2_;
  TVector3 p3Met_;
  float WmT_, invMass_, projectedMet_;
  int closestPV, hardestPV;

  // reco channel
  bool m_channelEE, m_channelMM, m_channelEM;
  bool m_channelEE_soft, m_channelMM_soft, m_channelEM_soft;

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;
  CutBasedEleIDSelector EgammaTightID_NC;

  //! be verbose during runtime
  bool _verbose;

  //! vector of isolation variables for reconstructed electrons
  float trackerIsolRel, ecalIsolRel, hcalIsolRel;
  std::vector<float> trackerIsolation, ecalIsolation, hcalIsolation;
  std::vector<float> combinedIsolation;

  //! vectors of cleaned jets
  std::vector<Jet> goodPFJets_;

  //! the CMS calorimeters
  enum { ecal=0, hcal=1 };
  
  //! the processes
  enum { wjets=0, zjets=1, wother=2, zother=3, all=4 };
  int m_signal;
  
  //! the antiKT jets
  std::vector<Jet> _theAK5PFJets;
  std::vector<BTagJet> _theAK5BTagPFJets;

  //! number of reco jets (with reco leptons removal)
  int howManyPFJets;

  //! number of b-tag jets with different thresholds
  int nTagWPFJets[5];

  //! name of rootfile with dataset
  std::string _datasetName;

  //! the prefix to add to the output files
  const char *m_prefix;
  char namefile[200];

  //! best (di)-leptons variables
  float pt_[2], eta_[2], charge_[2];

  bool isData_;
  McTruthEvent mcevent;
  int TopToEEDecay, TopToMMDecay, TopToEMDecay;

};
#endif
