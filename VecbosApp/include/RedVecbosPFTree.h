#ifndef RedVecbosPFTree_h
#define RedVecbosPFTree_h

#include <vector>

class TFile;
class TTree;

class G3EventProxy;

class RedVecbosPFTree {
public:
   RedVecbosPFTree(const char * filename = "vecbosTreePF.root");
  ~RedVecbosPFTree();

  //! add some MC truth infos
  void addMcTruthInfos();
  //! add the electron ID+iso variables for the selected best electrons
  void addElectronPFInfos();
  void fillJetMultiplicities(int nJets, int nPFJets, int nGenJets, int caloJetSelected, int PFJetSelected);
  void fillElectronsPF(float pt[2], float eta[2], float mva[2], float ci[2], float ni[2], float pi[2], float combIso[2], int charge[2]);
  void fillKinematicsPF(float calotmass, float tctmass, float pftmass, float calomet, float tcmet, float pfmet);
  void fillMcTruth(int promptDecay);
  void fillRunInfo(int run, int lumi, int event);
  // void fillArcs(float gpt);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  int   myRun, myLumi, myEvent;
  int   myNJets;       
  int   myNPFJets;
  int   myNGenJets;
  int   myNMuon;
  float myCaloTransvMass, myTCTransvMass, myPFTransvMass;
  float myCaloMet, myTcMet, myPFMet;
  double myWeight;
  int myWToENuDecay;
  int  myCaloJetSelected, myPFJetSelected;
  int myPromptDecay;
  // float myGenWpt;

  // electron variables
  float myPFPt[2], myPFEta[2];
  float myPFMva[2], myPFChargedIso[2];
  float myPFNeutralIso[2], myPFPhotonsIso[2];
  float myPFCombinedIso[2];
  int myPFCharge[2];

  TFile* myFile;
  TTree* myTree;

};

#endif // RedVecbosTree_h
