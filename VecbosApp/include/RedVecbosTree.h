#ifndef RedVecbosTree_h
#define RedVecbosTree_h

#include <vector>

class TFile;
class TTree;

class G3EventProxy;

class RedVecbosTree {
public:
   RedVecbosTree(const char * filename = "vecbosTree.root");
  ~RedVecbosTree();

  //! add some MC truth infos
  void addMcTruthInfos();
  //! add variables to compute PF ele efficiencies
  void addPFelectrons();
  //! add the EVT full btag output
  void addBTagEVTInfos();
  //! add the informations on MHTphiMET and MHTphiJet
  void addEventShapeInfos();
  //! add the electron ID+iso variables for the selected best electrons
  void addElectronInfos();
  //! event by event final dataset fill
  void fillJetMultiplicities(int nJetsHi, int nPFJetsHi, int nTrackJetsHi, int nGenJetsHi, int nJetsLo, int nPFJetsLo, int nTrackJetsLo, int nGenJetsLo, int caloJetSelectedHi, int PFJetSelectedHi, int trackJetSelectedHi, int caloJetSelectedLo, int PFJetSelectedLo, int trackJetSelectedLo, float leadingCaloPt, float leadingPFPt, float emFracCaloJet, float emFracPFJet);
  void fillPFelectrons(float pfPt, float pfEta, float pfPhi, float pfdR, float pfmva, int npfe);
  void fillElectrons(int recoflag[2], float pt[2], float eta[2], float phi[2], int classification[2], int nbrems[2], float deta[2], float dphi[2], float hoe[2], float see[2], float spp[2], float eop[2], float fbrem[2], float trackerIso[2], float hcalIso[2], float ecalJIso[2], float ecalGTIso[2], float combinedIso[2], int charge[2], int missHits[2], float dist[2], float dcot[2], float lh[2], float e9esc[2],float e25esc[2]);
  void fillMoreElectrons(float esc[2], float eseed[2], float pin[2]);
  void fillKinematics(float invmass, float tmass, float tctmass, float pftmass, float met, float tcmet, float pfmet, float ptW, float ptWrecoil, float ptZ, float dijetmass[2], float dijetpt[2], float deltaetajet[2], float deltaphimet[2], float etaW, float etaZ, float mTZMet, float mZJet);  
  void fillWTemplates(float metCorr, float metUncorr, float mTCorr, float mTUncorr);
  void fillEventShape(std::vector<float> mhtphijet, std::vector<float> mhtphiPFjet, float mhtphimet);
  //! fill the variable to distinguish W->enu prompt decays from the rest
  void fillMcTruth(int promptDecay, int nInteractions=1);
  //! fill the EVT full btag outputs
  void fillBTagEVT(float combinedSecondaryVertexBJetTags, float combinedSecondaryVertexMVABJetTags,
                   float jetBProbabilityBJetTags, float jetProbabilityBJetTags, float simpleSecondaryVertexBJetTags,
                   float softMuonBJetTags,
                   float trackCountingHighPurBJetTags, float trackCountingHighEffBJetTags, 
		   int fB, int fBm, int nB);
  //! fill the jet-wise btag output (for one algo: 5 WPs)
  void fillBTagJetMultiplicities(int nJetsBTagged[5]);
  //! check  MC acceptance
  void fillGenInfo(float p1, float p2, float e1, float e2, float mm);
  //! store run, event, lumi
  void fillRunInfo(int run, int lumi, int event);
  //! store event variables
  void fillEventInfo(int nvtx, float rho);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

 private:

  int   myRun, myLumi, myEvent;
  int   myNJetsHi, myNJetsLo;       
  int   myNTrackJetsHi, myNTrackJetsLo;       
  int   myNPFJetsHi, myNPFJetsLo;
  int   myNGenJetsHi, myNGenJetsLo;
  int   myNMuon;
  float myInvMass;
  float myTransvMass, myTCTransvMass, myPFTransvMass;
  float myMet, myTCMet, myPFMet;
  double myWeight;
  int myWToENuDecay;
  float myPtW, myPtWRecoil, myEtaW;
  float myPtZ, myEtaZ, myMTZMet, myMZJet;
  float myMHTphiJet[6];        
  float myMHTphiPFJet[6];        
  float myMHTphiMET;        
  float myTemplMetUncorr, myTemplMetCorr;
  float myTemplTrMassUncorr, myTemplTrMassCorr;
  float myLeadingCaloJetPt, myLeadingPFJetPt;
  float myEmFracCaloJet, myEmFracPFJet;
  int  myCaloJetSelectedHi, myPFJetSelectedHi, myTrackJetSelectedHi,
    myCaloJetSelectedLo, myPFJetSelectedLo, myTrackJetSelectedLo;
  int myPromptDecay;
  float myDiJetMass[2]; // hi,lo jet ET
  float myDiJetPt[2];
  float myDeltaEtaJet[2];
  float myDeltaPhiJetMet[2];
  int myNumInteractions;
  int myNumVtx;
  float myRho;

  // electron variables
  int myRecoflag[2];
  float myPt[2], myEta[2], myPhi[2];
  int myClassification[2], myNBrems[2];
  float myDeta[2], myDphi[2], myHoe[2], mySee[2], mySpp[2], myEop[2], myFbrem[2];
  float myTrackerIso[2], myHcalIso[2], myEcalJIso[2], myEcalGTIso[2], myCombinedIso[2];
  int myCharge[2];
  int myMissHits[2];
  float myDist[2], myDcot[2];
  float myLh[2];
  float myE9ESC[2],myE25ESC[2];
  float myESC[2], myEseed[2], myPin[2];

  // to compute PF electron efficiencies
  float myPFelePt, myPFeleEta, myPFelePhi;
  float myPFeleDeltaR, myPFeleMva;
  int myNPFele;

  float myCombinedSecondaryVertexBJetTags;
  float myCombinedSecondaryVertexMVABJetTags;
  float myJetBProbabilityBJetTags;
  float myJetProbabilityBJetTags;
  float mySimpleSecondaryVertexBJetTags;
  float mySoftMuonBJetTags;
  float myTrackCountingHighPurBJetTags;
  float myTrackCountingHighEffBJetTags;
  int myFoundMcB;
  int myFoundMcBmum;
  int myNumB;
  int myNJetsBTagged[5];

  float myPtGen1,  myPtGen2;
  float myEtaGen1, myEtaGen2;
  float myMeeGen;

  TFile* myFile;
  TTree* myTree;

};

#endif // RedVecbosTree_h
