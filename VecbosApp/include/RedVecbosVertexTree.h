#ifndef RedVecbosVertexTree_h
#define RedVecbosVertexTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedVecbosVertexTree {
public:
   RedVecbosVertexTree(const char * filename = "vecbosVertexTree.root");
  ~RedVecbosVertexTree();

  //! add the CSA07 processID and weight block
  void addCSA07Infos();
  //! add some MC truth infos
  void addMcTruthInfos();
  //! event by event final dataset fill
  void fillAll(int nj, float dxyPv, float dxyOr, float dxyErr, float dszPv, float dszOr, float dszErr, float vtxX, float vtxY, float vtxZ, float pvx, float pvy, float pvz, float dzErr, float trackerIso, float ecalIso, float hcalIso, float wtmass, float met, int chargeEle);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! fill the variable to distinguish W->enu prompt decays from the rest
  void fillMcTruth(int wtoenudecay);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  int   myNJets;       
  float myTrackDxyPV, myTrackDxyOr, myTrackDxyError; 
  float myTrackDszPV, myTrackDszOr, myTrackDszError; 
  float myTrackVertexX, myTrackVertexY, myTrackVertexZ; 
  float myPVx, myPVy, myPVz;
  float myTrackDzError;
  float myTrackerIso, myEcalIso, myHcalIso;
  float myWTransvMass, myMet;
  int myChargeEle;
  double myWeight, myProcesId;
  float myLumi;
  int myWToENuDecay;

  TFile* myFile;
  TTree* myTree;

};

#endif // RedVecbosVertexTree_h
