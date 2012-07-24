#ifndef RedTopTree_h
#define RedTopTree_h

#include <vector>

class TFile;
class TTree;

class G3EventProxy;

class RedTopTree {
public:

   RedTopTree(const char * filename = "vecbosTree.root");

  ~RedTopTree();
  
  // general infos
  void fillGeneral(int cEE, int cMM, int cEM, float p1, float p2, float e1, float e2, int c1, int c2, float im, int npfj, float met, float ml1, float ml2);

  // lepton compatibility with vtx
  void fillVertexComp(float de, float dm);
  
  //! fill the jet-wise btag output (for one algo: 5 WPs)
  void fillBTagJetMultiplicities(int nJetsBTagged[5]);
  void fillBTag(float etj1, float etj2, int fb, int fbh, int ma, int mj1, int mj2, int tj1, int tj2);
  
  //! store run, event, lumi
  void fillRunInfo(int run, int lumi, int event);

  //! effectively store the events in the tree
  void store();

  //! save in the ROOT file
  void save();

private:
  int   myRun, myLumi, myEvent;

  int myChannelEE, myChannelMM, myChannelEM;
  float myPt1, myPt2;
  float myEta1, myEta2;
  int myCharge1, myCharge2;
  float myInvMass, myMet;
  int myNPFJets;
  float myMcLeptDr1, myMcLeptDr2;

  int myNJetsBTagged[5];

  float myEtJet1, myEtJet2;
  int myFoundB, myFoundBhad;
  int myMatchedAll;
  int myMatchedJet1, myMatchedJet2;
  int myTaggedJet1, myTaggedJet2;

  float myDzEle, myDzMu;

  TFile* myFile;
  TTree* myTree;
};

#endif // RedTopTree_h
