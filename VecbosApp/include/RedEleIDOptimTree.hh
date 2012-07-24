#ifndef RedEleIDOptimTree_h
#define RedEleIDOptimTree_h
class TFile;
class TTree;

class RedEleIDOptimTree {
public:
  RedEleIDOptimTree(const char * filename = "eleID.root");
  ~RedEleIDOptimTree();

  //! event by event final dataset fill
  void fillAll(int classification, int recoflag, float eta, float pt, 
               float deta, float dphi, float hoe, float s9s25, float see, float eopout, float fbrem);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:

  float myDeltaEta;
  float myDeltaPhi;
  float myHoe;
  float myS9s25;
  float mySee;
  float myEopOut;
  float myFBrem;
  int myClassification;
  int myRecoFlag;
  float myEta;
  float myPt;

  TFile* myFile;
  TTree* myTree;
};

#endif // RedEleIDOptimTree_h
