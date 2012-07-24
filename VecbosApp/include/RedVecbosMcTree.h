#ifndef RedVecbosMcTree_h
#define RedVecbosMcTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedVecbosMcTree {
public:
   RedVecbosMcTree(const char * filename = "vecbosMcTree.root");
  ~RedVecbosMcTree();

  //! event by event final dataset fill
  void fillAll(int procId, float id1, float eta1, float phi1, float ene1, float id2, float eta2, float phi2, float ene2); 
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  int   myProcessId;
  float myLept1Id;  
  float myLept1Eta;  
  float myLept1Phi;  
  float myLept1Ene;  
  float myLept2Id;  
  float myLept2Eta;  
  float myLept2Phi;  
  float myLept2Ene;  

  TFile* myFile;
  TTree* myTree;

};

#endif // RedVecbosTree_h
