#ifndef RedEleIdIsolPFTree_h
#define RedEleIdIsolPFTree_h
class TFile;
class TTree;

class RedEleIdIsolPFTree {
public:
  RedEleIdIsolPFTree(const char * filename = "eleID.root");
  ~RedEleIdIsolPFTree();

  //! event by event final dataset fill
  void fillAll(float pt, float eta, float mva, float hoe, float deta, float dphi, float eop, float eopo, float see, float ch, float nh, float ph, float comb, int il, float tip, float tisol, float eisol, float hisol);

  void fillIsol(float cnv3, float cnv4, float cnv5, float cv3, float cv4, float cv5, 
		float nvcnv3, float nvcnv4, float nvcnv5, float nvcv3, float nvcv4, float nvcv5, 
		float nnv3, float nnv4, float nnv5, float nv3, float nv4, float nv5, 
		float pnv3, float pnv4, float pnv5, float pv3, float pv4, float pv5);

  void fillSetpsAfter(int eleId, int isol, int conv);
  
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:

  float myEta, myPt;
  float myMva, myHoE;
  float myDeltaEta, myDeltaPhi;
  float myEoP, myEoPout;
  float mySigmaIeIe;
  float myCharged, myNeutral, myPhoton, myCombined;
  float myTrackerIsol, myEcalIsol, myHcalIsol;
  int   myInnerLayers;
  float myTranImpPar;

  int myPassedEleId, myPassedIsol, myPassedConv;

  float myEleChargedIso03nV, myEleChargedIso04nV, myEleChargedIso05nV;
  float myEleChargedIso03v,  myEleChargedIso04v,  myEleChargedIso05v;
  
  float myEleChargedIsoNvc03nV, myEleChargedIsoNvc04nV, myEleChargedIsoNvc05nV;
  float myEleChargedIsoNvc03v,  myEleChargedIsoNvc04v,  myEleChargedIsoNvc05v;
  
  float myEleNeutralIso03nV, myEleNeutralIso04nV, myEleNeutralIso05nV;
  float myEleNeutralIso03v,  myEleNeutralIso04v,  myEleNeutralIso05v;

  float myElePhotonsIso03nV, myElePhotonsIso04nV, myElePhotonsIso05nV;
  float myElePhotonsIso03v,  myElePhotonsIso04v,  myElePhotonsIso05v;
        
  TFile* myFile;
  TTree* myTree;
};

#endif // RedEleIDOptimTree_h
