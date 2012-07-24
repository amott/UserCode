#ifndef RedIsolationVtxOptimTree_h
#define RedIsolationVtxOptimTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedIsolationVtxOptimTree {
public:
   RedIsolationVtxOptimTree(const char * filename = "isolVtxOptim.root");
  ~RedIsolationVtxOptimTree();

  //! add the CSA07 processID and weight block
  void addCSA07Infos();
  //! add the kinematics branches
  void addKinematicsInfos();
  //! fill electrons kinematics
  void fillKinematics(float pt, float eta);
  //! event by event final dataset fill
  void fillAll(int promptDecay, float ftracker, float fhcal, float fecalJ, float fecalGT, float fcombined, float dzvtx, float dxyvtx, float dxyErrvtx);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:

  int myPromptDecay;
  float myPt, myEta;
  float myTrackerIsol; 
  float myHcalIsol;
  float myEcalJurIsol, myEcalGTIsol;
  float myCombinedIsol;
  float myDzVtx;
  float myDxyVtx;
  float myDxyErrVtx;
  double myWeight;
  double myProcesId;
  float myLumi;

  TFile* myFile;
  TTree* myTree;
};

#endif // RedIsolationVtxOptimTree_h
