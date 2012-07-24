#ifndef RedVecbosMuTree_h
#define RedVecbosMuTree_h
class TFile;
class TTree;

class G3EventProxy;

class RedVecbosMuTree {
public:
   RedVecbosMuTree(const char * filename = "vecbosTree.root");
  ~RedVecbosMuTree();

  //! add the CSA07 processID and weight block
  void addCSA07Infos();
  //! event by event final dataset fill
  void fillAll(int nj, float eta1, float phi1, float pt1, float eta2, float phi2, float pt2, float helic, float imass, float tmass, float met);
  void fillAll(int nj, float eta1, float phi1, float pt1, float eta2, float phi2, float pt2, float helic, float imass, float tmass, float met,float metx,float mety);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  //void saveEfficiencies(std::vector<float>);

  //! save in the ROOT file
  void save();
  //reset all vars
  void newEvent();

  ///Generic
  float weight;

  ///Jets Var
  int   myNJets; 
  int   allNJets;  
  float energyJet[20];
  float etJet[20]; 
  float etaJet[20]; 
  float phiJet[20]; 
  float vertexXJet[20];
  float vertexYJet[20];
  float vertexZJet[20];
  float emFracJet[20];
  float hadFracJet[20];
  
  ///MuonVar
  int   nMuons;
  int   nGoodMuons;
  float  ZMassForVeto;  

  float bestMuEta[2];  
  float bestMuPhi[2];  
  float bestMuPt[2];  
  float bestMuvtxX[2];
  float bestMuvtxY[2];
  float bestMuvtxZ[2];
      
  float bestMucharge[2]; 
  float bestMuPx[2];   
  float bestMuPy[2];  
  float muTrackDxy[2];         
  float muTrackD0[2];	       
  float muTrackDsz[2];         
  float muTrackDz[2];	       
  float muTrackDxyError[2];  
  float muTrackD0Error[2];   
  float muTrackDszError[2];  
  float muTrackDzError[2];   
  float sumPt03[2];    
  float emEt03[2];     
  float hadEt03[2];    
  float hoEt03[2];     
  float nTrk03[2];     
  float nJets03[2];    
  float sumPt05[2];    
  float emEt05[2];     
  float hadEt05[2];    
  float hoEt05[2];     
  float nTrk05[2];     
  float nJets05[2];    
  float EcalExpDepo[2];        
  float HcalExpDepo[2];        
  float HoExpDepo[2];	       
  float emS9[2];       
  float hadS9[2];      
  float hoS9[2];       
  float CaloComp[2];	       
  float SumMuonsPx;
  float SumMuonsPy;
   
  ///MET var 
  float myMet;
  float myMetx;  
  float myMety;
  int nZmumuCand;

  ///VecBos Var
  float myTransvMass;
  float Z0Mass;
  float Z0Pt  ;
  float Z0Eta ;
  float Z0Phi ;

  ///Misc
  float myHelicity;
  double myWeight;
  double myProcesId;
  float  myLumi;
   int      nPVTX;
   float    PVTXxPV[10];   //[nPV]
   float    PVTXyPV[10];   //[nPV]
   float    PVTXzPV[10];   //[nPV]
   float    PVTXErrxPV[10];   //[nPV]
   float    PVTXErryPV[10];   //[nPV]
   float    PVTXErrzPV[10];   //[nPV]
   float    SumPtPVTX[10];   //[nPV]
   float    ndofPVTX[10];   //[nPV]
   float    chi2PVTX[10];   //[nPV]


private:

  TFile* myFile;
  TTree* myTree;

};

#endif // RedVecbosMuTree_h
