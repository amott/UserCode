 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef CaltechReducer_h
#define CaltechReducer_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

#include "VecbosEGObject.hh"
#include "VecbosJet.hh"
#include "VecbosMu.hh"
#include "VecbosVtx.hh"
#include "VecbosMET.hh"


#include "HggVertexingNew.hh"
#include "HggEGEnergyCorrector.hh"
#include "HggEnergyScale.hh"

class CaltechReducer : public Vecbos{
public:
  CaltechReducer(TTree *tree=0); /// Class Constructor
  CaltechReducer(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false,int mod=-1); /// Class Constructor

  void setConfig(string s){config = s;}

  void Loop(string outFileName, int start, int stop);

  void addTrigger(string s){triggerNames.push_back(s);}
  void setDoVertexingVMA(bool b){doVertexingMVA=b;}

private:
  TTree * outTree;

  //private member functions:
  void init(string); // do variable initialization 
  void clearAll();
  void setOutputBranches();

  //methods for filling objects
  void fillEventInfo();
  void fillTriggerInfo();
  void fillMuons();
  void fillElectrons();
  void fillPhotons();
  void fillJets();
  void fillMets();
  void fillVertexInfo();
  void fillGeneratorInfo();

  void setupCollections(ReadConfig&);
  std::vector<std::string> validCollections; 
  std::map<std::string,bool> collectionSwitches;
  std::map<std::string,float> ptThresholds;
  std::map<std::string,int>   reqObjects;

  //configuration files for the various steps
  std::string config;

  bool _goodRunLS; 
  bool _isData;

  //configuration file and MVA parameters
  
  //for the Hgg-like vertexing MVA
  bool doVertexingMVA;
  void doVertexing();
  HggVertexing *vertexer;
  //energy correction variables
  std::string correctionType;
  HggEGEnergyCorrector *phocorrector;
  HggEGEnergyCorrector *elecorrector;

  // define variables for the output tree:
  std::vector<std::string> triggerNames; // list of all the triggers to consider
  std::vector<std::vector<std::string> > masks;
  int * triggerBits;       // this will be an array of the trigger decision per event (for the output tree)

  // ...
  //Event info
  int lumiBlock_; 
  int runNumber_; 
  int evtNumber_; 
  int bunchX_; 
  int orbitNumber_; 
  int evtTime_; 

  int phyDeclared_; 
  float rho_; 
  float rhoEtaMax44_; 

  //main object collections for the reduced tree 
  int nPho_;
  PhoCollection  Photons_;
  int nMu_;
  MuCollection Muons_;
  int nEle_;
  EleCollection Electrons_;
  int nJet_;
  JetCollection Jets_;
  int nCaloMet_;
  CaloMETCollection CaloMet_;
  int nPFMet_;
  PFMETCollection   PFMet_;
  int nTCMet_;
  TCMETCollection   TCMet_;
  int nVtx_;
  VtxCollection Vertices_;
  //generator collections
  //gen-leve phton
  int nGenPho_;
  GenCollection GenPhotons_;
  int nGenMu_;
  GenCollection GenMuons_;
  int nGenEle_;
  GenCollection GenElectrons_;
  int nGenHiggs_;
  GenCollection GenHiggs_;
  int nGenOthers_;
  GenCollection GenOthers_;



  //this is the collection of the TMVA selected vertices for each photon pair
  //the format is std::pair< (iPho1 << 14) + iPho2,iVrt> 
  int nPair_;
  vector<pair<int,int> > ggVerticesPhotonIndices; 
  vector<pair<int,float> > ggVerticesVertexIndex;  // vertex, MVA score
  vector<float>            ggVerticesPerEvtMVA;
  //std::vector<int> ggVertices_

  int procID;
  float qScale;
  float nPu;

  float caloMet;
  float caloMetPhi;

  float pfMet;
  float pfMetPhi;

  float tcMet;
  float tcMetPhi;
};
#endif


