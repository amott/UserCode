 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef HggReducer_h
#define HggReducer_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"



#include "HggVertexing.hh"
#include "HggEGEnergyCorrector.hh"
#include "HggEnergyScale.hh"

using namespace std;

struct PreSelCuts{
  float scet1;
  float scet2;
  float maxeta;
  float ecaliso[2]; //eb,ee
  float trkiso[2]; //eb,ee
  float hcaliso[2]; //eb,ee 
  float sieie[2]; //eb,ee
  float hoe;
};

enum EnergyRegressionMethod{JoshV1, JoshV2,Yong};

class HggReducer : public Vecbos{
public:
  HggReducer(TTree *tree=0); /// Class Constructor
  HggReducer(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false,int mod=-1); /// Class Constructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void addTrigger(string s){triggerNames.push_back(s);}
  void setConfig(string s){config = s;}
  void setCorrectionType(int s){correctionType = s;}
  void setPreselectionSet(int s){preSelSet = s;}
  void setScaleSmear(int s){applyScaleSmear = s;}
  void SetVtxConfigFile(string s){vertexCFG = s;}
  void SetScaleConfigFile(string s){energyScaleCFG = s;}
  void SetSmearConfigFile(string s){energySmearCFG = s;}
  void SetMinPhoSelection(int n){minPhoSel = n;}
private:
  bool _goodRunLS; 
  bool _isData;
  float _weight;
  TTree *_treeCond;

  float minPhoSel;
  //configuration file and MVA parameters
  string config;
  string vertexCFG;
  string energyScaleCFG;
  string energySmearCFG;

  void init(string); // do variable initialization 
  void clearAll();
  void setOutputBranches();
  
  bool passPreselection(VecbosPho *);

  float computeTrackIso(int iPho, int iVtx,
			float ptMinTrkIso,
			float outerConeTrkIso,
			float innerConeTrkIso,
			float etaStripHalfWidthTrkIso,
			float dzMaxTrkIso,
			float dyMaxTrkIso);

  //photon preselection variables
  void setupPreSelection();
  void fillMuons();
  void fillElectrons();
  void fillJets();
  vector<PreSelCuts> preselections;
  int preSelSet;

  HggVertexing *vertexer;
  //energy correction variables
  string correctionType;
  HggEGEnergyCorrector *corrector;
  HggEGEnergyCorrector *elecorrector;

  //energy smearing
  int applyScaleSmear;
  HggEnergyScale *energyScale;
  HggEnergyScale *energySmear;

  TTree * outTree;
  // define variables for the output tree:
  vector<string> triggerNames; // list of all the triggers to consider
  int * triggerBits;       // this will be an array of the trigger decision per event (for the output tree)

  // ...
  //Event info
  int lumiBlockO; 
  int runNumberO; 
  int evtNumberO; 
  int bunchX; 
  int orbitNumber; 
  int evtTime; 

  int phyDeclared; 
  float rho; 
  float rhoEtaMax44; 

  vector<short> *pileupBunchX; 
  vector<short> *pileupNInteraction; 
  float pileupTrueNumInterations;

  //main object collections for the reduced tree 
  const static int maxPho=100;
  int nPho_;
  vector<VecbosPho>  Photons_;
  int nSC_;
  /*
  SCCollection   SuperClusters_; // this contains every SC matched to a photon
  int nPFSC_
  PFSCCollection   PFSuperClusters_; // this contains every PFSC matched to a photon
  int nBC_;
  BCCollection   BasicClusters_; // this contains up to 4 BCs for each above SC
  int nConv_;
  ConvCollection Conversions_;   // this contains all conversions match to a photon
  */  
  void matchPhotonsElectrons();
  bool photonMatchedElectron[maxPho];

  int nMu_;
  MuCollection Muons_;
  int nEle_;
  EleCollection Electrons_;

  //this is the collection of the TMVA selected vertices for each photon pair
  //the format is std::pair< (iPho1 << 14) + iPho2,iVrt> 
  int nPair_;
  vector<pair<int,int> > ggVerticesPhotonIndices; 
  vector<pair<int,float> > ggVerticesVertexIndex01;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndex02;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndex03;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndexOld01;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndexOld02;  // vertex, MVA score
  vector<pair<int,float> > ggVerticesVertexIndexOld03;  // vertex, MVA score
  //std::vector<int> ggVertices_

  //vertex information
  void fillVertexInfo();
  int nVtx; 
  static const int maxVtx = 100;
  float vtxX[maxVtx];
  float vtxY[maxVtx];
  float vtxZ[maxVtx];
  float vtxChi2[maxVtx];
  float vtxNdof[maxVtx];
  float vtxNormalizedChi2[maxVtx];
  int vtxTrackSize[maxVtx];
  int vtxIsFake[maxVtx];
  int vtxIsValid[maxVtx];


  //GENERATOR information
  virtual void fillGeneratorInfo();  // this will be overwritten in the MC class to actually do something
  static const int MAXGenSaved = 1000;
  //gen-leve phton
  int nGenPho;
  GenCollection GenPhotons;

  int nGenMu;
  GenCollection GenMuons;

  int nGenEle;
  GenCollection GenElectrons;

  int nGenHiggs;
  GenCollection GenHiggs;

  int nGenOthers;
  GenCollection GenOthers;

  int procID;
  float qScale;
  float nPu;

  const static int maxJets=50;
  int nJets;
  float ptJet[maxJets];
  float etaJet[maxJets];
  float phiJet[maxJets];
  float energyJet[maxJets];
  
  float caloMet;
  float caloMetPhi;

  float pfMet;
  float pfMetPhi;

  float tcMet;
  float tcMetPhi;
};
#endif


