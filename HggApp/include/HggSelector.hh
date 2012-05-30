

#include <VecbosEGObject.hh>
#include <HggMassResolution.hh>

#include <vector>
#include <iostream>
#include <string>

#include <TChain.h>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
using namespace std;

                                    
#include "HggVertexing.hh"

class HggSelector{
public:
  HggSelector();
  ~HggSelector();
  HggSelector(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  void setConfig(string s){configFile = s;}
  void setMassResConfig(string s){massResConfig = s;}
  bool isValid(){return valid;}
  void addTrigger(string s){triggers.push_back(s);}
  void setDoMuMuGamma(bool d=true){doMuMuGamma=d;}
  void setNSigma(int ns){nSigma=ns;}

  void suppressElectronVeto(){doElectronVeto=false;}
  void Loop();
private:
  bool valid;
  bool doMuMuGamma;
  TChain* fChain;
  TTree* outTree;
  TTree* outTreeMuMuG;
  string configFile;
  string outputFile;

  std::pair<int,int> getBestPair(float*,int,int);
  
  bool doElectronVeto;
  bool isData_;
  
  string massResConfig;
  HggMassResolution *massRes;

  const static float rhoFac    = 0.17;
  const static float rhoFacBad = 0.52;
  const static float isoSumConst = 0;//5;
  const static float isoSumConstBad = 0;//7;

  vector<string> triggers;
  int *triggerDec;
  bool requireTrigger();

  int getVertexIndex(int,int);
  float getVertexMVA(int,int);

  int init();
  void clear();
  void setDefaults();
  void setBranchAddresses();
  void setupOutputTree();
  void setupTMVA();

  void fillMuMuGamma();
  //TMVA stuff
  string weightFile_IdEB;
  string weightFile_IdEE;
  string weightFile_diPho;
  string methodName_Id;
  string methodName_diPho;
  TMVA::Reader *photonMVA_EB;
  TMVA::Reader *photonMVA_EE;
  TMVA::Reader *diPhotonMVA;
  float getIdMVA(int,int,bool,bool,int);
  float getDiPhoMVA(int,int,float,float,bool);

  float getMPair(int,int);
  //These are the variables that will be filled for the TMVA:
  //better to do it this way so we can do PFSC or regular SC without major changes

  
  //photon ID MVA:
  float hoe;
  float sigietaieta;
  float isosumoet;
  float trkisooet;
  float isosumoetbad;
  float r9;
  float ecalisodr03;
  float hcalisodr04;
  float nVertexf;
  float etasc;
  float scetawidth;
  float scphiwidth;
  //diPhotonMVA:
  float smearedMassErrByMass; // mass error smeared/mass
  float smearedMassErrByMassWrongVtx;
  float vtxprob;
  float pho1PtByMass;
  float pho2PtByMass;
  float pho1Eta;
  float pho2Eta;
  float cosDPhi;
  float pho1IdMVA;
  float pho2IdMVA;
  //-------------

  //Output Variables:
  int trigger_;
  float mPair_;
  float mPairNoCorr_;
  float mPairRes_;
  float mPairResWrongVtx_;
  float diPhoMVA_;
  float pho1MVA_;
  float pho2MVA_;
  int diPhoVtx_;
  float diPhoVtxX_;
  float diPhoVtxY_;
  float diPhoVtxZ_;
  float pfmPair_;
  float pfmPairRes_;
  float pfmPairResWrongVtx_;
  float pfdiPhoMVA_;
  float pfpho1MVA_;
  float pfpho2MVA_;
  int pfdiPhoVtx_;

  float energyPho1;
  float energyNoCorrPho1;
  float etaPho1;
  float pxPho1;
  float pyPho1;
  float pzPho1;
  float r9Pho1;
  float indexPho1;

  float energyPho2;
  float energyNoCorrPho2;
  float etaPho2;
  float pxPho2;
  float pyPho2;
  float pzPho2;
  float r9Pho2;
  float indexPho2;

  float energyPFPho1;
  float etaPFPho1;
  float pxPFPho1;
  float pyPFPho1;
  float pzPFPho1;
  float r9PFPho1;
  float indexPFPho1;

  float energyPFPho2;
  float etaPFPho2;
  float pxPFPho2;
  float pyPFPho2;
  float pzPFPho2;
  float r9PFPho2;
  float indexPFPho2;

  VecbosPho pho1_;
  VecbosPho pho2_;
  VecbosPho pfpho1_;
  VecbosPho pfpho2_;

  int nSigma;
  std::vector<float> mPairScale;
  std::vector<float> pho1MVAScale;
  std::vector<float> pho2MVAScale;
  std::vector<float> diPhoMVAScale;
  std::vector<float> mPairSmear;
  std::vector<float> pho1MVASmear;
  std::vector<float> pho2MVASmear;
  std::vector<float> diPhoMVASmear;

  //for mumuG
  const static int maxMuMuG = 500;
  int nMuMuG;
  float massMuMuGamma[maxMuMuG];
  float massMuMuRegGamma[maxMuMuG];
  float massMuMuScaleGamma[maxMuMuG];
  float massMuMuGenGamma[maxMuMuG];
  float massMuMu[maxMuMuG];
  float puWeight[maxMuMuG];
  
  float ptMu1[maxMuMuG];
  float etaMu1[maxMuMuG];
  float phiMu1[maxMuMuG];
  int chargeMu1[maxMuMuG];
  float isoMu1[maxMuMuG];
  float drPhoMu1[maxMuMuG];
  //gen muon, same charge, stat 1
  //dR < 0.5, ptRelDif < 0.5
  //two gen objects cannot match reco object
  int genMatchMu1[maxMuMuG];
  float genPtMu1[maxMuMuG];
  int genIdMu1[maxMuMuG];
  int genStatusMu1[maxMuMuG];
  int genIdMomMu1[maxMuMuG];
  int genStatusMomMu1[maxMuMuG];

  float ptMu2[maxMuMuG];
  float etaMu2[maxMuMuG];
  float phiMu2[maxMuMuG];
  int chargeMu2[maxMuMuG];
  float isoMu2[maxMuMuG];
  float drPhoMu2[maxMuMuG];
  int genMatchMu2[maxMuMuG];
  float genPtMu2[maxMuMuG];
  int genIdMu2[maxMuMuG];
  int genStatusMu2[maxMuMuG];
  int genIdMomMu2[maxMuMuG];
  int genStatusMomMu2[maxMuMuG];
  
  float defEnergyPho[maxMuMuG];
  float regEnergyPho[maxMuMuG];
  float scaleEnergyPho[maxMuMuG];
  float etaPho[maxMuMuG];
  float phiPho[maxMuMuG];
  int   eleMatchPho[maxMuMuG];
  float r9Pho[maxMuMuG];
  float hOePho[maxMuMuG];
  float dr03EcalIsoPho[maxMuMuG];
  float dr04HcalIsoPho[maxMuMuG];
  float isosumoetPho[maxMuMuG];
  float dr03TrkSumHollowConePho[maxMuMuG];  
  float mvaPho[maxMuMuG];
  // photon matching: dR < 0.2
  // |recoPt-truePt/truePt| < 1

  int   genMatchPho[maxMuMuG];
  float genEnergyPho[maxMuMuG];
  int   genIdPho[maxMuMuG];
  int   genStatusPho[maxMuMuG];
  int   genIdMomPho[maxMuMuG];
  int   genStatusMomPho[maxMuMuG];
  
  
  int getGenMatchPho(VecbosPho*);
  int getGenMatchMu(VecbosMu*);

  //----
  float genHiggsPt;
  float genHiggsVx;
  float genHiggsVy;
  float genHiggsVz;

  float ptGenPho1;
  float etaGenPho1;
  float phiGenPho1;
  float energyGenPho1;
  float ptGenPho2;
  float etaGenPho2;
  float phiGenPho2;
  float energyGenPho2;
  
  float nPU_;

  //member data
  const static int maxPho = 100;
  int nPho_;
  std::vector<VecbosPho> *Photons_; // this contains ALL photons
  
  bool photonMatchedElectron[maxPho];

  //this is the collection of the TMVA selected vertices for each photon pair
  //the format is std::pair< (iPho1 << 14) + iPho2,iVrt> 
  int nPair_;
  std::vector<std::pair<int,int> > *ggVerticesPhotonIndices;
  std::vector<std::pair<int, float> > *ggVerticesVertexIndex;

  int nMu_;
  std::vector<VecbosMu> *Muons_;

  
  // for each photon, a vector of floats giving the track iso from each ggVertex

  int nVtx; 
  static const int MAXVX = 100;
  float vtxX[MAXVX];
  float vtxY[MAXVX];
  float vtxZ[MAXVX];
  float vtxChi2[MAXVX];
  float vtxNdof[MAXVX];
  float vtxNormalizedChi2[MAXVX];
  int vtxTrackSize[MAXVX];
  int vtxIsFake[MAXVX];
  int vtxIsValid[MAXVX];

  float rho;
  float rhoEtaMax44;

  int lumiBlock;
  int runNumber;
  long evtNumber;
  bool _isData;

  float higgsPtIn;
  float higgsVxIn;
  float higgsVyIn;
  float higgsVzIn;
  static const int MAXGenSaved=1000;
  int nGenPho;
  float etaGenPho[MAXGenSaved];      
  float phiGenPho[MAXGenSaved];      
  float ptGenPho[MAXGenSaved];       
  float energyGenPho[MAXGenSaved];   
  int pidMomGenPho[MAXGenSaved];

  int nGenMu;
  float etaGenMu[MAXGenSaved];      
  float phiGenMu[MAXGenSaved];      
  float ptGenMu[MAXGenSaved];       
  float energyGenMu[MAXGenSaved];   
  int pidMomGenMu[MAXGenSaved];

  float inPU;
};
