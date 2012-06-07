#include <VecbosEGObject.hh>
#include <HggMassResolution.hh>
#include <HggPhotonID.hh>

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
  std::pair<int,int> getBestPairPFCiC(int,int);
  
  bool doElectronVeto;
  bool isData_;
  
  string massResConfig;
  HggMassResolution *massRes;

  HggPhotonID *PhotonID;

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
  void fillGenInfo();

  void fillMuMuGamma();
  //TMVA stuff
  string weightFile_diPho;
  string methodName_diPho;
  TMVA::Reader *diPhotonMVA;

  float getDiPhoMVA(int,int,float,float,bool);

  float getMPair(int,int);
  //These are the variables that will be filled for the TMVA:
  //better to do it this way so we can do PFSC or regular SC without major changes

  
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

  float mPairPFCiC_;
  float mPairNoCorrPFCiC_;
  float mPairResPFCiC_;
  float mPairResWrongVtxPFCiC_;
  int diPhoVtxPFCiC_;
  float diPhoVtxXPFCiC_;
  float diPhoVtxYPFCiC_;
  float diPhoVtxZPFCiC_;

  std::vector<VecbosPho> OutPhotons_;
  std::vector<VecbosPho> OutPhotonsPFCiC_;

  float energyPho1;
  float energyNoCorrPho1;
  float etaPho1;
  float pxPho1;
  float pyPho1;
  float pzPho1;
  float r9Pho1;
  int catPho1;
  int passPFCiCPho1;
  float indexPho1;

  float energyPho2;
  float energyNoCorrPho2;
  float etaPho2;
  float pxPho2;
  float pyPho2;
  float pzPho2;
  float r9Pho2;
  int catPho2;
  int passPFCiCPho2;
  float indexPho2;

  VecbosPho pho1_;
  VecbosPho pho2_;

  int nSigma;
  std::vector<float> mPairScale;
  std::vector<float> pho1MVAScale;
  std::vector<float> pho2MVAScale;
  std::vector<float> diPhoMVAScale;
  std::vector<float> mPairSmear;
  std::vector<float> pho1MVASmear;
  std::vector<float> pho2MVASmear;
  std::vector<float> diPhoMVASmear;

  std::vector<float> mPairScalePFCiC;
  std::vector<float> mPairSmearPFCiC;

  //for mumuG
  const static int maxMuMuG = 500;
  int nMuMuG;
  float massMuMuGamma[maxMuMuG];
  float massMuMuRegGamma[maxMuMuG];
  float massMuMuScaleGamma[maxMuMuG];
  float massMuMuGenGamma[maxMuMuG];
  float massMuMu[maxMuMuG];
  float puWeight[maxMuMuG];
  
  MuCollection MMG_Mu1;
  MuCollection MMG_Mu2;
  PhoCollection MMG_Pho;
  float mvaPho[maxMuMuG];
  float isosumoetPho[maxMuMuG];
  
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

  int nGenHiggs;
  std::vector<VecbosGen> *GenHiggs;

  int nGenPho;
  std::vector<VecbosGen> *GenPhotons;

  float inPU;
};
