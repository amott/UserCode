#ifndef VecbosEGObject_h
#define VecbosEGObject_h

#include "VecbosBase.hh"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "../src/HggPhysUtils.cc"

//object for a basic cluster
class VecbosBC{
public:
  VecbosBC();
  VecbosBC(VecbosBase*, int);
  void Init(VecbosBase*,int);
  int index;
  float energy;
  float eta;
  float phi;
  float e3x3;
  float e5x5;
  float eTop;
  float eLeft;
  float eRight;
  float eBottom;
  float eMax;
  float e2nd;
  float etaCrystal;
  float phiCrystal;
  int iEta;
  int iPhi;
  float thetaTilt;
  float phiTilt;
  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;

  //BCInfo getStruct();
};

class VecbosPFBC{
public:
  VecbosPFBC();
  VecbosPFBC(VecbosBase*, int);
  void Init(VecbosBase*,int);
  int index;
  float energy;
  float eta;
  float phi;
  float e3x3;
  float e5x5;
  float eTop;
  float eLeft;
  float eRight;
  float eBottom;
  float eMax;
  float e2nd;
  float etaCrystal;
  float phiCrystal;
  int iEta;
  int iPhi;
  float thetaTilt;
  float phiTilt;
  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;

  //BCInfo getStruct();
};

//super cluster
class VecbosSC{
public:
  VecbosSC();
  VecbosSC(VecbosBase*, int);
  virtual void Init(VecbosBase*, int);
  std::vector<VecbosBC> basicClusters; // basic clusters associated with this supercluser
  int index;
  float energy;
  float esEnergy; // in the trees, but not yet in VecbosBase
  float eta;
  float phi;
  float e3x3;
  float e5x5;
  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;
  
  TVector3 CaloPos;

  float rawE;
  float phiWidth;
  float etaWidth;
  float HoverE;

  //for Hgg corrector
  float r9Scale;

  float r9(){return e3x3/rawE*r9Scale;}
  //  SCInfo getStruct();
};

class VecbosPFSC : public VecbosSC{
public:
  VecbosPFSC();
  VecbosPFSC(VecbosBase*, float, float); // for eta/phi match to the SC position
  VecbosPFSC(VecbosBase*, int);
  void Init(VecbosBase*, int);
  void Init(VecbosBase*, float, float);
  std::vector<VecbosPFBC> pfClusters;
};

class VecbosConversion{
public:
  VecbosConversion();
  VecbosConversion(VecbosBase*,int);
  void Init(VecbosBase*, int);
  int index;
  TVector3 pPair;
  TVector3 pRefittedPair;
  TLorentzVector p4RefittedPair;

  TVector3 CaloPos;

  float eOverP;

  TVector3 vtx;
  float vtxChi2;
  float vtxChi2Prob;
  bool vtxIsValid;
  int vtxNTracks;
  float vtxMVA;

  //track info
  float trk1Dz;       
  float trk1DzError;  
  float trk1Charge;   
  float trk1Algo;     
  float trk1D0;       
  float trk1Pout;
  float trk1Pin;

  float trk2Dz;       
  float trk2DzError;  
  float trk2Charge;   
  float trk2Algo;     
  float trk2D0;       
  float trk2Pout;     
  float trk2Pin; 

  //ConversionInfo getStruct();
};

class VecbosPho{
public:
  VecbosPho();
  VecbosPho(VecbosBase*,int);
  void Init(VecbosBase*, int);
  void matchConversion(VecbosBase*,bool);
  int index;
  float energy;
  float eta;
  float phi;


  //Hgg Correction Variables
  float correctedEnergy;
  float correctedEnergyError;
  float scaledEnergy;
  float scaledEnergyError;
  //float smearedEnergy;
  //float smearedEnergyError;
  float finalEnergy;
  float finalEnergyError;

  float dEoE;
  float dEoEErr;
  VecbosSC SC;
  VecbosPFSC PFSC;
  float HoverE;

  int hasPixel;
  TVector3 CaloPos;

  std::vector<float> photonTrkIsoFromVtx;
  pair<float,int> photonWorstIsoDR03;
  pair<float,int> photonWorstIsoDR04;


  TLorentzVector p4FromVtx(TVector3 vtx,float E,bool pf=false);
  VecbosConversion conversion;
  //isolation variables
  float dr03EcalRecHitSumEtCone;
  float dr03HcalTowerSumEtCone;
  float dr03TrkSumPtCone;
  float dr03TrkSumPtHollowCone;
  float dr04EcalRecHitSumEtCone;
  float dr04HcalTowerSumEtCone;
  float dr04TrkSumPtCone;
  float dr04TrkSumPtHollowCone;

  bool isBarrel(){return (fabs(this->SC.eta) < 1.48);}
  int  getCategory(){ (SC.r9()>0.94)+2*(isBarrel()); } //get the category 0-3 of the photon
  //PhoInfo getStruct();
};
class VecbosEle{
public:
  VecbosEle();
  VecbosEle(VecbosBase*,int);
  void Init(VecbosBase*, int);
  int index;
  float energy;
  float eta;
  float phi;

  VecbosSC SC;
  float esEnergy;
  float HoverE;
  bool isEcalDriven;
  bool isTrackerDriven;

  TVector3 CaloPos;
  
  //EleInfo getStruct();
};

class VecbosMu{
public:
  VecbosMu();
  VecbosMu(VecbosBase*, int);
  void Init(VecbosBase*, int);
  int index;
  float energy;
  float pt;
  float eta;
  float phi;
  TLorentzVector p4;
  int charge;
  float combinedIso;
  float emIso;
  float hadIso;
  float trkIso;
  bool isGlobalMuon;
  bool isTrackerMuon;
  bool isPromptMuon;
  int nTrackHits;
  int nPixelHits;
  float trackImpactPar;

  bool isLooseMuon;
  bool isTightMuon;
};

typedef std::vector<VecbosPho> PhoCollection;
typedef std::vector<VecbosSC> SCCollection;
typedef std::vector<VecbosPFSC> PFSCCollection;
typedef std::vector<VecbosBC> BCCollection;
typedef std::vector<VecbosConversion> ConvCollection;
typedef std::vector<VecbosMu> MuCollection;
#endif