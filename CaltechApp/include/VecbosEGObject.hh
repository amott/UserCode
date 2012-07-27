#ifndef VecbosEGObject_h
#define VecbosEGObject_h

#include "VecbosBase.hh"
#include "VecbosGen.hh"
#include "VecbosBaseObject.hh"
#include "VecbosPhysicsObject.hh"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "../src/HggPhysUtils.cc"

//object for a basic cluster
class VecbosBC : public VecbosBaseObject{
public:
  VecbosBC();
  VecbosBC(VecbosBase*, int);
  virtual void Init(VecbosBase*,int);
  float e3x3;
  float e5x5;
  float eTop;
  float eLeft;
  float eRight;
  float eBottom;
  float eMax;
  float e2nd;

  float e2x5Max;
  float e2x5Left;
  float e2x5Right;
  float e2x5Top;
  float e2x5Bottom;

  float etaCrystal;
  float phiCrystal;
  int iEta;
  int iPhi;
  float thetaTilt;
  float phiTilt;
  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;
};

//super cluster
class VecbosSC : public VecbosBaseObject{
public:
  VecbosSC();
  VecbosSC(VecbosBase*, int);
  virtual void Init(VecbosBase*, int);
  VecbosBC BCSeed;
  int nBCs;
  std::vector<VecbosBC> basicClusters; // basic clusters associated with this supercluser
  float esEnergy;
  float e3x3;
  float e5x5;

  float e3x1;
  float e1x3;
  float e4x4;
  float eMax;
  float e2x2;
  float e2nd;
  float e1x5;
  float e2x5Max;
  float e2x5Left;
  float e2x5Right;
  float e2x5Top;
  float e2x5Bottom;
  
  float eLeft;
  float eRight;
  float eTop;
  float eBottom;

  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;

  float esEffSigRR;
  
  TVector3 CaloPos;

  float rawE;
  float phiWidth;
  float etaWidth;
  float HoverE;

  float r9;
};

class VecbosConversion : public VecbosBaseObject{
public:
  VecbosConversion();
  VecbosConversion(VecbosBase*,int);
  void Init(VecbosBase*, int);
  TVector3 pPair;
  TVector3 pRefittedPair;
  TLorentzVector p4RefittedPair;

  TVector3 CaloPos;

  float eOverP;

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


class VecbosPho : public VecbosPhysicsObject{
public:
  VecbosPho();
  VecbosPho(VecbosBase*,int);
  void Init(VecbosBase*, int);
  void matchConversion(VecbosBase*,bool);

  //Hgg Correction Variables
  float regressionEnergy;
  float regressionEnergyError;
  float scaledEnergy;
  float scaledEnergyError;
  //float smearedEnergy;
  //float smearedEnergyError;
  float finalEnergy;
  float finalEnergyError;

  float dEoE;
  float dEoEErr;

  VecbosSC SC;
  VecbosConversion conversion;

  float HoverE;
  float HTowOverE;

  int hasPixel;
  TVector3 CaloPos;

  std::vector<float> photonTrkIsoFromVtx;
  pair<float,int> photonWorstIsoDR03;
  pair<float,int> photonWorstIsoDR04;

  TLorentzVector p4FromVtx(TVector3 vtx,float E);
  //isolation variables
  float dr03EcalRecHitSumEtCone;
  float dr03HcalTowerSumEtCone;
  float dr03TrkSumPtCone;
  float dr03TrkSumPtHollowCone;
  float dr04EcalRecHitSumEtCone;
  float dr04HcalTowerSumEtCone;
  float dr04TrkSumPtCone;
  float dr04TrkSumPtHollowCone;

  // these are the special versions used in the Hgg analysis -- 1 for each vertex
  float computeTrackIso(VecbosBase*,int,float,float,float,float,float,float);
  float dr03HggTrackIso[100];
  float dr04HggTrackIso[100];

  //pfIsolation
  int           nPV;
  float         dr01ChargedHadronPFIso[100];
  float         dr02ChargedHadronPFIso[100];
  float         dr03ChargedHadronPFIso[100];
  float         dr04ChargedHadronPFIso[100];
  float         dr05ChargedHadronPFIso[100];
  float         dr06ChargedHadronPFIso[100];

  float                      dr01NeutralHadronPFIso;
  float                      dr02NeutralHadronPFIso;
  float                      dr03NeutralHadronPFIso;
  float                      dr04NeutralHadronPFIso;
  float                      dr05NeutralHadronPFIso;
  float                      dr06NeutralHadronPFIso;
  float                      dr01PhotonPFIso;
  float                      dr02PhotonPFIso;
  float                      dr03PhotonPFIso;
  float                      dr04PhotonPFIso;
  float                      dr05PhotonPFIso;
  float                      dr06PhotonPFIso;

  bool isBarrel(){return (fabs(this->SC.eta) < 1.48);}
  int  getCategory(){ (SC.r9>0.94)+2*(isBarrel()); } //get the category 0-3 of the photon
  //PhoInfo getStruct();

  void doGenMatch(VecbosBase*,int);

  int eleMatchIndex;
  void matchElectron(VecbosBase* o);
};

class VecbosEle : public VecbosPhysicsObject{
public:
  VecbosEle();
  VecbosEle(VecbosBase*,int);
  void Init(VecbosBase*, int);

  float regressionEnergy;
  float regressionEnergyError;

  VecbosSC SC;
  float esEnergy;
  float HoverE;
  bool isEcalDriven;
  bool isTrackerDriven;
  
  float EOverP;

  float dEtaSCTrack;
  float dPhiSCTrack;

  float dr03ChargedHadronPFIso;
  float dr03NeutralHadronPFIso;
  float dr03PhotonPFIso;

  float dr04ChargedHadronPFIso;
  float dr04NeutralHadronPFIso;
  float dr04PhotonPFIso;

  float dr03TkSumPt;
  float dr03EcalRecHitSumEt;
  float dr03HcalTowerSumEt;

  float dr04TkSumPt;
  float dr04EcalRecHitSumEt;
  float dr04HcalTowerSumEt;
};

typedef std::vector<VecbosPho> PhoCollection;
typedef std::vector<VecbosSC> SCCollection;
typedef std::vector<VecbosBC> BCCollection;
typedef std::vector<VecbosConversion> ConvCollection;
typedef std::vector<VecbosEle> EleCollection;
#endif
