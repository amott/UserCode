#ifndef VecbosEGObject_h
#define VecbosEGObject_h

#include "VecbosBase.hh"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "../src/HggPhysUtils.cc"

#define OLD

//object for a basic cluster
class VecbosBC{
public:
  VecbosBC(VecbosBase*, int);
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
  float etaTilt;
  float phiTilt;
  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;
};

VecbosBC::VecbosBC(VecbosBase* o, int i){
  if(i>o->nBC) return;
  energy  = o->energyBC[i];
  eta     = o->etaBC[i];
  phi     = o->phiBC[i];
  e3x3    = o->e3x3BC[i];
  e5x5    = o->e5x5BC[i];
  eTop    = 0; //temporary b/c Vecbos doesn't have these yet
  eLeft   = 0;
  eRight  = 0;
  eBottom = 0;
  eMax    = o->eMaxBC[i];
  e2nd    = o->e2ndBC[i];
  
#ifndef OLD
  etaCrystal = o->etaCrystalBC[i];
  phiCrystal = o->phiCrystalBC[i];
  iEta = o->iEtaBC[i];
  iPhi = o->iPhiBC[i];
  thetaTilt = o->thetaTiltBC[i];
  phiTilt   = o->phiTiltBC[i];
#endif
  sigmaIEtaIEta = o->covIEtaIEtaBC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiBC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiBC[i];
};


//super cluster
class VecbosSC{
public:
  VecbosSC(VecbosBase*, int);
  std::vector<VecbosBC> basicClusters; // basic clusters associated with this supercluser
  float energy;
  float esEnergy; // in the trees, but not yet in VecbosBase
 float eta;
  float phi;
  float e3x3;
  float e5x5;
  float sigmaIEtaIEta;
  float sigmaIEtaIPhi;
  float sigmaIPhiIPhi;

  float rawE;
  float phiWidth;
  float etaWidth;
  float HoverE;

  float r9(){return e3x3/rawE;}
};

struct sort_pred {
  bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
    return left.second < right.second;
  }
};

VecbosSC::VecbosSC(VecbosBase* o, int i){
  if(i>o->nSC) return;
  energy   = o->energySC[i];
#ifndef OLD
  esEnergy = o->esEnergySC[i];
#endif
  eta      = o->etaSC[i];
  phi      = o->phiSC[i];
  e3x3     = o->e3x3SC[i];
  e5x5     = o->e5x5SC[i];
  sigmaIEtaIEta = o->covIEtaIEtaSC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiSC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiSC[i];

  rawE     = o->rawEnergySC[i];
  phiWidth = o->phiWidthSC[i];
  etaWidth = o->etaWidthSC[i];
  HoverE   = o->hOverESC[i];

  //get the basic clusters
  std::vector< std::pair<int,float> > indexEnergyMap;
  for(int j=0; j<o->nBC; j++){
    if(o->indexSCBC[j] == i){ // the basic cluster points to this SC
      indexEnergyMap.push_back(std::pair<int,float>(j,o->energyBC[j]) );
    }
  }
  std::sort(indexEnergyMap.begin(),indexEnergyMap.end(),sort_pred()); // sort the BCs by energy
  for(int j=0;j<indexEnergyMap.size();j++){
    basicClusters.push_back(VecbosBC(o,indexEnergyMap.at(j).first));
  }

};

class VecbosConversion{
public:
  VecbosConversion(VecbosBase*,int);
  TVector3 pPair;
  TVector3 pRefittedPair;
  TLorentzVector p4RefittedPair;

  float eOverP;

  TVector3 vtx;
  float vtxChi2;
  float vtxChi2Prob;
  bool vtxIsValid;
  int vtxNTracks;
  float vtxMVA;
};

VecbosConversion::VecbosConversion(VecbosBase *o, int i)
{
  if(i > o->nConv || i<0){ //initialize to invalid values for the MVA
    pPair.SetXYZ( -999., -999.,-999.);
    pRefittedPair.SetXYZ( -999., -999.,-999.);
    vtx.SetXYZ( -999., -999.,-999.);
    eOverP = -999;
    vtxChi2 = 0;
    vtxChi2Prob = 0;
    vtxNTracks  = 0;
    vtxIsValid = 0;
    vtxMVA     = -999;
  }else{
    pPair.SetXYZ(o->pxPairConv[i],o->pyPairConv[i],o->pzPairConv[i]);
    pRefittedPair.SetXYZ(o->pxRefittedPairConv[i],o->pyRefittedPairConv[i],o->pzRefittedPairConv[i]);
    vtx.SetXYZ(o->xVtxConv[i],o->yVtxConv[i],o->zVtxConv[i]);

    p4RefittedPair.SetPtEtaPhiE(o->ptRefittedPairConv[i], 
				o->etaRefittedPairConv[i], 
				o->phiRefittedPairConv[i], 
				o->energyRefittedPairConv[i]);
    eOverP = o->eOverPRefittedPairConv[i];
    
    vtxChi2 = o->chi2VtxConv[i];
    vtxChi2Prob = o->chi2ProbVtxConv[i];
    vtxIsValid  = o->isValidVtxConv[i];
    vtxNTracks  = o->nTracksVtxConv[i];
    vtxMVA      = o->mvaOutVtxConv[i];
  }
}


class VecbosPho{
public:
  VecbosPho(VecbosBase*,int);
  ~VecbosPho();
  void matchConversion(VecbosBase*,bool);
  float energy;
  float eta;
  float phi;

  float correctedEnergy;
  float correctedEnergyError;

  float scaledEnergy;
  
  float dEoE;
  VecbosSC SC;
  float HoverE;

  int hasPixel;
  TVector3 CaloPos;

  VecbosConversion *conversion;

  bool isBarrel(){return (fabs(this->SC.eta) < 1.48);}
};

VecbosPho::VecbosPho(VecbosBase* o, int i):
  SC(o,o->superClusterIndexPho[i]),
  CaloPos(0.,0.,0.)
{
  if(i>o->nPho) return;
  energy = o->energyPho[i];
  eta    = o->etaPho[i];
  phi    = o->phiPho[i];

  HoverE = o->hOverEPho[i];
  hasPixel = o->hasPixelSeedPho[i];
};

VecbosPho::~VecbosPho(){
  delete conversion;
}

void VecbosPho::matchConversion(VecbosBase *o,bool dR){ //for some reason, the h2gglobe code doesn't always use dR
  float minDR   = 999; int iMinDR   = -1;
  float minDeta = 999; 
  float minDphi = 999; int iMinDetaDphi = -1;
  for(int i=0; i<o->nConv; i++){
    if(o->ptRefittedPairConv[i] < 1 ||
       !o->isValidVtxConv[i] ||
       o->nTracksVtxConv[i] < 2 ||
       o->chi2ProbVtxConv[i]< 1e-6) continue;  //skip bad conversions

    float convPhi = o->phiRefittedPairConv[i];
    //convert the eta based on the z or the PV
    float convEta = etaTransformation(o->etaRefittedPairConv[i], o->zOfPVFromTracksConv[i] ); 
    
    float Deta = fabs(convEta - this->eta);
    float Dphi = DeltaPhi(convPhi,this->phi);
    float DR   = DeltaR(convEta,this->eta,convPhi,this->phi);

    if(DR < minDR){
      minDR = DR;
      iMinDR = i;
    }
    if( Deta < minDeta && Dphi < minDphi){
      minDphi = Dphi; minDeta = Deta;
      iMinDetaDphi = i;
    }
  }
  int matchIndex = -1;
  if(dR){
    if(minDR < 0.1) matchIndex = iMinDR;
  }else{
    if(minDeta<0.1 && minDphi < 0.1) matchIndex = iMinDetaDphi;
  }
  if(matchIndex == -1) conversion = new VecbosConversion(o,-1);
  else conversion = new VecbosConversion(o,matchIndex);
}


class VecbosEle{
public:
  VecbosEle(VecbosBase*,int);
  float energy;
  float eta;
  float phi;

  VecbosSC SC;
  float esEnergy;
  float HoverE;
  bool isEcalDriven;
  bool isTrackerDriven;

  TVector3 CaloPos;
};

VecbosEle::VecbosEle(VecbosBase* o,int i):
  SC(o,o->superClusterIndexPho[i])
{
  if(i>o->nEle) return;
  energy = o->energyEle[i];
  eta    = o->etaEle[i];
  phi    = o->phiEle[i];
#ifndef OLD  
  esEnergy = o->esEnergyEle[i];
#endif
  HoverE   = o->hOverEPho[i];

  int tmp = o->recoFlagsEle[i];
  //extract the flag information into booleans
  isTrackerDriven = tmp & 1; 
  isEcalDriven    = (tmp >> 1) & 1;
};

#endif
