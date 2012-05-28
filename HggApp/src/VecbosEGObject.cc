#include <VecbosEGObject.hh>
#include "CommonTools/include/Utils.hh"

VecbosBC::VecbosBC(){

}

VecbosBC::VecbosBC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosBC::Init(VecbosBase* o, int i){
  if(i>o->nBC){
    index = -1;
    return;
  }
  index   = i;
  energy  = o->energyBC[i];
  eta     = o->etaBC[i];
  phi     = o->phiBC[i];
  e3x3    = o->e3x3BC[i];
  e5x5    = o->e5x5BC[i];
  eTop    = o->eTopBC[i]; 
  eLeft   = o->eLeftBC[i];
  eRight  = o->eRightBC[i];
  eBottom = o->eBottomBC[i];
  eMax    = o->eMaxBC[i];
  e2nd    = o->e2ndBC[i];
  
  etaCrystal = o->etaCrystalBC[i];
  phiCrystal = o->phiCrystalBC[i];
  iEta = o->iEtaBC[i];
  iPhi = o->iPhiBC[i];
  thetaTilt = o->thetaTiltBC[i];
  phiTilt   = o->phiTiltBC[i];

  sigmaIEtaIEta = o->covIEtaIEtaBC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiBC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiBC[i];
};

/*
BCInfo VecbosBC::getStruct(){
  BCInfo out = {
    index,
    energy,
    eta,
    phi,
    e3x3,
    e5x5,
    eTop,
    eLeft,
    eRight,
    eBottom,
    eMax,
    e2nd,
    etaCrystal,
    phiCrystal,
    iEta,
    iPhi,
    etaTilt,
    phiTilt,
    sigmaIEtaIEta,
    sigmaIEtaIPhi,
    sigmaIPhiIPhi    
  };
  return out;
}
*/
VecbosPFBC::VecbosPFBC(){

}

VecbosPFBC::VecbosPFBC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosPFBC::Init(VecbosBase* o, int i){
  if(i>o->nPFBC) return;
  index   = i;
  energy  = o->energyPFBC[i];
  eta     = o->etaPFBC[i];
  phi     = o->phiPFBC[i];
  e3x3    = o->e3x3PFBC[i];
  e5x5    = o->e5x5PFBC[i];
  eTop    = o->eTopPFBC[i]; 
  eLeft   = o->eLeftPFBC[i];
  eRight  = o->eRightPFBC[i];
  eBottom = o->eBottomPFBC[i];
  eMax    = o->eMaxPFBC[i];
  e2nd    = o->e2ndPFBC[i];
  
  etaCrystal = o->etaCrystalPFBC[i];
  phiCrystal = o->phiCrystalPFBC[i];
  iEta = o->iEtaPFBC[i];
  iPhi = o->iPhiPFBC[i];
  thetaTilt = o->thetaTiltPFBC[i];
  phiTilt   = o->phiTiltPFBC[i];

  sigmaIEtaIEta = o->covIEtaIEtaPFBC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiPFBC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiPFBC[i];
};


struct sort_pred {
  bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
    return left.second < right.second;
  }
};

VecbosSC::VecbosSC(){}

VecbosSC::VecbosSC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosSC::Init(VecbosBase* o, int i){
  if(i>o->nSC){
    index = -1;
    return;
  }    
  index    = i;
  energy   = o->energySC[i];
  esEnergy = o->esEnergySC[i];
  eta      = o->etaSC[i];
  phi      = o->phiSC[i];
  e3x3     = o->e3x3SC[i];
  e5x5     = o->e5x5SC[i];
  sigmaIEtaIEta = o->covIEtaIEtaSC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiSC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiSC[i];

  CaloPos.SetXYZ(o->xPosSC[i],o->yPosSC[i],o->zPosSC[i]);

  rawE     = o->rawEnergySC[i];
  phiWidth = o->phiWidthSC[i];
  etaWidth = o->etaWidthSC[i];
  HoverE   = o->hOverESC[i];
  r9Scale  = 1;
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
  
 }

/*
SCInfo VecbosSC::getStruct(){
  SCInfo out = {
  index,
  0,   // we set the BC pointer seperately
  energy,
  esEnergy, // in the trees, but not yet in VecbosBase
 eta,
  phi,
  e3x3,
  e5x5,
  sigmaIEtaIEta,
  sigmaIEtaIPhi,
  sigmaIPhiIPhi,
  rawE,
  phiWidth,
  etaWidth,
  HoverE,
  this->r9()
  };

  out.BCs[0]=-1.;   out.BCs[1]=-1.;   out.BCs[2]=-1.;   out.BCs[3]=-1.; //initialize
  if(basicClusters.size() >0) out.BCs[0] = basicClusters.at(0).index;
  if(basicClusters.size() >1) out.BCs[1] = basicClusters.at(1).index;
  if(basicClusters.size() >2) out.BCs[2] = basicClusters.at(basicClusters.size()-1).index;
  if(basicClusters.size() >0) out.BCs[0] = basicClusters.at(basicClusters.size()-2).index;

}
*/
VecbosPFSC::VecbosPFSC(){}

VecbosPFSC::VecbosPFSC(VecbosBase* o, int i){
  this->Init(o,i);
}
VecbosPFSC::VecbosPFSC(VecbosBase* o, float etaSC,float phiSC){
  this->Init(o,etaSC,phiSC);
}
void VecbosPFSC::Init(VecbosBase* o, float etaSC,float phiSC){
  const float maxDR = 0.3;

  float DR=9999;  int index = -1;
  for(int i=0;i<o->nPhoPFSC;i++){
    float thisDR = DeltaR(etaSC,o->etaPhoPFSC[i],phiSC,o->phiPhoPFSC[i]);
    if(thisDR < DR && thisDR < maxDR){
      DR = thisDR;
      index = i;      
    }
  }
  this->Init(o,index);

}
void VecbosPFSC::Init(VecbosBase* o, int i){
  if(i>o->nPhoPFSC || i<0){
    index = -1;
    return;
  }
  index    = i;
  energy   = o->energyPhoPFSC[i];
  esEnergy = o->esEnergyPhoPFSC[i];
  eta      = o->etaPhoPFSC[i];
  phi      = o->phiPhoPFSC[i];
  e3x3     = o->e3x3PhoPFSC[i];
  e5x5     = o->e5x5PhoPFSC[i];
  sigmaIEtaIEta = o->covIEtaIEtaPhoPFSC[i];
  sigmaIEtaIPhi = o->covIEtaIPhiPhoPFSC[i];
  sigmaIPhiIPhi = o->covIPhiIPhiPhoPFSC[i];

  CaloPos.SetXYZ(o->xPosPFSC[i],o->yPosPFSC[i],o->zPosPFSC[i]);


  rawE     = o->rawEnergyPhoPFSC[i];
  phiWidth = o->phiWidthPhoPFSC[i];
  etaWidth = o->etaWidthPhoPFSC[i];
  HoverE   = o->hOverEPhoPFSC[i];
  r9Scale  = 1;  

  //get the basic clusters
  std::vector< std::pair<int,float> > indexEnergyMap;
  for(int j=0; j<o->nPFBC; j++){
    if(o->indexSCPFBC[j] == i){ // the basic cluster points to this SC
      indexEnergyMap.push_back(std::pair<int,float>(j,o->energyPFBC[j]) );
    }
  }
  std::sort(indexEnergyMap.begin(),indexEnergyMap.end(),sort_pred()); // sort the BCs by energy
  for(int j=0;j<indexEnergyMap.size();j++){
    pfClusters.push_back(VecbosPFBC(o,indexEnergyMap.at(j).first));
  }

 }

 VecbosConversion::VecbosConversion(){}

VecbosConversion::VecbosConversion(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosConversion::Init(VecbosBase* o, int i){
  index = i;
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

    trk1Dz      = o->trk1DzConv[i];
    trk1DzError = o->trk1DzErrorConv[i];  
    trk1Charge  = o->trk1ChargeConv[i];   
    trk1Algo    = o->trk1AlgoConv[i];     
    trk1D0      = o->trk1D0Conv[i];       
    trk1Pout    = o->trk1PoutConv[i];
    trk1Pin     = o->trk1PinConv[i];
    
    trk2Dz      = o->trk2DzConv[i];
    trk2DzError = o->trk2DzErrorConv[i];  
    trk2Charge  = o->trk2ChargeConv[i];   
    trk2Algo    = o->trk2AlgoConv[i];     
    trk2D0      = o->trk2D0Conv[i];       
    trk2Pout    = o->trk2PoutConv[i];
    trk2Pin     = o->trk2PinConv[i];
    
  }
}
/*
ConversionInfo VecbosConversion::getStruct(){
  ConversionInfo out = {
    index,
    pPair.X(),
    pPair.Y(),
    pPair.Z(),
    pRefittedPair.X(),
    pRefittedPair.Y(),
    pRefittedPair.Z(),
    p4RefittedPair.Pt(),
    p4RefittedPair.Eta(),
    p4RefittedPair.Phi(),
    p4RefittedPair.Energy(),
    eOverP,
    vtx.X(),
    vtx.Y(),
    vtx.Z(),
    vtxChi2,
    vtxChi2Prob,
    vtxIsValid,
    vtxNTracks,
    vtxMVA,
    trk1Dz,       
    trk1DzError,  
    trk1Charge,   
    trk1Algo,     
    trk1D0,       
    trk1Pout,
    trk1Pin,
    trk2Dz,       
    trk2DzError,  
    trk2Charge,   
    trk2Algo,     
    trk2D0,       
    trk2Pout,     
    trk2Pin
  };
  return out;
}
*/

VecbosPho::VecbosPho():
  CaloPos(0.,0.,0.)
{
  
}

VecbosPho::VecbosPho(VecbosBase* o, int i):
  CaloPos(0.,0.,0.)
{
  this->Init(o,i);
}

void VecbosPho::Init(VecbosBase* o, int i){
  if(i>o->nPho){
    index  = -1;
    return;
  }
  SC.Init(o,o->superClusterIndexPho[i]);
  PFSC.Init(o,SC.eta,SC.phi);
  energy = o->energyPho[i];
  eta    = o->etaPho[i];
  phi    = o->phiPho[i];
  index  = i;

  HoverE = o->hOverEPho[i];
  hasPixel = o->hasPixelSeedPho[i];

  dr03EcalRecHitSumEtCone = o->dr03EcalRecHitSumEtPho[i];
  dr03HcalTowerSumEtCone  = o->dr03HcalTowerSumEtPho[i];
  dr03TrkSumPtCone        = o->dr03TkSumPtPho[i];
  dr03TrkSumPtHollowCone  = o->dr03HollowTkSumPtPho[i]; 

  dr04EcalRecHitSumEtCone = o->dr04EcalRecHitSumEtPho[i];
  dr04HcalTowerSumEtCone  = o->dr04HcalTowerSumEtPho[i];
  dr04TrkSumPtCone        = o->dr04TkSumPtPho[i];
  dr04TrkSumPtHollowCone  = o->dr04HollowTkSumPtPho[i]; 
  int SCI = o->superClusterIndexPho[i];
  CaloPos.SetXYZ(o->xPosSC[SCI],o->yPosSC[SCI],o->zPosSC[SCI]);
};

TLorentzVector VecbosPho::p4FromVtx(TVector3 vtx,float E,bool pf){
  TVector3 dir;
  if(pf) dir = PFSC.CaloPos-vtx;
  else   dir = SC.CaloPos-vtx;
  TVector3 p   = dir.Unit()*E;
  TLorentzVector p4(p.x(),p.y(),p.z(),E);
  return p4;
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

  VecbosConversion c(o,-1);
  conversion = c;
}

/*
PhoInfo VecbosPho::getStruct(){
  PhoInfo out={
    index,
    energy,
    eta,
    phi,
    correctedEnergy,
    correctedEnergyError,
    scaledEnergy,
    smearedEnergy,
    dEoE,
    HoverE,
    hasPixel,
    CaloPos.X(),
    CaloPos.Y(),
    CaloPos.Z(),
    dr03EcalRecHitSumEtCone,
    dr03HcalTowerSumEtCone,
    dr03TrkSumPtCone,
    dr03TrkSumPtHollowCone,
    dr04EcalRecHitSumEtCone,
    dr04HcalTowerSumEtCone,
    dr04TrkSumPtCone,
    dr04TrkSumPtHollowCone,
    conversion->index,
    SC.index
  };
  return out;
}
*/

VecbosEle::VecbosEle(){}

VecbosEle::VecbosEle(VecbosBase* o,int i){
  this->Init(o,i);
}
void VecbosEle::Init(VecbosBase* o,int i){
  if(i>o->nEle) return;
  SC.Init(o,o->superClusterIndexPho[i]);
  energy = o->energyEle[i];
  eta    = o->etaEle[i];
  phi    = o->phiEle[i];
  esEnergy = SC.esEnergy;
  HoverE   = o->hOverEPho[i];

  int tmp = o->recoFlagsEle[i];
  //extract the flag information into booleans
  isTrackerDriven = tmp & 1; 
  isEcalDriven    = (tmp >> 1) & 1;
};
/*
EleInfo VecbosEle::getStruct(){
  EleInfo out={
  index,
  energy,
  eta,
  phi,
  esEnergy,
  HoverE,
  isEcalDriven,
  isTrackerDriven,
  SC.index,
  CaloPos.X(),  
  CaloPos.Y(),  
  CaloPos.Z()  
  };
  return out;
}
*/

VecbosMu::VecbosMu(){}

VecbosMu::VecbosMu(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosMu::Init(VecbosBase* o, int i){
  if(i<0 || i > o->nMuon){
    index = -1;
  }else{
    index = i;            
    energy = o->energyMuon[i];         
    pt = TMath::Sqrt( TMath::Power(o->pxMuon[i],2) + TMath::Power(o->pyMuon[i],2) );
    eta = o->etaMuon[i];
    phi = o->phiMuon[i];
    p4.SetPtEtaPhiE(pt,eta,phi,energy);
    charge = o->chargeMuon[i];
    combinedIso = (o->emEt03Muon[i] + o->hadEt03Muon[i] + o->sumPt03Muon[i] - o->rhoFastjet * TMath::Pi()*0.3*0.3)/pt;;
    emIso = o->emEt03Muon[i];
    hadIso = o->hadEt03Muon[i];
    trkIso = o->sumPt03Muon[i];
    Utils AnalysisUtilities;
    isGlobalMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::AllGlobalMuons);
    isTrackerMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::AllTrackerMuons);
    isPromptMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::GlobalMuonPromptTight);
    int iTrack = o->trackIndexMuon[i];
    if(iTrack<0 || iTrack > o->nTrack) {
      nTrackHits=0;
      nPixelHits=0;
      trackImpactPar=0;
    }else{
      nTrackHits = o->numberOfValidStripTIBHitsTrack[iTrack]
	+ o->numberOfValidStripTIDHitsTrack[iTrack]
	+ o->numberOfValidStripTOBHitsTrack[iTrack]
	+ o->numberOfValidStripTECHitsTrack[iTrack];
      nPixelHits = o->numberOfValidPixelBarrelHitsTrack[iTrack] + o->numberOfValidPixelEndcapHitsTrack[iTrack];
      trackImpactPar = fabs(o->transvImpactParTrack[iTrack]);
    }
    isLooseMuon = true;
    if(!isGlobalMuon || nTrackHits<=10) isLooseMuon=false;
    isTightMuon = isLooseMuon;
    if(!isTrackerMuon || !isPromptMuon || combinedIso >=0.15 || nPixelHits == 0 || trackImpactPar >=0.2) isTightMuon = false;
  }
}