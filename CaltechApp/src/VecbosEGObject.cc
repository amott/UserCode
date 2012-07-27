#include <VecbosEGObject.hh>
#include "CommonTools/include/Utils.hh"
#include "VecbosVtx.hh"
#define debugEGObject 0

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

  e2x5Max  = o->e2x5MaxBC[i];
  e2x5Left = o->e2x5LeftBC[i];
  e2x5Right = o->e2x5RightBC[i];
  e2x5Top = o->e2x5TopBC[i];
  e2x5Bottom = o->e2x5BottomBC[i];
  
  etaCrystal = o->etaCrystalBC[i];
  phiCrystal = o->phiCrystalBC[i];
  iEta = o->iEtaBC[i];
  iPhi = o->iPhiBC[i];
  thetaTilt = o->thetaTiltBC[i];
  phiTilt   = o->phiTiltBC[i];

  sigmaIEtaIEta = sqrt(o->covIEtaIEtaBC[i]);
  sigmaIEtaIPhi = o->covIEtaIPhiBC[i];
  sigmaIPhiIPhi = sqrt(o->covIPhiIPhiBC[i]);
};

struct sort_pred {
  bool operator()(const std::pair<int,float> &left, const std::pair<int,float> &right) {
    return left.second > right.second;
  }
};

VecbosSC::VecbosSC(){}

VecbosSC::VecbosSC(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosSC::Init(VecbosBase* o, int i){
  if(i>o->nSC || i<0){
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

  e3x1     = o->e3x1SC[i];
  e1x3     = o->e1x3SC[i];
  e4x4     = o->e4x4SC[i];
  eMax     = o->eMaxSC[i];
  e2x2     = o->e2x2SC[i];
  e2nd     = o->e2ndSC[i];
  e1x5     = o->e1x5SC[i];
  e2x5Max  = o->e2x5MaxSC[i];
  e2x5Left = o->e2x5LeftSC[i];
  e2x5Right = o->e2x5RightSC[i];
  e2x5Top = o->e2x5TopSC[i];
  e2x5Bottom = o->e2x5BottomSC[i];
  
  eLeft    = o->eLeftSC[i];
  eRight   = o->eRightSC[i];
  eTop     = o->eTopSC[i];
  eBottom  = o->eBottomSC[i];


  sigmaIEtaIEta = sqrt(o->covIEtaIEtaSC[i]);
  sigmaIEtaIPhi = o->covIEtaIPhiSC[i];
  sigmaIPhiIPhi = sqrt(o->covIPhiIPhiSC[i]);

  esEffSigRR = TMath::Sqrt( TMath::Power(o->esEffsIxIxSC[i],2) +
			    TMath::Power(o->esEffsIyIySC[i],2) );

  CaloPos.SetXYZ(o->xPosSC[i],o->yPosSC[i],o->zPosSC[i]);

  rawE     = o->rawEnergySC[i];
  phiWidth = o->phiWidthSC[i];
  etaWidth = o->etaWidthSC[i];
  HoverE   = o->hOverESC[i];
  r9       = e3x3/rawE;
  //get the basic clusters


  nBCs = o->nBCSC[i];
  std::vector< std::pair<int,float> > indexEnergyMap;
  float maxE=-1;
  int maxI=-1;
  for(int j=0; j<o->nBC; j++){
    if(o->indexSCBC[j] == i){ // the basic cluster points to this SC
      indexEnergyMap.push_back(std::pair<int,float>(j,o->energyBC[j]) );
      if(o->energyBC[j] > maxE){
	maxE = o->energyBC[j];
	maxI = j;
      }
    }
  }
  
  BCSeed.Init(o,maxI);
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
    vtxX = -999; vtxY = -999; vtxZ = -999;
    eOverP = -999;
    vtxChi2 = 0;
    vtxChi2Prob = 0;
    vtxNTracks  = 0;
    vtxIsValid = 0;
    vtxMVA     = -999;
  }else{
    pPair.SetXYZ(o->pxPairConv[i],o->pyPairConv[i],o->pzPairConv[i]);
    pRefittedPair.SetXYZ(o->pxRefittedPairConv[i],o->pyRefittedPairConv[i],o->pzRefittedPairConv[i]);
    vtxX = o->xVtxConv[i];
    vtxY = o->yVtxConv[i];
    vtxZ = o->zVtxConv[i];

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
  if(i>o->nPho || i<0){
    index  = -1;
    return;
  }
  SC.Init(o,o->superClusterIndexPho[i]);
  energy = o->energyPho[i];
  eta    = o->etaPho[i];
  phi    = o->phiPho[i];
  pt     = TMath::Sqrt(TMath::Power(o->pxPho[i],2) + TMath::Power(o->pyPho[i],2));
  index  = i;

  HoverE = o->hOverEPho[i];
  HTowOverE = o->hTowOverEPho[i];
  hasPixel = o->hasPixelSeedPho[i];

  dr03EcalRecHitSumEtCone = o->dr03EcalRecHitSumEtPho[i];
  dr03HcalTowerSumEtCone  = o->dr03HcalTowerSumEtPho[i];
  dr03TrkSumPtCone        = o->dr03TkSumPtPho[i];
  dr03TrkSumPtHollowCone  = o->dr03HollowTkSumPtPho[i]; 

  dr04EcalRecHitSumEtCone = o->dr04EcalRecHitSumEtPho[i];
  dr04HcalTowerSumEtCone  = o->dr04HcalTowerSumEtPho[i];
  dr04TrkSumPtCone        = o->dr04TkSumPtPho[i];
  dr04TrkSumPtHollowCone  = o->dr04HollowTkSumPtPho[i]; 

  nPV=o->nPV;
  for(int iPV=0;iPV<o->nPV;iPV++){
    dr03HggTrackIso[iPV] = computeTrackIso(o,iPV,0.,0.3,0.2,0.0,1.0,0.1);
    dr04HggTrackIso[iPV] = computeTrackIso(o,iPV,0.,0.4,0.2,0.0,1.0,0.1);

    dr01ChargedHadronPFIso[iPV] = o->dr01ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr02ChargedHadronPFIso[iPV] = o->dr02ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr03ChargedHadronPFIso[iPV] = o->dr03ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr04ChargedHadronPFIso[iPV] = o->dr04ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr05ChargedHadronPFIso[iPV] = o->dr05ChargedHadronPFIsoPho[i*o->nPV+iPV];
    dr06ChargedHadronPFIso[iPV] = o->dr06ChargedHadronPFIsoPho[i*o->nPV+iPV];
  }
  dr01NeutralHadronPFIso  = o->dr01NeutralHadronPFIsoPho[i];
  dr02NeutralHadronPFIso  = o->dr02NeutralHadronPFIsoPho[i];
  dr03NeutralHadronPFIso  = o->dr03NeutralHadronPFIsoPho[i];
  dr04NeutralHadronPFIso  = o->dr04NeutralHadronPFIsoPho[i];
  dr05NeutralHadronPFIso  = o->dr05NeutralHadronPFIsoPho[i];
  dr06NeutralHadronPFIso  = o->dr06NeutralHadronPFIsoPho[i];

  dr01PhotonPFIso  = o->dr01PhotonPFIsoPho[i];
  dr02PhotonPFIso  = o->dr02PhotonPFIsoPho[i];
  dr03PhotonPFIso  = o->dr03PhotonPFIsoPho[i];
  dr04PhotonPFIso  = o->dr04PhotonPFIsoPho[i];
  dr05PhotonPFIso  = o->dr05PhotonPFIsoPho[i];
  dr06PhotonPFIso  = o->dr06PhotonPFIsoPho[i];

  int SCI = o->superClusterIndexPho[i];
  CaloPos.SetXYZ(o->xPosSC[SCI],o->yPosSC[SCI],o->zPosSC[SCI]);

  this->matchElectron(o);


  this->matchConversion(o,true);
  genMatch.Init(o,-1);
};

void VecbosPho::matchElectron(VecbosBase* o){
  eleMatchIndex=-1;
  for(int iEle=0; iEle < o->nEle;iEle++){
    if(o->superClusterIndexEle[iEle] != SC.index) continue;
    if(o->hasMatchedConversionEle[iEle]) continue;
    if(o->gsfTrackIndexEle[iEle]<0 || o->gsfTrackIndexEle[iEle] >=o->nGsfTrack) continue;
    if(o->expInnerLayersGsfTrack[o->gsfTrackIndexEle[iEle]]>0) continue;
    eleMatchIndex = iEle;
    break;
  }
}

TLorentzVector VecbosPho::p4FromVtx(TVector3 vtx,float E){
  TVector3 dir = SC.CaloPos-vtx;
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
       o->nTracksVtxConv[i] != 2 ||
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

  conversion.Init(o,matchIndex);
}

void VecbosPho::doGenMatch(VecbosBase* o, int phoID){
  const float maxDR = 0.2;
  float dEoEBest = 9999;
  int indexGen = -1;
  for(int i=0;i<o->nMc;i++){  
    if(!o->statusMc[i]==1) continue; //require status 1 particles
    if(!o->idMc[i] == phoID) continue; //gen photon
    if(o->energyMc[i] < 1.) continue;
    if(DeltaR(SC.eta,o->etaMc[i],SC.phi,o->phiMc[i]) > maxDR) continue;
    float dEoE = fabs(finalEnergy-o->energyMc[i])/o->energyMc[i];
    if(dEoE > 1.) continue;
    if(dEoE < dEoEBest){
      dEoEBest = dEoE;
      indexGen = i;
    }
  }
  genMatch.Init(o,indexGen);
}

float VecbosPho::computeTrackIso(VecbosBase*o,int vtxI,
				 float ptMin,
				 float outerCone,
				 float innerCone,
				 float etaStripHalfWidth,
				 float dzMax,
				 float dxyMax){
  if(vtxI<0 || vtxI >=o->nPV) return 0;
  float vX=0,vY=0,vZ=0;
  vX = o->PVxPV[vtxI];
  vY = o->PVyPV[vtxI];
  vZ = o->PVzPV[vtxI];


  float SumTrackPt=0;
  for(int iTrack=0;iTrack<o->nTrack;iTrack++){
    float ptTrack = sqrt( o->pxTrack[iTrack]*o->pxTrack[iTrack] + o->pyTrack[iTrack]*o->pyTrack[iTrack] );
    if(ptTrack < ptMin) continue;
    float dZ = fabs( (o->trackVzTrack[iTrack] - vZ) - 
		      ( (o->trackVxTrack[iTrack] - vX)*o->pxTrack[iTrack] 
			+ (o->trackVyTrack[iTrack] - vY)*o->pyTrack[iTrack])/ptTrack*o->pzTrack[iTrack]/ptTrack);
    if(dZ > dzMax) continue;
    float dXY = ( (o->trackVyTrack[iTrack] - vY)*o->pyTrack[iTrack] - (o->trackVxTrack[iTrack] - vX)*o->pxTrack[iTrack] )/ptTrack;
    if( fabs(dXY) > dxyMax ) continue;
    TVector3 trackP(o->pxTrack[iTrack],o->pyTrack[iTrack],o->pzTrack[iTrack]);
    float dEta = fabs(eta - trackP.Eta());
    float dR   = DeltaR(eta,trackP.Eta(),phi,trackP.Phi());
    if( dR < outerCone && dR >= innerCone && dEta >= etaStripHalfWidth) SumTrackPt+=ptTrack;
  }
  return SumTrackPt;
}

VecbosEle::VecbosEle(){}

VecbosEle::VecbosEle(VecbosBase* o,int i){
  this->Init(o,i);
}
void VecbosEle::Init(VecbosBase* o,int i){
  if(i>o->nEle) return;
  SC.Init(o,o->superClusterIndexPho[i]);
  energy = o->energyEle[i];
  regressionEnergy = 0;

  eta    = o->etaEle[i];
  phi    = o->phiEle[i];
  pt     = TMath::Sqrt(TMath::Power(o->pxEle[i],2)+TMath::Power(o->pyEle[i],2));

  esEnergy = SC.esEnergy;
  HoverE   = o->hOverEPho[i];

  vtxX = o->vertexXEle[i];
  vtxY = o->vertexYEle[i];
  vtxZ = o->vertexZEle[i];

  EOverP = o->eSuperClusterOverPEle[i];

  dr03ChargedHadronPFIso = o->pfCandChargedIso03Ele[i];
  dr03NeutralHadronPFIso = o->pfCandNeutralIso03Ele[i];
  dr03PhotonPFIso        = o->pfCandPhotonIso03Ele[i];

  dr04ChargedHadronPFIso = o->pfCandChargedIso04Ele[i];
  dr04NeutralHadronPFIso = o->pfCandNeutralIso04Ele[i];
  dr04PhotonPFIso        = o->pfCandPhotonIso04Ele[i];

  dr03TkSumPt          = o->dr03TkSumPtEle[i];
  dr03EcalRecHitSumEt  = o->dr03EcalRecHitSumEtEle[i];
  dr03HcalTowerSumEt   = o->dr03HcalTowerSumEtEle[i];

  dr04TkSumPt          = o->dr04TkSumPtEle[i];
  dr04EcalRecHitSumEt  = o->dr04EcalRecHitSumEtEle[i];
  dr04HcalTowerSumEt   = o->dr04HcalTowerSumEtEle[i];

  dEtaSCTrack = o->deltaEtaEleClusterTrackAtCaloEle[i];
  dPhiSCTrack = o->deltaPhiEleClusterTrackAtCaloEle[i];

  int tmp = o->recoFlagsEle[i];
  //extract the flag information into booleans
  isTrackerDriven = tmp & 1; 
  isEcalDriven    = (tmp >> 1) & 1;

  genMatch.Init(o,-1);
};

