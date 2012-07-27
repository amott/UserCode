#include "VecbosMET.hh"

VecbosMET::VecbosMET(){ index=-1; }

void VecbosMET::Init(VecbosBase *o, int i){
  index=-1;
}

VecbosCaloMET::VecbosCaloMET(){ index = -1; }

VecbosCaloMET::VecbosCaloMET(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosCaloMET::Init(VecbosBase* o, int i){
  if(i!=0){
    index=-1;
    return;
  }

  index=i;
  energy = o->energyMet[i];
  eta    = o->etaMet[i];
  phi    = o->phiMet[i];
  pt     = TMath::Sqrt(TMath::Power(o->pxMet[i],2) + TMath::Power(o->pyMet[i],2) );
  charge = o->chargeMet[i];

  vtxX   = o->vertexXMet[i];
  vtxY   = o->vertexYMet[i];
  vtxZ   = o->vertexZMet[i];

  sumEt  = o->sumEtMet[i];
  mEtSig = o-> mEtSigMet[i];
  significance = o->significanceMet[i];
}

VecbosPFMET::VecbosPFMET(){ index = -1; }

VecbosPFMET::VecbosPFMET(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosPFMET::Init(VecbosBase* o, int i){
  if(i!=0){
    index=-1;
    return;
  }

  index=i;
  energy = o->energyPFMet[i];
  eta    = o->etaPFMet[i];
  phi    = o->phiPFMet[i];
  pt     = TMath::Sqrt(TMath::Power(o->pxPFMet[i],2) + TMath::Power(o->pyPFMet[i],2) );
  charge = o->chargePFMet[i];

  vtxX   = o->vertexXPFMet[i];
  vtxY   = o->vertexYPFMet[i];
  vtxZ   = o->vertexZPFMet[i];

  sumEt  = o->sumEtPFMet[i];
  mEtSig = o-> mEtSigPFMet[i];
  significance = o->significancePFMet[i];
}

VecbosTCMET::VecbosTCMET(){ index = -1; }

VecbosTCMET::VecbosTCMET(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosTCMET::Init(VecbosBase* o, int i){
  if(i!=0){
    index=-1;
    return;
  }

  index=i;
  energy = o->energyTCMet[i];
  eta    = o->etaTCMet[i];
  phi    = o->phiTCMet[i];
  pt     = TMath::Sqrt(TMath::Power(o->pxTCMet[i],2) + TMath::Power(o->pyTCMet[i],2) );
  charge = o->chargeTCMet[i];

  vtxX   = o->vertexXTCMet[i];
  vtxY   = o->vertexYTCMet[i];
  vtxZ   = o->vertexZTCMet[i];

  sumEt  = o->sumEtTCMet[i];
  mEtSig = o-> mEtSigTCMet[i];
  significance = o->significanceTCMet[i];
}

