#include "VecbosVtx.hh"

VecbosVtx::VecbosVtx(){}

VecbosVtx::VecbosVtx(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosVtx::Init(VecbosBase* o, int i){
  if(i<0 || i >= o->nPV){
    index = -1;
    return;
  }

  x = o->PVxPV[i];
  y = o->PVyPV[i];
  z = o->PVzPV[i];

  xe = o->PVErrxPV[i];
  ye = o->PVErryPV[i];
  ze = o->PVErrzPV[i];

  pxChargedMet = o->pxChMetPV[i];
  pyChargedMet = o->pyChMetPV[i];
  pzChargedMet = o->pzChMetPV[i];
  
  SumPt = o->SumPtPV[i];
  ndof = o->ndofPV[i];
  rho  = o->rhoPV[i];
  chi2 = o->chi2PV[i];
  trackSize = o->trackSizePV[i];
}
