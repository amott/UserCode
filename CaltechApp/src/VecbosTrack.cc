#include "VecbosTrack.hh"
#include "TMath.h"
VecbosTrack::VecbosTrack(){index=-1;}

VecbosTrack::VecbosTrack(VecbosBase *o, int i){
  this->Init(o,i);
}

void VecbosTrack::Init(VecbosBase *o, int i){
  if(i<0 || i>o->nTrack){
    index=-1;
    return;
  }
  index=i;
  energy=TMath::Sqrt( TMath::Power(o->pxTrack[i],2) +
		       TMath::Power(o->pyTrack[i],2) +
		       TMath::Power(o->pzTrack[i],2) );
  pt = TMath::Sqrt( TMath::Power(o->pxTrack[i],2) +
		     TMath::Power(o->pyTrack[i],2) );
  eta = acosh(energy/pt);
  phi = acos(o->pyTrack[i]/o->pxTrack[i]);

  charge = o->chargeTrack[i];

  vtxX = o->trackVxTrack[i];
  vtxY = o->trackVyTrack[i];
  vtxZ = o->trackVzTrack[i];
		     
  vtxIndex = o->vtxIndexTrack[i];
  vtxWeight = o->vtxWeightTrack[i];

  d0 = o->d0Track[i];
  d0Error = o->d0ErrorTrack[i];
  dz = o->dzTrack[i];
  dzError = o->dzErrorTrack[i];
  
  expInnerLayersHits = o->expInnerLayersTrack[i];
  
  nValidBPIX = o->numberOfValidPixelBarrelHitsTrack[i];
  nValidFPIX = o->numberOfValidPixelEndcapHitsTrack[i];
  nValidTIB  = o->numberOfValidStripTIBHitsTrack[i];
  nValidTID  = o->numberOfValidStripTIDHitsTrack[i];
  nValidTOB  = o->numberOfValidStripTOBHitsTrack[i];
  nValidTEC  = o->numberOfValidStripTECHitsTrack[i];
  nValidMuons = o->numberOfValidMuonHitsTrack[i];
}

VecbosGsfTrack::VecbosGsfTrack(){index=-1;}

VecbosGsfTrack::VecbosGsfTrack(VecbosBase *o, int i){
  this->Init(o,i);
}

void VecbosGsfTrack::Init(VecbosBase *o, int i){
  if(i<0 || i>o->nGsfTrack){
    index=-1;
    return;
  }
  index=i;
  energy=TMath::Sqrt( TMath::Power(o->pxGsfTrack[i],2) +
		       TMath::Power(o->pyGsfTrack[i],2) +
		      TMath::Power(o->pzGsfTrack[i],2) );
  pt = TMath::Sqrt( TMath::Power(o->pxGsfTrack[i],2) +
		     TMath::Power(o->pyGsfTrack[i],2) );
  eta = acosh(energy/pt);
  phi = acos(o->pyGsfTrack[i]/o->pxGsfTrack[i]);

  charge = o->chargeGsfTrack[i];

  vtxX = o->trackVxGsfTrack[i];
  vtxY = o->trackVyGsfTrack[i];
  vtxZ = o->trackVzGsfTrack[i];
		     
  vtxIndex = o->vtxIndexGsfTrack[i];
  vtxWeight = o->vtxWeightGsfTrack[i];

  d0 = o->d0GsfTrack[i];
  d0Error = o->d0ErrorGsfTrack[i];
  dz = o->dzGsfTrack[i];
  dzError = o->dzErrorGsfTrack[i];
  
  expInnerLayersHits = o->expInnerLayersGsfTrack[i];
  
  nValidBPIX = o->numberOfValidPixelBarrelHitsGsfTrack[i];
  nValidFPIX = o->numberOfValidPixelEndcapHitsGsfTrack[i];
  nValidTIB  = o->numberOfValidStripTIBHitsGsfTrack[i];
  nValidTID  = o->numberOfValidStripTIDHitsGsfTrack[i];
  nValidTOB  = o->numberOfValidStripTOBHitsGsfTrack[i];
  nValidTEC  = o->numberOfValidStripTECHitsGsfTrack[i];
  nValidMuons = o->numberOfValidMuonHitsGsfTrack[i];
}

VecbosMuonTrack::VecbosMuonTrack(){index=-1;}

VecbosMuonTrack::VecbosMuonTrack(VecbosBase *o, int i){
  this->Init(o,i);
}

void VecbosMuonTrack::Init(VecbosBase *o, int i){
  if(i<0 || i>o->nGlobalMuonTrack){
    index=-1;
    return;
  }
  index=i;
  energy=TMath::Sqrt( TMath::Power(o->pxGlobalMuonTrack[i],2) +
		       TMath::Power(o->pyGlobalMuonTrack[i],2) +
		       TMath::Power(o->pzGlobalMuonTrack[i],2) );
  pt = TMath::Sqrt( TMath::Power(o->pxGlobalMuonTrack[i],2) +
		     TMath::Power(o->pyGlobalMuonTrack[i],2) );
  eta = acosh(energy/pt);
  phi = acos(o->pyGlobalMuonTrack[i]/o->pxGlobalMuonTrack[i]);

  charge = o->chargeGlobalMuonTrack[i];

  vtxX = o->trackVxGlobalMuonTrack[i];
  vtxY = o->trackVyGlobalMuonTrack[i];
  vtxZ = o->trackVzGlobalMuonTrack[i];
		     
  vtxIndex = o->vtxIndexGlobalMuonTrack[i];
  vtxWeight = o->vtxWeightGlobalMuonTrack[i];

  d0 = o->d0GlobalMuonTrack[i];
  d0Error = o->d0ErrorGlobalMuonTrack[i];
  dz = o->dzGlobalMuonTrack[i];
  dzError = o->dzErrorGlobalMuonTrack[i];
  
  expInnerLayersHits = o->expInnerLayersGlobalMuonTrack[i];
  
  nValidBPIX = o->numberOfValidPixelBarrelHitsGlobalMuonTrack[i];
  nValidFPIX = o->numberOfValidPixelEndcapHitsGlobalMuonTrack[i];
  nValidTIB  = o->numberOfValidStripTIBHitsGlobalMuonTrack[i];
  nValidTID  = o->numberOfValidStripTIDHitsGlobalMuonTrack[i];
  nValidTOB  = o->numberOfValidStripTOBHitsGlobalMuonTrack[i];
  nValidTEC  = o->numberOfValidStripTECHitsGlobalMuonTrack[i];
  nValidMuons = o->numberOfValidMuonHitsGlobalMuonTrack[i];
}

