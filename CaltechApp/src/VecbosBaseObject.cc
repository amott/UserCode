#include "VecbosBaseObject.hh"

VecbosBaseObject::VecbosBaseObject()
{
  index   = -1;
  energy  = 0.0001;
  pt      = 0.0001;
  eta     = 0.;
  phi     = 0.;
  charge  = 0;
}

TLorentzVector VecbosBaseObject::getP4(){
  TLorentzVector p4; p4.SetPtEtaPhiE(pt,eta,phi,energy);
  return p4;
}

