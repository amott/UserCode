#include "VecbosGen.hh"

VecbosGen::VecbosGen(){}

VecbosGen::VecbosGen(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosGen::Init(VecbosBase* o, int i){
  if(i < 0 || i > o->nMc){
    index = -1;
    return;
  }

  index = i;
  pt = o->pMc[i]/cosh(o->etaMc[i]);
  energy = o->energyMc[i];
  eta = o->etaMc[i];
  phi = o->phiMc[i];
  
  vtxX = o->vxMc[i];
  vtxY = o->vyMc[i];
  vtxZ = o->vzMc[i];

  status = o->statusMc[i];
  id     = o->idMc[i];
  indexMother = o->mothMc[i];
  if(indexMother >=0){
    statusMother = o->statusMc[indexMother];
    idMother     = o->idMc[indexMother];
  }else{
    statusMother = -1;
    idMother     =  0;
  }
}
