#include "VecbosPhysicsObject.hh"
#include "HggPhysUtils.cc"

VecbosPhysicsObject::VecbosPhysicsObject(){}

void VecbosPhysicsObject::doGenMatch(VecbosBase* o, int partId){
  const float maxDR = 0.2;
  float dEoEBest = 9999;
  int indexGen = -1;
  for(int i=0;i<o->nMc;i++){  
    if(!o->statusMc[i]==1) continue; //require status 1 particles
    if(abs(o->idMc[i]) != partId) continue; //gen mu
    if(o->energyMc[i] < 1.) continue;
    if(DeltaR(eta,o->etaMc[i],phi,o->phiMc[i]) > maxDR) continue;
    float dEoE = fabs(energy-o->energyMc[i])/o->energyMc[i];
    if(dEoE > 0.5) continue;
    if(dEoE < dEoEBest){
      dEoEBest = dEoE;
      indexGen = i;
    }
  }
  genMatch.Init(o,indexGen);

}
