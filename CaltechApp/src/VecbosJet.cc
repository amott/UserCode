#include "VecbosJet.hh"

VecbosJet::VecbosJet(){}

VecbosJet::VecbosJet(VecbosBase* o, int i,VecbosJet::JetType jtype = PFPUcorr){
  this->Init(o,i,jtype);
}

void VecbosJet::Init(VecbosBase* o, int i,VecbosJet::JetType jtype = PFPUcorr){
  type = jtype;
  switch(type){
  case PFPUcorr:
    if(i < 0 || i > o->nAK5PFPUcorrJet){
      index = -1;
      return;
    }
    index = i;
   energy = o->energyAK5PFPUcorrJet[i];
   uncorrEnergy = o->uncorrEnergyAK5PFPUcorrJet[i];
   eta = o->etaAK5PFPUcorrJet[i];
   phi = o->phiAK5PFPUcorrJet[i];
   pt = TMath::Sqrt(TMath::Power(o->pxAK5PFPUcorrJet[i],2) + TMath::Power(o->pyAK5PFPUcorrJet[i],2));
   
   charge = o->chargeAK5PFPUcorrJet[i];

   vtxX = o->vertexXAK5PFPUcorrJet[i];
   vtxY = o->vertexYAK5PFPUcorrJet[i];
   vtxZ = o->vertexZAK5PFPUcorrJet[i];

   area =  o->areaAK5PFPUcorrJet[i];
   chargedHadronFraction = o->chargedHadronEnergyAK5PFPUcorrJet[i];
   neutralHadronFraction = o->neutralHadronEnergyAK5PFPUcorrJet[i];

   jetIdMva = o->jetIdMvaAK5PFPUcorrJet[i];
  
   betaStar = o->betastarAK5PFPUcorrJet[i];
   betaStarIdMVA = o->betaIdMvaAK5PFPUcorrJet[i];
   betaStarClassicIdMVA = o->betastarclassicIdMvaAK5PFPUcorrJet[i];
   rmsCands = o->rmsCandAK5PFPUcorrJet[i];
   rmsCandsHand = o->rmsCandsHandAK5PFPUcorrJet[i];

   combinedSecondaryVertex = o->combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[i];
   simpleSecondaryVertexHighPur = o->simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[i];
   simpleSecondaryVertexHighEff = o->simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[i];
   break;

  case PFNoPU:
        if(i < 0 || i > o->nAK5PFNoPUJet){
      index = -1;
      return;
    }
    index = i;
   energy = o->energyAK5PFNoPUJet[i];
   uncorrEnergy = o->uncorrEnergyAK5PFNoPUJet[i];
   eta = o->etaAK5PFNoPUJet[i];
   phi = o->phiAK5PFNoPUJet[i];
   pt = TMath::Sqrt(TMath::Power(o->pxAK5PFNoPUJet[i],2) + TMath::Power(o->pyAK5PFNoPUJet[i],2));
   
   charge = o->chargeAK5PFNoPUJet[i];

   vtxX = o->vertexXAK5PFNoPUJet[i];
   vtxY = o->vertexYAK5PFNoPUJet[i];
   vtxZ = o->vertexZAK5PFNoPUJet[i];

   area =  o->areaAK5PFNoPUJet[i];
   chargedHadronFraction = o->chargedHadronEnergyAK5PFNoPUJet[i];
   neutralHadronFraction = o->neutralHadronEnergyAK5PFNoPUJet[i];

   jetIdMva = o->jetIdMvaAK5PFNoPUJet[i];
  
   betaStar = o->betastarAK5PFNoPUJet[i];
   betaStarIdMVA = o->betaIdMvaAK5PFNoPUJet[i];
   betaStarClassicIdMVA = o->betastarclassicIdMvaAK5PFNoPUJet[i];
   rmsCands = o->rmsCandAK5PFNoPUJet[i];
   rmsCandsHand = o->rmsCandsHandAK5PFNoPUJet[i];

   combinedSecondaryVertex = o->combinedSecondaryVertexBJetTagsAK5PFNoPUJet[i];
   simpleSecondaryVertexHighPur = o->simpleSecondaryVertexHighPurBJetTagsAK5PFNoPUJet[i];
   simpleSecondaryVertexHighEff = o->simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet[i];
   break;
 
 Default:
   index = -1;
   break;
  }

}

