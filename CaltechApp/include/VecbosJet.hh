#ifndef VecbosJet_hh
#define VecbosJet_hh

#include "VecbosPhysicsObject.hh"

class VecbosJet : public VecbosPhysicsObject{
public:
  enum JetType{PFPUcorr,PFNoPU};
  VecbosJet();
  VecbosJet(VecbosBase*, int,JetType);
  void Init(VecbosBase*, int,JetType);

  float uncorrEnergy;
  VecbosJet::JetType type;
  
  float area;
  float chargedHadronFraction;
  float neutralHadronFraction;

  float jetIdMva;
  
  float betaStar;
  float betaStarIdMVA;
  float betaStarClassicIdMVA;
  float rmsCands;
  float rmsCandsHand;

  float combinedSecondaryVertex;
  float simpleSecondaryVertexHighPur;
  float simpleSecondaryVertexHighEff;
};

typedef std::vector<VecbosJet> JetCollection;

#endif
