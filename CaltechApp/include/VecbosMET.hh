#ifndef VecbosMET_hh
#define VecbosMET_hh

#include "VecbosPhysicsObject.hh"

class VecbosMET : public VecbosPhysicsObject{ //pure virtual base class
public:
  VecbosMET();
  virtual void Init(VecbosBase*, int);
  
  float sumEt;
  float mEtSig;
  float significance;
};

class VecbosCaloMET : public VecbosMET{
public:
  VecbosCaloMET();
  VecbosCaloMET(VecbosBase*, int);
  void Init(VecbosBase*, int);
};

class VecbosPFMET : public VecbosMET{
public:
  VecbosPFMET();
  VecbosPFMET(VecbosBase*, int);
  void Init(VecbosBase*, int);
};

class VecbosTCMET : public VecbosMET{
public:
  VecbosTCMET();
  VecbosTCMET(VecbosBase*, int);
  void Init(VecbosBase*, int);
};

typedef std::vector<VecbosMET> METCollection;
typedef std::vector<VecbosCaloMET> CaloMETCollection;
typedef std::vector<VecbosPFMET> PFMETCollection;
typedef std::vector<VecbosTCMET> TCMETCollection;

#endif
