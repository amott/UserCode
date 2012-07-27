#ifndef VecbosVtx_hh
#define VecbosVtx_hh

#include "VecbosBase.hh"
#include "VecbosBaseObject.hh"

class VecbosVtx : public VecbosBaseObject{
public:
  VecbosVtx();
  VecbosVtx(VecbosBase*,int);
  void Init(VecbosBase*,int);
  float x,y,z;
  float xe,ye,ze;

  float pxChargedMet,pyChargedMet,pzChargedMet;

  float SumPt;
  float ndof;
  float rho;
  float chi2;
  int trackSize;
};

typedef std::vector<VecbosVtx> VtxCollection;
#endif
