#ifndef VecbosBaseObject_hh
#define VecbosBaseObject_hh

#include "TLorentzVector.h"

class VecbosBaseObject{
public:
  VecbosBaseObject();
  int index;
  float energy;
  float pt;
  float eta;
  float phi;

  int charge;

  float vtxX;
  float vtxY;
  float vtxZ;
  TVector3 getVertex(){return TVector3(vtxX,vtxY,vtxZ);}

  TLorentzVector getP4();

  bool isValid(){return (index!=-1);}
};


typedef std::vector<VecbosBaseObject> BaseObjectCollection;

#endif
