#ifndef VecbosGen_hh
#define VecbosGen_hh

#include "VecbosBase.hh"
#include "VecbosBaseObject.hh"

class VecbosGen : public VecbosBaseObject{
public:
  VecbosGen();
  VecbosGen(VecbosBase*, int);
  void Init(VecbosBase*, int);

  int status;
  int id;
  int statusMother;
  int idMother;
  int indexMother;
};

typedef std::vector<VecbosGen> GenCollection;

#endif
