#ifndef VecbosPhysicsObject_hh
#define VecbosPhysicsObject_hh

#include "VecbosBaseObject.hh"
#include "VecbosGen.hh"

class VecbosPhysicsObject : public VecbosBaseObject{
public:
  VecbosPhysicsObject();
  VecbosGen genMatch;
  virtual void doGenMatch(VecbosBase* o,int);

};
typedef std::vector<VecbosPhysicsObject> PhysicsObjectCollection;

#endif
