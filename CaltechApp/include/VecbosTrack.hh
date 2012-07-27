#ifndef VecbosTrack_hh
#define VecbosTrack_hh

#include "VecbosBaseObject.hh"

class VecbosTrack : public VecbosBaseObject{
public:
  enum TrackType{std,gsf,muon};
  VecbosTrack();
  VecbosTrack(VecbosBase*,int,TrackType);
  void init(VecbosBase*,int,TrackType);
  int vtxIndex;
  float vtxWeight;


}

#endif
