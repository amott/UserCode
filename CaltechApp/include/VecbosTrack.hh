#ifndef VecbosTrack_hh
#define VecbosTrack_hh

#include "VecbosBaseObject.hh"
#include "VecbosBase.hh"

class VecbosTrack : public VecbosBaseObject{
public:
  VecbosTrack();
  VecbosTrack(VecbosBase*,int);
  virtual void Init(VecbosBase*,int);
  int vtxIndex;
  float vtxWeight;

  float d0;
  float d0Error;
  float dz;
  float dzError;
  
  int expInnerLayersHits;

  int nValidBPIX;
  int nValidFPIX;
  int nValidTIB;
  int nValidTID;
  int nValidTOB;
  int nValidTEC;
  int nValidMuons;
};

class VecbosGsfTrack : public VecbosTrack{
public:
  VecbosGsfTrack();
  VecbosGsfTrack(VecbosBase*,int);
  virtual void Init(VecbosBase*, int);
};

class VecbosMuonTrack : public VecbosTrack{
public:
  VecbosMuonTrack();
  VecbosMuonTrack(VecbosBase*,int);
  virtual void Init(VecbosBase*, int);
};

typedef std::vector<VecbosTrack> TrackCollection;
typedef std::vector<VecbosGsfTrack> GsfTrackCollection;
typedef std::vector<VecbosMuonTrack> MuonTrackCollection;

#endif
