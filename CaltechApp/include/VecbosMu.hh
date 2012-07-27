#ifndef VecbosMu_hh
#define VecbosMu_hh

#include "VecbosPhysicsObject.hh"

class VecbosMu : public VecbosPhysicsObject{
public:
  VecbosMu();
  VecbosMu(VecbosBase*, int);
  void Init(VecbosBase*, int);

  float combinedIso;
  float emIso;
  float hadIso;
  float trkIso;
  bool isGlobalMuon;
  bool isTrackerMuon;
  bool isPromptMuon;
  int nTrackHits;
  int nPixelHits;
  float trackImpactPar;

  bool isLooseMuon;
  bool isTightMuon;
};

typedef std::vector<VecbosMu> MuCollection;
#endif
