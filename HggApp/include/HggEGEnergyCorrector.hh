#include "VecbosBase.hh"
#include <string>

#include "../MitPhysics/Utils/interface/GBRForest.h"
#include "../MitPhysics/Utils/interface/GBRTree.h"                         
//#include "h2gglobe/VertexAnalysis/interface/HggVertexAnalyzer.h"
//#include "h2gglobe/VertexAnalysis/interface/PhotonInfo.h"
//#include "h2gglobe/VertexAnalysis/interface/VertexAlgoParameters.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"

#include "../src/ecalGap.cc"
using namespace TMVA;

class HggEGEnergyCorrector{
 public:
  HggEGEnergyCorrector(VecbosBase*,int,Bool_t);
  std::pair<double,double> CorrectedEnergyWithError(int);
  std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithError(int);
  std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithErrorv2(int);
  std::pair<double,double> photonEnergyCorrector_CorrectedEnergyWithErrorv2(int);
 private:
  //private methods
  void Init(int,Bool_t);
  //vars
  
  ECAL_GEO ecalGeometry;

  VecbosBase *base;

  GBRForest *fReadereb;
  GBRForest *fReaderebvariance;
  GBRForest *fReaderee;
  GBRForest *fReadereevariance;
  Float_t *fVals;
};
