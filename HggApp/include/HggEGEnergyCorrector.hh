#include "VecbosBase.hh"
#include <string>

#include "GBRForest.h"
#include "GBRTree.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"

#include "VecbosEGObject.hh"

#include "../src/ecalGap.cc"
using namespace TMVA;

class HggEGEnergyCorrector{
 public:
  HggEGEnergyCorrector(VecbosBase*,string,Bool_t);
  std::pair<double,double> getPhotonEnergyCorrection(int);
  //std::pair<double,double> getElectronEnergyCorrection(int);

  std::pair<double,double> CorrectedEnergyWithError(int);
  std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithError(int);
  std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithErrorv2(int);
  std::pair<double,double> photonEnergyCorrector_CorrectedEnergyWithErrorv2(VecbosPho&);

  std::pair<double,double> photonEnergyCorrector_May2012(VecbosPho&);
  std::pair<double,double> electronEnergyCorrector_May2012(int);
 private:
  //private methods
  void Init();
  //vars
  string configFile;
  Bool_t isRealData;
  string version;
  ECAL_GEO ecalGeometry;

  VecbosBase *base;

  GBRForest *fReadereb;
  GBRForest *fReaderebvariance;
  GBRForest *fReaderee;
  GBRForest *fReadereevariance;
  Float_t *fVals;
};
