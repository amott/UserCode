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
  void getPhotonEnergyCorrection(VecbosPho&);
  void getElectronEnergyCorrection(VecbosEle&);
  //std::pair<double,double> getElectronEnergyCorrection(int);

  std::pair<double,double> CorrectedEnergyWithError(int);
  std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithError(int);
  std::pair<double,double> electronEnergyCorrector_CorrectedEnergyWithErrorv2(VecbosEle&);
  std::pair<double,double> photonEnergyCorrector_CorrectedEnergyWithErrorv2(VecbosPho&);

  std::pair<double,double> photonEnergyCorrector_May2012(VecbosPho&);
  std::pair<double,double> electronEnergyCorrector_May2012(VecbosEle&);

  void useElectronWeights(){usePhoton=false;}
 private:
  //private methods
  void Init();
  //vars
  string configFile;
  Bool_t isRealData;
  Bool_t usePhoton;
  string version;
  ECAL_GEO ecalGeometry;

  VecbosBase *base;

  GBRForest *fReadereb;
  GBRForest *fReaderebvariance;
  GBRForest *fReaderee;
  GBRForest *fReadereevariance;
  Float_t *fVals;
};
