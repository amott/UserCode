

#include <VecbosEGObject.hh>

#include <vector>
#include <iostream>
#include <string>

#include <TChain.h>
#include <TH1F.h>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TVector3.h"
using namespace std;

                                    
class HggPhotonID{
public:
  HggPhotonID();
  ~HggPhotonID();
  void setConfig(string s){configFile = s;}
  bool isValid(){return valid;}
  void Init();

  bool getPreSelection(VecbosPho*,int,float,TVector3,int);
  float getIdMVA(VecbosPho*,int,float, TVector3, int);
  bool getIdCiCPF(VecbosPho*,int,float, TVector3, int);
  bool getIdCiC(VecbosPho*,int,float, TVector3, int);
  
  int getCiCCat(VecbosPho*);

  bool getPreSelection2011(VecbosPho*,int,float,TVector3,int);
  bool getPreSelectionMay2012(VecbosPho*,int,float,TVector3,int);

  std::map<std::string,TH1F*>* getHists(){return &InputDists;}
private:
  string configFile;
  bool valid;

  string version;

  const static float isoSumConst = 0;
  const static float isoSumConstPF = 2.5;
  const static float rhoFac = 0.09;
  const static float rhoFacBad = 0.23;

  string weightFile_IdEB_2011;
  string weightFile_IdEE_2011;
  string weightFile_IdEB_2012;
  string weightFile_IdEE_2012;
  string methodName_Id;

  TMVA::Reader *photonMVA_EB_2011;
  TMVA::Reader *photonMVA_EE_2011;
  TMVA::Reader *photonMVA_EB_2012;
  TMVA::Reader *photonMVA_EE_2012;
  void setupTMVA();

  std::map<std::string,TH1F*> InputDists;

  void fillVariables(VecbosPho*,int,float,TVector3,int);
  //photon ID MVA
  float hoe;
  float sigietaieta;
  float isosumoet;
  float isosum;
  float trkisooet;
  float isosumoetbad;
  float isosumbad;
  float r9;
  float ecalisodr03;
  float ecalisodr04;
  float hcalisodr03;
  float hcalisodr04;
  float nVertexf;
  float etasc;
  float scetawidth;
  float scphiwidth;

  float pfChargedIsoGood03;
  float pfChargedIsoBad04;
  float pfChargedIsoGood03oet;
  float pfChargedIsoBad04oet;
  float pfPhotonIso03;
  float pfPhotonIso03oet;
  float pfPhotonIso04;
  float pfPhotonIso04oet;
  float sigietaiphi;
  float s4Ratio;
  float rho;
  float sigRR;
  float isosumoetPF;
  float isosumoetbadPF;
  float isosumPF;
  float isosumbadPF;
  
};
