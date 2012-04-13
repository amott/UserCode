 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef HggVertexing_h
#define HggVertexing_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<unsigned short> >+; 
#endif
//#include "TMVAGui.C"                         
//#if not defined(__CINT__) || defined(__MAKECINT__) 
#include "../h2gglobe/VertexAnalysis/interface/HggVertexAnalyzer.h"
#include "../h2gglobe/VertexAnalysis/interface/PhotonInfo.h"
#include "../h2gglobe/VertexAnalysis/interface/VertexAlgoParameters.h"
                 
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
//#endif 
using namespace TMVA; 
using namespace std;
#include "VecbosEGObject.hh"


class HggVertexing{
public:
  HggVertexing(VecbosBase*); /// Class Constructor
  virtual ~HggVertexing();     /// Class Destructor
  void setConfigFile(string s){configFilePath = s;}

  int vertex_tmva(int,int);
  void useConversions(bool b=true){useConversion = b;}
private:
  VecbosBase *base;
  void init(); // do variable initialization 
  
  TMVA::Reader * perVtxReader; 
  TMVA::Reader * perEvtReader;
  vector<string> rankVariables; 

  //read these from the config file
  string configFilePath;
  string perVtxMvaWeights;
  string perVtxMvaMethod;
  string perEvtMvaWeights; 
  string perEvtMvaMethod;

  vector<int> rankmethod;
  vector<string> varNameUsed;
  vector<short> *indvertexSelected_allpairpresel;
  vector<vector<float> > *photontrkisoselvtxdr03;

  VertexAlgoParameters vtxAlgoParams;
  HggVertexAnalyzer vAna;

  //options
  bool useConversion;
};
#endif


