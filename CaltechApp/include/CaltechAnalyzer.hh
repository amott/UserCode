#ifndef CaltechAnalyzer_hh
#define CaltechAnalyzer_hh

#include "VecbosEGObject.hh"
#include "VecbosJet.hh"
#include "VecbosMu.hh"
#include "VecbosVtx.hh"
#include "VecbosMET.hh"

#include "TTree.h"
#include "TString.h"
#include <vector>

#include "ReadConfig.hh"

#include <BaseAnalyzer.h>

class CaltechAnalyzer{
public:
  CaltechAnalyzer(TTree *tree,TString treeName,TString outputFileName);

  void process();

  bool    getIsValid(){return isValid;}
  TString getProblemDescription(){return reason;} //if we have a processing failure, this explains why

  
private:
  TString __outputFileName;
  TString __outputTreeName;

  TTree *__inputTree;

  void masterBegin(); //reserved: set branch addresses and open output file
  void masterEnd();   //reserved: close files, cleanup 
  void setBranchAddress(); //by default activates everything

  bool isValid; //if this is false, we bail out of the processing
  TString reason; //reason will be printed 
  
protected:
  TH1F  *__outputFile;
  TTree *__outputTree;

  ReadConfig *__config;
  void setConfigFile(TString);

  void masterClear();

  //slave functions run in process:
  virtual void slaveBegin()=0; // any pre-loop setup can be done here
  virtual void processEvent()=0; //this is the basic method to override to do your analysis
  virtual void slaveEnd()=0;   // any post-loop cleanup
  
  void makeInvalid(TString problem=""){isValid=false; reason.Append("; "+problem);}


  void disableCollection(TString colName){outputTree->SetBranchStatus(colName,0);}
  void enableCollection(TString colName){outputTree->SetBranchStatus(colName,1);}
  
  //member data
  int __runNumber;
  int __lumiBlock;
  int __evtNumber;

  float __rho;
  float __rhoEtaMax44;

  const static int maxTriggers=100;
  int __triggers[maxTriggers];
  
  //object collections:
  int __nPho;
  PhoCollection *__Photons;
  int __nMu;
  MuCollection *__Muons;
  int __nEle;
  EleCollection *__Electrons;
  int __nJet;
  JetCollection *__Jets;
  int __nCaloMet;
  CaloMETCollection *__CaloMet;
  int __nPFMet;
  PFMETCollection   *__PFMet;
  int __nTCMet;
  TCMETCollection   *__TCMet;
  int __nVtx;
  VtxCollection *__Vertices;
  //generator collections
  //gen-leve phton
  int __nGenPho;
  GenCollection *__GenPhotons;
  int __nGenMu;
  GenCollection *__GenMuons;
  int __nGenEle;
  GenCollection *__GenElectrons;
  int __nGenHiggs;
  GenCollection *__GenHiggs;
  int __nGenOthers;
  GenCollection *__GenOthers;

  int __nPair;
  std::vector<std::pair<int,int> > *__ggVerticesPhotonIndices;
  std::vector<std::pair<int,float> > *__ggVerticesVertexIndex;  // vertex, MVA score
  std::vector<std::float>            *__ggVerticesPerEvtMVA;
  
  
}

#endif
