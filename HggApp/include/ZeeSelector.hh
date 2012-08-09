#include <VecbosEGObject.hh>
#include <HggEGEnergyCorrector.hh>
#include <HggVertexing.hh>
#include <HggEnergyScale.hh>

#include <vector>
#include <iostream>
#include <string>

#include <TChain.h>
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TH1F.h"

using namespace std;
#include "HggVertexing.hh"

class ZeeSelector{
public:
  ZeeSelector();
  ~ZeeSelector();
  ZeeSelector(vector<string> fNames,string treeName,string outputFile);
  void loadChain(vector<string> fNames, string treeName);
  void setOutputFile(string s){outputFile = s;}
  bool isValid(){return valid;}
  void setIsData(bool d){isData_=d;}
  void Loop();
private:
  bool valid;;
  TChain* fChain;
  TTree* outTree;
  TTree* outTreeZee;
  string outputFile;

  bool isData_;

  int init();
  void setBranchAddresses();
  void setupOutputTree();
  

  // Mass Selection
  float DZmassref; 
  float DZmass;
  float Zeemass;
  int   nEle;
  int   lpass;
  int   tpass;
  int   mvapass;

  // Input variables
  int runNumber;
  int evtNumber;
  int nVtx;

  // Electron Selection
  static const int MAX = 100;
  int   charge[MAX];
  float eta[MAX];
  float phi[MAX];

  float energySC[MAX];
  float etaSC[MAX];
  float phiSC[MAX];
  float r9[MAX];

  float HoverE[MAX];
  float dEta[MAX];
  float dPhi[MAX];
  float sigmaIEtaIEta[MAX];
  float PFIsoOverPT1;
  float PFIsoOverPT2;
  int   d0[MAX];
  int   dz[MAX];

  float idMVA[MAX];

  float PT[MAX];
  float chargedPFiso[MAX];
  float neutralPFiso[MAX];
  float photonPFiso[MAX];
  float rho;

  bool  hasMatchedConversion[MAX];
  int   expInnerLayersHits[MAX];

  // Variables that will be outputted
  float mass;
  float Ele1mva;
  float Ele2mva;
  int   nEleOut;
  float Ele1eta;
  float Ele2eta;
  float Ele1r9;
  float Ele2r9;
  int   passtight;
  int   passmva;	
  int   passloose;
};
