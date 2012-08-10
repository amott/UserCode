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
  bool  isZmass;
  float rho;
  float PFIsoOverPT1;
  float PFIsoOverPT2;

  // Input variables
  int runNumber;
  int evtNumber;
  int nVtx;

  // Electron Selection
  std::vector<VecbosEle> *Electrons;

  // Variables that will be outputted
  float mass;
  int   nEleOut;
  float Ele1mva;
  float Ele2mva;
  float Ele1eta;
  float Ele2eta;
  float Ele1r9;
  float Ele2r9;
  int   passtight;
  int   passmva;	
  int   passloose;
};
