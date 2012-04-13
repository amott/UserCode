//read the Hgg Energy Scale as a csv

#ifndef HggEnergyScale_hh
#define HggEnergyScale_hh

#include "ReadConfig.hh"
#include <string>
#include <map>
#include <vector>
#include "VecbosEGObject.hh"
//#include "TString.h"


class HggEnergyScale{
public:
  HggEnergyScale(string);
  float getDEoE(VecbosPho,int);
  bool isValid();{return valid;}
private:
  bool valid;
   const static int nRegions = 8;
  const static string configNames[nRegions] = {
    "EBlowEtaBadDeltaE",
    "EBlowEtaGoldDeltaE",
    "EBhiEtaBadDeltaE",
    "EBhiEtaGoldDeltaE",
    "EElowEtaBadDeltaE",
    "EElowEtaGoldDeltaE"
    "EEhiEtaBadDeltaE",
    "EEhiEtaGoldDeltaE",
  };

  const static float r9Cut = 0.94;
  const static bool highR9 = {false,true,false,true,false,true,false,true};
  const static float minEta = {0.,0.,1.,1.,1.48,1.48,2.,2.};
  const static float maxEta = {1.,1.,1.48,1.48,2.,2.,3.,3.};

  std::vector<int> runs;
  std::vector<float> energyScales[nRegions]; 

};

HggEnergyScale::HggEnergyScale(string path){
  valid= false; // only s

  ReadConfig reader(path);

  if(!reader.is_init()) return;

  //std::vector<TString> strings;
  //for(int i=0;i<nRegions;i++) strings.push_back( (TString)reader.getParameter(configNames[i]) );
  
  std::string runString = reader.getParameter("RunUpper");
  //tokenize the string
  char *rs = strtok(runString.c_str(),",");
  while(rs){
    runs.push_back(atoi(rs));
    rs = strtok(NULL,",");
  }

  for(int i=0;i<nRegions;i++){
    std::string valString = reader.getParameter(configNames[i]);
    char *vs = strtok(valString.c_str(),",");
    while(vs){
      energyScales[i].push_back(atof(vs));
      vs = strtok(NULL,",");
    }
  }

  valid = true;
  for(int i=0;i<nRegions;i++){
    if(runs.size()!=energyScales[i].size()) valid=false;
  }

}

float HggEnergyScale::getDEoE(VecbosPho pho, int run){
  if(!valid) return -9999;

  int runIndex = 0;
  for(int i=0;i<runs.size();i++){
    if(run <= runs.at(i)) break;
    runIndex++;
  }

  int selectRegion=-1;
  for(int iReg = 0; iReg< nRegions; iReg++){
    if( std::fabs(pho.eta) > minEta[iReg] 
	&& std::fabs(pho.eta) <=maxEta[iReg] 
	&& ((pho.SC.r9() > r9Cut) == highR9[iReg]) ){
      selectRegion == iReg;
      break;
    }
  }

  if(selectRegion == -1) return -999;

  return energyScales[selectRegion].at(runIndex);
}
#endif
