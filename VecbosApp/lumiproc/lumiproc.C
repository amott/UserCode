#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "cajun/json/elements.h"
#include "cajun/json/reader.h"
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

typedef std::pair<unsigned int,unsigned int> aLSSegment;
typedef std::vector< std::pair<unsigned int,unsigned int> > LSSegments;
typedef unsigned int aRun;
typedef std::map< aRun, LSSegments > runsLSSegmentsMap;
typedef std::pair < aRun, LSSegments > aRunsLSSegmentsMapElement;

runsLSSegmentsMap* fillRunLSMap(char* jsonFile){
  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile);
  if (!jsonFileStream.is_open()){
    std::cout << "Unable to open file " << jsonFile << std::endl;
    return 0 ;
  }

  runsLSSegmentsMap* goodRunLS = new runsLSSegmentsMap;
  

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile,jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun){
    const json::Array& lsSegment = (*itRun).element;
    LSSegments thisRunSegments;
    for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator){
      json::Array lsSegment=(*lsIterator);
      json::Number lsStart=lsSegment[0];
      json::Number lsEnd=lsSegment[1];
      aLSSegment thisSegment;
      thisSegment.first=lsStart.Value();
      thisSegment.second=lsEnd.Value();
      thisRunSegments.push_back(thisSegment);
      //       std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]);
    }
    goodRunLS->insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
  }
  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS->size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS->begin(); itR!=goodRunLS->end(); ++itR){
    std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
    for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
      std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] ";
    std::cout << std::endl;
  }
  return goodRunLS;
}

bool isGoodRunLS( runsLSSegmentsMap&  goodRunLS, int runNumber, int lumiBlock)
{
  //std::cout << "GoodRunLS" << std::endl;
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
        return true;
    }
  return false;
}

int main(int argc, char* argv[]){

  bool processall = true;

  if(argc == 1){
    std::cout << "no lumi file given" << std::endl;
    std::cout << "usage: lumiproc input.root [goodrunls.json]" <<std::endl;
    return 1;
  }
  
  if(argc >3 ){
    std::cout << "too many parameters" << std::endl;
    std::cout << "usage: lumiproc input.root goodrunls.json" <<std::endl;
    return 1;
  }
  
  runsLSSegmentsMap* map=0;// = fillRunLSMap("/afs/cern.ch/user/m/mmozer/scratch0/VecBosApp/config/vecbos/json/dataset_eg_Nov4ReReco.json");


  if(argc == 3 ){ //json file available
    processall = false;
    map = fillRunLSMap(argv[2]);
  }
  
  float totallumi =0;

  TFile* infile = TFile::Open(argv[1]);
  TTree* tree = (TTree*)infile->Get("lumiAna/lumiInfo/lumiTree");

  Int_t           run;
  Int_t           lumi;
  Double_t           lumibyLS;
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("lumiSection", &lumi);
  tree->SetBranchAddress("totalLumiByLS", &lumibyLS);

  int nentriesOrig = tree->GetEntries();
  for(int i=0; i<nentriesOrig; i++) {
    tree->GetEntry(i);
    if(processall || isGoodRunLS(*map,run,lumi)){
      totallumi += lumibyLS;
    }
  }
  
  std::cout << "total luminosity = " << totallumi << " ub" << std::endl;

  return 0;

}
