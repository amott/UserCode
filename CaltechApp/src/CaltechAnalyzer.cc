#include "CaltechAnalyzer.hh"
#include <iostream>

CaltechAnalyzer::CaltechAnalyzer(TTree *tree, TString treeName,TString outputFileName):
  __Photons(0),
  __Muons(0),
  __Electrons(0),
  __Jets(0),
  __CaloMet(0),
  __PFMet(0),
  __TCMet(0),
  __Vertices(0),
  __GenPho(0),
  __GenMuons(0),
  __GenElectrons(0),
  __GenHiggs(0),
  __GenOthers(0),
  __ggVerticesPhotonIndices(0),
  __ggVerticesVertexIndex(0),
  __ggVerticesPerEvtMVA(0)
{
  isValid = true;

  __outputFileName=outputFileName;
  __outputTreeName=treeName;

  __inputTree = tree;
}

void CaltechAnalyzer::process(){
  //main function for processing events
  this->masterBegin();
  if( !this->getIsValid() ) return;
  this->slaveBegin();
  if( !this->getIsValid() ) return;

  Long64_t iEntry=-1;
  while(tree->GetEntry(++iEntry)){//Main Loop

    if(iEntry % 5000 == 0) std::cout << "Procesing Entry " << iEntry << std::endl;

    this->processEvent();
    if( !this->getIsValid() ) return;
    
  }//End Main Loop

  this->slaveEnd();
  if( !this->getIsValid() ) return;
  this->masterEnd();
}

void CaltechAnalyzer::masterBegin(){
  if( !this->getIsValid() ) return;
  this->setBranchAddress();
  
  __outputFile = new TFile(__outputFileName,"RECREATE");
  __outputTree = new TTree(__outputTreeName,"");
  
}

void CaltechAnalyzer::masterEnd(){
  if( !this->getIsValid() ) return;

  __outputTree->Write();
  __outputFile->Close();

  this->masterClear();
}


void CaltechAnalyzer::masterClear(){
  __Photons->clear();
  __Muons->clear();
  __Electrons->clear();
  __Jets->clear();
  __CaloMet->clear();
  __PFMet->clear();
  __TCMet->clear();
  __Vertices->clear();
  __GenPho->clear();
  __GenMuons->clear();
  __GenElectrons->clear();
  __GenHiggs->clear();
  __GenOthers->clear();
  __ggVerticesPhotonIndices->clear();
  __ggVerticesVertexIndex->clear();
  __ggVerticesPerEvtMVA->clear();

  __nPho=0;
  __nMu=0;
  __nEle=0;
  __nJets==0;
  __nCaloMet=0;
  __nPFMet=0;
  __nTCMet=0;
  __nVtx=0;
  __nGenPho=0;
  __nGenMu=0;
  __nGenEle=0;
  __nGenHiggs=0;
  __nGenOthers=0;
  __nPair=0;
}

void CaltechAnalyzer::setConfigFile(TString cfgName){
  if(__config != 0) delete __config;
  __config = new ReadConfig();
  
  if(__config.read(cfgName)!=0){
    this->makeInvalid("error parsing config");
    return;
  }
  // look for explicit collections to disable in the config file
  std::vector<TString> collectionsToDisable = __config.getTokens("disableCollections",",");

  std::vector<TString>::const_iterator colIt;
  for(colIt = collectionsToDisable.begin(); colIt != collectionsToDisable.end(); colIt++)
    this->disableCollection(*colIt);
  
}

void CaltechAnalyzer::setBranchAddress(){
  //sets the branch addresses for all collections defined in the Reduced trees
  //if any collections are not desired, they can be turned off by means of the 
  // disableCollection method

  __inputTree->SetBranchAddress("runNumber",&__runNumber);
  __inputTree->SetBranchAddress("evtNumber",&__evtNumber);
  __inputTree->SetBranchAddress("lumiBlock",&__lumiBlock);

  __inputTree->SetBranchAddress("rho",&__rho);
  __inputTree->SetBranchAddress("rhoEtaMax44",&__rhoEtaMax44);

  __inputTree->SetBranchAddress("nPho",&__nPho);
  __inputTree->SetBranchAddress("Photons",&__Photons);

  __inputTree->SetBranchAddress("nMu",&__nMu);
  __inputTree->SetBranchAddress("Muons",&__Muons);
  
  __inputTree->SetBranchAddress("nEle",&__nEle);
  __inputTree->SetBranchAddress("Electrons",&__Electrons);
  
  __inputTree->SetBranchAddress("nJets",&__nJets);
  __inputTree->SetBranchAddress("Jets",&__Jets);
  
  __inputTree->SetBranchAddress("nCaloMet",&__nCaloMet);
  __inputTree->SetBranchAddress("CaloMet",&__CaloMet);
  
  __inputTree->SetBranchAddress("nPFMet",&__nPFMet);
  __inputTree->SetBranchAddress("PFMet",&__PFMet);
  
  __inputTree->SetBranchAddress("nTCMet",&__nTCMet);
  __inputTree->SetBranchAddress("TCMet",&__TCMet);
  
  __inputTree->SetBranchAddress("nVtx",&__nVtx);
  __inputTree->SetBranchAddress("Vertices",&__Vertices);
  
  __inputTree->SetBranchAddress("nGenPho",&__nGenPho);
  __inputTree->SetBranchAddress("GenPhotons",&__GenPhotons);
  
  __inputTree->SetBranchAddress("nGenMu",&__nGenMu);
  __inputTree->SetBranchAddress("GenMuons",&__GenMuons);
  
  __inputTree->SetBranchAddress("nGenEle",&__nGenEle);
  __inputTree->SetBranchAddress("GenElectrons",&__GenElectrons);
  
  __inputTree->SetBranchAddress("nGenHiggs",&__nGenHiggs);
  __inputTree->SetBranchAddress("GenHiggs",&__GenHiggs);
  
  __inputTree->SetBranchAddress("nGenOthers",&__nGenOthers);
  __inputTree->SetBranchAddress("GenOthers",&__GenOthers);
  
  __inputTree->SetBranchAddress("nPair",&__nPair);
  __inputTree->SetBranchAddress("ggVerticesPhotonIndices",&__ggVerticesPhotonIndices);
  __inputTree->SetBranchAddress("ggVerticesVertexIndex",&__ggVerticesVertexIndex);
  __inputTree->SetBranchAddress("ggVerticesPerEvtMVA",&__ggVerticesPerEvtMVA);
}
