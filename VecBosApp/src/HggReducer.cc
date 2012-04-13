// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "ReadConfig.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#incclude "VecbosEGObject.hh"
#include "HggReducerSetBranches.cc" // define the setOutputBranches() method here, since it is very long
#include "HggReducer.hh"

HggReducer::HggReducer(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
  preSelSet = 1; // default
}

HggReducer::HggReducer(TTree *tree, string json, bool goodRunLS, bool isData,int mod) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
  preSelSet = 1;
}

HggReducer::~HggReducer() {
  delete vertexer;
  delete selector;
  delete energyScale;
}

void HggReducer::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void HggReducer::SetWeight(double weight) {
  _weight = weight;
}

void HggReducer::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  //do initializations (setup output tree, etc.)
  this->init();

  if(preSelSet < 0 || preSelSet >= preselections.size()){
    cout << "ERROR: invalid preselection set " << preSelSet << endl;
    return;
  }
 
  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // setup for triggers:
  const int nTrigs = triggerNames->size();
  vector<vector<string> > masks
  for(int iTrig=0;iTrigg<nTrigs;iTrig++){
    masks.at(i).push_back(triggerNames.at(i));// setup the trigger masks for Vecbos
  }

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
            // hadronic reload trigger masks and get HLT decisions
      for(int iTrig=0;iTrig<nTrigs;iTrig++){
	setRequiredTriggers(masks.at(i));
	reloadTriggerMask(true);
	triggerBits[i] = hasPassedHLT();
      }
    }

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }    

    std::vector<VecbosPho> photons;
    for(int iPho=0;iPho<nPho;iPho++){
      VecbosPho pho(this,iPho);

      //DO ENERGY CORRECTION
      std::pair<float,float> cor(pho.energy,0.);
      if(correctionType == 99){
	cor = corrector.photonEnergyCorrector_CorrectedEnergyWithError(j);
      }else if(correctionType == 96){
	cor = corrector.photonEnergyCorrector_CorrectedEnergyWithErrorv2(j);
      }
      pho.correctedEnergy = cor.first;
      pho.correctedEnergyError = cor.second;

      pho.dEoE = energyScale.getDEoE(pho,runNumber);
      pho.scaledEnergy = pho.correctedEnergy*(1-pho.dEoE);

      
    }//for(int iPho=0; ...


  }

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  //write output
 file->Close();
}




void HggReducer::init(){
  //read the config file to setup TMVA
  ReadConfig cfg(configFilePath);
  perVtxMvaWeights = cfg.getParameter("perVtxMvaWeights");
  perVtxMvaMethod  = cfg.getParameter("perVtxMvaMethod");
  perEvtMvaWeights = cfg.getParameter("perEvtMvaWeights");
  perEvtMvaMethod  = cfg.getParameter("perEvtMvaMethod");
  
  
  //define the tree
 outTree = new TTree("outTree","outTree");

 //setup the trigger objects
 triggerBits = new int[triggerNames->size()];

 //initialize pointers in the class def
 photonscbclusterenergy = new std::vector<std::vector<float> >; photonscbclusterenergy->clear();
 photonscbclusterposx = new std::vector<std::vector<float> >; photonscbclusterposx->clear();
 photonscbclusterposy = new std::vector<std::vector<float> >; photonscbclusterposy->clear();
 photonscbclusterposz = new std::vector<std::vector<float> >; photonscbclusterposz->clear();
 pileupBunchX = new std::vector<short>; pileupBunchX->clear();
 pileupNInteraction = new std::vector<short>; pileupNInteraction->clear();
 
 this->setOutputBranches(); // initialize the branches of the output tree
 this->setupPreSelection(); // define the photon preselection cuts

 vertexer  = new HggVertexing(this);
 vertexer.setConfigFile("../hgg/default.cfg");
 vertexer.useConversions();
 corrector = new HggEGEnergyCorrector(this,correctionType,_isData);
 if(correctionType == 99){
   energyScale = new HggEnergyScale("../hgg/energy_scale.default");
 }else if(correctionType == 96){
   energyScale = new HggEnergyScale("../hgg/energy_scale.default");
 }
}

void HggReducer::setupPreSelection(){
  //used to define the sets of prescale cuts we can use
  //maybe this should be totally configurable, but for now hard-code
  PreSelCuts set1;
  set1.scet1 = 35;                                                                                                                              
  set1.scet2 = 25;                                                                                                                              
  set1.maxeta = 2.5;                                                                                                                            
  set1.ecaliso_eb  = 10;                                                                                                                        
  set1.ecaliso_ee = 10;                                                                                                                         
  set1.hcaliso_eb = -1;                                                                                                                         
  set1.hcaliso_ee = -1;                                                                                                                         
  set1.sieie_eb = 0.013;                                                                                                                        
  set1.sieie_ee = 0.03;                                                                                                                         
  set1.hoe = 0.15;

  PreSelCuts set2;
  set2.scet1 = 20;                                                                                                                              
  set2.scet2 = 20;                                                                                                                              
  set2.maxeta = 2.5;                                                                                                                            
  set2.ecaliso_eb  = -1;                                                                                                                        
  set2.ecaliso_ee = -1;                                                                                                                         
  set2.hcaliso_eb = -1;                                                                                                                         
  set2.hcaliso_ee = -1;                                                                                                                         
  set2.sieie_eb = 0.013;                                                                                                                        
  set2.sieie_ee = 0.03;                                                                                                                         
  set2.hoe = 0.15; 

  PreSelCuts set3;
  set3.scet1 = 20;                                                                                                                            
  set3.scet2 = 20;                                                                                                                            
  set3.maxeta = -1;                                                                                                                           
  set3.ecaliso_eb  = -1;                                                                                                                      
  set3.ecaliso_ee = -1;                                                                                                                       
  set3.hcaliso_eb = -1;                                                                                                                       
  set3.hcaliso_ee = -1;                                                                                                                       
  set3.sieie_eb = -1;                                                                                                                         
  set3.sieie_ee = -1;                                                                                                                         
  set3.hoe = -1;       

  PreSelCuts set4;
  set4.scet1 = 28;                                                                                                                            
  set4.scet2 = 28;                                                                                                                            
  set4.maxeta = -1;                                                                                                                           
  set4.ecaliso_eb  = -1;                                                                                                                      
  set4.ecaliso_ee = -1;                                                                                                                       
  set4.hcaliso_eb = -1;                                                                                                                       
  set4.hcaliso_ee = -1;                                                                                                                       
  set4.sieie_eb = -1;                                                                                                                         
  set4.sieie_ee = -1;                                                                                                                         
  set4.hoe = -1;        

  PreSelCuts set5;
  set5.scet1 = 25;                                                                                                                            
  set5.scet2 = 25;                                                                                                                            
  set5.maxeta = -1;                                                                                                                           
  set5.ecaliso_eb  = -1;                                                                                                                      
  set5.ecaliso_ee = -1;                                                                                                                       
  set5.hcaliso_eb = -1;                                                                                                                       
  set5.hcaliso_ee = -1;                                                                                                                       
  set5.sieie_eb = -1;                                                                                                                         
  set5.sieie_ee = -1;                                                                                                                         
  set5.hoe = -1;    

  preselections.push_back(set1);
  preselections.push_back(set2);
  preselections.push_back(set3);
  preselections.push_back(set4);
  preselections.push_back(set5);

}


