// std includes

#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
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
#include "VecbosEGObject.hh"
#include "CaltechReducer.hh"
#include "../src/HggPhysUtils.cc"
#include "assert.h"

#define debugReducer 0

CaltechReducer::CaltechReducer(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  doVertexingMVA=true;
}

CaltechReducer::CaltechReducer(TTree *tree, string json, bool goodRunLS, bool isData,int mod) : Vecbos(tree) {
  _goodRunLS = goodRunLS;
  _isData = isData;
  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
  doVertexingMVA=true;
}

struct less_than_pt_photon{
  inline bool operator() (const VecbosPho& p1, const VecbosPho& p2){ 
    return p1.finalEnergy/cosh(p1.eta) < p2.finalEnergy/cosh(p2.eta);
  }
};

void CaltechReducer::Loop(string outFileName, int start, int stop) {
  if(fChain == 0){
    cout << "fChain Not defined! QUITTING" << endl;
    return;
  }
  this->init(outFileName);    


  //do initializations (setup output tree, etc.) 
  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // setup for triggers:
  const int nTrigs = triggerNames.size();
  for(int iTrig=0;iTrig<nTrigs;iTrig++){
    std::vector<string> tmp;
    tmp.push_back(triggerNames.at(iTrig));// setup the trigger masks for Vecbos
    masks.push_back(tmp);
  }

  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Starting with Entry: " << start << endl;
  cout << "Number of entries = " << stop << endl;
  
  Long64_t jentry = start-1;
  //****************************************************************************
  //           MAIN    LOOP
  //****************************************************************************

  while(fChain->GetEntry(++jentry)){
    if(jentry == stop) break;
    if (jentry%500 == 0) cout << ">>> Processing event # " << jentry << endl;

    //Good Run selection
    //APPLY THE JSON
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }    

    this->clearAll();
    this->fillEventInfo();
    if(_isData) this->fillTriggerInfo();

    if(collectionSwitches["Vertices"]) this->fillVertexInfo();
    if(collectionSwitches["Muons"]) this->fillMuons();
    if(collectionSwitches["Electrons"]) this->fillElectrons();
    if(collectionSwitches["Jets"]) this->fillJets();
    this->fillMets(); //selection done within the method
    if(collectionSwitches["Photons"]) this->fillPhotons();
    if(collectionSwitches["Gen"]) if(!_isData) this->fillGeneratorInfo();

    //vetoes
    if(debugReducer) {
      std::cout << "nVtx: " << nVtx_ << std::endl
		<< "nMu:  " << nMu_  << std::endl
		<< "nEle: " << nEle_ << std::endl
		<< "nJet: " << nJet_ << std::endl
		<< "nPho: " << nPho_ << std::endl
		<< "nMet: " << nCaloMet_ << "  " << nPFMet_ << "  " << nTCMet_ << std::endl;
     
	
    }
    if(nVtx_ < reqObjects["Vertices"]) continue;
    if(nMu_  < reqObjects["Muons"]) continue;
    if(nEle_ < reqObjects["Electrons"]) continue;
    if(nJet_ < reqObjects["Jets"]) continue;
    if(nPho_ < reqObjects["Photons"]) continue;
    if(nCaloMet_ < reqObjects["CaloMET"]) continue;
    if(nPFMet_ < reqObjects["PFMET"]) continue;
    if(nTCMet_ < reqObjects["TCMET"]) continue;


    //do the vertexing TMVA
    if(doVertexingMVA) this->doVertexing();
    
    outTree->Fill();
  } // end of main loop

  cout << "Writing Tree:" << endl;
  //open main file
  TFile *file = new TFile(outFileName.c_str(),"RECREATE");

  outTree->Write();

  //write output
 file->Close();
}

void CaltechReducer::clearAll(){
  //clear collections
  nPho_=0;
  Photons_.clear();
  nMu_=0;
  Muons_.clear();
  nEle_=0;
  Electrons_.clear();
  nJet_=0;
  Jets_.clear();
  CaloMet_.clear();
  PFMet_.clear();
  TCMet_.clear();

  nVtx_=0;
  Vertices_.clear();

  nGenPho_=0;
  nGenMu_=0;
  nGenEle_=0;
  nGenHiggs_=0;

  GenPhotons_.clear();
  GenMuons_.clear();
  GenElectrons_.clear();
  GenHiggs_.clear();
  GenOthers_.clear();

  ggVerticesPhotonIndices.clear();
  ggVerticesVertexIndex.clear();
  ggVerticesPerEvtMVA.clear();

}

void CaltechReducer::init(string outputFileName){
  //define the tree
  outTree = new TTree("CaltechReduce","CaltechReduce");

 //setup the trigger objects
 triggerBits = new int[triggerNames.size()];

 this->setOutputBranches(); // initialize the branches of the output tree

 ReadConfig cfg;
 int errorCode = cfg.read(config);
 if(errorCode){
   cout << "ERROR READING CONFIG FILE!!!! " << endl
	<< "ABORTING" << endl;
   throw 1;
   return;
 }
 // allow triggers to either be a comma-separated list or a list that looks like:
 //Triggers_1=
 //Triggers_2=
 //Triggers_3=
 //etc.  we will search for up to 1000 such entries (that should be enough ...)
 if(cfg.getParameter("Triggers").compare("")!=0){ //comma seperated list
   triggerNames  = cfg.getTokens("Triggers",",");
 }else{
   for(int i=0;i<1000;i++){
     if(cfg.getParameter( Form("Triggers_%d",i) ).compare("")!=0)
       triggerNames.push_back(cfg.getParameter( Form("Triggers_%d",i) ));
   }
 }
 
 std::vector<std::string>::const_iterator trigIt;
 std::cout << "Saving decisions for triggers: " << std::endl;
 for(trigIt = triggerNames.begin(); trigIt != triggerNames.end(); trigIt++){
   std::cout << *trigIt << std::endl;
 }


 //setup the vertexing configuration
 doVertexingMVA = cfg.getParameter("doVertexingMVA").compare("no")!=0;
 
 vertexer  = new HggVertexing(this);
 vertexer->setConfigFile(config);
 vertexer->useConversions();
 vertexer->init();

 phocorrector = new HggEGEnergyCorrector(this,config,_isData);
 elecorrector = new HggEGEnergyCorrector(this,config,_isData);
 elecorrector->useElectronWeights();
 
 this->setupCollections(cfg);
}

void CaltechReducer::setupCollections(ReadConfig &cfg){
  validCollections.push_back("Jets");
  validCollections.push_back("CaloMET");
  validCollections.push_back("PFMET");
  validCollections.push_back("TCMET");
  validCollections.push_back("Photons");
  validCollections.push_back("Electons");
  validCollections.push_back("Muons");
  validCollections.push_back("Vertices");
  validCollections.push_back("Gen");
  std::vector<std::string>::const_iterator colIt;

  std::cout << "Filling Collections:" <<std::endl;
  for(colIt = validCollections.begin(); colIt != validCollections.end(); colIt++){
    collectionSwitches[*colIt] = cfg.getParameter(Form("fill_%s",colIt->c_str())).compare("no")!=0;
    std::string ptS = cfg.getParameter(Form("minPt_%s",colIt->c_str()));
    if(ptS.compare("")!=0) ptThresholds[*colIt] = atof(ptS.c_str());
    else ptThresholds[*colIt] = 0;
    std::string reqS = cfg.getParameter(Form("required_%s",colIt->c_str()));
    if(reqS.compare("")!=0) reqObjects[*colIt] = atof(reqS.c_str());
    else reqObjects[*colIt] = 0;

    std::cout << *colIt << "  Fill: " << collectionSwitches[*colIt]
	      << "   min pT: " << ptThresholds[*colIt]
	      << "   min #:  " << reqObjects[*colIt] << std::endl;
  }
}
void CaltechReducer::fillPhotons(){
  nPho_=0;
  for(int iPho=0;iPho<nPho;iPho++){
    VecbosPho pho(this,iPho);
    if(debugReducer) std::cout << pho.pt <<std::endl;
    if(pho.pt < ptThresholds["Photons"]) continue;
    phocorrector->getPhotonEnergyCorrection(pho);
    pho.finalEnergy = pho.regressionEnergy;
    pho.finalEnergyError = pho.regressionEnergyError; //temporary
    Photons_.push_back(pho);
    nPho_++;
  }
  
}

void CaltechReducer::fillEventInfo(){
  runNumber_=runNumber;
  evtNumber_=eventNumber;
  lumiBlock_=lumiBlock;
  rho_=rhoJetsFastJet;  // rho from kt6PFJets
  nPu = nPU[1];
  
}

void CaltechReducer::fillTriggerInfo(){
  //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
  if(_isData) {
    // hadronic reload trigger masks and get HLT decisions
    for(int iTrig=0;iTrig<triggerNames.size();iTrig++){
      setRequiredTriggers(masks.at(iTrig));
      reloadTriggerMask(true);
      triggerBits[iTrig] = hasPassedHLT();
    }
  }
}

void CaltechReducer::fillMuons(){
  nMu_=0;
  for(int iMuon = 0; iMuon<nMuon;iMuon++){
    VecbosMu mu(this,iMuon);
    if(mu.pt < ptThresholds["Muons"]) continue;
    Muons_.push_back(mu);
    nMu_++;
  }
}
void CaltechReducer::fillElectrons(){
  nEle_=0;
  for(int iEle = 0; iEle<nEle;iEle++){
    VecbosEle ele(this,iEle);
    if(ele.pt < ptThresholds["Electron"]) continue;
    elecorrector->getElectronEnergyCorrection(ele);
    Electrons_.push_back(ele);
    nEle_++;
  }
}

void CaltechReducer::fillJets(){
  nJet_=0;
  for(int iJet = 0; iJet<nAK5PFPUcorrJet;iJet++){
    VecbosJet jet(this,iJet,VecbosJet::PFPUcorr);
    if(jet.pt < ptThresholds["Jets"]) continue;
    Jets_.push_back(jet);
    nJet_++;
  }
}

void CaltechReducer::fillMets(){
  if(collectionSwitches["CaloMET"]){
    VecbosCaloMET cmet(this,0);
    if(cmet.pt >= ptThresholds["CaloMET"]) CaloMet_.push_back(cmet);
  }
  if(collectionSwitches["PFMet"]){  
    VecbosPFMET   pmet(this,0);
    if(pmet.pt >= ptThresholds["PFMET"]) PFMet_.push_back(pmet);
  }
  if(collectionSwitches["TCMet"]){
    VecbosTCMET   tmet(this,0);
    if(tmet.pt >= ptThresholds["TCMET"]) TCMet_.push_back(tmet);
  }
  
  nCaloMet_ = CaloMet_.size();
  nPFMet_   = PFMet_.size();
  nTCMet_   = TCMet_.size();
}

void CaltechReducer::fillVertexInfo(){
  nVtx_ = 0;
  for(int i=0;i<nPV;i++){
    VecbosVtx vtx(this,i);
    Vertices_.push_back(vtx);
    nVtx_++;
  }
}

void CaltechReducer::fillGeneratorInfo(){
  //GENERATOR information
  nGenPho_=0;
  nGenMu_=0;
  nGenEle_=0;
  nGenHiggs_=0;

  for(int iGen=0;iGen<nMc;iGen++){
    VecbosGen part(this,iGen);
    if(part.pt < ptThresholds["Gen"]) continue;

    switch(abs(idMc[iGen])){
    case 25:  //higgs
      GenHiggs_.push_back(part); break;
    case 22: //Photons
      GenPhotons_.push_back(part); break;
    case 13: //Muons
      GenMuons_.push_back(part); break;
    case 11: //Electrons
      GenElectrons_.push_back(part); break;
    Default:
      GenOthers_.push_back(part); break;
    }
    nGenHiggs_ = GenHiggs_.size();
    nGenPho_   = GenPhotons_.size();
    nGenEle_   = GenElectrons_.size();
    nGenMu_    = GenMuons_.size();
    nGenOthers_= GenOthers_.size();
  }    
  procID = 0; //genProcessId;
  qScale;
  
  //match the reco objects to gen objects
  
  PhoCollection::iterator phoIt;
  for(phoIt = Photons_.begin(); phoIt !=Photons_.end(); phoIt++)
    phoIt->doGenMatch(this,22);
  MuCollection::iterator muIt;
  for(muIt = Muons_.begin(); muIt !=Muons_.end(); muIt++)
    muIt->doGenMatch(this,13);
  EleCollection::iterator eleIt;
  for(eleIt = Electrons_.begin(); eleIt !=Electrons_.end(); eleIt++)
    eleIt->doGenMatch(this,11);

}

void CaltechReducer::doVertexing(){
  std::vector<VecbosPho>::iterator iPho1;
  std::vector<VecbosPho>::iterator iPho2;
  nPair_=0;
  for(iPho1 = Photons_.begin(); iPho1 != Photons_.end(); iPho1++){
    //VecbosPho pho1 = *iPho1;
    for(iPho2 = iPho1+1; iPho2 != Photons_.end(); iPho2++){
      if(iPho2==iPho1) continue;
      //VecbosPho pho2 = *iPho2;
      float perEvt=-1;
      vector<pair<int,float> > vtxPair = vertexer->vertex_tmva(&*iPho2,&*iPho2,perEvt);	  
      const int nTop=1;
      pair<int,float> top[nTop];
      for(int i=0;i<nTop;i++) top[i] = pair<int,float>(0,-1);
      vector<pair<int,float> >::const_iterator vtxIt;
      for(vtxIt = vtxPair.begin(); vtxIt != vtxPair.end(); vtxIt++){
	pair<int,float> tmp = *vtxIt;
	for(int i=0;i<nTop;i++){
	  if(tmp.second > top[i].second){
	    swap(tmp,top[i]);
	  }
	}
      }
      ggVerticesPhotonIndices.push_back(pair<int,int>(iPho1->index,iPho2->index) );
      ggVerticesVertexIndex.push_back(top[0]);
      ggVerticesPerEvtMVA.push_back(perEvt);
      nPair_++;
    }
  }
  
}

void CaltechReducer::setOutputBranches(){

  //Event info
  outTree->Branch("lumiBlock",&lumiBlock_,"lumiBlock/I");
  outTree->Branch("runNumber",&runNumber_,"runNumber/I");
  outTree->Branch("evtNumber",&evtNumber_,"evtNumber/I");
  
 //physics declared -- should be set by the JSON, but can't hurt
 
 
  outTree->Branch("rho", &rho_,"rho/F");
  outTree->Branch("rhoEtaMax44", &rhoEtaMax44_,"rhoEtaMax44/F");
 
  //trigger
  for(int i=0;i<triggerNames.size();i++){
    outTree->Branch(triggerNames.at(i).c_str(), &(triggerBits[i]), Form("%s/I",triggerNames.at(i).c_str()) );  // this will produce 1 int per trigger in the output tree
  }
 

 //objects
  if(doVertexingMVA){
    outTree->Branch("nPair",&nPair_);
    outTree->Branch("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
    outTree->Branch("ggVerticesVertexIndex",&ggVerticesVertexIndex);
    outTree->Branch("ggVerticesPerEvtMVA",&ggVerticesPerEvtMVA);
  }
 ///information for the vertex
  outTree->Branch("nPho",&nPho_);
  outTree->Branch("Photons",&Photons_);
  
  outTree->Branch("nVtx",&nVtx_,"nVtx/I"); 
  outTree->Branch("Vertices",&Vertices_);
  
  outTree->Branch("nMu",&nMu_,"nMu/I");
  outTree->Branch("Muons",&Muons_);
  
  outTree->Branch("nEle",&nEle_,"nEle/I");
  outTree->Branch("Electrons",&Electrons_);
  
  outTree->Branch("nJet",&nJet_,"nJet/I");
  outTree->Branch("Jets",&Jets_);
  
  outTree->Branch("nCaloMET",&nCaloMet_,"nCaloMET/I");
  outTree->Branch("CaloMET",&CaloMet_);
  outTree->Branch("nPFMET",&nPFMet_,"nPFMET/I");
  outTree->Branch("PFMET",&PFMet_);
  outTree->Branch("nTCMET",&nTCMet_,"nTCMET/I");
  outTree->Branch("TCMET",&TCMet_);

 //FOR MONTE CARLO:
  if(!_isData){
    //generator level information
    outTree->Branch("procID",&procID,"procID/I");
    outTree->Branch("qScale",&qScale,"qScale/F");
    outTree->Branch("nPU",&nPu,"nPU/F");

    ///gen electron, muon,photon
    outTree->Branch("nGenPho",&nGenPho_,"nGenPho/I");
    outTree->Branch("GenPhotons",&GenPhotons_);

    outTree->Branch("nGenEle",&nGenEle_,"nGenEle/I");
    outTree->Branch("GenElectrons",&GenElectrons_);
    
    outTree->Branch("nGenMu",&nGenMu_,"nGenMu/I");
    outTree->Branch("GenMuons",&GenMuons_);
    //generator level higgs
    outTree->Branch("nGenHiggs",&nGenHiggs_,"nGenHiggs/I");
    outTree->Branch("GenHiggs",&GenHiggs_);
    //other particles
    outTree->Branch("nGenOthers",&nGenOthers_,"nGenOthers/I");
    outTree->Branch("GenOthers",&GenOthers_);
  }

}
