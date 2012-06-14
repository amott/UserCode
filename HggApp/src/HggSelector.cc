#include <HggSelector.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
using namespace std;
using namespace TMVA;

#define debugSelector 0

HggSelector::HggSelector():
  fChain(0),
  valid(false),
  doElectronVeto(true),
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  Photons_(0),
  nSigma(3),
  doMuMuGamma(false),
  isData_(true)
{
}

HggSelector::HggSelector(vector<string> fNames, string treeName,string outFName):
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  doElectronVeto(true),
  Photons_(0),
  nSigma(3),
  doMuMuGamma(false),
  isData_(true)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

HggSelector::~HggSelector(){
  delete fChain;
}


void HggSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int HggSelector::init(){
  if(!valid) return -1;

  mPair_ = -1;
  mPairNoCorr_ = -1;
  diPhoMVA_ = -999.;

  triggerDec = new int[triggers.size()];

  
  this->setBranchAddresses();
  this->setupOutputTree();
  ReadConfig cfg;
  if(cfg.read(configFile)!=0){
    cout << "Error reading configuration file!";
    valid = false;
    return -1;
  }
  weightFile_diPho = cfg.getParameter("weightFile_diPho");

  methodName_diPho  = cfg.getParameter("methodName_diPho");

  triggers = cfg.getTokens("Triggers",",");
  cout << "Parameters: " << endl
       << weightFile_diPho << endl
       << methodName_diPho << endl;

  string MassResConfig = cfg.getParameter("MassResolutionConfig");
  massRes = new HggMassResolution();
  massRes->setConfigFile(MassResConfig);
  massRes->init();


  PhotonID = new HggPhotonID();
  PhotonID->setConfig(configFile);
  PhotonID->Init();

  //setup MVA inputs histograms
  MVAInputs["massRes"] =  new TH1F("massRes","",200,0.,0.2);
  MVAInputs["massResWrongVtx"] = new TH1F("massResWrongVtx","",200,0.,0.2);
  MVAInputs["vtxProb"] = new TH1F("vtxProb","",100,0,1);
  MVAInputs["p1_EtByM"] = new TH1F("p1_EtByM","",200,0,1);
  MVAInputs["p2_EtByM"] = new TH1F("p2_EtByM","",200,0,1);
  MVAInputs["p1_Eta"] = new TH1F("p1_Eta","",120,-3,3);
  MVAInputs["p2_Eta"] = new TH1F("p2_Eta","",120,-3,3);
  MVAInputs["CosDPhi"] = new TH1F("cosDPhi","",200,-1,1);
  MVAInputs["p1_idMVA"] = new TH1F("p1_idMVA","",200,-1,1);
  MVAInputs["p2_idMVA"] = new TH1F("p2_idMVA","",200,-1,1);

  MjjDists["BeforeVBF"] = new TH1F("MjjBeforeVBF","",200,0,2000);
  MjjDists["AfterDEta"] = new TH1F("MjjAfterDEta","",200,0,2000);
  MjjDists["AfterZ"] = new TH1F("MjjAfterZ","",200,0,2000);
  MjjDists["Final"] = new TH1F("Final","",200,0,2000);


  this->setupTMVA();
}

void HggSelector::Loop(){
  if(!valid) return;

  this->init();
  this->setDefaults();
  cout << "Getting Entries ... "  << endl;
  Long64_t nEntries = fChain->GetEntries();
  Long64_t jentry=-1;
  int index1=-1,index2=-1;
  int index1PFCiC=-1,index2PFCiC=-1;
  int index1CiC=-1,index2CiC=-1;
  while(fChain->GetEntry(++jentry)){
    if(jentry%500==0) cout << ">> Processing Entry " << jentry << "/" << nEntries << endl;

    nPU_ = inPU;
    this->clear();
    if(!isData_) this->fillGenInfo();
    if(doMuMuGamma) this->fillMuMuGamma();
    
    if(debugSelector) cout << "requiring triggers ... " << flush;
    trigger_ = this->requireTrigger();
    //    if(!trigger_){
    //  outTree->Fill();
    //  continue; // don't bother doing the MVA 
    //}

    if(debugSelector) cout << "done" << endl;
    if(debugSelector) cout << "# Photons: " << nPho_ << endl;

    if(debugSelector) cout << "Starting Photon Loop" << endl;

    float mvaOut[3];
    std::pair<int,int> indices;
    
    if(nSigma>0 && !isData_){ //do the energy scale and energy resolution systematics
      for(int iSmear=-nSigma;iSmear<=nSigma;iSmear++){
	//Do the Smearing for the MVA analysis
	indices = getBestPair(mvaOut,iSmear,0);
	diPhoMVASmear.push_back(mvaOut[0]);
	pho1MVASmear.push_back(mvaOut[1]);
	pho2MVASmear.push_back(mvaOut[2]);

	if(indices.first == -1 || indices.second==-1){
	  mPairSmear.push_back(-1);
	}else{
	  mPairSmear.push_back(getMPair(indices.first,indices.second));
	}
	//do the smearing for the PFCiC Analysis
	indices = getBestPairCiC(iSmear,0,true);
	if(indices.first == -1 || indices.second==-1){
	  mPairSmearPFCiC.push_back(-1);
	}else{
	  mPairSmearPFCiC.push_back(getMPair(indices.first,indices.second));
	}
	//do the smearing for the CiC Analysis
	indices = getBestPairCiC(iSmear,0,false);
	if(indices.first == -1 || indices.second==-1){
	  mPairSmearCiC.push_back(-1);
	}else{
	  mPairSmearCiC.push_back(getMPair(indices.first,indices.second));
	}
      } // Done with smearing

      for(int iScale=-nSigma;iScale<=nSigma;iScale++){ //do the scaling systematic
	indices = getBestPair(mvaOut,-999,iScale);
	diPhoMVAScale.push_back(mvaOut[0]);
	pho1MVAScale.push_back(mvaOut[1]);
	pho2MVAScale.push_back(mvaOut[2]);

	if(indices.first == -1 || indices.second==-1){
	  mPairScale.push_back(-1);
	}else{
	  mPairScale.push_back(getMPair(indices.first,indices.second));
	}
		//do the smearing for the PFCiC Analysis
	indices = getBestPairCiC(-999,iScale,true);
	if(indices.first == -1 || indices.second==-1){
	  mPairScalePFCiC.push_back(-1);
	}else{
	  mPairScalePFCiC.push_back(getMPair(indices.first,indices.second));
	}

		//do the smearing for the CiC Analysis
	indices = getBestPairCiC(-999,iScale,false);
	if(indices.first == -1 || indices.second==-1){
	  mPairScaleCiC.push_back(-1);
	}else{
	  mPairScaleCiC.push_back(getMPair(indices.first,indices.second));
	}
      }
    }
    
    indices = getBestPair(&diPhoMVA_,0,0); // no scaling and default smearing
    index1 = indices.first;
    index2 = indices.second;

    //Do the PFCiC selection as well!
    indices = getBestPairCiC(0,0,true); // no scaling and default smearing
    index1PFCiC = indices.first;
    index2PFCiC = indices.second;

    //Do the CiC selection as well!
    indices = getBestPairCiC(0,0,false); // no scaling and default smearing
    index1CiC = indices.first;
    index2CiC = indices.second;

    if(debugSelector) cout << "LOOP DONE" << endl;	  
    

    if(debugSelector) cout << "indices: " << index1 << "  " << index2 << endl;
    if(debugSelector) cout << "indicesPFCiC: " << index1PFCiC << "  " << index2PFCiC << endl;
    if(debugSelector) cout << "indicesCiC: " << index1CiC << "  " << index2CiC << endl;

    if(index1 > -1 && index2 > -1){
      //fill MVA variables
      int selectedVertex = getVertexIndex(index1,index2);
      if(debugSelector) cout << "Final Selection: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1);
      pho2_ = Photons_->at(index2);
      OutPhotons_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
      OutPhotons_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));

      if(debugSelector) cout << "Photon Indices: " << pho1_.index << "  " << pho2_.index << endl;

      mPair_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorr_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairRes_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,false);
      mPairResWrongVtx_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtx_ = selectedVertex;
      diPhoVtxX_ = vtxX[selectedVertex];
      diPhoVtxY_ = vtxY[selectedVertex];
      diPhoVtxZ_ = vtxZ[selectedVertex];
      float jpt[2];
      Mjj_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
      ptJet1_ = jpt[0];
      ptJet2_ = jpt[1];
    }else{
      mPair_=-1;      
    }

    if(index1PFCiC > -1 && index2PFCiC > -1){
      //fill PFCiC variables
      int selectedVertex = getVertexIndex(index1PFCiC,index2PFCiC);
      if(debugSelector) cout << "Final Selection: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1);
      pho2_ = Photons_->at(index2);
      OutPhotonsPFCiC_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
      OutPhotonsPFCiC_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));

      mPairPFCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorrPFCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairResPFCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,false);
      mPairResWrongVtxPFCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtxPFCiC_ = selectedVertex;
      diPhoVtxXPFCiC_ = vtxX[selectedVertex];
      diPhoVtxYPFCiC_ = vtxY[selectedVertex];
      diPhoVtxZPFCiC_ = vtxZ[selectedVertex];
      float jpt[2];
      MjjPFCiC_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
      ptJet1PFCiC_ = jpt[0];
      ptJet2PFCiC_ = jpt[1];
    }else{
      mPairPFCiC_=-1;
    }

    if(index1CiC > -1 && index2CiC > -1){
      //fill CiC variables
      int selectedVertex = getVertexIndex(index1CiC,index2CiC);
      if(debugSelector) cout << "Final Selection: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1);
      pho2_ = Photons_->at(index2);
      OutPhotonsCiC_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
      OutPhotonsCiC_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));

      mPairCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorrCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairResCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,false);
      mPairResWrongVtxCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtxCiC_ = selectedVertex;
      diPhoVtxXCiC_ = vtxX[selectedVertex];
      diPhoVtxYCiC_ = vtxY[selectedVertex];
      diPhoVtxZCiC_ = vtxZ[selectedVertex];
      float jpt[2];
      MjjCiC_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
      ptJet1CiC_ = jpt[0];
      ptJet2CiC_ = jpt[1];
    }else{
      mPairCiC_=-1;
    }

    nOutPhotons_ = OutPhotons_.size();
    nOutPhotonsPFCiC_ = OutPhotonsPFCiC_.size();
    nOutPhotonsCiC_ = OutPhotonsCiC_.size();
    outTree->Fill();
  }//while(fChain...

  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  if(doMuMuGamma) outTreeMuMuG->Write();
  std::map<std::string,TH1F*>::iterator it;
  for(it = MVAInputs.begin(); it!=MVAInputs.end(); it++) (*it).second->Write();
  for(it = MjjDists.begin();  it!=MjjDists.end();  it++) (*it).second->Write();
  if(PhotonID->getHists()){
    for(it = PhotonID->getHists()->begin();  it!=PhotonID->getHists()->end();  it++) (*it).second->Write();

  }

  f->Close();
}

float HggSelector::getMPair(int i1, int i2){
  int selectedVertex = getVertexIndex(i1,i2);
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
  
  VecbosPho pho1 = Photons_->at(i1);
  VecbosPho pho2 = Photons_->at(i2);
  
  return (pho1.p4FromVtx(vtxPos,pho1.finalEnergy) + pho2.p4FromVtx(vtxPos,pho2.finalEnergy)).M();
}

std::pair<int,int> HggSelector::getBestPairCiC(int smearShift,int scaleShift,bool usePF=true){
  std::pair<int,int> indices(-1,-1);
  //int bestCat=-1;
  //float bestMass=-99;
  float highestPtSum=0; // is this the best way to select the pairs?
  TRandom3 rng(0);

    for(int iPho1=0; iPho1<nPho_;iPho1++){
      if(photonMatchedElectron[iPho1] && doElectronVeto) return std::pair<int,int>(-1,-1);
      for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
	if(iPho1==iPho2) continue;
	if(photonMatchedElectron[iPho2] && doElectronVeto) return std::pair<int,int>(-1,-1);
	if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
	//scale/smear the energy of the photon
	VecbosPho* pho1 = &(Photons_->at(iPho1));
	VecbosPho* pho2 = &(Photons_->at(iPho2));
	if(!isData_){
	  
	  //apply scale shift	  
	  pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
	  pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;

	  float smear = pho1->dEoE + smearShift*pho1->dEoEErr;
	  if(smear < 0) smear = 0;
	  if(smearShift<-100) smear = 0;
	  float rand=0;
	  if(smear > 0) rand = rng.Gaus(0,smear);
	  if(rand < -1) rand=-1;
	  if(rand > 1E3) rand = 1E3;
	  pho1->finalEnergy = pho1->finalEnergy*(1+rand);
	  if(smear > 0) rand = rng.Gaus(0,smear);
	  if(rand < -1) rand = -1;
	  if(rand > 1E3) rand = 1E3;
	  pho2->finalEnergy = pho2->finalEnergy*(1+rand);
	}
	int selVtxI = this->getVertexIndex(iPho1,iPho2);
	bool CiC1,CiC2;
	if(usePF){
	  CiC1 = PhotonID->getIdCiCPF(pho1,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
	  CiC2 = PhotonID->getIdCiCPF(pho2,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
	}else{
	  CiC1 = PhotonID->getIdCiC(pho1,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
	  CiC2 = PhotonID->getIdCiC(pho2,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
	}
	if(!CiC1 || !CiC2) continue;
	float thisPtSum = pho1->p4FromVtx(TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),pho1->finalEnergy).Pt()
	  + pho2->p4FromVtx(TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),pho2->finalEnergy).Pt();	
	if(thisPtSum > highestPtSum){
	  highestPtSum = thisPtSum;
	  indices.first = iPho1;
	  indices.second = iPho2;
	}
      }// for(iPho2...
    }// for(iPho1...
    return indices;
}

std::pair<int,int> HggSelector::getBestPair(float* mvaOut, int smearShift,int scaleShift){
  std::pair<int,int> indices(-1,-1);
  float diPhoMVAMax=-99;
  TRandom3 rng(0);
  *mvaOut = -999.;
  for(int iPho1=0; iPho1<nPho_;iPho1++){
    if(photonMatchedElectron[iPho1] && doElectronVeto) continue;
    for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
      if(iPho1==iPho2) continue;
      if(photonMatchedElectron[iPho2] && doElectronVeto) continue;
      if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
      //NON-PF block
      //scale/smear the energy of the photon
      VecbosPho* pho1 = &(Photons_->at(iPho1));
      VecbosPho* pho2 = &(Photons_->at(iPho2));
      if(!isData_){
	
	//apply scale shift	  
	pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
	pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;
	
	float smear = pho1->dEoE + smearShift*pho1->dEoEErr;
	if(smear < 0) smear = 0;
	if(smearShift<-100) smear = 0;
	float rand=0;
	if(smear > 0) rand = rng.Gaus(0,smear);
	if(rand < -1) rand=-1;
	if(rand > 1E3) rand = 1E3;
	pho1->finalEnergy = pho1->finalEnergy*(1+rand);
	if(smear > 0) rand = rng.Gaus(0,smear);
	if(rand < -1) rand = -1;
	if(rand > 1E3) rand = 1E3;
	pho2->finalEnergy = pho2->finalEnergy*(1+rand);
      }
      int selVtxI = this->getVertexIndex(iPho1,iPho2);
      float mva1 = PhotonID->getIdMVA(pho1,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
      float mva2 = PhotonID->getIdMVA(pho2,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
      float diPhoMVA =  getDiPhoMVA(iPho1,iPho2,mva1,mva2,false);
      if(debugSelector) cout << "\t\t" << mva1 << "  " << mva2 << "  " << diPhoMVA << endl;
      if(diPhoMVA > diPhoMVAMax && diPhoMVA>=-1){
	indices.first = iPho1;
	indices.second = iPho2;
	diPhoMVAMax = diPhoMVA;
	*mvaOut = diPhoMVA;
      }
    }// for(iPho2...
  }// for(iPho1...
  
    return indices;
}

void HggSelector::setDefaults(){
  mPair_=-1;
  mPairNoCorr_=-1;
  genHiggsPt = -1;
  nPU_ = inPU;
}
void HggSelector::clear(){
  mPair_=-1;
  mPairNoCorr_=-1;
  diPhoMVA_=-999;
  diPhoVtx_= -1;
  diPhoVtxX_=0;
  diPhoVtxY_=0;
  diPhoVtxZ_=0;

  if(debugSelector) std::cout << "Clearing Scale/Smear variables" << std::endl;
  mPairScale.clear();
  pho1MVAScale.clear();
  pho2MVAScale.clear();
  diPhoMVAScale.clear();
  mPairSmear.clear();
  pho1MVASmear.clear();
  pho2MVASmear.clear();
  diPhoMVASmear.clear();

  if(debugSelector) std::cout << "Clearing Output Variables" << std::endl;
  OutPhotons_.clear();
  OutPhotonsPFCiC_.clear();
  OutPhotonsCiC_.clear();

  if(debugSelector) std::cout << "Clearing MMG Output Variables" << std::endl;
  if(doMuMuGamma){
    MMG_Mu1.clear();
    MMG_Mu2.clear();
    MMG_Pho.clear();
  }
}

void HggSelector::fillGenInfo(){
  if(debugSelector) std::cout << "Setting Gen Higgs Info" << std::endl;
  if(nGenHiggs > 0){
    genHiggsPt = GenHiggs->front().pt;
    genHiggsVx = GenHiggs->front().Vx;
    genHiggsVy = GenHiggs->front().Vy;
    genHiggsVz = GenHiggs->front().Vz;
  }else{
    genHiggsPt = 0;
    genHiggsVx = 0;
    genHiggsVy = 0;
    genHiggsVz = 0;
  }

  if(debugSelector) std::cout << "Clearing Output Variables" << std::endl;
  int selGenPho=0;
  GenCollection::const_iterator genPho;
  for(genPho = GenPhotons->begin(); genPho != GenPhotons->end(); genPho++){
    if(genPho->idMother != 25) continue; // not a higgs
    if(selGenPho==0){
      etaGenPho1 = genPho->eta;
      phiGenPho1 = genPho->phi;
      ptGenPho1 =  genPho->pt;
      energyGenPho1 = genPho->energy;
    }
    else if(selGenPho==1){
      etaGenPho2 = genPho->eta;
      phiGenPho2 = genPho->phi;
      ptGenPho2 =  genPho->pt;
      energyGenPho2 = genPho->energy;
    }else{
      cout << "More than 2 photons from a higgs decay!!" << endl;
    }
    selGenPho++;
  }
  
}

ReducedPhotonData HggSelector::getReducedData(VecbosPho* pho,TVector3 selVtx,int selVtxI){
  ReducedPhotonData data;
  TLorentzVector p4 = pho->p4FromVtx(selVtx,pho->finalEnergy,false);
  TLorentzVector p4NoCorr = pho->p4FromVtx(selVtx,pho->energy,false);
  TLorentzVector p4Gen;
  if(pho->genMatch.index!=-1) {
    p4Gen.SetPtEtaPhiE(pho->genMatch.pt,pho->genMatch.eta,pho->genMatch.phi,pho->genMatch.energy);
    data.pt_Gen = p4Gen.Pt(); data.eta_Gen = p4Gen.Eta(); data.phi_Gen = p4Gen.Phi(); data.E_Gen = p4Gen.E();
  }
  else{
    p4Gen.SetPtEtaPhiE(0.,0.,0.,0.);
    data.pt_Gen = 0.; data.eta_Gen = 0.; data.phi_Gen = 0.; data.E_Gen = 0.;
  
  }
  data.pt = p4.Pt(); data.eta = p4.Eta(); data.phi = p4.Phi(); data.E = p4.E();
  data.pt_NoCorr = p4NoCorr.Pt(); data.eta_NoCorr = p4NoCorr.Eta(); data.phi_NoCorr = p4NoCorr.Phi(); data.E_NoCorr = p4NoCorr.E();
  data.index = pho->index;
  data.r9 = pho->SC.r9();
  data.passPFCiC = PhotonID->getIdCiCPF(pho,nVtx,rho,selVtx,selVtxI); 
  data.category = (data.r9 < 0.94)+2*(fabs(data.eta) > 1.48); 
  data.idMVA = PhotonID->getIdMVA(pho,nVtx,rho,selVtx,selVtxI);
  return data;
}

void HggSelector::setupTMVA(){
  diPhotonMVA = new TMVA::Reader( "!Color:!Silent" );; 

  // diphoton                
  diPhotonMVA->AddVariable("masserrsmeared/mass",&smearedMassErrByMass); 
  diPhotonMVA->AddVariable("masserrsmearedwrongvtx/mass",&smearedMassErrByMass);         
  diPhotonMVA->AddVariable("vtxprob",&vtxprob);             
  diPhotonMVA->AddVariable("ph1.pt/mass",&pho1PtByMass);         
  diPhotonMVA->AddVariable("ph2.pt/mass",&pho2PtByMass);         
  diPhotonMVA->AddVariable("ph1.eta",&pho1Eta);             
  diPhotonMVA->AddVariable("ph2.eta",&pho2Eta);             
  diPhotonMVA->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&cosDPhi);            
  diPhotonMVA->AddVariable("ph1.idmva",&pho1IdMVA);
  diPhotonMVA->AddVariable("ph2.idmva",&pho2IdMVA);

  //book MVAs:
  diPhotonMVA->BookMVA(  methodName_diPho, weightFile_diPho);
}

void HggSelector::setBranchAddresses(){
  if(!valid) return;
  //Event info
  fChain->SetBranchAddress("lumiBlock",&lumiBlock);
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("evtNumber",&evtNumber);
  fChain->SetBranchAddress("isRealData",&_isData);
  

 ///information for the vertex
  fChain->SetBranchAddress("nVtx",&nVtx);
  fChain->SetBranchAddress("vtxX",vtxX);
  fChain->SetBranchAddress("vtxY",vtxY);
  fChain->SetBranchAddress("vtxZ",vtxZ);
  fChain->SetBranchAddress("vtxChi2",vtxChi2);
  fChain->SetBranchAddress("vtxNdof",vtxNdof);
  fChain->SetBranchAddress("vtxNormalizedChi2",vtxNormalizedChi2);
  fChain->SetBranchAddress("vtxTrackSize",vtxTrackSize);
  fChain->SetBranchAddress("vtxIsFake",vtxIsFake);
  fChain->SetBranchAddress("vtxIsValid",vtxIsValid);
  
  fChain->SetBranchAddress("rho", &rho);
  fChain->SetBranchAddress("rhoEtaMax44", &rhoEtaMax44);
 
 //objects
 fChain->SetBranchAddress("nPho",&nPho_);
 fChain->SetBranchAddress("Photons",&Photons_);

 fChain->SetBranchAddress("photonMatchedElectron",photonMatchedElectron);
 fChain->SetBranchAddress("nPair",&nPair_); 
 fChain->SetBranchAddress("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
 fChain->SetBranchAddress("ggVerticesVertexIndex",&ggVerticesVertexIndex);

 //fChain->SetBranchAddress("nMu",&nMu_);
 //fChain->SetBranchAddress("Muons",&Muons_);

 fChain->SetBranchAddress("nGenHiggs",&nGenHiggs);
 fChain->SetBranchAddress("GenHiggs",&GenHiggs);

 fChain->SetBranchAddress("nGenPho",&nGenPho);
 fChain->SetBranchAddress("GenPhotons",&GenPhotons);

 fChain->SetBranchAddress("nPU",&inPU);

 fChain->SetBranchAddress("nJets",&nJets);
 fChain->SetBranchAddress("ptJet",ptJets);
 fChain->SetBranchAddress("etaJet",etaJets);
 fChain->SetBranchAddress("phiJet",phiJets);
 fChain->SetBranchAddress("energyJet",energyJets);

 vector<string>::const_iterator trigIt;
 int i=0;
 for(trigIt=triggers.begin();trigIt!=triggers.end();trigIt++,i++){
   fChain->SetBranchAddress(trigIt->c_str(),&(triggerDec[i]));
 }

}

void HggSelector::setupOutputTree(){
  outTree = new TTree("HggOutput","");
  outTree->Branch("trigger",&trigger_,"trigger/I");
  outTree->Branch("mPair",&mPair_,"mPair/F");
  outTree->Branch("mPairNoCorr",&mPairNoCorr_,"mPairNoCorr/F");
  outTree->Branch("mPairRes",&mPairRes_,"mPairRes/F");
  outTree->Branch("mPairResWrongVtx",&mPairResWrongVtx_,"mPairResWrongVtx/F");
  outTree->Branch("diPhotonMVA",&diPhoMVA_,"diPhotonMVA/F");
  outTree->Branch("diPhotonVtx",&diPhoVtx_,"diPhotonVtx/I");
  outTree->Branch("diPhotonVtxX",&diPhoVtxX_,"diPhotonVtxX/F");
  outTree->Branch("diPhotonVtxY",&diPhoVtxY_,"diPhotonVtxY/F");
  outTree->Branch("diPhotonVtxZ",&diPhoVtxZ_,"diPhotonVtxZ/F");
  outTree->Branch("Mjj",&Mjj_,"Mjj");
  outTree->Branch("ptJet1",&ptJet1_,"ptJet1");
  outTree->Branch("ptJet2",&ptJet2_,"ptJet2");

  outTree->Branch("mPairPFCiC",&mPairPFCiC_,"mPairPFCiC/F");
  outTree->Branch("mPairNoCorrPFCiC",&mPairNoCorrPFCiC_,"mPairNoCorrPFCiC/F");
  outTree->Branch("mPairResPFCiC",&mPairResPFCiC_,"mPairResPFCiC/F");
  outTree->Branch("mPairResWrongVtxPFCiC",&mPairResWrongVtxPFCiC_,"mPairResWrongVtxPFCiC/F");
  outTree->Branch("diPhotonVtxPFCiC",&diPhoVtxPFCiC_,"diPhotonVtxPFCiC/I");
  outTree->Branch("diPhotonVtxXPFCiC",&diPhoVtxXPFCiC_,"diPhotonVtxXPFCiC/F");
  outTree->Branch("diPhotonVtxYPFCiC",&diPhoVtxYPFCiC_,"diPhotonVtxYPFCiC/F");
  outTree->Branch("diPhotonVtxZPFCiC",&diPhoVtxZPFCiC_,"diPhotonVtxZPFCiC/F");
  outTree->Branch("MjjPFCiC",&MjjPFCiC_,"MjjPFCiC");
  outTree->Branch("ptJet1PFCiC",&ptJet1PFCiC_,"ptJet1PFCiC");
  outTree->Branch("ptJet2PFCiC",&ptJet2PFCiC_,"ptJet2PFCiC");
  
  outTree->Branch("mPairCiC",&mPairCiC_,"mPairCiC/F");
  outTree->Branch("mPairNoCorrCiC",&mPairNoCorrCiC_,"mPairNoCorrCiC/F");
  outTree->Branch("mPairResCiC",&mPairResCiC_,"mPairResCiC/F");
  outTree->Branch("mPairResWrongVtxCiC",&mPairResWrongVtxCiC_,"mPairResWrongVtxCiC/F");
  outTree->Branch("diPhotonVtxCiC",&diPhoVtxCiC_,"diPhotonVtxCiC/I");
  outTree->Branch("diPhotonVtxXCiC",&diPhoVtxXCiC_,"diPhotonVtxXCiC/F");
  outTree->Branch("diPhotonVtxYCiC",&diPhoVtxYCiC_,"diPhotonVtxYCiC/F");
  outTree->Branch("diPhotonVtxZCiC",&diPhoVtxZCiC_,"diPhotonVtxZCiC/F");
  outTree->Branch("MjjCiC",&MjjCiC_,"MjjCiC");
  outTree->Branch("ptJet1CiC",&ptJet1CiC_,"ptJet1CiC");
  outTree->Branch("ptJet2CiC",&ptJet2CiC_,"ptJet2CiC");
  
  outTree->Branch("nPhoton",nOutPhotons_);
  outTree->Branch("Photon",&OutPhotons_);
  outTree->Branch("nPhotonPFCiC",nOutPhotonsPFCiC_);
  outTree->Branch("PhotonPFCiC",&OutPhotonsPFCiC_);
  outTree->Branch("nPhotonCiC",nOutPhotonsCiC_);
  outTree->Branch("PhotonCiC",&OutPhotonsCiC_);

  outTree->Branch("genHiggsPt",&genHiggsPt);
  outTree->Branch("genHiggsVx",&genHiggsVx);
  outTree->Branch("genHiggsVy",&genHiggsVy);
  outTree->Branch("genHiggsVz",&genHiggsVz);
  

  outTree->Branch("ptGenPho1",&ptGenPho1);
  outTree->Branch("etaGenPho1",&etaGenPho1);
  outTree->Branch("phiGenPho1",&phiGenPho1);
  outTree->Branch("energyGenPho1",&energyGenPho1);

  outTree->Branch("ptGenPho2",&ptGenPho2);
  outTree->Branch("etaGenPho2",&etaGenPho2);
  outTree->Branch("phiGenPho2",&phiGenPho2);
  outTree->Branch("energyGenPho2",&energyGenPho2);

  outTree->Branch("mPairScale",&mPairScale);
  outTree->Branch("pho1MVAScale",&pho1MVAScale);
  outTree->Branch("pho2MVAScale",&pho2MVAScale);
  outTree->Branch("diPhoMVAScale",&diPhoMVAScale);
  outTree->Branch("mPairSmear",&mPairSmear);
  outTree->Branch("pho1MVASmear",&pho1MVASmear);
  outTree->Branch("pho2MVASmear",&pho2MVASmear);
  outTree->Branch("diPhoMVASmear",&diPhoMVASmear);

  outTree->Branch("nPU",&nPU_);
  if(doMuMuGamma){
    outTreeMuMuG = new TTree("MuMuGamma","");
    outTreeMuMuG->Branch("nMuMuG",&nMuMuG);

    outTreeMuMuG->Branch("massMuMuGamma",massMuMuGamma,"massMuMuGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMuRegGamma",massMuMuRegGamma,"massMuMuRegGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMuScaleGamma",massMuMuScaleGamma,"massMuMuScaleGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMuGenGamma",massMuMuGenGamma,"massMuMuGenGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMu",massMuMu,"massMuMu[nMuMuG]");
    outTreeMuMuG->Branch("puWeight",puWeight,"puWeight[nMuMuG]");

    outTreeMuMuG->Branch("Muon1",&MMG_Mu1);
    outTreeMuMuG->Branch("Muon2",&MMG_Mu2);
    outTreeMuMuG->Branch("Photon",&MMG_Pho);
  
    outTreeMuMuG->Branch("isosumoetPho",isosumoetPho,"isosumoetPho[nMuMuG]");
    outTreeMuMuG->Branch("mvaPho",mvaPho,"mvaPho[nMuMuG]");
  }
}

void HggSelector::fillMuMuGamma(){
  nMuMuG=0;
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);
  for(int iMu1=0;iMu1<nMu_;iMu1++){
    VecbosMu mu1 = Muons_->at(iMu1);    
    if(mu1.pt < 10) continue;
    if(!mu1.isGlobalMuon || !mu1.isTrackerMuon) continue;
    if(mu1.nTrackHits <= 10 | mu1.nPixelHits==0) continue;
    if(mu1.trackImpactPar >=0.2) continue;
    if(mu1.trkIso >= 3) continue;
    for(int iMu2=iMu1+1;iMu2<nMu_;iMu2++){
      VecbosMu mu2 = Muons_->at(iMu2);
      if(mu2.pt < 10) continue;
      if(mu1.pt<20 && mu2.pt<20) continue; // 20/10 selection
      if(mu1.charge*mu2.charge >=0) continue; //opposite charge
      if(!mu2.isGlobalMuon || !mu2.isTrackerMuon) continue;
      if(mu2.nTrackHits <= 10 | mu2.nPixelHits==0) continue;
      if(mu2.trackImpactPar >=0.2) continue;
      if(mu2.trkIso >= 3) continue;
      for(int iPho=0; iPho<nPho_;iPho++){
	VecbosPho pho = Photons_->at(iPho);

	//fill the different masses
	TLorentzVector p4Mu1; p4Mu1.SetPtEtaPhiM(mu1.pt,mu1.eta,mu1.phi,0.106);
	TLorentzVector p4Mu2; p4Mu2.SetPtEtaPhiM(mu2.pt,mu2.eta,mu2.phi,0.106);

	massMuMu[nMuMuG] = (p4Mu1+p4Mu2).M();
	massMuMuGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.energy,false)).M();
	massMuMuRegGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.correctedEnergy,false)).M();
	massMuMuScaleGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.scaledEnergy,false)).M();
	if(pho.genMatch.index>=0) massMuMuGenGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.genMatch.energy,false)).M();
	else massMuMuGenGamma[nMuMuG] = -1;

	MMG_Mu1.push_back(mu1);
	MMG_Mu2.push_back(mu2);
	MMG_Pho.push_back(pho);
	float eT = pho.p4FromVtx(vtx,pho.energy,false).Et();
	isosumoetPho[nMuMuG] = (pho.dr03EcalRecHitSumEtCone + pho.dr04HcalTowerSumEtCone + pho.dr03TrkSumPtHollowCone + isoSumConst - rho*rhoFac)/eT;
	mvaPho[nMuMuG] = PhotonID->getIdMVA(&pho,nVtx,rho,TVector3(vtxX[0],vtxY[0],vtxZ[0]),0); 
	
	nMuMuG++;
      }
    }
  }


  if(doMuMuGamma) outTreeMuMuG->Fill();
}

#define debugMjj 1
float HggSelector::getVBFMjj(VecbosPho* pho1, VecbosPho* pho2,TVector3 SelVtx,float *jetPts){
  TLorentzVector p1 = pho1->p4FromVtx(SelVtx,pho1->finalEnergy);
  TLorentzVector p2 = pho2->p4FromVtx(SelVtx,pho2->finalEnergy);
  jetPts[0]=0.; jetPts[1]=0.;
  if( max(p1.Pt(),p2.Pt())/(p1+p2).M() < 60./120. ) return 0;
  if( min(p1.Pt(),p2.Pt()) < 25. ) return 0;
  //if( (p1+p2).M() < 120. ) return 0; 
  if(nJets<2) return 0;
  if(debugMjj) cout << "NJets: " << nJets << endl;
  int i1=-1,i2=-1;
  float maxSumPt=0;
  for(int iJ1=0;iJ1<nJets;iJ1++){
    if(ptJets[iJ1] < 20.) continue;
    for(int iJ2=iJ1+1;iJ2<nJets;iJ2++){
      if(ptJets[iJ2] < 20.) continue;
      //if(ptJets[iJ1] < 30. && ptJets[iJ2] < 30.) continue; // require 1 30 GeV jet
      if(ptJets[iJ1]+ptJets[iJ2] > maxSumPt){
	maxSumPt = ptJets[iJ1]+ptJets[iJ2];
	i1=iJ1; i2=iJ2;
      }
    }
  }
  if(debugMjj) cout << "Selected Jet Indices: " << i1 << " " << i2 << endl;
  if(i1==-1 || i2==-1) return 0; //didn't find 2 30 GeV Jets
  TLorentzVector j1; j1.SetPtEtaPhiE(ptJets[i1],etaJets[i1],phiJets[i1],energyJets[i1]);
  TLorentzVector j2; j2.SetPtEtaPhiE(ptJets[i2],etaJets[i2],phiJets[i2],energyJets[i2]);
  MjjDists["BeforeVBF"]->Fill((j1+j2).M());
  float dEtaJ = fabs(etaJets[i1]-etaJets[i2]);
  if(debugMjj) cout << "dEtaJ: " << dEtaJ << endl;
  if(dEtaJ < 3.) return 0;
  MjjDists["AfterDEta"]->Fill( (j1+j2).M() );
  TLorentzVector ggSystem = p1+p2;
  float Z = ggSystem.Eta() - (etaJets[i1]+etaJets[i2])/2;
  if(debugMjj) cout << "Z: " << Z << endl;
  if(fabs(Z)>2.5) return 0;
  MjjDists["AfterZ"]->Fill( (j1+j2).M() );
  TLorentzVector jjSystem = j1+j2;
  if(debugMjj) cout << "dPhi jj gg: " << fabs(jjSystem.DeltaPhi(ggSystem)) << endl;
  if( fabs(jjSystem.DeltaPhi(ggSystem)) < 2.6 ) return 0;
  MjjDists["Final"]->Fill( (j1+j2).M() );
  jetPts[0] = j1.Pt();
  jetPts[1] = j2.Pt();
  if(debugMjj) cout << "jj Mass: " << jjSystem.M() <<endl;
  return jjSystem.M();
}

bool HggSelector::requireTrigger(){
  if(!_isData) return true; //no triggers on MC

  if(triggers.size()==0) return true; //no trigger selection
  
  for(int i=0;i<triggers.size();i++){
    if(triggerDec[i]) return true;
  }
  return false;
}


int HggSelector::getVertexIndex(int indexPho1,int indexPho2){
  //get the vertex selected by the di photon vertex MVA
  int selectedVertex = -1;
  if(indexPho1 > nPho_ || indexPho2 > nPho_) return -1;
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  if(debugSelector) cout << "Original Indices:" << origIndex1 << "  " << origIndex2 << endl;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      selectedVertex = ggVerticesVertexIndex->at(i).first;
      break;
    }
  }
  return selectedVertex;
}
float HggSelector::getVertexMVA(int indexPho1,int indexPho2){
  //get the vertex selected by the di photon vertex MVA
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  float selectedVertex = -1;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      selectedVertex = ggVerticesVertexIndex->at(i).second;
      break;
    }
  }
  return selectedVertex;
}

float HggSelector::getDiPhoMVA(int indexPho1, int indexPho2, float mva1, float mva2, bool usePFSC){
  if(debugSelector) cout << "getDiPhoMVA" <<endl;  
  int selectedVertex = getVertexIndex(indexPho1,indexPho2);
  if(selectedVertex<0 || selectedVertex >= nVtx){
    cout << "WARNING: Photons " << indexPho1 << " and " << indexPho2 
	 << " have no selected vertex!" << endl
	 << "Skipping this pair" << endl;
    return -9999.;
  }
  if(debugSelector) cout << mva1 << "  " << mva2 << endl;
  if(mva1 < -999 || mva2 < -999) return -9999.;
  if(debugSelector) cout << "getting photons" <<endl;  
  VecbosPho pho1 = Photons_->at(indexPho1);
  VecbosPho pho2 = Photons_->at(indexPho2);
  if(debugSelector) cout << "selected Vertex: " << selectedVertex <<endl;  
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);

  TLorentzVector p1 = pho1.p4FromVtx(vtxPos,pho1.finalEnergy);
  TLorentzVector p2 = pho2.p4FromVtx(vtxPos,pho2.finalEnergy);

  VecbosPho* phoLead = (p1.Et() > p2.Et() ? &pho1 : & pho2); //select the leading photon
  VecbosPho* phoSubLead = (p1.Et() > p2.Et() ? &pho2 : & pho1); //select the leading photon

  VecbosSC *scLead;
  VecbosSC *scSubLead;
  if(usePFSC){
    scLead = &(phoLead->PFSC);
    scSubLead = &(phoSubLead->PFSC);
  }else{
    scLead = &(phoLead->SC);
    scSubLead = &(phoSubLead->SC);
  }

  float mvaLead = (p1.Et() > p2.Et() ? mva1 : mva2);
  float mvaSubLead = (p1.Et() > p2.Et() ? mva2 : mva1);

  float mPair = (p1+p2).M();
  //fill variables
  smearedMassErrByMass = massRes->getMassResolutionEonly(&pho1,&pho2,vtxPos)/mPair;
  smearedMassErrByMassWrongVtx = massRes->getMassResolution(&pho1,&pho2,vtxPos,true)/mPair;
  if(debugSelector) cout << "Mass Erro resolved. Selected Vertex: " << selectedVertex << endl;
  vtxprob=1.-0.49*(getVertexMVA(indexPho1,indexPho2)+1.0);
  if(debugSelector) cout << "vtx prob: " << vtxprob << endl;
  pho1PtByMass = max(p1.Et(),p2.Et())/mPair;
  pho2PtByMass = min(p1.Et(),p2.Et())/mPair;
  pho1Eta = (p1.Et() > p2.Et() ? p1.Eta(): p2.Eta());
  pho2Eta = (p1.Et() > p2.Et() ? p2.Eta(): p1.Eta());
  cosDPhi = TMath::Cos(DeltaPhi(scLead->phi,scSubLead->phi));
  pho1IdMVA = mvaLead;
  pho2IdMVA = mvaSubLead;

  MVAInputs["massRes"]->Fill(smearedMassErrByMass);
  MVAInputs["massResWrongVtx"]->Fill(smearedMassErrByMassWrongVtx);
  MVAInputs["vtxProb"]->Fill(vtxprob);
  MVAInputs["p1_EtByM"]->Fill(pho1PtByMass);
  MVAInputs["p2_EtByM"]->Fill(pho2PtByMass);
  MVAInputs["p1_Eta"]->Fill(pho1Eta);
  MVAInputs["p2_Eta"]->Fill(pho2Eta);
  MVAInputs["CosDPhi"]->Fill(cosDPhi);
  MVAInputs["p1_idMVA"]->Fill(pho1IdMVA);
  MVAInputs["p2_idMVA"]->Fill(pho2IdMVA);

  return diPhotonMVA->EvaluateMVA(methodName_diPho);
}
