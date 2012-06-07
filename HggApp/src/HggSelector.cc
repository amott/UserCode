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
  doMuMuGamma(false)
{
}

HggSelector::HggSelector(vector<string> fNames, string treeName,string outFName):
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  doElectronVeto(true),
  Photons_(0),
  nSigma(3),
  doMuMuGamma(false)
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
  pho1MVA_  = -999.;
  pho2MVA_  = -999.;

  triggerDec = new int[triggers.size()];

  massRes = new HggMassResolution();
  massRes->setConfigFile(massResConfig);
  massRes->init();
  
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

  PhotonID = new HggPhotonID();
  PhotonID->setConfig(configFile);
  PhotonID->Init();

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
  while(fChain->GetEntry(++jentry)){
    if(jentry%500==0) cout << ">> Processing Entry " << jentry << "/" << nEntries << endl;

    nPU_ = inPU;
    this->clear();
    if(!isData_) this->fillGenInfo();
    if(doMuMuGamma) this->fillMuMuGamma();
    
    if(debugSelector) cout << "requiring triggers ... " << flush;
    trigger_ = this->requireTrigger();
    if(!trigger_){
      outTree->Fill();
      continue; // don't bother doing the MVA 
    }

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
	indices = getBestPairPFCiC(iSmear,0);
	if(indices.first == -1 || indices.second==-1){
	  mPairSmearPFCiC.push_back(-1);
	}else{
	  mPairSmearPFCiC.push_back(getMPair(indices.first,indices.second));
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
	indices = getBestPairPFCiC(-999,iScale);
	if(indices.first == -1 || indices.second==-1){
	  mPairScalePFCiC.push_back(-1);
	}else{
	  mPairScalePFCiC.push_back(getMPair(indices.first,indices.second));
	}

      }

    }
    
    indices = getBestPair(mvaOut,0,0); // no scaling and default smearing
    index1 = indices.first;
    index2 = indices.second;
    diPhoMVA_ = mvaOut[0];
    pho1MVA_  = mvaOut[1];
    pho2MVA_  = mvaOut[2];

    //Do the PFCiC selection as well!
    indices = getBestPairPFCiC(0,0); // no scaling and default smearing
    index1PFCiC = indices.first;
    index2PFCiC = indices.second;

    if(debugSelector) cout << "LOOP DONE" << endl;	  
    

    if(debugSelector) cout << "indices: " << index1 << "  " << index2 << endl;

    if(index1 > -1 && index2 > -1){
      //fill MVA variables
      int selectedVertex = getVertexIndex(index1,index2);
      if(debugSelector) cout << "Final Selection: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1);
      pho2_ = Photons_->at(index2);
      OutPhotons_.push_back(pho1_);
      OutPhotons_.push_back(pho2_);
      energyPho1 = pho1_.finalEnergy;
      energyNoCorrPho1 = pho1_.energy;
      etaPho1    = pho1_.SC.eta;
      pxPho1     = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy).Px();
      pyPho1     = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy).Py();
      pzPho1     = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy).Pz();
      r9Pho1     = pho1_.SC.r9();
      indexPho1  = pho1_.index;

      energyPho2 = pho2_.finalEnergy;
      energyNoCorrPho2 = pho2_.energy;      
      etaPho2    = pho2_.SC.eta;
      pxPho2     = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy).Px();
      pyPho2     = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy).Py();
      pzPho2     = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy).Pz();
      r9Pho2     = pho2_.SC.r9();
      indexPho2  = pho2_.index;
      if(debugSelector) cout << "Photon Indices: " << pho1_.index << "  " << pho2_.index << endl;

      mPair_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorr_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairRes_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,false);
      mPairResWrongVtx_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtx_ = selectedVertex;
      diPhoVtxX_ = vtxX[selectedVertex];
      diPhoVtxY_ = vtxY[selectedVertex];
      diPhoVtxZ_ = vtxZ[selectedVertex];
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
      OutPhotonsPFCiC_.push_back(pho1_);
      OutPhotonsPFCiC_.push_back(pho2_);

      mPairPFCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorrPFCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairResPFCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,false);
      mPairResWrongVtxPFCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtxPFCiC_ = selectedVertex;
      diPhoVtxXPFCiC_ = vtxX[selectedVertex];
      diPhoVtxYPFCiC_ = vtxY[selectedVertex];
      diPhoVtxZPFCiC_ = vtxZ[selectedVertex];
    }else{
      mPairPFCiC_=-1;
    }

    outTree->Fill();
  }//while(fChain...

  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  if(doMuMuGamma) outTreeMuMuG->Write();
  f->Close();
}

float HggSelector::getMPair(int i1, int i2){
  int selectedVertex = getVertexIndex(i1,i2);
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
  
  VecbosPho pho1 = Photons_->at(i1);
  VecbosPho pho2 = Photons_->at(i2);
  
  return (pho1.p4FromVtx(vtxPos,pho1.finalEnergy) + pho2.p4FromVtx(vtxPos,pho2.finalEnergy)).M();
}

std::pair<int,int> HggSelector::getBestPairPFCiC(int smearShift,int scaleShift){
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
	bool CiC1 = PhotonID->getIdCiCPF(pho1,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
	bool CiC2 = PhotonID->getIdCiCPF(pho2,nVtx,rho,TVector3(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]),selVtxI);
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
  mvaOut[0] = -999;
  mvaOut[1] = -999;
  mvaOut[2] = -999;
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
	mvaOut[0] = diPhoMVA;
	mvaOut[1] = mva1;
	mvaOut[2] = mva2;
      }
    }// for(iPho2...
  }// for(iPho1...
  
    return indices;
}

void HggSelector::setDefaults(){
  mPair_=-1;
  mPairNoCorr_=-1;
  indexPho1 = -1;
  indexPho2 = -1;
  genHiggsPt = -1;
  nPU_ = inPU;
}
void HggSelector::clear(){
  mPair_=-1;
  mPairNoCorr_=-1;
  diPhoMVA_=-999;
  pho1MVA_=-999;
  pho2MVA_=-999;
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
  outTree->Branch("Photon1MVA",&pho1MVA_,"Photon1MVA/F");
  outTree->Branch("Photon2MVA",&pho2MVA_,"Photon2MVA/F");
  outTree->Branch("diPhotonVtx",&diPhoVtx_,"diPhotonVtx/I");
  outTree->Branch("diPhotonVtxX",&diPhoVtxX_,"diPhotonVtxX/F");
  outTree->Branch("diPhotonVtxY",&diPhoVtxY_,"diPhotonVtxY/F");
  outTree->Branch("diPhotonVtxZ",&diPhoVtxZ_,"diPhotonVtxZ/F");

  outTree->Branch("mPairPFCiC",&mPairPFCiC_,"mPairPFCiC/F");
  outTree->Branch("mPairNoCorrPFCiC",&mPairNoCorrPFCiC_,"mPairNoCorrPFCiC/F");
  outTree->Branch("mPairResPFCiC",&mPairResPFCiC_,"mPairResPFCiC/F");
  outTree->Branch("mPairResWrongVtxPFCiC",&mPairResWrongVtxPFCiC_,"mPairResWrongVtxPFCiC/F");
  outTree->Branch("diPhotonVtxPFCiC",&diPhoVtxPFCiC_,"diPhotonVtxPFCiC/I");
  outTree->Branch("diPhotonVtxXPFCiC",&diPhoVtxXPFCiC_,"diPhotonVtxXPFCiC/F");
  outTree->Branch("diPhotonVtxYPFCiC",&diPhoVtxYPFCiC_,"diPhotonVtxYPFCiC/F");
  outTree->Branch("diPhotonVtxZPFCiC",&diPhoVtxZPFCiC_,"diPhotonVtxZPFCiC/F");

  
  outTree->Branch("Photon",&OutPhotons_);
  outTree->Branch("PhotonPFCiC",&OutPhotonsPFCiC_);

  outTree->Branch("energyPho1",&energyPho1);
  outTree->Branch("energyNoCorrPho1",&energyNoCorrPho1);
  outTree->Branch("etaPho1",&etaPho1);
  outTree->Branch("pxPho1",&pxPho1);
  outTree->Branch("pyPho1",&pyPho1);
  outTree->Branch("pzPho1",&pzPho1);
  outTree->Branch("r9Pho1",&r9Pho1);
  outTree->Branch("catPho1",&catPho1,"catPho1/I");
  outTree->Branch("passPFCiCPho1",&passPFCiCPho1,"passPFCiCPho1/I");  
  outTree->Branch("indexPho1",&indexPho1);

  outTree->Branch("energyPho2",&energyPho2);
  outTree->Branch("energyNoCorrPho2",&energyNoCorrPho2);
  outTree->Branch("etaPho2",&etaPho2);
  outTree->Branch("pxPho2",&pxPho2);
  outTree->Branch("pyPho2",&pyPho2);
  outTree->Branch("pzPho2",&pzPho2);
  outTree->Branch("r9Pho2",&r9Pho2);
  outTree->Branch("catPho2",&catPho2,"catPho2/I");
  outTree->Branch("passPFCiCPho2",&passPFCiCPho2,"passPFCiCPho2/I");  
  outTree->Branch("indexPho2",&indexPho2);

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
  smearedMassErrByMass = massRes->getMassResolution(&pho1,&pho2,vtxPos,false)/mPair;
  smearedMassErrByMassWrongVtx = massRes->getMassResolution(&pho1,&pho2,vtxPos,true)/mPair;
  if(debugSelector) cout << "Mass Erro resolved. Selected Vertex: " << selectedVertex << endl;
  vtxprob=1.-0.49*(getVertexMVA(indexPho1,indexPho2)+1.0);
  if(debugSelector) cout << "vtx prob: " << vtxprob << endl;
  pho1PtByMass = max(p1.Et(),p2.Et())/mPair;
  pho2PtByMass = min(p1.Et(),p2.Et())/mPair;
  pho1Eta = scLead->eta; 
  pho2Eta = scSubLead->eta;
  cosDPhi = TMath::Cos(DeltaPhi(scLead->phi,scSubLead->phi));
  pho1IdMVA = mvaLead;
  pho2IdMVA = mvaSubLead;

  return diPhotonMVA->EvaluateMVA(methodName_diPho);
}
