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
  weightFile_IdEB = cfg.getParameter("weightFile_IdEB");
  weightFile_IdEE = cfg.getParameter("weightFile_IdEE");
  weightFile_diPho = cfg.getParameter("weightFile_diPho");

  methodName_Id  = cfg.getParameter("methodName_Id");
  methodName_diPho  = cfg.getParameter("methodName_diPho");

  triggers = cfg.getTokens("Triggers",",");
  cout << "Parameters: " << endl
       << weightFile_IdEB << endl
       << weightFile_IdEB << endl
       << weightFile_diPho << endl
       << methodName_Id << endl
       << methodName_diPho << endl;

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
  int PFindex1=-1,PFindex2=-1;
  while(fChain->GetEntry(++jentry)){
    if(jentry%500==0) cout << ">> Processing Entry " << jentry << "/" << nEntries << endl;

    nPU_ = inPU;
    this->clear();
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
	indices = getBestPair(mvaOut,iSmear,0);
	diPhoMVASmear.push_back(mvaOut[0]);
	pho1MVASmear.push_back(mvaOut[1]);
	pho2MVASmear.push_back(mvaOut[2]);

	if(indices.first == -1 || indices.second==-1){
	  mPairSmear.push_back(-1);
	}else{
	  mPairSmear.push_back(getMPair(indices.first,indices.second));
	}
      }

      for(int iScale=-nSigma;iScale<=nSigma;iScale++){
	indices = getBestPair(mvaOut,-999,iScale);
	diPhoMVAScale.push_back(mvaOut[0]);
	pho1MVAScale.push_back(mvaOut[1]);
	pho2MVAScale.push_back(mvaOut[2]);

	if(indices.first == -1 || indices.second==-1){
	  mPairScale.push_back(-1);
	}else{
	  mPairScale.push_back(getMPair(indices.first,indices.second));
	}
      }

    }
    
    indices = getBestPair(mvaOut,0,0); // no scaling and default smearing
    index1 = indices.first;
    index2 = indices.second;
    diPhoMVA_ = mvaOut[0];
    pho1MVA_  = mvaOut[1];
    pho2MVA_  = mvaOut[2];


    if(debugSelector) cout << "LOOP DONE" << endl;	  
    if(index1 == -1 || index2==-1 || index1 >=nPho_ || index2 >=nPho_){
      //cout << "No selected photon pair!!!" << endl;
      outTree->Fill();
      continue;
    }

    if(debugSelector) cout << "indices: " << index1 << "  " << index2 << endl;
    int selectedVertex = getVertexIndex(index1,index2);
    if(debugSelector) cout << "Final Selection: " << selectedVertex << endl;
    
    TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
    if(debugSelector) cout << "Vtx Pos: " << vtxPos.X() << "  " << vtxPos.Y() << "  " << vtxPos.Z() << endl;
    if(index1 > -1 && index2 > -1){
      pho1_ = Photons_->at(index1);
      pho2_ = Photons_->at(index2);
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

    if(debugSelector) cout << "indices: " << PFindex1 << "  " << PFindex2 << endl;
    if(PFindex1 != -1 && PFindex2!=-1 && PFindex1 <nPho_ && PFindex2 < nPho_){
      int PFselectedVertex = getVertexIndex(PFindex1,PFindex2);
      if(debugSelector) cout << "PF Final Selection: " << PFselectedVertex << endl;
      TVector3 PFvtxPos(vtxX[PFselectedVertex],vtxY[PFselectedVertex],vtxZ[PFselectedVertex]);
      if(PFindex1 > -1 && PFindex2 > -1){
	pfpho1_ = Photons_->at(PFindex1);
	pfpho2_ = Photons_->at(PFindex2);
	energyPFPho1 = pfpho1_.finalEnergy;
	etaPFPho1    = pfpho1_.SC.eta;
	pxPFPho1     = pfpho1_.p4FromVtx(vtxPos,pfpho1_.finalEnergy,true).Px();
	pyPFPho1     = pfpho1_.p4FromVtx(vtxPos,pfpho1_.finalEnergy,true).Py();
	pzPFPho1     = pfpho1_.p4FromVtx(vtxPos,pfpho1_.finalEnergy,true).Pz();
	r9PFPho1     = pfpho1_.SC.r9();
	indexPFPho1  = pfpho1_.index;
	
	energyPFPho2 = pfpho2_.finalEnergy;
	etaPFPho2    = pfpho2_.SC.eta;
	pxPFPho2     = pfpho2_.p4FromVtx(vtxPos,pfpho2_.finalEnergy).Px();
	pyPFPho2     = pfpho2_.p4FromVtx(vtxPos,pfpho2_.finalEnergy).Py();
	pzPFPho2     = pfpho2_.p4FromVtx(vtxPos,pfpho2_.finalEnergy).Pz();
	r9PFPho2     = pfpho2_.SC.r9();
	indexPFPho2  = pfpho2_.index;

	pfmPair_ = (pfpho1_.p4FromVtx(PFvtxPos,pfpho1_.finalEnergy,true) + pfpho2_.p4FromVtx(PFvtxPos,pfpho2_.finalEnergy,true)).M();
	pfmPairRes_ = massRes->getMassResolution(&pfpho1_,&pfpho2_,vtxPos,false);
	pfmPairResWrongVtx_ = massRes->getMassResolution(&pfpho1_,&pfpho2_,vtxPos,true);
	pfdiPhoVtx_ = PFselectedVertex;
      }else{
	pfmPair_=-1;
      }
    }
    outTree->Fill();
  }//while(fChain...

  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  outTreeMuMuG->Write();
  f->Close();
}

float HggSelector::getMPair(int i1, int i2){
  int selectedVertex = getVertexIndex(i1,i2);
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
  
  VecbosPho pho1 = Photons_->at(i1);
  VecbosPho pho2 = Photons_->at(i2);
  
  return (pho1.p4FromVtx(vtxPos,pho1.finalEnergy) + pho2.p4FromVtx(vtxPos,pho2.finalEnergy)).M();
}

std::pair<int,int> HggSelector::getBestPair(float* mvaOut, int smearShift,int scaleShift){
  std::pair<int,int> indices(-1,-1);
  float diPhoMVAMax=-99;
  TRandom3 rng(0);

    for(int iPho1=0; iPho1<nPho_;iPho1++){
      for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
	if(iPho1==iPho2) continue;
	if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
	//NON-PF block
	//scale/smear the energy of the photon
	if(!isData_){
	  VecbosPho* pho1 = &(Photons_->at(iPho1));
	  VecbosPho* pho2 = &(Photons_->at(iPho2));
	  
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
	float mva1 = getIdMVA(iPho1,iPho2,true,false,-1);
	float mva2 = getIdMVA(iPho1,iPho2,false,false,-1);
	float diPhoMVA =  getDiPhoMVA(iPho1,iPho2,mva1,mva2,false);
	if(debugSelector) cout << "\t\t" << mva1 << "  " << mva2 << "  " << diPhoMVA << endl;
	if(diPhoMVA > diPhoMVAMax){
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
  pfmPair_=-1;
  indexPho1 = -1;
  indexPho2 = -1;
  indexPFPho1 = -1;
  indexPFPho2 = -1;      
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

  pfmPair_=-1;
  pfdiPhoMVA_=-999;
  pfpho1MVA_=-999;
  pfpho2MVA_=-999;  
  pfdiPhoVtx_=-1;

  genHiggsPt = higgsPtIn;
  genHiggsVx = higgsVxIn;
  genHiggsVy = higgsVyIn;
  genHiggsVz = higgsVzIn;

  mPairScale.clear();
  pho1MVAScale.clear();
  pho2MVAScale.clear();
  diPhoMVAScale.clear();
  mPairSmear.clear();
  pho1MVASmear.clear();
  pho2MVASmear.clear();
  diPhoMVASmear.clear();

  int selGenPho=0;
  for(int iPho=0;iPho<nGenPho;iPho++){
    if(pidMomGenPho[iPho] == 25){ //higgs
      if(selGenPho==0){
	etaGenPho1 = etaGenPho[iPho];
	phiGenPho1 = phiGenPho[iPho];
	ptGenPho1 = ptGenPho[iPho];
	energyGenPho1 = energyGenPho[iPho];
      }
      else if(selGenPho==1){
	etaGenPho2 = etaGenPho[iPho];
	phiGenPho2 = phiGenPho[iPho];
	ptGenPho2 = ptGenPho[iPho];
	energyGenPho2 = energyGenPho[iPho];
      }else{
	cout << "More than 2 photons from a higgs decay!!" << endl;
      }
      selGenPho++;
    }
  }
}

void HggSelector::setupTMVA(){
  photonMVA_EB = new TMVA::Reader( "!Color:!Silent" );;
  photonMVA_EE = new TMVA::Reader( "!Color:!Silent" );;
  diPhotonMVA = new TMVA::Reader( "!Color:!Silent" );; 

  //setup the ID MVAs:
  ///EB
  photonMVA_EB->AddVariable("HoE",&hoe);         
  photonMVA_EB->AddVariable("covIEtaIEta",&sigietaieta);         
  photonMVA_EB->AddVariable("tIso1abs",&isosumoet);              
  photonMVA_EB->AddVariable("tIso3abs",&trkisooet);              
  photonMVA_EB->AddVariable("tIso2abs",&isosumoetbad);           
  photonMVA_EB->AddVariable("R9",&r9);           
  photonMVA_EB->AddVariable("absIsoEcal",&ecalisodr03);          
  photonMVA_EB->AddVariable("absIsoHcal",&hcalisodr04);          
  photonMVA_EB->AddVariable("NVertexes",&nVertexf);              
  photonMVA_EB->AddVariable("ScEta",&etasc);     
  photonMVA_EB->AddVariable("EtaWidth",&scetawidth);             
  photonMVA_EB->AddVariable("PhiWidth",&scphiwidth);             
       
  ///EE
  photonMVA_EE->AddVariable("HoE",&hoe);         
  photonMVA_EE->AddVariable("covIEtaIEta",&sigietaieta);         
  photonMVA_EE->AddVariable("tIso1abs",&isosumoet);              
  photonMVA_EE->AddVariable("tIso3abs",&trkisooet);              
  photonMVA_EE->AddVariable("tIso2abs",&isosumoetbad);           
  photonMVA_EE->AddVariable("R9",&r9);           
  photonMVA_EE->AddVariable("absIsoEcal",&ecalisodr03);          
  photonMVA_EE->AddVariable("absIsoHcal",&hcalisodr04);          
  photonMVA_EE->AddVariable("NVertexes",&nVertexf);              
  photonMVA_EE->AddVariable("ScEta",&etasc);     
  photonMVA_EE->AddVariable("EtaWidth",&scetawidth);             
  photonMVA_EE->AddVariable("PhiWidth",&scphiwidth);             
  
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
  photonMVA_EB->BookMVA( methodName_Id, weightFile_IdEB);
  photonMVA_EE->BookMVA( methodName_Id, weightFile_IdEE);
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

 fChain->SetBranchAddress("nMu",&nMu_);
 fChain->SetBranchAddress("Muons",&Muons_);

 fChain->SetBranchAddress("higgsPt",&higgsPtIn);
 fChain->SetBranchAddress("higgsVx",&higgsVxIn);
 fChain->SetBranchAddress("higgsVy",&higgsVyIn);
 fChain->SetBranchAddress("higgsVz",&higgsVzIn);

 fChain->SetBranchAddress("nGenPho",&nGenPho);
 fChain->SetBranchAddress("etaGenPho",etaGenPho);
 fChain->SetBranchAddress("phiGenPho",phiGenPho);
 fChain->SetBranchAddress("ptGenPho",ptGenPho);
 fChain->SetBranchAddress("energyGenPho",energyGenPho);
 fChain->SetBranchAddress("pidMomGenPho",pidMomGenPho);

 fChain->SetBranchAddress("nGenMu",&nGenMu);
 fChain->SetBranchAddress("etaGenMu",etaGenMu);
 fChain->SetBranchAddress("phiGenMu",phiGenMu);
 fChain->SetBranchAddress("ptGenMu",ptGenMu);
 fChain->SetBranchAddress("energyGenMu",energyGenMu);
 fChain->SetBranchAddress("pidMomGenMu",pidMomGenMu);

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

  outTree->Branch("PFmPair",&pfmPair_,"PFmPair/F");
  outTree->Branch("PFmPairRes",&pfmPairRes_,"PFmPairRes/F");
  outTree->Branch("PFmPairResWrongVtx",&pfmPairResWrongVtx_,"PFmPairResWrongVtx/F");
  outTree->Branch("PFdiPhotonMVA",&pfdiPhoMVA_,"PFdiPhotonMVA/F");
  outTree->Branch("PFPhoton1MVA",&pfpho1MVA_,"PFPhoton1MVA/F");
  outTree->Branch("PFPhoton2MVA",&pfpho2MVA_,"PFPhoton2MVA/F");
  outTree->Branch("PFdiPhotonVtx",&pfdiPhoVtx_,"PFdiPhotonVtx/I");

  outTree->Branch("energyPho1",&energyPho1);
  outTree->Branch("energyNoCorrPho1",&energyNoCorrPho1);
  outTree->Branch("etaPho1",&etaPho1);
  outTree->Branch("pxPho1",&pxPho1);
  outTree->Branch("pyPho1",&pyPho1);
  outTree->Branch("pzPho1",&pzPho1);
  outTree->Branch("r9Pho1",&r9Pho1);
  outTree->Branch("indexPho1",&indexPho1);

  outTree->Branch("energyPho2",&energyPho2);
  outTree->Branch("energyNoCorrPho2",&energyNoCorrPho2);
  outTree->Branch("etaPho2",&etaPho2);
  outTree->Branch("pxPho2",&pxPho2);
  outTree->Branch("pyPho2",&pyPho2);
  outTree->Branch("pzPho2",&pzPho2);
  outTree->Branch("r9Pho2",&r9Pho2);
  outTree->Branch("indexPho2",&indexPho2);

  outTree->Branch("energyPFPho1",&energyPFPho1);
  outTree->Branch("etaPFPho1",&etaPFPho1);
  outTree->Branch("pxPFPho1",&pxPFPho1);
  outTree->Branch("pyPFPho1",&pyPFPho1);
  outTree->Branch("pzPFPho1",&pzPFPho1);
  outTree->Branch("r9PFPho1",&r9PFPho1);
  outTree->Branch("indexPFPho1",&indexPFPho1);

  outTree->Branch("energyPFPho2",&energyPFPho2);
  outTree->Branch("etaPFPho2",&etaPFPho2);
  outTree->Branch("pxPFPho2",&pxPFPho2);
  outTree->Branch("pyPFPho2",&pyPFPho2);
  outTree->Branch("pzPFPho2",&pzPFPho2);
  outTree->Branch("r9PFPho2",&r9PFPho2);
  outTree->Branch("indexPFPho2",&indexPFPho2);

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
  
    outTreeMuMuG->Branch("ptMu1",ptMu1,"ptMu1[nMuMuG]");
    outTreeMuMuG->Branch("etaMu1",etaMu1,"etaMu1[nMuMuG]");
    outTreeMuMuG->Branch("phiMu1",phiMu1,"phiMu1[nMuMuG]");
    outTreeMuMuG->Branch("chargeMu1",chargeMu1,"chargeMu1[nMuMuG]/I");
    outTreeMuMuG->Branch("isoMu1",isoMu1,"isoMu1[nMuMuG]");
    outTreeMuMuG->Branch("drPhoMu1",drPhoMu1,"drPhoMu1[nMuMuG]");
  //gen muon, same charge, stat 1
  //dR < 0.5, ptRelDif < 0.5
  //two gen objects cannot match reco object
    outTreeMuMuG->Branch("genMatchMu1",genMatchMu1,"genMatchMu1[nMuMuG]/I");
    outTreeMuMuG->Branch("genPtMu1",genPtMu1,"genPtMu1[nMuMuG]");
    outTreeMuMuG->Branch("genIdMu1",genIdMu1,"genIdMu1[nMuMuG]/I");
    outTreeMuMuG->Branch("genStatusMu1",genStatusMu1,"genStatusMu1[nMuMuG]/I");
    outTreeMuMuG->Branch("genIdMomMu1",genIdMomMu1,"genIdMomMu1[nMuMuG]/I");
    outTreeMuMuG->Branch("genStatusMomMu1",genStatusMomMu1,"genStatusMomMu1[nMuMuG]/I");

    outTreeMuMuG->Branch("ptMu2",ptMu2,"ptMu2[nMuMuG]");
    outTreeMuMuG->Branch("etaMu2",etaMu2,"etaMu2[nMuMuG]");
    outTreeMuMuG->Branch("phiMu2",phiMu2,"phiMu2[nMuMuG]");
    outTreeMuMuG->Branch("chargeMu2",chargeMu2,"chargeMu2[nMuMuG]/I");
    outTreeMuMuG->Branch("isoMu2",isoMu2,"isoMu2[nMuMuG]");
    outTreeMuMuG->Branch("drPhoMu2",drPhoMu2,"drPhoMu2[nMuMuG]");
  //gen muon, same charge, stat 1
  //dR < 0.5, ptRelDif < 0.5
  //two gen objects cannot match reco object
    outTreeMuMuG->Branch("genMatchMu2",genMatchMu2,"genMatchMu2[nMuMuG]/I");
    outTreeMuMuG->Branch("genPtMu2",genPtMu2,"genPtMu2[nMuMuG]");
    outTreeMuMuG->Branch("genIdMu2",genIdMu2,"genIdMu2[nMuMuG]/I");
    outTreeMuMuG->Branch("genStatusMu2",genStatusMu2,"genStatusMu2[nMuMuG]/I");
    outTreeMuMuG->Branch("genIdMomMu2",genIdMomMu2,"genIdMomMu2[nMuMuG]/I");
    outTreeMuMuG->Branch("genStatusMomMu2",genStatusMomMu2,"genStatusMomMu2[nMuMuG]/I");

  
    outTreeMuMuG->Branch("defEnergyPho",defEnergyPho,"defEnergyPho[nMuMuG]");
    outTreeMuMuG->Branch("regEnergyPho",regEnergyPho,"regEnergyPho[nMuMuG]");
    outTreeMuMuG->Branch("scaleEnergyPho",scaleEnergyPho,"scaleEnergyPho[nMuMuG]");
    outTreeMuMuG->Branch("etaPho",etaPho,"etaPho[nMuMuG]");
    outTreeMuMuG->Branch("phiPho",phiPho,"phiPho[nMuMuG]");
    outTreeMuMuG->Branch("eleMatchPho",eleMatchPho,"eleMatchPho[nMuMuG]/I");
    outTreeMuMuG->Branch("r9Pho",r9Pho,"r9Pho[nMuMuG]");
    outTreeMuMuG->Branch("hOePho",hOePho,"hOePho[nMuMuG]");
    outTreeMuMuG->Branch("dr03EcalIsoPho",dr03EcalIsoPho,"dr03EcalIsoPho[nMuMuG]");
    outTreeMuMuG->Branch("dr04HcalIsoPho",dr04HcalIsoPho,"dr04HcalIsoPho[nMuMuG]");
    outTreeMuMuG->Branch("isosumoetPho",isosumoetPho,"isosumoetPho[nMuMuG]");
    outTreeMuMuG->Branch("dr03TrkSumHollowConePho",dr03TrkSumHollowConePho,"dr03TrkSumHollowConePho[nMuMuG]");  
    outTreeMuMuG->Branch("mvaPho",mvaPho,"mvaPho[nMuMuG]");
  // photon matching: dR < 0.2
  // |recoPt-truePt/truePt| < 1

    outTreeMuMuG->Branch("genMatchPho",genMatchPho,"genMatchPho[nMuMuG]/I");
    outTreeMuMuG->Branch("genEnergyPho",genEnergyPho,"genEnergyPho[nMuMuG]");
    outTreeMuMuG->Branch("genIdPho",genIdPho,"genIdPho[nMuMuG]/I");
    outTreeMuMuG->Branch("genStatusPho",genStatusPho,"genStatusPho[nMuMuG]/I");
    outTreeMuMuG->Branch("genIdMomPho",genIdMomPho,"genIdMomPho[nMuMuG]/I");
    outTreeMuMuG->Branch("genStatusMomPho",genStatusMomPho,"genStatusMomPho[nMuMuG]/I");
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
    int genIndexMu1 = getGenMatchMu(&mu1);
    for(int iMu2=iMu1+1;iMu2<nMu_;iMu2++){
      VecbosMu mu2 = Muons_->at(iMu2);
      if(mu2.pt < 10) continue;
      if(mu1.pt<20 && mu2.pt<20) continue; // 20/10 selection
      if(mu1.charge*mu2.charge >=0) continue; //opposite charge
      if(!mu2.isGlobalMuon || !mu2.isTrackerMuon) continue;
      if(mu2.nTrackHits <= 10 | mu2.nPixelHits==0) continue;
      if(mu2.trackImpactPar >=0.2) continue;
      if(mu2.trkIso >= 3) continue;
      int genIndexMu2 = getGenMatchMu(&mu2);
      for(int iPho=0; iPho<nPho_;iPho++){
	VecbosPho pho = Photons_->at(iPho);
	int genIndexPho = getGenMatchPho(&pho);

	//fill the different masses
	TLorentzVector p4Mu1; p4Mu1.SetPtEtaPhiM(mu1.pt,mu1.eta,mu1.phi,0.106);
	TLorentzVector p4Mu2; p4Mu2.SetPtEtaPhiM(mu2.pt,mu2.eta,mu2.phi,0.106);

	massMuMu[nMuMuG] = (p4Mu1+p4Mu2).M();
	massMuMuGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.energy,false)).M();
	massMuMuRegGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.correctedEnergy,false)).M();
	massMuMuScaleGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.scaledEnergy,false)).M();
	if(genIndexPho>=0) massMuMuGenGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,energyGenPho[genIndexPho],false)).M();
	else massMuMuGenGamma[nMuMuG] = -1;

	//muon info
	ptMu1[nMuMuG] = mu1.pt;
	etaMu1[nMuMuG] = mu1.eta;
	phiMu1[nMuMuG] = mu1.phi;
	chargeMu1[nMuMuG] = mu1.charge;
	isoMu1[nMuMuG] = mu1.combinedIso;
	drPhoMu1[nMuMuG] = DeltaR(mu1.eta,pho.SC.eta,mu2.phi,pho.SC.phi);
	genMatchMu1[nMuMuG] = (genIndexMu1>=0);
	if(genMatchMu1[nMuMuG]){
	  genPtMu1[nMuMuG] = ptGenMu[genIndexMu1];
	  genIdMu1[nMuMuG] = mu1.charge*(-13);
	  genStatusMu1[nMuMuG] = 1;
	  genIdMomMu1[nMuMuG] = pidMomGenMu[genIndexMu1];
	  genStatusMomMu1[nMuMuG] = 3;
	}else{
	  genPtMu1[nMuMuG] = 0;
	  genIdMu1[nMuMuG] = 0;
	  genStatusMu1[nMuMuG] = 0;
	  genIdMomMu1[nMuMuG] = 0;
	  genStatusMomMu1[nMuMuG] = 0;
	}
	
	ptMu2[nMuMuG] = mu2.pt;
	etaMu2[nMuMuG] = mu2.eta;
	phiMu2[nMuMuG] = mu2.phi;
	chargeMu2[nMuMuG] = mu2.charge;
	isoMu2[nMuMuG] = mu2.combinedIso;
	drPhoMu2[nMuMuG] = DeltaR(mu2.eta,pho.SC.eta,mu2.phi,pho.SC.phi);
	genMatchMu2[nMuMuG] = (genIndexMu2>=0);
	if(genMatchMu2[nMuMuG]){
	  genPtMu2[nMuMuG] = ptGenMu[genIndexMu2];
	  genIdMu2[nMuMuG] = mu2.charge*(-13);
	  genStatusMu2[nMuMuG] = 1;
	  genIdMomMu2[nMuMuG] = pidMomGenMu[genIndexMu2];
	  genStatusMomMu2[nMuMuG] = 3;
	}else{
	  genPtMu2[nMuMuG] = 0;
	  genIdMu2[nMuMuG] = 0;
	  genStatusMu2[nMuMuG] = 0;
	  genIdMomMu2[nMuMuG] = 0;
	  genStatusMomMu2[nMuMuG] = 0;
	}

	//photon
	defEnergyPho[nMuMuG] = pho.energy;
	regEnergyPho[nMuMuG] = pho.correctedEnergy;
	scaleEnergyPho[nMuMuG] = pho.scaledEnergy;
	etaPho[nMuMuG] = pho.SC.eta;
	phiPho[nMuMuG] = pho.SC.phi;
	eleMatchPho[nMuMuG] = photonMatchedElectron[iPho];
	r9Pho[nMuMuG] = pho.SC.r9();
	hOePho[nMuMuG] = pho.HoverE;
	dr03EcalIsoPho[nMuMuG] = pho.dr03EcalRecHitSumEtCone;
	dr04HcalIsoPho[nMuMuG] = pho.dr04HcalTowerSumEtCone;
	dr03TrkSumHollowConePho[nMuMuG] = pho.dr03TrkSumPtHollowCone;
	float eT = pho.p4FromVtx(vtx,pho.energy,false).Et();
	isosumoetPho[nMuMuG] = (pho.dr03EcalRecHitSumEtCone + pho.dr04HcalTowerSumEtCone + pho.dr03TrkSumPtHollowCone + isoSumConst - rho*rhoFac)/eT;
	mvaPho[nMuMuG] = getIdMVA(iPho,0,0,false,0);
	genMatchPho[nMuMuG] = (genIndexPho>=0);
	if(genMatchPho[nMuMuG]){
	  genEnergyPho[nMuMuG] = energyGenPho[genIndexPho];
	  genIdPho[nMuMuG] = 22;
	  genStatusPho[nMuMuG] = 1;
	  genIdMomPho[nMuMuG] = pidMomGenPho[genIndexPho];
	  genStatusMomPho[nMuMuG] = 3;
	}else{
	  genEnergyPho[nMuMuG] = 0;
	  genIdPho[nMuMuG] = 0;
	  genStatusPho[nMuMuG] = 0;
	  genIdMomPho[nMuMuG] = 0;
	  genStatusMomPho[nMuMuG] = 0;
	}
	
	nMuMuG++;
      }
    }
  }


  if(doMuMuGamma) outTreeMuMuG->Fill();
}

int HggSelector::getGenMatchPho(VecbosPho *pho){
  const float maxDR = 0.2;
  float bestdE = 9999;
  int indexGenPho=-1;
  for(int iP=0;iP<nGenPho;iP++){
    if(energyGenPho[iP] < 1.) continue;
    if(DeltaR(pho->SC.eta,etaGenPho[iP],pho->SC.phi,phiGenPho[iP]) > maxDR) continue;
    if(fabs(pho->finalEnergy-energyGenPho[iP])/energyGenPho[iP] > 1.) continue;
    if( fabs(pho->finalEnergy-energyGenPho[iP])<bestdE ){
      bestdE = fabs(pho->finalEnergy-energyGenPho[iP]);
      indexGenPho = iP;
    }
  }
  return indexGenPho;
}

int HggSelector::getGenMatchMu(VecbosMu *mu){
  const float maxDR = 0.5;
  float bestdE = 9999;
  int indexGenMu=-1;
  for(int iM=0;iM<nGenMu;iM++){
    if(energyGenMu[iM] < 1.) continue;
    if(DeltaR(mu->eta,etaGenMu[iM],mu->phi,phiGenMu[iM]) > maxDR) continue;
    if(fabs(mu->energy-energyGenMu[iM])/energyGenMu[iM] > 0.5) continue;
    if( fabs(mu->energy-energyGenMu[iM])<bestdE ){
      bestdE = fabs(mu->energy-energyGenMu[iM]);
      indexGenMu = iM;
    }
  }
  return indexGenMu;
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

float HggSelector::getIdMVA(int indexPho1,int indexPho2, bool selPho1, bool usePFSC,int forceVertex=-1){
  
  if(debugSelector) cout << "getIdMVA" <<endl;
  if(debugSelector) cout << "get vertex index ... " << flush;
 
  int selectedVertex = forceVertex;
  if(forceVertex == -1) selectedVertex = getVertexIndex(indexPho1,indexPho2);
  if(debugSelector) cout << "Done" <<endl;
  if(selectedVertex<0 || selectedVertex >= nVtx){
    cout << "WARNING: Photons " << indexPho1 << " and " << indexPho2 
	 << " have no selected vertex!" << endl
	 << "Skipping this pair" << endl;
    return -9999.;
  }
  
  int index = (selPho1==true ? indexPho1: indexPho2);
  if(debugSelector) cout << "index: " << index << endl;
  if(debugSelector) cout << "selectedVertex: " << selectedVertex << " / " << nVtx << endl;
  VecbosPho pho = Photons_->at(index);
  //get the vtx w.r.t. which this has the worst isolation
  float worstVtxTrkIso    = pho.photonWorstIsoDR04.first;
  int worstVtxIndex       = pho.photonWorstIsoDR04.second;
  if(debugSelector) cout << "savedVertices: " << pho.photonTrkIsoFromVtx.size() << endl;
  float selectedVtxTrkIso = pho.photonTrkIsoFromVtx.at(selectedVertex);

  if(debugSelector) cout << worstVtxTrkIso << "  " << worstVtxIndex << "  " 
			 << selectedVtxTrkIso << endl;

  //vetos:
  if(photonMatchedElectron[index] && doElectronVeto) return -8999.;
  
  //set variables
  VecbosSC *SC;
  if(usePFSC) SC = &(pho.PFSC);
  else SC = &(pho.SC);
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);

  float eT = pho.p4FromVtx(vtxPos,pho.finalEnergy,false).Et();
  hoe = pho.HoverE;
  sigietaieta = SC->sigmaIEtaIEta;        
  isosumoet   = (selectedVtxTrkIso + pho.dr03EcalRecHitSumEtCone + pho.dr04HcalTowerSumEtCone + isoSumConst - rho*rhoFac)/eT;
  trkisooet = selectedVtxTrkIso;
  isosumoetbad = (worstVtxTrkIso + pho.dr04EcalRecHitSumEtCone + pho.dr04HcalTowerSumEtCone + isoSumConstBad - rho*rhoFacBad)/eT;
  r9 = SC->r9();          
  ecalisodr03 = pho.dr03EcalRecHitSumEtCone - rho*rhoFac;         
  hcalisodr04 = pho.dr04HcalTowerSumEtCone  - rho*rhoFac;         
  nVertexf  = nVtx;             
  etasc = SC->eta;    
  scetawidth = SC->etaWidth;            
  scphiwidth = SC->phiWidth;            

  if(debugSelector) cout << "et: " << eT << "  " << SC->r9() << "  " << pho.dr03EcalRecHitSumEtCone << "  " << pho.dr04HcalTowerSumEtCone << "  " << selectedVtxTrkIso << "  " 
			 << isosumoet << "  " << pho.dr03TrkSumPtHollowCone << endl;
  if(SC->r9() < 0.9){
    if( (pho.dr03EcalRecHitSumEtCone - 0.012*eT) > 4 
	|| (pho.dr04HcalTowerSumEtCone - 0.005*eT) > 4 
	|| (selectedVtxTrkIso - 0.002*eT) > 4
	|| pho.dr03EcalRecHitSumEtCone > 3 
	|| pho.dr04HcalTowerSumEtCone >  3
	|| isosumoet > 2.8 
	|| pho.dr03TrkSumPtHollowCone > 4) return -9998.;
    if( (pho.isBarrel() && (pho.HoverE > 0.075 || SC->sigmaIEtaIEta > 0.014) )
	|| (!pho.isBarrel() && (pho.HoverE > 0.075 || SC->sigmaIEtaIEta > 0.034) ) ) return -9997.;
  }else{ //(SC->r9() > 0.9
    if( (pho.dr03EcalRecHitSumEtCone - 0.012*eT) > 50 
	|| (pho.dr04HcalTowerSumEtCone - 0.005*eT) > 50 
	|| (selectedVtxTrkIso - 0.002*eT) > 50
	|| pho.dr03EcalRecHitSumEtCone > 3 
	|| pho.dr04HcalTowerSumEtCone >  3
	|| isosumoet > 2.8 
	|| pho.dr03TrkSumPtHollowCone > 4) return -9996.;
    if( (pho.isBarrel() && (pho.HoverE > 0.075 || SC->sigmaIEtaIEta > 0.014) )
	|| (!pho.isBarrel() && (pho.HoverE > 0.075 || SC->sigmaIEtaIEta > 0.034) ) ) return -9995.;
  }

  return (pho.isBarrel() ? photonMVA_EB->EvaluateMVA(methodName_Id) : photonMVA_EE->EvaluateMVA(methodName_Id) );	      
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
