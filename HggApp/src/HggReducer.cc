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
#include "HggReducer.hh"
#include "../src/HggPhysUtils.cc"

#define debugReducer 0

HggReducer::HggReducer(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
  preSelSet = 1; // default
  vertexCFG = "";
  energyScaleCFG = "";
  energySmearCFG = "";
  minPhoSel = 2;
}

HggReducer::HggReducer(TTree *tree, string json, bool goodRunLS, bool isData,int mod) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;
  vertexCFG = "";
  energyScaleCFG = "";
  energySmearCFG = "";
  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
  preSelSet = 1;
  minPhoSel = 2;
}

void HggReducer::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void HggReducer::SetWeight(double weight) {
  _weight = weight;
}

void HggReducer::Loop(string outFileName, int start, int stop) {
  if(fChain == 0){
    cout << "fChain Not defined! QUITTING" << endl;
    return;
  }
  this->init(outFileName);    


  //do initializations (setup output tree, etc.)

  if(debugReducer) cout << "Doing Initialization ... " << flush;
  if(debugReducer) cout << "Done" <<endl;
  
  TRandom3 rng(0);

  if(preSelSet < 0 || preSelSet >= preselections.size()){
    cout << "ERROR: invalid preselection set " << preSelSet << endl;
    return;
  }
 
  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // setup for triggers:
  if(debugReducer) cout << "setting Triggers: " << endl;
  const int nTrigs = triggerNames.size();
  vector<vector<string> > masks;
  for(int iTrig=0;iTrig<nTrigs;iTrig++){
    if(debugReducer) cout << ">> " << triggerNames.at(iTrig) << endl; 
    std::vector<string> tmp;
    tmp.push_back(triggerNames.at(iTrig));// setup the trigger masks for Vecbos
    masks.push_back(tmp);
  }

  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Starting with Entry: " << start << endl;
  cout << "Number of entries = " << stop << endl;
  
  Long64_t jentry = start-1;
  while(fChain->GetEntry(++jentry)){
    if(jentry == stop) break;
    if (jentry%500 == 0) cout << ">>> Processing event # " << jentry << endl;

    runNumberO=runNumber;
    evtNumberO=eventNumber;
    lumiBlockO=lumiBlock;

    //filters:
    if(nTrack > 600 || nPV > 40) {
      if(debugReducer) cout << "dropping event: too many tracks/PVs: " << nTrack << "  " << nPV << endl; 
      continue;
    }

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      // hadronic reload trigger masks and get HLT decisions
      for(int iTrig=0;iTrig<nTrigs;iTrig++){
	setRequiredTriggers(masks.at(iTrig));
	reloadTriggerMask(true);
	triggerBits[iTrig] = hasPassedHLT();
      }
    }

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      if(debugReducer) cout << "Dropping: " << runNumber << ":" << lumiBlock << endl;
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }    

    this->clearAll();
    this->fillVertexInfo();

    this->fillMuons();
    this->fillJets();
    // setup the photons
    std::vector<VecbosPho> tmpPhotons; //temporarily hold a collection of VecbosPhos
    int nSelectedPho=0; // number of selected photons
    if(debugReducer) cout << ">> " << nPho << " Photons" << endl;

    for(int iPho=0;iPho<nPho;iPho++){
      VecbosPho pho(this,iPho);
      if(debugReducer) cout << iPho << ": energy=" <<  pho.energy << "  eta=" << pho.SC.eta << endl;
      if(!this->passPreselection(&pho)) continue;
      if(debugReducer) cout << "pass" << endl;

      //DO ENERGY CORRECTION
      std::pair<float,float> cor = corrector->getPhotonEnergyCorrection(iPho);
      pho.correctedEnergy = cor.first;
      pho.correctedEnergyError = cor.second;
            if(debugReducer) cout << "Corrected Photon: E=" << pho.energy << "  CorE=" << pho.correctedEnergy << "  eta=" 
      	   << pho.SC.eta << "  r9=" << pho.SC.r9 << endl;
      


      if(_isData){ // get the energy scale for data
	if(debugReducer) cout << "Doing Energy Scale ... " << flush;
	std::pair<float,float> dE = energyScale->getDEoE(pho,runNumber);
	if(debugReducer) cout << "Done" << endl;;	
	pho.dEoE    = dE.first;
	pho.dEoEErr = 0; 
	pho.scaledEnergy = pho.correctedEnergy*(1-pho.dEoE);
	pho.scaledEnergyError = pho.correctedEnergyError*(1-(pho.dEoE+pho.dEoEErr));
      }

      if(debugReducer) cout << pho.dEoE <<"   " << pho.scaledEnergy << endl;

      if(!_isData){ // monte carlo, get energy smearing and scale error
	//first get the scale error
	pho.scaledEnergy = pho.correctedEnergy;
	pho.scaledEnergyError = pho.scaledEnergy*energySmear->getMCScaleErr(pho,applyScaleSmear);
	
	//now get the energy smearing
	std::pair<float,float> dE = energySmear->getDEoE(pho,applyScaleSmear);	
	pho.dEoE    = dE.first;
	pho.dEoEErr = dE.second;  // from these, we can generate the energy smearing later (its just a gaussian)
      }

      pho.finalEnergy = pho.scaledEnergy;
      pho.finalEnergyError = pho.scaledEnergyError;

      nSelectedPho++;
      
      // ADD things to the collections
      Photons_.push_back(pho);
      nPho_++;
      if(debugReducer) cout << "aa " << nPho_ << endl;
      //Photons_.push_back(pho);
      //nPho_++;
      /*
      VecbosSC sc = pho.SC;
      SuperClusters_.push_back(sc);
      nSC_++;
      VecbosPFSC pfsc = pho.PFSC;
      PFSuperClusters_.push_back(pfsc);      
      nPFSC++;
      //for(int i=0;i<4;i++){ // push back up to 4 BCs for each SC
      //	if(scSt.BCs[i] > -1) BasicClusters_.push_back(pho.SC.basicClusters[i].getStruct());
      //}
      if(pho.hasConv()){ //there is a valid conversion
	Conversions_.push_back(pho.conversions.at(0));
	nConv_++;
      }
      */
      //end of photon loop
    }//for(int iPho=0; ...

    if(nSelectedPho<minPhoSel){
      if(debugReducer) cout << "Fewer than " << minPhoSel << " selected photons" << endl;
      continue; // skip the event if there are fewer than 2 photons
    }
    if(debugReducer) cout << "More than 1 selected photons" << endl;
    //do the vertexing TMVA
    std::vector<VecbosPho>::iterator iPho1;
    std::vector<VecbosPho>::iterator iPho2;
    if(debugReducer) cout << "aa" << endl;
    nPair_=0;
    for(iPho1 = Photons_.begin(); iPho1 != Photons_.end(); iPho1++){
      if(debugReducer) cout << nPair_ << endl;
      //VecbosPho pho1 = *iPho1;
      for(iPho2 = iPho1+1; iPho2 != Photons_.end(); iPho2++){
	if(debugReducer) cout << ">> " << nPair_ << endl;
	if(iPho2==iPho1) continue;
	if(debugReducer) cout << ">>>>" << endl;
	//VecbosPho pho2 = *iPho2;
	if(debugReducer) cout << "Doing Vertexing ... " << endl;
	pair<int,float> vtxPair = vertexer->vertex_tmva(&*iPho2,&*iPho2);
	int vtx = vtxPair.first;
	if(debugReducer) cout << "Done!  Vtx Index: " << vtx << endl;
	std::pair<unsigned int,int> tmp;
	ggVerticesPhotonIndices.push_back(pair<int,int>(iPho1->index,iPho2->index) );
	ggVerticesVertexIndex.push_back(vtxPair);
	nPair_++;
      }
    }
    
    // fill the isolation variables for the photons
    if(debugReducer) cout << "Filling IsoVars" << endl;
    for(int iPho = 0; iPho < nPho_; iPho++){
      std::vector<float> thisPhoIso;
      float worstIso03 = 0; int worstIsoIndex03 = -1;
      float worstIso04 = 0; int worstIsoIndex04 = -1;
      for(int i=0;i<nVtx;i++){
	float thisIsoDR03 = computeTrackIso(iPho,i,0.,0.3,0.2,0.0,1.0,0.1);
	float thisIsoDR04 = computeTrackIso(iPho,i,0.,0.4,0.2,0.0,1.0,0.1);
	Photons_.at(iPho).photonTrkIsoFromVtx.push_back(thisIsoDR03);
	if(thisIsoDR03 > worstIso03){
	  worstIso03 = thisIsoDR03; worstIsoIndex03 = i;
	}
	if(thisIsoDR04 > worstIso04){
	  worstIso04 = thisIsoDR04; worstIsoIndex04 = i;
	}
      }
      Photons_.at(iPho).photonWorstIsoDR03 = std::pair<float,int>(worstIso03,worstIsoIndex03);
      Photons_.at(iPho).photonWorstIsoDR04 = std::pair<float,int>(worstIso04,worstIsoIndex04);
    } // end of photon isolation loop

    //fill remaining collections
    if(!_isData) this->fillGeneratorInfo();
    this->matchPhotonsElectrons();
   //cout << "Filling Collection" << endl;
    outTree->Fill();
  } // end of main loop

  cout << "Writing Tree:" << endl;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();

  //write output
 file->Close();
}


float HggReducer::computeTrackIso(int iPho, int iVtx,
				  float ptMin,
				  float outerCone,
				  float innerCone,
				  float etaStripHalfWidth,
				  float dzMax,
				  float dxyMax){

  float vX=0,vY=0,vZ=0;
  if(iVtx >= -0 && iVtx < nPV){
    vX = vtxX[iVtx];
    vY = vtxY[iVtx];
    vZ = vtxZ[iVtx];
  }

  double SumTrackPt=0;
  for(int iTrack=0;iTrack<nTrack;iTrack++){
    double ptTrack = sqrt( pxTrack[iTrack]*pxTrack[iTrack] + pyTrack[iTrack]*pyTrack[iTrack] );
    if(ptTrack < ptMin) continue;
    double dZ = fabs( (trackVzTrack[iTrack] - vZ) - 
		      ( (trackVxTrack[iTrack] - vX)*pxTrack[iTrack] 
			+ (trackVyTrack[iTrack] - vY)*pyTrack[iTrack])/ptTrack*pzTrack[iTrack]/ptTrack);
    if(dZ > dzMax) continue;
    double dXY = ( (trackVyTrack[iTrack] - vY)*pyTrack[iTrack] - (trackVxTrack[iTrack] - vX)*pxTrack[iTrack] )/ptTrack;
    if( fabs(dXY) > dxyMax ) continue;
    TVector3 trackP(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
    double dEta = fabs(etaPho[iPho] - trackP.Eta());
    double dR   = DeltaR<float>(etaPho[iPho],trackP.Eta(),phiPho[iPho],trackP.Phi());
    if( dR < outerCone && dR >= innerCone && dEta >= etaStripHalfWidth) SumTrackPt+=ptTrack;
  }
  return SumTrackPt;
}

void HggReducer::matchPhotonsElectrons(){
  vector<VecbosPho>::const_iterator phoIt;
  for(phoIt = Photons_.begin(); phoIt != Photons_.end();phoIt++){
    bool match = false;
    for(int i=0;i<nEle;i++){ 
      //copied from  ConversionTools::hasMatchedPromptElectron
      if(superClusterIndexEle[i] != phoIt->SC.index) continue;
      if(hasMatchedConversionEle[i]) continue;
      if(gsfTrackIndexEle[i]<0 || gsfTrackIndexEle[i] >=nGsfTrack) continue;
      if(expInnerLayersGsfTrack[gsfTrackIndexEle[i]]>0) continue;
      match = true;
      break;	
    }
    photonMatchedElectron[phoIt-Photons_.begin()] = match;
  }
}

void HggReducer::clearAll(){
//clear collections

  nPho_=0;
  Photons_.clear();
  
  nMu_=0;
  Muons_.clear();
  pileupBunchX->clear();
  pileupNInteraction->clear();
  
  ggVerticesPhotonIndices.clear();
  ggVerticesVertexIndex.clear();

  nGenPho=0;
  nGenMu =0;
  nGenEle=0;
  nGenHiggs=0;

  GenPhotons.clear();
  GenMuons.clear();
  GenElectrons.clear();
  GenHiggs.clear();

}

void HggReducer::init(string outputFileName){
  //define the tree
  outTree = new TTree("HggReduce","HggReduce");

 //setup the trigger objects
 triggerBits = new int[triggerNames.size()];

 pileupBunchX = new std::vector<short>; pileupBunchX->clear();
 pileupNInteraction = new std::vector<short>; pileupNInteraction->clear();

 this->setOutputBranches(); // initialize the branches of the output tree
 this->setupPreSelection(); // define the photon preselection cuts

 ReadConfig cfg;
 int errorCode = cfg.read(config);
 if(errorCode){
   cout << "ERROR READING CONFIG FILE!!!! " << endl
	<< "ABORTING" << endl;
   throw 1;
   return;
 }

 string EnergyScaleCFG = cfg.getParameter("EnergyScaleCFG");
 string EnergySmearCFG = cfg.getParameter("EnergySmearCFG");
 string sPreselection   = cfg.getParameter("Preselection");
 string sScaleSmear     = cfg.getParameter("ScaleSmear");
 string sMinPhoSel      = cfg.getParameter("MinPreselPhotons");
 triggerNames  = cfg.getTokens("Triggers",",");

 preSelSet = atoi(sPreselection.c_str());
 applyScaleSmear = atoi(sScaleSmear.c_str());
 if(sMinPhoSel.compare("")!=0) minPhoSel = atoi(sMinPhoSel.c_str());

 cout << "Config Parameters:" << endl
      << "EnergyScale: " << EnergyScaleCFG << endl
      << "EnergySmear: " << EnergySmearCFG << endl
      << "Preselection Set: " << sPreselection << endl
      << "ApplyScaleSmear: "  << sScaleSmear << endl;

 vertexer  = new HggVertexing(this);
 vertexer->setConfigFile(config);
 vertexer->useConversions();
 vertexer->init();
 corrector = new HggEGEnergyCorrector(this,config,_isData);
 energyScale = new HggEnergyScale(EnergyScaleCFG);
 energySmear = new HggEnergyScale(EnergySmearCFG);
 
}

bool HggReducer::passPreselection(VecbosPho *pho){
  if(pho->index == -1 || pho->SC.index == -1) return false;
  PreSelCuts cuts = preselections.at(preSelSet-1);
  if(debugReducer) cout << pho->SC.energy << "  " << pho->SC.eta << "  " << pho->SC.energy/cosh(pho->SC.eta) << endl;
  if( pho->SC.energy/cosh(pho->SC.eta) < cuts.scet2 || 
      (cuts.maxeta>0 && pho->SC.eta > cuts.maxeta) ) return false;  //baseline kinematics
  if(cuts.ecaliso[!pho->isBarrel()]>0 && pho->dr03EcalRecHitSumEtCone >= cuts.ecaliso[!pho->isBarrel()]) return false;
  if(cuts.hcaliso[!pho->isBarrel()]>0 && pho->dr03HcalTowerSumEtCone >= cuts.hcaliso[!pho->isBarrel()]) return false;
  if(cuts.sieie[!pho->isBarrel()]>0 && pho->SC.sigmaIEtaIEta >= cuts.sieie[!pho->isBarrel()]) return false;
  if(cuts.hoe>0 && pho->SC.HoverE >= cuts.hoe) return false;
  return true;

}

void HggReducer::fillMuons(){
  return;
  nMu_=0;
  for(int iMuon = 0; iMuon<nMuon;iMuon++){
    VecbosMu mu(this,iMuon);
    Muons_.push_back(mu);
    nMu_++;
  }
}

void HggReducer::fillJets(){
  const float minPt = 20;
  const int nJetCat=4;
  const float maxJetEta[nJetCat] = {2.5,2.75,3,4.7};
  const float betaStarSlope[nJetCat] = {0.2,0.3,999,999};
  const float rmsCut[nJetCat] = {0.06,0.05,0.05,0.055};

  nJets=0;
  for(int iJ=0;iJ<nAK5PFPUcorrJet; iJ++){
    int jetCat=0;
    for(; jetCat<4; jetCat++){ // get the jet category
      if(fabs(etaAK5PFPUcorrJet[iJ]) <maxJetEta[jetCat]) break;
    }
  if(jetCat>3) continue; //if its >4.7, reject the jet
  float pT = TMath::Sqrt(TMath::Power(pxAK5PFPUcorrJet[iJ],2)+TMath::Power(pyAK5PFPUcorrJet[iJ],2));
  if(pT < minPt) continue;
  if(betastarAK5PFPUcorrJet[iJ] > betaStarSlope[jetCat]*TMath::Log(nPV)-0.64) continue;
  if(rmsCandsHandAK5PFPUcorrJet[iJ] > rmsCut[jetCat]) continue; //jet ID variables

  ptJet[nJets] = pT;
  etaJet[nJets] = etaAK5PFPUcorrJet[iJ];
  phiJet[nJets] = phiAK5PFPUcorrJet[iJ];
  energyJet[nJets] = energyAK5PFPUcorrJet[iJ];
  nJets++;

  }
}
void HggReducer::setupPreSelection(){
  //used to define the sets of prescale cuts we can use
  //maybe this should be totally configurable, but for now hard-code
  PreSelCuts set1;
  set1.scet1 = 35;                                                                                                                              
  set1.scet2 = 25;                                                         
  set1.maxeta = 2.5;                                                       
  set1.ecaliso[0]  = 10;                                                   
  set1.ecaliso[1] = 10;                                                    
  set1.hcaliso[0] = -1;                                                    
  set1.hcaliso[1] = -1;                                                    
  set1.sieie[0] = 0.013;                                                   
  set1.sieie[1] = 0.03;                                                    
  set1.hoe = 0.15;

  PreSelCuts set2;
  set2.scet1 = 20;                                                         
  set2.scet2 = 20;                                                         
  set2.maxeta = 2.5;                                                       
  set2.ecaliso[0]  = -1;                                                   
  set2.ecaliso[1] = -1;                                                    
  set2.hcaliso[0] = -1;                                                    
  set2.hcaliso[1] = -1;                                                    
  set2.sieie[0] = 0.013;                                                   
  set2.sieie[1] = 0.03;                                                    
  set2.hoe = 0.15; 

  PreSelCuts set3;
  set3.scet1 = 20;                                                       
  set3.scet2 = 20;                                                       
  set3.maxeta = -1;                                                      
  set3.ecaliso[0]  = -1;                                                 
  set3.ecaliso[1] = -1;                                                  
  set3.hcaliso[0] = -1;                                                  
  set3.hcaliso[1] = -1;                                                  
  set3.sieie[0] = -1;                                                    
  set3.sieie[1] = -1;                                                    
  set3.hoe = -1;       

  PreSelCuts set4;
  set4.scet1 = 28;                                                       
  set4.scet2 = 28;                                                       
  set4.maxeta = -1;                                                      
  set4.ecaliso[0]  = -1;                                                 
  set4.ecaliso[1] = -1;                                                  
  set4.hcaliso[0] = -1;                                                  
  set4.hcaliso[1] = -1;                                                  
  set4.sieie[0] = -1;                                                    
  set4.sieie[1] = -1;                                                    
  set4.hoe = -1;        

  PreSelCuts set5;
  set5.scet1 = 25;                                                       
  set5.scet2 = 25;                                                       
  set5.maxeta = -1;                                                      
  set5.ecaliso[0]  = -1;                                                 
  set5.ecaliso[1] = -1;                                                  
  set5.hcaliso[0] = -1;                                                  
  set5.hcaliso[1] = -1;                                                  
  set5.sieie[0] = -1;                                                    
  set5.sieie[1] = -1;                                                    
  set5.hoe = -1;    

  PreSelCuts set6;
  set6.scet1 = 0;                                                       
  set6.scet2 = 0;                                                       
  set6.maxeta = -1;                                                      
  set6.ecaliso[0]  = -1;                                                 
  set6.ecaliso[1] = -1;                                                  
  set6.hcaliso[0] = -1;                                                  
  set6.hcaliso[1] = -1;                                                  
  set6.sieie[0] = -1;                                                    
  set6.sieie[1] = -1;                                                    
  set6.hoe = -1;    

  preselections.push_back(set1);
  preselections.push_back(set2);
  preselections.push_back(set3);
  preselections.push_back(set4);
  preselections.push_back(set5);
  preselections.push_back(set6);

}

void HggReducer::fillVertexInfo(){
  nVtx = nPV;
  for(int i=0;i<nPV;i++){
    vtxX[i] = PVxPV[i];             
    vtxY[i] = PVyPV[i];             
    vtxZ[i] = PVzPV[i];             
    vtxChi2[i] = chi2PV[i];          
    vtxNdof[i] = ndofPV[i];          
    vtxNormalizedChi2[i] = normChi2PV[i];
    vtxTrackSize[i] = trackSizePV[i];       
    vtxIsFake[i] = isFakePV[i];          
    vtxIsValid[i] = isValidPV[i];         
  }
  rho=rhoJetsFastJet;  // rho from kt6PFJets
}

void HggReducer::fillGeneratorInfo(){
  //GENERATOR information
  nGenPho=0;
  nGenMu =0;
  nGenEle=0;
  nGenHiggs=0;

  for(int iGen=0;iGen<nMc;iGen++){
    VecbosGen part(this,iGen);

    switch(abs(idMc[iGen])){
    case 25:  //higgs
      GenHiggs.push_back(part); break;
    case 22: //Photons
      GenPhotons.push_back(part); break;
    case 13: //Muons
      GenMuons.push_back(part); break;
    case 11: //Electrons
      GenElectrons.push_back(part); break;
    Default:
      GenOthers.push_back(part); break;
    }
    nGenHiggs = GenHiggs.size();
    nGenPho   = GenPhotons.size();
    nGenEle   = GenElectrons.size();
    nGenMu    = GenMuons.size();
    nGenOthers= GenOthers.size();
  }    
  procID = 0; //genProcessId;
  qScale;
  nPu = nPU[1];
  
  //match the reco objects to gen objects
  
  PhoCollection::iterator phoIt;
  for(phoIt = Photons_.begin(); phoIt !=Photons_.end(); phoIt++)
    phoIt->doGenMatch(this);
  MuCollection::iterator muIt;
  for(muIt = Muons_.begin(); muIt !=Muons_.end(); muIt++)
    muIt->doGenMatch(this);
  
}

void HggReducer::setOutputBranches(){

  //Event info
  outTree->Branch("lumiBlock",&lumiBlockO,"lumiBlock/I");
 outTree->Branch("runNumber",&runNumberO,"runNumber/I");
outTree->Branch("evtNumber",&evtNumberO,"evtNumber/I");
 outTree->Branch("bunchX",&bunchX,"bunchX/I");
 outTree->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
 outTree->Branch("evtTime",&evtTime,"evtTime/I");
 outTree->Branch("isRealData",&_isData,"isRealData/I");
  

 ///information for the vertex
 outTree->Branch("nVtx",&nVtx,"nVtx/I");
 outTree->Branch("vtxX",vtxX,"vtxX[nVtx]/F");
 outTree->Branch("vtxY",vtxY,"vtxY[nVtx]/F");
 outTree->Branch("vtxZ",vtxZ,"vtxZ[nVtx]/F");
 outTree->Branch("vtxChi2",vtxChi2,"vtxChi2[nVtx]/F");
 outTree->Branch("vtxNdof",vtxNdof,"vtxNdof[nVtx]/F");
 outTree->Branch("vtxNormalizedChi2",vtxNormalizedChi2,"vtxNormalizedChi2[nVtx]/F");
 outTree->Branch("vtxTrackSize",vtxTrackSize,"vtxTrackSize[nVtx]/I");
 outTree->Branch("vtxIsFake",vtxIsFake,"vtxIsFake[nVtx]/I");
 outTree->Branch("vtxIsValid",vtxIsValid,"vtxIsValid[nVtx]/I");
 
 //physics declared -- should be set by the JSON, but can't hurt
 outTree->Branch("phyDeclared",&phyDeclared,"phyDeclared/I");
 
 
 outTree->Branch("rho", &rho,"rho/F");
 outTree->Branch("rhoEtaMax44", &rhoEtaMax44,"rhoEtaMax44/F");
 
 
 
 outTree->Branch("pileupBunchX","std::vector<short>", &pileupBunchX);
 outTree->Branch("pileupNInteraction","std::vector<short>", &pileupNInteraction);
 outTree->Branch("pileupTrueNumInterations",&pileupTrueNumInterations);
 

  //trigger -- here we depart a bit from Yong's original code
 for(int i=0;i<triggerNames.size();i++){
   outTree->Branch(triggerNames.at(i).c_str(), &(triggerBits[i]), Form("%s/I",triggerNames.at(i).c_str()) );  // this will produce 1 int per trigger in the output tree
 }
 

 //objects
 /*
 outTree->Branch("nSC",&nSC_);
 outTree->Branch("nPFSC",&nPFSC_);
 outTree->Branch("nConv",&nConv_);
 outTree->Branch("SuperClusters",&SuperClusters_);                          
 outTree->Branch("PFSuperClusters",&PFSuperClusters_);                          
 //outTree->Branch("BasicClusters",&BasicClusters_);                          
 outTree->Branch("Conversions",&Conversions_);
 */
 outTree->Branch("nPho",&nPho_);
 outTree->Branch("nPair",&nPair_);
 outTree->Branch("Photons",&Photons_);
 outTree->Branch("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
 outTree->Branch("ggVerticesVertexIndex",&ggVerticesVertexIndex);
 outTree->Branch("photonMatchedElectron",photonMatchedElectron,"photonMatchedElectron[nPho]/O");

 outTree->Branch("nMu",&nMu_,"nMu/I");
 outTree->Branch("Muons",&Muons_);

 outTree->Branch("nJets",&nJets,"nJets/I");
 outTree->Branch("ptJet",ptJet,"ptJet[nJets]");
 outTree->Branch("etaJet",etaJet,"etaJet[nJets]");
 outTree->Branch("phiJet",phiJet,"phiJet[nJets]");
 outTree->Branch("energyJet",energyJet,"energyJet[nJets]");

 //FOR MONTE CARLO:
  if(!_isData){
    //generator level information
    outTree->Branch("procID",&procID,"procID/I");
    outTree->Branch("qScale",&qScale,"qScale/F");
    outTree->Branch("nPU",&nPu,"nPU/F");

    ///gen electron, muon,photon
    outTree->Branch("nGenPho",&nGenPho,"nGenPho/I");
    outTree->Branch("GenPhotons",&GenPhotons);

    outTree->Branch("nGenEle",&nGenEle,"nGenEle/I");
    outTree->Branch("GenElectrons",&GenElectrons);
    
    outTree->Branch("nGenMu",&nGenMu,"nGenMu/I");
    outTree->Branch("GenMuons",&GenMuons);
    //generator level higgs
    outTree->Branch("nGenHiggs",&nGenHiggs,"nGenHiggs/I");
    outTree->Branch("GenHiggs",&GenHiggs);
    //other particles
    outTree->Branch("nGenOthers",&nGenOthers,"nGenOthers/I");
    outTree->Branch("GenOthers",&GenOthers);
  }

}
