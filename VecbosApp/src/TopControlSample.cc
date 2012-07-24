#include <string>
#include <iostream>

#include <TTree.h>
#include <TObjString.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/BTagAlgoBits.h"
#include "EgammaAnalysisTools/include/ElectronTrackerIsolation.hh"
#include "EgammaAnalysisTools/include/ElectronCaloIsolation.hh"
#include "EgammaAnalysisTools/include/ElectronBestCandidateSelector.hh"
#include "include/CaloTower.hh"
#include "include/TopControlSample.hh"
#include "include/JetCounter.hh"


using namespace bits;

TopControlSample::TopControlSample(TTree *tree) 
  : Vecbos(tree) {
  
  std::string configDir;
  configDir="config/top/";
  
  // kinematic selection and counters
  std::string fileCuts     = configDir + "cuts.txt";
  std::string fileSwitches = configDir + "switches.txt";  
  _selection = new Selection(fileCuts,fileSwitches);  

  _pCounterEE = &_counterEE;
  _pCounterMM = &_counterMM;
  _pCounterEM = &_counterEM;

  ConfigureSelection(_selection, _pCounterEE);
  _pCounterEE->SetTitle("EE counter");
  ConfigureSelection(_selection, _pCounterMM);
  _pCounterMM->SetTitle("MM counter");
  ConfigureSelection(_selection, _pCounterEM);
  _pCounterEM->SetTitle("EM counter");


  // configuring eleID configuration
  TString selectionString(_selection->getStringParameter("electronIDType"));
  EgammaTightID_NC.ConfigureNoClass("config/vecbos/electronId/"+selectionString);  
  EgammaTightID_NC.ConfigureEcalCleaner("config/vecbos/electronId/");              
    
  // To read good run list!
  std::cout << "[GoodRunLS]::goodRunLS is " << _selection->getSwitch("goodRunLS") 
	    << " isData is " <<  _selection->getSwitch("isData") << std::endl;
  if (_selection->getSwitch("goodRunLS") && _selection->getSwitch("isData")) {
    std::string goodRunJsonFile = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }
  
  // MC truth
  if ( _selection->getSwitch("isData") ) mcevent.SetData(true);  
  TopToEEDecay = TopToMMDecay = TopToEMDecay = 0;  
}

TopControlSample::~TopControlSample(){ 

  myOutTree -> save();

  delete _selection;
}

// loop over events - real analysis
void TopControlSample::Loop() {
  
  _verbose=false;
  if(fChain == 0) return;

  // reduced tree
  sprintf(namefile,"%sOutTree.root",m_prefix);
  myOutTree = new RedTopTree(namefile);

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Total number of entries in the chain = " << nentries << std::endl;
  int maxEvents = nentries;

  // lumi
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  for (Long64_t jentry=0; jentry<maxEvents;jentry++) {
    
    bool eventToAnalyze = true;
    
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    if ( eventToAnalyze ) {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE 
      // WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
      bool newTriggerMask = true;
      reloadTriggerMask(newTriggerMask);
      
      // Good Run selection
      if (_selection->getSwitch("isData") && _selection->getSwitch("goodRunLS") 
	  && !isGoodRunLS()) {
	if ( lastRun!= runNumber || lastLumi != lumiBlock) {
	  lastRun  = runNumber;
	  lastLumi = lumiBlock;
	  std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi 
		    << " is rejected" << std::endl;
	}
	continue;
      }
      if (_selection->getSwitch("isData") && _selection->getSwitch("goodRunLS") 
	  && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
      }

      // all good events
      _counterEE.IncrVar("event",1.); 
      _counterMM.IncrVar("event",1.); 
      _counterEM.IncrVar("event",1.); 


      // for MC: to know if this is a tt->ele, mu, emu 
      /*
      if ( !_selection->getSwitch("isData") ) { 
	
	int indexMcTtoEle,     indexMcTtoPos     = -1;
	int indexMcTtoMu,      indexMcTtoAMu     = -1;
	int indexMcTtoEleTau,  indexMcTtoMuTau   = -1;
	int indexMcTtoPosATau, indexMcTtoAMuATau = -1;
	
	mcevent.LoadDecay(nMc,idMc,mothMc);
        mcevent.LoadMomentum(pMc,energyMc,thetaMc,phiMc);
	indexMcTtoEle     = mcevent.indexEleTopPrompt();
        indexMcTtoPos     = mcevent.indexPosTopPrompt();
        indexMcTtoMu      = mcevent.indexMuTopPrompt();
        indexMcTtoAMu     = mcevent.indexAMuTopPrompt();
	indexMcTtoEleTau  = mcevent.indexTauEleTopPrompt();
	indexMcTtoPosATau = mcevent.indexATauPosTopPrompt();
	indexMcTtoMuTau   = mcevent.indexTauMuTopPrompt();
	indexMcTtoAMuATau = mcevent.indexATauAMuTopPrompt();

	TopToEEDecay = TopToMMDecay = TopToEMDecay = 0;  

	if ( indexMcTtoEle>-1   && indexMcTtoPos>-1 )    TopToEEDecay = 1;
	if ( indexMcTtoEle>-1   && indexMcTtoPosATau>-1) TopToEEDecay = 1;
	if ( indexMcTtoPos>-1   && indexMcTtoEleTau>-1)  TopToEEDecay = 1;
	if ( indexMcTtoEle>-1   && indexMcTtoEleTau>-1)  TopToEEDecay = 1;
	if ( indexMcTtoPos>-1   && indexMcTtoPosATau>-1) TopToEEDecay = 1;

	if ( indexMcTtoMu>-1    && indexMcTtoAMu>-1 )    TopToMMDecay = 1;
	if ( indexMcTtoMu>-1    && indexMcTtoAMuATau>-1) TopToMMDecay = 1;
	if ( indexMcTtoAMu>-1   && indexMcTtoMuTau>-1)   TopToMMDecay = 1;
	if ( indexMcTtoMu>-1    && indexMcTtoMuTau>-1)   TopToMMDecay = 1;
	if ( indexMcTtoAMu>-1   && indexMcTtoAMuATau>-1) TopToMMDecay = 1;	

	if ( indexMcTtoEle>-1 && indexMcTtoAMu>-1 )        TopToEMDecay = 1;
	if ( indexMcTtoPos>-1 && indexMcTtoMu>-1 )         TopToEMDecay = 1;
	if ( indexMcTtoEle>-1 && indexMcTtoAMuATau>-1 )    TopToEMDecay = 1;
	if ( indexMcTtoPos>-1 && indexMcTtoMuTau>-1 )      TopToEMDecay = 1;
	if ( indexMcTtoAMu>-1 && indexMcTtoEleTau>-1)      TopToEMDecay = 1;
	if ( indexMcTtoMu>-1 && indexMcTtoPosATau>-1)      TopToEMDecay = 1;
	if ( indexMcTtoEleTau>-1 && indexMcTtoAMuATau>-1 ) TopToEMDecay = 1;
	if ( indexMcTtoPosATau>-1 && indexMcTtoMuTau>-1 )  TopToEMDecay = 1;
      */
	/*
	  if ( !TopToEEDecay && !TopToEMDecay && !TopToMMDecay) {
	  cout << endl;
	  cout << endl;
	  for(int imc=0; imc<nMc; imc++) {
	  int theid      = idMc[imc];
	  int themoth    = mothMc[imc];
	  int theGmoth   = mothMc[themoth];
	  int themothid  = idMc[themoth];
	  int theGmothid = idMc[theGmoth];
	  cout << "# = " << imc << ", id = " << theid << ", moth# = " << themoth << ", gMoth# = " << theGmoth 
	  << ", mothId = " << themothid << ", gMothId = " << theGmothid << endl; 
	  }
	  }*/
	
      /*
	if (_selection->getSwitch("mcTruth")) {
	  if (TopToEEDecay==0 && TopToMMDecay==0 && TopToEMDecay==0) continue;
	}
      }
    */
      _counterEE.IncrVar("mcTruth",1.); 
      _counterMM.IncrVar("mcTruth",1.); 
      _counterEM.IncrVar("mcTruth",1.); 
      // if( TopToEEDecay ) _counterEE.IncrVar("mcTruth",1.); 
      // if( TopToMMDecay ) _counterMM.IncrVar("mcTruth",1.); 
      // if( TopToEMDecay ) _counterEM.IncrVar("mcTruth",1.); 
      

      // trigger selection
      Utils anaUtils;
      bool passedHLT = hasPassedHLT();
      if (_selection->getSwitch("trigger") && !passedHLT) continue; 
      _counterEE.IncrVar("trigger",1.);
      _counterMM.IncrVar("trigger",1.);
      _counterEM.IncrVar("trigger",1.);

      // good events cleaning cuts: tracks quality
      bool okTracker = true;
      if (_selection->getSwitch("goodTracks")) {
	if (nTrack>10) {
	  int goodTracks=0;
	  for (int iTrack=0; iTrack<nTrack; iTrack++) {
	    bool goodTrack = (qualityMaskTrack[iTrack]>>2)%2;
	    if (goodTrack) goodTracks++;
	  }
	  float goodTracksF = (float)(goodTracks);
	  float totTracksF  = (float)(nTrack);
	  float tracksRatio = goodTracksF/totTracksF;
	  if(!_selection->passCut("goodTracksRatio",tracksRatio)) okTracker = false;
	}
      }
      if (!okTracker) continue;
      _counterEE.IncrVar("goodTracks",1.);
      _counterMM.IncrVar("goodTracks",1.);
      _counterEM.IncrVar("goodTracks",1.);

      // at least one good primary vertex required
      if (nPV<1) continue;
      bool isGoodPV = false;
      for (int iPV=0; iPV<nPV; iPV++) {
	float rhoVtx = sqrt(PVxPV[iPV]*PVxPV[iPV] + PVyPV[iPV]*PVyPV[iPV]);
	if(!_selection->passCut("pvNdf",ndofPV[iPV])) continue; 
	if(!_selection->passCut("pvDz", PVzPV[iPV]))  continue; 
	if(!_selection->passCut("pvRho", rhoVtx))     continue; 
	isGoodPV = true;  
      }
      if (!isGoodPV) continue;            
      _counterEE.IncrVar("primaryVertex",1.);
      _counterMM.IncrVar("primaryVertex",1.);
      _counterEM.IncrVar("primaryVertex",1.);

      // use the highest pT good PV in the event (here to check muonID)
      hardestPV = -1;
      float sumPtMax = 0.0;
      for (int iPV=0; iPV<nPV; iPV++) {
	float rhoVtx = sqrt(PVxPV[iPV]*PVxPV[iPV] + PVyPV[iPV]*PVyPV[iPV]);
	if(!_selection->passCut("pvNdf",ndofPV[iPV])) continue; 
	if(!_selection->passCut("pvDz", PVzPV[iPV]))  continue; 
	if(!_selection->passCut("pvRho", rhoVtx))     continue; 
	if(SumPtPV[iPV] > sumPtMax) {
	  sumPtMax = SumPtPV[iPV];
	  hardestPV = iPV;
	}
      }

      // electrons
      std::vector<int> acceptElectrons,  idElectrons; 
      std::vector<int> isolElectrons,    convRejElectrons;
      std::vector<int> furtherElectrons, vertexElectrons;
      trackerIsolation.clear();
      ecalIsolation.clear();
      hcalIsolation.clear();
      combinedIsolation.clear();

      // electrons in the acceptance
      for(int ii=0;ii<nEle;ii++) {             
	bool isGoodEle = true;
	TVector3 pLepton(pxEle[ii],pyEle[ii],pzEle[ii]);
	float thisPt=pLepton.Pt();
	if ( !_selection->passCut("etaElectronAcc",etaEle[ii])) isGoodEle = false;
	if ( !_selection->passCut("ptElectronAcc",thisPt) )     isGoodEle = false;
	if ( isGoodEle ) acceptElectrons.push_back(ii);
	
	// assign to each electron the three isolations
	trackerIsolation.push_back( dr03TkSumPtEle[ii] / pLepton.Pt() );
	ecalIsolation.push_back( dr03EcalRecHitSumEtEle[ii] / pLepton.Pt() );
	hcalIsolation.push_back(  dr03HcalTowerSumEtEle[ii] / pLepton.Pt() );

	// and the combined isolation (of relative isolations)
	trackerIsolRel = dr03TkSumPtEle[ii] / pLepton.Pt();
	ecalIsolRel    = dr03EcalRecHitSumEtEle[ii] / pLepton.Pt();
	hcalIsolRel    = dr03HcalTowerSumEtEle[ii] / pLepton.Pt();

	if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[ii],isEB)) {
	  combinedIsolation.push_back(  (dr03TkSumPtEle[ii] + (float) TMath::Max(float(dr03EcalRecHitSumEtEle[ii]-1.),float(0.)) + dr03HcalTowerSumEtEle[ii]) / (TMath::Max(pLepton.Pt(),20.)) );
	}
	else {
	  combinedIsolation.push_back(  (dr03TkSumPtEle[ii] + dr03EcalRecHitSumEtEle[ii] + dr03HcalTowerSumEtEle[ii]) / (TMath::Max(pLepton.Pt(),20.)) );
	}
      }
      int howManyAccEle = acceptElectrons.size();

      // 'tight' (optimized) identified electrons 
      for(int ii=0;ii<howManyAccEle;ii++) {      
        bool isGoodEle  = true;
	// id
        int accEleIndex = acceptElectrons[ii];
        bool eleId, isol, conv;
        eleId = isol = conv = false;
        isEleID(&EgammaTightID_NC,accEleIndex,&eleId,&isol,&conv);
        if (!eleId) isGoodEle = false;

	// not tracker driven only
	bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[accEleIndex], bits::isEcalDriven);
	if(!ecalDriven) isGoodEle = false;

	if (isGoodEle) idElectrons.push_back(acceptElectrons[ii]);            
      }
      int howManyIdEle = idElectrons.size();

      // tight optimized, tight isolated electrons
      for(int ii=0;ii<howManyIdEle;ii++) {      
        bool isGoodEle = true;
        int idEleIndex=idElectrons[ii];
	bool eleId, isol, conv;
        eleId = isol = conv = false;
        isEleID(&EgammaTightID_NC,idEleIndex,&eleId,&isol,&conv);
	if (!isol) isGoodEle = false;
        if (isGoodEle) isolElectrons.push_back(idEleIndex); 
      }
      int howManyIsolEle = isolElectrons.size();

      // tight conversion rejection electrons
      for(int ii=0;ii<howManyIsolEle;ii++) {      
        bool isGoodEle = true;
        int isolEleIndex=isolElectrons[ii];
	bool eleId, isol, conv;
        eleId = isol = conv = false;
        isEleID(&EgammaTightID_NC,isolEleIndex,&eleId,&isol,&conv);
	if (!conv) isGoodEle = false;
        if (isGoodEle) convRejElectrons.push_back(isolEleIndex); 
      }
      int howManyConvRejEle = convRejElectrons.size();
      
      // further requirements to match top analysis: 
      // ele far from global or tracker muons
      // SC eT > 10 GeV
      for(int iele=0; iele<howManyConvRejEle; iele++) {	
	bool isGoodEle = true;

	int thisEleIndex = convRejElectrons[iele];
	TVector3 pThisEle(pxEle[thisEleIndex],pyEle[thisEleIndex],pzEle[thisEleIndex]);
	float thisPt = pThisEle.Pt();

	int thisEleScIndex = superClusterIndexEle[thisEleIndex];
	float thisScEt = energySC[thisEleScIndex]*sin(thetaSC[thisEleScIndex]);
	if (thisScEt<15) isGoodEle = false;

	for(int iMu=0;iMu<nMuon;iMu++) {             
	  Utils anaUtilsMuon; 
	  bool muonOrTracker = anaUtilsMuon.muonIdVal(muonIdMuon[iMu],AllGlobalMuons) || anaUtilsMuon.muonIdVal(muonIdMuon[iMu],AllTrackerMuons);
	  if(!muonOrTracker) continue;

	  TVector3 pThisMu(pxMuon[iMu],pyMuon[iMu],pzMuon[iMu]);
	  float deltaREleMu = pThisMu.DeltaR(pThisEle);
	  if (deltaREleMu<0.1) isGoodEle = false;
	}
	
	if (isGoodEle) furtherElectrons.push_back(thisEleIndex); 
      }
      int howManyFurtherEle = furtherElectrons.size();

      // muons
      std::vector<int> acceptMuons, idMuons, isolMuons, vertexMuons;

      // muons in the acceptance 
      for(int ii=0;ii<nMuon;ii++) {             
	bool isGoodMu = true;
	TVector3 pLepton(pxMuon[ii],pyMuon[ii],pzMuon[ii]);
	float thisPt=pLepton.Pt();
	if ( !_selection->passCut("etaMuonAcc",etaMuon[ii])) isGoodMu = false;
	if ( !_selection->passCut("ptMuonAcc",thisPt) )      isGoodMu = false;
	if ( isGoodMu ) acceptMuons.push_back(ii);
      }
      int howManyAccMu = acceptMuons.size();
      
      // identified muons
      for(int ii=0;ii<howManyAccMu;ii++) {      
        bool isGoodMu = true;
        isMuonID(acceptMuons[ii],&isGoodMu);            
        if (isGoodMu) idMuons.push_back(acceptMuons[ii]);            
      }
      int howManyIdMu = idMuons.size();

      // isolated muons
      for(int ii=0;ii<howManyIdMu;ii++) {      
        bool isGoodMu = true;
        int idMuIndex = idMuons[ii];
	TVector3 pLepton(pxMuon[idMuIndex],pyMuon[idMuIndex],pzMuon[idMuIndex]);
	float thisPt=pLepton.Pt();

	float trackerIsolMuon = sumPt03Muon[idMuIndex];   
	float ecalIsolMuon    = emEt03Muon[idMuIndex];  
	float hcalIsolMuon    = hadEt03Muon[idMuIndex];
	float maxDen = thisPt;
	if (thisPt<20) maxDen = 20.;
	float combinedIsolMuonRel = (trackerIsolMuon+ecalIsolMuon+hcalIsolMuon)/maxDen;
	if ( !_selection->passCut("isolMuon",combinedIsolMuonRel) ) isGoodMu = false;
	if (isGoodMu) isolMuons.push_back(idMuIndex);  
      }
      int howManyIsolMu = isolMuons.size();

      // dummy
      int howManyConvRejMu = howManyIsolMu;
      int howManyFurtherMu = howManyIsolMu;


      // leptons consistency with primary vertex:
      // match within XX cm along the beamline 
      float z0 = PVzPV[hardestPV];       
      float dzEle = -100.;
      int nVtxInputEle = howManyFurtherEle;
      for(int iele=0; iele<nVtxInputEle; iele++) {	
	bool isGoodEle = true;
	int thisEle  = furtherElectrons[iele];
	int gsfTrack = gsfTrackIndexEle[thisEle];
	float dz = trackVzGsfTrack[gsfTrack]-z0;
	dzEle = dz;
	if(_selection->getSwitch("vertexCompatibility") && !_selection->passCut("dzVertexEle",dz)) isGoodEle = false;
	if(isGoodEle) vertexElectrons.push_back(thisEle); 
      }
      int howManyVertexEle = vertexElectrons.size();

      float dzMu = -100.;
      int nVtxInputMu = howManyIsolMu;
      for(int imu=0; imu<nVtxInputMu; imu++) {	
	bool isGoodMu = true;
	int thisMu  = isolMuons[imu];
	int kfTrack = trackIndexMuon[thisMu];
	float dz = trackVzTrack[kfTrack]-z0;
	dzMu = dz;
	if(_selection->getSwitch("vertexCompatibility") && !_selection->passCut("dzVertexMu",dz)) isGoodMu = false;
	if(isGoodMu) vertexMuons.push_back(thisMu); 
      }
      int howManyVertexMu = vertexMuons.size();

      
      // choose the reconstruction channel using only 
      // leptons in the acceptance (no isol/id/conv)
      // no choice applied here: one event can enter several counters
      // this is only to count for this part of selection!
      std::pair<int,int> bestElectronPair_soft = getBestGoodElePair(acceptElectrons);
      int theEle1_soft(bestElectronPair_soft.first);
      int theEle2_soft(bestElectronPair_soft.second);
      std::pair<int,int> bestMuonPair_soft = getBestGoodMuonPair(acceptMuons);
      int theMu1_soft(bestMuonPair_soft.first);
      int theMu2_soft(bestMuonPair_soft.second);

      // reconstructed channels
      m_channelEE_soft = false;     
      m_channelMM_soft = false;
      m_channelEM_soft = false;
      if (theEle1_soft > -1 && theEle2_soft > -1) m_channelEE_soft = true;
      if (theMu1_soft  > -1 && theMu2_soft  > -1) m_channelMM_soft = true;
      if (theEle1_soft > -1 && theMu1_soft  > -1) m_channelEM_soft = true;

      if(m_channelEE_soft) {

	bool isGood = true;
	
	if (_selection->getSwitch("leptonAcceptance") && howManyAccEle<2) isGood = false;
	if (isGood) _counterEE.IncrVar("leptonAcceptance",1.);  
	
	if (_selection->getSwitch("leptonId") && howManyIdEle<2) isGood = false;
	if (isGood) _counterEE.IncrVar("leptonId",1.);

	if (_selection->getSwitch("leptonIsol") && howManyIsolEle<2) isGood = false;
	if (isGood) _counterEE.IncrVar("leptonIsol",1.);

	if (_selection->getSwitch("convRej") && howManyConvRejEle<2) isGood = false;
	if (isGood) _counterEE.IncrVar("convRej",1.);

	if (_selection->getSwitch("further") && howManyFurtherEle<2) isGood = false;
	if (isGood) _counterEE.IncrVar("further",1.);

	if (_selection->getSwitch("vertexCompatibility") && howManyVertexEle<2) isGood = false;
	if (isGood) _counterEE.IncrVar("vertexCompatibility",1.);
      }


      if(m_channelMM_soft) {

	bool isGood = true;

	if (_selection->getSwitch("leptonAcceptance") && howManyAccMu<2) isGood = false;
	if (isGood) _counterMM.IncrVar("leptonAcceptance",1.);  
	      
	if (_selection->getSwitch("leptonId") && howManyIdMu<2) isGood = false;
	if (isGood) _counterMM.IncrVar("leptonId",1.);

	if (_selection->getSwitch("leptonIsol") && howManyIsolMu<2) isGood = false;
	if (isGood) _counterMM.IncrVar("leptonIsol",1.);

	if (_selection->getSwitch("convRej") && howManyConvRejMu<2) isGood = false;
	if (isGood) _counterMM.IncrVar("convRej",1.);

	if (_selection->getSwitch("further") && howManyFurtherMu<2) isGood = false;
	if (isGood) _counterMM.IncrVar("further",1.);

	if (_selection->getSwitch("vertexCompatibility") && howManyVertexMu<2) isGood = false;
	if (isGood) _counterMM.IncrVar("vertexCompatibility",1.);
      }

      if(m_channelEM_soft) {

	bool isGood = true;
	
	if (_selection->getSwitch("leptonAcceptance") && (howManyAccEle<1 || howManyAccMu<1)) isGood = false;
	if (isGood) _counterEM.IncrVar("leptonAcceptance",1.);  
	      
	if (_selection->getSwitch("leptonId") && (howManyIdEle<1 || howManyIdMu<1)) isGood = false;
	if (isGood) _counterEM.IncrVar("leptonId",1.);

	if (_selection->getSwitch("leptonIsol") && (howManyIsolEle<1 || howManyIsolMu<1)) isGood = false;
	if (isGood) _counterEM.IncrVar("leptonIsol",1.);

	if (_selection->getSwitch("convRej") && (howManyConvRejEle<1 || howManyConvRejMu<1)) isGood = false;
	if (isGood) _counterEM.IncrVar("convRej",1.);

	if (_selection->getSwitch("further") && (howManyFurtherEle<1 || howManyFurtherMu<1)) isGood = false;
	if (isGood) _counterEM.IncrVar("further",1.);

	if (_selection->getSwitch("vertexCompatibility") && (howManyVertexEle<1 || howManyVertexMu<1)) isGood = false;
	if (isGood) _counterEM.IncrVar("vertexCompatibility",1.);
      }
      
      // get the best electrons and muons - without charge requirement
      std::pair<int,int> bestElectronPair = getBestGoodElePair(vertexElectrons);
      int theEle1(bestElectronPair.first);
      int theEle2(bestElectronPair.second);
      std::pair<int,int> bestMuonPair = getBestGoodMuonPair(vertexMuons);
      int theMu1(bestMuonPair.first);
      int theMu2(bestMuonPair.second);

      // reconstructed channels
      m_channelEE = false;     
      m_channelMM = false;
      m_channelEM = false;
      if (theEle1 > -1 && theEle2 > -1) m_channelEE = true;
      if (theMu1  > -1 && theMu2  > -1) m_channelMM = true;
      if (theEle1 > -1 && theMu1  > -1) m_channelEM = true;

      if (_verbose) {
	std::cout << "nEle = " << nEle << "\tnMuon = " << nMuon << std::endl;
	std::cout << "indices: " << theEle1 << " " << theEle2 << " " << theMu1 << " " << theMu2 << std::endl;
	std::cout << "chargeEle1 = "    << chargeEle[theEle1] << "\tchargeEle2 = "     << chargeEle[theEle2] 
	  	  << "\tchargeMuon1 = " << chargeMuon[theMu1] << "\tchargeMuonPlus = " << chargeMuon[theMu2] 
		  << std::endl;
	std::cout << "ee = "    << m_channelEE 
		  << "\tmm = "  << m_channelMM 
		  << "\temu = " << m_channelEM << std::endl; 
      }

      // QUI DEVO SCEGLIERE: scelta .... 1=EM, 2=MM, 3=EE
      if( m_channelEM ) { 
	m_channelMM = false;
	m_channelEE = false;
      }	
      if( m_channelMM && !m_channelEM ) { 
	m_channelEE = false;
      }	

      if (!m_channelEE && !m_channelMM && !m_channelEM) continue;

      // starting here we split for different channels and one event enters only one class
      if (m_channelEE) _counterEE.IncrVar("recoChannel",1);      
      if (m_channelMM) _counterMM.IncrVar("recoChannel",1);      
      if (m_channelEM) _counterEM.IncrVar("recoChannel",1);      

      // set the event kinematics
      if(m_channelEE) setKinematics2Ele(theEle1, theEle2);
      if(m_channelMM) setKinematics2Mu (theMu1,  theMu2);
      if(m_channelEM) setKinematics1Ele1Mu(theEle1, theMu1);

      // to check if the selected reco lepton is the one matching MC truth (lepton from W)
      TVector3 lept1_3V, lept2_3V;
      if(m_channelEE) { 
	lept1_3V = p4Ele1_.Vect();
	lept2_3V = p4Ele2_.Vect();
      }
      if(m_channelMM) { 
	lept1_3V = p4Mu1_.Vect();
	lept2_3V = p4Mu2_.Vect();
      }
      if(m_channelEM) { 
	lept1_3V = p4Ele1_.Vect();
	lept2_3V = p4Mu1_.Vect();
      }

      std::vector<float> theDR;
      if ( !_selection->getSwitch("isData") ) { 
	
	// cases where we have MC W->munu / W->enu
	for (int iMc=0; iMc<nMc; iMc++) { 
	  int mother = mothMc[iMc];
	  if ( abs(idMc[mother])==24 && (abs(idMc[iMc])==11 || abs(idMc[iMc])==13) ) {
	    TVector3 p3mc;
	    p3mc.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);
	    float dR1 = p3mc.DeltaR(lept1_3V);
	    float dR2 = p3mc.DeltaR(lept2_3V);
	    float dR;
	    if (dR1<dR2) theDR.push_back(dR1);
	    if (dR2<dR1) theDR.push_back(dR2);
	  }
	}
	
	// cases where the MC W->munu / W->enu are not found
	if (theDR.size()==0) {
	  theDR.push_back(1000.);
	  theDR.push_back(1000.);
	} 
	if (theDR.size()==1) theDR.push_back(1000.);
	
      } else {
	theDR.push_back(1000.);
	theDR.push_back(1000.);
      }

      // opposite charge      
      if (_selection->getSwitch("oppositeCharge") && (charge_[0]*charge_[1])>0 ) continue;
      if (m_channelEE) _counterEE.IncrVar("oppositeCharge");
      if (m_channelMM) _counterMM.IncrVar("oppositeCharge");
      if (m_channelEM) _counterEM.IncrVar("oppositeCharge");

      // reject j/psi, upsilon etc
      if ( (_selection->getSwitch("lowInvMass")) && (!_selection->passCut("lowInvMass",invMass_)) ) continue;
      if (m_channelEE) _counterEE.IncrVar("lowInvMass");
      if (m_channelMM) _counterMM.IncrVar("lowInvMass");
      if (m_channelEM) _counterEM.IncrVar("lowInvMass");

      // reject Z
      if ( (m_channelEE || m_channelMM) && 
	   (_selection->getSwitch("highInvMass")) && (_selection->passCut("highInvMass",invMass_)) ) continue;
      if (m_channelEE) _counterEE.IncrVar("highInvMass");
      if (m_channelMM) _counterMM.IncrVar("highInvMass");
      if (m_channelEM) _counterEM.IncrVar("highInvMass");

      // count PF jets with Et > lowEtJetAcc
      _theAK5PFJets     = GetCorrPFJets();            
      _theAK5BTagPFJets = GetBTagCorrPFJets();
      howManyPFJets = 0;
      JetCounter goodPFJets_counter(_theAK5PFJets);
      goodPFJets_counter.SetThresholds(_selection->getLowerCut("lowEtJetAcc"), _selection->getUpperCut("etaJetAcc"));
      goodPFJets_counter.SetDistance(_selection->getUpperCut("jetConeWidth"));
      goodPFJets_counter.SetPFJetId(Jet::loose);     
      goodPFJets_counter.SetBTagJets(_theAK5BTagPFJets);

      if (m_channelEE) {
	std::vector<TVector3> Lepton;
	Lepton.push_back(p4Ele1_.Vect());
	Lepton.push_back(p4Ele2_.Vect());
	goodPFJets_counter.SetParticlesToRemove(Lepton);
      } 
      if (m_channelMM) {
	std::vector<TVector3> Lepton;
	Lepton.push_back(p4Mu1_.Vect());
	Lepton.push_back(p4Mu2_.Vect());
	goodPFJets_counter.SetParticlesToRemove(Lepton);
      }
      if (m_channelEM) {
	std::vector<TVector3> Lepton;
	Lepton.push_back(p4Ele1_.Vect());
	Lepton.push_back(p4Mu1_.Vect());
	goodPFJets_counter.SetParticlesToRemove(Lepton);
      }
      goodPFJets_   = goodPFJets_counter.getGoodJets();
      howManyPFJets = goodPFJets_.size();

      // at least 2 reconstructed PF jets:
      if (_selection->getSwitch("jetNumber")) {
	if (!_selection->passCut("jetNumber",howManyPFJets))   continue;
      }

      // count b-tagged PF jets with 5 threshold for one reference algorithm
      nTagWPFJets[0] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 5.0); // VVLoose
      nTagWPFJets[1] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 4.0); // VLoose
      nTagWPFJets[2] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 3.3); // Loose
      nTagWPFJets[3] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 2.1); // Tight
      nTagWPFJets[4] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 1.9); // VTight

      // order jet by ET
      float highestEtJet = goodPFJets_[0].pt();
      float secondEtJet  = goodPFJets_[1].pt();
      // at least 2 reconstructed PF jets:
      // so far we applied the lowest threshold
      // here we apply the highest threshold as well
      if (_selection->getSwitch("jetNumber")) {
	if (!_selection->passCut("jetNumber",howManyPFJets))   continue;
      	if (!_selection->passCut("highEtJetAcc",highestEtJet)) continue;
      }
      if (m_channelEE) _counterEE.IncrVar("jetNumber");
      if (m_channelMM) _counterMM.IncrVar("jetNumber");
      if (m_channelEM) _counterEM.IncrVar("jetNumber");

      // significative MET
      if (_selection->getSwitch("met") && !_selection->passCut("met",p3Met_.Mag())) continue; 
      if (m_channelEE) _counterEE.IncrVar("met");
      if (m_channelMM) _counterMM.IncrVar("met");
      if (m_channelEM) _counterEM.IncrVar("met");

      // projected MET (see WW selection)
      if (_selection->getSwitch("projmet") && !_selection->passCut("projmet",projectedMet_)) continue; 
      if (m_channelEE) _counterEE.IncrVar("projmet");
      if (m_channelMM) _counterMM.IncrVar("projmet");
      if (m_channelEM) _counterEM.IncrVar("projmet");

      // full selection
      if (m_channelEE) _counterEE.IncrVar("fullSelection");      
      if (m_channelMM) _counterMM.IncrVar("fullSelection");      
      if (m_channelEM) _counterEM.IncrVar("fullSelection");      




      // for the two highest ET jets in selected events check:
      // 
      // 1) if the match or not a MC truth b
      int matchedAll  = 1000;
      int matchedJet1 = 1000;
      int matchedJet2 = 1000;
      int foundBq     = 1000;
      int foundBh     = 1000;
      if(!_selection->getSwitch("isData")) {
	foundBq = foundBquarks();
	foundBh = foundBhadrons();
	matchedAll = countAllBTagJets(goodPFJets_);
	bool bmatchedJet1 = countBTagJets(goodPFJets_[0]);
	bool bmatchedJet2 = countBTagJets(goodPFJets_[1]);
	if (bmatchedJet1)  matchedJet1 = 1;
	if (bmatchedJet2)  matchedJet2 = 1;
	if (!bmatchedJet1) matchedJet1 = 0;
	if (!bmatchedJet2) matchedJet2 = 0;
      }
      // 
      // 2) if they passed or not the btag cut
      std::vector<Jet> btaggedJets = goodPFJets_counter.getGoodBTaggedJets(bits::trackCountingHighEffBJetTags, 3.3); 
      int tagged1 = 0;
      int tagged2 = 0;
      for(unsigned int j=0; j<btaggedJets.size(); ++j) {
	TVector3 p3Jet = btaggedJets[j].Get3Vector();
	float deta1 = fabs(p3Jet.Eta() - goodPFJets_[0].eta());
	float dpt1  = fabs(p3Jet.Pt() - goodPFJets_[0].pt());
	float deta2 = fabs(p3Jet.Eta() - goodPFJets_[1].eta());
	float dpt2  = fabs(p3Jet.Pt() - goodPFJets_[1].pt());
	if(deta1<0.001 && dpt1<0.001) tagged1=1;
	if(deta2<0.001 && dpt2<0.001) tagged2=1;
      }


      // filling reduced tree
      myOutTree -> fillRunInfo(runNumber, lumiBlock, eventNumber);
      
      myOutTree -> fillVertexComp(dzEle, dzMu);
      
      myOutTree -> fillGeneral( m_channelEE, m_channelMM, m_channelEM, pt_[0], pt_[1], eta_[0], eta_[1], charge_[0], charge_[1], invMass_, howManyPFJets, p3Met_.Mag(), theDR[0], theDR[1]);

      myOutTree -> fillBTagJetMultiplicities( nTagWPFJets );

      myOutTree -> fillBTag(highestEtJet,secondEtJet,foundBq,foundBh,matchedAll,matchedJet1,matchedJet2,tagged1,tagged2);

      myOutTree -> store();
      
    } // event to analyze      

  } // loop over entries

}

// two highest pT 'good' electrons
std::pair<int,int> TopControlSample::getBestGoodElePair(std::vector<int> goodElectrons) {
  int theEle1=-1;
  int theEle2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;

  // cout << endl;
  // cout << "dentro getBestGoodElePair" << endl;
  // for(int ii=0; ii<nEle;ii++)  cout << "ele " << ii << ": " << sqrt(pxEle[ii]*pxEle[ii]+pyEle[ii]*pyEle[ii]) << " " << etaEle[ii] << endl;

  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
    float thisPt=pEle.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theEle2 = theEle1; theEle1 = eleIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theEle2 = eleIndex; }
  }

  // cout << "scelti: " << theEle1 << " " << theEle2 << endl;

  return make_pair(theEle1,theEle2);
}

// two highest pT 'good' muons
std::pair<int,int> TopControlSample::getBestGoodMuonPair(std::vector<int> goodMuons) {
  int theMu1=-1;
  int theMu2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1000.;
  for(int iMu=0;iMu<goodMuons.size();iMu++) {
    int muonIndex=goodMuons[iMu];
    TVector3 pMu(pxMuon[muonIndex],pyMuon[muonIndex],pzMuon[muonIndex]);
    float thisPt=pMu.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theMu2 = theMu1; theMu1 = muonIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theMu2 = muonIndex; }
  }
  return make_pair(theMu1,theMu2);
}


void TopControlSample::setKinematics2Ele(int theEle1, int theEle2) {

  // 3-4 vectors
  TVector3 p1, p2;
  p1.SetXYZ(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1]);
  p2.SetXYZ(pxEle[theEle2],pyEle[theEle2],pzEle[theEle2]);

  pT3Ele1_.SetXYZ(pxEle[theEle1],pyEle[theEle1],0.0);
  pT3Ele2_.SetXYZ(pxEle[theEle2],pyEle[theEle2],0.0);
  
  p4Ele1_.SetXYZT(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1],energyEle[theEle1]);
  p4Ele2_.SetXYZT(pxEle[theEle2],pyEle[theEle2],pzEle[theEle2],energyEle[theEle2]);

  // charge
  charge_[0] = chargeEle[theEle1];
  charge_[1] = chargeEle[theEle2];

  // eta
  eta_[0] = etaEle[theEle1];
  eta_[1] = etaEle[theEle2];
  
  // pt
  pt_[0] = p4Ele1_.Pt();
  pt_[1] = p4Ele2_.Pt();

  // PF missing transverse energy 
  p3Met_.SetXYZ(pxPFMet[0],pyPFMet[0],0.0);

  // W transverse mass using PF met
  WmT_ = sqrt(2 * p4Ele1_.Pt() * p3Met_.Mag() * (1-cos(pT3Ele1_.Angle(p3Met_))) );   

  // invariant mass
  invMass_ = (p4Ele1_ + p4Ele2_).M();       

  // projected met
  projectedMet_ = GetProjectedMet(p4Ele1_.Vect(),p4Ele2_.Vect());
}

void TopControlSample::setKinematics2Mu(int theMu1, int theMu2) {

  // 3-4 vectors
  TVector3 p1, p2;
  p1.SetXYZ(pxMuon[theMu1],pyMuon[theMu1],pzMuon[theMu1]);
  p2.SetXYZ(pxMuon[theMu2],pyMuon[theMu2],pzMuon[theMu2]);

  pT3Mu1_.SetXYZ(pxMuon[theMu1],pyMuon[theMu1],0.0);
  pT3Mu2_.SetXYZ(pxMuon[theMu2],pyMuon[theMu2],0.0);
  
  p4Mu1_.SetXYZT(pxMuon[theMu1],pyMuon[theMu1],pzMuon[theMu1],energyMuon[theMu1]);
  p4Mu2_.SetXYZT(pxMuon[theMu2],pyMuon[theMu2],pzMuon[theMu2],energyMuon[theMu2]);

  // charge
  charge_[0] = chargeMuon[theMu1];
  charge_[1] = chargeMuon[theMu2];

  // eta
  eta_[0] = etaMuon[theMu1];
  eta_[1] = etaMuon[theMu2];
  
  // pt
  pt_[0] = p4Mu1_.Pt();
  pt_[1] = p4Mu2_.Pt();

  // PF missing transverse energy 
  p3Met_.SetXYZ(pxPFMet[0],pyPFMet[0],0.0);

  // W transverse mass using PF met
  WmT_ = sqrt(2 * p4Mu1_.Pt() * p3Met_.Mag() * (1-cos(pT3Mu1_.Angle(p3Met_))) );   

  // invariant mass
  invMass_ = (p4Mu1_ + p4Mu2_).M();       

  // projected met
  projectedMet_ = GetProjectedMet(p4Mu1_.Vect(),p4Mu2_.Vect());
}

void TopControlSample::setKinematics1Ele1Mu(int theEle1, int theMu1) {

  // 3-4 vectors
  TVector3 p1, p2;
  p1.SetXYZ(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1]);
  p2.SetXYZ(pxMuon[theMu1],pyMuon[theMu1],pzMuon[theMu1]);

  pT3Ele1_.SetXYZ(pxEle[theEle1],pyEle[theEle1],0.0);
  pT3Mu1_.SetXYZ(pxMuon[theMu1],pyMuon[theMu1],0.0);
  
  p4Ele1_.SetXYZT(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1],energyEle[theEle1]);
  p4Mu1_.SetXYZT(pxMuon[theMu1],pyMuon[theMu1],pzMuon[theMu1],energyMuon[theMu1]);

  // charge
  charge_[0] = chargeEle[theEle1];
  charge_[1] = chargeMuon[theMu1];

  // eta
  eta_[0] = etaEle[theEle1];
  eta_[1] = etaMuon[theMu1];
  
  // pt
  pt_[0] = p4Ele1_.Pt();
  pt_[1] = p4Mu1_.Pt();

  // PF missing transverse energy 
  p3Met_.SetXYZ(pxPFMet[0],pyPFMet[0],0.0);

  // W transverse mass using PF met and the highest Et lepton
  float elePt = p4Ele1_.Pt();
  float muPt  = p4Mu1_.Pt();

  if (elePt>muPt) 
    WmT_ = sqrt(2 * p4Ele1_.Pt() * p3Met_.Mag() * (1-cos(pT3Ele1_.Angle(p3Met_))) );   
  
  if (elePt<muPt) 
    WmT_ = sqrt(2 * p4Mu1_.Pt() * p3Met_.Mag() * (1-cos(pT3Mu1_.Angle(p3Met_))) );   

  // invariant mass
  invMass_ = (p4Ele1_ + p4Mu1_).M();       

  // projected met
  projectedMet_ = GetProjectedMet(p4Ele1_.Vect(),p4Mu1_.Vect());
}


void TopControlSample::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];

  deta    = deltaEtaAtVtxEle[eleIndex];
  dphiin  = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];

  fbrem  = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxSC[sc];
      e4SwissCross = e4SwissCrossSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagSC[sc];
      seedTime = timeSC[sc];
      seedChi2 = chi2SC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->applyElectronIDOnPFlowElectrons(true);
  selector->SetHOverE( HoE );
  selector->SetS9S25( s9s25 );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetDPhiOut( dphiout );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetSigmaPhiPhi( spp );
  selector->SetEOverPout( eopout );
  selector->SetEOverPin( eop );
  selector->SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  selector->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  selector->SetEcalIsolation( ecalIsolation[eleIndex] );
  selector->SetTrkIsolation( trackerIsolation[eleIndex] );
  selector->SetHcalIsolation( hcalIsolation[eleIndex] );
  selector->SetCombinedIsolation( combinedIsolation[eleIndex] );
  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
  selector->m_cleaner->SetE1(e1);
  selector->m_cleaner->SetE4SwissCross(e4SwissCross);
  selector->m_cleaner->SetFiducialFlag(fidFlagSC);
  selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
  selector->m_cleaner->SetSeedTime(seedTime);
  selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();
}

void TopControlSample::isMuonID(int muonIndex, bool *muonIdOutput) {

  *muonIdOutput = true;

  Utils anaUtils; 

  // tracker and global muon
  bool flag = anaUtils.muonIdVal(muonIdMuon[muonIndex],AllGlobalMuons) && anaUtils.muonIdVal(muonIdMuon[muonIndex],AllTrackerMuons);
  
  // the following cuts are based on KF and global muon track. So if the cut above has failed, return here
  if(!flag) {
    *muonIdOutput = false;
    return;
  }

  int track = trackIndexMuon[muonIndex];
  if(trackValidHitsTrack[track]<=10) *muonIdOutput = false;

  int globalMuonTrack = combinedTrackIndexMuon[muonIndex];
  if(trackNormalizedChi2GlobalMuonTrack[globalMuonTrack] >= 10) *muonIdOutput = false;
  if(trackValidHitsGlobalMuonTrack[globalMuonTrack] == 0) *muonIdOutput = false;

  float dxy = fabs(trackDxyPV(PVxPV[hardestPV], PVyPV[hardestPV], PVzPV[hardestPV], 
                              trackVxTrack[track], trackVyTrack[track], trackVzTrack[track], 
                              pxTrack[track], pyTrack[track], pzTrack[track]));
  if(dxy > 0.020) *muonIdOutput = false;
}

double TopControlSample::trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}

void TopControlSample::ConfigureSelection(Selection* _selection, Counters *_counter) {
  
  _selection->addSwitch("isData");
  _selection->addSwitch("goodRunLS");
  _selection->addSwitch("mcTruth");
  _selection->addSwitch("trigger");
  _selection->addSwitch("goodTracks");
  _selection->addSwitch("leptonAcceptance");  
  _selection->addSwitch("leptonId");
  _selection->addSwitch("leptonIsol");
  _selection->addSwitch("convRej");
  _selection->addSwitch("further");
  _selection->addSwitch("vertexCompatibility");
  _selection->addSwitch("oppositeCharge");
  _selection->addSwitch("lowInvMass");
  _selection->addSwitch("highInvMass");
  _selection->addSwitch("jetNumber");
  _selection->addSwitch("met");
  _selection->addSwitch("projmet");

  _selection->addCut("goodTracksRatio");
  _selection->addCut("pvNdf");
  _selection->addCut("pvDz");
  _selection->addCut("pvRho");
  _selection->addCut("etaElectronAcc");
  _selection->addCut("ptElectronAcc");
  _selection->addCut("etaMuonAcc");
  _selection->addCut("ptMuonAcc");
  _selection->addCut("isolMuon");
  _selection->addCut("dzVertexEle");
  _selection->addCut("dzVertexMu");
  _selection->addCut("lowEtJetAcc");
  _selection->addCut("highEtJetAcc");
  _selection->addCut("etaJetAcc");
  _selection->addCut("jetConeWidth");
  _selection->addCut("lowInvMass");
  _selection->addCut("highInvMass");
  _selection->addCut("jetNumber");
  _selection->addCut("met");
  _selection->addCut("projmet");

  _selection->addStringParameter("electronIDType");

  _selection->summary();

  _counter->SetTitle("COUNTER");
  _counter->AddVar("event");
  _counter->AddVar("mcTruth");
  _counter->AddVar("trigger");
  _counter->AddVar("goodTracks");
  _counter->AddVar("primaryVertex");
  _counter->AddVar("leptonAcceptance");  
  _counter->AddVar("leptonId");
  _counter->AddVar("leptonIsol");
  _counter->AddVar("convRej");
  _counter->AddVar("further");
  _counter->AddVar("vertexCompatibility");
  _counter->AddVar("recoChannel");
  _counter->AddVar("oppositeCharge");
  _counter->AddVar("lowInvMass");
  _counter->AddVar("highInvMass");
  _counter->AddVar("jetNumber");
  _counter->AddVar("met");
  _counter->AddVar("projmet");
  _counter->AddVar("fullSelection");
}

void TopControlSample::displayRecoEfficiencies(Counters _counter) {

  _counter.Draw();
  if(_selection->getSwitch("isData")) {
    _counter.Draw("trigger","event");
  } else {
    _counter.Draw("mcTruth","event");
    _counter.Draw("trigger","mcTruth");
  }
  _counter.Draw("goodTracks","trigger");
  _counter.Draw("primaryVertex","goodTracks");
  _counter.Draw("leptonAcceptance","primaryVertex");
  _counter.Draw("leptonId","leptonAcceptance");
  _counter.Draw("leptonIsol","leptonId");
  _counter.Draw("convRej","leptonIsol");
  _counter.Draw("further","convRej");
  _counter.Draw("vertexCompatibility","further");
  _counter.Draw("recoChannel","vertexCompatibility");
  _counter.Draw("oppositeCharge","recoChannel");
  _counter.Draw("lowInvMass","oppositeCharge");
  _counter.Draw("highInvMass","lowInvMass");
  _counter.Draw("jetNumber","highInvMass");
  _counter.Draw("met","jetNumber");
  _counter.Draw("projmet","met");
  _counter.Draw("fullSelection","event");
}

void TopControlSample::displayEfficiencies() {

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "EE: " << std::endl;
  displayRecoEfficiencies(_counterEE);

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "MM: " << std::endl;
  displayRecoEfficiencies(_counterMM);

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "EM: " << std::endl;
  displayRecoEfficiencies(_counterEM);

  std::cout << endl;
  std::cout << endl;
  std::cout << "ElectronID selection: " << std::endl;
  EgammaTightID_NC.displayEfficiencies();
  std::cout << endl;

  // save into ROOT files the events as a function of jet multiplicity 
  sprintf(namefile,"%sTopCounters.root",m_prefix);
  _counterEE.Save(namefile,"recreate");
  _counterMM.Save(namefile,"update");
  _counterEM.Save(namefile,"update");
}

float TopControlSample::GetProjectedMet(TVector3 p1, TVector3 p2) {

  float deltaPhi1 = fabs(p1.DeltaPhi(p3Met_));
  float deltaPhi2 = fabs(p2.DeltaPhi(p3Met_));
  float deltaphi = TMath::Min(deltaPhi1,deltaPhi2);
  if(deltaphi<TMath::Pi()/2.) return p3Met_.Mag() * sin(deltaphi);
  else return p3Met_.Mag();
}

// to check if the reco jet matched a mc truth b
bool TopControlSample::countBTagJets(Jet myJet) {

  bool matched = false;

  std::vector<TVector3> theBquarks;
  for (int iMc=0; iMc<nMc; iMc++) { 
    int mother = mothMc[iMc];
    if ( abs(idMc[iMc])==5 && abs(idMc[mother])==6 ) {
      TVector3 p3;
      p3.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);
      theBquarks.push_back(p3);
    }
  }
  
  TVector3 p3Jet = myJet.Get3Vector();
  for(unsigned int b=0; b<theBquarks.size(); b++) {
    float dR = p3Jet.DeltaR(theBquarks[b]);
    if(dR<0.3) {
      matched = true;
      break;
    }
  }
  
  return matched;
}

// find all the b's in the event MC truth
int TopControlSample::foundBquarks() {

  int foundBs = 0;
  std::vector<TVector3> theBquarks;
  for (int iMc=0; iMc<nMc; iMc++) { 
    int mother = mothMc[iMc];
    if ( abs(idMc[iMc])==5 && abs(idMc[mother])==6 ) {
      TVector3 p3;
      p3.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);
      theBquarks.push_back(p3);
    }
  }
  foundBs=theBquarks.size();
  return foundBs;
}

// find all the b hadrons in the event MC truth
int TopControlSample::foundBhadrons() {

  int foundBs = 0;
  std::vector<TVector3> theBhadrons;
  for (int iMc=0; iMc<nMc; iMc++) { 
    if ( (int(abs(idMc[iMc]))/1000) == 5 || // b-baryon
	 (int(abs(idMc[iMc]))/100)%100 == 5 || // b mesons
	 abs(idMc[iMc]) == 5 ) { // b quark
      TVector3 p3BHad;
      p3BHad.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);
      theBhadrons.push_back(p3BHad);
    }
  }

  foundBs=theBhadrons.size();
  return foundBs;
}

int TopControlSample::countAllBTagJets(std::vector<Jet> jets) {
  
  int matched = 0;
  
  std::vector<TVector3> theBquarks;
  for (int iMc=0; iMc<nMc; iMc++) { 
    int mother = mothMc[iMc];
    if ( abs(idMc[iMc])==5 && abs(idMc[mother])==6 ) {
      TVector3 p3;
      p3.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);
      theBquarks.push_back(p3);
    }
  }

  for(unsigned int j=0; j<jets.size(); ++j) {
    TVector3 p3Jet = jets[j].Get3Vector();
    for(unsigned int b=0; b<theBquarks.size(); b++) {
      float dR = p3Jet.DeltaR(theBquarks[b]);
      // cout << "jet " << j << ", pT = " << p3Jet.Pt() << ", eta = " << p3Jet.Eta() << ", b eta = " << theBquarks[b].Eta() << ", dR = " << dR << endl;
      if(dR<0.3) {
        matched++;
        break;
      }
    }
  }

  return matched;
}
