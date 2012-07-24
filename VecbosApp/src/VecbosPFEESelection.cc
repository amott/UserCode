#include <string>
#include <iostream>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/Utils.hh"
#include "include/VecbosPFEESelection.hh"
#include "include/JetCounter.hh"
   

using namespace bits;

VecbosPFEESelection::VecbosPFEESelection(TTree *tree) 
  : Vecbos(tree) {

  // default do not check on mc-truth
  m_signal = all;

  // common kinematic selections
  std::string theConfigDir       = "config/vecbos/";
  std::string fileCutsCommon     = theConfigDir + "CommonSelectionCutsPFElectrons.txt";
  std::string fileSwitchesCommon = theConfigDir + "CommonSelectionSwitchesPFElectrons.txt";
  
  _commonSel = new Selection(fileCutsCommon,fileSwitchesCommon);
  ConfigCommonSelections(_commonSel);

  // To read good run list!
  std::cout << "[GoodRunLS]::goodRunLS is " << _commonSel->getSwitch("goodRunLS") << " isData is " <<  _commonSel->getSwitch("isData") << std::endl;
  if (_commonSel->getSwitch("goodRunLS") && _commonSel->getSwitch("isData")) {
    std::string goodRunGiasoneFile       = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  // Z -> ee specific selections
  std::string fileCutsZee     = theConfigDir + "ZeeCuts.txt";
  std::string fileSwitchesZee = theConfigDir + "ZeeSwitches.txt";
  _zeeSel = new Selection(fileCutsZee,fileSwitchesZee);
  _zeePCounter      = &_zeeCounter;
  _zeePCounterZjets = &_zeeCounterZjets;
  _zeePCounterWjets = &_zeeCounterWjets;
  _zeePCounterTTbar = &_zeeCounterTTbar;
  ConfigZeeSelections(_zeeSel, _zeePCounter);
  _zeePCounter->SetTitle("Z+jets counter");
  ConfigZeeSelections(_zeeSel, _zeePCounterZjets); 
  _zeePCounterZjets->SetTitle("Z+jets counter Z events");
  ConfigZeeSelections(_zeeSel, _zeePCounterWjets); 
  _zeePCounterWjets->SetTitle("Z+jets counter W events");
  ConfigZeeSelections(_zeeSel, _zeePCounterTTbar); 
  _zeePCounterTTbar->SetTitle("Z+jets counter TTbar events");
  for(int jetbin=0; jetbin<6; jetbin++) {
    ConfigZeeSelections(_zeeSel, &m_zeeJetbinCounter[jetbin]);
    ConfigZeeSelections(_zeeSel, &m_zeeRecoJetbinCounter[jetbin]);
    ConfigZeeSelections(_zeeSel, &m_zeePFJetbinCounter[jetbin]);
    ConfigZeeSelections(_zeeSel, &m_zeeRecoPFJetbinCounter[jetbin]);

    char counterTitle[200];
    // with MC truth electrons subtraction
    sprintf(counterTitle,"Z+%djets counter",jetbin);
    m_zeeJetbinCounter[jetbin].SetTitle(counterTitle);
    // with reco electrons subtraction
    sprintf(counterTitle,"Z+%drecojets counter",jetbin);
    m_zeeRecoJetbinCounter[jetbin].SetTitle(counterTitle);
    // with MC truth electrons subtraction
    sprintf(counterTitle,"Z+%dPFjets counter",jetbin);
    m_zeePFJetbinCounter[jetbin].SetTitle(counterTitle);
    // with reco electrons subtraction
    sprintf(counterTitle,"Z+%drecoPFjets counter",jetbin);
    m_zeeRecoPFJetbinCounter[jetbin].SetTitle(counterTitle);
  }

  // W -> enu specific selections
  std::string fileCutsWenu     = theConfigDir + "WenuPFCuts.txt";
  std::string fileSwitchesWenu = theConfigDir + "WenuPFSwitches.txt";
  _wenuSel = new Selection(fileCutsWenu,fileSwitchesWenu);
  _wenuPCounter      = &_wenuCounter;
  _wenuPCounterZjets = &_wenuCounterZjets;
  _wenuPCounterWjets = &_wenuCounterWjets;
  _wenuPCounterTTbar = &_wenuCounterTTbar;
  ConfigWenuSelections(_wenuSel, _wenuPCounter);
  _wenuPCounter->SetTitle("W+jets counter");
  ConfigWenuSelections(_wenuSel, _wenuPCounterZjets); 
  _wenuPCounterZjets->SetTitle("W+jets counter Z events");
  ConfigWenuSelections(_wenuSel, _wenuPCounterWjets); 
  _wenuPCounterWjets->SetTitle("W+jets counter W events");
  ConfigWenuSelections(_wenuSel, _wenuPCounterTTbar); 
  _wenuPCounterTTbar->SetTitle("W+jets counter TTbar events");
  for(int jetbin=0; jetbin<6; jetbin++) {
    ConfigWenuSelections(_wenuSel, &m_wenuJetbinCounter[jetbin]);
    ConfigWenuSelections(_wenuSel, &m_wenuRecoJetbinCounter[jetbin]);
    ConfigWenuSelections(_wenuSel, &m_wenuPFJetbinCounter[jetbin]);
    ConfigWenuSelections(_wenuSel, &m_wenuRecoPFJetbinCounter[jetbin]);
    char counterTitle[200];
    // with MC truth electrons subtraction
    sprintf(counterTitle,"W+%djets counter",jetbin);
    m_wenuJetbinCounter[jetbin].SetTitle(counterTitle);
    // with reco electrons subtraction
    sprintf(counterTitle,"W+%drecojets counter",jetbin);
    m_wenuRecoJetbinCounter[jetbin].SetTitle(counterTitle);
    // with MC truth electrons subtraction
    sprintf(counterTitle,"W+%dPFjets counter",jetbin);
    m_wenuPFJetbinCounter[jetbin].SetTitle(counterTitle);
    // with reco electrons subtraction
    sprintf(counterTitle,"W+%drecoPFjets counter",jetbin);
    m_wenuRecoPFJetbinCounter[jetbin].SetTitle(counterTitle);
  }

  m_prefix = "";

  // counters to compute number of jets at the beginning/end of the selection
  for(int ii=0; ii<6; ii++){
    m_wInitJets[ii] = 0.;
    m_zInitJets[ii] = 0.;
    m_wEndJets[ii]  = 0.;
    m_zEndJets[ii]  = 0.;
  }

  WToENuDecay = ZToEEDecay = 0;
  isData_ = _commonSel->getSwitch("isData");
  if(isData_) mcevent.SetData(true);
}

VecbosPFEESelection::~VecbosPFEESelection(){
  
  // closing/saving reduced trees
  myOutTree_Zee   -> save();
  myOutTree_Wenu  -> save();
  if(_commonSel->getSwitch("dumpEleID")) myOutEleIdIsolPFTree_Wenu -> save();   

  // deleting selections
  delete _commonSel;
  delete _zeeSel;
  delete _wenuSel;  
}

// loop over events - real analysis
void VecbosPFEESelection::Loop() {

  _verbose=false;
  if(fChain == 0) return;

  // kinematics reduced tree
  sprintf(namefile,"%sVecBosOutput-out-Zee.root",m_prefix);
  myOutTree_Zee  = new RedVecbosPFTree(namefile);
  sprintf(namefile,"%sVecBosOutput-out-Wenu.root",m_prefix);
  myOutTree_Wenu = new RedVecbosPFTree(namefile);
  myOutTree_Zee ->addMcTruthInfos();
  myOutTree_Wenu->addMcTruthInfos();
  myOutTree_Zee ->addElectronPFInfos();
  myOutTree_Wenu->addElectronPFInfos();

  // reduced tree for electron ID and isolation studies
  if(_commonSel->getSwitch("dumpEleID")) {
    sprintf(namefile,"%sVecBosOutput-EleID.root",m_prefix);
    myOutEleIdIsolPFTree_Wenu = new RedEleIdIsolPFTree(namefile);
  }

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Total number of entries in the chain = " << nentries << std::endl;
  int maxEvents = nentries;
  if(_commonSel->getSwitch("eventRange")) {
    maxEvents = (int) _commonSel->getUpperCut("eventRange");
    cout << "WARNING! switch eventRange ON! Will run only on the first " << maxEvents << " events!" << endl;
  }

  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  for (Long64_t jentry=0; jentry<maxEvents;jentry++) {

    bool eventToAnalyze = true;

    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    if ( eventToAnalyze ) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE 
      reloadTriggerMask();

      // Good Run selection
      if (isData_ && _commonSel->getSwitch("goodRunLS") && !isGoodRunLS()) {
	if ( lastRun!= runNumber || lastLumi != lumiBlock) {
	  lastRun = runNumber;
	  lastLumi = lumiBlock;
	  std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
	}
	continue;
      }
      if (isData_ && _commonSel->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
      }

      // used in the past for the soups
      float weight  = 1;

      int indexMcEleWToENu = -1;
      if ( !isData_ ) { 

        mcevent.LoadDecay(nMc,idMc,mothMc);
        mcevent.LoadMomentum(pMc,energyMc,thetaMc,phiMc);

        // check if the decay is a W->enu prompt
        indexMcEleWToENu = mcevent.indexEleWPrompt();
        WToENuDecay = (indexMcEleWToENu >-1 ) ? 1 : 0;

        // check if the decay is a Z->ee prompt
        int idx1Z = mcevent.indicesEleZPrompt().first;
        int idx2Z = mcevent.indicesEleZPrompt().second;
        ZToEEDecay = (idx1Z > 0 && idx2Z > 0 );

        // evaluate the pthat of the photon + jet event (only for photon+j events)
        photonj_pthat = photonPt();
      }

      // trigger 
      Utils anaUtils;
      bool passedHLT = hasPassedHLT();

      // take from the event the gen jets
      if(!isData_) _theAK5GenJets = GetGenJets();

      // take from the event the reconstructed jets
      _theAK5CaloJets = GetUncorrJets();
      _theAK5BTagCaloJets = GetBTagCorrJets();

      if(_theAK5CaloJets.size() != _theAK5BTagCaloJets.size()) {
        cout << "NASTY ERROR: _theAK5CaloJets.size() != _theAK5BTagCaloJets.size()" << endl;
        return;
      }

      // count jets removing MC truth electrons, for efficiency normalization bin-by-bin
      nwjets_mc = nwpfjets_mc = nwgenjets_mc = nzjets_mc = nzpfjets_mc = nzgenjets_mc = 0;
      if(!isData_) {
        std::vector<TVector3> mcElectrons;

        if(mcevent.indexEleWPrompt() > -1) { 
          mcElectrons.push_back( mcevent.p4ElectronWenuPrompt().Vect() );

          // count calojets
          JetCounter goodJetsMC_counter(_theAK5CaloJets);
          goodJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
          goodJetsMC_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
          goodJetsMC_counter.SetParticlesToRemove(mcElectrons);
          nwjets_mc = goodJetsMC_counter.numGoodJets();

          // count PF jets
          JetCounter goodPFJetsMC_counter(_theAK5PFJets);
          goodPFJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
          goodPFJetsMC_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
          goodPFJetsMC_counter.SetParticlesToRemove(mcElectrons);
          nwpfjets_mc = goodJetsMC_counter.numGoodJets();

          if(!isData_) {
            // count GenJets
            JetCounter genJets_counter(_theAK5GenJets);
            genJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
            genJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
            genJets_counter.SetParticlesToRemove(mcElectrons);
            nwgenjets_mc = goodJetsMC_counter.numGoodJets();
          }
        }

        mcElectrons.clear();
        if(mcevent.indicesEleZPrompt().first > -1 && mcevent.indicesEleZPrompt().second > -1) { 
          mcElectrons.push_back( (mcevent.p4ElectronZeePrompt().first).Vect() );
          mcElectrons.push_back( (mcevent.p4ElectronZeePrompt().second).Vect() );

          // count calojets
          JetCounter goodJetsMC_counter(_theAK5CaloJets);
          goodJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
          goodJetsMC_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
          goodJetsMC_counter.SetParticlesToRemove(mcElectrons);
          nzjets_mc = goodJetsMC_counter.numGoodJets();

          // count PF jets
          JetCounter goodPFJetsMC_counter(_theAK5PFJets);
          goodPFJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
          goodPFJetsMC_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
          goodPFJetsMC_counter.SetParticlesToRemove(mcElectrons);
          nwpfjets_mc = goodJetsMC_counter.numGoodJets();

          if(!isData_) {
            // count GenJets
            JetCounter genJets_counter(_theAK5GenJets);
            genJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
            genJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
            genJets_counter.SetParticlesToRemove(mcElectrons);
            nzgenjets_mc = goodJetsMC_counter.numGoodJets();
          }
        }
      }

      // electrons
      std::vector<int> acceptElectrons; 
      std::vector<int> idTightElectrons,     idLooseElectrons; 
      std::vector<int> isolTightElectrons,   isolLooseElectrons; 
      std::vector<int> convTightElectrons,   convLooseElectrons; 
      std::vector<int> vertexTightElectrons, vertexLooseElectrons;
      
      // for isolation
      chargedIsolation.clear();
      photonIsolation.clear();
      neutralIsolation.clear();
      combinedIsolation.clear();
      
      // reconstructed PF electrons
      int howManyRecoEle = nPFEle;
      
      // PF electrons in the acceptance  
      for(int ii=0; ii<nPFEle; ii++) {             
        bool isGoodEle = true;
	
	// eta cut
	bool isInEcalFiducial = false;
        if(!_commonSel->passCut("etaElectronAcc",etaPFEle[ii])) isGoodEle = false;      
	if ( (fabs(etaPFEle[ii]) < 1.4442) || ( fabs(etaPFEle[ii])>1.560 && fabs(etaPFEle[ii])<2.5 ) ) isInEcalFiducial = true;
	if (!isInEcalFiducial) isGoodEle = false;

	// pt cut
        TVector3 pLepton(pxPFEle[ii],pyPFEle[ii],pzPFEle[ii]);
        float thisPt = pLepton.Pt();
        if(!_commonSel->passCut("ptElectronAcc",thisPt) ) isGoodEle = false; 
   
        if (isGoodEle) acceptElectrons.push_back(ii);
	
        // assign to each electron the three absolute PF isolations  

	// chiara: usato fino ad ora
	chargedIsolation.push_back( chIso04vetoPFEle[ii] ); // / pLepton.Pt() );
        photonIsolation.push_back ( phIso04vetoPFEle[ii] ); // / pLepton.Pt() );
        neutralIsolation.push_back( nhIso04vetoPFEle[ii] ); // / pLepton.Pt() );

	// chiara: dopo riottimizzazione con new samples
	// chargedIsolation.push_back( chIso03noVetoPFEle[ii] ); 
        // photonIsolation.push_back ( phIso03vetoPFEle[ii] );   
        // neutralIsolation.push_back( nhIso03noVetoPFEle[ii] ); 

	// chiara: test x confronto con florian
        // chargedIsolation.push_back( chIso04noVetoNVCPFEle[ii] ); // / pLepton.Pt() );
        // photonIsolation.push_back ( phIso04noVetoPFEle[ii] ); // / pLepton.Pt() );
        // neutralIsolation.push_back( nhIso04noVetoPFEle[ii] ); // / pLepton.Pt() );

        // and the combined isolation

	// relative
	// chargedIsolRel = chIso04vetoPFEle[ii] / pLepton.Pt();
	// photonIsolRel  = phIso04vetoPFEle[ii] / pLepton.Pt();
	// neutralIsolRel = nhIso04vetoPFEle[ii] / pLepton.Pt();

	// mix...
	// float myMaxPt = pLepton.Pt();
	// if (myMaxPt<30) myMaxPt = 30.;
        // chargedIsolRel = chIso04vetoPFEle[ii] / myMaxPt;
        // photonIsolRel  = phIso04vetoPFEle[ii] / myMaxPt;
        // neutralIsolRel = nhIso04vetoPFEle[ii] / myMaxPt;

	// absolute: chiara, usato fino ad ora
	chargedIsolRel = chIso04vetoPFEle[ii];
	photonIsolRel  = phIso04vetoPFEle[ii];
	neutralIsolRel = nhIso04vetoPFEle[ii];

	// absolute: chiara, dopo riottimizzazione con new samples
	// chargedIsolRel = chIso03noVetoPFEle[ii];
	// photonIsolRel  = phIso03vetoPFEle[ii];
	// neutralIsolRel = nhIso03noVetoPFEle[ii];

	// chiara: test x confronto con florian
	// chargedIsolRel = chIso04noVetoNVCPFEle[ii];
	// photonIsolRel  = phIso04noVetoPFEle[ii];
	// neutralIsolRel = nhIso04noVetoPFEle[ii];

	float bestCombination = chargedIsolRel + photonIsolRel + neutralIsolRel;
	combinedIsolation.push_back( bestCombination );
      }
      int howManyAccEle = acceptElectrons.size();
           
      // to study eleid/isolation/conversion optimization: distributions for the highest pt candidate and selection with the tight hypothesis 
      std::pair<int,int> bestElectronAccept = getBestGoodElePair(acceptElectrons);
      int theEleAccep1(bestElectronAccept.first);
      if (theEleAccep1>-1) {
	if (_commonSel->getSwitch("dumpEleID")) { 
	  bool isBarrel = false;
	  if( fabs(etaPFEle[theEleAccep1])<1.476 ) isBarrel = true;
	  int passedEleId       = 1;
	  int passedIsolation   = 1;
	  int passedConversions = 1;
	  if (!isEleID(theEleAccep1,1)) passedEleId = 0;
	  if ( isBarrel && !_commonSel->passCut("combinedIsolationTightEB",combinedIsolation[theEleAccep1]) ) passedIsolation = 0;
	  if (!isBarrel && !_commonSel->passCut("combinedIsolationTightEE",combinedIsolation[theEleAccep1]) ) passedIsolation = 0;
	  int thisMatchedTrack   = trackIndexPFEle[theEleAccep1];
	  int thisExpInnerLayers = expInnerLayersTrack[thisMatchedTrack];
	  if (!_commonSel->passCut("conversionRejection",thisExpInnerLayers)) passedConversions = 0;
	  FillEleIdIsolPFTree(theEleAccep1, passedEleId, passedIsolation, passedConversions);
	}
      }
            
      // 'tight' identified electrons 
      for(int ii=0;ii<howManyAccEle;ii++) {      
        bool isGoodEle = true;
        if (_commonSel->getSwitch("eleId") && !isEleID(acceptElectrons[ii],1)) isGoodEle = false;    
        if (isGoodEle) idTightElectrons.push_back(acceptElectrons[ii]);            
      }
      int howManyIdTightEle = idTightElectrons.size();
      
      // loose identified electrons, if needed 
      int howManyIdLooseEle = 0;
      if (_commonSel->getSwitch("asymmetricElectrons")) {   
        for(int ii=0;ii<howManyAccEle;ii++) {      
          int accIndex=acceptElectrons[ii];
          if ( (_commonSel->getSwitch("eleId") && isEleID(acceptElectrons[ii],0)) || (!_commonSel->getSwitch("eleId")) ) {       
	    bool alsoTight = false;
            for(int jj=0;jj<howManyIdTightEle;jj++) {      
              int tightIdEleIndex=idTightElectrons[jj];
              if (tightIdEleIndex==accIndex) alsoTight = true; 
            }
            if(!alsoTight) idLooseElectrons.push_back(accIndex);
          }
        }
        howManyIdLooseEle = idLooseElectrons.size();
      }
      int howManyIdEle = howManyIdTightEle + howManyIdLooseEle;

      // (combined) tight optimized, tight isolated electrons
      for(int ii=0;ii<howManyIdTightEle;ii++) {      
        bool isGoodEle = true;
        int idEleIndex = idTightElectrons[ii];
	if( _commonSel->getSwitch("combinedIsolationTightEB") ) {
	  if( fabs(etaPFEle[idEleIndex])<1.476 && !_commonSel->passCut("combinedIsolationTightEB",combinedIsolation[idEleIndex]) ) isGoodEle = false;
	}
	if( _commonSel->getSwitch("combinedIsolationTightEE") ) {
	  if( fabs(etaPFEle[idEleIndex])>1.476 && !_commonSel->passCut("combinedIsolationTightEE",combinedIsolation[idEleIndex]) ) isGoodEle = false;
	}
	if (isGoodEle) isolTightElectrons.push_back(idEleIndex); 
      }
      int howManyIsolTightEles = isolTightElectrons.size();

      // loose optimized, loose (combined) isolated electrons if needed
      int howManyIsolLooseEles = 0;
      if (_commonSel->getSwitch("asymmetricElectrons")) {   
        for(int ii=0;ii<howManyIdLooseEle;ii++) {      
          bool isGoodEle = true;
          int idEleIndex=idLooseElectrons[ii];
	  if( _commonSel->getSwitch("combinedIsolationTightEB") ) {
	    if( fabs(etaPFEle[idEleIndex])<1.476 &&
		!_commonSel->passCut("combinedIsolationLooseEB",combinedIsolation[idEleIndex]) ) isGoodEle = false;
	  }
	  if( _commonSel->getSwitch("combinedIsolationTightEE") ) {
	    if( fabs(etaPFEle[idEleIndex])>1.476 &&
		!_commonSel->passCut("combinedIsolationLooseEE",combinedIsolation[idEleIndex]) ) isGoodEle = false;
	  }
          if (isGoodEle) isolLooseElectrons.push_back(idEleIndex); 
        }
        howManyIsolLooseEles = isolLooseElectrons.size();
      }

      // conversion rejection: tight electrons    
      for(int ii=0;ii<howManyIsolTightEles;ii++) {      
        bool isGoodEle      = true;
        int convEleIndex    = isolTightElectrons[ii];
	int theMatchedTrack = trackIndexPFEle[convEleIndex];
	int expInnerLayers  = expInnerLayersTrack[theMatchedTrack]; 
	if( _commonSel->getSwitch("conversionRejection") ) {  
	  if(!_commonSel->passCut("conversionRejection",expInnerLayers)) isGoodEle = false;   
	}
	if (isGoodEle) convTightElectrons.push_back(convEleIndex);   
      }
      int howManyConvTightEles = convTightElectrons.size();
      
      // and loose electrons
      int howManyConvLooseEles = 0;
      if (_commonSel->getSwitch("asymmetricElectrons")) {   
        for(int ii=0;ii<howManyIsolLooseEles;ii++) {      
          bool isGoodEle = true;
          int convEleIndex=isolLooseElectrons[ii];
	  int theMatchedTrack = trackIndexPFEle[convEleIndex];
	  int expInnerLayers  = expInnerLayersTrack[theMatchedTrack]; 
	  if( _commonSel->getSwitch("conversionRejection") ) {
	    if(!_commonSel->passCut("conversionRejection",expInnerLayers)) isGoodEle = false;   
	  }
          if (isGoodEle) convLooseElectrons.push_back(convEleIndex); 
        }
        howManyConvLooseEles = convLooseElectrons.size();
      }
      
      // total number of electrons after conversion check
      int howManyConvEles = -999;
      if (!_commonSel->getSwitch("asymmetricElectrons")) howManyConvEles = howManyConvTightEles; 
      if ( _commonSel->getSwitch("asymmetricElectrons")) howManyConvEles = howManyConvTightEles + howManyConvLooseEles; 

      int howManyDzVertexEles      = 0; 
      int howManyDzVertexLooseEles = 0; 
      int howManyDzVertexTightEles = 0; 
 

      // here the final electrons are selected for the analysis -------------------------

      // the highest pt identified, isolated and not from conversions electron, used to define the EWK vertex      
      std::pair<int,int> bestElectronConvTightPair = getBestGoodElePair(convTightElectrons);
      int theConvEle1(bestElectronConvTightPair.first);
      int theConvEle1_GSFTrack = gsfTrackIndexPFEle[theConvEle1];
      float EWKz = trackVzGsfTrack[theConvEle1_GSFTrack];
      closestPV = -1;   
      float z0 = 0.0;
      if(nPV>0) {         
        float minDzPV = 999.;
        for(int v=0; v<nPV; v++) {
          if(fabs(PVzPV[v]-EWKz)<minDzPV) {
            minDzPV=fabs(PVzPV[v]-EWKz);
            closestPV = v;
          }
        }
        z0 = PVzPV[closestPV];    // if nPV=0 --> z0 = 0
	

        // loop on tight and loose list
        for(int eleList=0; eleList<2; eleList++) {
	  int nConvEles = (eleList==0) ? howManyConvTightEles : howManyConvLooseEles;
	  
          // electron consistency with primary vertex
          for(int iiconv=0; iiconv<nConvEles; iiconv++) {
            bool isGoodEle = true;
            int iele = (eleList==0) ? convTightElectrons[iiconv] : convLooseElectrons[iiconv];

            // zele - zPV
            int gsfTrack = gsfTrackIndexPFEle[iele];
            float dz = trackVzGsfTrack[gsfTrack]-PVzPV[closestPV];
            if(_commonSel->getSwitch("dzVertex") && !_commonSel->passCut("dzVertex",dz)) isGoodEle = false;
            else {
              if(eleList==0) howManyDzVertexTightEles++;
              if(eleList==1) howManyDzVertexLooseEles++;
            }

            // dxy
            float dxy = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxGsfTrack[gsfTrack], trackVyGsfTrack[gsfTrack], trackVzGsfTrack[gsfTrack], pxPFEle[iele], pyPFEle[iele], pzPFEle[iele]);
            if(_commonSel->getSwitch("dxyVertex") && !_commonSel->passCut("dxyVertex",dxy)) isGoodEle = false;
	    
            if(isGoodEle) {
              if(eleList==0) vertexTightElectrons.push_back(iele); 
              else if(eleList==1) vertexLooseElectrons.push_back(iele);
            }
          }
        }
      }

      // total number of 'good' electrons
      int howManyTightVertexEles = vertexTightElectrons.size();
      int howManyLooseVertexEles = vertexLooseElectrons.size();
      int howManyEles      = -999;
      int howManyTightEles = -999;
      if (!_commonSel->getSwitch("asymmetricElectrons")) { 
        howManyEles         = howManyTightVertexEles; 
        howManyTightEles    = howManyTightVertexEles; 
        howManyDzVertexEles = howManyDzVertexTightEles; 
      }
      if (_commonSel->getSwitch("asymmetricElectrons"))  { 
        howManyTightEles    = howManyTightVertexEles; 
        howManyEles         = howManyTightVertexEles+howManyLooseVertexEles; 
        howManyDzVertexEles = howManyDzVertexTightEles + howManyDzVertexLooseEles; 
      }   

      // two highest pt good electrons - no charge requirement
      std::pair<int,int> bestElectronTightPair = getBestGoodElePair(vertexTightElectrons);
      int theEle1(bestElectronTightPair.first);
      int theEle2(bestElectronTightPair.second);
      if (_commonSel->getSwitch("asymmetricElectrons") && theEle2<0 ){
        std::pair<int,int> bestElectronLoosePair = getBestGoodElePair(vertexLooseElectrons);
        int myEle2(bestElectronLoosePair.first);
        theEle2 = myEle2;
      }

      // set the event kinematics
      setKinematics(theEle1, theEle2);

      // take from the event the PF reconstructed jets 
      _theAK5PFJets = GetUncorrPFJets();
      
      howManyWJets = 0;
      howManyWPFJets = 0; 
      if(theEle1_ > -1) {
        std::vector<TVector3> WElectron;
        WElectron.push_back(p4Ele1_.Vect());
        
   	// count calo jets
        JetCounter goodJets_counter(_theAK5CaloJets);
        goodJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
        goodJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodJets_counter.SetParticlesToRemove(WElectron);
        goodWJets_ = goodJets_counter.getGoodJets();
        howManyWJets = goodWJets_.size();

        // get btagging low threshold calo jets 
        JetCounter goodLowThrJets_counter(_theAK5CaloJets);
        goodLowThrJets_counter.SetThresholds(_commonSel->getLowerCut("etJetBvetoAcc"), _commonSel->getUpperCut("etaJetAcc")); 
        goodLowThrJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodLowThrJets_counter.SetParticlesToRemove(WElectron);
        m_cleanedBvetoWBTagJets = goodLowThrJets_counter.getGoodBTagJets();
	
        // count PF jets
        JetCounter goodPFJets_counter(_theAK5PFJets);
        goodPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
        goodPFJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
        goodPFJets_counter.SetParticlesToRemove(WElectron);
        goodWPFJets_ = goodPFJets_counter.getGoodJets();
        howManyWPFJets = goodWPFJets_.size();

        // get btagging low threshold calo jets 
        JetCounter goodLowThrPFJets_counter(_theAK5PFJets);
        goodLowThrPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetBvetoAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
        goodLowThrJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
        goodLowThrJets_counter.SetParticlesToRemove(WElectron);
        m_cleanedBvetoWBTagPFJets = goodLowThrJets_counter.getGoodBTagJets();
      } 
      
      howManyZJets = 0;
      howManyZPFJets = 0;
      if(theEle1_ > -1 && theEle2_ > -1) {
        std::vector<TVector3> ZElectrons;
        ZElectrons.push_back(p4Ele1_.Vect());
        ZElectrons.push_back(p4Ele2_.Vect());

        // count calo jets
        JetCounter goodJets_counter(_theAK5CaloJets);
        goodJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc"));
        goodJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodJets_counter.SetParticlesToRemove(ZElectrons);
        goodZJets_ = goodJets_counter.getGoodJets();
        howManyZJets = goodZJets_.size();

        // get btagging low threshold calo jets 
        JetCounter goodLowThrJets_counter(_theAK5CaloJets);
        goodLowThrJets_counter.SetThresholds(_commonSel->getLowerCut("etJetBvetoAcc"), _commonSel->getUpperCut("etaJetAcc")); 
        goodLowThrJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodLowThrJets_counter.SetParticlesToRemove(ZElectrons);
        m_cleanedBvetoZBTagJets = goodLowThrJets_counter.getGoodBTagJets();

        // count PF jets
        JetCounter goodPFJets_counter(_theAK5PFJets);
        goodPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
        goodPFJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
        goodPFJets_counter.SetParticlesToRemove(ZElectrons);
        goodZPFJets_ = goodPFJets_counter.getGoodJets();
        howManyZPFJets = goodZPFJets_.size();

        // get btagging low threshold calo jets 
        JetCounter goodLowThrPFJets_counter(_theAK5PFJets);
        goodLowThrPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetBvetoAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
        goodLowThrJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
        goodLowThrJets_counter.SetParticlesToRemove(ZElectrons);
        m_cleanedBvetoZBTagPFJets = goodLowThrJets_counter.getGoodBTagJets();
      }

      //! B vetoes
      // remove muons from the second top->W(mu nu)b decay
      int nGoodMuons = CountMuons(2.4, 15.);

      //! proper B tagging
      BTagJet jetBTagEVTW, jetBTagEVTZ;
      if(theEle1_ > -1) {
        calcEventBVetoVariables(m_cleanedBvetoWJets,m_cleanedBvetoWBTagJets);
        jetBTagEVTW = m_maxBTagEvt;
        calcEventBVetoVariables(m_cleanedBvetoWPFJets,m_cleanedBvetoWBTagPFJets);
        jetBTagEVTZ = m_maxBTagEvt;
      }

      BTagJet pfjetBTagEVTW, pfjetBTagEVTZ;
      if(theEle1_ > -1 && theEle2_ > -1) {
        calcEventBVetoVariables(m_cleanedBvetoZJets,m_cleanedBvetoZBTagJets);
        pfjetBTagEVTW = m_maxBTagEvt;
        calcEventBVetoVariables(m_cleanedBvetoZPFJets,m_cleanedBvetoZBTagPFJets);
        pfjetBTagEVTZ = m_maxBTagEvt;
      }


      // --------------------------------------------------------------
      bool isCommonSel, isZeeUpToEnd, isWenuUpToEnd, isZeeUpToEndPFJet, isWenuUpToEndPFJet ;
      isCommonSel = false;
      isZeeUpToEnd = false;
      isWenuUpToEnd = false;
      isZeeUpToEndPFJet = false;
      isWenuUpToEndPFJet = false;

      // for the Z veto for W->enu selection
      bool ZcandFound = foundZCandidate(theEle1_,acceptElectrons);
      
      // common data for the selection
      SelectorData commonData;
      commonData.weight    = weight;
      commonData.passedHLT = passedHLT;
      commonData.nRecoEle  = howManyRecoEle;
      commonData.nAccEle   = howManyAccEle;
      commonData.nIdTightEle = howManyIdTightEle;
      commonData.nIdLooseEle = howManyIdLooseEle;
      commonData.nIsolTightEle = howManyIsolTightEles;
      commonData.nIsolLooseEle = howManyIsolLooseEles;
      commonData.nConvRejTightEle = howManyConvTightEles;
      commonData.nConvRejLooseEle = howManyConvLooseEles;
      commonData.nPV = nPV;
      commonData.nDzVertexEle  = howManyDzVertexEles;
      commonData.nDxyVertexEle = howManyEles;
      commonData.nMuons = nGoodMuons;

      // Z + calo jets selection
      SelectorData dataZeeJet(commonData);
      if(m_signal == zjets) dataZeeJet.signal = ZToEEDecay;
      else if(m_signal == zother) dataZeeJet.signal = !ZToEEDecay;
      else dataZeeJet.signal = true;
      dataZeeJet.njetsMC = nzjets_mc;
      dataZeeJet.njets   = howManyZJets;
      dataZeeJet.btagEVT = jetBTagEVTZ.combinedSecondaryVertexBJetTags;
      dataZeeJet.mInv    = mee_;
      // dataZeeJet.mhtMET  = mah_;

      CutBasedSelectorEE selectorZeeJet(dataZeeJet);
      selectorZeeJet.isMc(!isData_);
      isZeeUpToEnd  = selectorZeeJet.outputZ(_commonSel,_zeeSel,_zeePCounter,&m_zeeJetbinCounter,&m_zeeRecoJetbinCounter);

      // W + calo jets selection
      SelectorData dataWenuJet(commonData);
      if(m_signal == wjets) dataWenuJet.signal = WToENuDecay;
      else if(m_signal == wother) dataWenuJet.signal = !WToENuDecay; 
      else dataWenuJet.signal = true;
      dataWenuJet.njetsMC   = nwjets_mc;
      dataWenuJet.njets     = howManyWJets;
      dataWenuJet.btagEVT   = jetBTagEVTW.combinedSecondaryVertexBJetTags;
      dataWenuJet.foundAnyZ = ZcandFound;
      dataWenuJet.met       = p3Met_.Mag();
      dataWenuJet.mt        = WmT_;
      // dataWenuJet.mhtJet    = mah_;

      CutBasedSelectorEE selectorWenuJet(dataWenuJet);
      selectorWenuJet.isMc(!isData_);
      isWenuUpToEnd  = selectorWenuJet.outputW(_commonSel,_wenuSel,_wenuPCounter,&m_wenuJetbinCounter,&m_wenuRecoJetbinCounter);

      // Z + PF jets selection
      SelectorData dataZeePFJet(commonData);
      if(m_signal == zjets) dataZeePFJet.signal = ZToEEDecay;
      else if(m_signal == zother) dataZeePFJet.signal = !ZToEEDecay;
      else dataZeePFJet.signal = true;
      dataZeePFJet.njetsMC     = nzpfjets_mc;
      dataZeePFJet.njets       = howManyZPFJets;
      dataZeePFJet.btagEVT     = pfjetBTagEVTZ.combinedSecondaryVertexBJetTags;
      dataZeePFJet.mInv        = mee_;
      // dataZeePFJet.mhtMET      = -999.;

      CutBasedSelectorEE selectorZeePFJet(dataZeePFJet);
      selectorZeePFJet.isMc(!isData_);
      isZeeUpToEndPFJet  = selectorZeePFJet.outputZ(_commonSel,_zeeSel,0,&m_zeePFJetbinCounter,&m_zeeRecoPFJetbinCounter);

      // W + PF jets selection
      SelectorData dataWenuPFJet(commonData);
      if(m_signal == wjets) dataWenuPFJet.signal = WToENuDecay;
      else if(m_signal == wother) dataWenuPFJet.signal = !WToENuDecay; 
      else dataWenuPFJet.signal = true;
      dataWenuPFJet.njetsMC     = nwpfjets_mc;
      dataWenuPFJet.njets       = howManyWPFJets;
      dataWenuPFJet.btagEVT     = jetBTagEVTW.combinedSecondaryVertexBJetTags;
      dataWenuPFJet.foundAnyZ   = ZcandFound;
      dataWenuPFJet.met         = p3Met_.Mag();
      dataWenuPFJet.mt          = WmT_;
      // dataWenuPFJet.mhtJet      = -999.;

      CutBasedSelectorEE selectorWenuPFJet(dataWenuPFJet);
      selectorWenuPFJet.isMc(!isData_);
      isWenuUpToEndPFJet  = selectorWenuPFJet.outputW(_commonSel,_wenuSel,0,&m_wenuPFJetbinCounter,&m_wenuRecoPFJetbinCounter);

      // --------------------------------------------------------------
      if (isWenuUpToEnd || isWenuUpToEndPFJet) {       // filling reduced trees for the W->enu study
	
	
	/*
	// MC level W, chiara
	int indexGenNu=999, indexGenEle=999;
	for(int imc=0;imc<100;imc++) {
	  if( abs(idMc[imc])==11 && fabs(idMc[mothMc[imc]])==24 ) indexGenEle = imc;
	  if( abs(idMc[imc])==12 && fabs(idMc[mothMc[imc]])==24 ) indexGenNu  = imc;
	}
	if (indexGenNu < 0 || indexGenEle < 0) cout << "chiara, ma???" << endl;
	cout << "indexGenNu = " << indexGenNu << ", indexGenEle = " << indexGenEle << endl;
	TVector3 trueEleP, trueNuP;
	trueEleP.SetPtThetaPhi(pMc[indexGenEle]*fabs(sin(thetaMc[indexGenEle])),thetaMc[indexGenEle],phiMc[indexGenEle]);
	trueNuP.SetPtThetaPhi(pMc[indexGenNu]*fabs(sin(thetaMc[indexGenNu])),thetaMc[indexGenNu],phiMc[indexGenNu]);
	float wGenPt = (trueEleP+trueNuP).Pt();
	cout << wGenPt << endl;
	myOutTree_Wenu -> fillArcs(wGenPt);
	// chiara
	*/	  

        myOutTree_Wenu -> fillJetMultiplicities( howManyWJets, howManyWPFJets, nwgenjets_mc, isWenuUpToEnd, isWenuUpToEndPFJet);
        myOutTree_Wenu -> fillElectronsPF(pt_, eta_, mva_, chargedIso_, neutralIso_, photonIso_, combinedIso_,charge_);
        myOutTree_Wenu -> fillKinematicsPF( WCalomT_, WTCmT_, WPFmT_, GetPt(pxMet[0],pyMet[0]), GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]));
        if(!isData_) myOutTree_Wenu->fillMcTruth(WToENuDecay);
        myOutTree_Wenu -> fillRunInfo(runNumber, lumiBlock, eventNumber);
        myOutTree_Wenu -> store();
      }

      if (isZeeUpToEnd || isZeeUpToEndPFJet ) {        // filling reduced trees for the Z->ee study

        vector<float> dummyW;
        for(int i=0;i<6;i++) dummyW.push_back(-999.);

        myOutTree_Zee -> fillJetMultiplicities( howManyZJets, howManyZPFJets, nzgenjets_mc, isZeeUpToEnd, isZeeUpToEndPFJet);
        myOutTree_Zee -> fillElectronsPF(pt_, eta_,mva_,chargedIso_, neutralIso_, photonIso_, combinedIso_,charge_);
        myOutTree_Zee -> fillKinematicsPF( -1., -1., -1., -1., -1., -1. );
        if(!isData_) myOutTree_Zee->fillMcTruth(ZToEEDecay);
        myOutTree_Zee -> fillRunInfo(runNumber, lumiBlock, eventNumber);
        myOutTree_Zee -> store();
      }
    }
    
  } // loope over entries
}  


void VecbosPFEESelection::displayEfficiencies() {

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "+++ DETAILED EFFICIENCY FOR SIGNAL AS A FUNCTION OF JET MULTIPLICITY +++" << std::endl;
  std::cout << "+++ Z(ee)+jets +++" << std::endl;
  sprintf(namefile,"%sVecBosEfficiency.root",m_prefix);
  for(int ii=0; ii<6; ii++) {
    std::cout << "\tZ+" << ii << "jet" << std::endl;
    displayZeeEfficiencies(m_zeeJetbinCounter[ii]);
    const char* option = (ii==0) ? "recreate" : "update";
    m_zeeJetbinCounter[ii].Save(namefile,option);
    m_zeePFJetbinCounter[ii].Save(namefile,"update");
  }
  std::cout << std::endl << std::endl << std::endl << std::endl;
  std::cout << "+++ W(enu)+jets +++" << std::endl;
  for(int ii=0; ii<6; ii++) {
    std::cout << "\tW+" << ii << "jet" << std::endl;
    displayWenuEfficiencies(m_wenuJetbinCounter[ii]);
    m_wenuJetbinCounter[ii].Save(namefile,"update");
    m_wenuPFJetbinCounter[ii].Save(namefile,"update");
  }
  std::cout << "----------------------------------------" << std::endl;

  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Full Zee selections (inclusive): -->  " << std::endl;
  displayZeeEfficiencies(_zeeCounter);
  std::cout << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Full Wenu selections (inclusive): --> " << std::endl;
  displayWenuEfficiencies(_wenuCounter);
  std::cout << std::endl;
  
  // save into ROOT files the events as a function of jet multiplicity
  // they are different by ones stored in VecBosEfficiency.root because
  // the jet counting is done removing reco electrons, not with MC truth electrons
  sprintf(namefile,"%sVecBosCounters.root",m_prefix);
  for(int ii=0; ii<6; ii++) {
    const char* option = (ii==0) ? "recreate" : "update";
    m_zeeRecoJetbinCounter[ii].Save(namefile,option);
    m_zeeRecoPFJetbinCounter[ii].Save(namefile,"update");
  }
  for(int ii=0; ii<6; ii++) {
    m_wenuRecoJetbinCounter[ii].Save(namefile,"update");
    m_wenuRecoPFJetbinCounter[ii].Save(namefile,"update");
  }
}

// two highest pT 'good' electrons
std::pair<int,int> VecbosPFEESelection::getBestGoodElePair(std::vector<int> goodElectrons) {

  int theEle1 =-1;
  int theEle2 =-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;
  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    TVector3 pEle(pxPFEle[eleIndex],pyPFEle[eleIndex],pzPFEle[eleIndex]);
    float thisPt=pEle.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ 
      maxPt2  = maxPt1; 
      maxPt1  = thisPt; 
      theEle2 = theEle1; 
      theEle1 = eleIndex; 
    }
    if (thisPt<maxPt1 && thisPt>maxPt2){ 
      maxPt2  = thisPt; 
      theEle2 = eleIndex; 
    }
  }
  return make_pair(theEle1,theEle2);
}

void VecbosPFEESelection::setKinematics(int theEle1, int theEle2) {
  theEle1_ = theEle1;
  theEle2_ = theEle2;
  if (theEle1_ > -1) {
    TVector3 p1;
    p1.SetXYZ(pxPFEle[theEle1_],pyPFEle[theEle1_],pzPFEle[theEle1_]);
    p4Ele1_.SetVectM(p1,0.000511);
    pT3Ele1_.SetXYZ(pxPFEle[theEle1_],pyPFEle[theEle1_],0.0);
  } else {
    p4Ele1_.SetXYZT(0.,0.,0.,0.);
  }
  if (theEle2 > -1) {
    TVector3 p2;
    p2.SetXYZ(pxPFEle[theEle2],pyPFEle[theEle2],pzPFEle[theEle2]);
    p4Ele2_.SetVectM(p2,0.000511);
    pT3Ele2_.SetXYZ(pxPFEle[theEle2],pyPFEle[theEle2],0.0);
  } else {
    p4Ele2_.SetXYZT(0.,0.,0.,0.);
  }

  //! set best electron(s) variables
  setElectrons();

  // missing transverse energy (used for the cut, even if we dump also the other ones)
  if(_commonSel->getSwitch("useTCMet"))      p3Met_.SetXYZ(pxTCMet[0],pyTCMet[0],0.0);
  else if(_commonSel->getSwitch("usePFMet")) p3Met_.SetXYZ(pxPFMet[0],pyPFMet[0],0.0);
  else p3Met_.SetXYZ(pxMet[0],pyMet[0],0.0);
   p3CaloMet_.SetXYZ(pxMet[0],pyMet[0],0.0);
   p3TCMet_.SetXYZ(pxTCMet[0],pyTCMet[0],0.0);
   p3PFMet_.SetXYZ(pxPFMet[0],pyPFMet[0],0.0);

  // ----------------------------------------------------------------------------------
  // W SPECIFIC VARIABLES
  // ----------------------------------------------------------------------------------
  if(theEle1_ > -1) {

    // direct W Pt = Pt(ele)+MET
    pt3W_ = p4Ele1_.Perp() + p3Met_;    
    
    // W transverse mass (used in the selection)
    WmT_ = sqrt(2 * p4Ele1_.Pt() * p3Met_.Mag() * (1-cos(pT3Ele1_.Angle(p3Met_))) );
    
    // to dump in final tree
    WCalomT_ = sqrt(2 * p4Ele1_.Pt() * p3CaloMet_.Mag() * (1-cos(pT3Ele1_.Angle(p3CaloMet_))) );
    WTCmT_   = sqrt(2 * p4Ele1_.Pt() * p3TCMet_.Mag()   * (1-cos(pT3Ele1_.Angle(p3TCMet_))) );   
    WPFmT_   = sqrt(2 * p4Ele1_.Pt() * p3PFMet_.Mag()   * (1-cos(pT3Ele1_.Angle(p3PFMet_))) );   
  } else {
    WmT_ = -1;
  }
  
  // -------------------------------------------------------------------------------
  // Z SPECIFIC VARIABLES
  // -------------------------------------------------------------------------------
  if(theEle1_ > -1 && theEle2_ > -1) {
    mee_ = getMee(theEle1_, theEle2_);        // invariant mass
  } else {
    mee_ = -1.;
  }
}

void VecbosPFEESelection::setElectrons() {
  
  int bestElectrons[2];
  bestElectrons[0] = theEle1_;
  bestElectrons[1] = theEle2_;
  TLorentzVector p4[2];
  p4[0] = p4Ele1_;
  p4[1] = p4Ele2_;  

  // i=0: hard electron, i=1: slow electron
  for(int i=0; i<2; i++) {
    int ele = bestElectrons[i];
    if(ele > -1) {
      pt_[i]          = p4[i].Pt();
      eta_[i]         = p4[i].Eta();
      charge_[i]      = chargePFEle[ele];
      mva_[i]         = MvaOutputPFEle[ele];
      // chiara: usato fino ad ora
      chargedIso_[i]  = chIso04vetoPFEle[ele];
      neutralIso_[i]  = nhIso04vetoPFEle[ele];
      photonIso_[i]   = phIso04vetoPFEle[ele];
      // chiara: per confronto con florian
      // chargedIso_[i]  = chIso04noVetoNVCPFEle[ele];
      // neutralIso_[i]  = nhIso04noVetoPFEle[ele];
      // photonIso_[i]   = phIso04noVetoPFEle[ele];
      combinedIso_[i] = combinedIsolation[ele];
    } else { // the electron is not selected
      pt_[i]          = -1.;
      eta_[i]         = -1.;
      charge_[i]      = -1;
      mva_[i]         = -1.;
      chargedIso_[i]  = -1.;
      neutralIso_[i]  = -1.;
      photonIso_[i]   = -1.;
      combinedIso_[i] = -1.;
    }
  }
  
}

// invariant mass
float VecbosPFEESelection::getMee(int theEle1, int theEle2) {
  TLorentzVector *_p1 = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *_p2 = new TLorentzVector(0.,0.,0.,0.);
  _p1 ->SetXYZT(pxPFEle[theEle1],pyPFEle[theEle1],pzPFEle[theEle1],energyPFEle[theEle1]);
  _p2 ->SetXYZT(pxPFEle[theEle2],pyPFEle[theEle2],pzPFEle[theEle2],energyPFEle[theEle2]);
  double theInvMass = (*_p1+*_p2).M();       
  return theInvMass;
}

// electron identification
bool VecbosPFEESelection::isEleID(int eleIndex, bool wantTight) {

  bool isBarrel = true;
  if (fabs(etaPFEle[eleIndex])>1.476) isBarrel = false;

  bool isGoodEle = true;
  if( isBarrel ) {
    if(  wantTight && !_commonSel->passCut("tightMvaEB",MvaOutputPFEle[eleIndex]) ) isGoodEle = false;
    if( !wantTight && !_commonSel->passCut("looseMvaEB",MvaOutputPFEle[eleIndex]) ) isGoodEle = false;
  }
  if( !isBarrel ) {
    if(  wantTight && !_commonSel->passCut("tightMvaEE",MvaOutputPFEle[eleIndex]) ) isGoodEle = false;
    if( !wantTight && !_commonSel->passCut("looseMvaEE",MvaOutputPFEle[eleIndex]) ) isGoodEle = false;
  }

  return isGoodEle;
}

void VecbosPFEESelection::FillEleIdIsolPFTree(int eleIndex, int eleid, int isol, int conv) {

  TVector3 p3Ele(pxPFEle[eleIndex],pyPFEle[eleIndex],pzPFEle[eleIndex]);
  float elePT         = p3Ele.Pt(); 
  float eleETA        = etaPFEle[eleIndex];
  float eleMVA        = MvaOutputPFEle[eleIndex];

  int matchedTrack    = trackIndexPFEle[eleIndex];
  int expInnerLayers  = expInnerLayersTrack[matchedTrack];
  float transvImpPar  = transvImpactParTrack[matchedTrack]; 

  int gsfPF          = gsfTrackIndexPFEle[eleIndex];  
  // float transvImpPar = transvImpactParGsfTrack[gsfPF]; 

  float eleChargedIso03nV = chIso03noVetoPFEle[eleIndex];
  float eleChargedIso04nV = chIso04noVetoPFEle[eleIndex];
  float eleChargedIso05nV = chIso05noVetoPFEle[eleIndex];
  float eleChargedIso03v  = chIso03vetoPFEle[eleIndex];
  float eleChargedIso04v  = chIso04vetoPFEle[eleIndex];
  float eleChargedIso05v  = chIso05vetoPFEle[eleIndex];

  float eleChargedIso03nvcnv = chIso03noVetoNVCPFEle[eleIndex];
  float eleChargedIso04nvcnv = chIso04noVetoNVCPFEle[eleIndex];
  float eleChargedIso05nvcnv = chIso05noVetoNVCPFEle[eleIndex];
  float eleChargedIso03nvcv  = chIso03vetoNVCPFEle[eleIndex];
  float eleChargedIso04nvcv  = chIso04vetoNVCPFEle[eleIndex];
  float eleChargedIso05nvcv  = chIso05vetoNVCPFEle[eleIndex];

  float eleNeutralIso03nV = nhIso03noVetoPFEle[eleIndex];
  float eleNeutralIso04nV = nhIso04noVetoPFEle[eleIndex];
  float eleNeutralIso05nV = nhIso05noVetoPFEle[eleIndex];
  float eleNeutralIso03v  = nhIso03vetoPFEle[eleIndex];
  float eleNeutralIso04v  = nhIso04vetoPFEle[eleIndex];
  float eleNeutralIso05v  = nhIso05vetoPFEle[eleIndex];

  float elePhotonsIso03nV = phIso03noVetoPFEle[eleIndex];
  float elePhotonsIso04nV = phIso04noVetoPFEle[eleIndex];
  float elePhotonsIso05nV = phIso05noVetoPFEle[eleIndex];
  float elePhotonsIso03v  = phIso03vetoPFEle[eleIndex];
  float elePhotonsIso04v  = phIso04vetoPFEle[eleIndex];
  float elePhotonsIso05v  = phIso05vetoPFEle[eleIndex];

  float eleCombIso        = combinedIsolation[eleIndex];

  
  // search eleID variables from matched egamma ele
  float matchedHoE    = -100.; 
  float matchedDeta   = -100.;
  float matchedDphi   = -100.;
  float matchedEoP    = -100.;
  float matchedEoPout = -100.;
  for(int egEle=0; egEle<nEle; egEle++) {
    if (gsfTrackIndexEle[egEle]==gsfPF) { 
      matchedHoE    = hOverEEle[egEle];
      matchedDeta   = deltaEtaAtVtxEle[egEle];    
      matchedDphi   = deltaPhiAtVtxEle[egEle];    
      matchedEoP    = eSuperClusterOverPEle[egEle];   
      matchedEoPout = eSeedOverPoutEle[egEle];    
    }
  }
                                                                        
  // myOutEleIdIsolPFTree_Wenu->fillAll(elePT, eleETA, eleMVA, matchedHoE, matchedDeta, matchedDphi, matchedEoP, matchedEoPout, eleChargedIso04v, eleNeutralIso04v, elePhotonsIso04v, eleCombIso, expInnerLayers, transvImpPar);

  myOutEleIdIsolPFTree_Wenu->fillIsol(
				      eleChargedIso03nV,eleChargedIso04nV,eleChargedIso05nV,eleChargedIso03v,eleChargedIso04v,eleChargedIso05v,
				      eleChargedIso03nvcnv,eleChargedIso04nvcnv,eleChargedIso05nvcnv,eleChargedIso03nvcv,eleChargedIso04nvcv,eleChargedIso05nvcv,
				      eleNeutralIso03nV,eleNeutralIso04nV,eleNeutralIso05nV,eleNeutralIso03v,eleNeutralIso04v,eleNeutralIso05v,
				      elePhotonsIso03nV,elePhotonsIso04nV,elePhotonsIso05nV,elePhotonsIso03v,elePhotonsIso04v,elePhotonsIso05v
				      );

  myOutEleIdIsolPFTree_Wenu->fillSetpsAfter(eleid, isol, conv);

  myOutEleIdIsolPFTree_Wenu->store();
}

void VecbosPFEESelection::ConfigCommonSelections(Selection* _selection) {

  _selection->addSwitch("isData");
  _selection->addSwitch("goodRunLS");
  _selection->addSwitch("eleId");
  _selection->addSwitch("asymmetricElectrons");
  _selection->addSwitch("useTCMet");
  _selection->addSwitch("usePFMet");
  _selection->addSwitch("dumpEleID");
  _selection->addCut("eventRange");
  _selection->addCut("etaElectronAcc");
  _selection->addCut("ptElectronAcc");
  _selection->addCut("tightMvaEB");
  _selection->addCut("tightMvaEE");
  _selection->addCut("looseMvaEB");
  _selection->addCut("looseMvaEE");
  _selection->addCut("conversionRejection");
  _selection->addCut("combinedIsolationTightEB");
  _selection->addCut("combinedIsolationLooseEB");
  _selection->addCut("combinedIsolationTightEE");
  _selection->addCut("combinedIsolationLooseEE");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("etJetBvetoAcc");
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaPFJetAcc");
  _selection->addCut("etPFJetAcc");
  _selection->addCut("etPFJetBvetoAcc");
  _selection->addCut("pfJetConeWidth");
  _selection->addCut("nPrimaryVertices");
  _selection->addCut("dzVertexEles");
  _selection->addCut("dxyVertexEles");
  _selection->summary();
}

void VecbosPFEESelection::ConfigZeeSelections(Selection* _selection, Counters *_counter) {
  
  _selection->addSwitch("trigger");
  _selection->addCut("nRecoEles");
  _selection->addCut("nAccEles");
  _selection->addCut("nIdEles");
  _selection->addCut("nIsolEles");
  _selection->addCut("nConvRejEles");
  _selection->addCut("nDzVertexEles");
  _selection->addCut("nDxyVertexEles"); 
  _selection->addCut("nMuons");
  _selection->addCut("btagEVT");
  _selection->addCut("meeCut");
  _selection->addCut("MHTphiMET");
  _selection->summary();
  
  _counter->SetTitle("Z->EE EVENT COUNTER");
  _counter->AddVar("event");
  _counter->AddVar("mcTruth");
  _counter->AddVar("trigger");
  _counter->AddVar("nRecoEles");
  _counter->AddVar("nAccEles");
  _counter->AddVar("nIdEles");
  _counter->AddVar("nIsolEles");
  _counter->AddVar("nConvRejEles");
  _counter->AddVar("nPrimaryVertices");
  _counter->AddVar("nDzVertexEles");
  _counter->AddVar("nDxyVertexEles");
  _counter->AddVar("nMuons");
  _counter->AddVar("btagEVT");
  _counter->AddVar("meeCut");
  _counter->AddVar("MHTphiMET");
  _counter->AddVar("fullSelection");
}

void VecbosPFEESelection::ConfigWenuSelections(Selection* _selection, Counters *_counter) {

  _selection->addSwitch("trigger");
  
  _selection->addCut("nRecoEles");
  _selection->addCut("nAccEles");
  _selection->addCut("nIdEles");
  _selection->addCut("nIsolEles");
  _selection->addCut("nConvRejEles");
  _selection->addCut("nFinalEles");
  _selection->addSwitch("ZVeto");
  _selection->addCut("nDzVertexEles");
  _selection->addCut("nDxyVertexEles"); 
  _selection->addCut("nMuons");
  _selection->addCut("btagEVT");
  _selection->addCut("metCut");
  _selection->addCut("MHTphiJet");
  _selection->addCut("transvMassCut");
  _selection->summary();
  
  _counter->SetTitle("W->ENU EVENT COUNTER");
  _counter->AddVar("event");
  _counter->AddVar("mcTruth");
  _counter->AddVar("trigger");
  _counter->AddVar("nRecoEles");
  _counter->AddVar("nAccEles");
  _counter->AddVar("nIdEles");
  _counter->AddVar("nIsolEles");
  _counter->AddVar("nConvRejEles");
  _counter->AddVar("ZVeto");
  _counter->AddVar("nPrimaryVertices");
  _counter->AddVar("nDzVertexEles"); 
  _counter->AddVar("nDxyVertexEles");
  _counter->AddVar("nMuons");
  _counter->AddVar("btagEVT");
  _counter->AddVar("metCut");
  _counter->AddVar("transvMassCut");
  _counter->AddVar("MHTphiJet");
  _counter->AddVar("fullSelection");
}
 
void VecbosPFEESelection::displayZeeEfficiencies(Counters _counter){ 
  _counter.Draw();
  if(_commonSel->getSwitch("isData")) {
    _counter.Draw("trigger","event");
  } else {
    _counter.Draw("mcTruth","event");
    _counter.Draw("trigger","mcTruth");
  }
  _counter.Draw("nRecoEles","trigger");
  _counter.Draw("nAccEles","nRecoEles");
  _counter.Draw("nIdEles","nAccEles");
  _counter.Draw("nIsolEles","nIdEles");
  _counter.Draw("nConvRejEles","nIsolEles");
  _counter.Draw("nPrimaryVertices","nConvRejEles");
  _counter.Draw("nDzVertex","nPrimaryVertices");
  _counter.Draw("nDxyVertex","nDzVertex");
  _counter.Draw("nMuons","nDxyVertex");
  _counter.Draw("btagEVT","nMuons");
  _counter.Draw("meeCut","btagEVT");
  _counter.Draw("MHTphiMET","meeCut");
  _counter.Draw("fullSelection","event");
} 

void VecbosPFEESelection::displayWenuEfficiencies(Counters _counter){ 
  _counter.Draw();
  if(_commonSel->getSwitch("isData")) {
    _counter.Draw("trigger","event");
  } else {
    _counter.Draw("mcTruth","event");
    _counter.Draw("trigger","mcTruth");
  }
  _counter.Draw("nRecoEles","trigger");
  _counter.Draw("nAccEles","nRecoEles");
  _counter.Draw("nIdEles","nAccEles");
  _counter.Draw("nIsolEles","nIdEles");
  _counter.Draw("nConvRejEles","nIsolEles");
  _counter.Draw("ZVeto","nConvRejEles");
  _counter.Draw("nPrimaryVertices","ZVeto");
  _counter.Draw("nDzVertex","nPrimaryVertices");
  _counter.Draw("nDxyVertex","nDzVertex");
  _counter.Draw("nMuons","nDxyVertex");
  _counter.Draw("btagEVT","nMuons");
  _counter.Draw("metCut","btagEVT");
  _counter.Draw("transvMassCut","metCut");
  _counter.Draw("MHTphiJet","transvMassCut");
  _counter.Draw("fullSelection","event");
}

bool VecbosPFEESelection::foundZCandidate(int eleIndex, std::vector<int> accEles) {

  int howManyAccEle = accEles.size();
  for(int iAccEle=0; iAccEle<howManyAccEle; iAccEle++) {
    int iAccEleIndex = accEles[iAccEle];
    float mass = getMee(eleIndex,iAccEleIndex);
    if(_zeeSel->passCut("meeCut",mass)) return true;
  }
  
  return false;
}

int VecbosPFEESelection::CountMuons(float etaMax, float ptMin) {

  int nMuons = 0;
  for(int imu=0; imu<nMuon; imu++) {

    // acceptance
    if(fabs(etaMuon[imu])>etaMax) continue;
    if(GetPt(pxMuon[imu],pyMuon[imu])<ptMin) continue;
    
    nMuons++;

  }

  return nMuons;

}

bool VecbosPFEESelection::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 v(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = v.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(v.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

std::vector<float> VecbosPFEESelection::jetBTagVariables(Jet jet) {

  float nTracks   = 0;
  float sumNumDxy = 0.;
  float sumNumDsz = 0;
  float sumDen    = 0.;

  TVector3 p3Jet = jet.Get3Vector(); 
  TLorentzVector p4TracksInJet(0.,0.,0.,0.);

  for(int iTrack=0; iTrack<nTrack; iTrack++) {
    if (!isGoodTrack(iTrack,0.5,500,20,2.4,5)) continue;
    TVector3 p3Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
    TLorentzVector p4Track(p3Track,p3Track.Mag());      // assume mass=0...

    float deltaR = p3Jet.DeltaR(p3Track);
    if(fabs(deltaR)<0.5){

      nTracks++;

      float weight = p3Track.Pt() * p3Track.Pt() * p3Track.Pt() * p3Track.Pt();

      float dxy = fabs(eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxTrack[iTrack], trackVyTrack[iTrack], trackVzTrack[iTrack], pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack]));
      float w_dxy = dxy*weight;
      float dsz = fabs(eleDszPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxTrack[iTrack], trackVyTrack[iTrack], trackVzTrack[iTrack], pxTrack[iTrack], pyTrack[iTrack], pzTrack[iTrack]));
      float w_dsz = dsz*weight;

      sumNumDxy = sumNumDxy + w_dxy;
      sumNumDsz = sumNumDsz + w_dsz;
      sumDen    = sumDen + weight;
      
      p4TracksInJet += p4Track;
    }
  } 
  
  float DxyAverage = 0.0;
  float DszAverage = 0.0;

  if (sumDen != 0) {
    DxyAverage = sumNumDxy / sumDen;
    DszAverage = sumNumDsz / sumDen;
  }

  float jetMass = p4TracksInJet.M();

  std::vector<float> variables;
  variables.clear();
  variables.push_back(DxyAverage);
  variables.push_back(DszAverage);
  variables.push_back(jetMass);
  variables.push_back(nTracks);

  return variables;

}

void VecbosPFEESelection::calcEventBVetoVariables(std::vector<Jet> jets, std::vector<BTagJet> btags) {
  
  m_maxBTagEvt.combinedSecondaryVertexBJetTags = -1000.;
  m_maxBTagEvt.combinedSecondaryVertexMVABJetTags = -1000.;
  m_maxBTagEvt.jetBProbabilityBJetTags = -1000.;
  m_maxBTagEvt.jetProbabilityBJetTags = -1000.;
  m_maxBTagEvt.simpleSecondaryVertexBJetTags = -1000.;
  m_maxBTagEvt.softMuonBJetTags = -1000.;
  m_maxBTagEvt.trackCountingHighPurBJetTags = -1000.;
  m_maxBTagEvt.trackCountingHighEffBJetTags = -1000.;
  
  for (unsigned int iJet=0; iJet<jets.size(); iJet++){     
    Jet thisJet = jets[iJet];
    std::vector<float> variables = jetBTagVariables(thisJet);    
    float dxy = variables[0];
    float dsz = variables[1];
//     if (dxy > m_maxDxyEvt) m_maxDxyEvt=dxy;
//     if (dsz > m_maxDszEvt) m_maxDszEvt=dsz;
    
    // full b-tag output for each jet (only for calojets, so far)
    if ( btags.size() > 0 ) {
      if ( jets.size() != btags.size() ) cout << "VERY NASTY ERROR: jets.size() != btags.size()" << endl;
      BTagJet btag = btags[iJet];
      if (btag.combinedSecondaryVertexBJetTags > m_maxBTagEvt.combinedSecondaryVertexBJetTags) m_maxBTagEvt.combinedSecondaryVertexBJetTags=btag.combinedSecondaryVertexBJetTags;
      if (btag.combinedSecondaryVertexMVABJetTags > m_maxBTagEvt.combinedSecondaryVertexMVABJetTags) m_maxBTagEvt.combinedSecondaryVertexMVABJetTags=btag.combinedSecondaryVertexMVABJetTags;
      if (btag.jetBProbabilityBJetTags > m_maxBTagEvt.jetBProbabilityBJetTags) m_maxBTagEvt.jetBProbabilityBJetTags=btag.jetBProbabilityBJetTags;
      if (btag.jetProbabilityBJetTags > m_maxBTagEvt.jetProbabilityBJetTags) m_maxBTagEvt.jetProbabilityBJetTags=btag.jetProbabilityBJetTags;
      if (btag.simpleSecondaryVertexBJetTags > m_maxBTagEvt.simpleSecondaryVertexBJetTags) m_maxBTagEvt.simpleSecondaryVertexBJetTags=btag.simpleSecondaryVertexBJetTags;
      if (btag.softMuonBJetTags > m_maxBTagEvt.softMuonBJetTags) m_maxBTagEvt.softMuonBJetTags=btag.softMuonBJetTags;
      if (btag.trackCountingHighPurBJetTags > m_maxBTagEvt.trackCountingHighPurBJetTags) m_maxBTagEvt.trackCountingHighPurBJetTags=btag.trackCountingHighPurBJetTags;
      if (btag.trackCountingHighEffBJetTags > m_maxBTagEvt.trackCountingHighEffBJetTags) m_maxBTagEvt.trackCountingHighEffBJetTags=btag.trackCountingHighEffBJetTags;
    }
  }

}

