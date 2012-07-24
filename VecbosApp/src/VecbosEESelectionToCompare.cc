#include <string>
#include <iostream>
#include "TMath.h"

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/Skimmer.hh"
#include "EgammaAnalysisTools/include/ElectronTrackerIsolation.hh"
#include "EgammaAnalysisTools/include/ElectronCaloIsolation.hh"
#include "EgammaAnalysisTools/include/ElectronBestCandidateSelector.hh"
#include "include/CaloTower.hh"
#include "include/VecbosEESelectionToCompare.hh"
#include "include/JetCounter.hh"

   

using namespace bits;

VecbosEESelectionToCompare::VecbosEESelectionToCompare(TTree *tree) 
  : Vecbos(tree) {
  
  // default do not check on mc-truth
  m_signal = all;

  // single electron efficiency
  //  EgammaCutBasedID.Configure("config/vecbos/");
  cout << "=== CONFIGURING TIGHT ELECTRON ID ===" << endl;
  EgammaTightID_NC.ConfigureNoClass("config/vecbos/electronId/tightWP70");
  cout << "=== CONFIGURING LOOSE ELECTRON ID ===" << endl;
  EgammaLooseID_NC.ConfigureNoClass("config/vecbos/electronId/looseWP95");
  cout << "=== DONE ELECTRON ID CONFIGURATION." << endl;

  // common kinematic selections
  std::string theConfigDir       = "config/vecbos/";
  std::string fileCutsCommon     = theConfigDir + "CommonSelectionCutsElectronsDefaultSel.txt";
  std::string fileSwitchesCommon = theConfigDir + "CommonSelectionSwitchesElectronsDefaultSel.txt";

  _commonSel = new Selection(fileCutsCommon,fileSwitchesCommon);
  ConfigCommonSelections(_commonSel);


  std::cout << "[GoodRunLS]::goodRunLS is " << _commonSel->getSwitch("goodRunLS") << " isData is " <<  _commonSel->getSwitch("isData") << std::endl;

  //To read good run list!
  if (_commonSel->getSwitch("goodRunLS") && _commonSel->getSwitch("isData"))
    {
      std::string goodRunGiasoneFile       = "config/vecbos/json/goodRunLS.json";
      setJsonGoodRunList(goodRunGiasoneFile);
      fillRunLSMap();
    }

  // Z -> ee specific selections
  std::string fileCutsZee     = theConfigDir + "ZeeCutsDefaultSel.txt";
  std::string fileSwitchesZee = theConfigDir + "ZeeSwitchesDefaultSel.txt";
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
  std::string fileCutsWenu     = theConfigDir + "WenuCutsDefaultSel.txt";
  std::string fileSwitchesWenu = theConfigDir + "WenuSwitchesDefaultSel.txt";
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

  // configuration for the isolation Fisher discriminant
  // spectator variables used for the training (necessary to the reader, too, even if dummy)
  float dummy1, dummy2, dummy3;
  reader_EB = new TMVA::Reader( "!Color:Silent" );
  reader_EB->AddVariable( "trackerIsolRel := trackerIsol/pt", &trackerIsolRel );
  reader_EB->AddVariable( "ecalIsolRel := ecalJurIsol/pt", &ecalIsolRel );
  reader_EB->AddVariable( "hcalIsolRel := hcalIsol/pt", &hcalIsolRel );
  reader_EB->AddSpectator( "weight", &dummy1 );
  reader_EB->AddSpectator( "pt",     &dummy2 );
  reader_EB->AddSpectator( "eta",    &dummy3 );
  reader_EB->BookMVA( "Fisher", "config/vecbos/electronIsolation/EB_Fisher.weights.xml" );
  
  reader_EE = new TMVA::Reader( "!Color:Silent" );
  reader_EE->AddVariable( "trackerIsolRel := trackerIsol/pt", &trackerIsolRel );
  reader_EE->AddVariable( "ecalIsolRel := ecalJurIsol/pt", &ecalIsolRel );
  reader_EE->AddVariable( "hcalIsolRel := hcalIsol/pt", &hcalIsolRel );
  reader_EE->AddSpectator( "weight", &dummy1 );
  reader_EE->AddSpectator( "pt",     &dummy2 );
  reader_EE->AddSpectator( "eta",    &dummy3 );
  reader_EE->BookMVA( "Fisher", "config/vecbos/electronIsolation/EE_Fisher.weights.xml" );

  m_doBestElectronStudy = false;

  m_skimOutputFile = "";
  m_skimFile = "";
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

VecbosEESelectionToCompare::~VecbosEESelectionToCompare(){
  
  // closing/saving reduced trees
  myOutTree_Zee            -> save();
  myOutTree_Wenu           -> save();
//   myOutMcTree_Zee          -> save();
//   myOutMcTree_Wenu         -> save();
//   myOutVertexTree_FullWenu -> save();
//   myOutVertexTree_AccWenu  -> save();
//   myOutVertexTree_IdWenu   -> save();
  if(_commonSel->getSwitch("dumpIsolation")) myOutIsolationVtxOptimTree_Wenu -> save();
  if(_commonSel->getSwitch("dumpEleID")) myOutEleIDOptimTree_Wenu -> save();

  // deleting selections
  delete _commonSel;
  delete _zeeSel;
  delete _wenuSel;
  
}

// loop over events - real analysis
void VecbosEESelectionToCompare::Loop() {

  _verbose=false;
  if(fChain == 0) return;

  // kinematics reduced tree
  sprintf(namefile,"%sVecBosOutput-out-Zee.root",m_prefix);
  myOutTree_Zee  = new RedVecbosTree(namefile);
  sprintf(namefile,"%sVecBosOutput-out-Wenu.root",m_prefix);
  myOutTree_Wenu = new RedVecbosTree(namefile);
  myOutTree_Zee->addMcTruthInfos();
  myOutTree_Wenu->addMcTruthInfos();
  myOutTree_Zee->addBTagEVTInfos();
  myOutTree_Wenu->addBTagEVTInfos();
  myOutTree_Zee->addElectronInfos();
  myOutTree_Wenu->addElectronInfos();

  // mc truth tree for signal only
//   sprintf(namefile,"%sVecBosOutput-out-McZee.root",m_prefix);
//   myOutMcTree_Zee  = new RedVecbosMcTree(namefile);
//   sprintf(namefile,"%sVecBosOutput-out-McWenu.root",m_prefix);
//   myOutMcTree_Wenu = new RedVecbosMcTree(namefile);

  // study the best electron choice
  if ( m_doBestElectronStudy ) {
    bookBestCandidateHistos();
  }

  // reduced tree to study vertex variables
  //  sprintf(namefile,"%sVecBosOutput-out-Vtx-FullWenu.root",m_prefix);
//   myOutVertexTree_FullWenu = new RedVecbosVertexTree(namefile);
//   sprintf(namefile,"%sVecBosOutput-out-Vtx-AccWenu.root",m_prefix);
//   myOutVertexTree_AccWenu  = new RedVecbosVertexTree(namefile);
//   sprintf(namefile,"%sVecBosOutput-out-Vtx-IdWenu.root",m_prefix);
//   myOutVertexTree_IdWenu   = new RedVecbosVertexTree(namefile);
//   myOutVertexTree_FullWenu -> addMcTruthInfos();
//   myOutVertexTree_AccWenu  -> addMcTruthInfos();
//   myOutVertexTree_IdWenu   -> addMcTruthInfos();

  // reduced tree for isolation - vertexing optimization
  if(_commonSel->getSwitch("dumpIsolation")) {
    sprintf(namefile,"%sVecBosOutput-IsolVtx.root",m_prefix);
    myOutIsolationVtxOptimTree_Wenu = new RedIsolationVtxOptimTree(namefile);
    myOutIsolationVtxOptimTree_Wenu->addKinematicsInfos();
  }

  // reduced tree for electron ID
  if(_commonSel->getSwitch("dumpEleID")) {
    sprintf(namefile,"%sVecBosOutput-EleID.root",m_prefix);
    myOutEleIDOptimTree_Wenu = new RedEleIDOptimTree(namefile);
  }

  // the txtfile to write for the skim
  ofstream skimFile;
  bool writeSkim = false;
  if(strcmp(m_skimOutputFile,"") != 0 ) {
    writeSkim = true;
    skimFile.open(m_skimOutputFile);
  }

  // use certain skim specified by file
  bool SKIM = false;
  Skimmer skimmer(m_skimFile);
  if( strcmp(m_skimFile,"") != 0 ) {
    SKIM = true;
    skimmer.readFile();
  }

  // matrix NjReco - NjTrue
  TH2F *NWjetMatrix = new TH2F("NWjetMatrix","NWjetMatrix",6,0,6,6,0,6);
  NWjetMatrix->GetXaxis()->SetTitle("N_{jet}^{true}");
  NWjetMatrix->GetYaxis()->SetTitle("N_{jet}^{reco}");
  TH2F *NWPFjetMatrix = new TH2F("NWPFjetMatrix","NWPFjetMatrix",6,0,6,6,0,6);
  NWPFjetMatrix->GetXaxis()->SetTitle("N_{pfjet}^{true}");
  NWPFjetMatrix->GetYaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NWjetMatrixIncl = new TH2F("NWjetMatrixIncl","NWjetMatrixIncl",6,0,6,6,0,6);
  NWjetMatrixIncl->GetXaxis()->SetTitle("N_{jet}^{true}");
  NWjetMatrixIncl->GetYaxis()->SetTitle("N_{jet}^{reco}");
  TH2F *NWPFjetMatrixIncl = new TH2F("NWPFjetMatrixIncl","NWPFjetMatrixIncl",6,0,6,6,0,6);
  NWPFjetMatrixIncl->GetXaxis()->SetTitle("N_{PFjet}^{true}");
  NWPFjetMatrixIncl->GetYaxis()->SetTitle("N_{PFjet}^{reco}");

  TH2F *NZjetMatrix = new TH2F("NZjetMatrix","NZjetMatrix",6,0,6,6,0,6);
  NZjetMatrix->GetXaxis()->SetTitle("N_{jet}^{true}");
  NZjetMatrix->GetYaxis()->SetTitle("N_{jet}^{reco}");
  TH2F *NZjetMatrixIncl = new TH2F("NZjetMatrixIncl","NZjetMatrixIncl",6,0,6,6,0,6);
  NZjetMatrixIncl->GetXaxis()->SetTitle("N_{jet}^{true}");
  NZjetMatrixIncl->GetYaxis()->SetTitle("N_{jet}^{reco}");
  TH2F *NZPFjetMatrix = new TH2F("NZPFjetMatrix","NZPFjetMatrix",6,0,6,6,0,6);
  NZPFjetMatrix->GetXaxis()->SetTitle("N_{PFjet}^{true}");
  NZPFjetMatrix->GetYaxis()->SetTitle("N_{PFjet}^{reco}");
  TH2F *NZPFjetMatrixIncl = new TH2F("NZPFjetMatrixIncl","NZPFjetMatrixIncl",6,0,6,6,0,6);
  NZPFjetMatrixIncl->GetXaxis()->SetTitle("N_{PFjet}^{true}");
  NZPFjetMatrixIncl->GetYaxis()->SetTitle("N_{PFjet}^{reco}");

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

    if ( SKIM ) {
      if ( ! skimmer.output( jentry ) ) {
        eventToAnalyze = false;
      }
    }

    if ( eventToAnalyze ) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      

      //Good Run selection
      if (isData_ && _commonSel->getSwitch("goodRunLS") && !isGoodRunLS())
	{
	  if ( lastRun!= runNumber || lastLumi != lumiBlock)
	    {
	      lastRun = runNumber;
	      lastLumi = lumiBlock;
	      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
	    }
	  continue;
	}

      if (isData_ && _commonSel->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) )
	{
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
      //bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
      bool passedHLT = true;

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
      std::vector<int> acceptElectrons, acceptElectronsTmp; 
      std::vector<int> idTightElectrons, idLooseElectrons; 
      std::vector<int> isolTightElectrons, isolLooseElectrons; 
      std::vector<int> vertexTightElectrons, vertexLooseElectrons;

      trackerIsolation.clear();
      ecalIsolation.clear();
      hcalIsolation.clear();
      combinedIsolation.clear();
      // reconstructed electrons
      int howManyRecoEle = nEle;

     // electrons in the acceptance - with loose cut
      for(int ii=0;ii<howManyRecoEle;ii++) {
        bool isGoodEle = true;
        TVector3 pLepton(pxEle[ii],pyEle[ii],pzEle[ii]);
        float thisPt=pLepton.Pt();
        if(!_commonSel->passCut("etaElectronAcc",etaEle[ii]) &&
           anaUtils.isInECALFiducial(fiducialFlagsEle[ii]) ) isGoodEle = false;
        if(!_commonSel->passCut("ptElectronAcc",thisPt) ) isGoodEle = false;
        if (isGoodEle) acceptElectronsTmp.push_back(ii);
      }
      std::pair<int,int> bestElectronAcceptTmp = getBestGoodElePair(acceptElectronsTmp);
      int theEleAccepTmp(bestElectronAcceptTmp.first);
      TVector3 pLeptonTmp(pxEle[theEleAccepTmp],pyEle[theEleAccepTmp],pzEle[theEleAccepTmp]);
      float thisPtTmp=pLeptonTmp.Pt();

      // electrons in the acceptance - with tight cut on the first ele
      if (thisPtTmp>20.) {
  
	// electrons in the acceptance
	for(int ii=0;ii<howManyRecoEle;ii++) {             
	  bool isGoodEle = true;
	  TVector3 pLepton(pxEle[ii],pyEle[ii],pzEle[ii]);
	  float thisPt=pLepton.Pt();
	  if(!_commonSel->passCut("etaElectronAcc",etaEle[ii]) && 
	     anaUtils.isInECALFiducial(fiducialFlagsEle[ii]) ) isGoodEle = false;      
	  if(!_commonSel->passCut("ptElectronAcc",thisPt) ) isGoodEle = false; 
	  if (isGoodEle) acceptElectrons.push_back(ii);

	  // assign to each electron the three isolations
	  trackerIsolation.push_back( dr03TkSumPtEle[ii]  );
	  ecalIsolation.push_back( dr04EcalRecHitSumEtEle[ii]  );
	  hcalIsolation.push_back( dr04HcalTowerSumEtEle[ii] );
	  // and the combined isolation
	  trackerIsolRel = dr03TkSumPtEle[ii] / pLepton.Pt();
	  ecalIsolRel = dr04EcalRecHitSumEtEle[ii] / pLepton.Pt();
	  hcalIsolRel = dr04HcalTowerSumEtEle[ii] / pLepton.Pt();
	  float fisher = -999.;
	  if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[ii],isEB)) {
	    fisher = reader_EB->EvaluateMVA( "Fisher" );
	  } else {
	    fisher = reader_EE->EvaluateMVA( "Fisher" );
	  }
	  combinedIsolation.push_back( fisher );
	}
      }
      int howManyAccEle = acceptElectrons.size();

      // for vertex studies - to check before (eleId + isolation) and after (reco + acceptance)
      std::pair<int,int> bestElectronAccept = getBestGoodElePair(acceptElectrons);
      int theEleAccep1(bestElectronAccept.first);

      // 'tight' (optimized) identified electrons 
      for(int ii=0;ii<howManyAccEle;ii++) {      
        bool isGoodEle = true;
        // if (_commonSel->getSwitch("eleId") && !isEleID(&EgammaCutBasedID,acceptElectrons[ii])) isGoodEle = false; // class dependent eleID
        if (_commonSel->getSwitch("eleId") && !isEleID(&EgammaTightID_NC,acceptElectrons[ii])) isGoodEle = false;
        if (isGoodEle) idTightElectrons.push_back(acceptElectrons[ii]);            
        if(_commonSel->getSwitch("dumpEleID")) FillEleIDOptimTree(acceptElectrons[ii]);
      }
      int howManyIdTightEle = idTightElectrons.size();

      // egamma loose identified electrons, if needed 
      int howManyIdLooseEle = 0;
      if (_commonSel->getSwitch("asymmetricElectrons")) {   
        for(int ii=0;ii<howManyAccEle;ii++) {      
          int accIndex=acceptElectrons[ii];
          // if (anaUtils.electronIdVal(eleIdCutsEle[accIndex],eleIdLoose)) {
          if ((_commonSel->getSwitch("eleId") && isEleID(&EgammaLooseID_NC,acceptElectrons[ii])) || 
              (!_commonSel->getSwitch("eleId")) ) {
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

      // for vertex studies - to check before isolation and after (reco+accept+eleID)
      std::pair<int,int> bestElectronIdent = getBestGoodElePair(idTightElectrons);
      int theEleIdent1(bestElectronIdent.first);

      int howManyTkIsolEle=0,      howManyEcalIsolEle=0;      // for intermediate isolation efficiencies
      int howManyTkIsolTightEle=0, howManyEcalIsolTightEle=0; // for intermediate isolation efficiencies
      int howManyTkIsolLooseEle=0, howManyEcalIsolLooseEle=0; // for intermediate isolation efficiencies
      // tight optimized, tight isolated electrons
      for(int ii=0;ii<howManyIdTightEle;ii++) {      
        bool isGoodEle = true;
        int idEleIndex=idTightElectrons[ii];
        // combined isolated electrons
        if(_commonSel->getSwitch("combinedIsolationTightEB") && _commonSel->getSwitch("combinedIsolationTightEE")) {
          if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[idEleIndex],isEB) && 
             !_commonSel->passCut("combinedIsolationTightEB",combinedIsolation[idEleIndex])) isGoodEle = false;
          if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[idEleIndex],isEE) &&
             !_commonSel->passCut("combinedIsolationTightEE",combinedIsolation[idEleIndex])) isGoodEle = false;
        } else {

	  bool isBarrel = false;
          bool isEndcap = false;
          if(fabs(etaEle[idEleIndex])<1.476) isBarrel = true;
          if(fabs(etaEle[idEleIndex])>1.476) isEndcap = true;

          // tracker isolated electrons
	  bool isGoodTkIsoEle = true;
          float relSumPtTracks03 = trackerIsolation[idEleIndex];
	  if ( isBarrel && _commonSel->getSwitch("trackerPtSumTightEB") && !_commonSel->passCut("trackerPtSumTightEB", relSumPtTracks03) ) {
	    isGoodEle = false;
	    isGoodTkIsoEle = false;
	  }
	  if ( isEndcap && _commonSel->getSwitch("trackerPtSumTightEE") && !_commonSel->passCut("trackerPtSumTightEE", relSumPtTracks03) ) {
	    isGoodEle = false;
	    isGoodTkIsoEle = false;
	  }
          if(isGoodTkIsoEle) howManyTkIsolTightEle++;
          
          // ecal isolated electrons 
	  bool isGoodECALIsoEle = true;
          float relSumEtECAL04 = ecalIsolation[idEleIndex];
	  if ( isBarrel && _commonSel->getSwitch("ecalPtSumTightEB") && !_commonSel->passCut("ecalPtSumTightEB", relSumEtECAL04) ) {
	    isGoodEle = false;
	    isGoodECALIsoEle = false;
	  }
	  if ( isEndcap && _commonSel->getSwitch("ecalPtSumTightEE") && !_commonSel->passCut("ecalPtSumTightEE", relSumEtECAL04) ) {
	    isGoodEle = false;
	    isGoodECALIsoEle = false;
	  }
          if(isGoodECALIsoEle) howManyEcalIsolTightEle++;

          // hcal isolated electrons 
          float relSumEtHCAL04 = hcalIsolation[idEleIndex];
	  if ( isBarrel && _commonSel->getSwitch("hcalPtSumTightEB") && !_commonSel->passCut("hcalPtSumTightEB", relSumEtHCAL04) ) {
	    isGoodEle = false;
	  }
	  if ( isEndcap && _commonSel->getSwitch("hcalPtSumTightEE") && !_commonSel->passCut("hcalPtSumTightEE", relSumEtHCAL04) ) {
	    isGoodEle = false;
	  }
        }
        if (isGoodEle) isolTightElectrons.push_back(idEleIndex); 
      }
      int howManyIsolTightEles = isolTightElectrons.size();


      // loose optimized, loose isolated electrons if needed
      int howManyIsolLooseEles = 0;
      if (_commonSel->getSwitch("asymmetricElectrons")) {   
        for(int ii=0;ii<howManyIdLooseEle;ii++) {      
          bool isGoodEle = true;
          int idEleIndex=idLooseElectrons[ii];

          // combined isolated electrons
          if(_commonSel->getSwitch("combinedIsolationLooseEB") && _commonSel->getSwitch("combinedIsolationLooseEE")) {
            if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[idEleIndex],isEB) && 
               !_commonSel->passCut("combinedIsolationLooseEB",combinedIsolation[idEleIndex])) isGoodEle = false;
            if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[idEleIndex],isEE) &&
               !_commonSel->passCut("combinedIsolationLooseEE",combinedIsolation[idEleIndex])) isGoodEle = false;
          } else {

	    bool isBarrel = false;
	    bool isEndcap = false;
	    if(fabs(etaEle[idEleIndex])<1.476) isBarrel = true;
	    if(fabs(etaEle[idEleIndex])>1.476) isEndcap = true;
	    
	    // tracker isolated electrons
	    bool isGoodTkIsoEle = true;
	    float relSumPtTracks03 = trackerIsolation[idEleIndex];
	    if ( isBarrel && _commonSel->getSwitch("trackerPtSumLooseEB") && !_commonSel->passCut("trackerPtSumLooseEB", relSumPtTracks03) ) {
	      isGoodEle = false;
	      isGoodTkIsoEle = false;
	    }
	    if ( isEndcap && _commonSel->getSwitch("trackerPtSumLooseEE") && !_commonSel->passCut("trackerPtSumLooseEE", relSumPtTracks03) ) {
	      isGoodEle = false;
	      isGoodTkIsoEle = false;
	    }
	    if(isGoodTkIsoEle) howManyTkIsolLooseEle++;
	    
	    // ecal isolated electrons 
	    bool isGoodECALIsoEle = true;
	    float relSumEtECAL04 = ecalIsolation[idEleIndex];
	    if ( isBarrel && _commonSel->getSwitch("ecalPtSumLooseEB") && !_commonSel->passCut("ecalPtSumLooseEB", relSumEtECAL04) ) {
	      isGoodEle = false;
	      isGoodECALIsoEle = false;
	    }
	    if ( isEndcap && _commonSel->getSwitch("ecalPtSumLooseEE") && !_commonSel->passCut("ecalPtSumLooseEE", relSumEtECAL04) ) {
	      isGoodEle = false;
	      isGoodECALIsoEle = false;
	    }
	    if(isGoodECALIsoEle) howManyEcalIsolLooseEle++;
	    
	    // hcal isolated electrons 
	    float relSumEtHCAL04 = hcalIsolation[idEleIndex];
	    if ( isBarrel && _commonSel->getSwitch("hcalPtSumLooseEB") && !_commonSel->passCut("hcalPtSumLooseEB", relSumEtHCAL04) ) {
	      isGoodEle = false;
	    }
	    if ( isEndcap && _commonSel->getSwitch("hcalPtSumLooseEE") && !_commonSel->passCut("hcalPtSumLooseEE", relSumEtHCAL04) ) {
	      isGoodEle = false;
	    }
          }
          if (isGoodEle) isolLooseElectrons.push_back(idEleIndex); 
        }
        howManyIsolLooseEles = isolLooseElectrons.size();
      }

      // total number of isolated electrons
      int howManyIsolEles = -999;
      if (!_commonSel->getSwitch("asymmetricElectrons")) { 
        howManyIsolEles    = howManyIsolTightEles; 
        howManyTkIsolEle   = howManyTkIsolTightEle;
        howManyEcalIsolEle = howManyEcalIsolTightEle;
      }   
      if (_commonSel->getSwitch("asymmetricElectrons"))  { 
        howManyIsolEles    = howManyIsolTightEles    + howManyIsolLooseEles; 
        howManyTkIsolEle   = howManyTkIsolTightEle   + howManyTkIsolLooseEle;
        howManyEcalIsolEle = howManyEcalIsolTightEle + howManyEcalIsolLooseEle;
      }   

      // the highest pt isolated electron, used to define the EWK vertex
      std::pair<int,int> bestElectronIsolTightPair = getBestGoodElePair(isolTightElectrons);
      int theIsolEle1(bestElectronIsolTightPair.first);
      int theIsolEle1_GSFTrack = gsfTrackIndexEle[theIsolEle1];

      int howManyDzVertexEles = 0; 
      int howManyDzVertexLooseEles = 0; 
      int howManyDzVertexTightEles = 0; 
      float EWKz = trackVzGsfTrack[theIsolEle1_GSFTrack];

      // search for PV
      closestPV = -1;   // chiara: i mu hanno 2 def x W e Z
      float z0 = 0.0;
      if(nPV>0) {         // a primary vertex was found
        float minDzPV=999.;
        for(int v=0; v<nPV; v++) {
          if(fabs(PVzPV[v]-EWKz)<minDzPV) {
            minDzPV=fabs(PVzPV[v]-EWKz);
            closestPV = v;
          }
        }
        z0 = PVzPV[closestPV];    // if nPV=0 --> z0 = 0

        // loop on tight and loose list
        for(int eleList=0; eleList<2; eleList++) {

          int nIsolEles = (eleList==0) ? howManyIsolTightEles : howManyIsolLooseEles;

          // electron consistency with primary vertex
          for(int iisol=0; iisol<nIsolEles; iisol++) {
            bool isGoodEle = true;

            int iele = (eleList==0) ? isolTightElectrons[iisol] : isolLooseElectrons[iisol];
            int gsfTrack = gsfTrackIndexEle[iele];

            // zele - zPV
            float dz = trackVzGsfTrack[gsfTrack]-PVzPV[closestPV];
            if(_commonSel->getSwitch("dzVertex") && !_commonSel->passCut("dzVertex",dz)) isGoodEle = false;
            else {
              if(eleList==0) howManyDzVertexTightEles++;
              if(eleList==1) howManyDzVertexLooseEles++;
            }

            // dxy significance (used if dxy not used)
//             float dxySign = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], eleTrackVxEle[iele], eleTrackVyEle[iele], eleTrackVzEle[iele], pxEle[iele], pyEle[iele], pzEle[iele]) / eleTrackDxyErrorEle[iele];
//             if(_commonSel->getSwitch("dxySigVertex") && !_commonSel->passCut("dxySigVertex",dxySign)) isGoodEle = false;
            
            // dxy
            float dxy = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxGsfTrack[gsfTrack], trackVyGsfTrack[gsfTrack], trackVzGsfTrack[gsfTrack], pxEle[iele], pyEle[iele], pzEle[iele]);
            if(_commonSel->getSwitch("dxyVertex") && !_commonSel->passCut("dxyVertex",dxy)) isGoodEle = false;

            if(isGoodEle) {
              if(eleList==0) vertexTightElectrons.push_back(iele); 
              else if(eleList==1) vertexLooseElectrons.push_back(iele);
            }
          }
        }
      }

      int howManyTightVertexEles = vertexTightElectrons.size();
      int howManyLooseVertexEles = vertexLooseElectrons.size();

      // total number of 'good' electrons
      int howManyEles = -999;
      int howManyTightEles = -999;
      if (!_commonSel->getSwitch("asymmetricElectrons")) { 
        howManyEles      = howManyTightVertexEles; 
        howManyTightEles = howManyTightVertexEles; 
        howManyDzVertexEles = howManyDzVertexTightEles; 
      }
      if (_commonSel->getSwitch("asymmetricElectrons"))  { 
        howManyTightEles = howManyTightVertexEles; 
        howManyEles = howManyTightVertexEles+howManyLooseVertexEles; 
        howManyDzVertexEles = howManyDzVertexTightEles + howManyDzVertexLooseEles; 
      }   

      // two highest pt good electrons - no charge requirement
      std::pair<int,int> bestElectronTightPair = getBestGoodElePair(vertexTightElectrons);
      int theEle1(bestElectronTightPair.first);
      int theEle2(bestElectronTightPair.second);
      getBestGoodElePairFunny(isolTightElectrons);

      if (_commonSel->getSwitch("asymmetricElectrons") && theEle2<0 ){
        std::pair<int,int> bestElectronLoosePair = getBestGoodElePair(vertexLooseElectrons);
        int myEle2(bestElectronLoosePair.first);
        theEle2 = myEle2;
      }

      // set the event kinematics
      setKinematics(theEle1, theEle2);

      //PM Forming pfjets removing tracks in a veto cone around theEle1&2
//       _theCaloTowersForPFJets.clear();
//       getCaloTowersForPFJets(0.5,500.,20.,2.4,5,0.04, 0.1, 0.1, closestPV, theEle1, theEle2);  

//       // making jets
//       _theSISConePFJets.clear();

//       nfound = 0;
//       iCT = 0;
      
//       while(nfound == 0 && iCT < int(_theCaloTowersForPFJets.size())) {
//         if(_theCaloTowersForPFJets[iCT].et() > 0.5) nfound = 1;
//         iCT++;
//       }

//       if(nfound == 1) {
// 	_theSISConePFJets = SortJet(SISCone(_theCaloTowersForPFJets, 0.5, 0.0));
//       }

      // take from the event the PF reconstructed jets 
      _theAK5PFJets = GetUncorrPFJets();
      
      // the basic skim is that there is >=1 electron identified
      if( writeSkim ) {
        if ( howManyIdEle > 0 ) skimFile << 1 << std::endl;
        else skimFile << 0 << std::endl;
      }

      // for optimization studies
      if(_commonSel->getSwitch("dumpIsolation")) FillIsolationVertexTree(closestPV, theEleIdent1);
      
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

      //! event shape variables
      vector<float> MHTphiJet = MHTphiJetForW(theEle1_, goodWJets_);
      vector<float> MHTphiPFJet = MHTphiJetForW(theEle1_, goodWPFJets_);
      float MHTphiMET = MHTphiMETForZ(theEle1_, theEle2_);

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

      // study the best electron candidate selection for W
      if ( m_doBestElectronStudy ) {

        if ( howManyEles>0 ) {

          if( WToENuDecay ) {
            float ptMc = pMc[indexMcEleWToENu]*fabs(sin(thetaMc[indexMcEleWToENu]));

            Gene_eta[howManyWJets]->Fill(etaMc[indexMcEleWToENu]);
            Gene_pt[howManyWJets]->Fill(ptMc);

            int bestEleByPt = _bestPairByPt.first;
            float dr = deltaR_MCmatch(8,bestEleByPt);
            if ( electron_MCmatch_DeltaR(8,bestEleByPt,0.2) ) {
              H_etaBestByPt[howManyWJets]->Fill(etaMc[indexMcEleWToENu]);
              H_ptBestByPt[howManyWJets]->Fill(ptMc);
            }
            else {
              // these are inclusive
              H_etaMcNotMatched->Fill(etaMc[indexMcEleWToENu]);
              H_phiMcNotMatched->Fill(phiMc[indexMcEleWToENu]);
              H_ptNotMatched->Fill(ptMc);

              H_nRecoNotMatched->Fill(howManyAccEle);
              H_nIdNotMatched->Fill(howManyIdTightEle);
              H_nIsolNotMatched->Fill(howManyEles);
            }

            H_etaResolution->Fill(etaEle[bestEleByPt]-etaMc[indexMcEleWToENu]);
            H_phiResolution->Fill(phiEle[bestEleByPt]-phiMc[indexMcEleWToENu]);
            H_deltaR->Fill(dr);

          }

        }

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
      commonData.weight = weight;
      commonData.passedHLT = passedHLT;
      commonData.nRecoEle = howManyRecoEle;
      commonData.nAccEle = howManyAccEle;
      commonData.nIdTightEle = howManyIdTightEle;
      commonData.nIdLooseEle = howManyIdLooseEle;
      commonData.nTkIsolTightEle = howManyTkIsolTightEle;
      commonData.nTkIsolLooseEle = howManyTkIsolLooseEle;
      commonData.nEcalIsolTightEle = howManyEcalIsolTightEle;
      commonData.nEcalIsolLooseEle = howManyEcalIsolLooseEle;
      commonData.nHcalIsolTightEle = howManyIsolTightEles;
      commonData.nHcalIsolLooseEle = howManyIsolLooseEles;
      if(_commonSel->getSwitch("combinedIsolationTightEB") && _commonSel->getSwitch("combinedIsolationTightEE")) {
        commonData.nCombinedIsolTightEle = howManyIsolTightEles;
      } else commonData.nCombinedIsolTightEle = -1;
      if(_commonSel->getSwitch("combinedIsolationLooseEB") && _commonSel->getSwitch("combinedIsolationLooseEE")) {
        commonData.nCombinedIsolLooseEle = howManyIsolLooseEles;
      } else commonData.nCombinedIsolLooseEle = -1;
      commonData.nPV = nPV;
      commonData.nDzVertexEle = howManyDzVertexEles;
      commonData.nDxyVertexEle = howManyEles;
      commonData.nMuons = nGoodMuons;

      // Z + calo jets selection
      SelectorData dataZeeJet(commonData);
      if(m_signal == zjets) dataZeeJet.signal = ZToEEDecay;
      else if(m_signal == zother) dataZeeJet.signal = !ZToEEDecay;
      else dataZeeJet.signal = true;
      dataZeeJet.njetsMC = nzjets_mc;
      dataZeeJet.njets = howManyZJets;
      dataZeeJet.btagEVT = jetBTagEVTZ.combinedSecondaryVertexBJetTags;
      dataZeeJet.mInv = mee_;
      dataZeeJet.mhtMET = MHTphiMET;

      CutBasedSelectorEE selectorZeeJet(dataZeeJet);
      selectorZeeJet.isMc(!isData_);
      isZeeUpToEnd  = selectorZeeJet.outputZ(_commonSel,_zeeSel,_zeePCounter,&m_zeeJetbinCounter,&m_zeeRecoJetbinCounter);

      // W + calo jets selection
      SelectorData dataWenuJet(commonData);
      if(m_signal == wjets) dataWenuJet.signal = WToENuDecay;
      else if(m_signal == wother) dataWenuJet.signal = !WToENuDecay; 
      else dataWenuJet.signal = true;
      dataWenuJet.njetsMC = nwjets_mc;
      dataWenuJet.njets = howManyWJets;
      dataWenuJet.btagEVT = jetBTagEVTW.combinedSecondaryVertexBJetTags;
      dataWenuJet.foundAnyZ = ZcandFound;
      dataWenuJet.met = p3Met_.Mag();
      dataWenuJet.mt = WmT_;
      dataWenuJet.mhtJet = MHTphiJet;

      CutBasedSelectorEE selectorWenuJet(dataWenuJet);
      selectorWenuJet.isMc(!isData_);
      isWenuUpToEnd  = selectorWenuJet.outputW(_commonSel,_wenuSel,_wenuPCounter,&m_wenuJetbinCounter,&m_wenuRecoJetbinCounter);

      // Z + PF jets selection
      SelectorData dataZeePFJet(commonData);
      if(m_signal == zjets) dataZeePFJet.signal = ZToEEDecay;
      else if(m_signal == zother) dataZeePFJet.signal = !ZToEEDecay;
      else dataZeePFJet.signal = true;
      dataZeePFJet.njetsMC = nzpfjets_mc;
      dataZeePFJet.njets = howManyZPFJets;
      dataZeePFJet.btagEVT = pfjetBTagEVTZ.combinedSecondaryVertexBJetTags;
      dataZeePFJet.mInv = mee_;
      dataZeePFJet.mhtMET = MHTphiMET; // fix it with the PFjets one!

      CutBasedSelectorEE selectorZeePFJet(dataZeePFJet);
      selectorZeePFJet.isMc(!isData_);
      isZeeUpToEndPFJet  = selectorZeePFJet.outputZ(_commonSel,_zeeSel,0,&m_zeePFJetbinCounter,&m_zeeRecoPFJetbinCounter);

      // W + PF jets selection
      SelectorData dataWenuPFJet(commonData);
      if(m_signal == wjets) dataWenuPFJet.signal = WToENuDecay;
      else if(m_signal == wother) dataWenuPFJet.signal = !WToENuDecay; 
      else dataWenuPFJet.signal = true;
      dataWenuPFJet.njetsMC = nwpfjets_mc;
      dataWenuPFJet.njets = howManyWPFJets;
      dataWenuPFJet.btagEVT = jetBTagEVTW.combinedSecondaryVertexBJetTags;
      dataWenuPFJet.foundAnyZ = ZcandFound;
      dataWenuPFJet.met = p3Met_.Mag();
      dataWenuPFJet.mt = WmT_;
      dataWenuPFJet.mhtJet = MHTphiJet; // fix it with the PFjets one!

      CutBasedSelectorEE selectorWenuPFJet(dataWenuPFJet);
      selectorWenuPFJet.isMc(!isData_);
      isWenuUpToEndPFJet  = selectorWenuPFJet.outputW(_commonSel,_wenuSel,0,&m_wenuPFJetbinCounter,&m_wenuRecoPFJetbinCounter);

      // --------------------------------------------------------------
      if (isWenuUpToEnd || isWenuUpToEndPFJet) {       // filling reduced trees for the W->enu study
	
        myOutTree_Wenu -> fillJetMultiplicities( howManyWJets, howManyWPFJets, -1, nwgenjets_mc, isWenuUpToEnd, isWenuUpToEndPFJet, 0 );
        myOutTree_Wenu -> fillElectrons(recoflag_, pt_, eta_,
                                        classification_, deta_, dphi_, hoe_, see_, eop_, fbrem_,
                                        trackerIso_, hcalIso_, ecalJIso_, ecalGTIso_, combinedIso_);
        myOutTree_Wenu -> fillKinematics( -1., WCalomT_, WTCmT_, WPFmT_, GetPt(pxMet[0],pyMet[0]), GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), pt3W_.Mag(), U_.Mag() );
        myOutTree_Wenu -> fillEventShape( MHTphiJet, MHTphiPFJet, -1.);
        if(!isData_) myOutTree_Wenu->fillMcTruth(WToENuDecay);
        myOutTree_Wenu->fillBTagEVT(jetBTagEVTW.combinedSecondaryVertexBJetTags, jetBTagEVTW.combinedSecondaryVertexMVABJetTags,
                                    jetBTagEVTW.jetBProbabilityBJetTags, jetBTagEVTW.jetProbabilityBJetTags, jetBTagEVTW.simpleSecondaryVertexBJetTags,
                                    jetBTagEVTW.softMuonBJetTags,
                                    jetBTagEVTW.trackCountingHighPurBJetTags, jetBTagEVTW.trackCountingHighEffBJetTags);
        myOutTree_Wenu -> fillRunInfo(runNumber, lumiBlock, eventNumber);
        myOutTree_Wenu -> store();

        // fill the matrix
        if(WToENuDecay) {
          NWjetMatrix->Fill(nwjets_mc,howManyWJets);
          NWPFjetMatrix->Fill(nwpfjets_mc,howManyWPFJets);
          for(int jettrue=0; jettrue<nwjets_mc; jettrue++) {
            for(int jetreco=0; jetreco<howManyWJets; jetreco++) {
              NWjetMatrixIncl->Fill(jettrue,jetreco);
            }
          }
          for(int jettrue=0; jettrue<nwpfjets_mc; jettrue++) {
            for(int jetreco=0; jetreco<howManyWPFJets; jetreco++) {
              NWPFjetMatrixIncl->Fill(jettrue,jetreco);
            }
          }
        }
      }

      if (isZeeUpToEnd || isZeeUpToEndPFJet ) {        // filling reduced trees for the Z->ee study

        vector<float> dummyW;
        for(int i=0;i<6;i++) dummyW.push_back(-999.);

        myOutTree_Zee -> fillJetMultiplicities( howManyZJets, howManyZPFJets, -1, nzgenjets_mc, isZeeUpToEnd, isZeeUpToEndPFJet, 0 );
        myOutTree_Zee -> fillElectrons(recoflag_, pt_, eta_,
                                       classification_, deta_, dphi_, hoe_, see_, eop_, fbrem_,
                                       trackerIso_, hcalIso_, ecalJIso_, ecalGTIso_, combinedIso_);
        myOutTree_Zee -> fillKinematics( mee_, -1., -1., -1., GetPt(pxMet[0],pyMet[0]), GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), -1., -1.);
        myOutTree_Zee -> fillWTemplates( tempMETcorr_, tempMETuncorr_, WmTFromZCorr_, WmTFromZUncorr_);
        myOutTree_Zee -> fillEventShape( dummyW, dummyW, MHTphiMET);
        if(!isData_) myOutTree_Zee->fillMcTruth(ZToEEDecay);
        myOutTree_Zee->fillBTagEVT(jetBTagEVTZ.combinedSecondaryVertexBJetTags, jetBTagEVTZ.combinedSecondaryVertexMVABJetTags,
                                   jetBTagEVTZ.jetBProbabilityBJetTags, jetBTagEVTZ.jetProbabilityBJetTags, jetBTagEVTZ.simpleSecondaryVertexBJetTags,
                                   jetBTagEVTZ.softMuonBJetTags,
                                   jetBTagEVTZ.trackCountingHighPurBJetTags, jetBTagEVTZ.trackCountingHighEffBJetTags);
        myOutTree_Zee -> fillRunInfo(runNumber, lumiBlock, eventNumber);
        myOutTree_Zee -> store();

        // fill the matrix
        if(ZToEEDecay) {
          NZjetMatrix->Fill(nzpfjets_mc,howManyZJets);
          NZPFjetMatrix->Fill(nzjets_mc,howManyZPFJets);
          for(int jettrue=0; jettrue<nzjets_mc; jettrue++) {
            for(int jetreco=0; jetreco<howManyZJets; jetreco++) {
              NZjetMatrixIncl->Fill(jettrue,jetreco);
            }
          }
          for(int jettrue=0; jettrue<nzpfjets_mc; jettrue++) {
            for(int jetreco=0; jetreco<howManyZPFJets; jetreco++) {
              NZPFjetMatrixIncl->Fill(jettrue,jetreco);
            }
          }
        }
      }

    }

  } // loope over entries

  if( writeSkim ) {
    skimFile.close();
  }

  // save the Njet matrix
  char nameFileMatrix[200];
  sprintf(nameFileMatrix,"%sJetMatrix.root",m_prefix);
  TFile file(nameFileMatrix,"recreate");
  file.cd();
  NWjetMatrix->Write();
  NWjetMatrixIncl->Write();
  NZjetMatrix->Write();
  NZjetMatrixIncl->Write();
  NWPFjetMatrix->Write();
  NWPFjetMatrixIncl->Write();
  NZPFjetMatrix->Write();
  NZPFjetMatrixIncl->Write();
  file.Close();

}  


void VecbosEESelectionToCompare::displayEfficiencies() {

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
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Tight ElectronID selections: " << std::endl;
  EgammaTightID_NC.diplayEfficiencies();
  std::cout << endl;
  std::cout << endl;
  std::cout << "----------------------------------------"      << std::endl;
  std::cout << "Loose ElectronID selections: " << std::endl;
  EgammaLooseID_NC.diplayEfficiencies();
  std::cout << endl;
  std::cout << endl;
  std::cout << "----------------------------------------"      << std::endl;

  if ( m_doBestElectronStudy ) {

    for(int j=0; j<5; j++) {

      const char* opt = (j==0) ? "RECREATE" : "UPDATE";
      sprintf(namefile,"%sVecBosOutput-WenuJ_purity.root",m_prefix);
      EfficiencyEvaluator PurityEta(namefile,opt);
      PurityEta.AddNumerator(Gene_eta[j]);
      PurityEta.AddNumerator(H_etaBestByPt[j]);
      PurityEta.SetDenominator(Gene_eta[j]);
      PurityEta.ComputeEfficiencies(false);
      PurityEta.Write();

      EfficiencyEvaluator PurityPt(namefile,"UPDATE");
      PurityPt.AddNumerator(Gene_pt[j]);
      PurityPt.AddNumerator(H_ptBestByPt[j]);
      PurityPt.SetDenominator(Gene_pt[j]);
      PurityPt.ComputeEfficiencies(false);
      PurityPt.Write();

    }

    TFile *fileEff = TFile::Open(namefile,"UPDATE");
    fileEff->cd();
    H_etaResolution->Write();
    H_phiResolution->Write();
    H_deltaR->Write();
    H_etaMcNotMatched->Write();
    H_phiMcNotMatched->Write();
    H_ptNotMatched->Write();
    H_nRecoNotMatched->Write();
    H_nIdNotMatched->Write();
    H_nIsolNotMatched->Write();
    fileEff->Close();

  }

}

// two highest pT 'good' electrons
std::pair<int,int> VecbosEESelectionToCompare::getBestGoodElePair(std::vector<int> goodElectrons) {
  if (goodElectrons.size()<1)
    {
      return make_pair(-1,-1);
    }  
  else if (goodElectrons.size()<2)
    {
      return make_pair(goodElectrons[0],-1);
    }
  else
    {
      float ptEle[100];
      
      for(int iEle=0;iEle<goodElectrons.size();iEle++) {
	int eleIndex = goodElectrons[iEle];
	TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
	ptEle[iEle]=pEle.Pt();
      }
      int sortedPtIndex[100];
      for (int i=0;i<100;++i)
	sortedPtIndex[i]=-1;
      TMath::Sort<float,int>(goodElectrons.size(),ptEle,sortedPtIndex);
      return make_pair(goodElectrons[sortedPtIndex[0]],goodElectrons[sortedPtIndex[1]]);
    }
}

void VecbosEESelectionToCompare::setKinematics(int theEle1, int theEle2) {
  // electron candidates
  theEle1_ = theEle1;
  theEle2_ = theEle2;
  float massele = 0.000511;
  if (theEle1_ > -1) {
    TVector3 p1;
    p1.SetXYZ(pxEle[theEle1_],pyEle[theEle1_],pzEle[theEle1_]);
    p4Ele1_.SetVectM(p1,massele);
    pT3Ele1_.SetXYZ(pxEle[theEle1_],pyEle[theEle1_],0.0);
  } else {
    p4Ele1_.SetXYZT(0.,0.,0.,0.);
  }
  if (theEle2 > -1) {
    TVector3 p2;
    p2.SetXYZ(pxEle[theEle2],pyEle[theEle2],pzEle[theEle2]);
    p4Ele2_.SetVectM(p2,massele); 
    pT3Ele2_.SetXYZ(pxEle[theEle2],pyEle[theEle2],0.0);
  } else {
    p4Ele2_.SetXYZT(0.,0.,0.,0.);
  }

  //! set best electron(s) variables
  setElectrons();

  sumCalo_.SetXYZ(0,0,0);
  std::vector<CaloTower>::iterator tower;
  for(tower=_theCaloTowers.begin(); tower!=_theCaloTowers.end(); ++tower) {
    sumCalo_ += tower->Get3Vector();
  }

  // missing transverse energy (used for the cut, even if we dump also the other ones)
  if(_commonSel->getSwitch("useTCMet")) p3Met_.SetXYZ(pxTCMet[0],pyTCMet[0],0.0);
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
  
    // recoil W Pt = -[E_HCAL + E_ECAL - E_SC(ele)]
    TLorentzVector electronSCMomentum = getSC4Vector(theEle1_);
    TVector3 theSC3Vector = electronSCMomentum.Vect();
    TVector3 U_ = sumCalo_ - theSC3Vector; 

    // U_par (U_perp) = projections of U along (perp) to electron track
    U_par_ = pT3Ele1_ * U_ / pT3Ele1_.Mag();
    U_perp_ = U_.Perp(pT3Ele1_);
  
    // W transverse mass (used in the selection)
    WmT_ = sqrt(2 * p4Ele1_.Pt() * p3Met_.Mag() * (1-cos(pT3Ele1_.Angle(p3Met_))) );
    // to dump in final tree
    WCalomT_ = sqrt(2 * p4Ele1_.Pt() * p3CaloMet_.Mag() * (1-cos(pT3Ele1_.Angle(p3CaloMet_))) );
    WTCmT_ = sqrt(2 * p4Ele1_.Pt() * p3TCMet_.Mag() * (1-cos(pT3Ele1_.Angle(p3TCMet_))) );   
    WPFmT_ = sqrt(2 * p4Ele1_.Pt() * p3PFMet_.Mag() * (1-cos(pT3Ele1_.Angle(p3PFMet_))) );   
  } else {
    pt3W_.SetXYZ(0.,0.,0.);
    U_.SetXYZ(0.,0.,0.);
    U_par_ = U_perp_ = WmT_ = -1;
  }
  
  // -------------------------------------------------------------------------------
  // Z SPECIFIC VARIABLES
  // -------------------------------------------------------------------------------
  if(theEle1_ > -1 && theEle2_ > -1) {
    // invariant mass
    mee_ = getMee(theEle1_, theEle2_);

    //! --- template mT from Z ---
    TVector3 ele2SCp3 = getSC4Vector(theEle2_).Vect();
    TVector3 newUncorr = sumCalo_ - ele2SCp3; 
    TVector3 newUncorrPerp(newUncorr.x(), newUncorr.y(), 0.);
    tempMETuncorr_ = newUncorr.Perp();

    // Z direction for the correction
    pt3Z_.SetXYZ(p4Ele1_.Px()+p4Ele2_.Px(), p4Ele1_.Py()+p4Ele2_.Py(), 0.0);

    // correction for the mW/mZ ratio
    double zTrueMass = 91.1876;
    double wTrueMass = 80.403;	  
    double ratio = wTrueMass/zTrueMass;
    TVector3 newCorrPerp = ( (newUncorrPerp-pt3Z_)*ratio ) + pt3Z_;
    tempMETcorr_ = newCorrPerp.Mag();

    // W transverse mass from template MET
    WmTFromZUncorr_ = sqrt(2*(pT3Ele1_.Mag())*tempMETuncorr_*(1-cos(pT3Ele1_.Angle(newUncorrPerp))) );
    WmTFromZCorr_   = sqrt(2*(pT3Ele1_.Mag())*tempMETcorr_*(1-cos(pT3Ele1_.Angle(newCorrPerp))) );
  } else {
    pt3Z_.SetXYZ(0.,0.,0);
    mee_ = WmTFromZUncorr_ = WmTFromZCorr_ = -1;
  }

}

void VecbosEESelectionToCompare::setElectrons() {
  
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
      recoflag_[i] = recoFlagsEle[ele];
      pt_[i] = p4[i].Pt();
      eta_[i] = p4[i].Eta();
      classification_[i] = classificationEle[ele];
      deta_[i] = deltaEtaAtVtxEle[ele];
      dphi_[i] = deltaPhiAtVtxEle[ele];
      hoe_[i] = hOverEEle[ele];
      see_[i] = SigmaiEiE(ele);
      eop_[i] = eSuperClusterOverPEle[ele];
      fbrem_[i] = FBrem(ele);
      trackerIso_[i] = dr03TkSumPtEle[ele];
      hcalIso_[i] = dr04HcalTowerSumEtEle[ele];
      ecalJIso_[i] = dr04EcalRecHitSumEtEle[ele];
      ecalGTIso_[i] = scBasedEcalSum04Ele[ele];
      combinedIso_[i] = combinedIsolation[ele];
    } else { // the electron is not selected
      recoflag_[i] = -1;
      pt_[i] = -1;
      eta_[i] = -1;
      classification_[i] = -1;
      deta_[i] = -1;
      dphi_[i] = -1;
      hoe_[i] = -1;
      see_[i] = -1;
      eop_[i] = -1;
      fbrem_[i] = -1;
      trackerIso_[i] = -1;
      hcalIso_[i] = -1;
      ecalJIso_[i] = -1;
      ecalGTIso_[i] = -1;
      combinedIso_[i] = -1;
    }
  }

}

void VecbosEESelectionToCompare::getBestGoodElePairFunny(std::vector<int> goodElectrons) {

  std::vector<ElectronQualityData> electronQual;
  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    ElectronQualityData quality;
    quality.index = eleIndex;
    quality.ptAtVtx = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
    
    int sclu = superClusterIndexEle[eleIndex];
    int pfsclu = PFsuperClusterIndexEle[eleIndex];
    float scEnergy = -1;
    // if only tracker driven, take the PF supercluster, take the standard one otherwise
    Utils anaUtils;
    if(anaUtils.electronRecoType(recoFlagsEle[eleIndex],isEcalDriven)) {
      // if it is ECAL driven it must have the SC
      scEnergy = energySC[sclu];
    } else {
      // if it has the PF supercluster, use that, else use the tracker momentum at outermost layer
      if(pfsclu > -1) scEnergy = energyPFSC[pfsclu];
      else {
        int gsfTrack = gsfTrackIndexEle[eleIndex];
        TVector3 p3T(pxAtOuterTrack[gsfTrack],pyAtOuterTrack[gsfTrack],pzAtOuterTrack[gsfTrack]);
        scEnergy = p3T.Mag();
      }
    }

    quality.SCenergy = scEnergy;
    quality.trackerSumPt = dr03TkSumPtEle[eleIndex];
    quality.hcalSumEt = dr04HcalTowerSumEtEle[eleIndex];
    quality.electronIdLH = eleIdLikelihoodEle[eleIndex];
    electronQual.push_back(quality);
  }
  ElectronBestCandidateSelector selector(electronQual);

  _bestPairByPt = selector.bestByPt();
  _bestPairBySCenergy = selector.bestBySCenergy();
  _bestPairByTrackerIsolation = selector.bestByTrackerIsolation();
  _bestPairByHcalIsolation = selector.bestByEcalIsolation();
  _bestPairByElectronIdLH = selector.bestByElectronIdLH();

}

// invariant mass
float VecbosEESelectionToCompare::getMee(int theEle1, int theEle2) {
  TLorentzVector *_p1 = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *_p2 = new TLorentzVector(0.,0.,0.,0.);
  _p1 ->SetXYZT(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1],energyEle[theEle1]);
  _p2 ->SetXYZT(pxEle[theEle2],pyEle[theEle2],pzEle[theEle2],energyEle[theEle2]);
  double theInvMass = (*_p1+*_p2).M();       
  return theInvMass;
}

std::vector<Jet> VecbosEESelectionToCompare::buildSortedSISConeCaloJets() { 

  // calotowers and jets with respect to the origin
  _theCaloTowers.clear();
  _theCaloTowers = getCaloTowers(0.0, 0);
  
  // making jets
  _theSISConeCaloJets.clear();
  _theSISConeBTagCaloJets.clear();

  int nfound = 0;
  int iCT = 0;
  while(nfound == 0 && iCT < int(_theCaloTowers.size())) {
    if(_theCaloTowers[iCT].et() > 0.5) nfound = 1;
    iCT++;
  }

  if(nfound == 1) {
    _theSISConeCaloJets = SortJet(SISCone(_theCaloTowers, 0.5, 0.0));
  }

}

bool VecbosEESelectionToCompare::isEleID(CutBasedEleIDSelector *selector, int eleIndex) {
  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 pTrkAtOuter(pxAtOuterGsfTrack[gsf],pyAtOuterGsfTrack[gsf],pzAtOuterGsfTrack[gsf]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = FBrem(eleIndex);
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
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
  selector->SetElectronClass ( classificationEle[eleIndex] );
  selector->SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  selector->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  
  //  return selector->output(); // class dependent result
  return selector->outputNoClass();
}

void VecbosEESelectionToCompare::FillEleIDOptimTree(int eleIndex) {

  TVector3 p3Ele(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
  int gsftrack = gsfTrackIndexEle[eleIndex];
  TVector3 p3OutEle(pxAtOuterGsfTrack[gsftrack],pyAtOuterGsfTrack[gsftrack],pzAtOuterGsfTrack[gsftrack]);

  float fbrem = (p3Ele.Mag() - p3OutEle.Mag()) / p3Ele.Mag();

  float s9s25, see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      s9s25 = 999.;
      see = 999.;
    }
  }

  myOutEleIDOptimTree_Wenu->fillAll(classificationEle[eleIndex],
                                    recoFlagsEle[eleIndex],
                                    etaEle[eleIndex],
                                    GetPt(pxEle[eleIndex],pyEle[eleIndex]),
                                    deltaEtaAtVtxEle[eleIndex], 
                                    deltaPhiAtVtxEle[eleIndex], 
                                    hOverEEle[eleIndex],
                                    s9s25, 
                                    see,
                                    eSeedOverPoutEle[eleIndex],
                                    fbrem);
  
  myOutEleIDOptimTree_Wenu->store();

}

void VecbosEESelectionToCompare::FillIsolationVertexTree(int closestPV, int iele) {
  // isolation and vertex optimization (variables are for the hardest ID electron)
  if(closestPV > -1) { // a primary vertex was found
    if(iele > -1) {
      TVector3 p3Ele(pxEle[iele],pyEle[iele],pzEle[iele]);
      int gsfTrack = gsfTrackIndexEle[iele];
      float tkIsol = dr03TkSumPtEle[iele];
      float ecalJurIsol = dr04EcalRecHitSumEtEle[iele];
      float ecalGTIsol = scBasedEcalSum04Ele[iele];
      float hcalIsol = dr04HcalTowerSumEtEle[iele];
      float combinedIsol = combinedIsolation[iele];
      float dz = trackVzGsfTrack[gsfTrack]-PVzPV[closestPV];
      float dxy = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxGsfTrack[gsfTrack], trackVyGsfTrack[gsfTrack], trackVzGsfTrack[gsfTrack], pxEle[iele], pyEle[iele], pzEle[iele]);
      float dxyErr = transvImpactParErrorTrack[gsfTrack];
      myOutIsolationVtxOptimTree_Wenu->fillAll(WToENuDecay,tkIsol,hcalIsol,ecalJurIsol,ecalGTIsol,combinedIsol,dz,dxy,dxyErr);
      myOutIsolationVtxOptimTree_Wenu->fillKinematics(p3Ele.Pt(),p3Ele.Eta());
      myOutIsolationVtxOptimTree_Wenu->store();
    }
  }
}

std::vector<CaloTower> VecbosEESelectionToCompare::getCaloTowers(float zpv, int type) {

  // RecHit thresholds for CaloTowers
  vector<float> thresh;
  thresh.push_back(.9);             //HB
  thresh.push_back(99999999.);      //etc...
  thresh.push_back(1.4);   
  thresh.push_back(1.4);
  thresh.push_back(1.2);
  thresh.push_back(1.8);
  thresh.push_back(.09);
  thresh.push_back(.45);
  thresh.push_back(.2);
  thresh.push_back(.45);

  std::vector<CaloTower> calotowers = CreateCaloTowers(thresh, zpv, type);

  return calotowers;

}

void VecbosEESelectionToCompare::getCaloTowersForPFJets(float cptMin, float cptMax, float cchi2, float cetaMax, float cnHits, float cDxy, float cdZ, float cd3d, int iBestPV, int theEle1, int theEle2) 
{
      TVector3 vPV(PVxPV[iBestPV],PVyPV[iBestPV],PVzPV[iBestPV]);
      int i1=-1; if (theEle1>-1) i1=trackIndexEle[theEle1];
      int i2=-1; if (theEle2>-1) i2=trackIndexEle[theEle2];

      for(int i = 0; i < nTrack; i++){
	TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
        TVector3 v(pxTrack[i],pyTrack[i],pzTrack[i]);
	vT = vT-vPV;
	//remove reco electrons from track collection
	if(i == i1) continue;
	if(i == i2) continue;
	if(vT.Pt() > cDxy) continue;
	if(fabs(vT.z()) > cdZ) continue;
	if(fabs(vT.Mag()) > cd3d) continue;
	double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
	if(pt < cptMin) continue;
	if(pt > cptMax) continue;
	if(trackNormalizedChi2Track[i] > cchi2) continue;
	if(fabs(v.Eta()) > cetaMax) continue;
	if(trackValidHitsTrack[i] < cnHits) continue;
	_theCaloTowersForPFJets.push_back(CaloTower(pt*cosh(v.Eta()), 0.0, v, v, v));
      }
}

void VecbosEESelectionToCompare::ConfigCommonSelections(Selection* _selection) {

  _selection->addSwitch("isData");
  _selection->addSwitch("goodRunLS");
  _selection->addSwitch("eleId");
  _selection->addSwitch("asymmetricElectrons");
  _selection->addSwitch("useTCMet");
  _selection->addSwitch("usePFMet");
  _selection->addSwitch("dumpIsolation");
  _selection->addSwitch("dumpEleID");
  _selection->addCut("eventRange");
  _selection->addCut("etaElectronAcc");
  _selection->addCut("ptElectronAcc");
  _selection->addCut("trackerPtSumTightEB");
  _selection->addCut("trackerPtSumLooseEB");
  _selection->addCut("trackerPtSumTightEE");
  _selection->addCut("trackerPtSumLooseEE");
  _selection->addCut("hcalPtSumTightEB");
  _selection->addCut("hcalPtSumLooseEB");
  _selection->addCut("hcalPtSumTightEE");
  _selection->addCut("hcalPtSumLooseEE");
  _selection->addCut("ecalPtSumTightEB");
  _selection->addCut("ecalPtSumLooseEB");
  _selection->addCut("ecalPtSumTightEE");
  _selection->addCut("ecalPtSumLooseEE");
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

void VecbosEESelectionToCompare::ConfigZeeSelections(Selection* _selection, Counters *_counter) {
  
  _selection->addSwitch("trigger");

  _selection->addCut("nRecoEles");
  _selection->addCut("nAccEles");
  _selection->addCut("nIdEles");
  _selection->addCut("nTkIsolEles");
  _selection->addCut("nEcalIsolEles");
  _selection->addCut("nHcalIsolEles");
  _selection->addCut("nCombinedIsolEles");
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
  _counter->AddVar("nTkIsolEles");
  _counter->AddVar("nEcalIsolEles");
  _counter->AddVar("nHcalIsolEles");
  _counter->AddVar("nCombinedIsolEles");
  _counter->AddVar("nPrimaryVertices");
  _counter->AddVar("nDzVertexEles");
  _counter->AddVar("nDxyVertexEles");
  _counter->AddVar("nMuons");
  _counter->AddVar("btagEVT");
  _counter->AddVar("meeCut");
  _counter->AddVar("MHTphiMET");
  _counter->AddVar("fullSelection");
}

void VecbosEESelectionToCompare::ConfigWenuSelections(Selection* _selection, Counters *_counter) {

  _selection->addSwitch("trigger");
  
  _selection->addCut("nRecoEles");
  _selection->addCut("nAccEles");
  _selection->addCut("nIdEles");
  _selection->addCut("nTkIsolEles");
  _selection->addCut("nEcalIsolEles");
  _selection->addCut("nHcalIsolEles");
  _selection->addCut("nCombinedIsolEles");
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
  _counter->AddVar("nTkIsolEles");
  _counter->AddVar("nEcalIsolEles");
  _counter->AddVar("nHcalIsolEles");
  _counter->AddVar("nCombinedIsolEles");
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
 
void VecbosEESelectionToCompare::displayZeeEfficiencies(Counters _counter){ 
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
  if(_commonSel->getSwitch("combinedIsolationTightEB") && _commonSel->getSwitch("combinedIsolationTightEE") &&
     _commonSel->getSwitch("combinedIsolationLooseEB") && _commonSel->getSwitch("combinedIsolationLooseEE")) {
    _counter.Draw("nCombinedIsolEles","nIdEles");
    _counter.Draw("nPrimaryVertices","nCombinedIsolEles");
  } else {
    _counter.Draw("nTkIsolEles","nIdEles");
    _counter.Draw("nEcalIsolEles","nTkIsolEles");
    _counter.Draw("nHcalIsolEles","nEcalIsolEles");
    _counter.Draw("nPrimaryVertices","nHcalIsolEles");
  }
  _counter.Draw("nDzVertex","nPrimaryVertices");
  _counter.Draw("nDxyVertex","nDzVertex");
  _counter.Draw("nMuons","nDxyVertex");
  _counter.Draw("btagEVT","nMuons");
  _counter.Draw("meeCut","btagEVT");
  _counter.Draw("MHTphiMET","meeCut");
  _counter.Draw("fullSelection","event");
} 

void VecbosEESelectionToCompare::displayWenuEfficiencies(Counters _counter){ 
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
  if(_commonSel->getSwitch("combinedIsolationTightEB") && _commonSel->getSwitch("combinedIsolationTightEE") &&
     _commonSel->getSwitch("combinedIsolationLooseEB") && _commonSel->getSwitch("combinedIsolationLooseEE")) {
    _counter.Draw("nCombinedIsolEles","nIdEles");
    _counter.Draw("ZVeto","nCombinedIsolEles");
  } else {
    _counter.Draw("nTkIsolEles","nIdEles");
    _counter.Draw("nEcalIsolEles","nTkIsolEles");
    _counter.Draw("nHcalIsolEles","nEcalIsolEles");
    _counter.Draw("ZVeto","nHcalIsolEles");
  }
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

void VecbosEESelectionToCompare::bookBestCandidateHistos() {
  
  char histoname[200];
  
  for(int j=0; j<5; j++) {
    
    sprintf(histoname,"Gene_eta_%d",j);
    TH1F* Gene_eta_j      = new TH1F(histoname,      "generated #eta",                  40, -2.5,2.5);
    Gene_eta.push_back(Gene_eta_j);
    
    sprintf(histoname,"H_etaBestByPt_%d",j);
    TH1F* H_etaBestByPt_j = new TH1F(histoname, "#eta of best electrons by p_{T}", 40, -2.5,2.5);
    H_etaBestByPt.push_back(H_etaBestByPt_j);

    sprintf(histoname,"Gene_pt_%d",j);
    TH1F* Gene_pt_j      = new TH1F(histoname, "generated p_{T} (GeV)",  40, 15, 100);
    Gene_pt.push_back(Gene_pt_j);

    sprintf(histoname,"H_ptBestByPt_%d",j);
    TH1F* H_ptBestByPt_j = new TH1F(histoname, "p_{T} of best electrons by p_{T}", 40, 15, 100);
    H_ptBestByPt.push_back(H_ptBestByPt_j);

  }

  H_etaResolution = new TH1F("H_etaResolution", "#eta_{vtx} - #eta_{true}", 50, -0.006, 0.006);
  H_phiResolution = new TH1F("H_phiResolution", "#phi_{vtx} - #phi_{true}", 50, -0.01, 0.01);
  H_deltaR = new TH1F("H_deltaR","H_deltaR", 50, 0, 2*TMath::Pi());
  H_etaMcNotMatched = new TH1F("H_etaMcNotMatched", "#eta of generator electron not matched", 50, -2.5, 2.5);
  H_phiMcNotMatched = new TH1F("H_phiMcNotMatched", "#phi of generator electron not matched", 50, -TMath::Pi(), TMath::Pi());
  H_ptNotMatched = new TH1F("H_ptNotMatched", "p_{T} of generator electron not matched", 50, 0, 100);
  
  H_nRecoNotMatched = new TH1F("H_nRecoNotMatched", "# of reco electrons when not matched", 5, 0, 6);
  H_nIdNotMatched = new TH1F("H_nIdNotMatched", "# of id electrons when not matched", 5, 0, 6);
  H_nIsolNotMatched = new TH1F("H_nIsolNotMatched", "# of isol electrons when not matched", 5, 0, 6);

}

bool VecbosEESelectionToCompare::electron_MCmatch_DeltaR(int iMc, int iEle, float deltaR) {
  
  TVector3 pMcParticle;
  pMcParticle.SetMagThetaPhi(pMc[iMc], thetaMc[iMc], phiMc[iMc]);
  
  TVector3 pRecoElectron(pxEle[iEle],pyEle[iEle],pzEle[iEle]);

  float dr = pMcParticle.DeltaR(pRecoElectron);

  if ( dr < deltaR ) return true;
  else return false;

}

float VecbosEESelectionToCompare::deltaR_MCmatch(int iMc, int iEle) {
  
  TVector3 pMcParticle;
  pMcParticle.SetMagThetaPhi(pMc[iMc], thetaMc[iMc], phiMc[iMc]);
  
  TVector3 pRecoElectron(pxEle[iEle],pyEle[iEle],pzEle[iEle]);

  float dr = pMcParticle.DeltaR(pRecoElectron);

  return dr;

}

bool VecbosEESelectionToCompare::foundZCandidate(int eleIndex, std::vector<int> accEles) {


  int howManyAccEle = accEles.size();
  for(int iAccEle=0; iAccEle<howManyAccEle; iAccEle++) {
    int iAccEleIndex = accEles[iAccEle];
    float mass = getMee(eleIndex,iAccEleIndex);
    if(_zeeSel->passCut("meeCut",mass)) return true;
  }
  
  return false;
}

int VecbosEESelectionToCompare::CountMuons(float etaMax, float ptMin) {

  int nMuons = 0;
  for(int imu=0; imu<nMuon; imu++) {

    // acceptance
    if(fabs(etaMuon[imu])>etaMax) continue;
    if(GetPt(pxMuon[imu],pyMuon[imu])<ptMin) continue;
    
    nMuons++;

  }

  return nMuons;

}

float VecbosEESelectionToCompare::GetDiJetHeaviestMass(std::vector<Jet> cleanedJets) {

  float heaviestMass = -1.;
  // make the combinatorics
  for(unsigned int j1=0; j1<cleanedJets.size(); j1++) {
    TLorentzVector p4Jet1 = cleanedJets[j1].Get4Vector();
    for(unsigned int j2=j1+1; j2<cleanedJets.size(); j2++) {
      TLorentzVector p4Jet2 = cleanedJets[j2].Get4Vector();
      float mass = (p4Jet1 + p4Jet2).M();
      if ( mass > heaviestMass ) heaviestMass = mass;
    }
  }

  return heaviestMass;

}

bool VecbosEESelectionToCompare::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 v(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = v.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(v.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

std::vector<float> VecbosEESelectionToCompare::jetBTagVariables(Jet jet) {

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

void VecbosEESelectionToCompare::calcEventBVetoVariables(std::vector<Jet> jets, std::vector<BTagJet> btags) {
  
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


vector<float> VecbosEESelectionToCompare::MHTphiJetForW(int theEle1, std::vector<Jet>& cleanedWJets) {

  float MHTphiJetTemp[6];
  for(int i=0; i<6; i++) MHTphiJetTemp[i] = 0.0;
  

  LoadCorrections();

  int njets = (int)cleanedWJets.size();

  // shape variables for W
  // this variable is undefined for 0 jets. In this case put default=0 because it falls in the accepted range by the cut
  if (njets>0 && theEle1>-1){

    CoolTools* CTW = new CoolTools();
	
    // 1) correcting the met for electrons
    vector<TLorentzVector> CTinput_corr = CTW->Get4Vectors(CTW->CaloTowers2Jets(_theCaloTowers, 0));   
	
    TLorentzVector wEleLV;
    wEleLV.SetPxPyPzE(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1],energyEle[theEle1]);

    TLorentzVector minus_wEleLV;
    minus_wEleLV.SetPxPyPzE(-pxEle[theEle1],-pyEle[theEle1],-pzEle[theEle1],energyEle[theEle1]);
    CTinput_corr.push_back(minus_wEleLV);

    // this is the MET corrected for the electron-from-W contribution
    TLorentzVector METcorrForWele = (CreateMET(CTinput_corr)).metVector();
	
	
    // 2) use the candle correction to estimated the corrected MET direction (for the boost)
    double METpar  = METcorrForWele.Vect().Dot(wEleLV.Vect())/wEleLV.Pt();
    double METperp = sqrt(METcorrForWele.Pt()*METcorrForWele.Pt()-METpar*METpar);

    // do the corrections
    double METpar_corr  = METpar + METcorrForWele.Pt()*ParCorr(METpar, wEleLV.Pt());
    double METperp_corr = METperp + METcorrForWele.Pt()*PerpCorr(METperp, wEleLV.Pt());

    TVector3 vMETPar;
    TVector3 vMETPerp;
    TVector3 vMET = METcorrForWele.Vect();

    // get the signs of the two corrected components right (sometimes the corrections can flip the direction!)
    if(METpar >= 0){
      vMETPar.SetPtEtaPhi(fabs(METpar), 0.0, wEleLV.Phi());
    } else {
      vMETPar.SetPtEtaPhi(fabs(METpar), 0.0, wEleLV.Phi()+TMath::Pi());
    }
    vMETPerp = vMET-vMETPar;
	
    TVector3 vMETPar_corr;
    TVector3 vMETPerp_corr;
    TVector3 vMET_corr;
	
    if(METpar_corr >= 0.0){
      vMETPar_corr.SetPtEtaPhi(fabs(METpar_corr), 0.0, wEleLV.Phi());
    } else {
      vMETPar_corr.SetPtEtaPhi(fabs(METpar_corr), 0.0, wEleLV.Phi()+TMath::Pi());
    }
	
    if(METperp_corr >= 0.0){
      vMETPerp_corr.SetPtEtaPhi(fabs(METperp_corr), 0.0, vMETPerp.Phi());
    } else {
      vMETPerp_corr.SetPtEtaPhi(fabs(METperp_corr), 0.0, vMETPerp.Phi()+TMath::Pi());
    }

    // this is now the best estimate of the W momentum
    vMET_corr = vMETPerp_corr+vMETPar_corr;

    TLorentzVector CMET;
    CMET.SetPxPyPzE(vMET_corr.X(), vMET_corr.Y(), 0.0, vMET_corr.Mag());
	

    // 3) now we try to use the Mt constraint to estimate the 'W energy' so that we can get the boost (needs momentum and energy)	
    float Wtruemass = 80.403;
    double MW = Wtruemass;
	
    TLorentzVector Nu;
    Nu.SetPxPyPzE(CMET.Px()-wEleLV.Px(),CMET.Py()-wEleLV.Py(),0.0,
                  sqrt((CMET.Px()-wEleLV.Px())*(CMET.Px()-wEleLV.Px())+
                       (CMET.Py()-wEleLV.Py())*(CMET.Py()-wEleLV.Py())));
    double Mt = sqrt((wEleLV.Pt()+Nu.Pt())*(wEleLV.Pt()+Nu.Pt()) - 
                     (wEleLV.Px()+Nu.Px())*(wEleLV.Px()+Nu.Px()) - 
                     (wEleLV.Py()+Nu.Py())*(wEleLV.Py()+Nu.Py()));

    double Nu_z[2];
	
    TLorentzVector W_from_MET;
    if(Mt > MW){
      double E = sqrt(Wtruemass*Wtruemass+CMET.E()*CMET.E());
      W_from_MET.SetPxPyPzE(CMET.Px(), CMET.Py(), 0.0, E);
    } else {
      double dot = Nu.Px()*wEleLV.Px()+Nu.Py()*wEleLV.Py();
      Nu_z[0] = (MW*MW*wEleLV.Pz()+2*dot*wEleLV.Pz()+wEleLV.E()
                 *sqrt((MW*MW-Mt*Mt)*(MW*MW+2*(Nu.Pt()*wEleLV.Pt()+dot))))/(2*wEleLV.Pt()*wEleLV.Pt());
      Nu_z[1] = (MW*MW*wEleLV.Pz()+2*dot*wEleLV.Pz()-wEleLV.E()
                 *sqrt((MW*MW-Mt*Mt)*(MW*MW+2*(Nu.Pt()*wEleLV.Pt()+dot))))/(2*wEleLV.Pt()*wEleLV.Pt());

      double E1 = sqrt(Wtruemass*Wtruemass+CMET.E()*CMET.E() +(Nu_z[0]+wEleLV.Pz())*(Nu_z[0]+wEleLV.Pz()));			   
      double E2 = sqrt(Wtruemass*Wtruemass+CMET.E()*CMET.E() +(Nu_z[1]+wEleLV.Pz())*(Nu_z[1]+wEleLV.Pz()));
			   
      double E;
      if(E1 < E2){
        E = E1;
        W_from_MET.SetPxPyPzE(CMET.Px(), CMET.Py(), 0.0, E);
      } else {
        E = E2;
        W_from_MET.SetPxPyPzE(CMET.Px(), CMET.Py(), 0.0, E);
      }
    }


    // now we have the boost vector for calotowers
    TVector3 boostW = W_from_MET.BoostVector();
    boostW.SetZ(0.0);
    TLorentzVector MHT = CreateMET(CTW->BoostJets(CTW->CaloTowers2Jets(_theCaloTowers, 0), boostW)).metVector();	
	
    // we have to boost in the opposite direction the electron 
    TVector3 minus_boostW = -boostW;
    minus_boostW.SetZ(0.0);
    TLorentzVector minus_boostedLV_wEle = minus_wEleLV;
    minus_boostedLV_wEle.Boost(minus_boostW);
    
    // boosted mht
    TLorentzVector corrMHT = MHT + minus_boostedLV_wEle;
	
    // there is an angular variable for each jet multiplicity:
    // >=1j includes only the first 1, >=2j only the first 1,2 etc.
    for(int bin=1; bin<6; bin++) {
      if(bin > njets) {
        MHTphiJetTemp[bin] = 0.0;
      } else {
        TVector3 p3JetSum = (0.,0.,0.);
        for (int iJet=0; iJet<bin; iJet++){  
          
          TVector3 p3Jet = cleanedWJets[iJet].Get3Vector();
          
          p3JetSum = p3JetSum + p3Jet;
        }
        p3JetSum.SetZ(0.0);
	
        MHTphiJetTemp[bin] = DeltaPhi(p3JetSum.Phi(),corrMHT.Phi());	

      }
    }
    delete CTW;
  } 
  
  vector<float> MHTphiJet;
  for(int i=0; i<6; i++) MHTphiJet.push_back( MHTphiJetTemp[i] );

  return MHTphiJet;

}

float VecbosEESelectionToCompare::phiJetMETForW(std::vector<Jet>& cleanedWJets, std::vector<TLorentzVector>& tracksForPFJets)
{
  float phiJetMET  = 0.0;
  int njets = (int)cleanedWJets.size();
  if (njets>0) 
    {

      CoolTools* CTZ = new CoolTools();
      vector<float> thresh;
      thresh.push_back(.9);             //HB
      thresh.push_back(99999999.);      //etc...
      thresh.push_back(1.4);   
      thresh.push_back(1.4);
      thresh.push_back(99999999.);
      thresh.push_back(99999999.);
      thresh.push_back(.09);
      thresh.push_back(.45);
      thresh.push_back(.2);
      thresh.push_back(.45);
      std::vector<CaloTower> calotowers = CreateCaloTowers(thresh, 0. , 0);
      vector<TLorentzVector> CTinput_corr = CTZ->Get4Vectors(CTZ->CaloTowers2Jets(calotowers, 0));    
      //MET1
      TLorentzVector MET1 = (CreateMET(CTinput_corr)).metVector();
      //MET2
      TLorentzVector MET2 = (CreateMET(tracksForPFJets)).metVector();
      //Vectorial sum of jets
      TLorentzVector J(0.,0.,0.);
      std::vector<TLorentzVector> jetVector= CTZ->Get4Vectors(cleanedWJets);
      for (int i=0;i<njets;++i)
	J+=jetVector[i];
      phiJetMET = DeltaPhi(J.Phi(),  ( MET1*(1./MET1.Pt())+MET2*(1./MET2.Pt()) ).Phi() );
      delete CTZ;
    }

  return phiJetMET;
}

float VecbosEESelectionToCompare::MHTphiMETForZ(int theEle1, int theEle2) {

  float MHTphiMET      = 0.0;

  // shape variables for Z
  if (theEle1>=0 && theEle2>=0) {

    CoolTools* CTZ = new CoolTools();
	
    // 1) correcting the met for electrons
    vector<TLorentzVector> CTinput_corr = CTZ->Get4Vectors(CTZ->CaloTowers2Jets(_theCaloTowers, 0));   
	
    TLorentzVector zEle1LV;
    zEle1LV.SetPxPyPzE(pxEle[theEle1],pyEle[theEle1],pzEle[theEle1],energyEle[theEle1]);

    TLorentzVector zEle2LV;
    zEle2LV.SetPxPyPzE(pxEle[theEle2],pyEle[theEle2],pzEle[theEle2],energyEle[theEle2]);

    TLorentzVector minus_zEle1LV;
    minus_zEle1LV.SetPxPyPzE(-pxEle[theEle1],-pyEle[theEle1],-pzEle[theEle1],energyEle[theEle1]);
    CTinput_corr.push_back(minus_zEle1LV);

    TLorentzVector minus_zEle2LV;
    minus_zEle2LV.SetPxPyPzE(-pxEle[theEle2],-pyEle[theEle2],-pzEle[theEle2],energyEle[theEle2]);
    CTinput_corr.push_back(minus_zEle2LV);

    // this is the MET corrected for the electrons-from-Z contribution
    TLorentzVector METcorrForZele = (CreateMET(CTinput_corr)).metVector();
	
    // 2) use the Z 4momentum to boost the MET
    TLorentzVector theZLV = zEle1LV + zEle2LV;

    // 3) now we have the boost vector for calotowers
    TVector3 boostZ = theZLV.BoostVector();
    boostZ.SetZ(0.0);
    TLorentzVector MHT = CreateMET(CTZ->BoostJets(CTZ->CaloTowers2Jets(_theCaloTowers, 0), boostZ)).metVector();	

    // we have to boost in the opposite direction the electrons with opposite sign
    TVector3 minus_boostZ = -boostZ;
    minus_boostZ.SetZ(0.0);

    TLorentzVector minus_boostedLV_zEle1 = minus_zEle1LV;
    TLorentzVector minus_boostedLV_zEle2 = minus_zEle2LV;
    minus_boostedLV_zEle1.Boost(minus_boostZ);
    minus_boostedLV_zEle2.Boost(minus_boostZ);
	
    // boosted mht
    TLorentzVector corrMHT = MHT + minus_boostedLV_zEle1 + minus_boostedLV_zEle2;
    corrMHT.SetZ(0.0);
	
    MHTphiMET = DeltaPhi(METcorrForZele.Phi(),corrMHT.Phi());	
    delete CTZ;
  }

  return MHTphiMET;
}


void VecbosEESelectionToCompare::LoadCorrections(){

  TBranch *b_ParConst;
  TBranch *b_PerpConst;
  TBranch *b_Perpbin;
  TBranch *b_Parbin;
  TBranch *b_Mubin;

  TFile *infile = new TFile("Z_calibFall08.root", "READ");
  TTree *tree = (TTree*) infile->Get("data");

  tree->SetBranchAddress("ParConst", ParConst, &b_ParConst);
  tree->SetBranchAddress("PerpConst", PerpConst, &b_PerpConst);
  tree->SetBranchAddress("Perpbin", Perpbin, &b_Perpbin);
  tree->SetBranchAddress("Parbin", Parbin, &b_Parbin);
  tree->SetBranchAddress("Mubin", Mubin, &b_Mubin);

  tree->GetEntry(0);

  infile->Close();

  MubinCenter[0] = 0.0;
  PerpbinCenter[0] = 0.0;
  ParbinCenter[0] = 0.0;
  
  for(int i = 0; i < NMu-1; i++){
    MubinCenter[i] = (Mubin[i]+Mubin[i+1])/2.0;
  }
  MubinCenter[NMu-1] = 2.0*MubinCenter[NMu-2] - MubinCenter[NMu-3];

  for(int i = 0; i < NMET-1; i++){
    PerpbinCenter[i] = (Perpbin[i]+Perpbin[i+1])/2.0;
    ParbinCenter[i] = (Parbin[i]+Parbin[i+1])/2.0;
  }
  PerpbinCenter[NMET-1] = 2.0*PerpbinCenter[NMET-2] - PerpbinCenter[NMET-3];
  ParbinCenter[NMET-1] = 2.0*ParbinCenter[NMET-2] - ParbinCenter[NMET-3];
}
  
double VecbosEESelectionToCompare::ParCorr(double thePar, double theMu){
  int iMET_L = 0;
  int iMET_R = 0;
  int iMu_L = 0;
  int iMu_R = 0;

  if(thePar >= ParbinCenter[NMET-1] && theMu >= MubinCenter[NMu-1]){
    return ParConst[NMET-1][NMu-1];
  }

  if(thePar < ParbinCenter[0] && theMu < MubinCenter[0]){
    return ParConst[0][0];
   
  }

  if(thePar >= ParbinCenter[NMET-1] && theMu < MubinCenter[0]){
    return ParConst[NMET-1][0];
  }

  if(thePar < ParbinCenter[0] && theMu >= MubinCenter[NMu-1]){
    return ParConst[0][NMu-1];
  }

  if(thePar < ParbinCenter[0]){
    
    for(int i = 0; i < NMu-1; i++){
      if(theMu >= MubinCenter[i] && theMu < MubinCenter[i+1]){
	iMu_L = i;
	iMu_R = i+1;
	break;
      }
    }
    double y1 = ParConst[0][iMu_R];
    double y0 = ParConst[0][iMu_L];
    double x1 = MubinCenter[iMu_R];
    double x0 = MubinCenter[iMu_L];

    return (y0+(theMu-x0)*(y1-y0)/(x1-x0));
  }

  if(thePar >= ParbinCenter[NMET-1]){
    
    for(int i = 0; i < NMu-1; i++){
      if(theMu >= MubinCenter[i] && theMu < MubinCenter[i+1]){
	iMu_L = i;
	iMu_R = i+1;
	break;
      }
    }
    double y1 = ParConst[NMET-1][iMu_R];
    double y0 = ParConst[NMET-1][iMu_L];
    double x1 = MubinCenter[iMu_R];
    double x0 = MubinCenter[iMu_L];

    return (y0+(theMu-x0)*(y1-y0)/(x1-x0));
  }

  if(theMu < MubinCenter[0]){
    
    for(int i = 0; i < NMET-1; i++){
      if(thePar >= ParbinCenter[i] && thePar < ParbinCenter[i+1]){
	iMET_L = i;
	iMET_R = i+1;
	break;
      }
    }
    double y1 = ParConst[iMET_R][0];
    double y0 = ParConst[iMET_L][0];
    double x1 = ParbinCenter[iMET_R];
    double x0 = ParbinCenter[iMET_L];

    return (y0+(thePar-x0)*(y1-y0)/(x1-x0));
  }

  if(theMu >= MubinCenter[NMu-1]){
    
    for(int i = 0; i < NMET-1; i++){
      if(thePar >= ParbinCenter[i] && thePar < ParbinCenter[i+1]){
	iMET_L = i;
	iMET_R = i+1;
	break;
      }
    }
    double y1 = ParConst[iMET_R][NMu-1];
    double y0 = ParConst[iMET_L][NMu-1];
    double x1 = ParbinCenter[iMET_R];
    double x0 = ParbinCenter[iMET_L];

    return (y0+(thePar-x0)*(y1-y0)/(x1-x0));
  }

  //bilinear interpolation!!!!!
  
  for(int i = 0; i < NMET-1; i++){
    if(thePar >= ParbinCenter[i] && thePar < ParbinCenter[i+1]){
      iMET_L = i;
      iMET_R = i+1;
      break;
    }
  }
  for(int i = 0; i < NMu-1; i++){
    if(theMu >= MubinCenter[i] && theMu < MubinCenter[i+1]){
      iMu_L = i;
      iMu_R = i+1;
      break;
    }
  }

  double x = thePar;
  double y = theMu;
  double x1 = ParbinCenter[iMET_L];
  double x2 = ParbinCenter[iMET_R];
  double y1 = MubinCenter[iMu_L];
  double y2 = MubinCenter[iMu_R];
  double Q11 = ParConst[iMET_L][iMu_L];
  double Q12 = ParConst[iMET_L][iMu_R];
  double Q21 = ParConst[iMET_R][iMu_L];
  double Q22 = ParConst[iMET_R][iMu_R];
  
  double f = Q11*(x2-x)*(y2-y)/((x2-x1)*(y2-y1)) + 
    Q21*(x-x1)*(y2-y)/((x2-x1)*(y2-y1)) +
    Q12*(x2-x)*(y-y1)/((x2-x1)*(y2-y1)) + 
    Q22*(x-x1)*(y-y1)/((x2-x1)*(y2-y1));

  return f;

}
 
double VecbosEESelectionToCompare::PerpCorr(double thePerp, double theMu){
  int iMET_L = 0;
  int iMET_R = 0;
  int iMu_L = 0;
  int iMu_R = 0;

  if(thePerp >= PerpbinCenter[NMET-1] && theMu >= MubinCenter[NMu-1]){
    return PerpConst[NMET-1][NMu-1];
   
  }

  if(thePerp < PerpbinCenter[0] && theMu < MubinCenter[0]){
    return PerpConst[0][0];
  }

  if(thePerp >= PerpbinCenter[NMET-1] && theMu < MubinCenter[0]){
    return PerpConst[NMET-1][0];
  }

  if(thePerp < PerpbinCenter[0] && theMu >= MubinCenter[NMu-1]){
    return PerpConst[0][NMu-1];
  }

  if(thePerp < PerpbinCenter[0]){
    
    for(int i = 0; i < NMu-1; i++){
      if(theMu >= MubinCenter[i] && theMu < MubinCenter[i+1]){
	iMu_L = i;
	iMu_R = i+1;
	break;
      }
    }
    double y1 = PerpConst[0][iMu_R];
    double y0 = PerpConst[0][iMu_L];
    double x1 = MubinCenter[iMu_R];
    double x0 = MubinCenter[iMu_L];

    return (y0+(theMu-x0)*(y1-y0)/(x1-x0));
  }

  if(thePerp >= PerpbinCenter[NMET-1]){
    
    for(int i = 0; i < NMu-1; i++){
      if(theMu >= MubinCenter[i] && theMu < MubinCenter[i+1]){
	iMu_L = i;
	iMu_R = i+1;
	break;
      }
    }
    double y1 = PerpConst[NMET-1][iMu_R];
    double y0 = PerpConst[NMET-1][iMu_L];
    double x1 = MubinCenter[iMu_R];
    double x0 = MubinCenter[iMu_L];

    return (y0+(theMu-x0)*(y1-y0)/(x1-x0));
  }

  if(theMu < MubinCenter[0]){
    
    for(int i = 0; i < NMET-1; i++){
      if(thePerp >= PerpbinCenter[i] && thePerp < PerpbinCenter[i+1]){
	iMET_L = i;
	iMET_R = i+1;
	break;
      }
    }
    double y1 = PerpConst[iMET_R][0];
    double y0 = PerpConst[iMET_L][0];
    double x1 = PerpbinCenter[iMET_R];
    double x0 = PerpbinCenter[iMET_L];

    return (y0+(thePerp-x0)*(y1-y0)/(x1-x0));
  }

  if(theMu >= MubinCenter[NMu-1]){
    
    for(int i = 0; i < NMET-1; i++){
      if(thePerp >= PerpbinCenter[i] && thePerp < PerpbinCenter[i+1]){
	iMET_L = i;
	iMET_R = i+1;
	break;
      }
    }
    double y1 = PerpConst[iMET_R][NMu-1];
    double y0 = PerpConst[iMET_L][NMu-1];
    double x1 = PerpbinCenter[iMET_R];
    double x0 = PerpbinCenter[iMET_L];

    return (y0+(thePerp-x0)*(y1-y0)/(x1-x0));
  }

  //bilinear interpolation!!!!!
  
  for(int i = 0; i < NMET-1; i++){
    if(thePerp >= PerpbinCenter[i] && thePerp < PerpbinCenter[i+1]){
      iMET_L = i;
      iMET_R = i+1;
      break;
    }
  }
  for(int i = 0; i < NMu-1; i++){
    if(theMu >= MubinCenter[i] && theMu < MubinCenter[i+1]){
      iMu_L = i;
      iMu_R = i+1;
      break;
    }
  }

  double x = thePerp;
  double y = theMu;
  double x1 = PerpbinCenter[iMET_L];
  double x2 = PerpbinCenter[iMET_R];
  double y1 = MubinCenter[iMu_L];
  double y2 = MubinCenter[iMu_R];
  double Q11 = PerpConst[iMET_L][iMu_L];
  double Q12 = PerpConst[iMET_L][iMu_R];
  double Q21 = PerpConst[iMET_R][iMu_L];
  double Q22 = PerpConst[iMET_R][iMu_R];
  
 
  double f = Q11*(x2-x)*(y2-y)/((x2-x1)*(y2-y1)) + 
    Q21*(x-x1)*(y2-y)/((x2-x1)*(y2-y1)) +
    Q12*(x2-x)*(y-y1)/((x2-x1)*(y2-y1)) + 
    Q22*(x-x1)*(y-y1)/((x2-x1)*(y2-y1));



  return f;

}

TLorentzVector VecbosEESelectionToCompare::getSC4Vector(int eleindex) {
  Utils anaUtils;
  int sclu = superClusterIndexEle[eleindex];
  int pfsclu = PFsuperClusterIndexEle[eleindex];
  TLorentzVector electronSCMomentum;
  // if only tracker driven, take the PF supercluster, take the standard one otherwise
  if(anaUtils.electronRecoType(recoFlagsEle[eleindex],isEcalDriven)) {
    // if it is ECAL driven it must have the SC
    electronSCMomentum.SetPtEtaPhiE(energySC[sclu]*fabs(thetaSC[sclu]),etaSC[sclu],phiSC[sclu],energySC[sclu]);          
  } else {
    // if it has the PF supercluster, use that, else use the tracker momentum at outermost layer
    if(pfsclu > -1) electronSCMomentum.SetPtEtaPhiE(energyPFSC[pfsclu]*fabs(thetaPFSC[pfsclu]),
                                                    etaPFSC[pfsclu],phiPFSC[pfsclu],energyPFSC[pfsclu]);
    else {
      int gsfTrack = gsfTrackIndexEle[eleindex];
      TVector3 p3T(pxAtOuterTrack[gsfTrack],pyAtOuterTrack[gsfTrack],pzAtOuterTrack[gsfTrack]);
      float mass = 0.000511;
      electronSCMomentum.SetVectM(p3T,sqrt(p3T.Mag2()+mass*mass));
    }
  }
  return electronSCMomentum;
}

