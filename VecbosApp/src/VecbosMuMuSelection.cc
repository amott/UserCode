#include <string>
#include <iostream>

#include <TTree.h>
#include <TObjString.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/Skimmer.hh"
#include "CommonTools/include/BTagAlgoBits.h"
#include "EgammaAnalysisTools/include/ElectronTrackerIsolation.hh"
#include "EgammaAnalysisTools/include/ElectronCaloIsolation.hh"
#include "EgammaAnalysisTools/include/ElectronBestCandidateSelector.hh"
#include "include/CaloTower.hh"
#include "include/VecbosMuMuSelection.hh"
#include "include/JetCounter.hh"

// uncomment if the trees contain the PU informations
//#define PILEUPINFO

using namespace bits;

VecbosMuMuSelection::VecbosMuMuSelection(TTree *tree) 
  : Vecbos(tree) {
  
  // default do not check on mc-truth
  m_signal = all;

  // common kinematic selections
  std::string theConfigDir       = "config/vecbos/";
  std::string fileCutsCommon     = theConfigDir + "CommonSelectionCutsElectrons.txt";
  std::string fileSwitchesCommon = theConfigDir + "CommonSelectionSwitchesElectrons.txt";
//   std::string fileCutsCommon     = theConfigDir + "CommonSelectionCutsWCandle.txt";
//   std::string fileSwitchesCommon = theConfigDir + "CommonSelectionSwitchesWCandle.txt";

  _commonSel = new Selection(fileCutsCommon,fileSwitchesCommon);
  ConfigCommonSelections(_commonSel);

  // single electron efficiency
  //  EgammaCutBasedID.Configure("config/vecbos/");

  // Apply the JES uncertainties
  TString JESUncertainty(_commonSel->getStringParameter("JESUncertainty"));
  if(JESUncertainty==TString("Up")) {
    std::cout << "=== APPLYING SCALE UP JES UNCERTAINTY" << std::endl;
    jes_=1;
  } else if(JESUncertainty==TString("Down")) {
    std::cout << "=== APPLYING SCALE DOWN JES UNCERTAINTY" << std::endl;
    jes_=-1;
  } else jes_=0;

  //Configuring EleId+Iso WP
  TString selectionString(_commonSel->getStringParameter("electronIDType"));
  TObjArray* selectionTokens=selectionString.Tokenize("x");
  if (selectionTokens->GetEntries()!=2)
    {
      std::cout << "UNKNOWN ELECTRON IDENTIFICATION WORKING POINTS " << selectionTokens->GetEntries() << std::endl;
      exit(1);  
    }

  TString TightIdWP=((TObjString*)(*selectionTokens)[0])->GetString();
  TString LooseIdWP=((TObjString*)(*selectionTokens)[1])->GetString();
  cout << "=== CONFIGURING " << TightIdWP << " TIGHT ELECTRON ID ===" << endl;
  EgammaTightID_NC.ConfigureNoClass("config/vecbos/electronId/"+TightIdWP);
  cout << "=== CONFIGURING " << LooseIdWP << " LOOSE ELECTRON ID ===" << endl;
  EgammaLooseID_NC.ConfigureNoClass("config/vecbos/electronId/"+LooseIdWP);

  EgammaTightID_NC.ConfigureEcalCleaner("config/vecbos/electronId/");
  EgammaLooseID_NC.ConfigureEcalCleaner("config/vecbos/electronId/");

  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_MC.root");
  TDirectory *LHdir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useHoverE = false;
  LH = new ElectronLikelihood(&(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);


  //Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << _commonSel->getSwitch("goodRunLS") << " isData is " <<  _commonSel->getSwitch("isData") << std::endl;

  //To read good run list!
  if (_commonSel->getSwitch("goodRunLS") && _commonSel->getSwitch("isData"))
    {
      std::string goodRunGiasoneFile       = "config/vecbos/json/dataset_eg_Dec22ReReco.json";
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
  for(int ithr=0; ithr<2; ++ithr) {
    for(int jetbin=0; jetbin<6; jetbin++) {
      ConfigZeeSelections(_zeeSel, &m_zeeJetbinCounter[ithr][jetbin]);
      ConfigZeeSelections(_zeeSel, &m_zeeRecoJetbinCounter[ithr][jetbin]);
      ConfigZeeSelections(_zeeSel, &m_zeePFJetbinCounter[ithr][jetbin]);
      ConfigZeeSelections(_zeeSel, &m_zeeRecoPFJetbinCounter[ithr][jetbin]);
      
      char counterTitle[200];
      // with MC truth electrons subtraction
      sprintf(counterTitle,"Z+%djets thr%d counter",jetbin,ithr);
      m_zeeJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
      // with reco electrons subtraction
      sprintf(counterTitle,"Z+%drecojets thr%d counter",jetbin,ithr);
      m_zeeRecoJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
      // with MC truth electrons subtraction
      sprintf(counterTitle,"Z+%dPFjets thr%d counter",jetbin,ithr);
      m_zeePFJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
      // with reco electrons subtraction
      sprintf(counterTitle,"Z+%drecoPFjets thr%d counter",jetbin,ithr);
      m_zeeRecoPFJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
    }
  }

  // W -> enu specific selections
  std::string fileCutsWenu     = theConfigDir + "WenuCuts.txt";
  std::string fileSwitchesWenu = theConfigDir + "WenuSwitches.txt";
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
  for(int ithr=0; ithr<2; ++ithr) {
    for(int jetbin=0; jetbin<6; jetbin++) {
      ConfigWenuSelections(_wenuSel, &m_wenuJetbinCounter[ithr][jetbin]);
      ConfigWenuSelections(_wenuSel, &m_wenuRecoJetbinCounter[ithr][jetbin]);
      ConfigWenuSelections(_wenuSel, &m_wenuPFJetbinCounter[ithr][jetbin]);
      ConfigWenuSelections(_wenuSel, &m_wenuRecoPFJetbinCounter[ithr][jetbin]);

      char counterTitle[200];
      // with MC truth electrons subtraction
      sprintf(counterTitle,"W+%djets thr%d counter",jetbin,ithr);
      m_wenuJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
      // with reco electrons subtraction
      sprintf(counterTitle,"W+%drecojets thr%d counter",jetbin,ithr);
      m_wenuRecoJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
      // with MC truth electrons subtraction
      sprintf(counterTitle,"W+%dPFjets thr%d counter",jetbin,ithr);
      m_wenuPFJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
      // with reco electrons subtraction
      sprintf(counterTitle,"W+%drecoPFjets thr%d counter",jetbin,ithr);
      m_wenuRecoPFJetbinCounter[ithr][jetbin].SetTitle(counterTitle);
    }
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

  m_doBestElectronStudy = m_doBTagEfficiency = false;

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
  WinAccept   = ZinAccept  = 0;
  isData_ = _commonSel->getSwitch("isData");
  if(isData_) mcevent.SetData(true);

  // btagging efficiency counters
  numBJet_reco = numNoBJet_reco = numBJet_tag = numNoBJet_tag = 0;
}

VecbosMuMuSelection::~VecbosMuMuSelection(){
  
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
void VecbosMuMuSelection::Loop() {

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
  myOutTree_Wenu->addPFelectrons();

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
  TH2F *NWPFjetMatrix15 = new TH2F("NWPFjetMatrix15","NWPFjetMatrix15",5,-0.5,4.5,5,-0.5,4.5);
  NWPFjetMatrix15->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NWPFjetMatrix15->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NWPFjetMatrix30 = new TH2F("NWPFjetMatrix30","NWPFjetMatrix30",5,-0.5,4.5,5,-0.5,4.5);
  NWPFjetMatrix30->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NWPFjetMatrix30->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NWPFjetMatrixIncl15 = new TH2F("NWPFjetMatrixIncl15","NWPFjetMatrixIncl15",5,-0.5,4.5,5,-0.5,4.5);
  NWPFjetMatrixIncl15->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NWPFjetMatrixIncl15->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NWPFjetMatrixIncl30 = new TH2F("NWPFjetMatrixIncl30","NWPFjetMatrixIncl30",5,-0.5,4.5,5,-0.5,4.5);
  NWPFjetMatrixIncl30->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NWPFjetMatrixIncl30->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  // 
  TH2F *NZPFjetMatrix15 = new TH2F("NZPFjetMatrix15","NZPFjetMatrix15",5,-0.5,4.5,5,-0.5,4.5);
  NZPFjetMatrix15->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NZPFjetMatrix15->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NZPFjetMatrix30 = new TH2F("NZPFjetMatrix30","NZPFjetMatrix30",5,-0.5,4.5,5,-0.5,4.5);
  NZPFjetMatrix30->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NZPFjetMatrix30->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NZPFjetMatrixIncl15 = new TH2F("NZPFjetMatrixIncl15","NZPFjetMatrixIncl15",5,-0.5,4.5,5,-0.5,4.5);
  NZPFjetMatrixIncl15->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NZPFjetMatrixIncl15->GetXaxis()->SetTitle("N_{pfjet}^{reco}");
  TH2F *NZPFjetMatrixIncl30 = new TH2F("NZPFjetMatrixIncl30","NZPFjetMatrixIncl30",5,-0.5,4.5,5,-0.5,4.5);
  NZPFjetMatrixIncl30->GetYaxis()->SetTitle("N_{pfjet}^{true}");
  NZPFjetMatrixIncl30->GetXaxis()->SetTitle("N_{pfjet}^{reco}");


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

  bool goodLumi = true;

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
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        goodLumi = isGoodRunLS();
      }

      if (isData_ && _commonSel->getSwitch("goodRunLS") && !goodLumi) {
        if(lastRun!= runNumber || lastLumi != lumiBlock) {
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

      reloadTriggerMask(true);
      
      // Used in the past for the soups
      float weight  = 1;

      int indexMcEleWToENu = -1;
      int foundB    = 0;
      int foundBmum = -999;
      int numB = 0;

      if ( !isData_ ) { 

        mcevent.LoadDecay(nMc,idMc,mothMc);
        mcevent.LoadMomentum(pMc,energyMc,thetaMc,phiMc);

        // check if the decay is a W->enu prompt
        indexMcEleWToENu = mcevent.indexEleWPrompt(13);
        WToENuDecay = (indexMcEleWToENu >-1 ) ? 1 : 0;

        // check if the decay is a Z->ee prompt
        int idx1Z = mcevent.indicesEleZPrompt(13).first;
        int idx2Z = mcevent.indicesEleZPrompt(13).second;
        ZToEEDecay = (idx1Z > 0 && idx2Z > 0 );

	// check if the W mc electron is within the acceptance
	ptGen1  = 1000.;
	ptGen2  = 1000.;
	etaGen1 = 1000.;
	etaGen2 = 1000.;
	meeGen  = 1000.;
	if (WToENuDecay) {
	  float ptGenEle = pMc[indexMcEleWToENu]*sin(thetaMc[indexMcEleWToENu]);
	  int goodPt = 0;
	  if (ptGenEle>20.) goodPt=1;

	  float abseta = fabs(etaMc[indexMcEleWToENu]);
	  int goodEta=0;
	  if (abseta<1.442) goodEta=1;
	  if (abseta<2.5 && abseta>1.566) goodEta=1;

	  if (goodEta==1 && goodPt==1) WinAccept=1;
	  else WinAccept=0;

	  ptGen1  = pMc[indexMcEleWToENu]*sin(thetaMc[indexMcEleWToENu]);
	  etaGen1 = etaMc[indexMcEleWToENu];
	}

	// check if the Z mc electrons are within the acceptance
	if (ZToEEDecay) {

	  float ptGenEleH, ptGenEleL;
	  float ptGenEle1 = pMc[idx1Z]*sin(thetaMc[idx1Z]);
	  float ptGenEle2 = pMc[idx2Z]*sin(thetaMc[idx2Z]);
	  if (ptGenEle1>ptGenEle2) { ptGenEleH=ptGenEle1; ptGenEleL=ptGenEle2; }
	  if (ptGenEle1<ptGenEle2) { ptGenEleH=ptGenEle2; ptGenEleL=ptGenEle1; }

	  float abseta1 = fabs(etaMc[idx1Z]);
	  float abseta2 = fabs(etaMc[idx2Z]);
	  int goodEta1=0;
	  if (abseta1<1.442) goodEta1=1;
	  if (abseta1<2.5 && abseta1>1.566) goodEta1=1;
	  int goodEta2=0;
	  if (abseta2<1.442) goodEta2=1;
	  if (abseta2<2.5 && abseta2>1.566) goodEta2=1;

	  float goodInvMass=0;
	  TLorentzVector genEleTLV1 = mcevent.p4ElectronZeePrompt().first;
	  TLorentzVector genEleTLV2 = mcevent.p4ElectronZeePrompt().second;
	  double theInvMass = (genEleTLV1+genEleTLV2).M();       
	  if (theInvMass<120. && theInvMass>60.) goodInvMass=1;

	  if (goodEta1 && goodEta2 && ptGenEleH>20 && ptGenEleL>10 && goodInvMass) ZinAccept=1;
	  else ZinAccept=0;

	  ptGen1  = ptGenEleH;
	  ptGen2  = ptGenEleL;
	  etaGen1 = etaMc[idx1Z];
	  etaGen2 = etaMc[idx2Z];
	  meeGen  = theInvMass;
	}

        // evaluate the pthat of the photon + jet event (only for photon+j events)
        photonj_pthat = photonPt();

	// check if a b-quark is present at generator level
	for (int iMc=0; iMc<nMc; iMc++) { 
	  // if ( (abs(idMc[iMc])==5) && (idMc[mothMc[iMc]]>0) && (idMc[mothMc[iMc]]!=21) ) { 
	  if ( abs(idMc[iMc])==5 ) {
	    foundB = 1; 
	    foundBmum = idMc[mothMc[iMc]];
            //            if(abs(idMc[mothMc[iMc]])!=5) numB++;
	  }
	}
      }  // MC
      

      // trigger 
      Utils anaUtils;
      bool passedHLT = hasPassedHLT();
      //      bool passedHLT = true;
//       std::cout << "passed HLT" << std::endl; 
//       for(int ii=0;ii<nEle;ii++){
//  	if(triggerMatch(etaEle[ii],phiEle[ii],0.2))
//  	  std::cout << "electron "<<ii << " has Trigger match"<<std::endl;
//       }

      int nPileUp=-1;
      // take from the event the gen jets
      // and pileup informations from simulation
      if(!isData_) {
        _theAK5GenJets = GetGenJets();
#ifdef PILEUPINFO
        nPileUp = nPU;
#endif
      }

      // take from the event the reconstructed jets
      _theAK5CaloJets     = GetCorrJets(jes_);
      _theAK5BTagCaloJets = GetBTagCorrJets();  
      if(_commonSel->getSwitch("PUSubtraction")) {
        _theAK5PFJets       = GetCorrPUPFJets(jes_);
        _theAK5BTagPFJets   = GetBTagCorrPUPFJets();
      } else {
        _theAK5PFJets       = GetCorrPFJets(jes_);
        _theAK5BTagPFJets   = GetBTagCorrPFJets();
      }

      if(_theAK5CaloJets.size() != _theAK5BTagCaloJets.size()) {
        cout << "NASTY ERROR: _theAK5CaloJets.size() != _theAK5BTagCaloJets.size()" << endl;
        return;
      }

      if(_theAK5PFJets.size() != _theAK5BTagPFJets.size()) {
        cout << "NASTY ERROR: _theAK5PFJets.size() != _theAK5BTagPFJets.size()" << endl;
        return;
      }

      // count jets removing MC truth electrons, for efficiency normalization bin-by-bin
      for(int ithr=0; ithr<2; ++ithr)
        nwjets_mc[ithr] = nwpfjets_mc[ithr] = nwgenjets_mc[ithr] = nzjets_mc[ithr] = nzpfjets_mc[ithr] = nzgenjets_mc[ithr] = 0;
      if(!isData_) {
        std::vector<TVector3> mcElectrons;

        if(mcevent.indexEleWPrompt(13) > -1) { 
          mcElectrons.push_back( mcevent.p4ElectronWenuPrompt().Vect() );

          // count calojets 
          JetCounter goodJetsMC_counter(_theAK5CaloJets);
          goodJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
          goodJetsMC_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
          goodJetsMC_counter.SetParticlesToRemove(mcElectrons);
          goodJetsMC_counter.SetPFJetId(Jet::none); // not applied for calojets
          nwjets_mc[0] = goodJetsMC_counter.numGoodJets();

          goodJetsMC_counter.reset();
          goodJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etJetAccLow"), _commonSel->getUpperCut("etaJetAcc"));
          nwjets_mc[1] = goodJetsMC_counter.numGoodJets();

          // count PF jets
          JetCounter goodPFJetsMC_counter(_theAK5PFJets);
          goodPFJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
          goodPFJetsMC_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
          goodPFJetsMC_counter.SetParticlesToRemove(mcElectrons);
          goodPFJetsMC_counter.SetPFJetId(Jet::loose);
          nwpfjets_mc[0] = goodPFJetsMC_counter.numGoodJets();

          goodPFJetsMC_counter.reset();
          goodPFJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAccLow"), _commonSel->getUpperCut("etaPFJetAcc"));
          nwpfjets_mc[1] = goodPFJetsMC_counter.numGoodJets();

          if(!isData_) {
            // count GenJets
            JetCounter genJets_counter(_theAK5GenJets);
            genJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
            genJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
            genJets_counter.SetParticlesToRemove(mcElectrons);
            genJets_counter.SetPFJetId(Jet::none); // not applied for genjets
            nwgenjets_mc[0] = genJets_counter.numGoodJets();

            genJets_counter.reset();
            genJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAccLow"), _commonSel->getUpperCut("etaPFJetAcc"));
            nwgenjets_mc[1] = genJets_counter.numGoodJets();

          }
        }

        mcElectrons.clear();
        if(mcevent.indicesEleZPrompt(13).first > -1 && mcevent.indicesEleZPrompt(13).second > -1) { 
          mcElectrons.push_back( (mcevent.p4ElectronZeePrompt().first).Vect() );
          mcElectrons.push_back( (mcevent.p4ElectronZeePrompt().second).Vect() );

          // count calojets
          JetCounter goodJetsMC_counter(_theAK5CaloJets);
          goodJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
          goodJetsMC_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
          goodJetsMC_counter.SetParticlesToRemove(mcElectrons);
          goodJetsMC_counter.SetPFJetId(Jet::none); // not applied for calojets
          nzjets_mc[0] = goodJetsMC_counter.numGoodJets();

          goodJetsMC_counter.reset();
          goodJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etJetAccLow"), _commonSel->getUpperCut("etaJetAcc"));
          nzjets_mc[1] = goodJetsMC_counter.numGoodJets();
          

          // count PF jets
          JetCounter goodPFJetsMC_counter(_theAK5PFJets);
          goodPFJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
          goodPFJetsMC_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
          goodPFJetsMC_counter.SetParticlesToRemove(mcElectrons);
          goodPFJetsMC_counter.SetPFJetId(Jet::loose);
          nzpfjets_mc[0] = goodPFJetsMC_counter.numGoodJets();

          goodPFJetsMC_counter.reset();
          goodPFJetsMC_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAccLow"), _commonSel->getUpperCut("etaPFJetAcc"));
          nzpfjets_mc[1] = goodPFJetsMC_counter.numGoodJets();

          if(!isData_) {
            // count GenJets
            JetCounter genJets_counter(_theAK5GenJets);
            genJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
            genJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
            genJets_counter.SetParticlesToRemove(mcElectrons);
            genJets_counter.SetPFJetId(Jet::none); // not applied for genjets
            nzgenjets_mc[0] = genJets_counter.numGoodJets();

            genJets_counter.reset();
            genJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAccLow"), _commonSel->getUpperCut("etaPFJetAcc"));
            nzgenjets_mc[1] = genJets_counter.numGoodJets();
          }
        }
      }

      // electrons
      std::vector<int> acceptElectrons, acceptElectronsTmp; 
      std::vector<int> idTightElectrons, idLooseElectrons; 
      std::vector<int> isolTightElectrons, isolLooseElectrons; 
      std::vector<int> convRejTightElectrons, convRejLooseElectrons; 
      std::vector<int> vertexTightElectrons, vertexLooseElectrons;

      trackerIsolation.clear();
      ecalIsolation.clear();
      hcalIsolation.clear();
      combinedIsolation.clear();
      // reconstructed electrons
      int howManyRecoEle = nMuon;

      // electrons in the acceptance - with loose cut
      for(int ii=0;ii<howManyRecoEle;ii++) {
        bool isGoodEle = true;
        TVector3 pLepton(pxMuon[ii],pyMuon[ii],pzMuon[ii]);
        float thisPt=pLepton.Pt();
        if(!_commonSel->passCut("etaElectronAcc",etaMuon[ii])) isGoodEle = false;
        if(!_commonSel->passCut("ptSlowElectronAcc",thisPt) ) isGoodEle = false;
        if (isGoodEle) acceptElectronsTmp.push_back(ii);
      }
      std::pair<int,int> bestElectronAcceptTmp = getBestGoodElePair(acceptElectronsTmp);
      int theEleAccepTmp(bestElectronAcceptTmp.first);
      TVector3 pLeptonTmp(pxMuon[theEleAccepTmp],pyMuon[theEleAccepTmp],pzMuon[theEleAccepTmp]);
      float thisPtTmp=pLeptonTmp.Pt();

      // electrons in the acceptance - with tight cut on the first ele
      if (_commonSel->passCut("ptHardElectronAcc",thisPtTmp)) {

        // electrons in the acceptance
        for(int ii=0;ii<howManyRecoEle;ii++) {             
          bool isGoodEle = true;
          TVector3 pLepton(pxMuon[ii],pyMuon[ii],pzMuon[ii]);
          float thisPt=pLepton.Pt();
          if (!_commonSel->passCut("etaElectronAcc",etaMuon[ii]))  isGoodEle = false;
	  if (!_commonSel->passCut("ptSlowElectronAcc",thisPt) ) isGoodEle = false;
          if (isGoodEle) acceptElectrons.push_back(ii);

          // assign to each electron the three isolations
          trackerIsolation.push_back( sumPt03Muon[ii] / pLepton.Pt() );
          ecalIsolation.push_back( emEt03Muon[ii] / pLepton.Pt() );
          hcalIsolation.push_back(  hadEt03Muon[ii] / pLepton.Pt() );
          // and the combined isolation (of relative isolations)
          trackerIsolRel = sumPt03Muon[ii] / pLepton.Pt();
          ecalIsolRel = emEt03Muon[ii] / pLepton.Pt();
          hcalIsolRel = hadEt03Muon[ii] / pLepton.Pt();
          float fisher = -999.;
          if(anaUtils.fiducialFlagECAL(fiducialFlagsEle[ii],isEB)) {
            fisher = reader_EB->EvaluateMVA( "Fisher" );
          } else {
            fisher = reader_EE->EvaluateMVA( "Fisher" );
          }
          combinedIsolation.push_back( (sumPt03Muon[ii]+ emEt03Muon[ii] +hadEt03Muon[ii])/pLepton.Pt() );
        }
      }
      int howManyAccEle = acceptElectrons.size();

      // for vertex studies - to check before (eleId + isolation) and after (reco + acceptance)
      std::pair<int,int> bestElectronAccept = getBestGoodElePair(acceptElectrons);
      int theEleAccep1(bestElectronAccept.first);
        
      // 'tight' (optimized) identified electrons  --- DUMMY 
      for(int ii=0;ii<howManyAccEle;ii++) {      
        idTightElectrons.push_back(acceptElectrons[ii]);            
      }
      int howManyIdTightEle = idTightElectrons.size();

      int howManyIdEle = howManyIdTightEle;

      // for vertex studies - to check before isolation and after (reco+accept+eleID)
      std::pair<int,int> bestElectronIdent = getBestGoodElePair(idTightElectrons);
      int theEleIdent1(bestElectronIdent.first);
      
      // tight optimized, tight isolated electrons
      for(int ii=0;ii<howManyIdTightEle;ii++) {      

        bool isGoodEle = true;
        int idEleIndex=idTightElectrons[ii];
        TVector3 pLepton(pxMuon[idEleIndex],pyMuon[idEleIndex],pzMuon[idEleIndex]);
        if( _commonSel->getSwitch("isolation") && (sumPt03Muon[idEleIndex]+ emEt03Muon[idEleIndex] +hadEt03Muon[idEleIndex])/pLepton.Pt()>0.15) isGoodEle=false;
        if (isGoodEle) isolTightElectrons.push_back(idEleIndex); 
      }
      int howManyIsolTightEles = isolTightElectrons.size();
    
      // tight conversion rejection electrons - DUMMY
      for(int ii=0;ii<howManyIsolTightEles;ii++) {      
        bool isGoodEle = true;
        int isolEleIndex=isolTightElectrons[ii];
        convRejTightElectrons.push_back(isolEleIndex); 
      }
      int howManyConvRejTightEles = convRejTightElectrons.size();
    
      // the highest pt isolated electron, used to define the EWK vertex
      std::pair<int,int> bestElectronConvRejTightPair = getBestGoodElePair(convRejTightElectrons);
      int theConvRejEle1(bestElectronConvRejTightPair.first);
      int theConvRejEle1_Track = trackIndexMuon[theConvRejEle1];

      int howManyDzVertexEles = 0; 
      int howManyDzVertexLooseEles = 0; 
      int howManyDzVertexTightEles = 0; 
      float EWKz = trackVzTrack[theConvRejEle1_Track];

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

          int nConvRejEles = howManyConvRejTightEles;

          // electron consistency with primary vertex
          for(int iconv=0; iconv<nConvRejEles; iconv++) {
            bool isGoodEle = true;

            int iele = convRejTightElectrons[iconv];
            int ctfTrack = trackIndexMuon[iele];

            // zele - zPV
            float dz = trackVzTrack[ctfTrack]-PVzPV[closestPV];
            if(_commonSel->getSwitch("dzVertexEles") && !_commonSel->passCut("dzVertexEles",dz)) isGoodEle = false;
            else {
              howManyDzVertexTightEles++;
            }

            // dxy significance (used if dxy not used)
            // float dxySign = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], eleTrackVxEle[iele], eleTrackVyEle[iele], eleTrackVzEle[iele], pxEle[iele], pyEle[iele], pzEle[iele]) / eleTrackDxyErrorEle[iele];
            // if(_commonSel->getSwitch("dxySigVertex") && !_commonSel->passCut("dxySigVertex",dxySign)) isGoodEle = false;
            
            // dxy
            float dxy = eleDxyPV(PVxPV[closestPV], PVyPV[closestPV], PVzPV[closestPV], trackVxTrack[ctfTrack], trackVyTrack[ctfTrack], trackVzTrack[ctfTrack], pxMuon[iele], pyMuon[iele], pzMuon[iele]);
            if(_commonSel->getSwitch("dxyVertexEles") && !_commonSel->passCut("dxyVertexEles",dxy)) isGoodEle = false;

            if(isGoodEle) {
              vertexTightElectrons.push_back(iele); 
            }
          }
        }
      }

      int howManyTightVertexEles = vertexTightElectrons.size();

      // total number of 'good' electrons
      int howManyEles = -999;
      int howManyTightEles = -999;
      howManyEles      = howManyTightVertexEles; 
      howManyTightEles = howManyTightVertexEles; 
      howManyDzVertexEles = howManyDzVertexTightEles; 

      // two highest pt good electrons - no charge requirement
      std::pair<int,int> bestElectronTightPair = getBestGoodElePair(vertexTightElectrons);
      int theEle1(bestElectronTightPair.first);
      int theEle2(bestElectronTightPair.second);
      //      getBestGoodElePairFunny(convRejTightElectrons);

      // set the event kinematics
      setKinematics(theEle1, theEle2);
      
      // decide if the 1st electron matches with the HLT object that fired the single electron trigger
      // run only if requested in order not to spend CPU...
      bool HLTmatch = true;
      //      if(_zeeSel->getSwitch("matchedHLT") || _wenuSel->getSwitch("matchedHLT")) HLTmatch = triggerMatch(etaEle[theEle1],phiEle[theEle1],0.2);

      // the basic skim is that there is >=1 electron identified
      if( writeSkim ) {
        if ( howManyIdEle > 0 ) skimFile << 1 << std::endl;
        else skimFile << 0 << std::endl;
      }

      // for optimization studies
      if(_commonSel->getSwitch("dumpIsolation")) FillIsolationVertexTree(closestPV, theEleIdent1);
      
      howManyWJets[0] = howManyWJets[1] = 0;
      howManyWPFJets[0] = howManyWPFJets[1] = 0; 
      if(theEle1_ > -1) {
        std::vector<TVector3> WElectron;
        WElectron.push_back(p4Ele1_.Vect());
        
        // count calo jets
        JetCounter goodJets_counter(_theAK5CaloJets);
        goodJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc")); 
        goodJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodJets_counter.SetParticlesToRemove(WElectron);
        goodJets_counter.SetPFJetId(Jet::none); // not applied for calojets
        (goodWJets_)[0] = goodJets_counter.getGoodJets();
        howManyWJets[0] = (goodWJets_)[0].size();
	
        goodJets_counter.reset();
        goodJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAccLow"), _commonSel->getUpperCut("etaJetAcc"));
        goodWJets_[1] = goodJets_counter.getGoodJets();
        howManyWJets[1] = goodWJets_[1].size();

        // leading calo jet pt for W
        leadingCaloJetPtW = -1.;
        emFracCaloJetW = -999.;
        for (int theJet=0; theJet<howManyWJets[1]; theJet++) {
          float theJetPt = (goodWJets_[1])[theJet].pt();
          if ( theJetPt > leadingCaloJetPtW ) { 
            leadingCaloJetPtW = theJetPt; 
            emFracCaloJetW = (goodWJets_[1])[theJet].EmFrac();
          }
        }

        // get btagging low threshold calo jets 
        JetCounter goodLowThrJets_counter(_theAK5CaloJets);
        goodLowThrJets_counter.SetThresholds(_commonSel->getLowerCut("etJetBvetoAcc"), _commonSel->getUpperCut("etaJetBvetoAcc")); 
        goodLowThrJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodLowThrJets_counter.SetParticlesToRemove(WElectron);
	goodLowThrJets_counter.SetBTagJets(_theAK5BTagCaloJets);        
        goodLowThrJets_counter.SetPFJetId(Jet::none); // not applied for calojets
	m_cleanedBvetoWJets     = goodLowThrJets_counter.getGoodJets(); 
        m_cleanedBvetoWBTagJets = goodLowThrJets_counter.getGoodBTagJets();

        // count PF jets 
        JetCounter goodPFJets_counter(_theAK5PFJets);
        goodPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc"));
        goodPFJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
        goodPFJets_counter.SetParticlesToRemove(WElectron);
        goodPFJets_counter.SetBTagJets(_theAK5BTagPFJets);
        goodPFJets_counter.SetPFJetId(Jet::loose);
        goodWPFJets_[0] = goodPFJets_counter.getGoodJets();
        howManyWPFJets[0] = goodWPFJets_[0].size();
        // count b-tagged PF jets with 5 threshold for one reference algorithm
        nTagWPFJets[0] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 5.0); // VVLoose
        nTagWPFJets[1] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 4.0); // VLoose
        nTagWPFJets[2] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 3.3); // Loose
        nTagWPFJets[3] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 2.1); // Tight
        nTagWPFJets[4] = goodPFJets_counter.numBTaggedJets(bits::trackCountingHighEffBJetTags, 1.9); // VTight

        // estimate b-tag efficiency and mistag rate (needs MC truth to find the b)
        if(!isData_ && m_doBTagEfficiency) {
          std::vector<Jet> btaggedJets = goodPFJets_counter.getGoodBTaggedJets(bits::trackCountingHighEffBJetTags, 3.3);
          countBTagJets(goodWPFJets_[0], btaggedJets);
        }

        // look for b-hadrons inside the jet cone
        if(!isData_) {
          numB = countMatchingB(goodWPFJets_[0]);
        }

        goodPFJets_counter.reset();
        goodPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAccLow"), _commonSel->getUpperCut("etaPFJetAcc"));
        goodWPFJets_[1] = goodPFJets_counter.getGoodJets();
        howManyWPFJets[1] = goodWPFJets_[1].size();        

        // leading PF jet pt
        leadingPFJetPtW = -1;
        emFracPFJetW = -999.;
        for (int theJet=0; theJet<howManyWPFJets[1]; theJet++) {
          float theJetPt = (goodWPFJets_[1])[theJet].pt();
          if ( theJetPt > leadingPFJetPtW ) { 
            leadingPFJetPtW = theJetPt; 
            emFracPFJetW = (goodWPFJets_[1])[theJet].EmFrac();
          }
        }

        // invariant mass of multi-jets
        diPFJetMassW[0] = GetDiJetMass(goodWPFJets_[0]);
        diPFJetMassW[1] = GetDiJetMass(goodWPFJets_[1]);

        // get btagging low threshold PF jets  // still to be fully implemented, fixme
        JetCounter goodLowThrPFJets_counter(_theAK5PFJets);
        goodLowThrPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetBvetoAcc"), _commonSel->getUpperCut("etaPFJetBvetoAcc")); 
        goodLowThrPFJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));  
        goodLowThrPFJets_counter.SetParticlesToRemove(WElectron);                         
        goodLowThrPFJets_counter.SetPFJetId(Jet::loose);
	// goodLowThrPFJets_counter.SetBTagJets(_theAK5BTagPFJets);        
        m_cleanedBvetoWBTagPFJets = goodLowThrPFJets_counter.getGoodBTagJets();         
	m_cleanedBvetoWPFJets     = goodLowThrPFJets_counter.getGoodJets(); 
      } 
      
      howManyZJets[0] = howManyZJets[1] = 0;
      howManyZPFJets[0] = howManyZPFJets[1] = 0;
      if(theEle1_ > -1 && theEle2_ > -1) {
        std::vector<TVector3> ZElectrons;
        ZElectrons.push_back(p4Ele1_.Vect());
        ZElectrons.push_back(p4Ele2_.Vect());
	
        // count calo jets
        JetCounter goodJets_counter(_theAK5CaloJets);
        goodJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAcc"), _commonSel->getUpperCut("etaJetAcc"));
        goodJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodJets_counter.SetParticlesToRemove(ZElectrons);
        goodJets_counter.SetPFJetId(Jet::none); // not applied for calojets
        (goodZJets_)[0] = goodJets_counter.getGoodJets();
        howManyZJets[0] = (goodZJets_)[0].size();
	
        goodJets_counter.reset();
        goodJets_counter.SetThresholds(_commonSel->getLowerCut("etJetAccLow"), _commonSel->getUpperCut("etaJetAcc"));
        (goodZJets_)[1] = goodJets_counter.getGoodJets();
        howManyZJets[1] = (goodZJets_)[1].size();        

        // leading calo jet pt for Z
        leadingCaloJetPtZ = -1.;
        emFracCaloJetZ    = -999.;
        for (int theJet=0; theJet<howManyZJets[1]; theJet++) {
          float theJetPt = (goodZJets_[1])[theJet].pt();
          if ( theJetPt > leadingCaloJetPtZ ) { 
            leadingCaloJetPtZ = theJetPt; 
            emFracCaloJetZ    = (goodZJets_[1])[theJet].EmFrac();
          }
        }

        // invariant mass of multi-jets
        diPFJetMassZ[0] = GetDiJetMass(goodZPFJets_[0]);
        diPFJetMassZ[1] = GetDiJetMass(goodZPFJets_[1]);

        // set the kinematic variables that use the selected jets
        setKinematics2();

        // get btagging low threshold calo jets 
        JetCounter goodLowThrJets_counter(_theAK5CaloJets);
        goodLowThrJets_counter.SetThresholds(_commonSel->getLowerCut("etJetBvetoAcc"), _commonSel->getUpperCut("etaJetBvetoAcc")); 
        goodLowThrJets_counter.SetDistance(_commonSel->getUpperCut("jetConeWidth"));
        goodLowThrJets_counter.SetParticlesToRemove(ZElectrons);
	goodLowThrJets_counter.SetBTagJets(_theAK5BTagCaloJets);          
        goodLowThrJets_counter.SetPFJetId(Jet::none); // not applied for calojets
	m_cleanedBvetoZJets     = goodLowThrJets_counter.getGoodJets();   
        m_cleanedBvetoZBTagJets = goodLowThrJets_counter.getGoodBTagJets(); 

        // count PF jets
        JetCounter goodPFJets_counter(_theAK5PFJets);
        goodPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAcc"), _commonSel->getUpperCut("etaPFJetAcc")); 
        goodPFJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));
        goodPFJets_counter.SetParticlesToRemove(ZElectrons);
        goodPFJets_counter.SetPFJetId(Jet::loose);
        (goodZPFJets_)[0] = goodPFJets_counter.getGoodJets();
        howManyZPFJets[0] = (goodZPFJets_)[0].size();

        goodPFJets_counter.reset();
        goodPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetAccLow"), _commonSel->getUpperCut("etaPFJetAcc"));
        (goodZPFJets_)[1] = goodPFJets_counter.getGoodJets();
        howManyZPFJets[1] = (goodZPFJets_)[1].size();

        // leading PF jet pt
        leadingPFJetPtZ = -1;
        emFracPFJetZ    = -999.;
        for (int theJet=0; theJet<howManyZPFJets[1]; theJet++) {
          float theJetPt = (goodZPFJets_[1])[theJet].pt();
          if ( theJetPt > leadingPFJetPtZ ) { 
            leadingPFJetPtZ = theJetPt; 
            emFracPFJetZ = (goodZPFJets_[1])[theJet].EmFrac();
          }
        }

        // get btagging low threshold PF jets - btag still to be fully implemented, fixme
        JetCounter goodLowThrPFJets_counter(_theAK5PFJets);
        goodLowThrPFJets_counter.SetThresholds(_commonSel->getLowerCut("etPFJetBvetoAcc"), _commonSel->getUpperCut("etaPFJetBvetoAcc")); 
        goodLowThrPFJets_counter.SetDistance(_commonSel->getUpperCut("pfJetConeWidth"));   
        goodLowThrPFJets_counter.SetParticlesToRemove(ZElectrons);                         
	// goodLowThrPFJets_counter.SetBTagJets(_theAK5BTagPFJets);          
        goodLowThrPFJets_counter.SetPFJetId(Jet::loose);
        m_cleanedBvetoZBTagPFJets = goodLowThrPFJets_counter.getGoodBTagJets();   
	m_cleanedBvetoZPFJets     = goodLowThrPFJets_counter.getGoodJets();  
      }
      
      //! event shape variables
      vector<float> MHTphiJet, MHTphiPFJet;
      float MHTphiMET;

#ifdef USECALOTOWERS
      vector<float> MHTphiJet = MHTphiJetForW(theEle1_, goodWJets_[0]);
      vector<float> MHTphiPFJet = MHTphiJetForW(theEle1_, goodWPFJets_[0]);
      float MHTphiMET = MHTphiMETForZ(theEle1_, theEle2_);
#endif

      //! B vetoes
      // remove muons from the second top->W(mu nu)b decay
      int nGoodMuons = 0;

      //! proper B tagging 
      BTagJet jetBTagEVTW, jetBTagEVTZ;
      BTagJet pfjetBTagEVTW, pfjetBTagEVTZ;
//       if(theEle1_ > -1) {
//         calcEventBVetoVariables(m_cleanedBvetoWJets,m_cleanedBvetoWBTagJets);
//         jetBTagEVTW = m_maxBTagEvt;
//         calcEventBVetoVariables(m_cleanedBvetoWPFJets,m_cleanedBvetoWBTagPFJets);
// 	pfjetBTagEVTW = m_maxBTagEvt;
//       }
//       if(theEle1_ > -1 && theEle2_ > -1) {
//         calcEventBVetoVariables(m_cleanedBvetoZJets,m_cleanedBvetoZBTagJets);
// 	jetBTagEVTZ = m_maxBTagEvt;
//         calcEventBVetoVariables(m_cleanedBvetoZPFJets,m_cleanedBvetoZBTagPFJets);
//         pfjetBTagEVTZ = m_maxBTagEvt;
//       }

      // study the best electron candidate selection for W
      if ( m_doBestElectronStudy ) {

        if ( howManyEles>0 ) {

          if( WToENuDecay ) {
            float ptMc = pMc[indexMcEleWToENu]*fabs(sin(thetaMc[indexMcEleWToENu]));

            Gene_eta[howManyWJets[0]]->Fill(etaMc[indexMcEleWToENu]);
            Gene_pt[howManyWJets[0]]->Fill(ptMc);

            int bestEleByPt = _bestPairByPt.first;
            float dr = deltaR_MCmatch(8,bestEleByPt);
            if ( electron_MCmatch_DeltaR(8,bestEleByPt,0.2) ) {
              H_etaBestByPt[howManyWJets[0]]->Fill(etaMc[indexMcEleWToENu]);
              H_ptBestByPt[howManyWJets[0]]->Fill(ptMc);
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
      bool isCommonSel[2], isZeeUpToEnd[2], isWenuUpToEnd[2], isZeeUpToEndPFJet[2], isWenuUpToEndPFJet[2] ;
      for(int ithr=0; ithr<2; ++ithr) {
        isCommonSel[ithr] = false;
        isZeeUpToEnd[ithr] = false;
        isWenuUpToEnd[ithr] = false;
        isZeeUpToEndPFJet[ithr] = false;
        isWenuUpToEndPFJet[ithr] = false;

        // for the Z veto for W->enu selection
        bool ZcandFound = foundZCandidate(theEle1_,acceptElectrons);
	// bool ZcandFound = foundZCandidate(theEle1_,idLooseElectrons);

	// to apply the majority method for charge asymmetry --- DUMMY
	bool chargeMajMethod=true;

        // common data for the selection
        SelectorData commonData;
        commonData.weight = weight;
        commonData.passedHLT = passedHLT;
        commonData.nRecoEle = howManyRecoEle;
        commonData.nAccEle = howManyAccEle;
        commonData.nIdTightEle = howManyIdTightEle;
        commonData.nIdLooseEle = 0;
        commonData.nIsolTightEle = howManyIsolTightEles;
        commonData.nIsolLooseEle = 0;
        commonData.nConvRejTightEle = howManyConvRejTightEles;
        commonData.nConvRejLooseEle = 0;
        commonData.nPV = nPV;
        commonData.nDzVertexEle = howManyDzVertexEles;
        commonData.nDxyVertexEle = howManyEles;
        commonData.matchedHLT = HLTmatch;
        commonData.nMuons = nGoodMuons;

        // Z + calo jets selection
        SelectorData dataZeeJet(commonData);
        if(m_signal == zjets) dataZeeJet.signal = ZToEEDecay;
        else if(m_signal == zother) dataZeeJet.signal = !ZToEEDecay;
        else dataZeeJet.signal = true;
        if(m_signal == zjets) dataZeeJet.inacceptance = ZinAccept;
        else dataZeeJet.inacceptance = true;
        dataZeeJet.njetsMC = nzjets_mc[ithr];
        dataZeeJet.njets = howManyZJets[ithr];
        dataZeeJet.btagEVT = jetBTagEVTZ.combinedSecondaryVertexBJetTags;
        dataZeeJet.mInv = mee_;
        dataZeeJet.mhtMET = MHTphiMET;

        CutBasedSelectorEE selectorZeeJet(dataZeeJet);
        selectorZeeJet.isMc(!isData_);
        isZeeUpToEnd[ithr]  = selectorZeeJet.outputZ(_commonSel,_zeeSel,_zeePCounter,&(m_zeeJetbinCounter[ithr]),&(m_zeeRecoJetbinCounter[ithr]));

        // W + calo jets selection
        SelectorData dataWenuJet(commonData);
        if(m_signal == wjets) dataWenuJet.signal = WToENuDecay;
        else if(m_signal == wother) dataWenuJet.signal = !WToENuDecay; 
        else dataWenuJet.signal = true;
        if(m_signal == wjets) dataWenuJet.inacceptance = WinAccept;
        else dataWenuJet.inacceptance = true;
        dataWenuJet.njetsMC = nwjets_mc[ithr];
        dataWenuJet.njets = howManyWJets[ithr];
        dataWenuJet.btagEVT = jetBTagEVTW.combinedSecondaryVertexBJetTags;
        dataWenuJet.foundAnyZ = ZcandFound;
	dataWenuJet.chargeMajorityMethod = chargeMajMethod;
        dataWenuJet.met = p3Met_.Mag();
        dataWenuJet.mt = WmT_;
        dataWenuJet.mhtJet = MHTphiJet;

        CutBasedSelectorEE selectorWenuJet(dataWenuJet);
        selectorWenuJet.isMc(!isData_);
        isWenuUpToEnd[ithr]  = selectorWenuJet.outputW(_commonSel,_wenuSel,_wenuPCounter,&(m_wenuJetbinCounter[ithr]),&(m_wenuRecoJetbinCounter[ithr]));

        // Z + PF jets selection
        SelectorData dataZeePFJet(commonData);
        if(m_signal == zjets) dataZeePFJet.signal = ZToEEDecay;
        else if(m_signal == zother) dataZeePFJet.signal = !ZToEEDecay;
        else dataZeePFJet.signal = true;
        if(m_signal == zjets) dataZeePFJet.inacceptance = ZinAccept;
        else dataZeePFJet.inacceptance = true;
        dataZeePFJet.njetsMC = nzpfjets_mc[ithr];
        dataZeePFJet.njets = howManyZPFJets[ithr];
        dataZeePFJet.btagEVT = pfjetBTagEVTZ.combinedSecondaryVertexBJetTags;
        dataZeePFJet.mInv = mee_;
        dataZeePFJet.mhtMET = MHTphiMET; // fix it with the PFjets one!

        CutBasedSelectorEE selectorZeePFJet(dataZeePFJet);
        selectorZeePFJet.isMc(!isData_);
        isZeeUpToEndPFJet[ithr]  = selectorZeePFJet.outputZ(_commonSel,_zeeSel,0,&(m_zeePFJetbinCounter[ithr]),&(m_zeeRecoPFJetbinCounter[ithr]));

        // W + PF jets selection
        SelectorData dataWenuPFJet(commonData);
        if(m_signal == wjets) dataWenuPFJet.signal = WToENuDecay;
        else if(m_signal == wother) dataWenuPFJet.signal = !WToENuDecay; 
        else dataWenuPFJet.signal = true;
        if(m_signal == wjets) dataWenuPFJet.inacceptance = WinAccept;
        else dataWenuPFJet.inacceptance = true;
        dataWenuPFJet.njetsMC = nwpfjets_mc[ithr];
        dataWenuPFJet.njets = howManyWPFJets[ithr];
        dataWenuPFJet.btagEVT = pfjetBTagEVTW.combinedSecondaryVertexBJetTags;
        dataWenuPFJet.foundAnyZ = ZcandFound;
	dataWenuPFJet.chargeMajorityMethod = chargeMajMethod;
        dataWenuPFJet.met = p3Met_.Mag();
        dataWenuPFJet.mt = WmT_;
        dataWenuPFJet.mhtJet = MHTphiJet; // fix it with the PFjets one!

        CutBasedSelectorEE selectorWenuPFJet(dataWenuPFJet);
        selectorWenuPFJet.isMc(!isData_);
        isWenuUpToEndPFJet[ithr] = selectorWenuPFJet.outputW(_commonSel,_wenuSel,0,&(m_wenuPFJetbinCounter[ithr]),&(m_wenuRecoPFJetbinCounter[ithr]));
      }

      float dummy1[2], dummy2[2], dummy3[2];

      // --------------------------------------------------------------
      if (isWenuUpToEnd[0] || isWenuUpToEnd[1] || isWenuUpToEndPFJet[0] || isWenuUpToEndPFJet[1]) {       // filling reduced trees for the W->enu study
	
        myOutTree_Wenu -> fillJetMultiplicities( howManyWJets[0], howManyWPFJets[0], -1, nwgenjets_mc[0],  
                                                 howManyWJets[1], howManyWPFJets[1], -1, nwgenjets_mc[1],
                                                 isWenuUpToEnd[0], isWenuUpToEndPFJet[0], 0, 
                                                 isWenuUpToEnd[1], isWenuUpToEndPFJet[1], 0, 
                                                 leadingCaloJetPtW, leadingPFJetPtW, emFracCaloJetW, emFracPFJetW);
        myOutTree_Wenu -> fillBTagJetMultiplicities( nTagWPFJets );
        myOutTree_Wenu -> fillElectrons(recoflag_, pt_, eta_, phi_,
                                        classification_, nbrems_, deta_, dphi_, hoe_, see_, spp_, eop_, fbrem_,
                                        trackerIso_, hcalIso_, ecalJIso_, ecalGTIso_, combinedIso_, charge_, missHits_, dist_, dcot_, lh_, e9esc_, e25esc_);
        myOutTree_Wenu -> fillKinematics( -1., WCalomT_, WTCmT_, WPFmT_, GetPt(pxMet[0],pyMet[0]), GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), p3W_.Perp(), U_.Mag(), -1., diPFJetMassW,
                                          dummy1, dummy2, dummy3, p3W_.Eta(), -1., -1., -1.);
        //        myOutTree_Wenu -> fillEventShape( MHTphiJet, MHTphiPFJet, -1.);
        myOutTree_Wenu -> fillEventInfo(nPV,rhoFastjet);
        if(!isData_) myOutTree_Wenu->fillMcTruth(WToENuDecay,nPileUp);
        myOutTree_Wenu->fillBTagEVT(jetBTagEVTW.combinedSecondaryVertexBJetTags, jetBTagEVTW.combinedSecondaryVertexMVABJetTags,
                                    jetBTagEVTW.jetBProbabilityBJetTags, jetBTagEVTW.jetProbabilityBJetTags, jetBTagEVTW.simpleSecondaryVertexBJetTags,
                                    jetBTagEVTW.softMuonBJetTags,
                                    jetBTagEVTW.trackCountingHighPurBJetTags, jetBTagEVTW.trackCountingHighEffBJetTags, foundB, foundBmum, numB);


        myOutTree_Wenu -> fillGenInfo(ptGen1, ptGen2, etaGen1, etaGen2, meeGen); 

        myOutTree_Wenu -> fillRunInfo(runNumber, lumiBlock, eventNumber);

        myOutTree_Wenu -> store();
      }


      // fill the matrix for gen electrons in acceptance
      if(WToENuDecay && WinAccept==1) {    

	int numeroPFmc0, numeroPFmc1;
	int numeroGEN0,  numeroGEN1;

	if ( nwpfjets_mc[0]<=3 ) numeroPFmc0 = nwpfjets_mc[0];
	if ( nwpfjets_mc[0]>=4 ) numeroPFmc0 = 4;

	if ( nwpfjets_mc[1]<=3 ) numeroPFmc1 = nwpfjets_mc[1];
	if ( nwpfjets_mc[1]>=4 ) numeroPFmc1 = 4;

	if ( nwgenjets_mc[0]<=3 ) numeroGEN0 = nwgenjets_mc[0];
	if ( nwgenjets_mc[0]>=4 ) numeroGEN0 = 4;

	if ( nwgenjets_mc[1]<=3 ) numeroGEN1 = nwgenjets_mc[1];
	if ( nwgenjets_mc[1]>=4 ) numeroGEN1 = 4;

	NWPFjetMatrix30->Fill(numeroPFmc0, numeroGEN0);
	NWPFjetMatrix15->Fill(numeroPFmc1, numeroGEN1);

	// for(int jettrue=0; jettrue<nwgenjets_mc[0]; jettrue++) {
	//   for(int jetreco=0; jetreco<nwpfjets_mc[0]; jetreco++) {
	//    NWPFjetMatrixIncl30->Fill(jetreco,jettrue);
	//  }
	//	}
	// for(int jettrue=0; jettrue<nwgenjets_mc[1]; jettrue++) {
	// for(int jetreco=0; jetreco<nwpfjets_mc[1]; jetreco++) {
	//   NWPFjetMatrixIncl15->Fill(jetreco,jettrue);
	// }
	// }
      }


      if (isZeeUpToEnd[0] || isZeeUpToEnd[1] || isZeeUpToEndPFJet[0] || isZeeUpToEndPFJet[1] ) {        // filling reduced trees for the Z->ee study

        vector<float> dummyW;
        for(int i=0;i<6;i++) dummyW.push_back(-999.);

        myOutTree_Zee -> fillJetMultiplicities( howManyZJets[0], howManyZPFJets[0], -1, nzgenjets_mc[0], 
                                                howManyZJets[1], howManyZPFJets[1], -1, nzgenjets_mc[1],
                                                isZeeUpToEnd[0], isZeeUpToEndPFJet[0], 0, 
                                                isZeeUpToEnd[1], isZeeUpToEndPFJet[1], 0,
                                                leadingCaloJetPtZ, leadingPFJetPtZ, emFracCaloJetZ, emFracPFJetZ);
        myOutTree_Zee -> fillElectrons(recoflag_, pt_, eta_, phi_,
                                       classification_, nbrems_, deta_, dphi_, hoe_, see_, spp_, eop_, fbrem_,
                                       trackerIso_, hcalIso_, ecalJIso_, ecalGTIso_, combinedIso_, charge_, missHits_, dist_, dcot_, lh_, e9esc_,e25esc_);
        myOutTree_Zee -> fillKinematics( mee_, -1., -1., -1., GetPt(pxMet[0],pyMet[0]), GetPt(pxTCMet[0],pyTCMet[0]), GetPt(pxPFMet[0],pyPFMet[0]), -1., -1., pt3Z_.Mag(), diPFJetMassZ, 
                                         dummy1, dummy2, dummy3, -1., p3Z_.Eta(), ZMetMT_, ZJetM_ );
        myOutTree_Zee -> fillEventInfo(nPV,rhoFastjet);
        myOutTree_Zee -> fillWTemplates( tempMETcorr_, tempMETuncorr_, WmTFromZCorr_, WmTFromZUncorr_);
        //        myOutTree_Zee -> fillEventShape( dummyW, dummyW, MHTphiMET);
        if(!isData_) myOutTree_Zee->fillMcTruth(ZToEEDecay,nPileUp);
        myOutTree_Zee->fillBTagEVT(jetBTagEVTZ.combinedSecondaryVertexBJetTags, jetBTagEVTZ.combinedSecondaryVertexMVABJetTags,
                                   jetBTagEVTZ.jetBProbabilityBJetTags, jetBTagEVTZ.jetProbabilityBJetTags, jetBTagEVTZ.simpleSecondaryVertexBJetTags,
                                   jetBTagEVTZ.softMuonBJetTags,
                                   jetBTagEVTZ.trackCountingHighPurBJetTags, jetBTagEVTZ.trackCountingHighEffBJetTags, foundB, foundBmum, numB);

        myOutTree_Zee -> fillGenInfo(ptGen1, ptGen2, etaGen1, etaGen2, meeGen);   

        myOutTree_Zee -> fillRunInfo(runNumber, lumiBlock, eventNumber);
        myOutTree_Zee -> store();

      } // 

      // fill the matrix for electrons in the acceptance
      if(ZToEEDecay && ZinAccept==1) { 
	
	int numeroPFmc0, numeroPFmc1;
	int numeroGEN0,  numeroGEN1;

	if ( nzpfjets_mc[0]<=3 ) numeroPFmc0 = nzpfjets_mc[0];
	if ( nzpfjets_mc[0]>=4 ) numeroPFmc0 = 4;

	if ( nzpfjets_mc[1]<=3 ) numeroPFmc1 = nzpfjets_mc[1];
	if ( nzpfjets_mc[1]>=4 ) numeroPFmc1 = 4;

	if ( nzgenjets_mc[0]<=3 ) numeroGEN0 = nzgenjets_mc[0];
	if ( nzgenjets_mc[0]>=4 ) numeroGEN0 = 4;

	if ( nzgenjets_mc[1]<=3 ) numeroGEN1 = nzgenjets_mc[1];
	if ( nzgenjets_mc[1]>=4 ) numeroGEN1 = 4;
	
	NZPFjetMatrix30->Fill(numeroPFmc0, numeroGEN0);
	NZPFjetMatrix15->Fill(numeroPFmc1, numeroGEN1);
	
	//for(int jettrue=0; jettrue<nzgenjets_mc[0]; jettrue++) {
	// for(int jetreco=0; jetreco<nzpfjets_mc[0]; jetreco++) {
	//  NZPFjetMatrixIncl30->Fill(jetreco,jettrue);
	//}
	//}
	//for(int jettrue=0; jettrue<nzgenjets_mc[1]; jettrue++) {
	//for(int jetreco=0; jetreco<nzpfjets_mc[1]; jetreco++) {
	//  NZPFjetMatrixIncl15->Fill(jetreco,jettrue);
	//}
	//}
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
  NWPFjetMatrix15->Write();
  NWPFjetMatrix30->Write();
  NZPFjetMatrix15->Write();
  NZPFjetMatrix30->Write();
  // NWPFjetMatrixIncl15->Write();
  // NWPFjetMatrixIncl30->Write();
  // NZPFjetMatrixIncl15->Write();
  // NZPFjetMatrixIncl30->Write();
  file.Close();
}  


void VecbosMuMuSelection::displayEfficiencies() {

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "+++ DETAILED EFFICIENCY FOR SIGNAL AS A FUNCTION OF JET MULTIPLICITY +++" << std::endl;
  std::cout << "+++ Z(ee)+jets +++" << std::endl;
  sprintf(namefile,"%sVecBosEfficiency.root",m_prefix);
  for(int ithr=0; ithr<2; ++ithr) {
    if(ithr==0) std::cout << "\t ---> HIGH ET THRESHOLD JETS <---" << std::endl;
    else std::cout << "\t ---> LOW ET THRESHOLD JETS <---" << std::endl;

    for(int ii=0; ii<6; ii++) {
      std::cout << "\tZ+" << ii << "jet" << std::endl;
      displayZeeEfficiencies(m_zeeJetbinCounter[ithr][ii]);
      const char* option = (ii==0 && ithr==0) ? "recreate" : "update";
      m_zeeJetbinCounter[ithr][ii].Save(namefile,option);
      m_zeePFJetbinCounter[ithr][ii].Save(namefile,"update");
    }
    std::cout << std::endl << std::endl << std::endl << std::endl;
    std::cout << "+++ W(enu)+jets +++" << std::endl;
    for(int ii=0; ii<6; ii++) {
      std::cout << "\tW+" << ii << "jet" << std::endl;
      displayWenuEfficiencies(m_wenuJetbinCounter[ithr][ii]);
      m_wenuJetbinCounter[ithr][ii].Save(namefile,"update");
      m_wenuPFJetbinCounter[ithr][ii].Save(namefile,"update");
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
  }

  // save into ROOT files the events as a function of jet multiplicity
  // they are different by ones stored in VecBosEfficiency.root because
  // the jet counting is done removing reco electrons, not with MC truth electrons
  sprintf(namefile,"%sVecBosCounters.root",m_prefix);
  for(int ithr=0; ithr<2; ++ithr) {
    for(int ii=0; ii<6; ii++) {
      const char* option = (ii==0 && ithr==0) ? "recreate" : "update";
      m_zeeRecoJetbinCounter[ithr][ii].Save(namefile,option);
      m_zeeRecoPFJetbinCounter[ithr][ii].Save(namefile,"update");
    }
    for(int ii=0; ii<6; ii++) {
      m_wenuRecoJetbinCounter[ithr][ii].Save(namefile,"update");
      m_wenuRecoPFJetbinCounter[ithr][ii].Save(namefile,"update");
    }
  }

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Tight ElectronID selections: " << std::endl;
  EgammaTightID_NC.displayEfficiencies();
  std::cout << endl;
  std::cout << endl;
  std::cout << "----------------------------------------"      << std::endl;
  std::cout << "Loose ElectronID selections: " << std::endl;
  EgammaLooseID_NC.displayEfficiencies();
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

  if ( m_doBTagEfficiency && !isData_ ) {
    std::cout << "<======= BTag efficiency =======>" << std::endl;
    std::cout << "efficiency for b-tagging true b-jets = " << float(numBJet_tag) / float(numBJet_reco) << endl;
    std::cout << "efficiency for b-tagging  non b-jets = " << float(numNoBJet_tag) / float(numNoBJet_reco) << endl;
    std::cout << "<======= End BTag efficiency =======>" << std::endl;
  }

}

// two highest pT 'good' electrons
std::pair<int,int> VecbosMuMuSelection::getBestGoodElePair(std::vector<int> goodElectrons) {
  int theEle1=-1;
  int theEle2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;
  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    TVector3 pEle(pxMuon[eleIndex],pyMuon[eleIndex],pzMuon[eleIndex]);
    float thisPt=pEle.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theEle2 = theEle1; theEle1 = eleIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theEle2 = eleIndex; }
  }
  return make_pair(theEle1,theEle2);
}

void VecbosMuMuSelection::setKinematics(int theEle1, int theEle2) {
  // electron candidates
  theEle1_ = theEle1;
  theEle2_ = theEle2;
  float massele = 0.000511;
  if (theEle1_ > -1) {
    TVector3 p1;
    p1.SetXYZ(pxMuon[theEle1_],pyMuon[theEle1_],pzMuon[theEle1_]);
    p4Ele1_.SetVectM(p1,massele);
    pT3Ele1_.SetXYZ(pxMuon[theEle1_],pyMuon[theEle1_],0.0);
  } else {
    p4Ele1_.SetXYZT(0.,0.,0.,0.);
  }
  if (theEle2 > -1) {
    TVector3 p2;
    p2.SetXYZ(pxMuon[theEle2],pyMuon[theEle2],pzMuon[theEle2]);
    p4Ele2_.SetVectM(p2,massele); 
    pT3Ele2_.SetXYZ(pxMuon[theEle2],pyMuon[theEle2],0.0);
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
    p3W_ = p4Ele1_.Vect() + p3Met_;
  
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
    p3W_.SetXYZ(0.,0.,0.);
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
    p3Z_.SetXYZ(p4Ele1_.Px()+p4Ele2_.Px(), p4Ele1_.Py()+p4Ele2_.Py(), p4Ele1_.Pz()+p4Ele2_.Pz());

    // correction for the mW/mZ ratio
    double zTrueMass = 91.1876;
    double wTrueMass = 80.403;	  
    double ratio = wTrueMass/zTrueMass;
    TVector3 newCorrPerp = ( (newUncorrPerp-pt3Z_)*ratio ) + pt3Z_;
    tempMETcorr_ = newCorrPerp.Mag();

    // W transverse mass from template MET
    WmTFromZUncorr_ = sqrt(2*(pT3Ele1_.Mag())*tempMETuncorr_*(1-cos(pT3Ele1_.Angle(newUncorrPerp))) );
    WmTFromZCorr_   = sqrt(2*(pT3Ele1_.Mag())*tempMETcorr_*(1-cos(pT3Ele1_.Angle(newCorrPerp))) );
    
    // Z+MET transverse mass (negletting Z mass for the moment)
    ZMetMT_ = sqrt(2 * pt3Z_.Mag() * p3Met_.Mag() * (1-cos(pt3Z_.Angle(p3Met_))) );
  } else {
    pt3Z_.SetXYZ(0.,0.,0);
    p3Z_.SetXYZ(0.,0.,0.);
    mee_ = WmTFromZUncorr_ = WmTFromZCorr_ = -1;
    ZMetMT_ = -1.;
  }

}

void VecbosMuMuSelection::setKinematics2() {

  // jets are pT-ordered
  if((goodZPFJets_[0]).size()>0) {
    TLorentzVector p4LeadingJet = (goodZPFJets_[0])[0].Get4Vector();
    TLorentzVector p4Z = p4Ele1_ + p4Ele2_;
    ZJetM_ = (p4Z+p4LeadingJet).M();
  } else ZJetM_ = -1.;

}

void VecbosMuMuSelection::setElectrons() {
  
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
      recoflag_[i] = 1;
      pt_[i] = p4[i].Pt();
      eta_[i] = p4[i].Eta();
      phi_[i] = p4[i].Phi();
      classification_[i] = 1;
      nbrems_[i] = 1;
//       if(_commonSel->getSwitch("isData")) 
//  	{
//  	  deta_[i] = deltaEtaAtVtxEle[ele] - detaCorrections(eta_[i],phi_[i]);
//  	  dphi_[i] = deltaPhiAtVtxEle[ele] - dphiCorrections(eta_[i],phi_[i]);
//  	}
//       else
// 	{
      deta_[i] = 1;
      dphi_[i] = 1;
          //        }
      hoe_[i] = 1;
      see_[i] = 1;
      spp_[i] = 1;
      eop_[i] = 1;
      fbrem_[i] = 1;
      lh_[i] = 1;
      e9esc_[i] = 1;
      e25esc_[i] = 1;
      if(_commonSel->getSwitch("PUSubtractionLeptons")) {
        trackerIso_[i] = sumPt03Muon[ele] - rhoFastjet*TMath::Pi()*0.3*0.3;
        hcalIso_[i] = hadEt03Muon[ele] - rhoFastjet*TMath::Pi()*0.3*0.3;
        ecalJIso_[i] = emEt03Muon[ele] - rhoFastjet*TMath::Pi()*0.3*0.3;
        ecalGTIso_[i] = emEt03Muon[ele] - rhoFastjet*TMath::Pi()*0.3*0.3;
        combinedIso_[i] = combinedIsolation[ele] - rhoFastjet*TMath::Pi()*0.3*0.3;
      } else {
        trackerIso_[i] = sumPt03Muon[ele];
        hcalIso_[i] = hadEt03Muon[ele];
        ecalJIso_[i] = emEt03Muon[ele];
        ecalGTIso_[i] = emEt03Muon[ele];
        combinedIso_[i] = combinedIsolation[ele];
      }
      charge_[i] =  chargeMuon[ele];
      int ctf = trackIndexEle[ele];
      missHits_[i] = 0;
      dist_[i] = 999;
      dcot_[i] = 999;
    } else { // the electron is not selected
      recoflag_[i] = -1;
      pt_[i] = -1;
      eta_[i] = -1;
      phi_[i] = -1;
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
      charge_[i] = -1;
      missHits_[i] = -1;
      dist_[i] = -1;
      dcot_[i] = -1;
    }
  }

}


// invariant mass
float VecbosMuMuSelection::getMee(int theEle1, int theEle2) {
  TLorentzVector *_p1 = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *_p2 = new TLorentzVector(0.,0.,0.,0.);
  _p1 ->SetXYZT(pxMuon[theEle1],pyMuon[theEle1],pzMuon[theEle1],energyMuon[theEle1]);
  _p2 ->SetXYZT(pxMuon[theEle2],pyMuon[theEle2],pzMuon[theEle2],energyMuon[theEle2]);
  double theInvMass = (*_p1+*_p2).M();       
  return theInvMass;
}

void VecbosMuMuSelection::FillEleIDOptimTree(int eleIndex) {

  TVector3 p3Ele(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
  int gsftrack = gsfTrackIndexEle[eleIndex];
  TVector3 p3OutEle(pxAtOuterGsfTrack[gsftrack],pyAtOuterGsfTrack[gsftrack],pzAtOuterGsfTrack[gsftrack]);

  float fbrem = fbremEle[eleIndex];

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

  float deta,dphi;

//   if(_commonSel->getSwitch("isData")) 
//     {
//       deta = deltaEtaAtVtxEle[eleIndex] - detaCorrections(etaEle[eleIndex],phiEle[eleIndex]);
//       dphi = deltaPhiAtVtxEle[eleIndex] - dphiCorrections(etaEle[eleIndex],phiEle[eleIndex]);
//     }
//   else
//     {
      deta = deltaEtaAtVtxEle[eleIndex];
      dphi = deltaPhiAtVtxEle[eleIndex];
//     }

  myOutEleIDOptimTree_Wenu->fillAll(classificationEle[eleIndex],
                                    recoFlagsEle[eleIndex],
                                    etaEle[eleIndex],
                                    GetPt(pxEle[eleIndex],pyEle[eleIndex]),
				    deta,
				    dphi,
                                    hOverEEle[eleIndex],
                                    s9s25, 
                                    see,
                                    eSeedOverPoutEle[eleIndex],
                                    fbrem);
  
  myOutEleIDOptimTree_Wenu->store();

}

void VecbosMuMuSelection::FillIsolationVertexTree(int closestPV, int iele) {
  // isolation and vertex optimization (variables are for the hardest ID electron)
  if(closestPV > -1) { // a primary vertex was found
    if(iele > -1) {
      TVector3 p3Ele(pxEle[iele],pyEle[iele],pzEle[iele]);
      int gsfTrack = gsfTrackIndexEle[iele];
      float tkIsol = dr04TkSumPtEle[iele];
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

#ifdef USECALOTOWERS
std::vector<CaloTower> VecbosMuMuSelection::getCaloTowers(float zpv, int type) {

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

void VecbosMuMuSelection::getCaloTowersForPFJets(float cptMin, float cptMax, float cchi2, float cetaMax, float cnHits, float cDxy, float cdZ, float cd3d, int iBestPV, int theEle1, int theEle2) 
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

std::vector<Jet> VecbosMuMuSelection::buildSortedSISConeCaloJets() { 

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

vector<float> VecbosMuMuSelection::MHTphiJetForW(int theEle1, std::vector<Jet>& cleanedWJets) {

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

float VecbosMuMuSelection::phiJetMETForW(std::vector<Jet>& cleanedWJets, std::vector<TLorentzVector>& tracksForPFJets)
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

float VecbosMuMuSelection::MHTphiMETForZ(int theEle1, int theEle2) {

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


void VecbosMuMuSelection::LoadCorrections(){

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
  
double VecbosMuMuSelection::ParCorr(double thePar, double theMu){
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
 
double VecbosMuMuSelection::PerpCorr(double thePerp, double theMu){
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

#endif

void VecbosMuMuSelection::ConfigCommonSelections(Selection* _selection) {

  _selection->addSwitch("isData");
  _selection->addSwitch("inacceptance");
  _selection->addSwitch("inclusiveCounting");
  _selection->addSwitch("goodRunLS");
  _selection->addSwitch("eleId");
  _selection->addSwitch("isolation");
  _selection->addSwitch("conversionRejection");
  _selection->addSwitch("asymmetricElectrons");
  _selection->addStringParameter("electronIDType");
  _selection->addStringParameter("JESUncertainty");
  _selection->addSwitch("useTCMet");
  _selection->addSwitch("usePFMet");
  _selection->addSwitch("dumpIsolation");
  _selection->addSwitch("dumpEleID");
  _selection->addSwitch("PUSubtraction");
  _selection->addSwitch("PUSubtractionLeptons");
  _selection->addCut("eventRange");
  _selection->addCut("etaElectronAcc");
  _selection->addCut("ptHardElectronAcc");
  _selection->addCut("ptSlowElectronAcc");
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("etJetAccLow");
  _selection->addCut("etJetBvetoAcc");
  _selection->addCut("etaJetBvetoAcc");
  _selection->addCut("jetConeWidth");
  _selection->addCut("etaPFJetAcc");
  _selection->addCut("etPFJetAcc");
  _selection->addCut("etPFJetAccLow");
  _selection->addCut("etPFJetBvetoAcc");
  _selection->addCut("etaPFJetBvetoAcc");
  _selection->addCut("pfJetConeWidth");
  _selection->addCut("nPrimaryVertices");
  _selection->addCut("dzVertexEles");
  _selection->addCut("dxyVertexEles");
  _selection->summary();
}

void VecbosMuMuSelection::ConfigZeeSelections(Selection* _selection, Counters *_counter) {
  
  _selection->addSwitch("trigger");

  _selection->addCut("nRecoEles");
  _selection->addCut("nAccEles");
  _selection->addCut("nIdEles");
  _selection->addCut("nIsolEles");
  _selection->addCut("nConvRejEles");
  _selection->addCut("nDzVertexEles");
  _selection->addCut("nDxyVertexEles"); 
  _selection->addSwitch("matchedHLT"); 
  _selection->addCut("nMuons");
  _selection->addCut("btagEVT");
  _selection->addCut("meeCut");
  _selection->addCut("MHTphiMET");
  _selection->summary();
  
  _counter->SetTitle("Z->EE EVENT COUNTER");
  _counter->AddVar("event");
  _counter->AddVar("mcTruth");
  _counter->AddVar("inacceptance");
  _counter->AddVar("trigger");
  _counter->AddVar("nRecoEles");
  _counter->AddVar("nAccEles");
  _counter->AddVar("nIdEles");
  _counter->AddVar("nIsolEles");
  _counter->AddVar("nConvRejEles");
  _counter->AddVar("nPrimaryVertices");
  _counter->AddVar("nDzVertexEles");
  _counter->AddVar("nDxyVertexEles");
  _counter->AddVar("matchedHLT"); 
  _counter->AddVar("nMuons");
  _counter->AddVar("btagEVT");
  _counter->AddVar("meeCut");
  _counter->AddVar("MHTphiMET");
  _counter->AddVar("fullSelection");
}

void VecbosMuMuSelection::ConfigWenuSelections(Selection* _selection, Counters *_counter) {

  _selection->addSwitch("trigger");
  
  _selection->addCut("nRecoEles");
  _selection->addCut("nAccEles");
  _selection->addCut("nIdEles");
  _selection->addCut("nIsolEles");
  _selection->addSwitch("charge");
  _selection->addCut("nConvRejEles");
  _selection->addCut("nFinalEles");
  _selection->addSwitch("ZVeto");
  _selection->addCut("nDzVertexEles");
  _selection->addCut("nDxyVertexEles"); 
  _selection->addSwitch("matchedHLT"); 
  _selection->addCut("nMuons");
  _selection->addCut("btagEVT");
  _selection->addCut("metCut");
  _selection->addCut("MHTphiJet");
  _selection->addCut("transvMassCut");
  _selection->summary();
  
  _counter->SetTitle("W->ENU EVENT COUNTER");
  _counter->AddVar("event");
  _counter->AddVar("mcTruth");
  _counter->AddVar("inacceptance");
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
  _counter->AddVar("charge");
  _counter->AddVar("matchedHLT"); 
  _counter->AddVar("nMuons");
  _counter->AddVar("btagEVT");
  _counter->AddVar("metCut");
  _counter->AddVar("transvMassCut");
  _counter->AddVar("MHTphiJet");
  _counter->AddVar("fullSelection");
}
 
void VecbosMuMuSelection::displayZeeEfficiencies(Counters _counter){ 
  _counter.Draw();
  if(_commonSel->getSwitch("isData")) {
    _counter.Draw("trigger","event");
  } else {
    _counter.Draw("mcTruth","event");
    _counter.Draw("inacceptance","mcTruth");
    _counter.Draw("trigger","inacceptance");
  }
  _counter.Draw("nRecoEles","trigger");
  _counter.Draw("nAccEles","nRecoEles");
  _counter.Draw("nIdEles","nAccEles");
  _counter.Draw("nIsolEles","nIdEles");
  _counter.Draw("nConvRejEles","nIsolEles");
  _counter.Draw("nPrimaryVertices","nConvRejEles");
  _counter.Draw("nDzVertexEles","nPrimaryVertices");
  _counter.Draw("nDxyVertexEles","nDzVertexEles");
  _counter.Draw("matchedHLT","nDxyVertexEles");
  _counter.Draw("nMuons","matchedHLT");
  _counter.Draw("btagEVT","nMuons");
  _counter.Draw("meeCut","btagEVT");
  _counter.Draw("MHTphiMET","meeCut");
  _counter.Draw("fullSelection","event");
} 

void VecbosMuMuSelection::displayWenuEfficiencies(Counters _counter){ 
  _counter.Draw();
  if(_commonSel->getSwitch("isData")) {
    _counter.Draw("trigger","event");
  } else {
    _counter.Draw("mcTruth","event");
    _counter.Draw("inacceptance","mcTruth");
    _counter.Draw("trigger","inacceptance");
  }
  _counter.Draw("nRecoEles","trigger");
  _counter.Draw("nAccEles","nRecoEles");
  _counter.Draw("nIdEles","nAccEles");
  _counter.Draw("nIsolEles","nIdEles");
  _counter.Draw("nConvRejEles","nIsolEles");
  _counter.Draw("ZVeto","nConvRejEles");
  _counter.Draw("nPrimaryVertices","ZVeto");
  _counter.Draw("nDzVertexEles","nPrimaryVertices");
  _counter.Draw("nDxyVertexEles","nDzVertexEles");
  _counter.Draw("charge","nDxyVertexEles");
  _counter.Draw("matchedHLT","charge");
  _counter.Draw("nMuons","matchedHLT");
  _counter.Draw("btagEVT","nMuons");
  _counter.Draw("metCut","btagEVT");
  _counter.Draw("transvMassCut","metCut");
  _counter.Draw("MHTphiJet","transvMassCut");
  _counter.Draw("fullSelection","event");
}

void VecbosMuMuSelection::bookBestCandidateHistos() {
  
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

bool VecbosMuMuSelection::electron_MCmatch_DeltaR(int iMc, int iEle, float deltaR) {
  
  TVector3 pMcParticle;
  pMcParticle.SetMagThetaPhi(pMc[iMc], thetaMc[iMc], phiMc[iMc]);
  
  TVector3 pRecoElectron(pxEle[iEle],pyEle[iEle],pzEle[iEle]);

  float dr = pMcParticle.DeltaR(pRecoElectron);

  if ( dr < deltaR ) return true;
  else return false;

}

float VecbosMuMuSelection::deltaR_MCmatch(int iMc, int iEle) {
  
  TVector3 pMcParticle;
  pMcParticle.SetMagThetaPhi(pMc[iMc], thetaMc[iMc], phiMc[iMc]);
  
  TVector3 pRecoElectron(pxEle[iEle],pyEle[iEle],pzEle[iEle]);

  float dr = pMcParticle.DeltaR(pRecoElectron);

  return dr;

}

bool VecbosMuMuSelection::foundZCandidate(int eleIndex, std::vector<int> accEles) {


  int howManyAccEle = accEles.size();
  for(int iAccEle=0; iAccEle<howManyAccEle; iAccEle++) {
    int iAccEleIndex = accEles[iAccEle];

    // veto a Tight x Loose combination making Z

//     bool eleId, isol, conv;
//     eleId = isol = conv = false;
//     isEleID(&EgammaLooseID_NC,iAccEleIndex,&eleId,&isol,&conv);
//     if(!eleId || !isol || !conv) continue;

    float mass = getMee(eleIndex,iAccEleIndex);
    if(_zeeSel->passCut("meeCut",mass)) return true;
  }
  
  return false;
}

int VecbosMuMuSelection::CountMuons(float etaMax, float ptMin) {

  int nMuons = 0;
  for(int imu=0; imu<nMuon; imu++) {

    // acceptance
    if(fabs(etaMuon[imu])>etaMax) continue;
    if(GetPt(pxMuon[imu],pyMuon[imu])<ptMin) continue;
    
    nMuons++;

  }

  return nMuons;

}

float VecbosMuMuSelection::GetDiJetHeaviestMass(std::vector<Jet> cleanedJets) {

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

bool VecbosMuMuSelection::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 v(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = v.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(v.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

std::vector<float> VecbosMuMuSelection::jetBTagVariables(Jet jet) {

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

void VecbosMuMuSelection::calcEventBVetoVariables(std::vector<Jet> jets, std::vector<BTagJet> btags) {

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

TLorentzVector VecbosMuMuSelection::getSC4Vector(int eleindex) {
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

void VecbosMuMuSelection::countBTagJets(std::vector<Jet> jets, std::vector<Jet> btagjets) {

  std::vector<TVector3> theBquarks;
  // find all the b's in the event MC truth
  for (int iMc=0; iMc<nMc; iMc++) { 
    if ( abs(idMc[iMc])==5 ) {
      TVector3 p3;
      p3.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);
      theBquarks.push_back(p3);
    }
  }

  //  cout << "---> found " << theBquarks.size() << " b's in this event " << endl;

  // denominator: look for MC truth b inside a reco jet
  for(unsigned int j=0; j<jets.size(); ++j) {
    bool matched = false;
    TVector3 p3Jet = jets[j].Get3Vector();
    for(unsigned int b=0; b<theBquarks.size(); b++) {
      float dR = p3Jet.DeltaR(theBquarks[b]);
      if(dR<0.3) {
        numBJet_reco++;
        matched = true;
        break;
      }
      if(!matched) numNoBJet_reco++;
    }
  }

  //  cout << "this event has " << jets.size() << " good jets, num matched = " << numBJet_reco << endl;

  // numerator: look for MC truth b inside a tagged jet 
  for(unsigned int j=0; j<btagjets.size(); ++j) {
    bool matched = false;
    TVector3 p3Jet = btagjets[j].Get3Vector();
    for(unsigned int b=0; b<theBquarks.size(); b++) {
      float dR = p3Jet.DeltaR(theBquarks[b]);
      if(dR<0.3) {
        numBJet_tag++;
        matched = true;
        break;
      }
      if(!matched) numNoBJet_tag++;
    }
  }

  //  cout << "this event has " << btagjets.size() << " good btag-jets, num matched = " << numBJet_tag << endl;

}

int VecbosMuMuSelection::countMatchingB(std::vector<Jet> jets) {
  int nb=0;

  for(unsigned int j=0; j<jets.size(); j++) {
    TVector3 p3Jet = jets[j].Get3Vector();

    // find all the b's in the event MC truth
    for (int iMc=0; iMc<nMc && iMc<100; iMc++) { 
      if ( (int(abs(idMc[iMc]))/1000) == 5 || // b-baryon
           (int(abs(idMc[iMc]))/100)%100 == 5 || // b mesons
           abs(idMc[iMc]) == 5 ) // b quark
        {
          TVector3 p3BHad;
          p3BHad.SetMagThetaPhi(pMc[iMc],thetaMc[iMc],phiMc[iMc]);

          float dR = p3Jet.DeltaR(p3BHad);
          if(dR<0.3) {
            nb++;
            break;
          }
        }
    }
  }
  return nb;
}
