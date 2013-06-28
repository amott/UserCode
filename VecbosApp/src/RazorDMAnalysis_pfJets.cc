// std includes
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
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
#include "AnalysisSelector.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "RazorDMAnalysis.hh"

//Float_t Jet_Min_Pt = 70.0;
Float_t Jet_Min_Pt = 80.0;//at least 2 jest PT>80 GeV

RazorDMAnalysis::RazorDMAnalysis(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

RazorDMAnalysis::RazorDMAnalysis(TTree *tree, string jsonFile,bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
  }
  
}

RazorDMAnalysis::~RazorDMAnalysis() {}

void RazorDMAnalysis::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorDMAnalysis::SetWeight(double weight){
  _weight=weight;
}

void RazorDMAnalysis::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;
  
  int HLT_Razor;
  int HLT_Razor_prescaled;
  int passedHLT;
  
  bool ECALTPFilterFlag;
  bool drBoundary;
  bool drDead;
  bool CSCHaloFilterFlag;
  bool trackerFailureFilterFlag;
  bool BEECALFlag; 
  bool eeBadScFilterFlag;//New filter, Cristian.(number 8 )
  bool HBHENoiseFilterResultFlag;//New filter, Cristian.(number 6 )

  // PF block
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double R[4];
  double RSQ[4];
  double MR[4];
  double MRT[4];
  double pTHem1;
  double etaHem1;
  double phiHem1;
  double pTHem2;
  double etaHem2;
  double phiHem2;
  double Jet_PT[20];
  double Jet_Eta[20];
  double Jet_Phi[20];
  double CSV[20];
  int    nBtag;
  int    nBtagTight;
  double W = _weight;
  double mst, mchi;
  int    BOX_NUM;
  Double_t Mu_Px_[2], Mu_Py_[2], Mu_Pz_[2], Mu_E_[2];
  
  double metX[4], metY[4], metCorrX[4], metCorrY[4], ht;
  double mht[3];//xyz
  
  //Cristian MC information
  
  //Int_t           nMC;
  //Float_t         pMC[2001];   //[nMC]
  //Float_t         thetaMC[2001];   //[nMC]
  //Float_t         etaMC[2001];   //[nMC]
  //Float_t         phiMC[2001];   //[nMC]
  //Float_t         energyMC[2001];   //[nMC]
  //Int_t           idMC[2001];   //[nMC]
  //Int_t           mothMC[2001];   //[nMC]
  //Int_t           statusMC[2001];   //[nMC]
  
  float mssm[3];
  //int    ss;
  int    nPV;
  int Jet_Multiplicity;
  int N_Jets;
  // gen level info
  double pT1, pT2, eta1, eta2, phi1, phi2;
  int i1, i2;
  // ttbar decay: 0 = nolep, 1 = semilep; 2 = fully lep
  int iTopDecay;
  
  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
  //outTree->Branch("ss", &ss, "ss/I");
  
  // HLT bits
  outTree->Branch("HLT_Razor", &HLT_Razor, "HLT_Razor/I"); 
  outTree->Branch("HLT_Razor_prescaled", &HLT_Razor_prescaled, "HLT_Razor_prescaled/I"); 
  outTree->Branch("passedHLT", &passedHLT, "passedHLT/I");

  //  block
  outTree->Branch("R", R, "R[4]/D");
  outTree->Branch("RSQ", RSQ, "RSQ[4]/D");
  outTree->Branch("MR", MR, "MR[4]/D");
  outTree->Branch("MRT", MRT, "MRT[4]/D");
  outTree->Branch("pTHem1", &pTHem1, "pTHem1/D");
  outTree->Branch("etaHem1", &etaHem1, "etaHem1/D");
  outTree->Branch("phiHem1", &phiHem1, "phiHem1/D");
  outTree->Branch("pTHem2", &pTHem2, "pTHem2/D");
  outTree->Branch("etaHem2", &etaHem2, "etaHem2/D");
  outTree->Branch("phiHem2", &phiHem2, "phiHem2/D");
  outTree->Branch("N_Jets", &N_Jets, "N_Jets/I");
  outTree->Branch("Jet_PT", Jet_PT, "Jet_PT[N_Jets]/D");
  outTree->Branch("Jet_Eta", Jet_Eta, "Jet_Eta[N_Jets]/D");
  outTree->Branch("Jet_Phi", Jet_Phi, "Jet_Phi[N_Jets]/D");
  outTree->Branch("CSV", CSV, "CSV[N_Jets]/D");
  outTree->Branch("nBtag", &nBtag, "nBtag/I");
  outTree->Branch("nBtagTight", &nBtagTight, "nBtagTight/I");
  outTree->Branch("W", &W, "W/D");
  outTree->Branch("mst", &mst, "mst/D");
  outTree->Branch("mchi", &mchi, "mchi/D");
  outTree->Branch("nPV", &nPV, "nPV/I");

  outTree->Branch("i1", &i1, "i1/I");
  outTree->Branch("pT1", &pT1, "pT1/D");
  outTree->Branch("eta1", &eta1, "eta1/D");
  outTree->Branch("phi1", &phi1, "phi1/D");
  outTree->Branch("i2", &i2, "i2/I");
  outTree->Branch("pT2", &pT2, "pT2/D");
  outTree->Branch("eta2", &eta2, "eta2/D");
  outTree->Branch("phi2", &phi2, "phi2/D");
  
  //outTree->Branch("N_Jets", &N_Jets, "N_Jets/I");
  
  outTree->Branch("Mu_Px", Mu_Px_,"Mu_Px_[2]/D");
  outTree->Branch("Mu_Py", Mu_Py_,"Mu_Py_[2]/D");
  outTree->Branch("Mu_Pz", Mu_Pz_,"Mu_Pz_[2]/D");
  outTree->Branch("Mu_E", Mu_E_,"Mu_E_[2]/D");
  
  outTree->Branch("iTopDecay", &iTopDecay, "iTopDecay/I");
  outTree->Branch("mssm", mssm, "mssm[3]/F");
  
  //MET Info
  outTree->Branch("metX", metX, "metX[4]/D");
  outTree->Branch("metY", metY, "metY[4]/D");
  outTree->Branch("metCorrX", metCorrX, "metCorrX[4]/D");
  outTree->Branch("metCorrY", metCorrY, "metCorrY[4]/D");
  outTree->Branch("ht", &ht, "ht/D");
  outTree->Branch("mht", mht, "mht[3]/D");//xyz->[0,1,2]
  
  //MC GEN LEVEL INFO
  
  if( _isData == 0 ){
    
    outTree->Branch("nMC", &nMc, "nMC/I");
    outTree->Branch("pMC", pMc, "pMC[101]/F");
    outTree->Branch("thetaMC", thetaMc, "thetaMC[101]/F");
    outTree->Branch("etaMC", etaMc, "etaMC[101]/F");
    outTree->Branch("phiMC", phiMc, "phiMC[101]/F");
    outTree->Branch("energyMC", energyMc, "energyMC[101]/F");
    outTree->Branch("vxMC", vxMc, "vxMC[101]/F");
    outTree->Branch("vyMC", vyMc, "vyMC[101]/F");
    outTree->Branch("vzMC", vzMc, "vzMC[101]/F");
    outTree->Branch("idMC", idMc, "idMC[101]/I");
    outTree->Branch("mothMC", mothMc, "mothMC[101]/I");
    outTree->Branch("statusMC", statusMc, "statusMC[101]/I");
    
  }
  
  double Npassed_In = 0;
  double Npassed_PV = 0;
  //Jets
  double Npassed_2Jet = 0;
  //Leptons
  double Npassed_LepVeto=0;
  //B-tag
  double Npassed_0btag=0;

  double weightII = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  std::vector<std::string> maskHLT_Razor; 
  maskHLT_Razor.push_back("HLT_RsqMR55_Rsq0p09_MR150");
  maskHLT_Razor.push_back("HLT_RsqMR60_Rsq0p09_MR150");
  maskHLT_Razor.push_back("HLT_RsqMR65_Rsq0p09_MR150");
  
  std::vector<std::string> maskHLT_Razor_prescaled; 
  maskHLT_Razor_prescaled.push_back("HLT_RsqMR40_Rsq0p04");
  //  maskHLT_Razor_prescaled.push_back("HLT_RsqMR45_Rsq0p09");
  
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
      setRequiredTriggers(maskHLT_Razor); reloadTriggerMask(true); HLT_Razor = hasPassedHLT();
      setRequiredTriggers(maskHLT_Razor_prescaled); reloadTriggerMask(true); HLT_Razor_prescaled = hasPassedHLT();
      
      ECALTPFilterFlag = (METFlags >> 0)%2;
      drBoundary = (METFlags >> 1)%2;
      drDead = (METFlags >> 2)%2;
      CSCHaloFilterFlag = (METFlags >> 3)%2;
      trackerFailureFilterFlag = (METFlags >> 4)%2;
      BEECALFlag = (METFlags >> 5)%2; 
      HBHENoiseFilterResultFlag =  (METFlags >> 6)%2;
      eeBadScFilterFlag = (METFlags >> 8)%2;
    }
    
    //Good Run selection
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
    
    Npassed_In += weightII;
    
    double m0=9999, m12=9999, mc=9999;
    
    _isSMS = false;//Chage to true when running on SMS sample
    
    if( _isData == 0  && _isSMS ){
      //double m0=9999, m12=9999, mc=9999;
      
      int ctr_s = 0;
      
      //Read and store MSSM paramenters
      mssm[0] = mssm[1] = mssm[2] = 0;
      for (std::vector<string>::iterator it = commentLHE->begin() ; it != commentLHE->end(); ++it){
	//std::cout << ctr_s << "=======STRING=====" << *it << std::endl;
	istringstream iss (*it,istringstream::in);
	string val;
	for (int n = 0; n < 5; n++){
	  iss >> val;
	  if( n == 1){
	    if( val.compare( "model" ) )break;
	  }else if( n == 2 ){
	    //std::cout << ctr_s << "=======STRING=====" << *it << std::endl;  
	    cout << val.substr( 0,val.find("_") ) << endl;
	    string aa1 = val.substr( val.find("_")+ 1 );
	    string sm0 = aa1.substr( 0, val.find("_")-1 );
	    string sm1 = aa1.substr( val.find("_") );
	    std::cout << "sm0: " << sm0 << " sm1:  " << sm1 << std::endl;
	    mssm[0] = atof(sm0.c_str());
	    mssm[1] = atof(sm1.c_str());
	  }else if( n == 3 )mssm[2] = atof( val.c_str() );
	  //cout << n << " " << val << endl;
	  
	}
	ctr_s++;
      }
      
    }
    
    
    // to integrate with others
    /*
      if(!_isData && _isSMS){
      //find the simplified model parameters for T1tttt                                                                                              
      std::vector<std::string>::const_iterator c_begin = commentLHE->begin();
      std::vector<std::string>::const_iterator c_end = commentLHE->end();
      for(std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
      size_t found = (*cit).find("T1bbbb");
      if( found != std::string::npos) {
      size_t foundLength = (*cit).size();
      found = (*cit).find("=");
      std::string smaller = (*cit).substr(found+1,foundLength);
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      
      std::istringstream iss(smaller);
      iss >> m0;
      iss.clear();
      
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      iss.str(smaller);
      iss >> m12;
      iss.clear();
      
      found = smaller.find("_");
      smaller = smaller.substr(found+1,smaller.size());
      iss.str(smaller);
      iss >> mc;
      iss.clear();
	  
      }
      }
      }
    */
    mst=m0;
    mchi=mc;
    
    //HLT and Data Filter
    passedHLT = HLT_Razor + HLT_Razor_prescaled;

    if ( _isData == true ) {
      if ( passedHLT == 0 ) continue;//Comment for getting trigger turn-ons
      if ((ECALTPFilterFlag==0) || (drBoundary==0) || (drDead==0) || (CSCHaloFilterFlag==0) || (trackerFailureFilterFlag==0) || (BEECALFlag==0) || ( HBHENoiseFilterResultFlag ==0 ) || ( eeBadScFilterFlag == 0 ) ) continue;
    }
    
    // find highest-pT PV [replace with Sagar's code]
    int iPV = passPV();
    if( iPV < 0 ) continue;//If negative no PV found to pass the cuts
    Npassed_PV += weightII;
    nPV = N_PV_EVENT;
    
    //////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    ////////////////////////// Calo Jets + JetID /////////////////////// 
    ////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    
    /*
    vector<TLorentzVector> Jet;
    vector <int> iJet;
    bool badjet = false;
    std::cout << "calo jets number " << nAK5Jet << std::endl; 
    for( int i = 0; i < nAK5Jet; i++ ) {
      TLorentzVector jet;
      double px = pxAK5Jet[i];
      double py = pyAK5Jet[i];
      double pz = pzAK5Jet[i];
      double E = sqrt( px*px + py*py + pz*pz );//massless Jet
      jet.SetPxPyPzE( px, py, pz, E );
      double pt = sqrt( px*px + py*py );
      if( fabs( jet.Eta() ) >= 3.0 ) continue;
      if( nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98 ){//What are this variables for???
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) < 2.55 && emFracAK5Jet[i] <= 0.01 ){
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) >= 2.55 && jet.Pt() > 80. && emFracAK5Jet[i] >= 1. ){
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) >= 2.55 && emFracAK5Jet[i] <= -0.9 ){
	badjet = true;
	break;
      }
      
      if ( jet.Pt() > 40. && fabs( jet.Eta() ) < 3.0 ) {
	Jet.push_back(jet);
	iJet.push_back(i);
      }
    }
    
    // Number of Jets                                                                                             
    Jet_Multiplicity[0] = Jet.size();
    
    // jet ID                                                                                                 
    if (badjet == true) continue;// If any Jet is bad (see loop before) event is rejected                    
    
    // >= 2Jets with pT> 70 GeV                                                       
    if( int( Jet.size() ) < 2 ) continue;//At least 2 Jets                                                     
    int iFirstJet = HighestPt(Jet, -99);
    int iSecondJet = HighestPt(Jet, iFirstJet);
    if( Jet[iSecondJet].Pt() < 70. ) continue;//First and second most energetic Jets Must have Pt > 70 GeV        
    Npassed_2Jet+=weightII;
    
    //count btagged jets                                                                                          
    nBtag[0] = 0;
    nBtagTight[0] = 0;
    
    for( int b = 0; b < iJet.size(); b++ ){
      if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5Jet[iJet[b]] ) ) nBtag[0]++; //Loose              
      if( pfJetPassCSVT( combinedSecondaryVertexBJetTagsAK5Jet[ iJet[b] ] ) ) nBtagTight[0]++; //Tight        
    }
    //if(nBtag>0) continue;                                                                                       
    Npassed_0btag+=weightII;//What is the meaning of this variable  
    
    /////////////////////////////                                                                                 
    ////////////HT///////////////                                                                                 
    /////////////////////////////                                                                            

    ht =  mht[0] = mht[1] = mht[2] = 0;//initialize HT and MHT                                                  
    for(  std::vector<TLorentzVector>::iterator it = Jet.begin(); it != Jet.end(); ++it){
      ht += (*it).Pt();//Compute ht, scalar sum of Pt
      mht[0] += (*it).Px();//Vectorial Sum of vec(P) of the jets                                                  
      mht[1] += (*it).Py();//Vectorial Sum of vec(P) of the jets                                                  
      mht[2] += (*it).Pz();//Vectorial Sum of vec(P) of the jets                                                  
    }

    */

    /////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //////////////////////// PF JETS + JetID ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    
    vector<TLorentzVector> pfJets;
    vector<int> i_pfJets;
    vector<double> pfJets_f_photon, pfJets_f_electron, pfJets_f_muon, pfJets_f_neutralhad, pfJets_f_chargedhad,	\
      pfJets_f_hfhad, pfJets_f_hfem;
    vector<double> pfJets_mult_photon, pfJets_mult_electron, pfJets_mult_muon, pfJets_mult_neutralhad, \
      pfJets_mult_chargedhad, pfJets_mult_hfhad, pfJets_mult_hfem;
    bool bad_pfJet = false;
    int N_pfJets = 0, pfBtags = 0;
    double pfHT, pfMHTx, pfMHTy;
    
    bool good_pfjet = false;
    //std::cout << "======== pf Jets ======== " << nAK5PFNoPUJet << std::endl;
    for(int i = 0; i < nAK5PFNoPUJet; i++){
      TLorentzVector jet;
      double px = pxAK5PFNoPUJet[i];
      double py = pyAK5PFNoPUJet[i];
      double pz = pzAK5PFNoPUJet[i];
      double E = sqrt(px*px+py*py+pz*pz);
      double scale = 1.;
      jet.SetPxPyPzE(scale*px,scale*py,scale*pz,scale*E);
      
      good_pfjet = false;
      double EU = uncorrEnergyAK5PFNoPUJet[i];
      
      if(jet.Pt() > 40.0 && fabs(jet.Eta()) < 3.0){
	
	double fHAD = (neutralHadronEnergyAK5PFNoPUJet[i]+chargedHadronEnergyAK5PFNoPUJet[i])/EU;
	
	if(fHAD > 0.99){
	  N_pfJets = 0;
	  // clean NOISY event
	  break;
	}
	
	int nConstituents = chargedHadronMultiplicityAK5PFNoPUJet[i]+neutralHadronMultiplicityAK5PFNoPUJet[i]+photonMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i]+HFHadronMultiplicityAK5PFNoPUJet[i]+HFEMMultiplicityAK5PFNoPUJet[i];
	int chargedMult = chargedHadronMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i];
	
	float photonFrac = photonEnergyAK5PFNoPUJet[i]/EU;
	float electronFrac = electronEnergyAK5PFNoPUJet[i]/EU;
	float muonFrac = muonEnergyAK5PFNoPUJet[i]/EU;
	float neutralHadFrac = neutralHadronEnergyAK5PFNoPUJet[i]/EU;
	float chargedHadFrac = chargedHadronEnergyAK5PFNoPUJet[i]/EU;
	float HFHadFrac = HFHadronEnergyAK5PFNoPUJet[i]/EU;
	float HFEMFrac = HFEMEnergyAK5PFNoPUJet[i]/EU;
	
	int photonMult = photonMultiplicityAK5PFNoPUJet[i];
	int electronMult = electronMultiplicityAK5PFNoPUJet[i];
	int muonMult = muonMultiplicityAK5PFNoPUJet[i];
	int neutralHadMult = neutralHadronMultiplicityAK5PFNoPUJet[i];
	int chargedHadMult = chargedHadronMultiplicityAK5PFNoPUJet[i];
	int HFHadMult = HFHadronMultiplicityAK5PFNoPUJet[i];
	int HFEMMult = HFEMMultiplicityAK5PFNoPUJet[i];
	
	//std::cout << "nfrac: " << neutralHadFrac << "  ph frac: " <<  photonFrac << \
	  //  "  nConst: " << nConstituents << std::endl;
	
	  //std::cout << " eta: " << fabs(jet.Eta())  << "  Char had frac: " << chargedHadFrac << \
          //"  chargedMult: " << chargedMult << "  electronFrac: " << electronFrac <<  std::endl;
	  
	if((neutralHadFrac < 0.99) && (photonFrac < 0.99) && (nConstituents > 1)) {
	  //outside of tracker acceptance, these are the only requirements
	  if (fabs(jet.Eta())>=2.4) good_pfjet = true;
	  //inside of the tracker acceptance, there are extra requirements     
	  else {
	    if ((chargedHadFrac > 0.0) && (chargedMult > 0) && (electronFrac < 0.99)) good_pfjet = true;
	  }
	}
	
	if(good_pfjet){
	  //std::cout << "Good pfjet " << i << std::endl;
	  N_pfJets++;
	  pfJets.push_back(jet);
	  i_pfJets.push_back(i);
	  pfHT += jet.Pt();
	  pfMHTx -= jet.Px();
	  pfMHTy -= jet.Py();
	  
	  //pfJets_btag.push_back(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[i]);
	  if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5PFNoPUJet[i] ) )pfBtags++;
	  pfJets_f_photon.push_back(photonFrac);
	  pfJets_f_electron.push_back(electronFrac);
	  pfJets_f_muon.push_back(muonFrac);
	  pfJets_f_neutralhad.push_back(neutralHadFrac);
	  pfJets_f_chargedhad.push_back(chargedHadFrac);
	  pfJets_f_hfhad.push_back(HFHadFrac);
	  pfJets_f_hfem.push_back(HFEMFrac);
	  
	  pfJets_mult_photon.push_back(photonMult);
	  pfJets_mult_electron.push_back(electronMult);
	  pfJets_mult_muon.push_back(muonMult);
	  pfJets_mult_neutralhad.push_back(neutralHadMult);
	  pfJets_mult_chargedhad.push_back(chargedHadMult);
	  pfJets_mult_hfhad.push_back(HFHadMult);
	  pfJets_mult_hfem.push_back(HFEMMult);
	  
	}
	else {
	  cout << "clean NOISY event" << endl;
	  N_pfJets = 0;
	  break;//Only takes out the pfJets loop! But good_jet = false
	}
      }
    }
    
    
    // jet ID                                                                 
    if (N_pfJets <= 0 )  continue;// If any Jet is bad (see loop before) event is rejected
    
    
    //////////////////////////////////////////////////////////////
    /////////////////////Create Muon Collection///////////////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    vector<int> iMuLoose;
    vector<TLorentzVector> MuLoose;
    vector<int> iMuTight;
    vector<TLorentzVector> MuTight;

    for( int i = 0; i < nMuon; i++ ) {
      TLorentzVector thisMu( pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i] );
      if( ( isTightMuon(i, true) ) && ( thisMu.Pt() > 15. ) ) {
        iMuTight.push_back(i);
        MuTight.push_back(thisMu);
	iMuLoose.push_back(i);
        MuLoose.push_back(thisMu);
      }else if( ( isLooseMuon(i, true) ) && ( thisMu.Pt() > 15. ) ) {
        iMuLoose.push_back(i);
        MuLoose.push_back(thisMu);
      }
      

    }
    

    ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    ///////////////// Create Collection  /////////////////    
    ////////////////  pfJets muon subtracted ////////////
    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    
    vector<TLorentzVector> pfJets_noMu;
    int ctr = 0;
    vector<int> i_pfJets_noMu;
    for (std::vector<TLorentzVector>::iterator pfJet_it = pfJets.begin(); pfJet_it != pfJets.end(); ++pfJet_it){
      bool isMuon = false;
      for (std::vector<TLorentzVector>::iterator Mu_it = MuLoose.begin(); Mu_it != MuLoose.end(); ++Mu_it){
	
	if ( (*Mu_it).Pt() < 15. ) continue;
	
	if( (*pfJet_it).DeltaR( *Mu_it ) <= 0.3){
	  isMuon = true;
	  break;
	}
      }
      
      if (!isMuon){
	pfJets_noMu.push_back( *pfJet_it );//If no muon found push it back into the collection  
	i_pfJets_noMu.push_back(i_pfJets[ctr]);
      }
      
      ctr++;
    }
    
    
    // Number of Jets                                                                                             
    Jet_Multiplicity = pfJets_noMu.size();
    N_Jets = pfJets_noMu.size();
    std::map<double, double> JetPTMap;
    std::map<double, double> JetEtaMap;
    std::map<double, double> JetPhiMap;
    std::vector<double> JetPT;
    for(int j = 0; j < pfJets_noMu.size(); j++){
      JetPTMap[pfJets_noMu[j].Pt()] = combinedSecondaryVertexBJetTagsAK5Jet[  i_pfJets_noMu[j]  ];
      JetEtaMap[pfJets_noMu[j].Pt()] = pfJets_noMu[j].Eta();
      JetPhiMap[pfJets_noMu[j].Pt()] = pfJets_noMu[j].Phi();
      JetPT.push_back( pfJets_noMu[j].Pt() );
    }
    
    std::sort(JetPT.begin(), JetPT.end());
    std::reverse(JetPT.begin(), JetPT.end());
    for(int j = 0; j < pfJets_noMu.size(); j++){
      Jet_PT[j] = JetPT[j];
      Jet_Eta[j] = JetEtaMap[JetPT[j]];
      Jet_Phi[j] = JetPhiMap[JetPT[j]];
      CSV[j] = JetPTMap[JetPT[j]];
    }
    
    // >= 2Jets with pT > 80 GeV                                                            
    if( int( pfJets_noMu.size() ) < 2 ) continue;//At least 2 Jets                                         
    //std::cout << "debug -1" << std::endl;
    int iFirst_pfJet = HighestPt(pfJets_noMu, -99);
    int iSecond_pfJet = HighestPt(pfJets_noMu, iFirst_pfJet);
    if( pfJets_noMu[iSecond_pfJet].Pt() < Jet_Min_Pt ) continue;//First and second most energetic Jets \
    //Must have Pt > 80 GeV  
    
    Npassed_2Jet+=weightII;
    
    //count btagged jets                                                                                          
    nBtag = 0;
    nBtagTight = 0;
    
    for( int b = 0; b < i_pfJets_noMu.size(); b++ ){
      if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5Jet[ i_pfJets_noMu[b] ] ) ) nBtag++;//Loose
      if( pfJetPassCSVT( combinedSecondaryVertexBJetTagsAK5Jet[ i_pfJets_noMu[b] ] ) ) nBtagTight++;
      //Tight  
    }
    
    
    
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ///// Muon MET Correct pfJet_noMuon ///////
    ///////////////////////////////////////////
    //////////////////////////////////////////
    
    TVector3 Muon_MET_Correction(0,0,0);

    // Boxes
    Mu_Px_[0] = Mu_Py_[0] = Mu_Pz_[0] = Mu_E_[0] = Mu_Px_[1] = Mu_Py_[1] = Mu_Pz_[1] = Mu_E_[1] =  -9999;

    // Correct MET for Presence of Loose Muons with pT>15 GeV
    int iFirstMuon = HighestPt(MuLoose, -99);
    if(iFirstMuon >= 0) {
      //Muon_MET_Correction.SetX(Muon_MET_Correction.X() + pxMuon[iFirstMuon]);
      //Muon_MET_Correction.SetY(Muon_MET_Correction.Y() + pyMuon[iFirstMuon]);

      Muon_MET_Correction.SetX(Muon_MET_Correction.X() + MuLoose[iFirstMuon].Px());
      Muon_MET_Correction.SetY(Muon_MET_Correction.Y() + MuLoose[iFirstMuon].Py());

      /*
	std::cout << "pxMuon[iFirstMuon]: " <<				\
	sqrt(pxMuon[iFirstMuon]*pxMuon[iFirstMuon]+pyMuon[iFirstMuon]*pyMuon[iFirstMuon]) << \
	" MuLoose[iFirstMuon]: " << MuLoose[iFirstMuon].Px() << std::endl;
      */
      
      /*Mu_Px_[0] = pxMuon[iFirstMuon];
      Mu_Py_[0] = pyMuon[iFirstMuon];
      Mu_Pz_[0] = pzMuon[iFirstMuon];
      Mu_E_[0] = energyMuon[iFirstMuon];
      */
      Mu_Px_[0] =  MuLoose[iFirstMuon].Px();
      Mu_Py_[0] =  MuLoose[iFirstMuon].Py();
      Mu_Pz_[0] =  MuLoose[iFirstMuon].Pz();
      Mu_E_[0] =  MuLoose[iFirstMuon].E();
    }

    int iSecondMuon = HighestPt(MuLoose, iFirstMuon);
    if(iSecondMuon >= 0) {
      //Muon_MET_Correction.SetX(Muon_MET_Correction.X() + pxMuon[iSecondMuon]); 
      //Muon_MET_Correction.SetY(Muon_MET_Correction.Y() + pyMuon[iSecondMuon]);
      
      Muon_MET_Correction.SetX(Muon_MET_Correction.X() +  MuLoose[iSecondMuon].Px()); 
      Muon_MET_Correction.SetY(Muon_MET_Correction.Y() +  MuLoose[iSecondMuon].Py());
      /*std::cout << "pxMuon[iSecondMuon]: " <<				\
	sqrt( pxMuon[iSecondMuon]*pxMuon[iSecondMuon] +  pyMuon[iSecondMuon]*pyMuon[iSecondMuon] ) << \
	" MuLoose[iSecondMuon]: " << MuLoose[iSecondMuon].Pt() << std::endl;
      */
      
      /*
      Mu_Px_[1] = pxMuon[iSecondMuon];
      Mu_Py_[1] = pyMuon[iSecondMuon];
      Mu_Pz_[1] = pzMuon[iSecondMuon];
      Mu_E_[1] = energyMuon[iSecondMuon];
      */
      
      Mu_Px_[1] = MuLoose[iSecondMuon].Px();
      Mu_Py_[1] = MuLoose[iSecondMuon].Py();
      Mu_Pz_[1] = MuLoose[iSecondMuon].Pz();
      Mu_E_[1] = MuLoose[iSecondMuon].E();
      
      
    }

    // Tight Electron ID
    vector<int> iEleTight;
    for(int i=0; i < nEle; i++) {
      TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
      if(isTrigElectron(i) && thisEle.Pt() > 15.) {
	iEleTight.push_back(i);//Fill if the electron is tight and Pt > 15 GeV
      }
    }
    
    // Tight Tau ID
    vector<int> iTauTight;
    for ( int i = 0; i < nPFTau; i++) {
      TLorentzVector thisTau( pxPFTau[i], pyPFTau[i], pzPFTau[i], energyPFTau[i] );
      if ( isTightTau(i) && thisTau.Pt() > 20. ) {
	iTauTight.push_back(i);//Fill if the tau is tight and Pt > 20. GeV
      }
    }

    // Electron VETO
    if(iEleTight.size()>0) continue;

    
    // TAU VETO
    
    //if (iTauTight.size()>0) continue;//removed after a bug was found in the tau ID

    //////////////////////////////////////////////                                                                      
    //////////Indirect Lepton Veto (taus)/////////                                                                      
    //////////////////////////////////////////////  
    bool IsoPF = false;
    for(int kk = 0; kk < nPFCand; kk++){
      TLorentzVector pfCand4V(pxPFCand[kk], pyPFCand[kk], pzPFCand[kk], energyPFCand[kk]);
      bool isMuon = false;
      for (std::vector<TLorentzVector>::iterator Mu_it = MuLoose.begin(); Mu_it != MuLoose.end(); ++Mu_it){
	
        if ( (*Mu_it).Pt() < 15. ) continue;
	
        if( (*Mu_it).DeltaR( pfCand4V ) <= 0.3){
	  //std::cout << "Mu Energy: " << (*Mu_it).E() << " pfCand En: " << pfCand4V.E() << std::endl;
	  isMuon = true;
          break;//out of the muon iterator loop
        }
      }
      
      if( ILV( kk ) < 0.15 && isMuon == false){
        IsoPF = true;
        break;//out of the PFCand loop. This means we found an isolated tau or electron 
      }
    }
    
    if(IsoPF)continue;
    
    
    Npassed_LepVeto += weightII;//Record how many of th
    
    //cout << "NBTAG IS " << nBtag << endl;
    
    // BOX NUMBER
    BOX_NUM = -99;
    if(iMuLoose.size() == 0){
      BOX_NUM = 0;
    } else if(iMuLoose.size() == 1) {
      BOX_NUM = 1;
    } else if(iMuLoose.size() >= 2) {
      BOX_NUM = 2;
    }
    
    /*
    // dummy values                                          
    pTHem1[0] = -9999.;
    etaHem1[0] = -9999.;
    phiHem1[0] = -9999.;
    pTHem2[0] = -9999.;
    etaHem2[0] = -9999.;
    phiHem2[0] = -9999.;
    pTHem1[1] = -9999.;
    etaHem1[1] = -9999.;
    phiHem1[1] = -9999.;
    pTHem2[1] = -9999.;
    etaHem2[1] = -9999.;
    phiHem2[1] = -9999.;
    */
    for(int l = 0; l < 4; l++){
      R[l] = -99999.;
      RSQ[l] = -99999.;
      MR[l] = -99999.;
      MRT[l] = -99999.;
    }
    
    /*
    // hemispheres CALO JETS
    //std::cout << pxPFMet[0] << " " << pxPFMet[1] << " " << pxPFMet[2] << " " << pxPFMet[3] << std::endl;
    vector<TLorentzVector> tmp_Jet = CombineJets(Jet);
    if( tmp_Jet.size() >= 2 ) {
      TLorentzVector Hem1 = tmp_Jet[0];
      TLorentzVector Hem2 = tmp_Jet[1];
      
      // PFMET + Correction
      TVector3 MET[4];
      for(int k = 0; k < 4; k++ ){
	MET[k] = TVector3(pxPFMet[k], pyPFMet[k], 0.);
	
	metX[k] = pxPFMet[k];
	metY[k] = pyPFMet[k];
	
	metCorrX[k] = pxPFMet[k] + Muon_MET_Correction.Px();
	metCorrY[k] = pyPFMet[k] + Muon_MET_Correction.Py();
	
	MRT[k] = CalcMTR(Hem1, Hem2, MET[k] + Muon_MET_Correction);
	//MRT = CalcMTR(Hem1, Hem2, MET);
	double variable = -999999.;
	double Rvariable = -999999.;
	variable = CalcGammaMRstar(Hem1, Hem2);
	if(variable >0) Rvariable = MRT[k]/variable;
	
	// fill the R and hem part of the output tree
	
	R[k] = Rvariable;
	RSQ[k] = Rvariable*Rvariable;
	MR[k] = variable;
      }
      
      pTHem1[0] = Hem1.Pt();
      etaHem1[0] = Hem1.Eta();
      phiHem1[0] = Hem1.Phi();
      pTHem2[0] = Hem2.Pt();
      etaHem2[0] = Hem2.Eta();
      phiHem2[0] = Hem2.Phi();
    }
    */
    
    pTHem1 = -9999.;                                                                                           
    etaHem1 = -9999.;                                                                                          
    phiHem1 = -9999.;                                                                                          
    pTHem2 = -9999.;                                                                                           
    etaHem2 = -9999.;                                                                                          
    phiHem2 = -9999.;
    
    // hemispheres PfJets Muons subtracted                                                                       
    //std::cout << pxPFMet[0] << " " << pxPFMet[1] << " " << pxPFMet[2] << " " << pxPFMet[3] << std::endl;       
    vector<TLorentzVector> tmp_pfJet = CombineJets(pfJets_noMu);
    if( tmp_pfJet.size() >= 2 ) {
      TLorentzVector Hem1 = tmp_pfJet[0];
      TLorentzVector Hem2 = tmp_pfJet[1];
      
      // PFMET + Correction                                                                                     
      TVector3 MET[4];
      for(int k = 0; k < 4; k++ ){
        MET[k] = TVector3(pxPFMet[k], pyPFMet[k], 0.);
	
        metX[k] = pxPFMet[k];
        metY[k] = pyPFMet[k];
	
	metCorrX[k] = pxPFMet[k] + Muon_MET_Correction.Px();
        metCorrY[k] = pyPFMet[k] + Muon_MET_Correction.Py();
	
        MRT[k] = CalcMTR(Hem1, Hem2, MET[k] + Muon_MET_Correction);
        //MRT = CalcMTR(Hem1, Hem2, MET);                                                                        
        double variable = -999999.;
        double Rvariable = -999999.;
        variable = CalcGammaMRstar(Hem1, Hem2);
        if(variable >0) Rvariable = MRT[k]/variable;
	
        // fill the R and hem part of the output tree                                                           
	
        R[k] = Rvariable;
        RSQ[k] = Rvariable*Rvariable;
        MR[k] = variable;
      }
      
      pTHem1 = Hem1.Pt();
      etaHem1 = Hem1.Eta();
      phiHem1 = Hem1.Phi();
      pTHem2 = Hem2.Pt();
      etaHem2 = Hem2.Eta();
      phiHem2 = Hem2.Phi();
    }
    
    
    //gen-level info
    pT1 = -999;
    eta1 = -999;
    phi1 = -999;     
    pT2 = -999;
    eta2 = -999;
    phi2 = -999;
    i1 = -99;
    i2 = -99;
    if(!_isData) {
      iTopDecay = 0;
      int iL1 = -99;
      int iL2 = -99;
      for(int i=0; i<nMc; i++) {
	// TT final state
	if(abs(idMc[mothMc[i]]) == 24) {
	  if(idMc[i] >= 11 &&
	     idMc[i] <= 18) {
	    iTopDecay ++;
	  }
	}
	// Z daughters
	if(idMc[mothMc[i]] == 23) {
	  if(iL1 <0) iL1 = i;
	  else if(iL2<0) iL2 = i;
	}
      }
      iTopDecay = iTopDecay/2;
      
      if(iL1>=0) {
	pT1 = pMc[iL1]*sin(thetaMc[iL1]);
	eta1 = etaMc[iL1];
	phi1 = phiMc[iL1];
	i1 = idMc[iL1];
      } 
      if(iL2>=0) {
	pT2 = pMc[iL2]*sin(thetaMc[iL2]);
	eta2 = etaMc[iL2];
	phi2 = phiMc[iL2];
	i2 = idMc[iL2];
      } 
    }
    
    // fill output tree
    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;
    //std::cout << "metX: " << metX << " metY: " << metY << " metCorrX: " << metCorrX << " metCorrY: " << metCorrY << " Ht: " << ht << std::endl; 
    //std::cout << "=========================nMc: " << nMc << "=======================" << std::endl;
    outTree->Fill();
  }
  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("Npassed_2Jet",  &Npassed_2Jet,  "Npassed_2Jet/D");
  effTree->Branch("Npassed_0btag",  &Npassed_0btag,  "Npassed_0btag/D");
  effTree->Branch("Npassed_LepVeto",  &Npassed_LepVeto,  "Npassed_LepVeto/D");

  effTree->Fill();
  
  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  effTree->Write();
  //  if(isSMS)FullSMSTree->Write();
  file->Close();
}

int RazorDMAnalysis::HighestPt(vector<TLorentzVector> p, int iHIGHEST) {
  
  int iH = -99;
  double highestPT = 0.;
  for(int i=0; i<p.size();i++) {
    if((p[i].Pt()>= highestPT) && (i != iHIGHEST)) {
      iH = i;
      highestPT = p[i].Pt();
    }
  }
  return iH;
}
