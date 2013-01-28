// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom3.h>

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
#include "CommonTools/include/LumiReWeightingStandAlone.h"
#include "SUSYMultiTop.hh"

// to get PWD
#include <unistd.h>

SUSYMultiTop::SUSYMultiTop(TTree *tree, double weight, int tWfile) : Vecbos(tree) {

  char tmpSTRING[516];
  getcwd(tmpSTRING,sizeof(tmpSTRING));
  std::string pwd(tmpSTRING);

  std::string fileWeight;
  if(AnalysisSelector==1 || AnalysisSelector==6){
    fileWeight=pwd+"/BTagging/8TeVNew/WeightsMu.txt";
    if(tWfile==1)fileWeight=pwd+"/BTagging/8TeVNew/WeightstWMu.txt";
    if(tWfile==-1)fileWeight=pwd+"/BTagging/8TeVNew/WeightstbarWMu.txt";
  }
  if(AnalysisSelector==3){
    fileWeight=pwd+"/BTagging/8TeVNew/WeightsEle.txt";
    if(tWfile==1)fileWeight=pwd+"/BTagging/8TeVNew/WeightstWEle.txt";
    if(tWfile==-1)fileWeight=pwd+"/BTagging/8TeVNew/WeightstbarWEle.txt";
  }

  _goodRunLS = false;
  _isData = false;
  _weight=FindWeight(fileWeight, weight);  

  _sample=FindSample(fileWeight, weight);

  //Read SF and b-tagging eff
  bool isVerbose=false;

  std::string filepTB=pwd+"/BTagging/8TeVNew/pTBinsBNew.txt";
  std::string filepTL=pwd+"/BTagging/8TeVNew/pTBinsLNew.txt";
  std::string fileEtaB=pwd+"/BTagging/8TeVNew/EtaBinsBNew.txt";
  std::string fileEtaL=pwd+"/BTagging/8TeVNew/EtaBinsL.txt";
  std::string fileSFB=pwd+"/BTagging/8TeVNew/BEff_SF_CSVM.txt";
  std::string fileSFL=pwd+"/BTagging/8TeVNew/Mistag_SF_CSVM.txt";

  std::string fileEffL=pwd+"/BTagging/8TeVNew/"+_sample+"_EffL.txt";
  std::string fileEffB=pwd+"/BTagging/8TeVNew/"+_sample+"_EffB.txt";
  std::string fileEffC=pwd+"/BTagging/8TeVNew/"+_sample+"_EffC.txt";

  std::string fileEffLFast=pwd+"/BTagging/8TeVNew/"+_sample+"_EffLFast.txt";
  std::string fileEffBFast=pwd+"/BTagging/8TeVNew/"+_sample+"_EffBFast.txt";
  std::string fileEffCFast=pwd+"/BTagging/8TeVNew/"+_sample+"_EffCFast.txt";

  std::string fileEtaBFast=pwd+"/BTagging/8TeVNew/EtaBinsBFastSim.txt";
  std::string filepTBFast=pwd+"/BTagging/8TeVNew/pTBinsBFastSim.txt";
  std::string filepTLFast=pwd+"/BTagging/8TeVNew/pTBinsLFastSim.txt";
  std::string fileEtaLFast=pwd+"/BTagging/8TeVNew/EtaBinsLFastSim.txt";
  std::string fileSFBFast=pwd+"/BTagging/8TeVNew/BEff_SF_FastSim.txt";
  std::string fileSFCFast=pwd+"/BTagging/8TeVNew/CEff_SF_FastSim.txt";
  std::string fileSFLFast=pwd+"/BTagging/8TeVNew/Mistag_SF_FastSim.txt";

  //pT and Eta bins
  pTB=ReadBins(filepTB, isVerbose);
  EtaB=ReadBins(fileEtaB, isVerbose);
  pTL=ReadBins(filepTL, isVerbose);
  EtaL=ReadBins(fileEtaL, isVerbose);

  pTBFast=ReadBins(filepTBFast, isVerbose);
  EtaBFast=ReadBins(fileEtaBFast, isVerbose);
  pTLFast=ReadBins(filepTLFast, isVerbose);
  EtaLFast=ReadBins(fileEtaLFast, isVerbose);

  //SF
  SFBv=GetSF(fileSFB, 0, false);
  SFBvUp=GetSF(fileSFB, 1, false);
  SFBvDown=GetSF(fileSFB, -1, false);;
  SFCvUp=GetSF(fileSFB, 1, true);
  SFCvDown=GetSF(fileSFB, -1, true); 
  SFLv=GetSFLight(fileSFL, 0);
  SFLvUp=GetSFLight(fileSFL, 1);
  SFLvDown=GetSFLight(fileSFL, -1);

  if(isSMS){
    SFBvFast=GetSF(fileSFBFast, 0, false);
    SFBvUpFast=GetSF(fileSFBFast, 1, false);
    SFBvDownFast=GetSF(fileSFBFast, -1, false);
    SFCvFast=GetSF(fileSFCFast, 1, true);
    SFCvUpFast=GetSF(fileSFCFast, 1, true);
    SFCvDownFast=GetSF(fileSFCFast, -1, true); 
    SFLvFast=GetSF(fileSFLFast, 0, false);
    SFLvUpFast=GetSF(fileSFLFast, 1, false);
    SFLvDownFast=GetSF(fileSFLFast, -1, false);  
  }

  //Eff
  if(!isSMS && !isEffOnly){
    EffL=GetEff(fileEffL);
    EffB=GetEff(fileEffB);
    EffC=GetEff(fileEffC);
  }else if(!isEffOnly){
    EffL=GetEff(fileEffL);
    EffB=GetEff(fileEffB);
    EffC=GetEff(fileEffC);
    EffLFast=GetEff(fileEffLFast);
    EffBFast=GetEff(fileEffBFast);
    EffCFast=GetEff(fileEffCFast);
  }

  // load tag&probe histograms
  if(!_isData) LoadTagAndProbe();

}

double SUSYMultiTop::FindWeight(string fileWeight, double weight){

  //Get vent weight
  ifstream scan;
  char content[512];
  double EvWeight,tmp, xs;
  string label;
  scan.open(fileWeight.c_str(), ifstream::in);
  if(scan.good()) {
    while (!scan.eof()) {
      scan.getline(content,512);
      scan >> label >> xs >> tmp;
      if((xs+0.00005 > weight) && (xs-0.00005 < weight))EvWeight=tmp;
    }
  }
  return EvWeight;
}

std::string SUSYMultiTop::FindSample(string fileWeight, double weight){

  //Get sample                          
  ifstream scan;
  char content[512];
  double EvWeight, xs;
  string label, tmp;
  scan.open(fileWeight.c_str(), ifstream::in);
  if(scan.good()) {
    while (!scan.eof()) {
      scan.getline(content,512);
      scan >> tmp >> xs >> EvWeight;
      if((xs+0.00005 > weight) && (xs-0.00005 < weight))label=tmp;
    }
  }
  return label;
}


SUSYMultiTop::SUSYMultiTop(TTree *tree, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    std::string goodRunGiasoneFile = "json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

}

SUSYMultiTop::~SUSYMultiTop() {}

void SUSYMultiTop::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void SUSYMultiTop::SetWeight(double weight){
  _weight=weight;
}

int SUSYMultiTop::HighestPtJet(vector<TLorentzVector> Jet, int firstJet) {

  int index=-99;
  double pT=-999999.;
  for(int i=0; i<Jet.size(); i++) {
    if(i == firstJet) continue;
    if(Jet[i].Pt()>pT) {
      pT = Jet[i].Pt();
      index = i;
    }
  }
  return index;
}

bool SUSYMultiTop::isLooseJetMva(float pt, float eta, float id) {

  bool isOk = true;

  if (pt<10) {
    if (fabs(eta)<=2.5 && id<0.0)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<0.0) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<0.0) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.2)  isOk = false;
  }

  if (pt<20 && pt>=10) {
    if (fabs(eta)<=2.5 && id<-0.4)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<-0.4) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<-0.4) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.4)   isOk = false;
  }

  if (pt<30 && pt>=20) {
    if (fabs(eta)<=2.5 && id<0.0)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<0.0) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<0.2) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.6)  isOk = false;
  }

  if (pt<50 && pt>=30) {
    if (fabs(eta)<=2.5 && id<0.0)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<0.0) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<0.6) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.2)  isOk = false;
  }

  return isOk;
}


// bool SUSYMultiTop::isLooseMuon(int iMu){

//   bool ret = false;
//   Utils anaUtils;
//   bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);

//   if(isMuGlobal){

//     int iTrack = trackIndexMuon[iMu];
//     if(numberOfValidStripTIBHitsTrack[iTrack]+
//        numberOfValidStripTIDHitsTrack[iTrack]+
//        numberOfValidStripTOBHitsTrack[iTrack]+
//        numberOfValidStripTECHitsTrack[iTrack] > 10){

//       ret = true;

//     }
//   }
//   return ret;
// }

// bool SUSYMultiTop::isGlobalMuon(int iMu){
//   bool ret = false;
//   Utils anaUtils;
//   bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
//   bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
//   bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);

//   if(isMuGlobal && isMuGlobalPrompt && isMuTracker)ret=true;

//   return ret;
// }


// bool SUSYMultiTop::isTightMuon(int iMu){
  
//   bool ret = false;
//   Utils anaUtils;
//   bool isMuGlobal = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllGlobalMuons);
//   bool isMuGlobalPrompt = anaUtils.muonIdVal(muonIdMuon[iMu], bits::GlobalMuonPromptTight);
//   bool isMuTracker = anaUtils.muonIdVal(muonIdMuon[iMu], bits::AllTrackerMuons);
  
//   if(isMuGlobal && isMuGlobalPrompt && isMuTracker){

//     double pt = sqrt(pxMuon[iMu]*pxMuon[iMu]+pyMuon[iMu]*pyMuon[iMu]);
    
//     int iTrack = trackIndexMuon[iMu];
//     if(numberOfValidStripTIBHitsTrack[iTrack]+
//        numberOfValidStripTIDHitsTrack[iTrack]+
//        numberOfValidStripTOBHitsTrack[iTrack]+
//        numberOfValidStripTECHitsTrack[iTrack] > 10){
      
//       if(numberOfValidPixelBarrelHitsTrack[iTrack] > 0 ||
// 	 numberOfValidPixelEndcapHitsTrack[iTrack] > 0){
// 	if(fabs(transvImpactParTrack[iTrack]) < 0.2){
// 	  //if(trackNormalizedChi2GlobalMuonTrack[iTrack] < 10){                                                                                           

// 	  double IECAL = emEt03Muon[iMu];
// 	  double IHCAL = hadEt03Muon[iMu];
// 	  double ITRK  = sumPt03Muon[iMu];

// 	  double CombinedIso= (IECAL+IHCAL+ITRK) - rhoFastjet * TMath::Pi() * 0.3 * 0.3;

// 	  if(CombinedIso/pt < 0.3){
// 	    ret = true;
// 	  }
// 	}
//       }
//     }
//   }

//   return ret;
// }


// bool SUSYMultiTop::is80Electron(int iEle){

//   Utils anaUtils;

//   //is an ECAL driven electron
//   bool isECALdriven = anaUtils.electronRecoType(recoFlagsEle[iEle], bits::isEcalDriven);

//   if(isECALdriven == false) return false;

//   bool isBarrel = true;
//   bool isEndcap = true;
//   int iSC = superClusterIndexEle[iEle];
//   double ETA = etaSC[iSC];
//   double ET = energySC[iSC]/cosh(etaSC[iSC]);

//   if(ET <= 20.0) return false;

//   //is the electron in the ECAL fiducial region?                                                                                                             
//   if(fabs(ETA) > 1.4442) isBarrel = false;
//   if(fabs(ETA) > 2.5 || fabs(ETA) < 1.566) isEndcap = false;

//   if(isBarrel == false && isEndcap == false) return false;

//   double trackISO = dr03TkSumPtEle[iEle];
//   double ECALISO = dr03EcalRecHitSumEtEle[iEle];
//   double HCALISO = dr03HcalTowerSumEtEle[iEle];

//   double pt = sqrt(pxEle[iEle]*pxEle[iEle]+pyEle[iEle]*pyEle[iEle]);

//   double sigietaieta = covIEtaIEtaSC[iSC];
//   double dphi = deltaPhiAtVtxEle[iEle];
//   double deta = deltaEtaAtVtxEle[iEle];
//   double HoE = hOverEEle[iEle];
//   double convDist = convDistEle[iEle];
//   double convDcot = convDcotEle[iEle];

//   if(fabs(convDist) <= 0.02 && fabs(convDcot) <= 0.02) return false;

//   int iTrack = gsfTrackIndexEle[iEle];

//   if(expInnerLayersGsfTrack[iTrack] > 0) return false;

//   trackISO /= ET;
//   ECALISO /= ET;
//   HCALISO /= ET;

//   if(isBarrel){
//     //ISO WP80                                                                                                                                               
//     if(trackISO > 0.09) return false;
//     if(ECALISO > 0.07) return false;
//     if(HCALISO > 0.10) return false;

//     //ID WP80                                                                                                                                                
//     if(fabs(dphi) > 0.06) return false;
//     if(fabs(deta) > 0.004) return false;
//     if(HoE > 0.04) return false;
//     if(sigietaieta > 0.01) return false;
//   } else {
//     //ISO WP80                                                                                                                                               
//     if(trackISO > 0.04) return false;
//     if(ECALISO > 0.05) return false;
//     if(HCALISO > 0.025) return false;

//     //ID WP80                                                                                                                                                
//     if(fabs(dphi) > 0.03) return false;
//     //if(fabs(deta) > 0.007) return false;                                                                                                                   
//     if(HoE > 0.025) return false;
//     if(sigietaieta > 0.03) return false;
//   }

//   return true;
// }


// bool SUSYMultiTop::is95Electron(int iEle){

//   Utils anaUtils;

//   //is an ECAL driven electron                                                                                                                               
//   bool isECALdriven = anaUtils.electronRecoType(recoFlagsEle[iEle], bits::isEcalDriven);
//   if(isECALdriven == false) return false;

//   bool isBarrel = true;
//   bool isEndcap = true;
//   int iSC = superClusterIndexEle[iEle];
//   double ETA = etaSC[iSC];
//   double ET = energySC[iSC]/cosh(etaSC[iSC]);

//   if(ET <= 20.0) return false;

//   //is the electron in the ECAL fiducial region?                                                                                                             
//   if(fabs(ETA) > 1.4442) isBarrel = false;
//   if(fabs(ETA) > 2.5 || fabs(ETA) < 1.566) isEndcap = false;

//   if(isBarrel == false && isEndcap == false) return false;

//   double trackISO = dr03TkSumPtEle[iEle];
//   double ECALISO = dr03EcalRecHitSumEtEle[iEle];
//   double HCALISO = dr03HcalTowerSumEtEle[iEle];

//   double pt = sqrt(pxEle[iEle]*pxEle[iEle]+pyEle[iEle]*pyEle[iEle]);

//   double sigietaieta = covIEtaIEtaSC[iSC];
//   double deta = deltaEtaAtVtxEle[iEle];
//   double HoE = hOverEEle[iEle];
//   int iTrack = gsfTrackIndexEle[iEle];

//   if(expInnerLayersGsfTrack[iTrack] > 1) return false;

//   trackISO /= ET;
//   ECALISO /= ET;
//   HCALISO /= ET;


//   if(isBarrel){
//     //ISO WP95                                                                                                                                               
//     if(trackISO > 0.15) return false;
//     if(ECALISO > 2.0) return false;
//     if(HCALISO > 0.12) return false;

//     //ID WP95                                                                                                                                                

//     if(fabs(deta) > 0.007) return false;
//     if(HoE > 0.15) return false;
//     if(sigietaieta > 0.01) return false;
//   } else {
//     //ISO WP95                                                                                                                                               
//     if(trackISO > 0.08) return false;
//     if(ECALISO > 0.06) return false;
//     if(HCALISO > 0.05) return false;

//     //ID WP95                                                                                                                                                

//     if(HoE > 0.07) return false;
//     if(sigietaieta > 0.03) return false;
//   }

//   return true;
// }


int SUSYMultiTop::JetFlavorFull(TLorentzVector myJet){
  //Returns 2 for b-jets 1 for c-jets and 0 for light jets
  //If there's a b in the cone it returns 2 (even if there's also a c in it)
  int flavor=99;
  double R;
  double DR=0.5;
  for(int i=0; i<nMc; i++){
    if(flavor==99 || flavor==1){
      R=sqrt((myJet.Eta()-etaMc[i])*(myJet.Eta()-etaMc[i])+(myJet.Phi()-phiMc[i])*(myJet.Phi()-phiMc[i]));
      if(R < DR){
	if(idMc[i]==5 || idMc[i]==-5)flavor=2;
	else if(idMc[i]==4 || idMc[i]==-4)flavor=1;
      }
    }
  }
  if(flavor==99)flavor=0;
  return flavor;
}

int SUSYMultiTop::JetFlavorEasy(TLorentzVector myJet){
  //Returns 2 for b-jets and 0 for c and light jets                                                                                             
  int flavor=99;
  double R;
  double DR=0.5;
  for(int i=0; i<nMc; i++){
    if(flavor==99){
      R=sqrt((myJet.Eta()-etaMc[i])*(myJet.Eta()-etaMc[i])+(myJet.Phi()-phiMc[i])*(myJet.Phi()-phiMc[i]));
      if(R < DR){
	if(idMc[i]==5 || idMc[i]==-5)flavor=2;
      }
    }
  }
  if(flavor==99)flavor=0;
  return flavor;
}

vector <double> SUSYMultiTop::GetSF(string _file="", int Error=0, bool isC=false){
  vector <double> SFs;
  double SF, SFErr;
  ifstream inFile;
  char content[512];
  double pT, eta;
  inFile.open(_file.c_str(), ifstream::in);
  if(inFile.good()) {
    while (!inFile.eof()) {
      inFile.getline(content,512);
      inFile >> pT >> eta >> SF >> SFErr;
      if(Error==-1 && !isC)SF=SF-SFErr;
      if(Error==1 && !isC)SF=SF+SFErr;
      if(Error==-1 && isC)SF=SF-(2.0*SFErr);
      if(Error==1 && isC)SF=SF+(2.0*SFErr);
      if(SF >0)SFs.push_back(SF);
      else SFs.push_back(0);
    }
  }
  inFile.close();
  return SFs;
}

vector <double> SUSYMultiTop::GetEff(string _file=""){
  vector <double> SFs;
  double SF;
  ifstream inFile;
  char content[512];
  inFile.open(_file.c_str(), ifstream::in);
  if(inFile.good()) {
    while (!inFile.eof()) {
      inFile.getline(content,512);
      inFile >> SF;
      SFs.push_back(SF);
    }
  }
  inFile.close();
  return SFs;
}

vector <double> SUSYMultiTop::GetSFLight(string _file="", int Error=0){
  vector <double> SFs;
  double SF, SFUp, SFDown;
  ifstream inFile;
  char content[512];
  double pT, eta;
  inFile.open(_file.c_str(), ifstream::in);
  if(inFile.good()) {
    while (!inFile.eof()) {
      inFile.getline(content,512);
      inFile >> pT >> eta >> SF >> SFUp >> SFDown;
      if(Error==-1)SFs.push_back(SFDown);
      if(Error==1)SFs.push_back(SFUp);
      if(Error==0)SFs.push_back(SF);
    }
  }
  inFile.close();
  return SFs;
}


vector <double> SUSYMultiTop::ReadBins(string _file="", bool isVerbose=false){
  vector <double> bins;
  ifstream inFile;
  char content[512];
  double bin;
  if(isVerbose)cout<<"Reading from file: "<<_file<<endl<<"Bins"<<endl;
  inFile.open(_file.c_str(), ifstream::in);
  if(inFile.good()) {
    while (!inFile.eof()) {
      inFile.getline(content,512);
      inFile >> bin;
      if(isVerbose)cout<<bin<<endl;
      bins.push_back(bin);
    }
  }
  inFile.close();
  return bins;
}

int SUSYMultiTop::NearestInt(double num){
  int n=0;
  if(num > 0.0){
    if((num+0.5) >= (int(num)+1.0)){
      n=int(num)+1;
    }else{
      n=int(num);
    }
  }else if(num < 0.0){
    if((num-0.5) >= (int(num)-1.0)){
      n=int(num)-1;
    }else{
      n=int(num);
    }
  }
  return n;
}

bool SUSYMultiTop::SwitchTagPOG(bool isBTagged, double Btag_SF, double Btag_eff){
  
  bool newBTag = isBTagged;

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw the dice
  int seed=time(NULL);
  TRandom3 * rand_ = new TRandom3(seed);
  float coin = rand_->Uniform(1.);    
  
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF)/ (1.0 - (Btag_SF/Btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }

  delete rand_;

  return newBTag;
}


bool SUSYMultiTop::SwitchTag(bool isBTagged, double Btag_SF, double Btag_eff){
  
  bool newBTag = isBTagged;

  if (Btag_SF == 1) return newBTag; //no correction needed 
  if (Btag_eff == 0) return newBTag;

  //throw the dice
  int seed=time(NULL);
  TRandom3 * rand_ = new TRandom3(seed);
  float coin = rand_->Uniform(1.);    
  
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      //fraction of jets that need to be upgraded
      float mistagPercent = (-1.0 + Btag_SF)*Btag_eff;

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin < (1.0-Btag_SF)*Btag_eff ) {newBTag = false;}

  }

  delete rand_; 

  return newBTag;
}

int SUSYMultiTop::NewTag(bool isBTagged, int pdgIdPart, double Btag_SF, double Btag_eff, double Bmistag_SF, double Bmistag_eff){

  bool newBTag = isBTagged;

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 2 ||  abs( pdgIdPart ) == 1) { 

    double bctag_eff = Btag_eff;
    //    if ( abs(pdgIdPart)== 1 )  bctag_eff = Btag_eff/5.0; // take ctag eff as one 5th of Btag eff
    newBTag = SwitchTag(isBTagged, Btag_SF, bctag_eff);

  // light quarks:
  } else if( abs( pdgIdPart )== 0 ) {

    newBTag = SwitchTag(isBTagged, Bmistag_SF, Bmistag_eff);
    
  }

  isBTagged = newBTag;
  int Tag=0;
  if(isBTagged)Tag=1;
  return Tag;
}

int SUSYMultiTop::NewTagPOG(bool isBTagged, int pdgIdPart, double Btag_SF, double Btag_eff, double Bmistag_SF, double Bmistag_eff){

  bool newBTag = isBTagged;

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 2 ||  abs( pdgIdPart ) == 1) {

    double bctag_eff = Btag_eff;
    //    if ( abs(pdgIdPart)== 1 )  bctag_eff = Btag_eff/5.0; // take ctag eff as one 5th of Btag eff                                        
    newBTag = SwitchTagPOG(isBTagged, Btag_SF, bctag_eff);

    // light quarks: 
  } else if( abs( pdgIdPart )== 0 ) {

    newBTag = SwitchTagPOG(isBTagged, Bmistag_SF, Bmistag_eff);

  }

  isBTagged = newBTag;
  int Tag=0;
  if(isBTagged)Tag=1;
  return Tag;

}

void SUSYMultiTop::SetTags(int JetFlavor, double JetpT, double JetEta, int JetTag, int& JetTagNew, int& JetTagNewBUp, int& JetTagNewBDown, int& JetTagNewLUp, int& JetTagNewLDown, bool isPOG){

  //Put SF and Effs into arrays
  double SFB[pTB.size()-1][EtaB.size()-1];
  double SFC[pTB.size()-1][EtaB.size()-1];
  double SFL[pTL.size()-1][EtaL.size()-1];
  
  double SFBUp[pTB.size()-1][EtaB.size()-1];
  double SFCUp[pTB.size()-1][EtaB.size()-1];
  double SFLUp[pTL.size()-1][EtaL.size()-1];
  
  double SFBDown[pTB.size()-1][EtaB.size()-1];
  double SFCDown[pTB.size()-1][EtaB.size()-1];
  double SFLDown[pTL.size()-1][EtaL.size()-1];
  
  double EffLd[pTL.size()-1][EtaL.size()-1];
  double EffBd[pTB.size()-1][EtaB.size()-1];
  double EffCd[pTB.size()-1][EtaB.size()-1];  
       
  int NpTL=pTL.size()-1;
  int NpTB=pTB.size()-1;

  for(int r = 0; r<(pTL.size()-1)*(EtaL.size()-1); ++r){
    SFL[(r%NpTL)][int(r/NpTL)]=SFLv[r];
    SFLUp[(r%NpTL)][int(r/NpTL)]=SFLvUp[r];
    SFLDown[(r%NpTL)][int(r/NpTL)]=SFLvDown[r];
    EffLd[(r%NpTL)][int(r/NpTL)]=EffL[r];
  }
  
  for(int t = 0; t<(pTB.size()-1)*(EtaB.size()-1); ++t){
    SFB[(t%NpTB)][int(t/NpTB)]=SFBv[t];
    SFC[(t%NpTB)][int(t/NpTB)]=SFBv[t];
    SFBUp[(t%NpTB)][int(t/NpTB)]=SFBvUp[t];
    SFCUp[(t%NpTB)][int(t/NpTB)]=SFCvUp[t];
    SFBDown[(t%NpTB)][int(t/NpTB)]=SFBvDown[t];
    SFCDown[(t%NpTB)][int(t/NpTB)]=SFCvDown[t];
    EffBd[(t%NpTB)][int(t/NpTB)]=EffB[t];
    EffCd[(t%NpTB)][int(t/NpTB)]=EffC[t];
  }
  
  int npT, nEta;
  if(JetFlavor==0){
    for(int p = 0; p<(pTL.size()-1); ++p){
      if(JetpT < pTL[p+1] && JetpT >= pTL[p])npT=p;
    }
    if(JetpT < pTL[0] || JetpT > pTL[(pTL.size()-1)])npT=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaL.size()-1); ++q){
      if(fabs(JetEta) < EtaL[q+1] && fabs(JetEta) >= EtaL[q])nEta=q;
    }
    if(fabs(JetEta) < EtaL[0] || fabs(JetEta) > EtaL[(EtaL.size()-1)])nEta=-99; //For Jets outside the measured SF range
    if(npT!=-99 && nEta!=-99){
      bool isBTagged=true;
      if(JetTag==0)isBTagged=false;
      if(!isPOG){
	JetTagNew=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBUp=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBDown=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLUp=NewTag(isBTagged, 0, 1, 1, SFLUp[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLDown=NewTag(isBTagged, 0, 1, 1, SFLDown[npT][nEta], EffLd[npT][nEta]);
      }else{
	JetTagNew=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBUp=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBDown=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLUp=NewTagPOG(isBTagged, 0, 1, 1, SFLUp[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLDown=NewTagPOG(isBTagged, 0, 1, 1, SFLDown[npT][nEta], EffLd[npT][nEta]);
      }
    }
  }else if(JetFlavor==2){
    for(int p = 0; p<(pTB.size()-1); ++p){
      if(JetpT < pTB[p+1] && JetpT >= pTB[p])npT=p;
    }
    if(JetpT < pTB[0] || JetpT > pTB[(pTB.size()-1)])npT=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaB.size()-1); ++q){
      if(fabs(JetEta) < EtaB[q+1] && fabs(JetEta) >= EtaB[q])nEta=q;
    }
    if(fabs(JetEta) < EtaB[0] || fabs(JetEta) > EtaB[(EtaB.size()-1)])nEta=-99; //For Jets outside the measured SF range
    if(npT!=-99 && nEta!=-99){
      bool isBTagged=true;
      if(JetTag==0)isBTagged=false;
      if(!isPOG){
	JetTagNew=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTag(isBTagged, 2, SFBUp[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTag(isBTagged, 2, SFBDown[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
      }else{
	JetTagNew=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTagPOG(isBTagged, 2, SFBUp[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTagPOG(isBTagged, 2, SFBDown[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
      }
    }
  }else if(JetFlavor==1){
    for(int p = 0; p<(pTB.size()-1); ++p){
      if(JetpT < pTB[p+1] && JetpT >= pTB[p])npT=p;
    }
    if(JetpT < pTB[0] || JetpT > pTB[(pTB.size()-1)])npT=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaB.size()-1); ++q){
      if(fabs(JetEta) < EtaB[q+1] && fabs(JetEta) >= EtaB[q])nEta=q;
    }
    if(fabs(JetEta) < EtaB[0] || fabs(JetEta) > EtaB[(EtaB.size()-1)])nEta=-99; //For Jets outside the measured SF range          
    if(npT!=-99 && nEta!=-99){
      bool isBTagged=true;
      if(JetTag==0)isBTagged=false;
      if(!isPOG){
	JetTagNew=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTag(isBTagged, 1, SFCUp[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTag(isBTagged, 1, SFCDown[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
      }else{
	JetTagNew=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTagPOG(isBTagged, 1, SFCUp[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTagPOG(isBTagged, 1, SFCDown[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
      }
    }
  }

  return;
  
}


vector<int> SUSYMultiTop::SetTagsSMS(int JetFlavor, double JetpT, double JetEta, int JetTag, bool isPOG){

  vector <int> Tags;
  int JetTagNew;
  int JetTagNewBUp;
  int JetTagNewBDown;
  int JetTagNewLUp;
  int JetTagNewLDown;
  int JetTagNewBUpFast;
  int JetTagNewBDownFast;
  int JetTagNewCUpFast;
  int JetTagNewCDownFast;
  int JetTagNewLUpFast;
  int JetTagNewLDownFast;
  int JetTagNewFast;

  //Put SF and Effs into arrays
  //Full Data
  double SFB[pTB.size()-1][EtaB.size()-1];
  double SFC[pTB.size()-1][EtaB.size()-1];
  double SFL[pTL.size()-1][EtaL.size()-1];
  
  double SFBUp[pTB.size()-1][EtaB.size()-1];
  double SFCUp[pTB.size()-1][EtaB.size()-1];
  double SFLUp[pTL.size()-1][EtaL.size()-1];
  
  double SFBDown[pTB.size()-1][EtaB.size()-1];
  double SFCDown[pTB.size()-1][EtaB.size()-1];
  double SFLDown[pTL.size()-1][EtaL.size()-1];
  
  double EffLd[pTL.size()-1][EtaL.size()-1];
  double EffBd[pTB.size()-1][EtaB.size()-1];
  double EffCd[pTB.size()-1][EtaB.size()-1];
       
  int NpTL=pTL.size()-1;
  int NpTB=pTB.size()-1;
     
  for(int r = 0; r<(pTL.size()-1)*(EtaL.size()-1); ++r){
    SFL[(r%NpTL)][int(r/NpTL)]=SFLv[r];
    SFLUp[(r%NpTL)][int(r/NpTL)]=SFLvUp[r];
    SFLDown[(r%NpTL)][int(r/NpTL)]=SFLvDown[r];
    EffLd[(r%NpTL)][int(r/NpTL)]=EffL[r];
  }
  
  for(int t = 0; t<(pTB.size()-1)*(EtaB.size()-1); ++t){
    SFB[(t%NpTB)][int(t/NpTB)]=SFBv[t];
    SFC[(t%NpTB)][int(t/NpTB)]=SFBv[t];
    SFBUp[(t%NpTB)][int(t/NpTB)]=SFBvUp[t];
    SFCUp[(t%NpTB)][int(t/NpTB)]=SFCvUp[t];
    SFBDown[(t%NpTB)][int(t/NpTB)]=SFBvDown[t];
    SFCDown[(t%NpTB)][int(t/NpTB)]=SFCvDown[t];
    EffBd[(t%NpTB)][int(t/NpTB)]=EffB[t];
    EffCd[(t%NpTB)][int(t/NpTB)]=EffC[t];
  }

  //Fast Full
  double SFBFast[pTBFast.size()-1][EtaBFast.size()-1];
  double SFCFast[pTBFast.size()-1][EtaBFast.size()-1];
  double SFLFast[pTLFast.size()-1][EtaLFast.size()-1];
  
  double SFBUpFast[pTBFast.size()-1][EtaBFast.size()-1];
  double SFCUpFast[pTBFast.size()-1][EtaBFast.size()-1];
  double SFLUpFast[pTLFast.size()-1][EtaLFast.size()-1];
  
  double SFBDownFast[pTBFast.size()-1][EtaBFast.size()-1];
  double SFCDownFast[pTBFast.size()-1][EtaBFast.size()-1];
  double SFLDownFast[pTLFast.size()-1][EtaLFast.size()-1];
  
  double EffLdFast[pTLFast.size()-1][EtaLFast.size()-1];
  double EffBdFast[pTBFast.size()-1][EtaBFast.size()-1];
  double EffCdFast[pTBFast.size()-1][EtaBFast.size()-1];
  
  int NpTBFast=pTBFast.size()-1;
  int NpTLFast=pTLFast.size()-1;
  
  for(int r = 0; r<(pTLFast.size()-1)*(EtaLFast.size()-1); ++r){
    SFLFast[(r%NpTLFast)][int(r/NpTLFast)]=SFLvFast[r];
    SFLUpFast[(r%NpTLFast)][int(r/NpTLFast)]=SFLvUpFast[r];
    SFLDownFast[(r%NpTLFast)][int(r/NpTLFast)]=SFLvDownFast[r];
    EffLdFast[(r%NpTLFast)][int(r/NpTLFast)]=EffLFast[r];
  }
  
  for(int t = 0; t<(pTBFast.size()-1)*(EtaBFast.size()-1); ++t){
    SFBFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvFast[t];
    SFCFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvFast[t];
    SFBUpFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvUpFast[t];
    SFCUpFast[(t%NpTBFast)][int(t/NpTBFast)]=SFCvUpFast[t];
    SFBDownFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvDownFast[t];
    SFCDownFast[(t%NpTBFast)][int(t/NpTBFast)]=SFCvDownFast[t];
    EffBdFast[(t%NpTBFast)][int(t/NpTBFast)]=EffBFast[t];
    EffCdFast[(t%NpTBFast)][int(t/NpTBFast)]=EffCFast[t];
  }
  
  int npT, nEta;
  int npTFast, nEtaFast;

  if(JetFlavor==0){
    for(int p = 0; p<(pTL.size()-1); ++p){
      if(JetpT < pTL[p+1] && JetpT >= pTL[p])npT=p;
    }
    if(JetpT < pTL[0] || JetpT > pTL[(pTL.size()-1)])npT=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaL.size()-1); ++q){
      if(fabs(JetEta) < EtaL[q+1] && fabs(JetEta) >= EtaL[q])nEta=q;
    }
    if(fabs(JetEta) < EtaL[0] || fabs(JetEta) > EtaL[(EtaL.size()-1)])nEta=-99; //For Jets outside the measured SF range

    for(int p = 0; p<(pTLFast.size()-1); ++p){
      if(JetpT < pTLFast[p+1] && JetpT >= pTLFast[p])npTFast=p;
    }
    if(JetpT < pTLFast[0] || JetpT > pTLFast[(pTLFast.size()-1)])npTFast=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaLFast.size()-1); ++q){
      if(fabs(JetEta) < EtaLFast[q+1] && fabs(JetEta) >= EtaLFast[q])nEtaFast=q;
    }
    if(fabs(JetEta) < EtaLFast[0] || fabs(JetEta) > EtaLFast[(EtaLFast.size()-1)])nEtaFast=-99; //For Jets outside the measured SF range

    if(npTFast!=-99 && nEtaFast!=-99){
      bool isBTagged=true;
      if(JetTag==0)isBTagged=false;
      if(!isPOG){
        JetTagNew=NewTag(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
	JetTagNewFast=NewTag(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewBUpFast=NewTag(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewBDownFast=NewTag(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
	JetTagNewCUpFast=NewTag(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewCDownFast=NewTag(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewLUpFast=NewTag(isBTagged, 0, 1, 1, SFLUpFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewLDownFast=NewTag(isBTagged, 0, 1, 1, SFLDownFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
      }else{
        JetTagNew=NewTagPOG(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
	JetTagNewFast=NewTagPOG(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewBUpFast=NewTagPOG(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewBDownFast=NewTagPOG(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
	JetTagNewCUpFast=NewTagPOG(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewCDownFast=NewTagPOG(isBTagged, 0, 1, 1, SFLFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewLUpFast=NewTagPOG(isBTagged, 0, 1, 1, SFLUpFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
        JetTagNewLDownFast=NewTagPOG(isBTagged, 0, 1, 1, SFLDownFast[npTFast][nEtaFast], EffLdFast[npTFast][nEtaFast]);
      }
    }

    if(npT!=-99 && nEta!=-99){
      bool isBTagged=true;
      if(!isPOG){
	if(JetTagNew==0)isBTagged=false;
	JetTagNew=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBUp=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBDown=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLUp=NewTag(isBTagged, 0, 1, 1, SFLUp[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLDown=NewTag(isBTagged, 0, 1, 1, SFLDown[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewBUpFast==0)isBTagged=false;
	JetTagNewBUpFast=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewBDownFast==0)isBTagged=false;
	JetTagNewBDownFast=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewCUpFast==0)isBTagged=false;
	JetTagNewCUpFast=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewCDownFast==0)isBTagged=false;
        JetTagNewCDownFast=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewLUpFast==0)isBTagged=false;
	JetTagNewLUpFast=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewLDownFast==0)isBTagged=false;
        JetTagNewLDownFast=NewTag(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
      }else{
	if(JetTagNew==0)isBTagged=false;
	JetTagNew=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBUp=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewBDown=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLUp=NewTagPOG(isBTagged, 0, 1, 1, SFLUp[npT][nEta], EffLd[npT][nEta]);
	JetTagNewLDown=NewTagPOG(isBTagged, 0, 1, 1, SFLDown[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewBUpFast==0)isBTagged=false;
	JetTagNewBUpFast=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewBDownFast==0)isBTagged=false;
        JetTagNewBDownFast=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
        if(JetTagNewCUpFast==0)isBTagged=false;
        JetTagNewCUpFast=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
        if(JetTagNewCDownFast==0)isBTagged=false;
        JetTagNewCDownFast=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
	if(JetTagNewLUpFast==0)isBTagged=false;
	JetTagNewLUpFast=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
	isBTagged=true;
        if(JetTagNewLDownFast==0)isBTagged=false;
        JetTagNewLDownFast=NewTagPOG(isBTagged, 0, 1, 1, SFL[npT][nEta], EffLd[npT][nEta]);
      }
    }

  }else if(JetFlavor==2){
    for(int p = 0; p<(pTB.size()-1); ++p){
      if(JetpT < pTB[p+1] && JetpT >= pTB[p])npT=p;
    }
    if(JetpT < pTB[0] || JetpT > pTB[(pTB.size()-1)])npT=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaB.size()-1); ++q){
      if(fabs(JetEta) < EtaB[q+1] && fabs(JetEta) >= EtaB[q])nEta=q;
    }
    if(fabs(JetEta) < EtaB[0] || fabs(JetEta) > EtaB[(EtaB.size()-1)])nEta=-99; //For Jets outside the measured SF range

    for(int p = 0; p<(pTBFast.size()-1); ++p){
      if(JetpT < pTBFast[p+1] && JetpT >= pTBFast[p])npTFast=p;
    }
    if(JetpT < pTBFast[0] || JetpT > pTBFast[(pTBFast.size()-1)])npTFast=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaBFast.size()-1); ++q){
      if(fabs(JetEta) < EtaBFast[q+1] && fabs(JetEta) >= EtaBFast[q])nEtaFast=q;
    }
    if(fabs(JetEta) < EtaBFast[0] || fabs(JetEta) > EtaBFast[(EtaBFast.size()-1)])nEtaFast=-99; //For Jets outside the measured SF range        

    if(npTFast!=-99 && nEtaFast!=-99){
      bool isBTagged=true;
      if(JetTag==0)isBTagged=false;
      if(!isPOG){
        JetTagNew=NewTag(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewFast=NewTag(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewBUpFast=NewTag(isBTagged, 2, SFBUpFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewBDownFast=NewTag(isBTagged, 2, SFBDownFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewCUpFast=NewTag(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewCDownFast=NewTag(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLUpFast=NewTag(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLDownFast=NewTag(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
      }else{
        JetTagNew=NewTagPOG(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewFast=NewTagPOG(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewBUpFast=NewTagPOG(isBTagged, 2, SFBUpFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewBDownFast=NewTagPOG(isBTagged, 2, SFBDownFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewCUpFast=NewTagPOG(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewCDownFast=NewTagPOG(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLUpFast=NewTagPOG(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLDownFast=NewTagPOG(isBTagged, 2, SFBFast[npTFast][nEtaFast], EffBdFast[npTFast][nEtaFast], 1, 1);
      }
    }

    if(npT!=-99 && nEta!=-99){
      bool isBTagged=true;
      if(!isPOG){
	if(JetTagNew==0)isBTagged=false;
	JetTagNew=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTag(isBTagged, 2, SFBUp[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTag(isBTagged, 2, SFBDown[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewBUpFast==0)isBTagged=false;
	JetTagNewBUpFast=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewBDownFast==0)isBTagged=false;
	JetTagNewBDownFast=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewCUpFast==0)isBTagged=false;
	JetTagNewCUpFast=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewCDownFast==0)isBTagged=false;
        JetTagNewCDownFast=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewLUpFast==0)isBTagged=false;
	JetTagNewLUpFast=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewLDownFast==0)isBTagged=false;
        JetTagNewLDownFast=NewTag(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
      }else{
	if(JetTagNew==0)isBTagged=false;
	JetTagNew=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTagPOG(isBTagged, 2, SFBUp[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTagPOG(isBTagged, 2, SFBDown[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewBUpFast==0)isBTagged=false;
	JetTagNewBUpFast=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewBDownFast==0)isBTagged=false;
        JetTagNewBDownFast=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewCUpFast==0)isBTagged=false;
        JetTagNewCUpFast=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewCDownFast==0)isBTagged=false;
        JetTagNewCDownFast=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewLUpFast==0)isBTagged=false;
        JetTagNewLUpFast=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewLDownFast==0)isBTagged=false;
        JetTagNewLDownFast=NewTagPOG(isBTagged, 2, SFB[npT][nEta], EffBd[npT][nEta], 1, 1);
      }
    }

  }else if(JetFlavor==1){
    for(int p = 0; p<(pTB.size()-1); ++p){
      if(JetpT < pTB[p+1] && JetpT >= pTB[p])npT=p;
    }
    if(JetpT < pTB[0] || JetpT > pTB[(pTB.size()-1)])npT=-99; //For Jets outside the measured SF range
    for(int q = 0; q<(EtaB.size()-1); ++q){
      if(fabs(JetEta) < EtaB[q+1] && fabs(JetEta) >= EtaB[q])nEta=q;
    }
    if(fabs(JetEta) < EtaB[0] || fabs(JetEta) > EtaB[(EtaB.size()-1)])nEta=-99; //For Jets outside the measured SF range          

    for(int p = 0; p<(pTBFast.size()-1); ++p){
      if(JetpT < pTBFast[p+1] && JetpT >= pTBFast[p])npTFast=p;
    }
    if(JetpT < pTBFast[0] || JetpT > pTBFast[(pTBFast.size()-1)])npTFast=-99; //For Jets outside the measured SF range                      
    for(int q = 0; q<(EtaBFast.size()-1); ++q){
      if(fabs(JetEta) < EtaBFast[q+1] && fabs(JetEta) >= EtaBFast[q])nEtaFast=q;
    }
    if(fabs(JetEta) < EtaBFast[0] || fabs(JetEta) > EtaBFast[(EtaBFast.size()-1)])nEtaFast=-99; //For Jets outside the measured SF range      

    if(npTFast!=-99 && nEtaFast!=-99){
      bool isBTagged=true;
      if(JetTag==0)isBTagged=false;
      if(!isPOG){
        JetTagNew=NewTag(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewFast=NewTag(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewBUpFast=NewTag(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewBDownFast=NewTag(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewCUpFast=NewTag(isBTagged, 1, SFCUpFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewCDownFast=NewTag(isBTagged, 1, SFCDownFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLUpFast=NewTag(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLDownFast=NewTag(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
      }else{
        JetTagNew=NewTagPOG(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewFast=NewTagPOG(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewBUpFast=NewTagPOG(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewBDownFast=NewTagPOG(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewCUpFast=NewTagPOG(isBTagged, 1, SFCUpFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
	JetTagNewCDownFast=NewTagPOG(isBTagged, 1, SFCDownFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLUpFast=NewTagPOG(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
        JetTagNewLDownFast=NewTagPOG(isBTagged, 1, SFCFast[npTFast][nEtaFast], EffCdFast[npTFast][nEtaFast], 1, 1);
      }
    }

    if(npT!=-99 && nEta!=-99){
      bool isBTagged=true;
      if(!isPOG){
	if(JetTagNew==0)isBTagged=false;
	JetTagNew=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTag(isBTagged, 1, SFCUp[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTag(isBTagged, 1, SFCDown[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
	if(JetTagNewBUpFast==0)isBTagged=false;
	JetTagNewBUpFast=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewBDownFast==0)isBTagged=false;
	JetTagNewBDownFast=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewCUpFast==0)isBTagged=false;
	JetTagNewCUpFast=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewCDownFast==0)isBTagged=false;
        JetTagNewCDownFast=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewLUpFast==0)isBTagged=false;
	JetTagNewLUpFast=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewLDownFast==0)isBTagged=false;
        JetTagNewLDownFast=NewTag(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
      }else{
	if(JetTagNew==0)isBTagged=false;
	JetTagNew=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBUp=NewTagPOG(isBTagged, 1, SFCUp[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewBDown=NewTagPOG(isBTagged, 1, SFCDown[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLUp=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	JetTagNewLDown=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewBUpFast==0)isBTagged=false;
	JetTagNewBUpFast=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
	isBTagged=true;
        if(JetTagNewBDownFast==0)isBTagged=false;
        JetTagNewBDownFast=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
        isBTagged=true;
        if(JetTagNewCUpFast==0)isBTagged=false;
        JetTagNewCUpFast=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
        isBTagged=true;
        if(JetTagNewCDownFast==0)isBTagged=false;
        JetTagNewCDownFast=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
        isBTagged=true;
        if(JetTagNewLUpFast==0)isBTagged=false;
        JetTagNewLUpFast=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
        isBTagged=true;
        if(JetTagNewLDownFast==0)isBTagged=false;
        JetTagNewLDownFast=NewTagPOG(isBTagged, 1, SFB[npT][nEta], EffCd[npT][nEta], 1, 1);
      }
    }


  }

  Tags.push_back(JetTagNew);
  Tags.push_back(JetTagNewBUp);
  Tags.push_back(JetTagNewBDown);
  Tags.push_back(JetTagNewLUp);
  Tags.push_back(JetTagNewLDown);
  Tags.push_back(JetTagNewBUpFast);
  Tags.push_back(JetTagNewBDownFast);
  Tags.push_back(JetTagNewCUpFast);
  Tags.push_back(JetTagNewCDownFast);
  Tags.push_back(JetTagNewLUpFast);
  Tags.push_back(JetTagNewLDownFast);
  Tags.push_back(JetTagNewFast);

  return Tags;
  
}


void SUSYMultiTop::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  // PF block
  double MET;
  double ST;
  double SLT;
  int nJets;
  int JetMultiplicity;
  int EleFakeJetMultiplicity;
  int MuFakeJetMultiplicity;
  int PUJetMultiplicity;
  int NLooseMuons, NTightMuons; 
  int N80Eles, N95Eles;
  int nEleJets;
  int nBTagJets, nB, nBTagJetsCSV, nBCSV, nBTagJetsLoose;
  double PUmva[20];
  double PUpTPull[20];
  int PUindex[20];
  int PUmatchedjets;
  int PUunmatchedjets;
  double DRmin[20];
  double PUunpT[20];
  double PUunEta[20];
  double PUunmva[20];
  int PUunIndex[20];
  double JetpT[20];
  double JetEta[20];
  int JetTag[20];
  int JetTagNewBUp[20];
  int JetTagNewBDown[20];
  int JetTagNewLUp[20];
  int JetTagNewLDown[20];
  int JetTagNewBUpFast[20];
  int JetTagNewBDownFast[20];
  int JetTagNewCUpFast[20];
  int JetTagNewCDownFast[20];
  int JetTagNewLUpFast[20];
  int JetTagNewLDownFast[20];
  int JetTagLoose[20];
  int JetTagCSV[20];
  int JetFlavor[20];
  int JetTagNew[20];
  int JetTagNewFast[20];
  double HT;
  double MRpT;
  double RpT;
  double EleJetDR[1000];
  int PFgoodR;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double massPFHem;
  double PFR;
  double PFMR;
  double PFMT;
  double phiMETLept;
  double etaMETLept;
  double wLep;
  double wLepMu, wLepEle;
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double MuIso[2];
  double EleTrackIso[2], EleECALIso[2], EleHCALIso[2];
  int n_PV;
  double Weight=_weight;
  double mg, mst, mchi;
  double puW, puWFull;
  double EvW;
  double EvWFull;
  double wLepFull;

  //Prepare counters for efficiency
  double Npassed_In = 0;
  double Npassed_InPU = 0;
  //HLT
  double Npassed_HLT = 0;
  double Npassed_PV = 0;
  double Npassed_MET= 0;
  //Jets
  double NpassedPF_4Jet = 0;
  double NpassedPF_CenJet=0;
  double NpassedPF_8Jet=0;
  double NpassedPF_1Jet=0;
  double NpassedPF_2Jet=0;
  double NpassedPF_3Jet=0;
  double NpassedPF_5Jet=0;
  double NpassedPF_6Jet=0;
  double NpassedPF_7Jet=0;
  //Hemispheres
  double NpassedPF_beta=0;
  double NpassedPF_DeltaPhi=0;
  //Leptons
  double Npassed_TightMu=0;
  double Npassed_TwoMu=0;
  double Npassed_ThreeMu=0;
  double Npassed_FourMu=0;
  double Npassed_Ele=0;
  double Npassed_TwoEle=0;
  double Npassed_EleMu=0;
  //B-tag
  double Npassed_Btag=0;
  double Npassed_1b=0;
  double Npassed_2b=0;
  double Npassed_3b=0;
  double Npassed_4b=0;

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");

  //prepare tree for SMS eff
  TTree* FullSMSTree = new TTree("FullSMSTree", "FullSMSTree");
  FullSMSTree->Branch("mg", &mg, "mg/D");
  FullSMSTree->Branch("mchi", &mchi, "mchi/D");
  FullSMSTree->Branch("puW", &puWFull, "puW/D");
  FullSMSTree->Branch("wLep", &wLepFull, "wLep/D");
  FullSMSTree->Branch("EventWeight", &EvWFull, "EventWeight/D");

  //prepare b-tag tree
  TTree* TagTree= new TTree("TagTree", "TagTree");
  TagTree->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
  TagTree->Branch("JetpT", JetpT, "JetpT[JetMultiplicity]/D");
  TagTree->Branch("JetEta", JetEta, "JetEta[JetMultiplicity]/D");
  TagTree->Branch("JetTag", JetTagNew, "JetTag[JetMultiplicity]/I");
  TagTree->Branch("JetTagBUp", JetTagNewBUp, "JetTagBUp[JetMultiplicity]/I");
  TagTree->Branch("JetTagBDown", JetTagNewBDown, "JetTagBDown[JetMultiplicity]/I");
  TagTree->Branch("JetTagLUp", JetTagNewLUp, "JetTagLUp[JetMultiplicity]/I");
  TagTree->Branch("JetTagLDown", JetTagNewLDown, "JetTagLDown[JetMultiplicity]/I");
  TagTree->Branch("JetFlavor", JetFlavor, "JetFlavor[JetMultiplicity]/I");
  TagTree->Branch("puW", &puW, "puW/D");
  TagTree->Branch("wLep", &wLep, "wLep/D");
  TagTree->Branch("EventWeight", &EvW, "EventWeight/D");

  //prepare b-tag tree (POG Algorithm)
  TTree* TagTreePOG= new TTree("TagTreePOG", "TagTreePOG");
  TagTreePOG->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
  TagTreePOG->Branch("JetpT", JetpT, "JetpT[JetMultiplicity]/D");
  TagTreePOG->Branch("JetEta", JetEta, "JetEta[JetMultiplicity]/D");
  TagTreePOG->Branch("JetTag", JetTagNew, "JetTag[JetMultiplicity]/I");
  TagTreePOG->Branch("JetTagBUp", JetTagNewBUp, "JetTagBUp[JetMultiplicity]/I");
  TagTreePOG->Branch("JetTagBDown", JetTagNewBDown, "JetTagBDown[JetMultiplicity]/I");
  TagTreePOG->Branch("JetTagLUp", JetTagNewLUp, "JetTagLUp[JetMultiplicity]/I");
  TagTreePOG->Branch("JetTagLDown", JetTagNewLDown, "JetTagLDown[JetMultiplicity]/I");
  TagTreePOG->Branch("JetFlavor", JetFlavor, "JetFlavor[JetMultiplicity]/I");
  TagTreePOG->Branch("puW", &puW, "puW/D");
  TagTreePOG->Branch("wLep", &wLep, "wLep/D");
  TagTreePOG->Branch("EventWeight", &EvW, "EventWeight/D");

  //prepare b-tag tree SMS
  TTree* TagTreeSMS= new TTree("TagTreeSMS", "TagTreeSMS");
  TagTreeSMS->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
  TagTreeSMS->Branch("JetpT", JetpT, "JetpT[JetMultiplicity]/D");
  TagTreeSMS->Branch("JetEta", JetEta, "JetEta[JetMultiplicity]/D");
  TagTreeSMS->Branch("JetTag", JetTagNew, "JetTag[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagFast", JetTagNewFast, "JetTagFast[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagBUp", JetTagNewBUp, "JetTagBUp[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagBDown", JetTagNewBDown, "JetTagBDown[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagLUp", JetTagNewLUp, "JetTagLUp[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagLDown", JetTagNewLDown, "JetTagLDown[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagBUpFast", JetTagNewBUpFast, "JetTagBUpFast[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagBDownFast", JetTagNewBDownFast, "JetTagBDownFast[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagCUpFast", JetTagNewCUpFast, "JetTagCUpFast[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagCDownFast", JetTagNewCDownFast, "JetTagCDownFast[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagLUpFast", JetTagNewLUpFast, "JetTagLUpFast[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetTagLDownFast", JetTagNewLDownFast, "JetTagLDown[JetMultiplicity]/I");
  TagTreeSMS->Branch("JetFlavor", JetFlavor, "JetFlavor[JetMultiplicity]/I");
  TagTreeSMS->Branch("puW", &puW, "puW/D");
  TagTreeSMS->Branch("wLep", &wLep, "wLep/D");
  TagTreeSMS->Branch("EventWeight", &EvW, "EventWeight/D");
  TagTreeSMS->Branch("mg", &mg, "mg/D");
  TagTreeSMS->Branch("mchi", &mchi, "mchi/D");

  //prepare b-tag tree SMS (POG Algorithm)
  TTree* TagTreeSMSPOG= new TTree("TagTreeSMSPOG", "TagTreeSMSPOG");
  TagTreeSMSPOG->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
  TagTreeSMSPOG->Branch("JetpT", JetpT, "JetpT[JetMultiplicity]/D");
  TagTreeSMSPOG->Branch("JetEta", JetEta, "JetEta[JetMultiplicity]/D");
  TagTreeSMSPOG->Branch("JetTag", JetTagNew, "JetTag[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagFast", JetTagNewFast, "JetTagFast[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagBUp", JetTagNewBUp, "JetTagBUp[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagBDown", JetTagNewBDown, "JetTagBDown[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagLUp", JetTagNewLUp, "JetTagLUp[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagLDown", JetTagNewLDown, "JetTagLDown[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagBUpFast", JetTagNewBUpFast, "JetTagBUpFast[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagBDownFast", JetTagNewBDownFast, "JetTagBDownFast[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagCUpFast", JetTagNewCUpFast, "JetTagCUpFast[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagCDownFast", JetTagNewCDownFast, "JetTagCDownFast[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagLUpFast", JetTagNewLUpFast, "JetTagLUpFast[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetTagLDownFast", JetTagNewLDownFast, "JetTagLDown[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("JetFlavor", JetFlavor, "JetFlavor[JetMultiplicity]/I");
  TagTreeSMSPOG->Branch("puW", &puW, "puW/D");
  TagTreeSMSPOG->Branch("wLep", &wLep, "wLep/D");
  TagTreeSMSPOG->Branch("EventWeight", &EvW, "EventWeight/D");
  TagTreeSMSPOG->Branch("mg", &mg, "mg/D");
  TagTreeSMSPOG->Branch("mchi", &mchi, "mchi/D");

  //prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("nPV", &nPV, "nPV/I");

  //PF block
  outTree->Branch("puW", &puW, "puW/D");
  outTree->Branch("MET", &MET, "MET/D");
  outTree->Branch("wLep", &wLep, "wLep/D");
  outTree->Branch("ST", &ST, "ST/D");
  outTree->Branch("SLT", &SLT, "SLT/D");
  outTree->Branch("HT", &HT, "HT/D");
  outTree->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity/I");
  outTree->Branch("nBTagJets", &nBTagJets, "nBTagJets/I");
  outTree->Branch("nBTagJetsLoose", &nBTagJetsLoose, "nBTagJetsLoose/I");
  outTree->Branch("nB", &nB, "nB/I");
  outTree->Branch("nBTagJetsCSV", &nBTagJetsCSV, "nBTagJetsCSV/I");
  outTree->Branch("nBCSV", &nBCSV, "nBCSV/I");
  outTree->Branch("N80Eles", &N80Eles, "N80Eles/I");
  outTree->Branch("N95Eles", &N95Eles, "N95Eles/I");
  outTree->Branch("NLooseMuons", &NLooseMuons, "NLooseMuons/I");
  outTree->Branch("NTightMuons", &NTightMuons, "NTightMuons/I");
  outTree->Branch("EleFakeJetMultiplicity", &EleFakeJetMultiplicity, "EleFakeJetMultiplicity/I");
  outTree->Branch("MuFakeJetMultiplicity", &MuFakeJetMultiplicity, "MuFakeJetMultiplicity/I");
  outTree->Branch("PFMT", &PFMT, "PFMT/D");
  outTree->Branch("MuIso", MuIso, "MuIso[NLooseMuons]/D");
  outTree->Branch("EleTrackIso", EleTrackIso, "EleTrackIso[N95Eles]/D");
  outTree->Branch("EleECALIso", EleECALIso, "EleECALIso[N95Eles]/D");
  outTree->Branch("EleHCALIso", EleHCALIso, "EleHCALIso[N95Eles]/D");
  outTree->Branch("JetpT", JetpT, "JetpT[JetMultiplicity]/D");
  outTree->Branch("JetEta", JetEta, "JetEta[JetMultiplicity]/D");
  outTree->Branch("JetTag", JetTag, "JetTag[JetMultiplicity]/I");
  outTree->Branch("JetTagLoose", JetTagLoose, "JetTagLoose[JetMultiplicity]/I");
  outTree->Branch("JetTagCSV", JetTagCSV, "JetTagCSV[JetMultiplicity]/I");
  outTree->Branch("JetFlavor", JetFlavor, "JetFlavor[JetMultiplicity]/I");
  outTree->Branch("weight", &Weight, "weight/D");
  outTree->Branch("EventWeight", &EvW, "EventWeight/D");
  outTree->Branch("PUmatchedjets", &PUmatchedjets, "PUmatchedjets/I");
  outTree->Branch("PUunmatchedjets", &PUunmatchedjets, "PUunmatchedjets/I");
  outTree->Branch("PUJetMultiplicity", &PUJetMultiplicity, "PUJetMultiplicity/I");
  outTree->Branch("PUmva", PUmva, "PUmva[PUmatchedjets]/D");
  outTree->Branch("PUunmva", PUunmva, "PUunmva[PUunmatchedjets]/D");
  outTree->Branch("PUunpT", PUunpT, "PUunpT[PUunmatchedjets]/D");
  outTree->Branch("PUunEta", PUunEta, "PUunEta[PUunmatchedjets]/D");
  outTree->Branch("PUunIndex", PUunIndex, "PUunIndex[PUunmatchedjets]/I");
  outTree->Branch("DRmin", DRmin, "DRmin[JetMultiplicity]/D");
  outTree->Branch("PUindex", PUindex, "PUindex[PUmatchedjets]/I");
  outTree->Branch("PUpTPull", PUpTPull, "PUpTPull[PUmatchedjets]/D");
  //SMS
  outTree->Branch("mg", &mg, "mg/D");
  outTree->Branch("mchi", &mchi, "mchi/D");

  double weightII = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  Long64_t nbytes = 0;
  Long64_t nb = 0;

  //Prepare efftree
  TTree* effTree = new TTree("effTree", "effTree");
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_InPU",      &Npassed_InPU,      "Npassed_InPU/D");
  effTree->Branch("Npassed_HLT",      &Npassed_HLT,      "Npassed_HLT/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("NpassedPF_1Jet",      &NpassedPF_1Jet,      "NpassedPF_1Jet/D");
  effTree->Branch("NpassedPF_2Jet",      &NpassedPF_2Jet,      "NpassedPF_2Jet/D");
  effTree->Branch("NpassedPF_3Jet",      &NpassedPF_3Jet,      "NpassedPF_3Jet/D");
  effTree->Branch("NpassedPF_4Jet",      &NpassedPF_4Jet,      "NpassedPF_4Jet/D");
  effTree->Branch("NpassedPF_5Jet",      &NpassedPF_5Jet,      "NpassedPF_5Jet/D");
  effTree->Branch("NpassedPF_6Jet",      &NpassedPF_6Jet,      "NpassedPF_6Jet/D");
  effTree->Branch("NpassedPF_7Jet",      &NpassedPF_7Jet,      "NpassedPF_7Jet/D");
  effTree->Branch("NpassedPF_CenJet",      &NpassedPF_CenJet,      "NpassedPF_CenJet/D");
  effTree->Branch("NpassedPF_8Jet",      &NpassedPF_8Jet,      "NpassedPF_8Jet/D");
  effTree->Branch("Npassed_TightMu",      &Npassed_TightMu,      "Npassed_TightMu/D");
  effTree->Branch("Npassed_TwoMu",      &Npassed_TwoMu,      "Npassed_TwoMu/D");
  effTree->Branch("Npassed_ThreeMu",      &Npassed_ThreeMu,      "Npassed_ThreeMu/D");
  effTree->Branch("Npassed_FourMu",      &Npassed_FourMu,      "Npassed_FourMu/D");
  effTree->Branch("Npassed_Ele",      &Npassed_Ele,      "Npassed_Ele/D");
  effTree->Branch("Npassed_TwoEle",      &Npassed_TwoEle,      "Npassed_TwoEle/D");
  effTree->Branch("Npassed_EleMu",      &Npassed_EleMu,      "Npassed_EleMu/D");
  effTree->Branch("Npassed_MET",      &Npassed_MET,      "Npassed_MET/D");
  effTree->Branch("Npassed_Btag",      &Npassed_Btag,      "Npassed_Btag/D");
  effTree->Branch("Npassed_1b",      &Npassed_1b,      "Npassed_1b/D");
  effTree->Branch("Npassed_2b",      &Npassed_2b,      "Npassed_2b/D");
  effTree->Branch("Npassed_3b",      &Npassed_3b,      "Npassed_3b/D");
  effTree->Branch("Npassed_4b",      &Npassed_4b,      "Npassed_4b/D");

  //Prepare the LumiTree
  TTree * LumiTree = new TTree("LumiTree", "LumiTree");
  LumiTree->Branch("run", &run, "run/D");
  LumiTree->Branch("ls", &ls, "ls/D");

  cout << "Number of entries = " << stop << endl;

  //Build PU reweighting class
  //2012
  reweight::LumiReWeighting LumiWeights( "/afs/cern.ch/user/e/emanuele/workspace/public/pileup/s7pileup200.root",
					 "/afs/cern.ch/user/e/emanuele/workspace/public/pileup/puRun2012_5100ipb_71.root",
					 "hNPU","pileup");

  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;
    if (jentry%10000 == 0) cout << " Fraction of events analyzed "<< 100*(double(jentry)/double(stop))<<"%"<<endl;
    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE

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
    
    if(_isData)reloadTriggerMask(runNumber);
    
    //PU reweighting                                                                                                                                 
    if(!_isData)puW = LumiWeights.weight(nPU[1]);
    else puW=1;

    int iTightMuons[10], tightMuonCounter=0;
    int iLooseMuons[10], looseMuonCounter=0;
    int i95Eles[10], Ele95Counter=0;
    int i80Eles[10], Ele80Counter=0;

    //Muons
    //N.B. We use the 8TeV Loose Muon definition to have a looser isolation (0.4)
    SLT=0.;
    for(int j=0; j<nMuon;j++){
      double EtaMuon=etaMuon[j];
      double muonPt=sqrt(pxMuon[j]*pxMuon[j]+pyMuon[j]*pyMuon[j]);      
      if(Vecbos::isLooseMuon(j) && muonPt>35. && fabs(EtaMuon)< 2.1){
	iTightMuons[tightMuonCounter]=j;
	tightMuonCounter++;
	if(!_isData)wLepMu = Vecbos::getOfflineEff(muonPt, etaMuon[j], "Muon");
	else wLepMu=1.0;
      }
    }
    
    double IECAL; 
    double IHCAL;
    double ITRK;

    for(int j=0; j<nMuon;j++){
      double muonPt=sqrt(pxMuon[j]*pxMuon[j]+pyMuon[j]*pyMuon[j]);      
      if(Vecbos::isLooseMuon(j) && muonPt>20.){
	iLooseMuons[looseMuonCounter]=j;
	IECAL = emEt03Muon[j];
	IHCAL = hadEt03Muon[j];
	ITRK  = sumPt03Muon[j];
	MuIso[looseMuonCounter]=((IECAL+IHCAL+ITRK)/muonPt);
	looseMuonCounter++;
	SLT+=muonPt;
      }
    }
    
    //Electrons
    for(int k=0; k< nEle; k++){
      double ElePt=sqrt(pxEle[k]*pxEle[k]+pyEle[k]*pyEle[k]);
      if(Vecbos::isTightElectron(k) && ElePt>35. && fabs(etaEle[k])<2.5){
	i80Eles[Ele80Counter]=k;
	Ele80Counter++;
	if(!_isData)wLepEle = Vecbos::getOfflineEff(ElePt, etaEle[k], "Electron");
	else wLepEle = 1.0;
      }
    }
    
    double ptLept;
    if(AnalysisSelector==1 && tightMuonCounter > 0){
      wLep=wLepMu;
      ptLept=sqrt(pxMuon[iTightMuons[0]]*pxMuon[iTightMuons[0]]+pyMuon[iTightMuons[0]]*pyMuon[iTightMuons[0]]);
    }else if(AnalysisSelector==3 && Ele80Counter > 0){
      wLep=wLepEle;
      ptLept=sqrt(pxEle[i80Eles[0]]*pxEle[i80Eles[0]]+pyEle[i80Eles[0]]*pyEle[i80Eles[0]]);
    }else{
      wLep=1;
      ptLept=0;
    }
    
    int iSC;
    double ET;
    double trackISO;
    double ECALISO;
    double HCALISO;
    
    for(int k=0; k< nEle; k++){
      double ElePt=sqrt(pxEle[k]*pxEle[k]+pyEle[k]*pyEle[k]);
      if(Vecbos::isLooseElectron(k) && ElePt>20.){
	i95Eles[Ele95Counter]=k;
	iSC=superClusterIndexEle[k];
	ET=energySC[iSC]/cosh(etaSC[iSC]);
	trackISO = (dr03TkSumPtEle[k])/ET;
	ECALISO = (dr03EcalRecHitSumEtEle[k])/ET;
	HCALISO = (dr03HcalTowerSumEtEle[k])/ET;
	EleTrackIso[Ele95Counter]=trackISO;
	EleECALIso[Ele95Counter]=ECALISO;
	EleHCALIso[Ele95Counter]=HCALISO;
	Ele95Counter++;
	SLT+=ElePt;
      }
    }
    
    Npassed_In += weightII;
    
    Npassed_InPU += puW*wLep;


    double m0=9999, m12=9999, mc=9999;
     
    if(!_isData && isSMS && isT1tttt){

      //find the simplified model parameters for T1tttt                                                                                              
      std::vector<std::string>::const_iterator c_begin = commentLHE->begin();
      std::vector<std::string>::const_iterator c_end = commentLHE->end();
      for(std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
	size_t found = (*cit).find("T1tttt"); 
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
    }else if(!isT1tttt && isSMS && !_isData){
      string mass=_sample.substr(4, _sample.size());
      m12=(double)atof(mass.c_str());
    }

    mg=m12;
    mst=m0;
    mchi=mc;

    //HLT
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();

    if(!_isData)passedHLT = 1;
    if(!passedHLT)continue;

    Npassed_HLT += weightII;

    if(_isData)LumiTree->Fill();
    
    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    Npassed_PV += weightII;

    //Jets 
    vector<TLorentzVector> PFJet;
    HT=0.;
    double Pxtot=0., Pytot=0.;
    vector <int> iPFJet;
    for(int i=0; i< nAK5PFNoPUJet; i++) {
      TLorentzVector myJet(pxAK5PFNoPUJet[i], pyAK5PFNoPUJet[i], pzAK5PFNoPUJet[i], energyAK5PFNoPUJet[i]);
      TLorentzVector myJetCorr;
      if(!isJESSystematic)myJetCorr=myJet;
      if(isJESSystematic && JES==1) myJetCorr=GetJESCorrected(myJet, "Up");
      else if(isJESSystematic && JES==-1) myJetCorr=GetJESCorrected(myJet, "Down");
      if(myJetCorr.Pt()>30. && fabs(myJetCorr.Eta())< 2.4) {
	HT+=myJetCorr.Pt();
	Pxtot+=myJetCorr.Px();
	Pytot+=myJetCorr.Py();
	PFJet.push_back(myJetCorr);
	iPFJet.push_back(i);
      }
    }
    
    double DREleJet;
    vector<int> BadPFEleJet;
    int counter, double_check;
    for(int l=0; l<Ele80Counter; l++){
      for(int n=0; n < PFJet.size(); n++){
	TLorentzVector myJet = PFJet.at(n);
	DREleJet=sqrt((etaEle[i80Eles[l]]-myJet.Eta())*(etaEle[i80Eles[l]]-myJet.Eta())+(phiEle[i80Eles[l]]-myJet.Phi())*(phiEle[i80Eles[l]]-myJet.Phi()));
	double_check=0;
	counter=0;
	if(DREleJet < 0.5){
	  if(BadPFEleJet.size()!=0){
	    for(int f=0; f<BadPFEleJet.size(); f++){
	       if(iPFJet.at(n)!=BadPFEleJet.at(f))counter++;
	    }
	    if(counter == BadPFEleJet.size())double_check=1;
	    if(double_check!=0)BadPFEleJet.push_back(iPFJet.at(n));
	  }else{
	    BadPFEleJet.push_back(iPFJet.at(n));
	   }
	}
      }
    }
    
    double DRMuJet;
    vector<int> BadPFMuJet;
     int counterII, double_checkII, counter_Ele, double_checkEle;
     for(int l=0; l<looseMuonCounter; l++){
       for(int n=0; n < PFJet.size(); n++){
	 TLorentzVector myJet = PFJet.at(n);
	 DRMuJet=sqrt((etaMuon[iLooseMuons[l]]-myJet.Eta())*(etaMuon[iLooseMuons[l]]-myJet.Eta())+(phiMuon[iLooseMuons[l]]-myJet.Phi())*(phiMuon[iLooseMuons[l]]-myJet.Phi()));
	 double_checkII=0;
	 counterII=0;
	 counter_Ele=0;
	 double_checkEle=0;
	 if(DRMuJet < 0.5){
	   for(int z=0; z<BadPFEleJet.size(); z++){
             if(iPFJet.at(n)!=BadPFEleJet.at(z))counter_Ele++;
           }
	   if(counter_Ele==BadPFEleJet.size())double_checkEle=1;
	   if(BadPFMuJet.size()!=0){
	   for(int f=0; f<BadPFMuJet.size(); f++){
	     if(iPFJet.at(n)!=BadPFMuJet.at(f))counterII++;
	   }
	   if(counterII == BadPFMuJet.size())double_checkII=1;
	   if(double_checkII!=0 && double_checkEle!=0)BadPFMuJet.push_back(iPFJet.at(n));
	   }else{
	     if(double_checkEle!=0)BadPFMuJet.push_back(iPFJet.at(n));
	   }
	 }
       }
     }

     //Use MVA to determine whether a jet comes from PU or not
     bool goodJet;
     vector<double> BadPUJet;
     int igoodJet;
     double pTNoPU;
     int indexNoPU;
     PUmatchedjets=0;
     for(int j=0; j<nAK5PFPUcorrJet; j++){
       TLorentzVector myJet(pxAK5PFPUcorrJet[j], pyAK5PFPUcorrJet[j], pzAK5PFPUcorrJet[j], energyAK5PFPUcorrJet[j]);
       goodJet=false;
       for(int b=0; b < iPFJet.size(); b++){
         int n=iPFJet.at(b);
         int b_counter=0;
         for(int f=0; f<BadPFMuJet.size(); f++){
           if(n!=BadPFMuJet.at(f))b_counter++;
         }
         for(int d=0; d<BadPFEleJet.size(); d++){
           if(n!=BadPFEleJet.at(d))b_counter++;
         }
         if(b_counter==(BadPFEleJet.size()+BadPFMuJet.size())){
           double phi=PFJet[b].Phi();
           double eta=PFJet[b].Eta();
           double etaPUcorr=myJet.Eta();
           double phiPUcorr=myJet.Phi();
           double DR=sqrt((eta-etaPUcorr)*(eta-etaPUcorr)+(phi-phiPUcorr)*(phi-phiPUcorr));
	   if(DR<0.5){
	     goodJet=true;
	     igoodJet=b;
	     pTNoPU=PFJet[b].Pt();
	     indexNoPU=iPFJet.at(b);
	   }
         }
       }
       if(goodJet){
	 PUmva[PUmatchedjets]=jetIdMvaPhilV1AK5PFPUcorrJet[j];
	 double pull=(pTNoPU-myJet.Pt())/pTNoPU;
	 PUpTPull[PUmatchedjets]=pull;
	 PUindex[PUmatchedjets]=indexNoPU;
	 PUmatchedjets++;
       }
       if(goodJet && !isLooseJetMva(myJet.Pt(), myJet.Eta(), jetIdMvaPhilV1AK5PFPUcorrJet[j]))BadPUJet.push_back(igoodJet);
     }

     PUunmatchedjets=0;
     for(int j=0;j<iPFJet.size();j++){
       int m=iPFJet.at(j);
       TLorentzVector myJet(pxAK5PFNoPUJet[m], pyAK5PFNoPUJet[m], pzAK5PFNoPUJet[m], energyAK5PFNoPUJet[m]);
       int b_counter=0;
       for(int n=0; n< PUmatchedjets; n++){
	 if(PUindex[n]==m)b_counter++;
       }
       for(int f=0; f<BadPFMuJet.size(); f++){
	 if(m==BadPFMuJet.at(f))b_counter++;
       }
       for(int d=0; d<BadPFEleJet.size(); d++){
	 if(m==BadPFEleJet.at(d))b_counter++;
       }
       if(b_counter==0){
	 PUunmva[PUunmatchedjets]=jetIdMvaPhilV1AK5PFPUcorrJet[m];
	 PUunpT[PUunmatchedjets]=myJet.Pt();
	 PUunEta[PUunmatchedjets]=myJet.Eta();
	 PUunIndex[PUunmatchedjets]=m;
	 PUunmatchedjets++;
       }
     }

     for(int j=0;j<iPFJet.size();j++){
       for(int k=0; k<nAK5PFPUcorrJet; k++){
	 TLorentzVector myJet(pxAK5PFPUcorrJet[k], pyAK5PFPUcorrJet[k], pzAK5PFPUcorrJet[k], energyAK5PFPUcorrJet[k]);
	 double phi=PFJet[j].Phi();
         double eta=PFJet[j].Eta();
         double etaPUcorr=myJet.Eta();
         double phiPUcorr=myJet.Phi();
         double DR=sqrt((eta-etaPUcorr)*(eta-etaPUcorr)+(phi-phiPUcorr)*(phi-phiPUcorr));
	 if(k==0)DRmin[j]=DR;
	 else if(DR<DRmin[j])DRmin[j]=DR;
       }
     }

     //Subtract fake jets pT from HT and MRpT
     int c;
     for(int z=0; z< BadPFEleJet.size(); z++){
       c=BadPFEleJet.at(z);
       TLorentzVector myJet(pxAK5PFNoPUJet[c], pyAK5PFNoPUJet[c], pzAK5PFNoPUJet[c], energyAK5PFNoPUJet[c]);
       TLorentzVector myJetCorr;
       if(!isJESSystematic)myJetCorr=myJet;
       if(isJESSystematic && JES==1) myJetCorr=GetJESCorrected(myJet, "Up");
       else if(isJESSystematic && JES==-1) myJetCorr=GetJESCorrected(myJet, "Down");
       HT-=myJetCorr.Pt();
       Pxtot-=myJetCorr.Px();
       Pytot-=myJetCorr.Py();
     }
     for(int a=0; a < BadPFMuJet.size(); a++){
       c=BadPFMuJet.at(a);
       TLorentzVector myJet(pxAK5PFNoPUJet[c], pyAK5PFNoPUJet[c], pzAK5PFNoPUJet[c], energyAK5PFNoPUJet[c]);
       TLorentzVector myJetCorr;
       if(!isJESSystematic)myJetCorr=myJet;
       if(isJESSystematic && JES==1) myJetCorr=GetJESCorrected(myJet, "Up");
       else if(isJESSystematic && JES==-1) myJetCorr=GetJESCorrected(myJet, "Down");
       HT-=myJetCorr.Pt();
       Pxtot-=myJetCorr.Px();
       Pytot-=myJetCorr.Py();
     }
     for(int a=0; a < BadPUJet.size(); a++){
       c=iPFJet.at(BadPUJet.at(a));
       TLorentzVector myJet(pxAK5PFNoPUJet[c], pyAK5PFNoPUJet[c], pzAK5PFNoPUJet[c], energyAK5PFNoPUJet[c]);
       TLorentzVector myJetCorr;
       if(!isJESSystematic)myJetCorr=myJet;
       if(isJESSystematic && JES==1) myJetCorr=GetJESCorrected(myJet, "Up");
       else if(isJESSystematic && JES==-1) myJetCorr=GetJESCorrected(myJet, "Down");
       HT-=myJetCorr.Pt();
       Pxtot-=myJetCorr.Px();
       Pytot-=myJetCorr.Py();
     }


     //Add muon pT to compute MRpT and RpT
     for(int r=0; r<looseMuonCounter;r++){
       Pxtot+=pxMuon[iLooseMuons[r]];
       Pytot+=pyMuon[iLooseMuons[r]];
     }

     //Create arrays with the b-discriminators of the jets, in descending order
     //And save properties of the real jets in the tree after having given a dummy value to their entries
     for(int h=0; h < 20; h++){
       JetTagCSV[h]=-99;
       JetTag[h]=-99;
       JetEta[h]=-99;
       JetFlavor[h]=-99;
       JetpT[h]=-99;
     }
     int n;
     nBTagJets=0;
     nBTagJetsLoose=0;
     nB=0;
     nBTagJetsCSV=0;
     nBCSV=0;
     int nGoodJets=0;
     for(int b=0; b < iPFJet.size(); b++){
       n=iPFJet.at(b);
       int b_counter=0;
       for(int f=0; f<BadPFMuJet.size(); f++){
	 if(n!=BadPFMuJet.at(f))b_counter++;
       }
       for(int d=0; d<BadPFEleJet.size(); d++){
	 if(n!=BadPFEleJet.at(d))b_counter++;
       }
       for(int d=0; d<BadPUJet.size(); d++){
         if(b!=BadPUJet.at(d))b_counter++;
       }
       if(b_counter==(BadPFEleJet.size()+BadPFMuJet.size()+BadPUJet.size())){
	 if(trackCountingHighEffBJetTagsAK5PFNoPUJet[n] > 1.7){
           nBTagJetsLoose++;
           if(nGoodJets<20)JetTagLoose[nGoodJets]=1;
         }else{
           if(nGoodJets<20)JetTagLoose[nGoodJets]=0;
         }
	 if(trackCountingHighEffBJetTagsAK5PFNoPUJet[n] > 3.3){
	   nBTagJets++;
	   if(nGoodJets<20)JetTag[nGoodJets]=1;
	 }else{
	   if(nGoodJets<20)JetTag[nGoodJets]=0;
	 }
	 if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) {
           nBTagJetsCSV++;
           if(nGoodJets<20)JetTagCSV[nGoodJets]=1;
         }else{
           if(nGoodJets<20)JetTagCSV[nGoodJets]=0;
         }
	 TLorentzVector myJet(pxAK5PFNoPUJet[n], pyAK5PFNoPUJet[n], pzAK5PFNoPUJet[n], energyAK5PFNoPUJet[n]);
	 TLorentzVector myJetCorr;
	 if(!isJESSystematic)myJetCorr=myJet;
	 if(isJESSystematic && JES==1) myJetCorr=GetJESCorrected(myJet, "Up");
	 else if(isJESSystematic && JES==-1) myJetCorr=GetJESCorrected(myJet, "Down");
	 if(nGoodJets<20){
	   JetpT[nGoodJets]=myJetCorr.Pt();
	   JetEta[nGoodJets]=myJetCorr.Eta();
	   if(!_isData){
	     if(nGoodJets<20)JetFlavor[nGoodJets]=JetFlavorFull(myJetCorr);
	   }else{
	     if(nGoodJets<20)JetFlavor[nGoodJets]=99;
	   }
	 }
	 nGoodJets++;
       }
     }


     puWFull=puW;
     wLepFull=wLep;
     EvWFull=_weight*puW*wLep;
     if(isSMS)FullSMSTree->Fill();

     for(int x=0; x< nGoodJets; x++)if(JetTag[x]==1)nB++;
     for(int x=0; x< nGoodJets; x++)if(JetTagCSV[x]==1)nBCSV++;

     //Selection

#if AnalysisSelector == 1
     if(tightMuonCounter < 1)continue;
     if(Ele80Counter!=0)continue;
     Npassed_TightMu+=weightII;
#endif     

#if AnalysisSelector == 2
     if(looseMuonCounter!=2)continue;
     if(Ele80Counter!=0)continue;
     if(Ele95Counter!=0)continue;
     Npassed_TwoMu+=weightII;
#endif

#if AnalysisSelector == 3
     if(tightMuonCounter!=0)continue;
     if(Ele80Counter < 1)continue;
     //For an ele+3jets trigger
     //     if(PFJet.size()<3)continue;
     //     if((PFJet[0].Pt()< 40.0) || (PFJet[1].Pt()< 40.0) || (PFJet[2].Pt()< 40.0))continue;
     Npassed_Ele+=weightII;
#endif

#if AnalysisSelector == 4
     if(tightMuonCounter!=0)continue;
     if(looseMuonCounter!=0)continue;
     if(Ele95Counter!=2)continue;
     if(chargeEle[i95Eles[0]]!=chargeEle[i95Eles[1]])continue;
     Npassed_TwoEle+=weightII;
#endif

#if AnalysisSelector == 5
     if(looseMuonCounter!=1)continue;
     if(Ele80Counter!=1)continue;
     if(Ele95Counter!=1)continue;
     if(chargeEle[i80Eles[0]]!=chargeMuon[iLooseMuons[0]])continue;
     Npassed_EleMu+=weightII;
#endif

#if AnalysisSelector == 6
     bool ret=true;
     if(!_isData){
       for(int i=0; i<nMc; i++){
	 if(!ret)continue;
	 if((abs(idMc[i])==11 || abs(idMc[i])==13) && abs(idMc[mothMc[i]])==24)ret=false;
       }
     }
     if(ret)continue;
     //     if((tightMuonCounter < 1) && (Ele80Counter < 1))continue;
#endif

     EleFakeJetMultiplicity=BadPFEleJet.size();
     MuFakeJetMultiplicity=BadPFMuJet.size();
     PUJetMultiplicity=BadPUJet.size();
     JetMultiplicity=PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size();

     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()) == 1)NpassedPF_1Jet+=weightII;
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()) == 2)NpassedPF_2Jet+=weightII;
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()) == 3)NpassedPF_3Jet+=weightII;

#if isJustScaling == false
     if(PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size() < 4.) continue;
#endif

     Npassed_Btag+=weightII;
     
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size()) >= 4)NpassedPF_4Jet+= weightII;
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size()) == 5)NpassedPF_5Jet+=weightII;
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size()) >= 6)NpassedPF_6Jet+=weightII;
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size()) == 7)NpassedPF_7Jet+=weightII;

#if isBaseline == false
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size()) < 6.) continue;
#endif
     
     if((PFJet.size()-BadPFEleJet.size()-BadPFMuJet.size()-BadPUJet.size()) >= 8.)NpassedPF_8Jet+=weightII;
     
     if(sqrt(pxPFMet[0]*pxPFMet[0]+pyPFMet[0]*pyPFMet[0])>20.)Npassed_MET+=weightII;

     //fill output tree
     n_PV=nPV;

     MET=sqrt(pxPFMet[0]*pxPFMet[0]+pyPFMet[0]*pyPFMet[0]);
     PFMT=sqrt(2*MET*ptLept*(1-cos(phiMETLept)));
     ST=MET+SLT+HT;
     nJets=nAK5PFNoPUJet;
     nEleJets=0.;
     N80Eles = Ele80Counter;
     N95Eles = Ele95Counter;
     NLooseMuons = looseMuonCounter;
     NTightMuons = tightMuonCounter;

     run = runNumber;
     evNum = eventNumber;
     bx = eventNumber;
     ls = lumiBlock;
     orbit = orbitNumber;

     outTree->Fill();

     //Duplicate events based on the MC weight, correct tagging efficiency and save in the new tree 
     double SFBFast[pTBFast.size()-1][EtaBFast.size()-1];
     double SFCFast[pTBFast.size()-1][EtaBFast.size()-1];
     double SFLFast[pTLFast.size()-1][EtaLFast.size()-1];

     double SFBUpFast[pTBFast.size()-1][EtaBFast.size()-1];
     double SFCUpFast[pTBFast.size()-1][EtaBFast.size()-1];
     double SFLUpFast[pTLFast.size()-1][EtaLFast.size()-1];

     double SFBDownFast[pTBFast.size()-1][EtaBFast.size()-1];
     double SFCDownFast[pTBFast.size()-1][EtaBFast.size()-1];
     double SFLDownFast[pTLFast.size()-1][EtaLFast.size()-1];
     
     double EffLdFast[pTLFast.size()-1][EtaLFast.size()-1];
     double EffBdFast[pTBFast.size()-1][EtaBFast.size()-1];
     double EffCdFast[pTBFast.size()-1][EtaBFast.size()-1];  

     int NpTBFast=pTBFast.size()-1;
     int NpTLFast=pTLFast.size()-1;

     if(isSMS){

       for(int r = 0; r<(pTLFast.size()-1)*(EtaLFast.size()-1); ++r){
	 SFLFast[(r%NpTLFast)][int(r/NpTLFast)]=SFLvFast[r];
	 SFLUpFast[(r%NpTLFast)][int(r/NpTLFast)]=SFLvUpFast[r];
	 SFLDownFast[(r%NpTLFast)][int(r/NpTLFast)]=SFLvDownFast[r];
	 EffLdFast[(r%NpTLFast)][int(r/NpTLFast)]=EffLFast[r];
       }

       for(int t = 0; t<(pTBFast.size()-1)*(EtaBFast.size()-1); ++t){
	 SFBFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvFast[t];
	 SFCFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvFast[t];
	 SFBUpFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvUpFast[t];
	 SFCUpFast[(t%NpTBFast)][int(t/NpTBFast)]=SFCvUpFast[t];
	 SFBDownFast[(t%NpTBFast)][int(t/NpTBFast)]=SFBvDownFast[t];
	 SFCDownFast[(t%NpTBFast)][int(t/NpTBFast)]=SFCvDownFast[t];
	 EffBdFast[(t%NpTBFast)][int(t/NpTBFast)]=EffBFast[t];
	 EffCdFast[(t%NpTBFast)][int(t/NpTBFast)]=EffCFast[t];
       }

     }

     int Ne=NearestInt(_weight*puW*wLep);
     if(Ne < 1){
       Ne=1;
       EvW=_weight*puW*wLep;
     }else{
       EvW=((_weight*puW*wLep)/double(Ne));
     }

     bool isPOG=true;
     vector<int> Tags;

     if(!_isData){

       if(!isSMS && !isEffOnly){
	 for(int y=0; y<Ne; y++){
	   for(int o=0; o<JetMultiplicity; o++){
	     SetTags(JetFlavor[o], JetpT[o], JetEta[o], JetTagCSV[o], JetTagNew[o], JetTagNewBUp[o], JetTagNewBDown[o], JetTagNewLUp[o], JetTagNewLDown[o],false);
	   }
	   TagTree->Fill();
	 }
	 
	 for(int y=0; y<Ne; y++){
	   for(int o=0; o<JetMultiplicity; o++){
	     SetTags(JetFlavor[o], JetpT[o], JetEta[o], JetTagCSV[o], JetTagNew[o], JetTagNewBUp[o], JetTagNewBDown[o], JetTagNewLUp[o], JetTagNewLDown[o],true);
	   }
	   TagTreePOG->Fill();
	 }
       }else{
	 if(!isPOG){
	   for(int o=0; o<JetMultiplicity; o++){
	     
	     Tags=SetTagsSMS(JetFlavor[o], JetpT[o], JetEta[o], JetTagCSV[o], false);
	     
	     JetTagNew[o]=Tags[0];
	     JetTagNewBUp[o]=Tags[1];
	     JetTagNewBDown[o]=Tags[2];
	     JetTagNewLUp[o]=Tags[3];
	     JetTagNewLDown[o]=Tags[4];
	     JetTagNewBUpFast[o]=Tags[5];
	     JetTagNewBDownFast[o]=Tags[6];
	     JetTagNewCUpFast[o]=Tags[7];
	     JetTagNewCDownFast[o]=Tags[8];
	     JetTagNewLUpFast[o]=Tags[9];
	     JetTagNewLDownFast[o]=Tags[10];
	     JetTagNewFast[o]=Tags[11];
	     
	   }
	   
	   TagTreeSMS->Fill();
	   
	 }else if(!isEffOnly){
	   for(int o=0; o<JetMultiplicity; o++){
	     
	     Tags=SetTagsSMS(JetFlavor[o], JetpT[o], JetEta[o], JetTagCSV[o], true);
	     
	     JetTagNew[o]=Tags[0];
	     JetTagNewBUp[o]=Tags[1];
	     JetTagNewBDown[o]=Tags[2];
	     JetTagNewLUp[o]=Tags[3];
	     JetTagNewLDown[o]=Tags[4];
	     JetTagNewBUpFast[o]=Tags[5];
	     JetTagNewBDownFast[o]=Tags[6];
	     JetTagNewCUpFast[o]=Tags[7];
	     JetTagNewCDownFast[o]=Tags[8];
	     JetTagNewLUpFast[o]=Tags[9];
	     JetTagNewLDownFast[o]=Tags[10];
	     JetTagNewFast[o]=Tags[11];
	     
	   }
	   
	   TagTreeSMSPOG->Fill();
	 }
       }
     }
  }
  
  //Fill efficiency tree
  effTree->Fill();

  file->cd();

  outTree->Write();
  effTree->Write();
  if(_isData)LumiTree->Write();

  bool isPOG=true; 

  if(!isSMS && !_isData && !isEffOnly)TagTree->Write();
  if(!isSMS && !_isData && !isEffOnly)TagTreePOG->Write();
  if(isSMS && !isPOG)TagTreeSMS->Write();
  if(isSMS && isPOG)TagTreeSMSPOG->Write();
  if(isSMS)FullSMSTree->Write();
  file->Close();

}

  
vector<TLorentzVector> SUSYMultiTop::CombineJets(vector<TLorentzVector> myjets){

  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;

  int N_comb = 1;
  for(int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }
  
  double M_min = 9999999999.0;
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
	j_temp1 += myjets[count];
      } else {
	j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }    
    double M_temp = j_temp1.M2()+j_temp2.M2();
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }
  
  j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }

  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;
  
}

void SUSYMultiTop::BubbleSort(vector <float> &num)
{
  int i, j, flag = 1;    // set flag to 1 to start first pass
  float temp;             // holding variable
  int numLength = num.size(); 
  for(i = 1; (i <= numLength) && flag; i++)
    {
      flag = 0;
      for (j=0; j < (numLength -1); j++)
	{
	  if (num[j+1] > num[j])      // ascending order simply changes to <
	    { 
	      temp = num[j];             // swap elements
	      num[j] = num[j+1];
	      num[j+1] = temp;
	      flag = 1;               // indicates that a swap occurred.
	    }
	}
    }
  return;   //arrays are passed to functions by address; nothing is returned
}
