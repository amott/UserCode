 //-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef HggReducer_h
#define HggReducer_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"

#include "HggVertexing.hh"
#include "HggEGEnergyCorrector.hh"
#include "HggEnergyScale.hh"


#ifdef __MAKECINT__                                                                                                                                 
#pragma link C++ class vector<vector<float> >+;                                                                                                     
#pragma link C++ class vector<vector<unsigned short> >+;                                                                                            
#endif                                                                                                                                              

using namespace std;

struct PreSelCuts{
  float scet1;
  float scet2;
  float maxeta;
  float ecaliso_eb;
  float ecaliso_ee;
  float trkiso_eb;
  float trkiso_ee;
  float hcaliso_eb;
  float hcaliso_ee;
  float sieie_eb;
  float sieie_ee;
  float hoe;
};

enum EnergyRegressionMethod{JoshV1, JoshV2,Yong};

class HggReducer : public Vecbos{
public:
  HggReducer(TTree *tree=0); /// Class Constructor
  HggReducer(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false,int mod=-1); /// Class Constructor
  virtual ~HggReducer();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void addTrigger(string s){triggerNames.push_back(s);}
  void setCorrectionType(int s){correctionType = s;}
  void setPreselectionSet(int s){preSelSet = s;}
  void setScaleSmear(float s){applyScaleSmear = s;}
private:
  void init(); // do variable initialization 
  void setOutputBranches();
  

  //photon preselection variables
  void setupPreSelection();
  vector<PreSelCuts> preselections;
  int preSelSet;

  HggVertexing *vertexer;
  //energy correction variables
  int correctionType;
  HggEGEnergyCorrector *corrector;

  //energy smearing
  float applyScaleSmear;
  HggEnergyScale *energyScale;


  TTree * outTree;
  // define variables for the output tree:
  vector<string> triggerNames; // list of all the triggers to consider
  int * triggerBits;       // this will be an array of the trigger decision per event (for the output tree)

  // ...
  //Event info
  int lumiBlock; 
  int runNumber; 
  int evtNumber; 
  int bunchX; 
  int orbitNumber; 
  int evtTime; 

  int phyDeclared; 
  float rho; 
  float rhoEtaMax44; 

  vector<short> *pileupBunchX; 
  vector<short> *pileupNInteraction; 
  float pileupTrueNumInterations;

  static const int MAXSC= 300; 
  int nSC; 
  float ptSC[MAXSC];
  float eSC[MAXSC];
  float eRawSC[MAXSC];
  float etaSC[MAXSC];
  float phiSC[MAXSC];
  int flagSC[MAXSC];

  static const int nPhotonMAX = 100;
  int nPhoton;
  float photonsigmaIetaIeta[nPhotonMAX];
  int photonhasPixelSeed[nPhotonMAX];
  float photonenergy[nPhotonMAX];
  float photonpt[nPhotonMAX];
  float photoneta[nPhotonMAX];
  float photonphi[nPhotonMAX];
  float photonvertexx[nPhotonMAX];
  float photonvertexy[nPhotonMAX];
  float photonvertexz[nPhotonMAX];
  int photonhasConversionTracks[nPhotonMAX];
  float photonscrawEnergy[nPhotonMAX];
  float photonsceta[nPhotonMAX];
  float photonscphi[nPhotonMAX];
  float photoncaloPositionx[nPhotonMAX];
  float photoncaloPositiony[nPhotonMAX];
  float photoncaloPositionz[nPhotonMAX];
  float photonsccaloPositionx[nPhotonMAX];
  float photonsccaloPositiony[nPhotonMAX];
  float photonsccaloPositionz[nPhotonMAX];
  float photonscenergy[nPhotonMAX];
  float photone3x3[nPhotonMAX]; 
  float photone1x5[nPhotonMAX]; 
  float photone2x5[nPhotonMAX]; 
  float    photone5x5[nPhotonMAX]; 
  float   photonmaxEnergyXtal[nPhotonMAX]; 
  float   photonr9[nPhotonMAX]; 
  float photone2nd[nPhotonMAX];
  int photonInEB[nPhotonMAX];
  
  //more 
  float photonscpreshowerEnergy[nPhotonMAX]; 
  float photonscphiWidth[nPhotonMAX]; 
  float photonscetaWidth[nPhotonMAX]; 
  int photonscclusterSize[nPhotonMAX]; 
  std::vector<std::vector<float> >* photonscbclusterenergy;
  std::vector<std::vector<float> >* photonscbclusterposx;
  std::vector<std::vector<float> >* photonscbclusterposy;
  std::vector<std::vector<float> >* photonscbclusterposz;

  
  
  int photonscbcseedind[nPhotonMAX];
  int photonscbc2ind[nPhotonMAX];
  int photonscbclastind[nPhotonMAX];
  int photonscbclast2ind[nPhotonMAX];
  float photonbcseedenergy[nPhotonMAX];
  float photonbcseedeta[nPhotonMAX];
  float photonbcseedphi[nPhotonMAX];
  float photonbcseedeMax[nPhotonMAX];
  float photonbcseede2nd[nPhotonMAX];
  float photonbcseedeLeft[nPhotonMAX];
  float photonbcseedeRight[nPhotonMAX];
  float photonbcseedeTop[nPhotonMAX];
  float photonbcseedeBottom[nPhotonMAX];
  float photonbcseede3x3[nPhotonMAX];
  float photonbcseede5x5[nPhotonMAX];
  float photonbcseedCovIEtaIEta[nPhotonMAX];
  float photonbcseedCovIEtaIPhi[nPhotonMAX];
  float photonbcseedCovIPhiIPhi[nPhotonMAX];
  float photonbc2energy[nPhotonMAX];
  float photonbc2eta[nPhotonMAX];
  float photonbc2phi[nPhotonMAX];
  float photonbc2eMax[nPhotonMAX];
  float photonbc2e2nd[nPhotonMAX];
  float photonbc2eLeft[nPhotonMAX];
  float photonbc2eRight[nPhotonMAX];
  float photonbc2eTop[nPhotonMAX];
  float photonbc2eBottom[nPhotonMAX];
  float photonbc2e3x3[nPhotonMAX];
  float photonbc2e5x5[nPhotonMAX];
  float photonbc2CovIEtaIEta[nPhotonMAX];
  float photonbc2CovIEtaIPhi[nPhotonMAX];
  float photonbc2CovIPhiIPhi[nPhotonMAX];
  float photonbclastenergy[nPhotonMAX];
  float photonbclasteta[nPhotonMAX];
  float photonbclastphi[nPhotonMAX];
  float photonbclasteMax[nPhotonMAX];
  float photonbclaste2nd[nPhotonMAX];
  float photonbclasteLeft[nPhotonMAX];
  float photonbclasteRight[nPhotonMAX];
  float photonbclasteTop[nPhotonMAX];
  float photonbclasteBottom[nPhotonMAX];
  float photonbclaste3x3[nPhotonMAX];
  float photonbclaste5x5[nPhotonMAX];
  float photonbclastCovIEtaIEta[nPhotonMAX];
  float photonbclastCovIEtaIPhi[nPhotonMAX];
  float photonbclastCovIPhiIPhi[nPhotonMAX];
  float photonbclast2energy[nPhotonMAX];
  float photonbclast2eta[nPhotonMAX];
  float photonbclast2phi[nPhotonMAX];
  float photonbclast2eMax[nPhotonMAX];
  float photonbclast2e2nd[nPhotonMAX];
  float photonbclast2eLeft[nPhotonMAX];
  float photonbclast2eRight[nPhotonMAX];
  float photonbclast2eTop[nPhotonMAX];
  float photonbclast2eBottom[nPhotonMAX];
  float photonbclast2e3x3[nPhotonMAX];
  float photonbclast2e5x5[nPhotonMAX];
  float photonbclast2CovIEtaIEta[nPhotonMAX];
  float photonbclast2CovIEtaIPhi[nPhotonMAX];
  float photonbclast2CovIPhiIPhi[nPhotonMAX];
  int photonbcseedieta[nPhotonMAX];
  int photonbcseediphi[nPhotonMAX];
  float photonbcseedetacry[nPhotonMAX];
  float photonbcseedphicry[nPhotonMAX];
  int photonbc2ieta[nPhotonMAX];
  int photonbc2iphi[nPhotonMAX];
  float photonbc2etacry[nPhotonMAX];
  float photonbc2phicry[nPhotonMAX];
  
  
  int photonscnhits[nPhotonMAX]; 
  float photoneLeft[nPhotonMAX];
  float photoneRight[nPhotonMAX];
  float photoneBottom[nPhotonMAX];
  float photoneTop[nPhotonMAX];

  float photone1x3[nPhotonMAX];
  float photone3x1[nPhotonMAX];
  float photone2x2[nPhotonMAX];
  float photone3x2[nPhotonMAX];
  float photone4x4[nPhotonMAX];
  float photone2x5Right[nPhotonMAX];
  float photone2x5Left[nPhotonMAX];
  float photone2x5Top[nPhotonMAX];
  float photone2x5Bottom[nPhotonMAX];
  float photone2x5Max[nPhotonMAX];

  float photonlat[nPhotonMAX][3];
  float photonCovEtaEta[nPhotonMAX];
  float photonCovEtaPhi[nPhotonMAX];
  float photonCovPhiPhi[nPhotonMAX];
  
  float photonCovIEtaIEta[nPhotonMAX];
  float photonCovIEtaIPhi[nPhotonMAX];
  float photonCovIPhiIPhi[nPhotonMAX];
  float photonscCovIEtaIEta[nPhotonMAX];
  float photonscCovIEtaIPhi[nPhotonMAX];
  float photonscCovIPhiIPhi[nPhotonMAX];
  float photonzernike20[nPhotonMAX];
  float photonzernike42[nPhotonMAX];

  //isolation
  float    photonhadronicOverEm[nPhotonMAX]; 
  float    photonecalRecHitSumEtConeDR03[nPhotonMAX]; 
  float    photonhcalDepth1TowerSumEtConeDR03[nPhotonMAX]; 
  float    photonhcalDepth2TowerSumEtConeDR03[nPhotonMAX]; 
  float    photonhcalTowerSumEtConeDR03[nPhotonMAX]; 
  float   photontrkSumPtHollowConeDR03[nPhotonMAX]; 
  float   photontrkSumPtSolidConeDR03[nPhotonMAX]; 
  int    photonnTrkHollowConeDR03[nPhotonMAX]; 
  int   photonnTrkSolidConeDR03[nPhotonMAX]; 
  float   photonecalRecHitSumEtConeDR04[nPhotonMAX]; 
  float   photonhcalDepth1TowerSumEtConeDR04[nPhotonMAX]; 
  float   photonhcalDepth2TowerSumEtConeDR04[nPhotonMAX]; 
  float   photonhcalTowerSumEtConeDR04[nPhotonMAX]; 
  float photontrkSumPtHollowConeDR04[nPhotonMAX]; 
  float   photontrkSumPtSolidConeDR04[nPhotonMAX]; 
  int   photonnTrkHollowConeDR04[nPhotonMAX]; 
  int   photonnTrkSolidConeDR04[nPhotonMAX]; 
  
  int photonconversionsize[nPhotonMAX]; 
  float photonconversionVertexx[nPhotonMAX]; 
  float photonconversionVertexy[nPhotonMAX]; 
  float photonconversionVertexz[nPhotonMAX]; 
  float photonconversionrefittedPairMomentumx[nPhotonMAX];
  float photonconversionrefittedPairMomentumy[nPhotonMAX];
  float photonconversionrefittedPairMomentumz[nPhotonMAX];

  float photonconversionpairInvariantMass[nPhotonMAX];
  float photonconversionpairCotThetaSeparation[nPhotonMAX];
  float photonconversionEoverPrefittedTracks[nPhotonMAX];
  float photonconversionzOfPrimaryVertexFromTracks[nPhotonMAX];
  float photonconversiondistOfMinimumApproach[nPhotonMAX];
  float photonconversiondPhiTracksAtVtx[nPhotonMAX];
  float photonconversiondPhiTracksAtEcal[nPhotonMAX];
  float photonconversiondEtaTracksAtEcal[nPhotonMAX];
  int photonconversionnTracks[nPhotonMAX];
  float photonconversionMVAout[nPhotonMAX];
  int photonconversionVertexisValid[nPhotonMAX];
  float photonconversionVertexchi2[nPhotonMAX];
  float photonconversionChiSquaredProbability[nPhotonMAX];
  float photonconversion_track1_dz[nPhotonMAX];
  float photonconversion_track1_dzError[nPhotonMAX];
  int photonconversion_track1_charge[nPhotonMAX];
  float photonconversion_track1_d0[nPhotonMAX];
  float photonconversion_track1_tracksPout[nPhotonMAX];
  float photonconversion_track1_tracksPin[nPhotonMAX];
  int photonconversion_track1_algo[nPhotonMAX];
  float photonconversion_track2_dz[nPhotonMAX];
  float photonconversion_track2_dzError[nPhotonMAX];
  int photonconversion_track2_charge[nPhotonMAX];
  int photonconversion_track2_algo[nPhotonMAX];
  float photonconversion_track2_d0[nPhotonMAX];
  float photonconversion_track2_tracksPout[nPhotonMAX];
  float photonconversion_track2_tracksPin[nPhotonMAX];



  float photonseedtime[nPhotonMAX];
  float photonseedoutOfTimeChi2[nPhotonMAX]; 
  float photonseedchi2[nPhotonMAX]; 
  int photonseedrecoFlag[nPhotonMAX]; 
  int photonseedseverityLevel[nPhotonMAX]; 
  int photonfiducialFlag[nPhotonMAX]; 
  int photonscindex[nPhotonMAX]; 
  float photonswissCross[nPhotonMAX];
  int photonieta[nPhotonMAX];
  int photoniphi[nPhotonMAX];
  
  float photonE2overE9[nPhotonMAX];
  int photonmatchToallConv[nPhotonMAX];
  
  //photon-generator matching
  float photongenphtmatch[nPhotonMAX][3];
  float photongenelematch[nPhotonMAX][3];
  float photongenphtconv[nPhotonMAX][4];
  float photonPartonmatch[nPhotonMAX][6];
  float photonPartonmatchp[nPhotonMAX][3];
  
  int photonhasMatchedPromptElectron[nPhotonMAX];

  //vertex information
  int nVertex; 
  static const int MAXVX = 100;
  float vertexx[MAXVX];
  float vertexy[MAXVX];
  float vertexz[MAXVX];
  float vertexchi2[MAXVX];
  float vertexndof[MAXVX];
  float vertexnormalizedChi2[MAXVX];
  int vertextrackSize[MAXVX];
  int vertexisFake[MAXVX];
  int vertexisValid[MAXVX];


  //GENERATOR information
  static const int MAXGenSaved = 1000;
  //gen-leve phton
  int nGenPht;
  float etaGenPht[MAXGenSaved];
  float phiGenPht[MAXGenSaved];
  float ptGenPht[MAXGenSaved];
  int pidmomGenPht[MAXGenSaved];
  int pidmom2GenPht[MAXGenSaved];
  int indmom2GenPht[MAXGenSaved];
  int pidmom3GenPht[MAXGenSaved];
  int statusGenPht[MAXGenSaved];
  float vxGenPht[MAXGenSaved];
  float vyGenPht[MAXGenSaved];
  float vzGenPht[MAXGenSaved];
  
  //gen Muon 
  int nGenMu; 
  int chaGenMu[MAXGenSaved];
  float etaGenMu[MAXGenSaved];
  float phiGenMu[MAXGenSaved];
  float ptGenMu[MAXGenSaved];
  int pidmomGenMu[MAXGenSaved];
  int pidmom2GenMu[MAXGenSaved];
  int pidmom3GenMu[MAXGenSaved];
  int statusGenMu[MAXGenSaved];
  int pidGenMu[MAXGenSaved];
  float vxGenMu[MAXGenSaved];
  float vyGenMu[MAXGenSaved];
  float vzGenMu[MAXGenSaved];
  
  //gen electron
  int nGenEle; 
  int chaGenEle[MAXGenSaved];
  float etaGenEle[MAXGenSaved];
  float phiGenEle[MAXGenSaved];
  float ptGenEle[MAXGenSaved];
  int pidmomGenEle[MAXGenSaved];
  int pidmom2GenEle[MAXGenSaved];
  int pidmom3GenEle[MAXGenSaved];
  int statusGenEle[MAXGenSaved];
  int pidGenEle[MAXGenSaved];
  float vxGenEle[MAXGenSaved];
  float vyGenEle[MAXGenSaved];
  float vzGenEle[MAXGenSaved];

};
#endif


