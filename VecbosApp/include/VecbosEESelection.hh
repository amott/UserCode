//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed W-Z + jets
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef VecbosEESelection_h
#define VecbosEESelection_h

#include <vector>
#include <TVector3.h>
#include <TMath.h>
#include <TF2.h>
#include <TLorentzVector.h>
//#include "TMVA/Reader.h"

#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "include/Vecbos.hh"
#include "VecbosBase.hh"
#include "include/RedVecbosTree.h"
#include "include/RedVecbosVertexTree.h"
#include "include/RedVecbosMcTree.h"
#include "include/RedIsolationVtxOptimTree.hh"
#include "include/RedEleIDOptimTree.hh"
#include "include/McTruthEvent.hh"
#include "include/CutBasedSelectorEE.hh"

class VecbosEESelection : public Vecbos {
public:

  //! constructor
  VecbosEESelection(TTree *tree=0);
  //! destructor
  virtual ~VecbosEESelection();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies();
  //! do the best electron study
  void doBestElectronStudy(bool dobestelectronstudy) { m_doBestElectronStudy = dobestelectronstudy; }
  //! do the b-tag efficiency study
  void doBTagEfficiencyStudy(bool what) { m_doBTagEfficiency = what; }
  //! set the name of the skim file; if not set, do not write the file
  void setSkimOutputFile(const char *filename) { m_skimOutputFile = filename; }
  //! use the skim given in the file
  void useTheSkim(const char *filename) { m_skimFile = filename; }
  //! use a prefix for the output 
  void setPrefix(const char *prefix) { m_prefix = prefix; }
  //! set the signal type to select based on MC-truth
  void setSignal(int signal) { m_signal = signal; }

private:

  //! set the MC truth tree
  void setMcTruthTree();
  //! return the two highest pt electrons  
  std::pair<int, int> getBestGoodElePair(std::vector<int> goodElectrons);
  //!return the two best electrons with different criteria
  void getBestGoodElePairFunny(std::vector<int> goodElectrons);
  //! set the event kinematics
  void setKinematics(int theEle1, int theEle2);
  //! set the event kinematics after having chosen the jets
  void setKinematics2();
  //! set the best (di)-electron variables
  void setElectrons();
  //! compute invariant mass
  float getMee(int theEle1, int theEle2);
  
  float detaCorrections(float etaEle,float phiEle)
  {
    if (fabs(etaEle)<1.479)
      return 0.;
    
    float alignPar[3]={0.,0.,0.};
    
    if ((etaEle)>=1.479)
      {
	alignPar[0]= 0.52;
	alignPar[1]= -0.81;
	alignPar[2]= 0.81;
      }
    else if ((etaEle)<=-1.479)
      {
        alignPar[0]= -0.02;
        alignPar[1]= -0.81;
        alignPar[2]= -0.94;
      }
    
    int zIndex=etaEle>0 ? 1: -1;
    float xSCNew =(zIndex* 330. * TMath::Cos(phiEle) / TMath::SinH(etaEle)) + alignPar[0];
    float ySCNew =(zIndex* 330. * TMath::Sin(phiEle) / TMath::SinH(etaEle)) + alignPar[1];
    float zSCNew = (zIndex * 330.) + alignPar[2];
    
    //   std::cout << "Eta new:  " << TVector3(xSCNew,ySCNew,zSCNew).Eta() << "Eta old: " << etaEle << std::endl;
    
    return etaEle-TVector3(xSCNew,ySCNew,zSCNew).Eta();
    
//     TF2 f12_fit2("f12_fit2","[0]+(TMath::TanH(y)/325.)*([1]-([2]*TMath::SinH(y)*TMath::Cos(x-[3])))",-TMath::Pi(),TMath::Pi(),-10.,10.);
//     if (etaEle>1.479)
//       {
// 	f12_fit2.SetParameter(0,0.0013);
// 	f12_fit2.SetParameter(1,-0.06);
// 	f12_fit2.SetParameter(2,0.52);
// 	f12_fit2.SetParameter(3,2.17);
// 	return f12_fit2.Eval(phiEle,etaEle);
//       }
//     else if (etaEle<-1.479)
//       {
// 	f12_fit2.SetParameter(0,-0.0013);
// 	f12_fit2.SetParameter(1,-0.32);
// 	f12_fit2.SetParameter(2,0.45);
// 	f12_fit2.SetParameter(3,-1.58);
// 	return f12_fit2.Eval(phiEle,etaEle);
//       }

//     return 0.;
  }

  float dphiCorrections(float etaEle,float phiEle)
  {
    if (fabs(etaEle)<1.479)
      return 0.;
    
    float alignPar[3]={0.,0.,0.};
    
    if ((etaEle)>=1.479)
      {
	alignPar[0]=0.52;
	alignPar[1]=-0.81;
	alignPar[2]=0.81;
      }
    else if ((etaEle)<=-1.479)
      {
	alignPar[0]=-0.02;
	alignPar[1]=-0.81;
	alignPar[2]=-0.94;
      }
    
    
    int zIndex=etaEle>0 ? 1: -1;
    float xSCNew =(zIndex*330.) * TMath::Cos(phiEle) / (TMath::SinH(etaEle)) + alignPar[0];
    float ySCNew =(zIndex*330.) * TMath::Sin(phiEle) / (TMath::SinH(etaEle)) + alignPar[1];
    float zSCNew = (zIndex*330.) + alignPar[2];
    
    return phiEle-TVector3(xSCNew,ySCNew,zSCNew).Phi();

//     TF2 f12_fit2("f12_fit2","[0]+[1]*(TMath::SinH(y)/325.)*(TMath::Sin([2]-x))",-TMath::Pi(),TMath::Pi(),-10.,10.);
    
//     if (etaEle>1.479)
//       {
// 	f12_fit2.SetParameter(0,0.);
// 	f12_fit2.SetParameter(1,0.52);
// 	f12_fit2.SetParameter(2,2.17);
// 	return f12_fit2.Eval(phiEle,etaEle);
//       }
//     else if (etaEle<-1.479)
//       {
// 	f12_fit2.FixParameter(0,0.);
// 	f12_fit2.FixParameter(1,0.45);
// 	f12_fit2.FixParameter(2,-1.58);
// 	return f12_fit2.Eval(phiEle,etaEle);
//       }
//     return 0.;
  }
  
 //! returns the output of the custom cut electron ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! do the counting of b-tagged jets
  void countBTagJets(std::vector<Jet> jets, std::vector<Jet> btagjets);
  //! do the counting of B-hadrons (MC pnly)
  int countMatchingB(std::vector<Jet> jets);

#ifdef USECALOTOWERS
  //! returns the list of the calotowers 
  std::vector<CaloTower> getCaloTowers(float zpv, int type);
  //! returns the list of the calotowers for track jets removing electrons selected in the analysis 
  void getCaloTowersForPFJets(float ptMin, float ptMax, float chi2, float etaMax, float nHits, float Dxy, float dZ, float d3d, int closestPV, int theEle1, int theEle2); 
  //! make a list of SIS cone jets from calotowers sorted in Pt
  std::vector<Jet> buildSortedSISConeCaloJets();
  //! for shape variables
  vector<float> MHTphiJetForW(int ele1, std::vector<Jet>& cleanedWJets);
  float MHTphiMETForZ(int ele1, int ele2);
  float phiJetMETForW(std::vector<Jet>& cleanedWJets, std::vector<TLorentzVector>& tracksForPFJets);
#endif

  //! returns true if the W electron combined with a reco electron makes Z mass
  bool foundZCandidate(int eleIndex, std::vector<int> accEles);

  //! fill the tree for electron ID oprimization
  void FillEleIDOptimTree(int eleIndex);
  //! fill the tree for vertex and isolation optimization
  void FillIsolationVertexTree(int closestPV, int iele);

  //! selections common to all analyses
  void ConfigCommonSelections(Selection* _selection);

  //! selections for Z->ee analysis
  void ConfigZeeSelections(Selection* _selection, Counters *_counter);
  void displayZeeEfficiencies(Counters _counter);

  //! selections for W->enu analysis
  void ConfigWenuSelections(Selection* _selection, Counters *_counter);
  void displayWenuEfficiencies(Counters _counter);

  //! book the hostograms to study the best electron choice
  void bookBestCandidateHistos();


  //! DR matching between MC particle and reco electron
  bool electron_MCmatch_DeltaR(int iMc, int iEle, float deltaR);
  float deltaR_MCmatch(int iMc, int iEle);

  //! methods for B veto
  bool isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits);
  vector<float> jetBTagVariables(Jet jet);
  float jetDxy(Jet jet);
  float jetDsz(Jet jet);
  float invMassTracks(Jet jet);
  int tracksMultiplicity(Jet jet);
  float fractionalEnergy(Jet jet);
  void calcEventBVetoVariables(std::vector<Jet> jets, std::vector<BTagJet> btags);
  
  //! get the supercluster 4-vector of one electron
  TLorentzVector getSC4Vector(int eleindex);

  //! extra functions for shape variables
  void LoadCorrections();
  double ParCorr(double, double);
  double PerpCorr(double, double);

  float ParConst[50][50];
  float PerpConst[50][50];
  float Mubin[50];
  float Parbin[50];
  float Perpbin[50];
  float MubinCenter[50];
  float ParbinCenter[50];
  float PerpbinCenter[50];
  static const int NMET = 50;
  static const int NMu = 50;

  //! the members
  int theEle1_, theEle2_;
  TLorentzVector p4Ele1_, p4Ele2_;
  TVector3 pT3Ele1_, pT3Ele2_, p3Met_, U_, p3W_, pt3Z_, p3Z_, sumCalo_;
  TVector3 p3CaloMet_, p3TCMet_, p3PFMet_;
  float WCalomT_, WTCmT_, WPFmT_, ZMetMT_, ZJetM_;
  float U_par_, U_perp_, WmT_, WmTFromZUncorr_, WmTFromZCorr_, tempMETuncorr_, tempMETcorr_, mee_;
  int closestPV;

  //! get the mass of the heaviest di-jet pair
  float GetDiJetHeaviestMass(std::vector<Jet> cleanedJets);
  //! count muons in the event
  int CountMuons(float etaMax, float ptMin);

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;
  CutBasedEleIDSelector EgammaTightID_NC, EgammaLooseID_NC;
  ElectronLikelihood *LH;
  //! be verbose during runtime
  bool _verbose;
  //! to compute final quantities
  int theBestZee; 
  //! vectors to store indices of best candidates
  std::vector<int> goodElectrons;

  //! vector of isolation variables for reconstructed electrons
  std::vector<float> trackerIsolation, ecalIsolation, hcalIsolation;
  //! input variables to isolation Fisher 
  float trackerIsolRel, ecalIsolRel, hcalIsolRel;
  //! vector of combined isolation for reconstructed electrons
  std::vector<float> combinedIsolation;
  //  TMVA::Reader *reader_EB, *reader_EE;

  //! vectors of cleaned jets
  std::vector<Jet> goodWJets_[2], goodZJets_[2], goodWPFJets_[2], goodZPFJets_[2];
  //! vectors of cleaned jets for b-veto
  std::vector<Jet> m_cleanedBvetoWJets, m_cleanedBvetoZJets, m_cleanedBvetoWPFJets, m_cleanedBvetoZPFJets;
  std::vector<BTagJet> m_cleanedBvetoWBTagJets, m_cleanedBvetoZBTagJets, m_cleanedBvetoWBTagPFJets, m_cleanedBvetoZBTagPFJets;
  //! EVT b-tagging
  BTagJet m_maxBTagEvt;

  //! the CMS calorimeters
  enum { ecal=0, hcal=1 };

  //! the processes
  enum { wjets=0, zjets=1, wother=2, zother=3, all=4 };
  int m_signal;

  //! the calotowers of the event
  std::vector<CaloTower> _theCaloTowers;
  std::vector<CaloTower> _theCaloTowersForShape;
  std::vector<CaloTower> _theCaloTowersForPFJets;
  
  //! the SISCone jets of the event
  std::vector<Jet> _theSISConeCaloJets;
  std::vector<Jet> _theSISConePFJets;

  //! the antiKT jets of the event
  std::vector<Jet> _theAK5CaloJets;
  std::vector<Jet> _theAK5PFJets;
  std::vector<Jet> _theAK5GenJets;

  //! the btag collection
  std::vector<BTagJet> _theSISConeBTagCaloJets;
  std::vector<BTagJet> _theAK5BTagCaloJets;
  std::vector<BTagJet> _theAK5BTagPFJets;

  //! to compute shape variables
  CoolTools *cTools[2];
  TLorentzVector MHT;
 
  //! best electron based on different criteria
  std::pair<int,int> _bestPairByPt, _bestPairBySCenergy, _bestPairByTrackerIsolation, _bestPairByHcalIsolation, _bestPairByElectronIdLH;

  //!counters and selections
  Selection *_commonSel,       *_zeeSel,         *_wenuSel; 
  Counters _zeeCounter,      _wenuCounter;    
  Counters _zeeCounterZjets, _wenuCounterZjets; 
  Counters _zeeCounterWjets, _wenuCounterWjets; 
  Counters _zeeCounterTTbar, _wenuCounterTTbar; 
  Counters *_zeePCounter,      *_wenuPCounter;    
  Counters *_zeePCounterZjets, *_wenuPCounterZjets; 
  Counters *_zeePCounterWjets, *_wenuPCounterWjets; 
  Counters *_zeePCounterTTbar, *_wenuPCounterTTbar; 

  //! number of reco jets (with reco electrons removal)
  int howManyWJets[2], howManyZJets[2];
  int howManyForEffWJets[2], howManyForEffZJets[2];
  int howManyWPFJets[2], howManyZPFJets[2];
  int howManyForEffWPFJets[2], howManyForEffZPFJets[2];
  
  //! number of b-tag jets with different thresholds
  int nTagWPFJets[5];

  //! leading jets pt for W studies
  float leadingCaloJetPtW, leadingPFJetPtW;
  float emFracCaloJetW, emFracPFJetW;

  //! leading jets pt for Z studies
  float leadingCaloJetPtZ, leadingPFJetPtZ;
  float emFracCaloJetZ, emFracPFJetZ;

  //! di-jet mass
  float diPFJetMassW[2], diPFJetMassZ[2];

  //! di-jet pT
  float diPFJetPtW[2], diPFJetPtZ[2];

  //! delta eta between the jets
  float deltaEtaPFJetW[2], deltaEtaPFJetZ[2];

  //! delta phi between leading jet and met
  float deltaPhiJetMetW[2], deltaPhiJetMetZ[2];

  //! counters by jet multiplicity for signal efficiency (with MC truth electrons removal)
  //! array: [inclusive multiplicity][threshold:0=15GeV/1=30GeV]
  Counters m_wenuJetbinCounter[2][6], m_zeeJetbinCounter[2][6];
  Counters m_wenuRecoJetbinCounter[2][6], m_zeeRecoJetbinCounter[2][6];
  Counters m_wenuPFJetbinCounter[2][6], m_zeePFJetbinCounter[2][6];
  Counters m_wenuRecoPFJetbinCounter[2][6], m_zeeRecoPFJetbinCounter[2][6];
  
  //! reduced trees
  RedVecbosTree   *myOutTree_Zee;
  RedVecbosTree   *myOutTree_Wenu;
  RedVecbosMcTree *myOutMcTree_Zee;
  RedVecbosMcTree *myOutMcTree_Wenu;
  RedVecbosVertexTree *myOutVertexTree_FullWenu;
  RedVecbosVertexTree *myOutVertexTree_AccWenu;
  RedVecbosVertexTree *myOutVertexTree_IdWenu;
  RedIsolationVtxOptimTree *myOutIsolationVtxOptimTree_Wenu;
  RedEleIDOptimTree *myOutEleIDOptimTree_Wenu;

  //! name of rootfile with dataset
  std::string _datasetName;

  //! do best elctron study or b-tag efficiency study
  bool m_doBestElectronStudy, m_doBTagEfficiency;
  //! histograms to study the purity of the best electron selection
  vector<TH1F*> Gene_eta, H_etaBestByPt, H_etaBestBySCenergy, H_etaBestByTkIsol, H_etaBestByHcalIsol, H_etaBestByLH;
  vector<TH1F*> Gene_pt, H_ptBestByPt;
  TH1F *H_etaResolution, *H_phiResolution, *H_deltaR, *H_etaMcNotMatched, *H_phiMcNotMatched, *H_ptNotMatched, 
    *H_nRecoNotMatched, *H_nIdNotMatched, *H_nIsolNotMatched;
  //! skim output file
  const char *m_skimOutputFile, *m_skimFile;
  //! WToENuDecay = 1 if W->e nu promptly (exclude other W decays, as W->tau nu)
  int WToENuDecay;
  int ZToEEDecay;
  int WinAccept, ZinAccept;
  //! pthat for the case of photon+jet
  float photonj_pthat;
  //! the prefix to add to the output files
  const char *m_prefix;
  char namefile[200];
  //! shower shape variables 
  vector<TLorentzVector> CTinput;

  //! counters to compute number of jets at the beginning of the selection
  double m_wInitJets[6];
  double m_zInitJets[6];
  double m_wEndJets[6];
  double m_zEndJets[6];
  int nwjets_mc[2], nzjets_mc[2], nwpfjets_mc[2], nzpfjets_mc[2], nwgenjets_mc[2], nzgenjets_mc[2];

  //! best (di)-electrons variables
  int recoflag_[2];
  float pt_[2], eta_[2], phi_[2];
  int classification_[2], nbrems_[2];
  float deta_[2], dphi_[2], hoe_[2], see_[2], spp_[2], eop_[2], fbrem_[2], lh_[2], e9esc_[2],e25esc_[2];
  float trackerIso_[2], hcalIso_[2], ecalJIso_[2], ecalGTIso_[2], combinedIso_[2];
  float dist_[2], dcot_[2];
  float esc_[2], eseed_[2], pin_[2], pout_[2];
  int charge_[2], missHits_[2];

  //! to estimate b-tag efficiency and mistag rate
  int numBJet_reco, numNoBJet_reco, numBJet_tag, numNoBJet_tag;

  //! jes uncertainty (+1/-1/0)
  int jes_;
  
  //! to check MC acceptance
  float ptGen1, ptGen2, etaGen1, etaGen2, meeGen;

  int isData_;
  McTruthEvent mcevent;

};
#endif
