//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed Vecbos+jets
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini 
//    CERN
//    Thiago Tomei 
//    SPRACE, Sao Paulo, Brazil
//    Ilaria Segoni
//    CERN
//-------------------------------------------------------

/// The VecbosAnalysis class is contains the analysis algorithm for the
/// study of W+jets and Z+jets.  It is a port of Ilaria's original code
/// and supersedes the original VecbosAnalysis class.

#ifndef VecbosAnalysis_h
#define VecbosAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// FINAL class. Please don't derive from this, derive from Vecbos instead.
class VecbosAnalysis : public Vecbos{
public:

  VecbosAnalysis(TTree *tree=0, string namefile="output.root"); ///< Class Constructor
  virtual ~VecbosAnalysis();     ///< Class Destructor

  /// This method implements the selection algorithm of the W+jets and Z+jets Vecbos analyses 
  void Loop();
  
protected:
  
  void bookCaloJetHistoes(std::string jetAlgoName);
  void bookProcessIDHistoes();

  /// Functions that play the same role as beginJob() and endJob() in CMSSW
  void beginJob();
  void endJob();
  
  // These functions had a (const edm::Event& iEvent) arg in CMSSW 
  void PlotsBeforeSelection();
  bool doZed();
  bool doW();
  bool doMET();
  bool doRecoJets();

  void fillHistoes();
  //  void fillTriggerHist(TH2D* hltPlot,TH1D* hltPlotProfile, unsigned int eventcount);

  //CONFIGURABLE MEMBER DATA
  
  std::string soupType_;
  
  bool checkOtherRecoTypes;
  //  std::string Zsrc_string;
  bool useGlobalWeight;
  bool UseBestCandidate;
  
  //   edm::InputTag  DaughterIso1_;
  //   edm::InputTag  DaughterIso2_;
  
  std::string    outputHistoFileName;
  bool WSelection; bool ZSelection; bool METSelection; bool JetSelection;
  
  uint32_t jetAlgos;
  bool jetExclusive;
  unsigned int jetMultiplicityCut;
  unsigned int maxCandidatesCut;
  
  double metMinCut;
  double metMaxCut;

    
  // STANDARD MEMBER DATA AND HISTOGRAMS

  std::vector<int> ZCandidates; 
  unsigned int NumberoOfZedCand;
  static bool foundZCandidate[3];
  static float TypeOfZCounter[3];
  uint32_t ZRecoFlagNumber;
  std::vector<double> ZMassRegion;

  std::vector<float> TransWMass;
  unsigned int TotalEventsProcessed;
  unsigned int EventsAfterJetSelection;
  unsigned int EventsAfterZedSelection;
  
  //  edm::TriggerNames trigNames;
  //  edm::Handle<edm::TriggerResults> hltresults;
  unsigned int numTriggers;
 
  /// Muons
  //  edm::Handle<CandDoubleAssociations> DaughterIso[2];
  //edm::Handle<CandDoubleAssociations> DaughterIso1;
  //edm::Handle<CandDoubleAssociations> DaughterIso2;
  std::vector<double> LeptonLegsPtCut;
  std::vector<double> LeptonLegsEtaCut;
  
  ////JETs 
  std::vector<std::string> src_jet;
  double jetPtCut;
  double jetEtaMinCut;
  double jetEtaMaxCut;
  std::map<std::string, std::vector<int> > highPtJets;
  
  double dauIso[2];
  
  double EventWeight; float EventWeightFloat;
  int processIDInt;
    
  // STANDARD MEMBER DATA & HISTOGRAMS
  // The tree that Ilaria was filling.

  TFile*    fOutputFile ;
  TTree*    treeForFit;
  TTree*    treeBeforeCuts;
  Float_t   WNumberofCandidates;
  Float_t   MassForFit;
  Float_t   JetMultiplicity;
  Float_t   HltBits[90];
  Float_t   processIDIntOriginal;
  Float_t   caloMET,caloMETx,caloMETy;
  Float_t   WMassTransverse[5];
  Float_t   WMuonPt[5];
  Float_t   WMuonPx[5];
  Float_t   WMuonPy[5];
  Float_t   WMuonEta[5];
  Float_t   WMuonPhi[5];
  Float_t   MuDauIsolation[5];
  Float_t   ZMuonPt[2];
  Float_t   ZMuonPx[2];
  Float_t   ZMuonPy[2];
  Float_t   ZMuonEta[2];
  Float_t   ZMuonPhi[5];
  
  
  Float_t   BestWMuonPt;
  Float_t   BestWMuonPx;
  Float_t   BestWMuonPy;
  Float_t   BestWMuonEta;
  Float_t   BestWMuonPhi;
  
  Float_t   JetsPt[15];
  Float_t   JetsEta[15];
  Float_t   JetsPhi[15];
  
  Float_t   metAftercuts;
  Float_t   metMuCorrAftercuts;
  Float_t   metXAftercuts;
  Float_t   metXMuCorrAftercuts;
  Float_t   metYAftercuts;
  Float_t   metYMuCorrAftercuts;
  
  Float_t   metBeforecuts;
  Float_t   metMuCorrBeforecuts;
  Float_t   metXBeforecuts;
  Float_t   metXMuCorrBeforecuts;
  Float_t   metYBeforecuts;
  Float_t   metYMuCorrBeforecuts;
 
  // HISTOGRAMS
  ///Before Selection

  std::map<int, TH1D*> AllEventsMET;
  std::map<int, TH1D*> AllMuonsIso1;
  std::map<int, TH1D*> AllMuonsIso2;
  std::map<int, TH2D*> AllEventsMET_vs_Iso1;
  std::map<int, TH2D*> AllEventsMET_vs_Iso2;
  TH1D*       ProcessIDBeforeCuts;  
  TH1D*       ProcessIDAfterJets;	 
  TH1D*       ProcessIDAfterZed;	 
  TH1D*       ProcessIDInZedRegion; 
  
  TH1D*       Weight_vs_ProcessID;  
  
  std::map<int, TH2D*>        TriggersBeforeSelection;  
  std::map<int, TH2D*>        TriggersAfterEWKSelection;  
  std::map<int, TH2D*>        TriggersAfterJETSelection;  
  std::map<int, TH2D*>        TriggersInZMassRegion;  
  std::map<int, TH1D*>        TriggersBeforeSelection1D;  
  std::map<int, TH1D*>        TriggersAfterEWKSelection1D;  
  std::map<int, TH1D*>        TriggersAfterJETSelection1D;  
  std::map<int, TH1D*>        TriggersInZMassRegion1D;  

  ///After Selection
  
  TH1D*       ProcessIDAfterCuts; 
  TH1D*       NumberOfEvents; 
  
  std::map<int, TH1D*>	    NumberOfWCandidates;
  std::map<int, TH1D*>      WTransvMass;  
  
  std::map<int, TH1D*>      ZCandidatesMass[5];
  std::map<int, TH1D*>	    ZCandidatesPt;
  std::map<int, TH1D*>	    NumberOfZCandidates;
  
  std::map<int, TH1D*>	    muP[3];
  std::map<int, TH1D*>	    muPt[3];
  std::map<int, TH1D*>	    muPx[3];
  std::map<int, TH1D*>	    muPy[3];
  std::map<int, TH1D*>	    muPz[3];
  std::map<int, TH1D*>	    muEta[3];
  std::map<int, TH1D*>	    muPhi[3];
  std::map<int, TH2D*>      muEtaPt[3]; 
  std::map<int, TH1D*>	    muIso[3];
  std::map<int, TH2D*>	    ZCandidatesMassvsMuonEta[3];
  std::map<int, TH1D*>	    muEtaCorrelation;
  std::map<int, TH1D*>	    muPhiCorrelation;
  
  std::map<int, TH1D*>	    muP_SB[3];
  std::map<int, TH1D*>	    muPt_SB[3];
  std::map<int, TH1D*>	    muPx_SB[3];
  std::map<int, TH1D*>	    muPy_SB[3];
  std::map<int, TH1D*>	    muPz_SB[3];
  std::map<int, TH1D*>	    muEta_SB[3];
  std::map<int, TH1D*>	    muPhi_SB[3];
  std::map<int, TH2D*>      muEtaPt_SB[3]; 
  std::map<int, TH1D*>	    muIso_SB[3];
  std::map<int, TH2D*>	    ZCandidatesMassvsMuonEta_SB[3];
  
  std::map<int, TH1D*>	    TypeOfZCandidates;
  
  std::map<int, TH1D*>	     MET;
  std::map<int, TH1D*>	     METx;
  std::map<int, TH1D*>	     METy;
  
  std::map<int, TH1D*>	     GenMET;
  std::map<int, TH1D*>	     GenMETx;
  std::map<int, TH1D*>	     GenMETy;
  
  std::map<int, TH1D*>	     GenNoNuMET;
  std::map<int, TH1D*>	     GenNoNuMETx;
  std::map<int, TH1D*>	     GenNoNuMETy;
  
  std::map<int, TH1D*>	     NuPtGen;
  std::map<int, TH1D*>	     NuPtResidual;
  
  std::map<int, TH1D*>  JetMult;     	 //std::map<std::string,TH1D*>     JetMult;
  std::map<int, TH1D*>  JetPTAll;    	 //std::map<std::string,TH1D*>     JetPTAll;
  std::map<int, TH1D*>  JetEtaAll;   	 //std::map<std::string,TH1D*>     JetEtaAll;
  std::map<int, TH1D*>  JetPhiAll;   	 //std::map<std::string,TH1D*>     JetPhiAll;
  std::map<int, TH2D*>  JetEtaVSPTAll;	 //std::map<std::string,TH2D*>     JetEtaVSPTAll;

  /// Initialize the analysis parameters to their
  /// default values. All the values can be changed by
  /// the user specifying their value in a config file,
  /// to pass as input to the ReadParameters() function.
  virtual void InitParameters(string namefile);
  
  /// Taking as input a map containing the name and the
  /// variables of the setting inputs, it initializes the
  /// analysis parameters. The analysis parameters
  /// are private members of the Vecbos class.
  virtual void AssignParameters(map<string, double>);

  // Useful algorithms
  double transverse(double x, double y); ///<Calculates sqrt(x*x+y*y);

};

bool VecbosAnalysis::foundZCandidate[3];
float VecbosAnalysis::TypeOfZCounter[3];

#endif
