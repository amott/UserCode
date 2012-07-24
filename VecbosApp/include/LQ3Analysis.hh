//---------------------------------------------------------------------------
// Description:
//    Class for exporting light-weight tree for LQ studies
//    And also make histograms
// Authors:
//    Y. Chen
//---------------------------------------------------------------------------
#ifndef LQ3Analysis_h
#define LQ3Analysis_h
//---------------------------------------------------------------------------
#include "TTree.h"
//---------------------------------------------------------------------------
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/TriggerMask.hh"
#include "LQ3Helper.hh"
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
class LQ3Analysis;
struct MuonCandidate;
struct ElectronCandidate;
class LQ3OutputRecord;
//---------------------------------------------------------------------------
class LQ3Analysis : public Vecbos
{
public:
   LQ3Analysis(TTree *tree = 0); /// Class Constructor
   LQ3Analysis(TTree *tree = 0, bool isData = false, string JSONFile = ""); /// Class Constructor
   virtual ~LQ3Analysis();     /// Class Destructor
   void Loop(string outFileName, int start, int stop);
   void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
   void SetConditions(TTree* treeCond);

private:
   double SumPt(int iMu, int iZ);
   int BestPV();
   double pTMuon(int imu);
   double DeltaPhi_PiHalf(double phi1, double phi2);

   double weight;
   vector<int> m_requiredTriggers;

   bool _isData;
   bool _isNewMC;
   bool _goodRunLS;
   int _selectZMuMu;   // if signal MC, go into MC information and find Z and its daughters
                       // -1: don't care, 1: require one, 0: no ZMuMu allowed
   TTree* _treeCond;

private:
   map<string, TH1D *> GenerateQM1DHistograms();
   map<string, TH2D *> GenerateQM2DHistograms();
   void WriteQMHistograms(map<string, TH1D *> Histograms);
   void WriteQMHistograms(map<string, TH2D *> Histograms);
   void DeleteQMHistograms(map<string, TH1D *> Histograms);
   void DeleteQMHistograms(map<string, TH2D *> Histograms);

private:
   void QMFillBasicInformation(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillCaloJet(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillPFJet(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillCaloMET(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillPFMET(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillCaloBTag(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillPFBTag(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms);
   void QMFillCaloHemisphere(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
      vector<FourVector> &Groups);
   void QMFillPFHemisphere(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
      vector<FourVector> &Groups);
   void QMFillCaloRazor(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
      vector<FourVector> &Groups, int TCHELCount, int SSVHEMCount, int TCHETCount);
   void QMFillPFRazor(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
      vector<FourVector> &Groups, int TCHELCount, int SSVHEMCount, int TCHETCount);

private:
   bool CheckCaloJetID(int Index);
   bool CheckPFJetID(int Index);
   bool CheckVertex(int Index);
   bool CheckHLTBit(string PathName);

   vector<MuonCandidate> MakeMuonCandidates();
   vector<ElectronCandidate> MakeElectronCandidates();

   void FillTopMuons(vector<MuonCandidate> &Candidates, LQ3OutputRecord &M);
   void FillTopElectrons(vector<ElectronCandidate> &Candidates, LQ3OutputRecord &M);
   void FillCaloJet(LQ3OutputRecord &M);
   void FillPFJet(LQ3OutputRecord &M);
   void FillPV(LQ3OutputRecord &M);
   void FillMET(LQ3OutputRecord &M);
};
//---------------------------------------------------------------------------
struct MuonCandidate
{
   bool IsGlobal;
   bool IsPromptTight;
   bool IsTracker;
   int PixelHit;
   int StripHit;
   double Chi2;
   int ValidMuonHit;
   int MuonStations;
   double Dxy;
   double Isolation;
   double CombinedIsolation;
   double Phi;
   double Eta;
   double PT;
   double P;
   double Px;
   double Py;
   double Pz;
   int Charge;

   bool PassMuonID;
   bool PassMuonTight;
   bool PassMuonLoose;
};
//---------------------------------------------------------------------------
struct ElectronCandidate
{
   double SuperClusterEnergy;
   double SuperClusterEta;
   double SuperClusterEtaWidth;
   double SuperClusterPhi;
   double SuperClusterPhiWidth;
   double SuperClusterHcalTowerSumEt03;
   double SuperClusterSigmaIEtaIEta;
   int MissingHits;
   double ConversionDistance;
   double ConversionDeltaCotTheta;
   double DeltaEtaAtCalo;
   double DeltaPhiAtCalo;
   double HcalIsolation;
   double TrackIsolation;
   double EcalIsolation;
   double HOverE;
   double Phi;
   double Eta;
   double PT;
   double P;
   double Px;
   double Py;
   double Pz;
   int Charge;

   bool PassWP80;
   bool PassWP85;
   bool PassWP90;
   bool PassWP95;
};
//---------------------------------------------------------------------------
class LQ3OutputRecord
{
public:
   double RunNumber;
   double EventNumber;
   double BunchCrossing;
   double LumiSection;
   double Orbit;
   int PileUp;

   int CaloJetCount;
   int CaloJetCount30;
   int CaloJetCount50;
   int CaloJetCount100;
   double CaloJetE[100];
   double CaloJetPT[100];
   double CaloJetEta[100];
   double CaloJetPhi[100];
   double CaloJetCSVTag[100];
   double CaloJetCSVMTag[100];
   double CaloJetTCHPTag[100];
   double CaloJetTCHETag[100];
   double CaloJetProbabilityTag[100];
   double CaloJetBProbabilityTag[100];
   double CaloJetSSVHETag[100];
   double CaloJetSSVHPTag[100];
   
   int PFJetCount;
   int PFJetCount30;
   int PFJetCount50;
   int PFJetCount100;
   double PFJetE[100];
   double PFJetPT[100];
   double PFJetEta[100];
   double PFJetPhi[100];
   double PFJetCSVTag[100];
   double PFJetCSVMTag[100];
   double PFJetTCHPTag[100];
   double PFJetTCHETag[100];
   double PFJetProbabilityTag[100];
   double PFJetBProbabilityTag[100];
   double PFJetSSVHETag[100];
   double PFJetSSVHPTag[100];

   int PrimaryVertexCount;
   double PrimaryVertexMaxSumPT;

   double CaloMET[2];
   double PFMET[2];

   int GoodMuonCount;
   int GoodMuonCountTight;
   int GoodMuonCountLoose;
   MuonCandidate Muons[10];

   int GoodElectronCount;
   int GoodElectronCount80;
   int GoodElectronCount85;
   int GoodElectronCount90;
   int GoodElectronCount95;
   ElectronCandidate Electrons[10];

   bool PassHLT;
   bool PassNoiseFilter;
   bool PassCaloJetID;
   bool PassCaloJet60;
   bool PassPFJetID;
   bool PassPFJet60;

   bool PassDiJetAve30;
   bool PassDiJetAve60;
   bool PassDiJetAve80;
   bool PassDiJetAve110;
   bool PassDiJetAve150;
   bool PassDiJetAve190;
   bool PassDiJetAve240;
   bool PassDiJetAve300;
   bool PassDiJetAve370;

   bool PassDiJetAve15U;
   bool PassDiJetAve30U;
   bool PassDiJetAve50U;
   bool PassDiJetAve70U;
   bool PassDiJetAve100U;
   bool PassDiJetAve140U;
   bool PassDiJetAve180U;
   bool PassDiJetAve300U;

   bool PassSingleMu3;
   bool PassSingleMu5;
   bool PassSingleMu8;
   bool PassSingleMu12;
   bool PassSingleMu15;
   bool PassSingleMu20;
   bool PassSingleMu24;
   bool PassSingleMu30;
   bool PassSingleMu40;
   bool PassSingleMu100;

   bool PassR020MR500;
   bool PassR020MR550;
   bool PassR023MR550;
   bool PassR025MR400;
   bool PassR025MR450;
   bool PassR029MR450;
   bool PassR033MR300;
   bool PassR033MR350;
   bool PassR036MR350;
   bool PassR038MR200;
   bool PassR038MR250;
   bool PassR042MR250;

   bool PassEleL10R020MR200;
   bool PassEleL10R025MR200;
   bool PassEleL10R029MR200;
   bool PassEleT10R020MR200;

   bool PassMu8R020MR200;
   bool PassMu8R025MR200;
   bool PassMu8R029MR200;

   bool PassMu17Mu8;
   bool PassMu13Mu8;

   bool PassMu8EleL17;
   bool PassMu17EleL8;

public:
   void Clear()
   {
      RunNumber = 0;
      EventNumber = 0;
      LumiSection = 0;
      BunchCrossing = 0;
      Orbit = 0;
      PileUp = 0;

      CaloJetCount = 0;
      CaloJetCount30 = 0;
      CaloJetCount50 = 0;
      CaloJetCount100 = 0;
      for(int i = 0; i < 50; i++)
      {
         CaloJetE[i] = -1;
         CaloJetPT[i] = -1;
         CaloJetEta[i] = -1;
         CaloJetPhi[i] = -1;
         CaloJetCSVTag[i] = -1;
         CaloJetCSVMTag[i] = -1;
         CaloJetTCHPTag[i] = -1;
         CaloJetTCHETag[i] = -1;
      }

      PFJetCount = 0;
      PFJetCount30 = 0;
      PFJetCount50 = 0;
      PFJetCount100 = 0;
      for(int i = 0; i < 50; i++)
      {
         PFJetE[i] = -1;
         PFJetPT[i] = -1;
         PFJetEta[i] = -1;
         PFJetPhi[i] = -1;
         PFJetCSVTag[i] = -1;
         PFJetCSVMTag[i] = -1;
         PFJetTCHPTag[i] = -1;
         PFJetTCHETag[i] = -1;
      }

      PrimaryVertexCount = 0;
      PrimaryVertexMaxSumPT = 0;
   }
   
   void MakeBranches(TTree *Tree)
   {
      // General
      Tree->Branch("RunNumber", &RunNumber, "RunNumber/D");
      Tree->Branch("EventNumber", &EventNumber, "EventNumber/D");
      Tree->Branch("BunchCrossing", &BunchCrossing, "BunchCrossing/D");
      Tree->Branch("LumiSection", &LumiSection, "LumiSection/D");
      Tree->Branch("Orbit", &Orbit, "Orbit/D");
      Tree->Branch("PileUp", &PileUp, "PileUp/I");

      // calo jets
      Tree->Branch("CaloJetCount", &CaloJetCount, "CaloJetCount/I");
      Tree->Branch("CaloJetCount30", &CaloJetCount30, "CaloJetCount30/I");
      Tree->Branch("CaloJetCount50", &CaloJetCount50, "CaloJetCount50/I");
      Tree->Branch("CaloJetCount100", &CaloJetCount100, "CaloJetCount100/I");
      Tree->Branch("CaloJetE", CaloJetE, "CaloJetE[100]/D");
      Tree->Branch("CaloJetPT", CaloJetPT, "CaloJetPT[100]/D");
      Tree->Branch("CaloJetEta", CaloJetEta, "CaloJetEta[100]/D");
      Tree->Branch("CaloJetPhi", CaloJetPhi, "CaloJetPhi[100]/D");
      Tree->Branch("CaloJetCSVTag", CaloJetCSVTag, "CaloJetCSVTag[100]/D");
      Tree->Branch("CaloJetCSVMTag", CaloJetCSVMTag, "CaloJetCSVMTag[100]/D");
      Tree->Branch("CaloJetSSVHETag", CaloJetSSVHETag, "CaloJetSSVHETag[100]/D");
      Tree->Branch("CaloJetSSVHPTag", CaloJetSSVHPTag, "CaloJetSSVHPTag[100]/D");
      Tree->Branch("CaloJetTCHPTag", CaloJetTCHPTag, "CaloJetTCHPTag[100]/D");
      Tree->Branch("CaloJetTCHETag", CaloJetTCHETag, "CaloJetTCHETag[100]/D");

      // PF jets
      Tree->Branch("PFJetCount", &PFJetCount, "PFJetCount/I");
      Tree->Branch("PFJetCount30", &PFJetCount30, "PFJetCount30/I");
      Tree->Branch("PFJetCount50", &PFJetCount50, "PFJetCount50/I");
      Tree->Branch("PFJetCount100", &PFJetCount100, "PFJetCount100/I");
      Tree->Branch("PFJetE", PFJetE, "PFJetE[100]/D");
      Tree->Branch("PFJetPT", PFJetPT, "PFJetPT[100]/D");
      Tree->Branch("PFJetEta", PFJetEta, "PFJetEta[100]/D");
      Tree->Branch("PFJetPhi", PFJetPhi, "PFJetPhi[100]/D");
      Tree->Branch("PFJetCSVTag", PFJetCSVTag, "PFJetCSVTag[100]/D");
      Tree->Branch("PFJetCSVMTag", PFJetCSVMTag, "PFJetCSVMTag[100]/D");
      Tree->Branch("PFJetSSVHETag", PFJetSSVHETag, "PFJetSSVHETag[100]/D");
      Tree->Branch("PFJetSSVHPMTag", PFJetSSVHPTag, "PFJetSSVHPTag[100]/D");
      Tree->Branch("PFJetTCHPTag", PFJetTCHPTag, "PFJetTCHPTag[100]/D");
      Tree->Branch("PFJetTCHETag", PFJetTCHETag, "PFJetTCHETag[100]/D");

      // Primary vertex
      Tree->Branch("PrimaryVertexCount", &PrimaryVertexCount,
            "PrimaryVertexCount/I");
      Tree->Branch("PrimaryVertexMaxSumPT", &PrimaryVertexMaxSumPT,
            "PrimaryVertexMaxSumPT/D");

      // Missing energy
      Tree->Branch("CaloMET", CaloMET, "CaloMET[2]/D");
      Tree->Branch("PFMET", PFMET, "PFMET[2]/D");

      // Muons
      Tree->Branch("GoodMuonCount", &GoodMuonCount, "GoodMuonCount/I");
      Tree->Branch("GoodMuonCountTight", &GoodMuonCountTight, "GoodMuonCountTight/I");
      Tree->Branch("GoodMuonCountLoose", &GoodMuonCountLoose, "GoodMuonCountLoose/I");

      for(int i = 0; i < 10; i++)
      {
         Tree->Branch(Form("Muon%dPT", i + 1), &Muons[i].PT, Form("Muon%dPT/D", i + 1));
         Tree->Branch(Form("Muon%dEta", i + 1), &Muons[i].Eta, Form("Muon%dEta/D", i + 1));
         Tree->Branch(Form("Muon%dPhi", i + 1), &Muons[i].Phi, Form("Muon%dPhi/D", i + 1));
         Tree->Branch(Form("Muon%dPassTight", i + 1), &Muons[i].PassMuonTight, Form("Muon%dPassTight/O", i + 1));
         Tree->Branch(Form("Muon%dPassLoose", i + 1), &Muons[i].PassMuonLoose, Form("Muon%dPassLoose/O", i + 1));
      }

      // Electrons
      Tree->Branch("GoodElectronCount", &GoodElectronCount, "GoodElectronCount/I");
      Tree->Branch("GoodElectronCount80", &GoodElectronCount80, "GoodElectronCount80/I");
      Tree->Branch("GoodElectronCount85", &GoodElectronCount85, "GoodElectronCount85/I");
      Tree->Branch("GoodElectronCount90", &GoodElectronCount90, "GoodElectronCount90/I");
      Tree->Branch("GoodElectronCount95", &GoodElectronCount95, "GoodElectronCount95/I");

      for(int i = 0; i < 10; i++)
      {
         Tree->Branch(Form("Electron%dPT", i + 1), &Electrons[i].PT, Form("Electron%dPT/D", i + 1));
         Tree->Branch(Form("Electron%dEta", i + 1), &Electrons[i].Eta, Form("Electron%dEta/D", i + 1));
         Tree->Branch(Form("Electron%dPhi", i + 1), &Electrons[i].Phi, Form("Electron%dPhi/D", i + 1));
         Tree->Branch(Form("Electron%dPass95", i + 1), &Electrons[i].PassWP95, Form("Electron%dPass95/O", i + 1));
         Tree->Branch(Form("Electron%dPass90", i + 1), &Electrons[i].PassWP90, Form("Electron%dPass90/O", i + 1));
         Tree->Branch(Form("Electron%dPass85", i + 1), &Electrons[i].PassWP85, Form("Electron%dPass85/O", i + 1));
         Tree->Branch(Form("Electron%dPass80", i + 1), &Electrons[i].PassWP80, Form("Electron%dPass80/O", i + 1));
      }

      // baseline event selection
      Tree->Branch("PassHLT", &PassHLT, "PassHLT/O");
      Tree->Branch("PassNoiseFilter", &PassNoiseFilter, "PassNoiseFilter/O");
      Tree->Branch("PassCaloJetID", &PassCaloJetID, "PassCaloJetID/O");
      Tree->Branch("PassCaloJet60", &PassCaloJet60, "PassCaloJet60/O");
      Tree->Branch("PassPFJetID", &PassPFJetID, "PassPFJetID/O");
      Tree->Branch("PassPFJet60", &PassPFJet60, "PassPFJet60/O");

      // auxiliary HLT bits
      Tree->Branch("PassDiJetAve30", &PassDiJetAve30, "PassDiJetAve30/O");
      Tree->Branch("PassDiJetAve60", &PassDiJetAve60, "PassDiJetAve60/O");
      Tree->Branch("PassDiJetAve80", &PassDiJetAve80, "PassDiJetAve80/O");
      Tree->Branch("PassDiJetAve110", &PassDiJetAve110, "PassDiJetAve110/O");
      Tree->Branch("PassDiJetAve150", &PassDiJetAve150, "PassDiJetAve150/O");
      Tree->Branch("PassDiJetAve190", &PassDiJetAve190, "PassDiJetAve190/O");
      Tree->Branch("PassDiJetAve240", &PassDiJetAve240, "PassDiJetAve240/O");
      Tree->Branch("PassDiJetAve300", &PassDiJetAve300, "PassDiJetAve300/O");
      Tree->Branch("PassDiJetAve370", &PassDiJetAve370, "PassDiJetAve370/O");

      Tree->Branch("PassDiJetAve15U", &PassDiJetAve15U, "PassDiJetAve15U/O");
      Tree->Branch("PassDiJetAve30U", &PassDiJetAve30U, "PassDiJetAve30U/O");
      Tree->Branch("PassDiJetAve50U", &PassDiJetAve50U, "PassDiJetAve50U/O");
      Tree->Branch("PassDiJetAve70U", &PassDiJetAve70U, "PassDiJetAve70U/O");
      Tree->Branch("PassDiJetAve100U", &PassDiJetAve100U, "PassDiJetAve100U/O");
      Tree->Branch("PassDiJetAve140U", &PassDiJetAve140U, "PassDiJetAve140U/O");
      Tree->Branch("PassDiJetAve180U", &PassDiJetAve180U, "PassDiJetAve180U/O");
      Tree->Branch("PassDiJetAve300U", &PassDiJetAve300U, "PassDiJetAve300U/O");

      Tree->Branch("PassSingleMu3", &PassSingleMu3, "PassSingleMu3/O");
      Tree->Branch("PassSingleMu5", &PassSingleMu5, "PassSingleMu5/O");
      Tree->Branch("PassSingleMu8", &PassSingleMu8, "PassSingleMu8/O");
      Tree->Branch("PassSingleMu12", &PassSingleMu12, "PassSingleMu12/O");
      Tree->Branch("PassSingleMu15", &PassSingleMu15, "PassSingleMu15/O");
      Tree->Branch("PassSingleMu20", &PassSingleMu20, "PassSingleMu20/O");
      Tree->Branch("PassSingleMu24", &PassSingleMu24, "PassSingleMu24/O");
      Tree->Branch("PassSingleMu30", &PassSingleMu30, "PassSingleMu30/O");
      Tree->Branch("PassSingleMu40", &PassSingleMu40, "PassSingleMu40/O");
      Tree->Branch("PassSingleMu100", &PassSingleMu100, "PassSingleMu100/O");

      Tree->Branch("PassR020MR500", &PassR020MR500, "PassR020MR500/O");
      Tree->Branch("PassR020MR550", &PassR020MR550, "PassR020MR550/O");
      Tree->Branch("PassR023MR550", &PassR023MR550, "PassR023MR550/O");
      Tree->Branch("PassR025MR400", &PassR025MR400, "PassR025MR400/O");
      Tree->Branch("PassR025MR450", &PassR025MR450, "PassR025MR450/O");
      Tree->Branch("PassR029MR450", &PassR029MR450, "PassR029MR450/O");
      Tree->Branch("PassR033MR300", &PassR033MR300, "PassR033MR300/O");
      Tree->Branch("PassR033MR350", &PassR033MR350, "PassR033MR350/O");
      Tree->Branch("PassR036MR350", &PassR036MR350, "PassR036MR350/O");
      Tree->Branch("PassR038MR200", &PassR038MR200, "PassR038MR200/O");
      Tree->Branch("PassR038MR250", &PassR038MR250, "PassR038MR250/O");
      Tree->Branch("PassR042MR250", &PassR042MR250, "PassR042MR250/O");

      Tree->Branch("PassEleL10R020MR200", &PassEleL10R020MR200, "PassEleL10R020MR200/O");
      Tree->Branch("PassEleL10R025MR200", &PassEleL10R025MR200, "PassEleL10R025MR200/O");
      Tree->Branch("PassEleL10R029MR200", &PassEleL10R029MR200, "PassEleL10R029MR200/O");
      Tree->Branch("PassEleT10R020MR200", &PassEleT10R020MR200, "PassEleT10R020MR200/O");

      Tree->Branch("PassMu8R020MR200", &PassMu8R020MR200, "PassMu8R020MR200/O");
      Tree->Branch("PassMu8R025MR200", &PassMu8R025MR200, "PassMu8R025MR200/O");
      Tree->Branch("PassMu8R029MR200", &PassMu8R029MR200, "PassMu8R029MR200/O");

      Tree->Branch("PassMu17Mu8", &PassMu17Mu8, "PassMu17Mu8/O");
      Tree->Branch("PassMu13Mu8", &PassMu13Mu8, "PassMu13Mu8/O");

      Tree->Branch("PassMu8EleL17", &PassMu8EleL17, "PassMu8EleL17/O");
      Tree->Branch("PassMu17EleL8", &PassMu17EleL8, "PassMu17EleL8/O");
   }
};
//---------------------------------------------------------------------------
#endif

