//---------------------------------------------------------------------------
// Description:
//    Class for exporting light-weight tree for LQ studies
//    And also make histograms
// Authors:
//    Y. Chen
//---------------------------------------------------------------------------
#ifndef LQ3NewVariableBrainstorm_h
#define LQ3NewVariableBrainstorm_h
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
class LQ3NewVariableBrainstorm;
struct MuonCandidate;
struct ElectronCandidate;
class LQ3NewVariabeOutputRecord;
//---------------------------------------------------------------------------
class LQ3NewVariableBrainstorm : public Vecbos
{
public:
   LQ3NewVariableBrainstorm(TTree *tree = 0); /// Class Constructor
   LQ3NewVariableBrainstorm(TTree *tree = 0, bool isData = false, string JSONFile = ""); /// Class Constructor
   virtual ~LQ3NewVariableBrainstorm();     /// Class Destructor
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
   bool CheckCaloJetID(int Index);
   bool CheckPFJetID(int Index);
   bool CheckVertex(int Index);
   bool CheckHLTBit(string PathName);

   vector<MuonCandidate> MakeMuonCandidates();
   vector<ElectronCandidate> MakeElectronCandidates();
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
   double CombinedIsolation;
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
class LQ3NewVariableOutputRecord
{
public:
   double MR0, R0;   // this is the original MR/R
   double MR1, R1;
   double MR2, R2;
   double MR3, R3;
   double MR4, R4;
   double MR5, R5;
   double MR6, R6;
   double MR7, R7;
   double MR8, R8;
   double MR9, R9;
   double MR11, R11a, R11b, R11c, R11d, R11e, R11f, R11g;
   int WToENuCount, WToMuNuCount, WToTauNuCount;
public:
   void MakeBranches(TTree *Tree)
   {
      Tree->Branch("MR0", &MR0, "MR0/D");
      Tree->Branch("R0", &R0, "R0/D");
      Tree->Branch("MR1", &MR1, "MR1/D");
      Tree->Branch("MR2", &MR2, "MR2/D");
      Tree->Branch("MR3", &MR3, "MR3/D");
      Tree->Branch("MR4", &MR4, "MR4/D");
      Tree->Branch("MR5", &MR5, "MR5/D");
      Tree->Branch("MR6", &MR6, "MR6/D");
      Tree->Branch("MR7", &MR7, "MR7/D");
      Tree->Branch("MR8", &MR8, "MR8/D");
      Tree->Branch("MR9", &MR9, "MR9/D");
      Tree->Branch("MR11", &MR11, "MR11/D");
      Tree->Branch("R1", &R1, "R1/D");
      Tree->Branch("R2", &R2, "R2/D");
      Tree->Branch("R3", &R3, "R3/D");
      Tree->Branch("R4", &R4, "R4/D");
      Tree->Branch("R5", &R5, "R5/D");
      Tree->Branch("R6", &R6, "R6/D");
      Tree->Branch("R7", &R7, "R7/D");
      Tree->Branch("R8", &R8, "R8/D");
      Tree->Branch("R9", &R9, "R9/D");
      Tree->Branch("R11a", &R11a, "R11a/D");
      Tree->Branch("R11b", &R11b, "R11b/D");
      Tree->Branch("R11c", &R11c, "R11c/D");
      Tree->Branch("R11d", &R11d, "R11d/D");
      Tree->Branch("R11e", &R11e, "R11e/D");
      Tree->Branch("R11f", &R11f, "R11f/D");
      Tree->Branch("R11g", &R11g, "R11g/D");
      Tree->Branch("WToENuCount", &WToENuCount, "WToENuCount/I");
      Tree->Branch("WToMuNuCount", &WToMuNuCount, "WToMuNuCount/I");
      Tree->Branch("WToTauNuCount", &WToTauNuCount, "WToTauNuCount/I");
   }
};
//---------------------------------------------------------------------------
#endif

