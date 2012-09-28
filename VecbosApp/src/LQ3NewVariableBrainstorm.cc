//---------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;
//---------------------------------------------------------------------------
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
//---------------------------------------------------------------------------
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "LQ3NewVariableBrainstorm.hh"
//---------------------------------------------------------------------------
void LQ3NewVariableBrainstorm::SetConditions(TTree* treeCond)
{
   _treeCond = treeCond;
}
//---------------------------------------------------------------------------
LQ3NewVariableBrainstorm::LQ3NewVariableBrainstorm(TTree *tree, bool isData, string JSONFile) : Vecbos(tree)
{
   _isData = isData;
   _goodRunLS = false;

   // To read good run list!
   if(_isData == true && JSONFile != "")
   {
      cout << "JSON File Used: " << JSONFile << endl;
      setJsonGoodRunList(JSONFile);
      fillRunLSMap();
      _goodRunLS = true;
   }
}
//---------------------------------------------------------------------------
LQ3NewVariableBrainstorm::~LQ3NewVariableBrainstorm()
{
}
//---------------------------------------------------------------------------
void LQ3NewVariableBrainstorm::Loop(string outFileName, int start, int stop)
{
   cerr << "[LQ3NewVariableBrainstorm] You've entered the LQ3 analysis main function." << endl;

   if(fChain == 0)
      return;

   if(stop <= start)
      return;

   if(fChain->GetEntries() == 0)
      return;

   TFile *file = new TFile(outFileName.c_str(), "RECREATE");

   // output record!
   LQ3NewVariableOutputRecord Messenger;

   // prepare output trees!!!
   TTree* outTree = new TTree("LQ3Tree", "");
   Messenger.MakeBranches(outTree);

   // start looping!!!
   unsigned int lastLumi = 0;
   unsigned int lastRun = 0;

   Long64_t nbytes = 0;
   Long64_t nentries = fChain->GetEntries();
   cout << "Number of entries = " << nentries << endl;

   for(Long64_t jentry = start; jentry < stop; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if(ientry < 0)
         break;

      if(jentry % 2000 == 0)
         cout << "Processing entry " << jentry << "/" << nentries << endl;

      Long64_t nb = fChain->GetEntry(jentry);
      nbytes += nb;

      Utils anaUtils;
     
      // Good Run selection - for data!!!
      if(_isData == true && _goodRunLS == true && isGoodRunLS() == false)
      {
         if(lastRun != runNumber || lastLumi != lumiBlock)
         {
            lastRun = runNumber;
            lastLumi = lumiBlock;
            std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected!!!" << std::endl;
         }
         continue;
      }

      vector<MuonCandidate> MuonCandidates = MakeMuonCandidates();
      vector<ElectronCandidate> ElectronCandidates = MakeElectronCandidates();

      // Check HLT bits for the event!!
      reloadTriggerMask(true);

      bool PassingHLTrigger = hasPassedHLT();   // not sure which bits to use yet!!!!!
      PassingHLTrigger = true;

      bool failHcal2011Filter = false;    // not in 41X MC!!!
      if(_isData == true)
         failHcal2011Filter = fail2011Filter;

      // checking jet id!!!!
      bool AllCaloJetsPassingID = true;
      for(int i = 0; i < nAK5Jet; i++)
         if(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i] > 30 * 30 && fabs(etaAK5Jet[i]) < 3)
            if(CheckCaloJetID(i) == false)
               AllCaloJetsPassingID = false;

      bool AllPFJetsPassingID = true;
      for(int i = 0; i < nAK5PFPUcorrJet; i++)
         if(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i] > 30 * 30
            && fabs(etaAK5PFPUcorrJet[i]) < 2.4)
            if(CheckPFJetID(i) == false)
               AllPFJetsPassingID = false;

      // requiring 2 jets above 60 GeV!!!!
      bool HasTwoCaloJets60 = false;
      if(PassingHLTrigger == true && failHcal2011Filter == false && AllCaloJetsPassingID == true)
      {
         int JetCount60 = 0;
         for(int i = 0; i < nAK5Jet; i++)
         {
            if(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i] > 60 * 60 && fabs(etaAK5Jet[i]) < 3)
               JetCount60 = JetCount60 + 1;
         }
         if(JetCount60 >= 2)
            HasTwoCaloJets60 = true;
      }

      bool HasTwoPFJets60 = false;
      if(PassingHLTrigger == true && failHcal2011Filter == false && AllPFJetsPassingID == true)
      {
         int JetCount60 = 0;
         for(int i = 0; i < nAK5PFPUcorrJet; i++)
            if(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i] > 60 * 60
               && fabs(etaAK5PFPUcorrJet[i]) < 2.4)
               JetCount60 = JetCount60 + 1;
         if(JetCount60 >= 2)
            HasTwoPFJets60 = true;
      }

      bool PassBaselineCalo = PassingHLTrigger && !failHcal2011Filter && AllCaloJetsPassingID && HasTwoCaloJets60;
      bool PassBaselinePF = PassingHLTrigger && !failHcal2011Filter && AllPFJetsPassingID && HasTwoPFJets60;

      // b-tagging counts!!!!
      int NumberOfCaloJetTCHEL = 0;
      int NumberOfCaloJetTCHEM = 0;    vector<int> CaloJetTCHEMIndex;
      int NumberOfCaloJetTCHET = 0;
      int NumberOfCaloJetSSVHEM = 0;
      for(int i = 0; i < nAK5Jet; i++)
      {
         double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
         if(PT < 40)                                 continue;
         if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)   continue;
         if(CheckCaloJetID(i) == false)              continue;

         if(caloJetPassTCHEL(i) == true)
            NumberOfCaloJetTCHEL = NumberOfCaloJetTCHEL + 1;
         if(caloJetPassTCHEM(i) == true)
         {
            NumberOfCaloJetTCHEM = NumberOfCaloJetTCHEM + 1;
            CaloJetTCHEMIndex.push_back(i);
         }
         if(caloJetPassTCHET(i) == true)
            NumberOfCaloJetTCHET = NumberOfCaloJetTCHET + 1;
         if(caloJetPassSSVHEM(i) == true)
            NumberOfCaloJetSSVHEM = NumberOfCaloJetSSVHEM + 1;
      }

      if(NumberOfCaloJetTCHEM != 2 || PassBaselineCalo == false)
         continue;

      vector<FourVector> GoodJets;
      for(int i = 0; i < nAK5Jet; i++)
      {
         double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
         if(PT < 40)                                 continue;
         if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)   continue;
         if(CheckCaloJetID(i) == false)              continue;

         FourVector JetP(energyAK5Jet[i], pxAK5Jet[i], pyAK5Jet[i], pzAK5Jet[i]);
         GoodJets.push_back(JetP);
      }
      vector<FourVector> Groups = SplitIntoGroups(GoodJets, true);

      FourVector Jet1(energyAK5Jet[CaloJetTCHEMIndex[0]], pxAK5Jet[CaloJetTCHEMIndex[0]],
         pyAK5Jet[CaloJetTCHEMIndex[0]], pzAK5Jet[CaloJetTCHEMIndex[0]]);
      FourVector Jet2(energyAK5Jet[CaloJetTCHEMIndex[1]], pxAK5Jet[CaloJetTCHEMIndex[1]],
         pyAK5Jet[CaloJetTCHEMIndex[1]], pzAK5Jet[CaloJetTCHEMIndex[1]]);
      FourVector ME(energyMet[0], pxMet[0], pyMet[0], pzMet[0]);

      FourVector ISR(0, 0, 0, 0);
      for(int i = 0; i < nAK5Jet; i++)
      {
         if(CaloJetTCHEMIndex[0] == i || CaloJetTCHEMIndex[1] == i)
            continue;
         
         double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
         if(PT < 40)                                 continue;
         if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)   continue;
         if(CheckCaloJetID(i) == false)              continue;

         FourVector JetP(energyAK5Jet[i], pxAK5Jet[i], pyAK5Jet[i], pzAK5Jet[i]);
         ISR = ISR + JetP;
      }

      int WToENuCount = 0;
      int WToMuNuCount = 0;
      int WToTauNuCount = 0;

      for(int i = 0; i < nMc; i++)
      {
         if(fabs(idMc[i]) == 11 && mothMc[i] >= 0 && fabs(idMc[mothMc[i]]) == 24)
            WToENuCount = WToENuCount + 1;
         if(fabs(idMc[i]) == 13 && mothMc[i] >= 0 && fabs(idMc[mothMc[i]]) == 24)
            WToMuNuCount = WToMuNuCount + 1;
         if(fabs(idMc[i]) == 15 && mothMc[i] >= 0 && fabs(idMc[mothMc[i]]) == 24)
            WToTauNuCount = WToTauNuCount + 1;
      }

      Messenger.WToENuCount = WToENuCount;
      Messenger.WToMuNuCount = WToMuNuCount;
      Messenger.WToTauNuCount = WToTauNuCount;

      Messenger.MR0 = Get2011MR(Groups[0], Groups[1]);
      Messenger.R0 = Get2011R(Groups[0], Groups[1], ME);
      Messenger.MR1 = GetISR2011MR(Jet1, Jet2, ME, 1);
      Messenger.R1 = GetISR2011R(Jet1, Jet2, ME, 1);
      Messenger.MR2 = GetISR2011MR(Jet1, Jet2, ME, 2);
      Messenger.R2 = GetISR2011R(Jet1, Jet2, ME, 2);
      Messenger.MR3 = GetISR2011MR(Jet1, Jet2, ME, 3);
      Messenger.R3 = GetISR2011R(Jet1, Jet2, ME, 3);
      Messenger.MR4 = GetISR2011MR(Jet1, Jet2, ME, 4);
      Messenger.R4 = GetISR2011R(Jet1, Jet2, ME, 4);
      Messenger.MR1 = GetISR2011MR(Jet1, Jet2, ME, 1);
      Messenger.R1 = GetISR2011R(Jet1, Jet2, ME, 1);
      Messenger.MR2 = GetISR2011MR(Jet1, Jet2, ME, 2);
      Messenger.R2 = GetISR2011R(Jet1, Jet2, ME, 2);
      Messenger.MR3 = GetISR2011MR(Jet1, Jet2, ME, 3);
      Messenger.R3 = GetISR2011R(Jet1, Jet2, ME, 3);
      Messenger.MR4 = GetISR2011MR(Jet1, Jet2, ME, 4);
      Messenger.R4 = GetISR2011R(Jet1, Jet2, ME, 4);
      Messenger.MR5 = GetISR2011MR(Jet1, Jet2, ME, 5);
      Messenger.R5 = GetISR2011R(Jet1, Jet2, ME, 5);
      Messenger.MR6 = GetISR2011MR(Jet1, Jet2, ME, 6);
      Messenger.R6 = GetISR2011R(Jet1, Jet2, ME, 6);
      Messenger.MR7 = GetISR2011MR(Jet1, Jet2, ME, 7);
      Messenger.R7 = GetISR2011R(Jet1, Jet2, ME, 7);
      Messenger.MR8 = GetISR2011MR(Jet1, Jet2, ME, 8);
      Messenger.R8 = GetISR2011R(Jet1, Jet2, ME, 8);
      Messenger.MR9 = GetISR2011MR(Jet1, Jet2, ME, 9);
      Messenger.R9 = GetISR2011R(Jet1, Jet2, ME, 9);
      Messenger.MR11 = GetISR2011MR(Jet1, Jet2, ME, 11, ISR);
      Messenger.R11a = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'a');
      Messenger.R11b = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'b');
      Messenger.R11c = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'c');
      Messenger.R11d = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'd');
      Messenger.R11e = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'e');
      Messenger.R11f = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'f');
      Messenger.R11g = GetISR2011R(Jet1, Jet2, ME, 11, ISR, 'g');

      outTree->Fill();
   }

   file->cd();
   outTree->Write();

   file->Close();
}
//---------------------------------------------------------------------------
int LQ3NewVariableBrainstorm::BestPV()
{
   // find the highestpT PV
   double maxpT = -9999.;
   for(int i = 0; i < nPV; i++)
   {
      if(SumPtPV[i] > maxpT)
      {
         iPV = i;
         maxpT = SumPtPV[i];
      }
   }

   return iPV;
}
//---------------------------------------------------------------------------
double LQ3NewVariableBrainstorm::pTMuon(int i)
{
   return sqrt(pxMuon[i] * pxMuon[i] + pyMuon[i] * pyMuon[i]);
}
//---------------------------------------------------------------------------
double LQ3NewVariableBrainstorm::SumPt(int iMu, int iZ)
{
   double eta0 = etaMuon[iMu];
   double phi0 = phiMuon[iMu];
   double sumPt_tmp = 0;
   for(int i = 0; i < nTrack; i++)
   {
      if(i == trackIndexMuon[iMu])
         continue; // take out the muon

      if(trackValidHitsTrack[i] < 5)
         continue;                                     // minimum number of hits  XXX

      if(fabs(transvImpactParTrack[i] / transvImpactParErrorTrack[i]) > 5.)
         continue;    // track incompatible with the vertex on (x,y)

      if(fabs(PVzPV[BestPV()] - trackVzTrack[i]) > 0.1)
         continue;              // track incompatible with the vertex on z

      TVector3 v(pxTrack[i], pyTrack[i], pzTrack[i]);
      if(sqrt(pow(v.Eta() - eta0, 2.) + pow(v.Phi() - phi0, 2.)) > 0.5)
         continue; // track outside the cone

      if(v.Pt() < 0.500)
         continue;     // minimum pT

      if(v.Pt() > 500.)
         continue;     // maximum pT

      sumPt_tmp += v.Pt();
   }

   return sumPt_tmp;
}
//---------------------------------------------------------------------------
double LQ3NewVariableBrainstorm::DeltaPhi_PiHalf(double phi1, double phi2)
{
   double dp = fabs(DeltaPhi(phi1, phi2));
   if(dp > asin(1.))
      dp = asin(1.) * 2. - dp;
   return dp;
}
//---------------------------------------------------------------------------
bool LQ3NewVariableBrainstorm::CheckCaloJetID(int Index)
{
   return true;   // not in the tree....

   // pure09, loose, 2011 Mar. 9
   if(etaAK5Jet[Index] < 3 && etaAK5Jet[Index] > -3)
   {
      if(nHit90AK5Jet[Index] <= 1)
         return false;
      if(emFracAK5Jet[Index] <= 0.01)
         return false;
      if(fHPDAK5Jet[Index] > 0.98)
         return false;
   }

   return true;
}
//---------------------------------------------------------------------------
bool LQ3NewVariableBrainstorm::CheckPFJetID(int Index)
{
   return true;   // not there....

   // loose jet ID, 2011 Mar. 9
   if(neutralHadronEnergyAK5PFPUcorrJet[Index] / energyAK5PFPUcorrJet[Index] >= 0.99)
      return false;
   if(neutralEmEnergyAK5PFPUcorrJet[Index] / energyAK5PFPUcorrJet[Index] >= 0.99)
      return false;
   if(chargedHadronMultiplicityAK5PFPUcorrJet[Index]
         + neutralHadronMultiplicityAK5PFPUcorrJet[Index]
         + photonMultiplicityAK5PFPUcorrJet[Index]
         + electronMultiplicityAK5PFPUcorrJet[Index]
         + muonMultiplicityAK5PFPUcorrJet[Index]
         + HFHadronMultiplicityAK5PFPUcorrJet[Index]
         + HFEMMultiplicityAK5PFPUcorrJet[Index] <= 1)
      return false;
   if(etaAK5PFPUcorrJet[Index] < 2.4 && etaAK5PFPUcorrJet[Index] > -2.4)
   {
      if(chargedHadronEnergyAK5PFPUcorrJet[Index] / energyAK5PFPUcorrJet[Index] <= 0)
         return false;
      if(chargedHadronMultiplicityAK5PFPUcorrJet[Index] + electronMultiplicityAK5PFPUcorrJet[Index]
            + muonMultiplicityAK5PFPUcorrJet[Index] <= 0)
         return false;
      if(chargedEmEnergyAK5PFPUcorrJet[Index] / energyAK5PFPUcorrJet[Index] >= 0.99)
         return false;
   }
   return true;
}
//---------------------------------------------------------------------------
bool LQ3NewVariableBrainstorm::CheckVertex(int Index)
{
   if(Index >= nPV || Index < 0)
      return false;
   if(Index >= 20)
   {
      cout << "Dude...there are more than 20 vertices.... (" << nPV << " to be exact)" << endl;
      return false;
   }

   if(ndofPV[Index] <= 4)
      return false;
   if(PVzPV[Index] > 15)
      return false;
   if(PVxPV[Index] * PVxPV[Index] + PVyPV[Index] * PVyPV[Index] > 2 * 2)
      return false;
   // if(isFakePV[Index] == true)
   //    return false;

   return true;
}
//---------------------------------------------------------------------------
bool LQ3NewVariableBrainstorm::CheckHLTBit(string PathName)
{
   // Check if a bit is set in the event.
   // Return false if nothing matches the input

   vector<int> TriggerBits;
   for(int i = 0; i < (int)nameHLT->size(); i++)
      if((*nameHLT)[i].find(PathName) != string::npos)
         TriggerBits.push_back(indexHLT[i]);

   if(TriggerBits.size() == 0)
      return false;

   Utils AnalysisUtilities;
   return AnalysisUtilities.getTriggersOR(TriggerBits, firedTrg);
}
//---------------------------------------------------------------------------
vector<MuonCandidate> LQ3NewVariableBrainstorm::MakeMuonCandidates()
{
   if(nMuon > 1000)
   {
      cout << "Something wrong!  Muon count = " << nMuon << endl;
      return vector<MuonCandidate>();
   }

   vector<MuonCandidate> MuonCandidates;
   MuonCandidates.resize(nMuon);

   Utils anaUtils;

   for(int i = 0; i < nMuon; i++)
   {
      int iTrack = trackIndexMuon[i];
      int iGlobalMuonTrack = combinedTrackIndexMuon[i];
      int iSTAMuonTrack = standAloneTrackIndexMuon[i];

      MuonCandidate &Candidate = MuonCandidates[i];

      Candidate.IsGlobal =
         (int)anaUtils.muonIdVal(muonIdMuon[i], bits::AllGlobalMuons);
      Candidate.IsPromptTight =
         (int)anaUtils.muonIdVal(muonIdMuon[i], bits::GlobalMuonPromptTight);
      Candidate.IsTracker =
         (int)anaUtils.muonIdVal(muonIdMuon[i], bits::AllTrackerMuons);

      Candidate.PixelHit = numberOfValidPixelBarrelHitsTrack[iTrack]
         + numberOfValidPixelEndcapHitsTrack[iTrack];
      Candidate.StripHit = numberOfValidStripTIBHitsTrack[iTrack]
         + numberOfValidStripTOBHitsTrack[iTrack]
         + numberOfValidStripTIDHitsTrack[iTrack]
         + numberOfValidStripTECHitsTrack[iTrack];
      Candidate.Chi2 = trackNormalizedChi2GlobalMuonTrack[iGlobalMuonTrack];
      Candidate.ValidMuonHit = trackValidHitsSTAMuonTrack[iSTAMuonTrack];
      Candidate.MuonStations = 10000;

      Candidate.Charge = chargeMuon[i];

      Candidate.Dxy = eleDxyPV(PVxPV[0], PVyPV[0], PVzPV[0],
         vertexXMuon[i], vertexYMuon[i], vertexZMuon[i],
         pxMuon[i], pyMuon[i], pzMuon[i]);
      Candidate.Isolation = sumPt03Muon[i];
      Candidate.CombinedIsolation = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i]
         - rhoFastjet * TMath::Pi() * 0.3 * 0.3;

      Candidate.Eta = etaMuon[i];
      Candidate.Phi = phiMuon[i];
      Candidate.PT = pTMuon(i);
      Candidate.P = energyMuon[i];
      Candidate.Px = pxMuon[i];
      Candidate.Py = pyMuon[i];
      Candidate.Pz = pzMuon[i];

      Candidate.PassMuonID = muonPassTight(i);
      Candidate.PassMuonTight = muonPassTight(i);
      Candidate.PassMuonLoose = muonPassLoose(i);
   }

   return MuonCandidates;
}
//---------------------------------------------------------------------------
vector<ElectronCandidate> LQ3NewVariableBrainstorm::MakeElectronCandidates()
{
   if(nEle > 1000)
   {
      cout << "Something wrong!  Electron count = " << nEle << endl;
      return vector<ElectronCandidate>();
   }

   vector<ElectronCandidate> ElectronCandidates;
   ElectronCandidates.resize(nEle);

   for(int i = 0; i < nEle; i++)
   {
      ElectronCandidate &Candidate = ElectronCandidates[i];

      int SuperClusterIndex = superClusterIndexEle[i];
      int TrackIndex = trackIndexEle[i];

      Candidate.SuperClusterEnergy = energySC[SuperClusterIndex];
      Candidate.SuperClusterEta = etaSC[SuperClusterIndex];
      Candidate.SuperClusterEtaWidth = etaWidthSC[SuperClusterIndex];
      Candidate.SuperClusterPhi = phiSC[SuperClusterIndex];
      Candidate.SuperClusterPhiWidth = phiWidthSC[SuperClusterIndex];
      Candidate.SuperClusterHcalTowerSumEt03 =
         hcalTowerSumEtConeDR03SC[SuperClusterIndex];
      Candidate.SuperClusterSigmaIEtaIEta = covIEtaIEtaSC[SuperClusterIndex];

      Candidate.MissingHits = trackLostHitsTrack[TrackIndex];

      Candidate.ConversionDistance = convDistEle[i];
      Candidate.ConversionDeltaCotTheta = convDcotEle[i];
      Candidate.DeltaEtaAtCalo = deltaEtaAtCaloEle[i];
      Candidate.DeltaPhiAtCalo = deltaPhiAtCaloEle[i];
      Candidate.HcalIsolation = dr03HcalTowerSumEtEle[i];
      Candidate.TrackIsolation = dr03TkSumPtEle[i];
      Candidate.EcalIsolation = dr03EcalRecHitSumEtEle[i];
      Candidate.CombinedIsolation = dr03HcalTowerSumEtEle[i] + dr03TkSumPtEle[i] + dr03EcalRecHitSumEtEle[i]
         - rhoFastjet * TMath::Pi() * 0.3 * 0.3;
      Candidate.HOverE = hOverEEle[i];
      Candidate.Phi = phiEle[i];
      Candidate.Eta = etaEle[i];
      Candidate.PT = sqrt(pxEle[i] * pxEle[i] + pyEle[i] * pyEle[i]);
      Candidate.P = energyEle[i];
      Candidate.Px = pxEle[i];
      Candidate.Py = pyEle[i];
      Candidate.Pz = pzEle[i];
      Candidate.Charge = chargeEle[i];

      double CombinedRelativeIsolation = (Candidate.HcalIsolation
         + Candidate.EcalIsolation + Candidate.TrackIsolation) / Candidate.PT;

      Candidate.PassWP80 = electronPassWP80(i);
      Candidate.PassWP85 = electronPassWP85(i);
      Candidate.PassWP90 = electronPassWP90(i);
      Candidate.PassWP95 = electronPassWP95(i);
   }

   return ElectronCandidates;
}
//---------------------------------------------------------------------------












