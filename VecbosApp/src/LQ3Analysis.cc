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
#include "LQ3Analysis.hh"
//---------------------------------------------------------------------------
void LQ3Analysis::SetConditions(TTree* treeCond)
{
   _treeCond = treeCond;
}
//---------------------------------------------------------------------------
LQ3Analysis::LQ3Analysis(TTree *tree, bool isData, string JSONFile) : Vecbos(tree)
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
LQ3Analysis::~LQ3Analysis()
{
}
//---------------------------------------------------------------------------
void LQ3Analysis::Loop(string outFileName, int start, int stop)
{
   cerr << "[LQ3Analysis] You've entered the LQ3 analysis main function." << endl;

   if(fChain == 0)
      return;

   if(stop <= start)
      return;

   if(fChain->GetEntries() == 0)
      return;

   TFile *file = new TFile(outFileName.c_str(), "RECREATE");

   // output record!
   LQ3OutputRecord Messenger;

   // prepare output trees!!!
   TTree* outTree = new TTree("LQ3Tree", "LQ3Tree");
   Messenger.MakeBranches(outTree);

   // prepare quality monitoring histograms!!!!
   map<string, TH1D *> QualityMonitoring1DHistograms = GenerateQM1DHistograms();
   map<string, TH2D *> QualityMonitoring2DHistograms = GenerateQM2DHistograms();

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

      QualityMonitoring2DHistograms["HProcessedEvents"]->Fill(0.0, 0.0);   // count events!
      QualityMonitoring1DHistograms["HProcessedEventsPU"]->Fill(nPU[0]);   // count events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(1);   // all events!
      QualityMonitoring1DHistograms["HCountingHistogramPFJet"]->Fill(1);   // all events!

      vector<MuonCandidate> MuonCandidates = MakeMuonCandidates();
      vector<ElectronCandidate> ElectronCandidates = MakeElectronCandidates();

      // Check HLT bits for the event!!
      reloadTriggerMask(true);

      bool PassingHLTrigger = hasPassedHLT();   // not sure which bits to use yet!!!!!

      bool failHcal2011Filter = false;    // not in 41X MC!!!

      if(PassingHLTrigger == true)
      {
         QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(2);   // passed HLT!
         QualityMonitoring1DHistograms["HCountingHistogramPFJet"]->Fill(2);   // passed HLT!
      }
      
      if(PassingHLTrigger == true && failHcal2011Filter == false)
      {
         QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(3);   // passed filter!
         QualityMonitoring1DHistograms["HCountingHistogramPFJet"]->Fill(3);   // passed filter!
      }

      // checking jet id!!!!
      bool AllCaloJetsPassingID = true;
      for(int i = 0; i < nAK5Jet; i++)
         if(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i] > 30 * 30 && fabs(etaAK5Jet[i]) < 3)
            if(CheckCaloJetID(i) == false)
               AllCaloJetsPassingID = false;
      if(PassingHLTrigger == true && failHcal2011Filter == false && AllCaloJetsPassingID == true)
      {
         QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(4);   // passed calo jet id?
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(4);   // passed calo jet id?
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(4);   // passed calo jet id?
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(4);   // passed calo jet id?
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(4);   // passed calo jet id?
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(4);   // passed calo jet id?
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(4);   // passed calo jet id?
      }

      bool AllPFJetsPassingID = true;
      for(int i = 0; i < nAK5PFPUcorrJet; i++)
         if(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i] > 30 * 30
            && fabs(etaAK5PFPUcorrJet[i]) < 2.4)
            if(CheckPFJetID(i) == false)
               AllPFJetsPassingID = false;
      if(PassingHLTrigger == true && failHcal2011Filter == false && AllPFJetsPassingID == true)
         QualityMonitoring1DHistograms["HCountingHistogramPFJet"]->Fill(4);   // passed PF jet id!!

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
      if(PassingHLTrigger == true && failHcal2011Filter == false && AllCaloJetsPassingID == true
         && HasTwoCaloJets60 == true)
      {
         QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(5);   // two jets above 60 GeV!!!!
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(5);
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(5);
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(5);
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(5);
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(5);
         QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(5);
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
      if(PassingHLTrigger == true && failHcal2011Filter == false && AllPFJetsPassingID == true
         && HasTwoPFJets60 == true)
         QualityMonitoring1DHistograms["HCountingHistogramPFJet"]->Fill(5);   // two jets above 60 GeV!!!!

      // fill output tree!!!   for all events!!!!!!
      Messenger.Clear();
      Messenger.RunNumber = runNumber;
      Messenger.EventNumber = eventNumber;
      Messenger.BunchCrossing = bunchCrossing;
      Messenger.LumiSection = lumiBlock;
      Messenger.Orbit = orbitNumber;
      Messenger.PileUp = nPU[0];
      FillTopMuons(MuonCandidates, Messenger);
      FillTopElectrons(ElectronCandidates, Messenger);
      FillCaloJet(Messenger);
      FillPFJet(Messenger);
      FillPV(Messenger);
      FillMET(Messenger);
      Messenger.PassHLT = PassingHLTrigger;
      Messenger.PassNoiseFilter = !failHcal2011Filter;
      Messenger.PassCaloJetID = AllCaloJetsPassingID;
      Messenger.PassCaloJet60 = HasTwoCaloJets60;
      Messenger.PassPFJetID = AllPFJetsPassingID;
      Messenger.PassPFJet60 = HasTwoPFJets60;
      Messenger.PassDiJetAve30 = CheckHLTBit("HLT_DiJetAve30_");
      Messenger.PassDiJetAve60 = CheckHLTBit("HLT_DiJetAve60_");
      Messenger.PassDiJetAve80 = CheckHLTBit("HLT_DiJetAve80_");
      Messenger.PassDiJetAve110 = CheckHLTBit("HLT_DiJetAve110_");
      Messenger.PassDiJetAve150 = CheckHLTBit("HLT_DiJetAve150_");
      Messenger.PassDiJetAve190 = CheckHLTBit("HLT_DiJetAve190_");
      Messenger.PassDiJetAve240 = CheckHLTBit("HLT_DiJetAve240_");
      Messenger.PassDiJetAve300 = CheckHLTBit("HLT_DiJetAve300_");
      Messenger.PassDiJetAve370 = CheckHLTBit("HLT_DiJetAve370_");
      Messenger.PassDiJetAve15U = CheckHLTBit("HLT_DiJetAve15U_");
      Messenger.PassDiJetAve30U = CheckHLTBit("HLT_DiJetAve30U_");
      Messenger.PassDiJetAve50U = CheckHLTBit("HLT_DiJetAve50U_");
      Messenger.PassDiJetAve70U = CheckHLTBit("HLT_DiJetAve70U_");
      Messenger.PassDiJetAve100U = CheckHLTBit("HLT_DiJetAve100U_");
      Messenger.PassDiJetAve140U = CheckHLTBit("HLT_DiJetAve140U_");
      Messenger.PassDiJetAve180U = CheckHLTBit("HLT_DiJetAve180U_");
      Messenger.PassDiJetAve300U = CheckHLTBit("HLT_DiJetAve300U_");
      Messenger.PassSingleMu3 = CheckHLTBit("HLT_Mu3_");
      Messenger.PassSingleMu5 = CheckHLTBit("HLT_Mu5_");
      Messenger.PassSingleMu8 = CheckHLTBit("HLT_Mu8_");
      Messenger.PassSingleMu12 = CheckHLTBit("HLT_Mu12_");
      Messenger.PassSingleMu15 = CheckHLTBit("HLT_Mu15_");
      Messenger.PassSingleMu20 = CheckHLTBit("HLT_Mu20_");
      Messenger.PassSingleMu24 = CheckHLTBit("HLT_Mu24_");
      Messenger.PassSingleMu30 = CheckHLTBit("HLT_Mu30_");
      Messenger.PassSingleMu40 = CheckHLTBit("HLT_Mu40_");
      Messenger.PassSingleMu100 = CheckHLTBit("HLT_Mu100_");
      Messenger.PassR020MR500 = CheckHLTBit("HLT_R020_MR500_");
      Messenger.PassR020MR550 = CheckHLTBit("HLT_R020_MR550_");
      Messenger.PassR023MR550 = CheckHLTBit("HLT_R023_MR550_");
      Messenger.PassR025MR400 = CheckHLTBit("HLT_R025_MR400_");
      Messenger.PassR025MR450 = CheckHLTBit("HLT_R025_MR450_");
      Messenger.PassR029MR450 = CheckHLTBit("HLT_R029_MR450_");
      Messenger.PassR033MR300 = CheckHLTBit("HLT_R033_MR300_");
      Messenger.PassR033MR350 = CheckHLTBit("HLT_R033_MR300_");
      Messenger.PassR036MR350 = CheckHLTBit("HLT_R036_MR350_");
      Messenger.PassR038MR200 = CheckHLTBit("HLT_R038_MR200_");
      Messenger.PassR038MR250 = CheckHLTBit("HLT_R038_MR250_");
      Messenger.PassR042MR250 = CheckHLTBit("HLT_R042_MR250_");
      Messenger.PassEleL10R020MR200 = CheckHLTBit("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200_");
      Messenger.PassEleL10R025MR200 = CheckHLTBit("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200_")
         || CheckHLTBit("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R025_MR200_");
      Messenger.PassEleL10R029MR200 = CheckHLTBit("HLT_Ele10_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_R029_MR200_");
      Messenger.PassEleT10R020MR200 = CheckHLTBit("HLT_Ele10_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200_");
      Messenger.PassMu8R020MR200 = CheckHLTBit("HLT_Mu8_R020_MR200_");
      Messenger.PassMu8R025MR200 = CheckHLTBit("HLT_Mu8_R025_MR200_");
      Messenger.PassMu8R029MR200 = CheckHLTBit("HLT_Mu8_R029_MR200_");
      Messenger.PassMu17Mu8 = CheckHLTBit("HLT_Mu17_Mu8_");
      Messenger.PassMu13Mu8 = CheckHLTBit("HLT_Mu13_Mu8_");
      Messenger.PassMu8EleL17 = CheckHLTBit("HLT_Mu8_Ele17_CaloIdL_");
      Messenger.PassMu17EleL8 = CheckHLTBit("HLT_Mu17_Ele8_CaloIdL_");
      outTree->Fill();

      bool PassBaselineCalo = PassingHLTrigger && !failHcal2011Filter && AllCaloJetsPassingID && HasTwoCaloJets60;
      bool PassBaselinePF = PassingHLTrigger && !failHcal2011Filter && AllPFJetsPassingID && HasTwoPFJets60;

      // Data quality-monitoring plots --- basic quantities!
      QMFillBasicInformation(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);

      // Basic jet histograms!!!
      if(PassBaselineCalo == true)
         QMFillCaloJet(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);
      if(PassBaselineCalo)
         QMFillPFJet(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);

      // Basic MET histograms!!!!!
      QMFillCaloMET(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);
      QMFillPFMET(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);

      // B-tagging histograms!!!!!!!!!!!!!!!
      QMFillCaloBTag(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);
      QMFillPFBTag(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms);

      // b-tagging counts!!!!
      int NumberOfCaloJetTCHEL = 0;
      int NumberOfCaloJetTCHET = 0;
      int NumberOfCaloJetSSVHEM = 0;
      for(int i = 0; i < nAK5Jet; i++)
      {
         double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
         if(PT < 30)                                 continue;
         if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)   continue;
         if(CheckCaloJetID(i) == false)              continue;

         if(caloJetPassTCHEL(i) == true)
            NumberOfCaloJetTCHEL = NumberOfCaloJetTCHEL + 1;
         if(caloJetPassTCHET(i) == true)
            NumberOfCaloJetTCHET = NumberOfCaloJetTCHET + 1;
         if(caloJetPassSSVHEM(i) == true)
            NumberOfCaloJetSSVHEM = NumberOfCaloJetSSVHEM + 1;
      }

      if(PassBaselineCalo == true)
      {
         QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(7);
         if(NumberOfCaloJetTCHEL >= 1)
            QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(7);
         if(NumberOfCaloJetTCHEL >= 2)
            QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(7);
         if(NumberOfCaloJetTCHET >= 1)
            QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(7);
         if(NumberOfCaloJetTCHET >= 2)
            QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(7);
         if(NumberOfCaloJetSSVHEM >= 1)
            QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(7);
         if(NumberOfCaloJetSSVHEM >= 2)
            QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(7);
      }

      // Group into hemispheres!!!
      vector<FourVector> CaloHemispheres(2);
      if(PassBaselineCalo == true)
      {
         vector<FourVector> CaloJetInput;
         for(int i = 0; i < nAK5Jet; i++)
         {
            double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
            if(PT < 30)                                 continue;
            if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)   continue;
            // if(CheckCaloJetID(i) == false)              continue;

            CaloJetInput.push_back(FourVector(energyAK5Jet[i], pxAK5Jet[i], pyAK5Jet[i], pzAK5Jet[i]));
         }
         if(CaloJetInput.size() <= 20)   // > <
            CaloHemispheres = SplitIntoGroups(CaloJetInput, false);
         else
         {
            CaloHemispheres.push_back(FourVector(0, 0, 0, 0));
            CaloHemispheres.push_back(FourVector(0, 0, 0, 0));
         }

         QMFillCaloHemisphere(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms, CaloHemispheres);

         CaloHemispheres[0][0] = CaloHemispheres[0].GetP();
         CaloHemispheres[1][0] = CaloHemispheres[1].GetP();
      }

      // Calculate MR/R (starred version)
      double CaloMRStar = GetMRStar(CaloHemispheres[0], CaloHemispheres[1]);
      double CaloRStar = GetRStar(CaloHemispheres[0], CaloHemispheres[1], FourVector(0, pxMet[0], pyMet[0], 0));

      if(PassBaselineCalo == true)
         QMFillCaloRazor(QualityMonitoring1DHistograms, QualityMonitoring2DHistograms, CaloHemispheres,
            NumberOfCaloJetTCHEL, NumberOfCaloJetSSVHEM, NumberOfCaloJetTCHET);

      // Simple MR/R cut at 200, 0.20 to get a feel (we won't go lower than this)
      if(PassBaselineCalo == true)
      {
         if(CaloMRStar > 200 && CaloRStar > 0.2)
         {
            QualityMonitoring1DHistograms["HCountingHistogramCaloJet"]->Fill(8);
            if(NumberOfCaloJetTCHEL >= 1)
               QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL"]->Fill(8);
            if(NumberOfCaloJetTCHEL >= 2)
               QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHEL2"]->Fill(8);
            if(NumberOfCaloJetTCHET >= 1)
               QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET"]->Fill(8);
            if(NumberOfCaloJetTCHET >= 2)
               QualityMonitoring1DHistograms["HCountingHistogramCaloJetTCHET2"]->Fill(8);
            if(NumberOfCaloJetSSVHEM >= 1)
               QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM"]->Fill(8);
            if(NumberOfCaloJetSSVHEM >= 2)
               QualityMonitoring1DHistograms["HCountingHistogramCaloJetSSVHEM2"]->Fill(8);
         }
      }
   }

   file->cd();
   outTree->Write();

   // Write QM histograms!!!
   file->cd();
   WriteQMHistograms(QualityMonitoring1DHistograms);
   WriteQMHistograms(QualityMonitoring2DHistograms);

   // Delete Histograms!!!!!!
   DeleteQMHistograms(QualityMonitoring1DHistograms);
   DeleteQMHistograms(QualityMonitoring2DHistograms);

   file->Close();
}
//---------------------------------------------------------------------------
int LQ3Analysis::BestPV()
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
double LQ3Analysis::pTMuon(int i)
{
   return sqrt(pxMuon[i] * pxMuon[i] + pyMuon[i] * pyMuon[i]);
}
//---------------------------------------------------------------------------
double LQ3Analysis::SumPt(int iMu, int iZ)
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
double LQ3Analysis::DeltaPhi_PiHalf(double phi1, double phi2)
{
   double dp = fabs(DeltaPhi(phi1, phi2));
   if(dp > asin(1.))
      dp = asin(1.) * 2. - dp;
   return dp;
}
//---------------------------------------------------------------------------
map<string, TH1D *> LQ3Analysis::GenerateQM1DHistograms()
{
   map<string, TH1D *> Histograms;

   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJet",
      new TH1D("HCountingHistogramCaloJet",
      "Counting number of events passing each step (calo)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(7, "(B-tagging)");
   Histograms["HCountingHistogramCaloJet"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");
   
   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJetTCHEL",
      new TH1D("HCountingHistogramCaloJetTCHEL",
      "Counting number of events passing each step (calo, TCHEL x1)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramCaloJetTCHEL"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");
   
   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJetTCHEL2",
      new TH1D("HCountingHistogramCaloJetTCHEL2",
      "Counting number of events passing each step (calo, TCHEL x2)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramCaloJetTCHEL2"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");

   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJetTCHET",
      new TH1D("HCountingHistogramCaloJetTCHET",
      "Counting number of events passing each step (calo, TCHET x1)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramCaloJetTCHET"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");
   
   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJetTCHET2",
      new TH1D("HCountingHistogramCaloJetTCHET2",
      "Counting number of events passing each step (calo, TCHET x2)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramCaloJetTCHET2"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");

   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJetSSVHEM",
      new TH1D("HCountingHistogramCaloJetSSVHEM",
      "Counting number of events passing each step (calo, SSVHEM x1)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramCaloJetSSVHEM"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");
   
   Histograms.insert(pair<string, TH1D *>("HCountingHistogramCaloJetSSVHEM2",
      new TH1D("HCountingHistogramCaloJetSSVHEM2",
      "Counting number of events passing each step (calo, SSVHEM x2)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramCaloJetSSVHEM2"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");

   Histograms.insert(pair<string, TH1D *>("HCountingHistogramPFJet",
      new TH1D("HCountingHistogramPFJet",
      "Counting number of events passing each step (PF)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(7, "(B-tagging)");
   Histograms["HCountingHistogramPFJet"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");

   /*
   Histograms.insert(pair<string, TH1D *>("HCountingHistogramPFJetTCHEL",
      new TH1D("HCountingHistogramPFJetTCHEL",
      "Counting number of events passing each step (PF, TCHEL)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramPFJetTCHEL"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");

   Histograms.insert(pair<string, TH1D *>("HCountingHistogramPFJetTCHEL2",
      new TH1D("HCountingHistogramPFJetTCHEL2",
      "Counting number of events passing each step (PF, TCHEL x2)", 20, 0.5, 20.5)));
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(1, "All events");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(2, "Pass trigger");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(3, "Hcal baseline filter");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(4, "All jets passing ID");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(5, "Two jets above 60 GeV");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(6, "(MET)");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(7, "B-tagging");
   Histograms["HCountingHistogramPFJetTCHEL2"]->GetXaxis()->SetBinLabel(8, "MRStar 200, RStar 0.2");
   */

   Histograms.insert(pair<string, TH1D *>("HBunchCrossing",
      new TH1D("HBunchCrossing", "Bunch crossing of events", 3502, -1.5, 3500.5)));

   Histograms.insert(pair<string, TH1D *>("HProcessedEventsPU",
      new TH1D("HProcessedEventsPU", "PU distribution of processed events", 51, -0.5, 50.5)));

   Histograms.insert(pair<string, TH1D *>("HGoodPV",
      new TH1D("HGoodPV", "Number of good primary vertices", 20, -0.5, 19.5)));

   Histograms.insert(pair<string, TH1D *>("HNumberOfBadCaloJet30",
      new TH1D("HNumberOfBadCaloJet30", "Number of bad calo jets (30 GeV)", 20, -0.5, 19.5)));
   Histograms.insert(pair<string, TH1D *>("HNumberOfGoodCaloJet30",
      new TH1D("HNumberOfGoodCaloJet30", "Number of good calo jets (30 GeV)", 20, -0.5, 19.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloJetPT",
      new TH1D("HCaloJetPT", "PT of all (good) calo jets!", 100, 0, 1000)));
   Histograms.insert(pair<string, TH1D *>("HCaloJetEta",
      new TH1D("HCaloJetEta", "Eta of all (good) calo jets!", 100, -5, 5)));
   Histograms.insert(pair<string, TH1D *>("HCaloJetPhi",
      new TH1D("HCaloJetPhi", "Phi of all (good) calo jets!", 100, -3.1415926535, 3.1415926535)));

   Histograms.insert(pair<string, TH1D *>("HNumberOfBadPFJet30",
      new TH1D("HNumberOfBadPFJet30", "Number of bad PF jets (30 GeV)", 20, -0.5, 19.5)));
   Histograms.insert(pair<string, TH1D *>("HNumberOfGoodPFJet30",
      new TH1D("HNumberOfGoodPFJet30", "Number of good PF jets (30 GeV)", 20, -0.5, 19.5)));
   Histograms.insert(pair<string, TH1D *>("HPFJetPT",
      new TH1D("HPFJetPT", "PT of all (good) PF jets!", 100, 0, 1000)));
   Histograms.insert(pair<string, TH1D *>("HPFJetEta",
      new TH1D("HPFJetEta", "Eta of all (good) PF jets!", 100, -5, 5)));
   Histograms.insert(pair<string, TH1D *>("HPFJetPhi",
      new TH1D("HPFJetPhi", "Phi of all (good) PF jets!", 100, -3.1415926535, 3.1415926535)));

   Histograms.insert(pair<string, TH1D *>("HCaloTCHE",
      new TH1D("HCaloTCHE", "Calo jets, TCHE", 100, 0, 15)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHP",
      new TH1D("HCaloTCHP", "Calo jets, TCHP", 100, 0, 15)));
   Histograms.insert(pair<string, TH1D *>("HCaloSSVHE",
      new TH1D("HCaloSSVHE", "Calo jets, SSVHE", 100, 0, 15)));
   Histograms.insert(pair<string, TH1D *>("HCaloSSVHP",
      new TH1D("HCaloSSVHP", "Calo jets, SSVHP", 100, 0, 15)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHELCount",
      new TH1D("HCaloTCHELCount", "Number of calo jets tagged by TCHE(L)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHEMCount",
      new TH1D("HCaloTCHEMCount", "Number of calo jets tagged by TCHE(M)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHETCount",
      new TH1D("HCaloTCHETCount", "Number of calo jets tagged by TCHE(T)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHPLCount",
      new TH1D("HCaloTCHPLCount", "Number of calo jets tagged by TCHP(L)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHPMCount",
      new TH1D("HCaloTCHPMCount", "Number of calo jets tagged by TCHP(M)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloTCHPTCount",
      new TH1D("HCaloTCHPTCount", "Number of calo jets tagged by TCHP(T)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloSSVHEMCount",
      new TH1D("HCaloSSVHEMCount", "Number of calo jets tagged by SSVHE(M)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloSSVHETCount",
      new TH1D("HCaloSSVHETCount", "Number of calo jets tagged by SSVHE(T)", 10, -0.5, 9.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloSSVHPTCount",
      new TH1D("HCaloSSVHPTCount", "Number of calo jets tagged by SSVHP(T)", 10, -0.5, 9.5)));

   Histograms.insert(pair<string, TH1D *>("HCaloHemispherePT",
      new TH1D("HCaloHemispherePT", "Calo-hemisphere PT", 100, 0, 600)));
   Histograms.insert(pair<string, TH1D *>("HCaloHemisphereMass",
      new TH1D("HCaloHemisphereMass", "Calo-hemisphere mass", 100, 0, 1000)));
   Histograms.insert(pair<string, TH1D *>("HCaloHemisphereEta",
      new TH1D("HCaloHemisphereEta", "Calo-hemisphere eta", 100, -5, 5)));
   Histograms.insert(pair<string, TH1D *>("HCaloHemispherePhi",
      new TH1D("HCaloHemispherePhi", "Calo-hemisphere phi", 100, -3.1415926535, 3.1415926535)));

   Histograms.insert(pair<string, TH1D *>("HCaloMRStar",
      new TH1D("HCaloMRStar", "MRStar from calo system", 100, 0, 1500)));
   Histograms.insert(pair<string, TH1D *>("HCaloMRStarTCHEL",
      new TH1D("HCaloMRStarTCHEL", "MRStar from calo system (TCHEL x1)", 100, 0, 1500)));
   Histograms.insert(pair<string, TH1D *>("HCaloMRStarTCHEL2",
      new TH1D("HCaloMRStarTCHEL2", "MRStar from calo system (TCHEL x2)", 100, 0, 1500)));
   Histograms.insert(pair<string, TH1D *>("HCaloMRStarTCHET",
      new TH1D("HCaloMRStarTCHET", "MRStar from calo system (TCHET x1)", 100, 0, 1500)));
   Histograms.insert(pair<string, TH1D *>("HCaloMRStarTCHET2",
      new TH1D("HCaloMRStarTCHET2", "MRStar from calo system (TCHET x2)", 100, 0, 1500)));
   Histograms.insert(pair<string, TH1D *>("HCaloMRStarSSVHEM",
      new TH1D("HCaloMRStarSSVHEM", "MRStar from calo system (SSVHEM x1)", 100, 0, 1500)));
   Histograms.insert(pair<string, TH1D *>("HCaloMRStarSSVHEM2",
      new TH1D("HCaloMRStarSSVHEM2", "MRStar from calo system (SSVHEM x2)", 100, 0, 1500)));
   
   Histograms.insert(pair<string, TH1D *>("HCaloRStar",
      new TH1D("HCaloRStar", "RStar from calo system", 100, 0, 1.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloRStarTCHEL",
      new TH1D("HCaloRStarTCHEL", "RStar from calo system (TCHEL x1)", 100, 0, 1.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloRStarTCHEL2",
      new TH1D("HCaloRStarTCHEL2", "RStar from calo system (TCHEL x2)", 100, 0, 1.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloRStarTCHET",
      new TH1D("HCaloRStarTCHET", "RStar from calo system (TCHET x1)", 100, 0, 1.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloRStarTCHET2",
      new TH1D("HCaloRStarTCHET2", "RStar from calo system (TCHET x2)", 100, 0, 1.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloRStarSSVHEM",
      new TH1D("HCaloRStarSSVHEM", "RStar from calo system (SSVHEM x1)", 100, 0, 1.5)));
   Histograms.insert(pair<string, TH1D *>("HCaloRStarSSVHEM2",
      new TH1D("HCaloRStarSSVHEM2", "RStar from calo system (SSVHEM x2)", 100, 0, 1.5)));
   
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStar",
      new TH1D("HCaloGammaRStar", "GammaRStar from calo system", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarTCHEL",
      new TH1D("HCaloGammaRStarTCHEL", "GammaRStar from calo system (TCHEL x1)", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarTCHEL2",
      new TH1D("HCaloGammaRStarTCHEL2", "GammaRStar from calo system (TCHEL x2)", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarTCHET",
      new TH1D("HCaloGammaRStarTCHET", "GammaRStar from calo system (TCHET x1)", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarTCHET2",
      new TH1D("HCaloGammaRStarTCHET2", "GammaRStar from calo system (TCHET x2)", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarSSVHEM",
      new TH1D("HCaloGammaRStarSSVHEM", "GammaRStar from calo system (SSVHEM x1)", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarSSVHEM2",
      new TH1D("HCaloGammaRStarSSVHEM2", "GammaRStar from calo system (SSVHEM x2)", 100, 0, 10)));
   
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarWithMETCut",
      new TH1D("HCaloGammaRStarWithMETCut", "GammaRStar from calo system (with MET cut (60))", 100, 0, 10)));
   Histograms.insert(pair<string, TH1D *>("HCaloGammaRStarWithHigherMETCut",
      new TH1D("HCaloGammaRStarWithHigherMETCut", "GammaRStar from calo system (with higher MET cut (120))",
      100, 0, 10)));

   return Histograms;
}
//---------------------------------------------------------------------------
map<string, TH2D *> LQ3Analysis::GenerateQM2DHistograms()
{
   map<string, TH2D *> Histograms;

   Histograms.insert(pair<string, TH2D *>("HProcessedEvents",
      new TH2D("HProcessedEvents", "Number of events processed", 1, -0.5, 0.5, 1, -0.5, 0.5)));

   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStar",
      new TH2D("HCaloMRStarVsRStar", "MRStar vs. RStar from calo system", 100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarTCHEL",
      new TH2D("HCaloMRStarVsRStarTCHEL", "MRStar vs. RStar from calo system (TCHEL x1)",
      100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarTCHEL2",
      new TH2D("HCaloMRStarVsRStarTCHEL2", "MRStar vs. RStar from calo system (TCHEL x2)",
      100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarTCHET",
      new TH2D("HCaloMRStarVsRStarTCHET", "MRStar vs. RStar from calo system (TCHET x1)",
      100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarTCHET2",
      new TH2D("HCaloMRStarVsRStarTCHET2", "MRStar vs. RStar from calo system (TCHET x2)",
      100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarSSVHEM",
      new TH2D("HCaloMRStarVsRStarSSVHEM", "MRStar vs. RStar from calo system (SSVHEM x1)",
      100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarSSVHEM2",
      new TH2D("HCaloMRStarVsRStarSSVHEM2", "MRStar vs. RStar from calo system (SSVHEM x2)",
      100, 0, 1500, 100, 0, 1.5)));
   
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarWithMETCut",
      new TH2D("HCaloMRStarVsRStarWithMETCut", "MRStar vs. RStar from calo system (with MET cut (60))",
      100, 0, 1500, 100, 0, 1.5)));
   Histograms.insert(pair<string, TH2D *>("HCaloMRStarVsRStarWithHigherMETCut",
      new TH2D("HCaloMRStarVsRStarWithHigherMETCut",
      "MRStar vs. RStar from calo system (with higher MET cut (120))", 100, 0, 1500, 100, 0, 1.5)));

   return Histograms;
}
//---------------------------------------------------------------------------
void LQ3Analysis::WriteQMHistograms(map<string, TH1D *> Histograms)
{
   for(map<string, TH1D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
   {
      if(iter->second != NULL)
         iter->second->Write();
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::WriteQMHistograms(map<string, TH2D *> Histograms)
{
   for(map<string, TH2D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
   {
      if(iter->second != NULL)
         iter->second->Write();
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::DeleteQMHistograms(map<string, TH1D *> Histograms)
{
   for(map<string, TH1D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
   {
      if(iter->second == NULL)
         continue;

      delete iter->second;
      iter->second = NULL;
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::DeleteQMHistograms(map<string, TH2D *> Histograms)
{
   for(map<string, TH2D *>::iterator iter = Histograms.begin(); iter != Histograms.end(); iter++)
   {
      if(iter->second == NULL)
         continue;

      delete iter->second;
      iter->second = NULL;
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillBasicInformation(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   QM1DHistograms["HBunchCrossing"]->Fill(bunchCrossing);

   if(nPV < 1000)
   {
      int GoodPVCount = 0;
      for(int i = 0; i < nPV && i < 20; i++)
      {
         if(CheckVertex(i) == false)
            continue;
         GoodPVCount = GoodPVCount + 1;
      }
      QM1DHistograms["HGoodPV"]->Fill(GoodPVCount);
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillCaloJet(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   int NumberOfBadJet30 = 0;
   int NumberOfGoodJet30 = 0;

   for(int i = 0; i < nAK5Jet; i++)
   {
      double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
      if(PT < 30)
         continue;
      if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)
         continue;

      if(CheckCaloJetID(i) == true)
         NumberOfGoodJet30 = NumberOfGoodJet30 + 1;
      else
         NumberOfBadJet30 = NumberOfBadJet30 + 1;
   }

   QM1DHistograms["HNumberOfGoodCaloJet30"]->Fill(NumberOfGoodJet30);
   QM1DHistograms["HNumberOfBadCaloJet30"]->Fill(NumberOfBadJet30);

   for(int i = 0; i < nAK5Jet; i++)
   {
      double PT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
      if(PT < 30)
         continue;
      if(etaAK5Jet[i] > 3 || etaAK5Jet[i] < -3)
         continue;
      if(CheckCaloJetID(i) == false)
         continue;

      QM1DHistograms["HCaloJetPT"]->Fill(PT);
      QM1DHistograms["HCaloJetEta"]->Fill(etaAK5Jet[i]);
      QM1DHistograms["HCaloJetPhi"]->Fill(phiAK5Jet[i]);
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillPFJet(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   int NumberOfBadJet30 = 0;
   int NumberOfGoodJet30 = 0;

   for(int i = 0; i < nAK5PFPUcorrJet; i++)
   {
      double PT = sqrt(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i]);
      if(PT < 30)
         continue;
      if(etaAK5PFPUcorrJet[i] > 2.4 || etaAK5PFPUcorrJet[i] < -2.4)
         continue;

      if(CheckPFJetID(i) == true)
         NumberOfGoodJet30 = NumberOfGoodJet30 + 1;
      else
         NumberOfBadJet30 = NumberOfBadJet30 + 1;
   }

   QM1DHistograms["HNumberOfGoodPFJet30"]->Fill(NumberOfGoodJet30);
   QM1DHistograms["HNumberOfBadPFJet30"]->Fill(NumberOfBadJet30);

   for(int i = 0; i < nAK5PFPUcorrJet; i++)
   {
      double PT = sqrt(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i]);
      if(PT < 30)
         continue;
      if(etaAK5PFPUcorrJet[i] > 2.4 || etaAK5PFPUcorrJet[i] < -2.4)
         continue;
      if(CheckPFJetID(i) == false)
         continue;

      QM1DHistograms["HPFJetPT"]->Fill(PT);
      QM1DHistograms["HPFJetEta"]->Fill(etaAK5PFJet[i]);
      QM1DHistograms["HPFJetPhi"]->Fill(phiAK5PFJet[i]);
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillCaloMET(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   // TODO
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillPFMET(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   // TODO
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillCaloBTag(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   int TCHELCount = 0;
   int TCHEMCount = 0;
   int TCHETCount = 0;
   int TCHPLCount = 0;
   int TCHPMCount = 0;
   int TCHPTCount = 0;
   int SSVHEMCount = 0;
   int SSVHETCount = 0;
   int SSVHPTCount = 0;

   for(int i = 0; i < nAK5Jet; i++)
   {
      double PT = sqrt(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i]);
      if(PT < 30)
         continue;
      if(etaAK5PFPUcorrJet[i] > 2.4 || etaAK5PFPUcorrJet[i] < -2.4)
         continue;
      if(CheckCaloJetID(i) == false)
         continue;
      
      QM1DHistograms["HCaloTCHE"]->Fill(trackCountingHighEffBJetTagsAK5Jet[i]);
      QM1DHistograms["HCaloTCHP"]->Fill(trackCountingHighPurBJetTagsAK5Jet[i]);
      QM1DHistograms["HCaloSSVHE"]->Fill(simpleSecondaryVertexHighEffBJetTagsAK5Jet[i]);
      QM1DHistograms["HCaloSSVHP"]->Fill(simpleSecondaryVertexHighPurBJetTagsAK5Jet[i]);

      if(trackCountingHighEffBJetTagsAK5Jet[i] > 1.7)
         TCHELCount = TCHELCount + 1;
      if(trackCountingHighEffBJetTagsAK5Jet[i] > 3.3)
         TCHEMCount = TCHEMCount + 1;
      if(trackCountingHighEffBJetTagsAK5Jet[i] > 10.2)
         TCHETCount = TCHETCount + 1;
      if(trackCountingHighPurBJetTagsAK5Jet[i] > 1.19)
         TCHPLCount = TCHPLCount + 1;
      if(trackCountingHighPurBJetTagsAK5Jet[i] > 1.93)
         TCHPMCount = TCHPMCount + 1;
      if(trackCountingHighPurBJetTagsAK5Jet[i] > 3.41)
         TCHPTCount = TCHPTCount + 1;
      if(simpleSecondaryVertexHighEffBJetTagsAK5Jet[i] > 1.74)
         SSVHEMCount = SSVHEMCount + 1;
      if(simpleSecondaryVertexHighEffBJetTagsAK5Jet[i] > 3.05)
         SSVHETCount = SSVHETCount + 1;
      if(simpleSecondaryVertexHighPurBJetTagsAK5Jet[i] > 2.00)
         SSVHPTCount = SSVHPTCount + 1;
   }

   QM1DHistograms["HCaloTCHELCount"]->Fill(TCHELCount);
   QM1DHistograms["HCaloTCHEMCount"]->Fill(TCHEMCount);
   QM1DHistograms["HCaloTCHETCount"]->Fill(TCHETCount);
   QM1DHistograms["HCaloTCHPLCount"]->Fill(TCHPLCount);
   QM1DHistograms["HCaloTCHPMCount"]->Fill(TCHPMCount);
   QM1DHistograms["HCaloTCHPTCount"]->Fill(TCHPTCount);
   QM1DHistograms["HCaloSSVHEMCount"]->Fill(SSVHEMCount);
   QM1DHistograms["HCaloSSVHETCount"]->Fill(SSVHETCount);
   QM1DHistograms["HCaloSSVHPTCount"]->Fill(SSVHPTCount);
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillPFBTag(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms)
{
   // TODO
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillCaloHemisphere(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
   vector<FourVector> &Groups)
{
   for(int i = 0; i < (int)Groups.size(); i++)
   {
      QM1DHistograms["HCaloHemispherePT"]->Fill(Groups[i].GetPT());
      QM1DHistograms["HCaloHemisphereMass"]->Fill(Groups[i].GetMass());
      QM1DHistograms["HCaloHemisphereEta"]->Fill(Groups[i].GetEta());
      QM1DHistograms["HCaloHemispherePhi"]->Fill(Groups[i].GetPhi());
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillPFHemisphere(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
   vector<FourVector> &Groups)
{
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillCaloRazor(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
   vector<FourVector> &Groups, int TCHELCount, int SSVHEMCount, int TCHETCount)
{
   FourVector CaloMET(0, pxMet[0], pyMet[0], 0);
   
   double MRStar = GetMRStar(Groups[0], Groups[1]);
   double RStar = GetRStar(Groups[0], Groups[1], CaloMET);
   double GammaRStar = GetGammaRStar(Groups[0], Groups[1]);

   QM1DHistograms["HCaloMRStar"]->Fill(MRStar);
   QM1DHistograms["HCaloRStar"]->Fill(RStar);
   QM2DHistograms["HCaloMRStarVsRStar"]->Fill(MRStar, RStar);
   QM1DHistograms["HCaloGammaRStar"]->Fill(GammaRStar);
   
   if(CaloMET.GetPT() > 60)
   {
      QM2DHistograms["HCaloMRStarVsRStarWithMETCut"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarWithMETCut"]->Fill(GammaRStar);
   }
   if(CaloMET.GetPT() > 120)
   {
      QM2DHistograms["HCaloMRStarVsRStarWithHigherMETCut"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarWithHigherMETCut"]->Fill(GammaRStar);
   }

   if(TCHELCount >= 1)
   {
      QM1DHistograms["HCaloMRStarTCHEL"]->Fill(MRStar);
      QM1DHistograms["HCaloRStarTCHEL"]->Fill(RStar);
      QM2DHistograms["HCaloMRStarVsRStarTCHEL"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarTCHEL"]->Fill(GammaRStar);
   }
   
   if(TCHELCount >= 2)
   {
      QM1DHistograms["HCaloMRStarTCHEL2"]->Fill(MRStar);
      QM1DHistograms["HCaloRStarTCHEL2"]->Fill(RStar);
      QM2DHistograms["HCaloMRStarVsRStarTCHEL2"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarTCHEL2"]->Fill(GammaRStar);
   }
   
   if(TCHETCount >= 1)
   {
      QM1DHistograms["HCaloMRStarTCHET"]->Fill(MRStar);
      QM1DHistograms["HCaloRStarTCHET"]->Fill(RStar);
      QM2DHistograms["HCaloMRStarVsRStarTCHET"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarTCHET"]->Fill(GammaRStar);
   }
   
   if(TCHETCount >= 2)
   {
      QM1DHistograms["HCaloMRStarTCHET2"]->Fill(MRStar);
      QM1DHistograms["HCaloRStarTCHET2"]->Fill(RStar);
      QM2DHistograms["HCaloMRStarVsRStarTCHET2"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarTCHET2"]->Fill(GammaRStar);
   }
   
   if(SSVHEMCount >= 1)
   {
      QM1DHistograms["HCaloMRStarSSVHEM"]->Fill(MRStar);
      QM1DHistograms["HCaloRStarSSVHEM"]->Fill(RStar);
      QM2DHistograms["HCaloMRStarVsRStarSSVHEM"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarSSVHEM"]->Fill(GammaRStar);
   }
   
   if(SSVHEMCount >= 2)
   {
      QM1DHistograms["HCaloMRStarSSVHEM2"]->Fill(MRStar);
      QM1DHistograms["HCaloRStarSSVHEM2"]->Fill(RStar);
      QM2DHistograms["HCaloMRStarVsRStarSSVHEM2"]->Fill(MRStar, RStar);
      QM1DHistograms["HCaloGammaRStarSSVHEM2"]->Fill(GammaRStar);
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::QMFillPFRazor(map<string, TH1D *> &QM1DHistograms, map<string, TH2D *> &QM2DHistograms,
   vector<FourVector> &Groups, int TCHELCount, int SSVHEMCount, int TCHETCount)
{
}
//---------------------------------------------------------------------------
bool LQ3Analysis::CheckCaloJetID(int Index)
{
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
bool LQ3Analysis::CheckPFJetID(int Index)
{
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
bool LQ3Analysis::CheckVertex(int Index)
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
bool LQ3Analysis::CheckHLTBit(string PathName)
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
vector<MuonCandidate> LQ3Analysis::MakeMuonCandidates()
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
      Candidate.CombinedIsolation = sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i];

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
vector<ElectronCandidate> LQ3Analysis::MakeElectronCandidates()
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
void LQ3Analysis::FillTopMuons(vector<MuonCandidate> &Candidates, LQ3OutputRecord &M)
{
   M.GoodMuonCountTight = 0;
   M.GoodMuonCountLoose = 0;

   multimap<double, int, greater<double> > GoodMuonIndices;
   for(int i = 0; i < (int)Candidates.size(); i++)
   {
      if(Candidates[i].PT < 10)
         continue;
   
      if(Candidates[i].PassMuonTight == true)
         M.GoodMuonCountTight = M.GoodMuonCountTight + 1;
      if(Candidates[i].PassMuonLoose == true)
         M.GoodMuonCountLoose = M.GoodMuonCountLoose + 1;

      if(Candidates[i].PassMuonLoose == false)
         continue;
      GoodMuonIndices.insert(pair<double, int>(Candidates[i].PT, i));
   }

   vector<int> GoodMuons;
   for(multimap<double, int, greater<double> >::iterator iter = GoodMuonIndices.begin();
      iter!= GoodMuonIndices.end(); iter++)
      GoodMuons.push_back(iter->second);

   M.GoodMuonCount = GoodMuons.size();
   if(GoodMuons.size() > 0)   M.Muons[0] = Candidates[GoodMuons[0]];
   if(GoodMuons.size() > 1)   M.Muons[1] = Candidates[GoodMuons[1]];
   if(GoodMuons.size() > 2)   M.Muons[2] = Candidates[GoodMuons[2]];
   if(GoodMuons.size() > 3)   M.Muons[3] = Candidates[GoodMuons[3]];
   if(GoodMuons.size() > 4)   M.Muons[4] = Candidates[GoodMuons[4]];
   if(GoodMuons.size() > 5)   M.Muons[5] = Candidates[GoodMuons[5]];
   if(GoodMuons.size() > 6)   M.Muons[6] = Candidates[GoodMuons[6]];
   if(GoodMuons.size() > 7)   M.Muons[7] = Candidates[GoodMuons[7]];
   if(GoodMuons.size() > 8)   M.Muons[8] = Candidates[GoodMuons[8]];
   if(GoodMuons.size() > 9)   M.Muons[9] = Candidates[GoodMuons[9]];
}
//---------------------------------------------------------------------------
void LQ3Analysis::FillTopElectrons(vector<ElectronCandidate> &Candidates, LQ3OutputRecord &M)
{
   M.GoodElectronCount95 = 0;
   M.GoodElectronCount90 = 0;
   M.GoodElectronCount85 = 0;
   M.GoodElectronCount80 = 0;

   multimap<double, int, greater<double> > GoodElectronIndices;
   for(int i = 0; i < (int)Candidates.size(); i++)
   {
      if(Candidates[i].PT < 20)
         continue;
      if(Candidates[i].SuperClusterEnergy / cosh(Candidates[i].SuperClusterEta) < 20)
         continue;

      if(Candidates[i].PassWP95 == true)
         M.GoodElectronCount95 = M.GoodElectronCount95 + 1;
      if(Candidates[i].PassWP90 == true)
         M.GoodElectronCount90 = M.GoodElectronCount90 + 1;
      if(Candidates[i].PassWP85 == true)
         M.GoodElectronCount85 = M.GoodElectronCount85 + 1;
      if(Candidates[i].PassWP80 == true)
         M.GoodElectronCount80 = M.GoodElectronCount80 + 1;

      if(Candidates[i].PassWP95 == false)
         continue;
      GoodElectronIndices.insert(pair<double, int>(Candidates[i].PT, i));
   }

   vector<int> GoodElectrons;
   for(multimap<double, int>::iterator iter = GoodElectronIndices.begin();
      iter!= GoodElectronIndices.end(); iter++)
      GoodElectrons.push_back(iter->second);

   M.GoodElectronCount = GoodElectrons.size();
   if(GoodElectrons.size() > 0)   M.Electrons[0] = Candidates[GoodElectrons[0]];
   if(GoodElectrons.size() > 1)   M.Electrons[1] = Candidates[GoodElectrons[1]];
   if(GoodElectrons.size() > 2)   M.Electrons[2] = Candidates[GoodElectrons[2]];
   if(GoodElectrons.size() > 3)   M.Electrons[3] = Candidates[GoodElectrons[3]];
   if(GoodElectrons.size() > 4)   M.Electrons[4] = Candidates[GoodElectrons[4]];
   if(GoodElectrons.size() > 5)   M.Electrons[5] = Candidates[GoodElectrons[5]];
   if(GoodElectrons.size() > 6)   M.Electrons[6] = Candidates[GoodElectrons[6]];
   if(GoodElectrons.size() > 7)   M.Electrons[7] = Candidates[GoodElectrons[7]];
   if(GoodElectrons.size() > 8)   M.Electrons[8] = Candidates[GoodElectrons[8]];
   if(GoodElectrons.size() > 9)   M.Electrons[9] = Candidates[GoodElectrons[9]];
}
//---------------------------------------------------------------------------
void LQ3Analysis::FillCaloJet(LQ3OutputRecord &M)
{
   multimap<double, int, greater<double> > SortCaloJets;   // sort by pt
   for(int i = 0; i < nAK5Jet; i++)
   {
      // if(CheckCaloJetID(i) == false)   // bad jet, ignore
      //    continue;
      if(etaAK5Jet[i] < -3 || etaAK5Jet[i] > 3)
         continue;

      double JetPT = sqrt(pxAK5Jet[i] * pxAK5Jet[i] + pyAK5Jet[i] * pyAK5Jet[i]);
      SortCaloJets.insert(pair<double, int>(JetPT, i));

      M.CaloJetCount = M.CaloJetCount + 1;
      if(JetPT > 30)    M.CaloJetCount30 = M.CaloJetCount30 + 1;
      if(JetPT > 50)    M.CaloJetCount50 = M.CaloJetCount50 + 1;
      if(JetPT > 100)   M.CaloJetCount100 = M.CaloJetCount100 + 1;
   }
   int CurrentCount = 0;
   for(multimap<double, int>::iterator iter = SortCaloJets.begin();
         iter != SortCaloJets.end(); iter++)
   {
      if(CurrentCount >= 50)
         break;
      M.CaloJetE[CurrentCount] = energyAK5Jet[iter->second];
      M.CaloJetPT[CurrentCount] = iter->first;
      M.CaloJetEta[CurrentCount] = etaAK5Jet[iter->second];
      M.CaloJetPhi[CurrentCount] = phiAK5Jet[iter->second];
      M.CaloJetCSVTag[CurrentCount] =
         combinedSecondaryVertexBJetTagsAK5Jet[iter->second];
      M.CaloJetCSVMTag[CurrentCount] =
         combinedSecondaryVertexMVABJetTagsAK5Jet[iter->second];
      M.CaloJetTCHPTag[CurrentCount] =
         trackCountingHighPurBJetTagsAK5Jet[iter->second];
      M.CaloJetTCHETag[CurrentCount] =
         trackCountingHighEffBJetTagsAK5Jet[iter->second];
      M.CaloJetProbabilityTag[CurrentCount] =
         jetProbabilityBJetTagsAK5Jet[iter->second];
      M.CaloJetBProbabilityTag[CurrentCount] =
         jetBProbabilityBJetTagsAK5Jet[iter->second];
      M.CaloJetSSVHETag[CurrentCount] =
         simpleSecondaryVertexHighEffBJetTagsAK5Jet[iter->second];
      M.CaloJetSSVHPTag[CurrentCount] =
         simpleSecondaryVertexHighPurBJetTagsAK5Jet[iter->second];
      CurrentCount = CurrentCount + 1;
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::FillPFJet(LQ3OutputRecord &M)
{
   multimap<double, int, greater<double> > SortPFJets;   // sort by pt
   for(int i = 0; i < nAK5PFPUcorrJet; i++)
   {
      // if(CheckPFJetID(i) == false)   // bad jet, ignore
      //    continue;
      if(etaAK5PFPUcorrJet[i] < -2.4 || etaAK5PFPUcorrJet[i] > 2.4)
         continue;

      double JetPT =
         sqrt(pxAK5PFPUcorrJet[i] * pxAK5PFPUcorrJet[i] + pyAK5PFPUcorrJet[i] * pyAK5PFPUcorrJet[i]);
      SortPFJets.insert(pair<double, int>(JetPT, i));

      M.PFJetCount = M.PFJetCount + 1;
      if(JetPT > 30)    M.PFJetCount30 = M.PFJetCount30 + 1;
      if(JetPT > 50)    M.PFJetCount50 = M.PFJetCount50 + 1;
      if(JetPT > 100)   M.PFJetCount100 = M.PFJetCount100 + 1;
   }
   int CurrentCount = 0;
   for(multimap<double, int>::iterator iter = SortPFJets.begin();
         iter != SortPFJets.end(); iter++)
   {
      if(CurrentCount >= 50)
         break;
      M.PFJetE[CurrentCount] = energyAK5PFPUcorrJet[iter->second];
      M.PFJetPT[CurrentCount] = iter->first;
      M.PFJetEta[CurrentCount] = etaAK5PFPUcorrJet[iter->second];
      M.PFJetPhi[CurrentCount] = phiAK5PFPUcorrJet[iter->second];
      M.PFJetCSVTag[CurrentCount] =
         combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetCSVMTag[CurrentCount] =
         combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetTCHPTag[CurrentCount] =
         trackCountingHighPurBJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetTCHETag[CurrentCount] =
         trackCountingHighEffBJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetProbabilityTag[CurrentCount] =
         jetProbabilityBJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetBProbabilityTag[CurrentCount] =
         jetBProbabilityBJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetSSVHETag[CurrentCount] =
         simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[iter->second];
      M.PFJetSSVHPTag[CurrentCount] =
         simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[iter->second];
      CurrentCount = CurrentCount + 1;
   }
}
//---------------------------------------------------------------------------
void LQ3Analysis::FillPV(LQ3OutputRecord &M)
{
   M.PrimaryVertexCount = 0;
   M.PrimaryVertexMaxSumPT = 0;
   
   if(nPV > 1000)
      return;

   int MaxSumPTPVIndex = -1;
   for(int i = 0; i < nPV && i < 20; i++)
   {
      if(CheckVertex(i) == false)
         continue;
      if(MaxSumPTPVIndex < 0 || SumPtPV[MaxSumPTPVIndex] < SumPtPV[i])
         MaxSumPTPVIndex = i;
      M.PrimaryVertexCount = M.PrimaryVertexCount + 1;
   }
   if(MaxSumPTPVIndex >= 0)
      M.PrimaryVertexMaxSumPT = SumPtPV[MaxSumPTPVIndex];
}
//---------------------------------------------------------------------------
void LQ3Analysis::FillMET(LQ3OutputRecord &M)
{
   M.CaloMET[0] = pxMet[0];
   M.CaloMET[1] = pyMet[0];
   M.PFMET[0] = pxPFMet[0];
   M.PFMET[1] = pyPFMet[0];
}
//---------------------------------------------------------------------------












