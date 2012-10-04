#ifndef DataScoutingAnalyzer_h
#define DataScoutingAnalyzer_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <TTree.h>
#include <TFile.h>

//
// class declaration
//
template <typename jettype, typename mettype>
class DataScoutingAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DataScoutingAnalyzer(const edm::ParameterSet&);
      ~DataScoutingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  
  edm::InputTag tag_recoJet;
  std::string s_recoJetCorrector;
  edm::InputTag tag_recoRho;
  edm::InputTag tag_recoMet;
  edm::InputTag tag_recoElectrons;
  edm::InputTag tag_recoMuons;
  edm::InputTag tag_hcalNoise;
  std::string s_outputFile;

  //Optional MET filters
  // 1 means we KEEP the event
  bool ECALTPFilterFlag;
  bool HBHENoiseFilterResultFlag;
  bool hcalLaserEventFilterFlag;
  bool eeBadScFilterFlag;
  int  ECALDeadDRFilterFlag;
  int  ECALBoundaryDRFilterFlag;
  
  double jetThreshold;

  TFile *outputFile;
  TTree *outputTree;

  static const int _kMaxJets = 200;

  int nDSJets;
  float dsJetPt[_kMaxJets];
  float dsJetEta[_kMaxJets];
  float dsJetPhi[_kMaxJets];
  float dsJetE[_kMaxJets];
  float dsJetFracHad[_kMaxJets];
  int dsJetMatchIndex[_kMaxJets];


  float dsRho;

  float dsMetPt;
  float dsMetPhi;
  float dsMetCleanPt;
  float dsMetCleanPhi;
  
  int nRECOJets;
  float recoJetPt[_kMaxJets];
  float recoJetEta[_kMaxJets];
  float recoJetPhi[_kMaxJets];
  float recoJetE[_kMaxJets];
  //  float recoJetFracHad[_kMaxJets];

  float recoRho;

  float recoMetPt;
  float recoMetPhi;
  float recoMetCleanPt;
  float recoMetCleanPhi;
  

};


#endif
