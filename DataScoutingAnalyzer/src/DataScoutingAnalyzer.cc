// -*- C++ -*-
//
// Package:    DataScoutingAnalyzer
// Class:      DataScoutingAnalyzer
// 
/**\class DataScoutingAnalyzer DataScoutingAnalyzer.cc tmp/DataScoutingAnalyzer/src/DataScoutingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alex Mott
//         Created:  Wed Sep  5 16:25:41 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataScoutingAnalyzer.h"

//objects
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloMetCollection.h"
#include "DataFormats/JetReco/interface/CaloMet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/Common/interface/AssociationMap.h"

#incluce "DataFormats/Math/interface/deltaR.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DataScoutingAnalyzer::DataScoutingAnalyzer(const edm::ParameterSet& iConfig):
  tag_recoJet(iConfig.getParameter<edm::InputTag>("jets")),
  tag_recoRho(iConfig.getParameter<edm::InputTag>("rho")),  
  jetMatchThreshold(iConfig.getParameter<double>("jetMatchingThreshold")),
  tag_recoMet(iConfig.getParameter<edm::InputTag>("met")),
  tag_recoElectrons(iConfig.getParameter<edm::InputTag>("electrons")),
  tag_recoMuons(iConfig.getParameter<edm::InputTag>("muons")),
  tag_hcalNoise(iConfig.getParameter<edm::InputTag>("noise")),
  s_outputFile(iConfig.getParameter<string>("outputFile"))
{
   //now do what ever initialization is needed

}


DataScoutingAnalyzer::~DataScoutingAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DataScoutingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::JetCollection> h_recoJet;
  Handle<reco::MetCollection> h_recoMet;
  Handle<double> h_recoRho;

  //Handle<reco::Electron> h_recoElectrons;
  //Handle<reco::Muon> h_recoMuons;
  
  iEvent.getByLabel(tag_recoJet,h_recoJet);
  iEvent.getByLabel(tag_recoMet,h_recoMet);
  iEvent.getByLabel(tag_recoRho,h_recoRho);
  //iEvent.getByLabel(tag_recoElectrons,h_recoElectrons);
  //iEvent.getByLabel(tag_recoMuons,h_recoMuons);

  Handle<reco::CaloJetCollection> h_dsJet;
  Handle<reco::CaloMetCollection> h_dsMet;
  Handle<reco::CaloMetCollection> h_dsMetClean;
  Handle<double> h_rho;
  /*
  Handle<reco::Electron> h_dsElectrons;
  Handle<reco::SuperCluster> h_dsSC;
  Handle<reco::Muon> h_dsMuons;

  Handle<reco::ElectronIsolationMap> trackIsoMap;
  Handle<reco::ElectronIsolationMap> trackdEtaMap;
  Handle<reco::ElectronIsolationMap> trackdPhiMap;

  Handle<reco::RecoEcalCandidateIsolationMap> sieieMap;
  Handle<reco::RecoEcalCandidateIsolationMap> ecalIsoMap;
  Handle<reco::RecoEcalCandidateIsolationMap> hcalIsoMap;
  Handle<reco::RecoEcalCandidateIsolationMap> hforHEMap;
  */  
  iEvent.getByLabel("hltCaloJetIDPassed",h_dsJet);
  iEvent.getByLabel("hltMet",h_dsMet);
  iEvent.getByLabel("hltMetClean",h_dsMetClean);
  iEvent.getByLabel("hltKT6CaloJets","rhos",h_dsRho);
  /*
  iEvent.getByLabel("hltPixelMatchElectronsActivity",h_dsElectrons);
  iEvent.getByLabel("hltRecoEcalSuperClusterActivityCandidate",h_dsSC);
  iEvent.getByLabel("hltL3MuonCandidates",h_dsMuons);
  
  iEvent.getByLabel("hltHitElectronActivityTrackIsol",trackIsoMap);
  iEvent.getByLabel("hltHitElectronActivityDetaDphi","Deta",trackEtaMap);
  iEvent.getByLabel("hltHitElectronActivityTrackIsol","Dphi",trackPhiMap);

  iEvent.getByLabel("hltActivityPhotonClusterShape",sieieMap);
  iEvent.getByLabel("hltActivityPhotonEcalIso",ecalIsoMap);
  iEvent.getByLabel("hltActivityPhotonHcalIso",hcalIsoMap);
  iEvent.getByLabel("hltActivityPhotonHcalForHE",hforHEMap);
  */
  
  //fill the tree
  
  //MET
  dsMetPt = h_dsMet->Front().pt();
  dsMetPhi = h_dsMet->Front().phi();
  dsMetCleanPt = h_dsMetClean->Front().pt();
  dsMetCleanPhi = h_dsMetClean->Front().phi();
  recoMetPt = h_recoMet->Front().pt();
  recoMetPhi = h_recoMet->Front().phi();
  
  dsRho = *h_dsRho;
  recoRho = *h_recoRho;
  
  reco::JetCollection::const_iterator i_recoJet;
  nRECOJets=0;
  for(i_recoJet = h_recoJet->begin(); i_recoJet != h_recoJet->end(); irecoJet++){
    recoJetPt[nRECOJets] = i_recoJet->pt();
    recoJetEta[nRECOJets] = i_recoJet->eta();
    recoJetPhi[nRECOJets] = i_recoJet->phi();
    recoJetE[nRECOJets] = i_recoJet->energy();
    nRECOJets++;
  }

  reco::CaloJetCollection::const_iterator i_dsJet;
  nDSJets=0;
  for(i_dsJet = h_dsJet->begin(); i_dsJet != h_dsJet->end(); idsJet++){
    dsJetPt[nDSJets] = i_dsJet->pt();
    dsJetEta[nDSJets] = i_dsJet->eta();
    dsJetPhi[nDSJets] = i_dsJet->phi();
    dsJetE[nDSJets] = i_dsJet->energy();
    dsJetFracHad[nDSJets] = i_dsJet->energyFractionHadronic();

    //do jet matching
    float bestdEoE = 9999;
    int bestIndex=-1;
    for( int iRECOJet=0; iRECOJet < nRECOJets; iRECOJet++){
      if( reco::deltaR(i_dsJet->eta(),i_dsJet->phi(),recoJetEta[iRECOJet],recoJetPhi[iRECOJet])
	  > 0.5) continue; //require DR match
      
      float dEoE = (i_dsJet->energy() - recoJetE[iRECOJet])/recoJetE[iRECOJet];
      if(dEoE < bestdEoE){
	bestdEoE = dEoE;
	bestIndex = iRECOJet;
      }
    }
    dsJetMatchIndex[nDSJets] = bestIndex;
    nDSJets++;
  }

  outputTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
DataScoutingAnalyzer::beginJob()
{
  outputFile = new TFile(s_outputFile,"RECREATE");
  outputTree = new TTree("DSComp","");

  outputTree->Branch("nDSJets",&nDSJets,"nDSJets/I");
  outputTree->Branch("dsJetPt",dsJetPt,"dsJetPt[nDSJets]");
  outputTree->Branch("dsJetEta",dsJetEta,"dsJetEta[nDSJets]");
  outputTree->Branch("dsJetPhi",dsJetPhi,"dsJetPhi[nDSJets]");
  outputTree->Branch("dsJetE",dsJetE,"dsJetE[nDSJets]");
  outputTree->Branch("dsJetFracHad",dsJetFracHad,"dsJetFracHad[nDSJets]");
  outputTree->Branch("dsJetMatchIndex",dsJetMatchIndex,"dsJetMatchIndex[nDSJets]/I");

  outputTree->Branch("dsRho",&dsRho);
  outputTree->Branch("dsMetPt",&dsMetPt);
  outputTree->Branch("dsMetPhi",&dsMetPhi);
  outputTree->Branch("dsMetCleanPt",&dsMetCleanPt);
  outputTree->Branch("dsMetCleanPhi",&dsMetCleanPhi);

  outputTree->Branch("nRECOJets",&nRECOJets,"nRECOJets/I");
  outputTree->Branch("recoJetPt",recoJetPt,"recoJetPt[nRECOJets]");
  outputTree->Branch("recoJetEta",recoJetEta,"recoJetEta[nRECOJets]");
  outputTree->Branch("recoJetPhi",recoJetPhi,"recoJetPhi[nRECOJets]");
  outputTree->Branch("recoJetE",recoJetE,"recoJetE[nRECOJets]");
  outputTree->Branch("recoJetE",recoJetE,"recoJetE[nRECOJets]");

  outputTree->Branch("recoRho",&recoRho);
  outputTree->Branch("recoMetPt",&recoMetPt);
  outputTree->Branch("recoMetPhi",&recoMetPhi);
  outputTree->Branch("recoMetCleanPt",&recoMetCleanPt);
  outputTree->Branch("recoMetCleanPhi",&recoMetCleanPhi);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DataScoutingAnalyzer::endJob() 
{
  outputTree->Write();
  outputFile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
DataScoutingAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DataScoutingAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DataScoutingAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DataScoutingAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DataScoutingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DataScoutingAnalyzer);
