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
// $Id: DataScoutingAnalyzer.cc,v 1.5 2012/10/04 13:16:44 amott Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "amott/DataScoutingAnalyzer/interface/DataScoutingAnalyzer.h"

//objects
//#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/METReco/interface/CaloMetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
template <typename jettype, typename mettype>
DataScoutingAnalyzer<jettype,mettype>::DataScoutingAnalyzer(const edm::ParameterSet& iConfig):
  tag_recoJet(iConfig.getParameter<edm::InputTag>("jets")),
  s_recoJetCorrector(iConfig.getParameter<std::string>("jetCorrections")),
  tag_recoRho(iConfig.getParameter<edm::InputTag>("rho")),  
  jetThreshold(iConfig.getParameter<double>("jetThreshold")),
  tag_recoMet(iConfig.getParameter<edm::InputTag>("met")),
  tag_recoElectrons(iConfig.getParameter<edm::InputTag>("electrons")),
  tag_recoMuons(iConfig.getParameter<edm::InputTag>("muons")),
  tag_hcalNoise(iConfig.getParameter<edm::InputTag>("noise")),
  s_outputFile(iConfig.getParameter<std::string>("outputFile"))
{
   //now do what ever initialization is needed

}



template <typename jettype, typename mettype>
DataScoutingAnalyzer<jettype,mettype>::~DataScoutingAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
template <typename jettype, typename mettype>
void
DataScoutingAnalyzer<jettype,mettype>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<std::vector<jettype> > h_recoJet;
  const JetCorrector* corrector = JetCorrector::getJetCorrector (s_recoJetCorrector, iSetup);  
  //Handle<std::vector<reco::Jet> > h_recoJet;
  //Handle<std::vector<reco::MET> > h_recoMet;
  Handle<std::vector<mettype> > h_recoMet;
  Handle<double> h_recoRho;

  //Handle<reco::Electron> h_recoElectrons;
  //Handle<reco::Muon> h_recoMuons;
  
  iEvent.getByLabel(tag_recoJet,h_recoJet);
  iEvent.getByLabel(tag_recoMet,h_recoMet);
  iEvent.getByLabel(tag_recoRho,h_recoRho);
  //iEvent.getByLabel(tag_recoElectrons,h_recoElectrons);
  //iEvent.getByLabel(tag_recoMuons,h_recoMuons);

  Handle<std::vector<reco::CaloJet> > h_dsJet;
  Handle<std::vector<reco::CaloMET> > h_dsMet;
  Handle<std::vector<reco::CaloMET> > h_dsMetClean;
  Handle<double> h_dsRho;
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
  iEvent.getByLabel("hltKT6CaloJets","rho",h_dsRho);
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

  //MET FILTERS
  edm::Handle< bool > ECALTPFilter;
  try { iEvent.getByLabel("EcalDeadCellEventFlagProducer", ECALTPFilter); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << ECALTPFilter; }
  ECALTPFilterFlag = *ECALTPFilter;

  edm::Handle< bool > HBHENoiseFilterResult;
  try { iEvent.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), HBHENoiseFilterResult); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << HBHENoiseFilterResult; }
  HBHENoiseFilterResultFlag = *HBHENoiseFilterResult;

  edm::Handle< bool > hcalLaserEventFilter;
  try { iEvent.getByLabel("hcalLaserEventFilter", hcalLaserEventFilter); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << hcalLaserEventFilter; }
  hcalLaserEventFilterFlag = *hcalLaserEventFilter;

  edm::Handle< bool > eeBadScFilter;
  try { iEvent.getByLabel("eeBadScFilter", eeBadScFilter); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get bool: " << eeBadScFilter; }
  eeBadScFilterFlag = *eeBadScFilter;

  edm::Handle< int > ECALDeadDRFilter;
  try { iEvent.getByLabel("simpleDRFlagProducer","deadCellStatus", ECALDeadDRFilter); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get int: " << ECALDeadDRFilter; }
  ECALDeadDRFilterFlag = *ECALDeadDRFilter;
    
  edm::Handle< int > ECALBoundaryDRFilter;
  try { iEvent.getByLabel("simpleDRFlagProducer","boundaryStatus", ECALBoundaryDRFilter); }
  catch ( cms::Exception& ex ) { edm::LogWarning("CmsMetFiller") << "Can't get int: " << ECALBoundaryDRFilter; }
  ECALBoundaryDRFilterFlag = *ECALBoundaryDRFilter;
  
  //fill the tree
  
  //MET
  dsMetPt = h_dsMet->front().pt();
  dsMetPhi = h_dsMet->front().phi();
  dsMetCleanPt = h_dsMetClean->front().pt();
  dsMetCleanPhi = h_dsMetClean->front().phi();
  recoMetPt = h_recoMet->front().pt();
  recoMetPhi = h_recoMet->front().phi();
  
  dsRho = *h_dsRho;
  recoRho = *h_recoRho;
  
  typename std::vector<jettype>::const_iterator i_recoJet;
  //std::vector<reco::Jet>::const_iterator i_recoJet;
  nRECOJets=0;
  for(i_recoJet = h_recoJet->begin(); i_recoJet != h_recoJet->end(); i_recoJet++){
    double scale = corrector->correction(*i_recoJet,iEvent,iSetup);
    ((reco::Jet)*i_recoJet).scaleEnergy(scale);
    if(i_recoJet->pt() < jetThreshold - 10) continue;
    recoJetPt[nRECOJets] = i_recoJet->pt();
    recoJetEta[nRECOJets] = i_recoJet->eta();
    recoJetPhi[nRECOJets] = i_recoJet->phi();
    recoJetE[nRECOJets] = i_recoJet->energy();
    nRECOJets++;
  }

  std::vector<reco::CaloJet>::const_iterator i_dsJet;
  nDSJets=0;
  for(i_dsJet = h_dsJet->begin(); i_dsJet != h_dsJet->end(); i_dsJet++){
    //apply pileup correction
    float pileupCorr = i_dsJet->jetArea()*dsRho;

    if(i_dsJet->pt()-pileupCorr < jetThreshold) continue;
    dsJetPt[nDSJets] = i_dsJet->pt()-pileupCorr;
    dsJetEta[nDSJets] = i_dsJet->eta();
    dsJetPhi[nDSJets] = i_dsJet->phi();
    dsJetE[nDSJets] = i_dsJet->energy()-pileupCorr*cosh(i_dsJet->eta());
    dsJetFracHad[nDSJets] = i_dsJet->energyFractionHadronic();

    //do jet matching
    float bestdEoE = 9999;
    int bestIndex=-1;
    for( int iRECOJet=0; iRECOJet < nRECOJets; iRECOJet++){
      //std::cout << "DeltaR: " << reco::deltaR(i_dsJet->eta(),i_dsJet->phi(),recoJetEta[iRECOJet],recoJetPhi[iRECOJet]) << std::endl;
      if( reco::deltaR(i_dsJet->eta(),i_dsJet->phi(),recoJetEta[iRECOJet],recoJetPhi[iRECOJet])> 0.5) continue; //require DR match
      
      float dEoE = fabs(i_dsJet->energy() - recoJetE[iRECOJet])/recoJetE[iRECOJet];
      //std::cout << "dEoE: " << dEoE << std::endl;
      if(dEoE < 0.5 && dEoE < bestdEoE){
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
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::beginJob()
{
  outputFile = new TFile(s_outputFile.c_str(),"RECREATE");
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

  outputTree->Branch("ECALTPFilterFlag",&ECALTPFilterFlag,"ECALTPFilterFlag/B");
  outputTree->Branch("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlag,"HBHENoiseFilterResultFlag/B");
  outputTree->Branch("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlag,"hcalLaserEventFilterFlag/B");
  outputTree->Branch("eeBadScFilterFlag",&eeBadScFilterFlag,"eeBadScFilterFlag/B");
  outputTree->Branch("ECALDeadDRFilterFlag",&ECALDeadDRFilterFlag,"ECALDeadDRFilterFlag/I");
  outputTree->Branch("ECALBoundaryDRFilterFlag",&ECALBoundaryDRFilterFlag,"ECALBoundaryDRFilterFlag/I");

}

// ------------ method called once each job just after ending the event loop  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::endJob() 
{
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();
}

// ------------ method called when starting to processes a run  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
template <typename jettype, typename mettype>
void 
DataScoutingAnalyzer<jettype,mettype>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <typename jettype, typename mettype>
void
DataScoutingAnalyzer<jettype,mettype>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


typedef DataScoutingAnalyzer<reco::CaloJet,reco::CaloMET> CaloScoutingAnalyzer;
typedef DataScoutingAnalyzer<reco::PFJet,reco::CaloMET> PFJetScoutingAnalyzer;
typedef DataScoutingAnalyzer<reco::PFJet,reco::PFMET> PFScoutingAnalyzer;
//define this as a plug-in
DEFINE_FWK_MODULE(CaloScoutingAnalyzer);
DEFINE_FWK_MODULE(PFJetScoutingAnalyzer);
DEFINE_FWK_MODULE(PFScoutingAnalyzer);
