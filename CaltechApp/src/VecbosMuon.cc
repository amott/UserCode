#include "VecbosMu.hh"

VecbosMu::VecbosMu(){}

VecbosMu::VecbosMu(VecbosBase* o, int i){
  this->Init(o,i);
}

void VecbosMu::Init(VecbosBase* o, int i){
  if(i<0 || i > o->nMuon){
    index = -1;
  }else{
    index = i;            
    energy = o->energyMuon[i];         
    pt = TMath::Sqrt( TMath::Power(o->pxMuon[i],2) + TMath::Power(o->pyMuon[i],2) );
    eta = o->etaMuon[i];
    phi = o->phiMuon[i];
    //p4.SetPtEtaPhiE(pt,eta,phi,energy);
    charge = o->chargeMuon[i];
    combinedIso = (o->emEt03Muon[i] + o->hadEt03Muon[i] + o->sumPt03Muon[i] - o->rhoFastjet * TMath::Pi()*0.3*0.3)/pt;;
    emIso = o->emEt03Muon[i];
    hadIso = o->hadEt03Muon[i];
    trkIso = o->sumPt03Muon[i];
    Utils AnalysisUtilities;
    isGlobalMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::AllGlobalMuons);
    isTrackerMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::AllTrackerMuons);
    isPromptMuon = AnalysisUtilities.muonIdVal(o->muonIdMuon[i], bits::GlobalMuonPromptTight);
    int iTrack = o->trackIndexMuon[i];
    if(iTrack<0 || iTrack > o->nTrack) {
      nTrackHits=0;
      nPixelHits=0;
      trackImpactPar=0;
    }else{
      nTrackHits = o->numberOfValidStripTIBHitsTrack[iTrack]
	+ o->numberOfValidStripTIDHitsTrack[iTrack]
	+ o->numberOfValidStripTOBHitsTrack[iTrack]
	+ o->numberOfValidStripTECHitsTrack[iTrack];
      nPixelHits = o->numberOfValidPixelBarrelHitsTrack[iTrack] + o->numberOfValidPixelEndcapHitsTrack[iTrack];
      trackImpactPar = fabs(o->transvImpactParTrack[iTrack]);
    }
    isLooseMuon = true;
    if(!isGlobalMuon || nTrackHits<=10) isLooseMuon=false;
    isTightMuon = isLooseMuon;
    if(!isTrackerMuon || !isPromptMuon || combinedIso >=0.15 || nPixelHits == 0 || trackImpactPar >=0.2) isTightMuon = false;
  }

  vtxX = o->vertexXMuon[i];
  vtxY = o->vertexYMuon[i];
  vtxZ = o->vertexZMuon[i];

  genMatch.Init(o,-1);
}
