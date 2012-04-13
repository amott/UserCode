#include "HggEGEnergyCorrector.hh"
#include "VecbosEGObject.hh"
#include "HggPhysUtils.cc"
#include <string>

HggEGEnergyCorrector::HggEGEnergyCorrector(VecbosBase *r,int phCorrs,Bool_t isRealData):
  base(r)
{
  this->Init(phCorrs,isRealData);
}


void HggEGEnergyCorrector::Init(int testPhotonCorr,Bool_t isRealData){
  //  fVals = new Float_t[18];
  fVals = new Float_t[73];
  
  std::string regweights; 
  if( testPhotonCorr == 99){
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrph.root";
  }else if( testPhotonCorr == 98) {
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrele.root";
  }else if( testPhotonCorr == 97){
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrv2ele.root";
  }else if( testPhotonCorr == 96){
    regweights = "/afs/cern.ch/cms/cit/yongy/regweights/gbrv2ph.root";
  }
  
  TFile *fgbr = new TFile(regweights.c_str(),"READ");
  fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
  fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");
  
  fReaderee = (GBRForest*)fgbr->Get("EECorrection");
  fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");
  fgbr->Close();

  //load the ECAL geometry information
  ecalGeometry = loadecalGapCoordinates(isRealData);
}



std::pair<double,double> HggEGEnergyCorrector::CorrectedEnergyWithError(int j){
  VecbosPho pho(base,j); // initialize the photon with index j

  int iSC = base->superClusterIndexPho[j]; // supercluster index for this photon

  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry,pho.SC.eta,pho.SC.phi);
  
  fVals[0]  = pho.SC.rawE;
  fVals[1]  = pho.SC.r9();
  fVals[2]  = pho.SC.eta;
  fVals[3]  = pho.SC.phi;
  fVals[4]  = pho.SC.e5x5/pho.SC.rawE;
  
  Bool_t isbarrel = fabs(pho.SC.eta) < 1.48; 
  
  if( isbarrel){
    fVals[5]  = thisGeometry._aC; 
    fVals[6]  = thisGeometry._aS;
    fVals[7]  = thisGeometry._aM;
    fVals[8]  = thisGeometry._bC;
    fVals[9]  = thisGeometry._bS;
    fVals[10] = thisGeometry._bM;
    fVals[11] = pho.HoverE;
    fVals[12] = pho.SC.etaWidth;
    fVals[13] = pho.SC.phiWidth;
    fVals[14] = pho.SC.sigmaIEtaIEta;
  }else{
    fVals[5]  = pho.SC.rawE; //photonscpreshowerEnergy[j] / photonscrawEnergy[j]; /// DO WE STORE THIS????
    fVals[6]  = thisGeometry._xZ;
    fVals[7]  = thisGeometry._aC;
    fVals[8]  = thisGeometry._aS;
    fVals[9]  = thisGeometry._aM;
    fVals[10] = thisGeometry._yZ;
    fVals[11] = thisGeometry._bC;
    fVals[12] = thisGeometry._bS;
    fVals[13] = thisGeometry._bM;
    fVals[14] = pho.SC.HoverE;
    fVals[15] = pho.SC.etaWidth;
    fVals[16] = pho.SC.phiWidth;
    fVals[17] = pho.SC.sigmaIEtaIEta;
  }

  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = pho.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = pho.SC.rawE; //NEEED PRESHOWER!!! photonscrawEnergy[j] + photonscpreshowerEnergy[j]; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}



std::pair<double,double> HggEGEnergyCorrector::electronEnergyCorrector_CorrectedEnergyWithError(int j){
  VecbosEle ele(base,j);
    
  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry,ele.SC.eta,ele.SC.phi);
  
  fVals[0]  = ele.SC.rawE;
  fVals[1]  = ele.SC.r9();
  fVals[2]  = ele.SC.eta;
  fVals[3]  = ele.SC.phi;
  fVals[4]  = ele.SC.e5x5/ele.SC.rawE;
  
  Bool_t isbarrel = fabs(ele.SC.eta) < 1.48; 
  
  if( isbarrel){
    fVals[5]  = thisGeometry._aC; 
    fVals[6]  = thisGeometry._aS;
    fVals[7]  = thisGeometry._aM;
    fVals[8]  = thisGeometry._bC;
    fVals[9]  = thisGeometry._bS;
    fVals[10] = thisGeometry._bM;
    fVals[11] = ele.SC.HoverE;
    fVals[12] = ele.SC.etaWidth;
    fVals[13] = ele.SC.phiWidth;
    fVals[14] = ele.SC.sigmaIEtaIEta;

  }else{
    fVals[5]  = ele.esEnergy/ele.SC.rawE;
    fVals[6]  = thisGeometry._xZ;
    fVals[7]  = thisGeometry._aC;
    fVals[8]  = thisGeometry._aS;
    fVals[9]  = thisGeometry._aM;
    fVals[10] = thisGeometry._yZ;
    fVals[11] = thisGeometry._bC;
    fVals[12] = thisGeometry._bS;
    fVals[13] = thisGeometry._bM;
    fVals[14] = ele.SC.HoverE;
    fVals[15] = ele.SC.etaWidth;
    fVals[16] = ele.SC.phiWidth;
    fVals[17] = ele.SC.sigmaIEtaIEta;
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = ele.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = ele.SC.rawE + ele.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}


std::pair<double,double> HggEGEnergyCorrector::electronEnergyCorrector_CorrectedEnergyWithErrorv2(int j){
  VecbosEle ele(base,j);

  ///no correction for electron of tracker-deriven seed only ( not used in analysis anyway)
  if( !ele.isEcalDriven ){
    return std::pair<double,double>(ele.SC.energy,0);
  }

  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry, ele.SC.eta, ele.SC.phi);
  
 

  Bool_t isbarrel = fabs(ele.SC.eta) < 1.48; 
  
  if( isbarrel){
    

    fVals[0]  = ele.SC.rawE;
    fVals[1]  = ele.SC.r9();
    fVals[2]  = ele.SC.eta;
    fVals[3]  = ele.SC.phi;
    fVals[4]  = ele.SC.e5x5/ele.SC.rawE;
    fVals[5] = ele.HoverE;
    fVals[6] = ele.SC.etaWidth;
    fVals[7] = ele.SC.phiWidth;
    ///bc
    VecbosBC BC1 = ele.SC.basicClusters[0];
    
    fVals[8]  = BC1.eta - ele.SC.eta;
    fVals[9]  = DeltaPhi(BC1.phi,ele.SC.phi);
    fVals[10]  = BC1.energy/ele.SC.rawE;
    fVals[11]  = BC1.e3x3/BC1.energy;
    fVals[12]  = BC1.e5x5/BC1.energy;
    fVals[13] = sqrt(BC1.sigmaIEtaIEta);
    fVals[14] = sqrt(BC1.sigmaIPhiIPhi);
    fVals[15] = BC1.sigmaIEtaIPhi;
    fVals[16] = BC1.eMax/BC1.energy;
    fVals[17] = log(BC1.e2nd/BC1.eMax);
    fVals[18] = log(BC1.eTop/BC1.eMax);
    fVals[19] = log(BC1.eBottom/BC1.eMax);
    fVals[20] = log(BC1.eLeft/BC1.eMax);
    fVals[21] = log(BC1.eRight/BC1.eMax);
    fVals[22] = (BC1.eTop-BC1.eBottom)/(BC1.eTop + BC1.eBottom);
    fVals[23] = (BC1.eLeft-BC1.eRight)/(BC1.eLeft+BC1.eRight);

    if(ele.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = ele.SC.basicClusters[1];
      fVals[24]  = BC2.eta - ele.SC.eta;
      fVals[25]  = DeltaPhi(BC2.phi,ele.SC.phi);
      fVals[26]  = BC2.energy/ele.SC.rawE;
      fVals[27]  = BC2.e3x3/BC2.energy;
      fVals[28]  = BC2.e5x5/BC2.energy;
      fVals[29] = sqrt(BC2.sigmaIEtaIEta);
      fVals[30] = sqrt(BC2.sigmaIPhiIPhi);
      fVals[31] = BC2.sigmaIEtaIPhi;
      fVals[32] = BC2.eMax/BC2.energy;
      fVals[33] = log(BC2.e2nd/BC2.eMax);
      fVals[34] = log(BC2.eTop/BC2.eMax);
      fVals[35] = log(BC2.eBottom/BC2.eMax);
      fVals[36] = log(BC2.eLeft/BC2.eMax);
      fVals[37] = log(BC2.eRight/BC2.eMax);
      fVals[38] = (BC2.eTop-BC2.eBottom)/(BC2.eTop + BC2.eBottom);
      fVals[39] = (BC2.eLeft-BC2.eRight)/(BC2.eLeft+BC2.eRight);
    }else{ // only one basic cluster
      for(int i=24;i<40;i++) fVals[i] = 0.;
    }
    
    if(ele.SC.basicClusters.size() > 2){ // now look for the lowest energy basic cluster (for pileup mitigation)
      VecbosBC BCL = ele.SC.basicClusters.back(); // last element in the energy sorted list
      fVals[40] = BCL.eta-ele.SC.eta;
      fVals[41] = DeltaPhi(BCL.phi,ele.SC.phi);
      fVals[42] = BCL.energy/ele.SC.rawE;
      fVals[43] = BCL.e3x3/BCL.energy;
      fVals[44] = BCL.e5x5/BCL.energy;
      fVals[45] = sqrt(BCL.sigmaIEtaIEta);
      fVals[46] = sqrt(BCL.sigmaIPhiIPhi);
      fVals[47] = BCL.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=40;i<48;i++) fVals[i] = 0;
    }
    
    if(ele.SC.basicClusters.size() > 3){ // now look for the 2nd lowest energy basic cluster (for pileup mitigation)
      VecbosBC BC2L = ele.SC.basicClusters.at(ele.SC.basicClusters.size()-2); // last element in the energy sorted list
      fVals[48] = BC2L.eta-ele.SC.eta;
      fVals[49] = DeltaPhi(BC2L.phi,ele.SC.phi);
      fVals[50] = BC2L.energy/ele.SC.rawE;
      fVals[51] = BC2L.e3x3/BC2L.energy;
      fVals[52] = BC2L.e5x5/BC2L.energy;
      fVals[53] = sqrt(BC2L.sigmaIEtaIEta);
      fVals[54] = sqrt(BC2L.sigmaIPhiIPhi);
      fVals[55] = BC2L.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=48;i<56;i++) fVals[i] = 0;
    }
    
   
 //local coordinates and crystal indices
    //seed cluster

    fVals[56] = BC1.iEta; //crystal ieta
    fVals[57] = BC1.iPhi; //crystal iphi
    fVals[58] = BC1.iEta%5; //submodule boundary eta symmetry
    fVals[59] = BC1.iPhi%2; //submodule boundary phi symmetry    
    fVals[60] = (TMath::Abs(BC1.iEta)<=25)*(BC1.iEta%25) + (TMath::Abs(BC1.iEta)>25)*((BC1.iEta-25*TMath::Abs(BC1.iEta)/BC1.iEta)%20);  //module boundary eta approximate symmetry            
    fVals[61] = BC1.iPhi%20; //module boundary phi symmetry                                                                                                                    
    fVals[62] = BC1.etaCrystal; //local coordinates with respect to closest crystal center at nominal shower depth                                                                 
    fVals[63] = BC1.phiCrystal;

    //2nd cluster (meaningful gap corrections for converted photons)                                                                                                        
    if(ele.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = ele.SC.basicClusters[1];
      fVals[64] = BC2.iEta;
      fVals[65] = BC2.iPhi;
      fVals[66] = BC2.iEta%5;
      fVals[67] = BC2.iPhi%2;
      fVals[68] = (TMath::Abs(BC2.iEta)<=25)*(BC2.iEta%25) + (TMath::Abs(BC2.iEta)>25)*((BC2.iEta-25*TMath::Abs(BC2.iEta)/BC2.iEta)%20);
      fVals[69] = BC2.iPhi%20;
      fVals[70] = BC2.etaCrystal;
      fVals[71] = BC2.phiCrystal;
    }else{
      for(int i=64;i<72;i++) fVals[i] =0.;
    }
    //Nb. of vertex
    fVals[72] = base->nPV;
    
    
    
  }else{
    
    fVals[0]  = ele.SC.rawE;
    fVals[1]  = ele.SC.r9();
    fVals[2]  = ele.SC.eta;
    fVals[3]  = ele.SC.phi;
    fVals[4]  = ele.SC.e5x5/ele.SC.rawE;
    fVals[5]  = ele.SC.etaWidth;
    fVals[6]  = ele.SC.phiWidth;
    //Nb. of vertex
    fVals[7] = base->nPV; 
    
  }
  
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = ele.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = ele.SC.rawE + ele.esEnergy;
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  
  return std::pair<double,double>(ecor,ecorerr);
  
  
}




std::pair<double,double> HggEGEnergyCorrector::photonEnergyCorrector_CorrectedEnergyWithErrorv2(int j){
    
  
  VecbosPho pho(base,j);

  THIS_ECAL_GEO thisGeometry = getGapCoordinates(ecalGeometry, pho.SC.eta, pho.SC.phi);
  
 

  Bool_t isbarrel = fabs(pho.SC.eta) < 1.48; 
  
  if( isbarrel){
        
    fVals[0]  = pho.SC.rawE;
    fVals[1]  = pho.SC.r9();
    fVals[2]  = pho.SC.eta;
    fVals[3]  = pho.SC.phi;
    fVals[4]  = pho.SC.e5x5/pho.SC.rawE;
    fVals[5] = pho.HoverE;
    fVals[6] = pho.SC.etaWidth;
    fVals[7] = pho.SC.phiWidth;
    ///bc
    VecbosBC BC1 = pho.SC.basicClusters[0];
    
    fVals[8]  = BC1.eta - pho.SC.eta;
    fVals[9]  = DeltaPhi(BC1.phi,pho.SC.phi);
    fVals[10]  = BC1.energy/pho.SC.rawE;
    fVals[11]  = BC1.e3x3/BC1.energy;
    fVals[12]  = BC1.e5x5/BC1.energy;
    fVals[13] = sqrt(BC1.sigmaIEtaIEta);
    fVals[14] = sqrt(BC1.sigmaIPhiIPhi);
    fVals[15] = BC1.sigmaIEtaIPhi;
    fVals[16] = BC1.eMax/BC1.energy;
    fVals[17] = log(BC1.e2nd/BC1.eMax);
    fVals[18] = log(BC1.eTop/BC1.eMax);
    fVals[19] = log(BC1.eBottom/BC1.eMax);
    fVals[20] = log(BC1.eLeft/BC1.eMax);
    fVals[21] = log(BC1.eRight/BC1.eMax);
    fVals[22] = (BC1.eTop-BC1.eBottom)/(BC1.eTop + BC1.eBottom);
    fVals[23] = (BC1.eLeft-BC1.eRight)/(BC1.eLeft+BC1.eRight);

    if(pho.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = pho.SC.basicClusters[1];
      fVals[24]  = BC2.eta - pho.SC.eta;
      fVals[25]  = DeltaPhi(BC2.phi,pho.SC.phi);
      fVals[26]  = BC2.energy/pho.SC.rawE;
      fVals[27]  = BC2.e3x3/BC2.energy;
      fVals[28]  = BC2.e5x5/BC2.energy;
      fVals[29] = sqrt(BC2.sigmaIEtaIEta);
      fVals[30] = sqrt(BC2.sigmaIPhiIPhi);
      fVals[31] = BC2.sigmaIEtaIPhi;
      fVals[32] = BC2.eMax/BC2.energy;
      fVals[33] = log(BC2.e2nd/BC2.eMax);
      fVals[34] = log(BC2.eTop/BC2.eMax);
      fVals[35] = log(BC2.eBottom/BC2.eMax);
      fVals[36] = log(BC2.eLeft/BC2.eMax);
      fVals[37] = log(BC2.eRight/BC2.eMax);
      fVals[38] = (BC2.eTop-BC2.eBottom)/(BC2.eTop + BC2.eBottom);
      fVals[39] = (BC2.eLeft-BC2.eRight)/(BC2.eLeft+BC2.eRight);
    }else{ // only one basic cluster
      for(int i=24;i<40;i++) fVals[i] = 0.;
    }
    
    if(pho.SC.basicClusters.size() > 2){ // now look for the lowest energy basic cluster (for pileup mitigation)
      VecbosBC BCL = pho.SC.basicClusters.back(); // last element in the energy sorted list
      fVals[40] = BCL.eta-pho.SC.eta;
      fVals[41] = DeltaPhi(BCL.phi,pho.SC.phi);
      fVals[42] = BCL.energy/pho.SC.rawE;
      fVals[43] = BCL.e3x3/BCL.energy;
      fVals[44] = BCL.e5x5/BCL.energy;
      fVals[45] = sqrt(BCL.sigmaIEtaIEta);
      fVals[46] = sqrt(BCL.sigmaIPhiIPhi);
      fVals[47] = BCL.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=40;i<48;i++) fVals[i] = 0;
    }
    
    if(pho.SC.basicClusters.size() > 3){ // now look for the 2nd lowest energy basic cluster (for pileup mitigation)
      VecbosBC BC2L = pho.SC.basicClusters.at(pho.SC.basicClusters.size()-2); // last element in the energy sorted list
      fVals[48] = BC2L.eta-pho.SC.eta;
      fVals[49] = DeltaPhi(BC2L.phi,pho.SC.phi);
      fVals[50] = BC2L.energy/pho.SC.rawE;
      fVals[51] = BC2L.e3x3/BC2L.energy;
      fVals[52] = BC2L.e5x5/BC2L.energy;
      fVals[53] = sqrt(BC2L.sigmaIEtaIEta);
      fVals[54] = sqrt(BC2L.sigmaIPhiIPhi);
      fVals[55] = BC2L.sigmaIEtaIPhi;
    }else{ //only 2 or fewer clusters
      for(int i=48;i<56;i++) fVals[i] = 0;
    }
    
    fVals[56] = BC1.iEta; //crystal ieta
    fVals[57] = BC1.iPhi; //crystal iphi
    fVals[58] = BC1.iEta%5; //submodule boundary eta symmetry
    fVals[59] = BC1.iPhi%2; //submodule boundary phi symmetry    
    fVals[60] = (TMath::Abs(BC1.iEta)<=25)*(BC1.iEta%25) + (TMath::Abs(BC1.iEta)>25)*((BC1.iEta-25*TMath::Abs(BC1.iEta)/BC1.iEta)%20);  //module boundary eta approximate symmetry            
    fVals[61] = BC1.iPhi%20; //module boundary phi symmetry                                                                                                                    
    fVals[62] = BC1.etaCrystal; //local coordinates with respect to closest crystal center at nominal shower depth                                                                 
    fVals[63] = BC1.phiCrystal;

    //2nd cluster (meaningful gap corrections for converted photons)                                                                                                        
    if(pho.SC.basicClusters.size() > 1){ // second highest energy cluster
      VecbosBC BC2 = pho.SC.basicClusters[1];
      fVals[64] = BC2.iEta;
      fVals[65] = BC2.iPhi;
      fVals[66] = BC2.iEta%5;
      fVals[67] = BC2.iPhi%2;
      fVals[68] = (TMath::Abs(BC2.iEta)<=25)*(BC2.iEta%25) + (TMath::Abs(BC2.iEta)>25)*((BC2.iEta-25*TMath::Abs(BC2.iEta)/BC2.iEta)%20);
      fVals[69] = BC2.iPhi%20;
      fVals[70] = BC2.etaCrystal;
      fVals[71] = BC2.phiCrystal;
    }else{
      for(int i=64;i<72;i++) fVals[i] =0.;
    }
    //Nb. of vertex
    fVals[72] = base->nPV;
    
    
    
  }else{
    
    fVals[0]  = pho.SC.rawE;
    fVals[1]  = pho.SC.r9();
    fVals[2]  = pho.SC.eta;
    fVals[3]  = pho.SC.phi;
    fVals[4]  = pho.SC.e5x5/pho.SC.rawE;
    fVals[5]  = pho.SC.etaWidth;
    fVals[6]  = pho.SC.phiWidth;
    //Nb. of vertex
    fVals[7] = base->nPV; 
    
  }

  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *greader;
  const GBRForest *greadervar;

  if (isbarrel) {
    den = pho.SC.rawE;
    greader = fReadereb;
    greadervar = fReaderebvariance;
  }
  else {
    den = pho.SC.rawE + pho.SC.esEnergy; 
    greader = fReaderee;
    greadervar = fReadereevariance;
  }
  
  Double_t ecor = greader->GetResponse(fVals)*den;
  Double_t ecorerr = greadervar->GetResponse(fVals)*den*varscale;
  
  
  return std::pair<double,double>(ecor,ecorerr);
  
  
}
