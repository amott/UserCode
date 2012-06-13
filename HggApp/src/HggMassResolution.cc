#include <HggMassResolution.hh>
#include <TLorentzVector.h>
#include "ReadConfig.hh"
#include "TMath.h"

#include <iostream>
using namespace std;
#define debugMassRes 0
HggMassResolution::HggMassResolution(){
  this->clear();

  Categories.push_back("EBlowEtaGoldGAP"); highR9.push_back(true);
  Categories.push_back("EBlowEtaGoldCM"); highR9.push_back(true);
  Categories.push_back("EBlowEtaBad"); highR9.push_back(false);
  Categories.push_back("EBhighEtaGold"); highR9.push_back(true);
  Categories.push_back("EBhighEtaBad"); highR9.push_back(false);
  Categories.push_back("EElowEtaGold"); highR9.push_back(true);
  Categories.push_back("EElowEtaBad"); highR9.push_back(false);
  Categories.push_back("EEhighEtaGold"); highR9.push_back(true);
  Categories.push_back("EEhighEtaBad"); highR9.push_back(false);

  minEta.push_back(0.); maxEta.push_back(1.);
  minEta.push_back(0.); maxEta.push_back(1.);
  minEta.push_back(0.); maxEta.push_back(1.);
  minEta.push_back(1.); maxEta.push_back(1.5);
  minEta.push_back(1.); maxEta.push_back(1.5);
  minEta.push_back(1.5); maxEta.push_back(2.);
  minEta.push_back(1.5); maxEta.push_back(2.);
  minEta.push_back(2.); maxEta.push_back(3.);
  minEta.push_back(2.); maxEta.push_back(3.);

  dzRes.push_back(0.1); dzRes.push_back(TMath::Sqrt(2.)*5.8);

}

void HggMassResolution::clear(){
  
}

double HggMassResolution::getMassResolution(VecbosPho *leadPho,VecbosPho *subleadPho, TVector3 vtx,bool isWrongVtx){
  TLorentzVector p4Pho1 = leadPho->p4FromVtx(vtx,leadPho->finalEnergy);
  TLorentzVector p4Pho2 = subleadPho->p4FromVtx(vtx,leadPho->finalEnergy);
  double angle   = p4Pho1.Angle(p4Pho2.Vect());
  if(debugMassRes) cout << ">> >> doing getResolution" << endl;
  double resPho1 = this->getResolution(leadPho);
  double resPho2 = this->getResolution(subleadPho);
  if(debugMassRes) cout << ">> >> doing getAngleResolution" << endl;
  double angleRes = this->getAngleResolution(leadPho,subleadPho,vtx,isWrongVtx);
  double higgsMass = (p4Pho1+p4Pho2).M();
  
  if(isWrongVtx){
    angleRes*=0.5*higgsMass;
    double massResEOnly = 0.5*higgsMass*TMath::Sqrt( (resPho1*resPho1)/(leadPho->finalEnergy*leadPho->finalEnergy) + (resPho2*resPho2)/(subleadPho->finalEnergy*subleadPho->finalEnergy) );
    return TMath::Sqrt((massResEOnly*massResEOnly)+(angleRes*angleRes));
  }
  
  return 0.5*higgsMass*TMath::Sqrt( (resPho1*resPho1)/(p4Pho1.E()*p4Pho1.E()) + (resPho2*resPho2)/(p4Pho2.E()*p4Pho2.E()) +
				    + ((angleRes*angleRes)*(TMath::Sin(angle)/(1.-TMath::Cos(angle))) * (TMath::Sin(angle)/(1.-TMath::Cos(angle)))) );
}

double HggMassResolution::getMassResolutionEonly(VecbosPho *leadPho,VecbosPho *subleadPho,TVector3 vtx){
  TLorentzVector p4Pho1 = leadPho->p4FromVtx(vtx,leadPho->finalEnergy);
  TLorentzVector p4Pho2 = subleadPho->p4FromVtx(vtx,leadPho->finalEnergy);
  double resPho1 = this->getResolution(leadPho);
  double resPho2 = this->getResolution(subleadPho);
  double higgsMass = (p4Pho1+p4Pho2).M();
  
  return 0.5*higgsMass*TMath::Sqrt( (resPho1*resPho1)/(leadPho->finalEnergy*leadPho->finalEnergy) + (resPho2*resPho2)/(subleadPho->finalEnergy*subleadPho->finalEnergy) );
}

void HggMassResolution::init(){
  ReadConfig cfg;
  cfg.read(config);

  char valString[400];
  for(int i=0;i<nCategories;i++){
    std::vector<string> thisSmear =  cfg.getTokens(Categories[i],",");
    if(thisSmear.size()!=2){
      cout << "ERROR: INVALID smearing configuration" << endl;
      cout << i << "    " << Categories[i] << "   " << thisSmear[0] << endl;
      throw -100;
      return;
    }
    smear[i] = pair<float,float>(atof(thisSmear[0].c_str()),atof(thisSmear[1].c_str()));
  }
}

double HggMassResolution::getAngleResolution(VecbosPho* pho1,VecbosPho* pho2, TVector3 vtx, bool wrongVtx){
  TVector3 pho1Pos= pho1->CaloPos - vtx;
  TVector3 pho2Pos = pho2->CaloPos - vtx;

  double r1 = pho1Pos.Mag();
  double r2 = pho2Pos.Mag();
  double cos = TMath::Cos(pho1Pos.Phi()-pho2Pos.Phi());
  double sech1 = 1./TMath::CosH(pho1Pos.Eta());
  double sech2 = 1./TMath::CosH(pho2Pos.Eta());
  double tanh1 = TMath::TanH(pho1Pos.Eta());
  double tanh2 = TMath::TanH(pho2Pos.Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos);
  double denominator = 1. - tanh1*tanh2 - sech1*sech2*cos;
  if(debugMassRes) cout << ">> >> got Angle Resolution: " << (-1.*dzRes[wrongVtx]/denominator)*(numerator1/r1 + numerator2/r2) << endl;
  return (-1.*dzRes[wrongVtx]/denominator)*(numerator1/r1 + numerator2/r2);
}

float HggMassResolution::getResolution(VecbosPho* pho){
  pair<float,float> catRes = smear[this->getCategory(pho)];
  if(debugMassRes) cout << "got Res" << endl;
  return TMath::Sqrt(pho->correctedEnergyError*pho->correctedEnergyError+catRes.second*catRes.second);
}

int HggMassResolution::getCategory(VecbosPho* pho){
  if(debugMassRes) cout << "getCategory" << endl;
  if(this->isSphericalPhoton(pho->SC.BCSeed.iEta,pho->SC.BCSeed.iPhi)){ // these are the "special" photons that have to be treated differently
    if(pho->SC.r9() > r9Cut) return sphericalIndex;
  }
  if(debugMassRes) cout << "not a special photon" << endl;

  for(int i=0;i<nCategories;i++){
    if(i == sphericalIndex) continue;

    if( minEta[i] <= fabs(pho->SC.eta)
	&& fabs(pho->SC.eta) < maxEta[i]
	&& (pho->SC.r9() > r9Cut) == highR9[i]) return i;
  }
  cout << "Photon with no category??? r9: " << pho->SC.r9() << "  eta: " << pho->SC.eta << "  phi: " << pho->SC.phi << endl;
  throw -99;
  return -1;
}

bool HggMassResolution::isSphericalPhoton(int ieta, int iphi)
{                                                             
  if(debugMassRes) cout << "isSphericalPhoton" << endl;
  if ((iphi %20)<=5 || (iphi%20)>=16){                      
    return false;                                         
  }                                                         
  
  int ietaTT=(std::abs(ieta)-1)/5+1;                        
    if                                                        
      (                                                     
       (ietaTT>= 2&&     ietaTT<    5 ) ||               
       (ietaTT>= 7&&     ietaTT<    9 ) ||               
       (ietaTT>= 11&&    ietaTT<    13) ||               
       (ietaTT>= 15&&    ietaTT<    17)                  
       ){                                                
      return true;                                          
    }                                                         
                                                              
                                                              
    return false;                                             
}                                                             
