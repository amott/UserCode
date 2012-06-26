#include "TTreeFormula.h"
#include "TString.h"
#include "TChain.h"

typedef std::map<TString,TTreeFormula*> CatMap;

void updateFormulas(CatMap* m){
  CatMap::iterator it;
  for(it=m->begin(); it!=m->end(); it++) (*it).second->UpdateFormulaLeaves();
}

int getCat(CatMap* cats, TString type, const int nCat){
  if( (*cats)[ Form("%s_R",type.Data()) ]->EvalInstance() ){
    return -1; // rejections
  }
  if( (*cats)[ Form("%s_PreSelEta",type.Data()) ]->EvalInstance() ){
    return -1; // rejection photons between EB and EE
  }
  if( (*cats)[ Form("%s_PreSelEt",type.Data()) ]->EvalInstance() ){
    return -1; // final pT requirements on the photons
  }
  

  for(int iCat=0;iCat<nCat;iCat++){
    if( (*cats)[ Form("%s_%d",type.Data(),iCat) ]->EvalInstance() ) return iCat;
  }
  return -1;
}

CatMap getCategoryCuts(TChain* fChain){
  std::map<TString,TTreeFormula*> categories;


  //categories["multiplicity"] = new TTreeFormula("Mult","Photon>=2",fChain);
  categories["MVA_PreSelEta"] = new TTreeFormula("MVA_PreSelEta","(abs(Photon[0].eta) > 1.4442 && abs(Photon[0].eta) < 1.566) || (abs(Photon[1].eta) > 1.4442 && abs(Photon[1].eta) < 1.566) || abs(Photon[0].eta) > 2.5 || abs(Photon[1].eta) > 2.5",fChain);
  categories["MVA_PreSelEt"] = new TTreeFormula("MVA_PreSelEt","Photon[0].pt < 30 || Photon[1].pt < 30 || (Photon[0].pt < 40 && Photon[1].pt < 40)",fChain);
  categories["MVA_R"] = new TTreeFormula("MVA_0","mPair<0 || diPhotonMVA<0.05 || nPhoton==0",fChain); // rejection criterion
  categories["MVA_0"] = new TTreeFormula("MVA_0"," (Mjj >= 500) && (diPhotonMVA>0.05)",fChain);  
  categories["MVA_1"] = new TTreeFormula("MVA_1"," (Mjj < 500) && (Mjj >= 250) && (diPhotonMVA>0.05)",fChain);  
  categories["MVA_2"] = new TTreeFormula("MVA_2"," (Mjj < 250) && (diPhotonMVA>=0.88)",fChain);  
  categories["MVA_3"] = new TTreeFormula("MVA_3"," (Mjj < 250) && (diPhotonMVA>=0.71) && (diPhotonMVA<0.88)",fChain);  
  categories["MVA_4"] = new TTreeFormula("MVA_4"," (Mjj < 250) && (diPhotonMVA>=0.50) && (diPhotonMVA<0.71)",fChain);  
  categories["MVA_5"] = new TTreeFormula("MVA_5"," (Mjj < 250) && (diPhotonMVA>=-0.05) && (diPhotonMVA<0.50)",fChain);  
  
  categories["PFCiC_PreSelEta"] = new TTreeFormula("PFCiC_PreSelEta","(abs(PhotonPFCiC[0].eta) > 1.4442 && abs(PhotonPFCiC[0].eta) < 1.566) || (abs(PhotonPFCiC[1].eta) > 1.4442 && abs(PhotonPFCiC[1].eta) < 1.566) || abs(PhotonPFCiC[0].eta) > 2.5 || abs(PhotonPFCiC[1].eta) > 2.5",fChain);
  categories["PFCiC_PreSelEt"] = new TTreeFormula("PFCiC_PreSelEt","PhotonPFCiC[0].pt < 30 || PhotonPFCiC[1].pt < 30 || (PhotonPFCiC[0].pt < 40 && PhotonPFCiC[1].pt < 40)",fChain);
  categories["PFCiC_R"] = new TTreeFormula("PFCiC_R","mPairPFCiC==-1 || nPhotonPFCiC==0",fChain); // rejection criterion
  categories["PFCiC_0"] = new TTreeFormula("PFCiC_0","(MjjPFCiC >= 500)",fChain);
  categories["PFCiC_1"] = new TTreeFormula("PFCiC_1","(MjjPFCiC < 500) && (MjjPFCiC >= 250)",fChain);
  categories["PFCiC_2"] = new TTreeFormula("PFCiC_2","(MjjPFCiC < 250) && (PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && (abs(PhotonPFCiC[0].eta) < 1.48 && abs(PhotonPFCiC[1].eta) < 1.48)",fChain);
  categories["PFCiC_3"] = new TTreeFormula("PFCiC_3","(MjjPFCiC < 250) && !(PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && (abs(PhotonPFCiC[0].eta) < 1.48 && abs(PhotonPFCiC[1].eta) < 1.48)",fChain);
  categories["PFCiC_4"] = new TTreeFormula("PFCiC_4","(MjjPFCiC < 250) && (PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && !(abs(PhotonPFCiC[0].eta) < 1.48 && abs(PhotonPFCiC[1].eta) < 1.48)",fChain);
  categories["PFCiC_5"] = new TTreeFormula("PFCiC_5","(MjjPFCiC < 250) && !(PhotonPFCiC[0].r9 > 0.94 && PhotonPFCiC[1].r9 > 0.94) && !(abs(PhotonPFCiC[0].eta) < 1.48 && abs(PhotonPFCiC[1].eta) < 1.48)",fChain);

  categories["CiC_PreSelEta"] = new TTreeFormula("CiC_PreSelEta","(abs(PhotonCiC[0].eta) > 1.4442 && abs(PhotonCiC[0].eta) < 1.566) || (abs(PhotonCiC[1].eta) > 1.4442 && abs(PhotonCiC[1].eta) < 1.566) || abs(PhotonCiC[0].eta) > 2.5 || abs(PhotonCiC[1].eta) > 2.5",fChain);
  categories["CiC_PreSelEt"] = new TTreeFormula("CiC_PreSelEt","PhotonCiC[0].pt < 30 || PhotonCiC[1].pt < 30 || (PhotonCiC[0].pt < 40 && PhotonCiC[1].pt < 40)",fChain);
  categories["CiC_R"] = new TTreeFormula("CiC_R","mPairCiC==-1 || nPhotonCiC==0",fChain); // rejection criterion
  categories["CiC_0"] = new TTreeFormula("CiC_0","(MjjCiC >= 500)",fChain);
  categories["CiC_1"] = new TTreeFormula("CiC_1","(MjjCiC < 500) && (MjjCiC >= 250)",fChain);
  categories["CiC_2"] = new TTreeFormula("CiC_2","(MjjCiC < 250) && (PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && (abs(PhotonCiC[0].eta) < 1.48 && abs(PhotonCiC[1].eta) < 1.48)",fChain);
  categories["CiC_3"] = new TTreeFormula("CiC_3","(MjjCiC < 250) && !(PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && (abs(PhotonCiC[0].eta) < 1.48 && abs(PhotonCiC[1].eta) < 1.48)",fChain);
  categories["CiC_4"] = new TTreeFormula("CiC_4","(MjjCiC < 250) && (PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && !(abs(PhotonCiC[0].eta) < 1.48 && abs(PhotonCiC[1].eta) < 1.48)",fChain);
  categories["CiC_5"] = new TTreeFormula("CiC_5","(MjjCiC < 250) && !(PhotonCiC[0].r9 > 0.94 && PhotonCiC[1].r9 > 0.94) && !(abs(PhotonCiC[0].eta) < 1.48 && abs(PhotonCiC[1].eta) < 1.48)",fChain);  

  CatMap::iterator it;
  for(it=categories.begin(); it!=categories.end(); it++) (*it).second->SetQuickLoad(true);  


  return categories;
}
