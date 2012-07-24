#include "include/RedVecbosMuTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedVecbosMuTree::RedVecbosMuTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myFile->cd();
  myTree = new TTree("VecBosMu","vecbos analysis");

  myTree->Branch("goodNJets",	&myNJets,	"goodNJets/I");  
  myTree->Branch("allNJets",	&allNJets,	"allNJets/I");  

  myTree->Branch("weight",	&weight,	"weight/F");  

  myTree->Branch("energyJet",	 energyJet,    "etenergyJet[20]/F");
  myTree->Branch("etJet", 	 etJet,        "etJet[20]/F"); 	     
  myTree->Branch("etaJet", 	 etaJet,       "etaJet[20]/F");       
  myTree->Branch("phiJet", 	 phiJet,       "phiJet[20]/F");       
  myTree->Branch("vertexXJet",	 vertexXJet,   "vertexXJet[20]/F");   
  myTree->Branch("vertexYJet",	 vertexYJet,   "vertexYJet[20]/F");   
  myTree->Branch("vertexZJet",	 vertexZJet,   "vertexZJet[20]/F");   
  myTree->Branch("emFracJet",	 emFracJet,    "emFracJet[20]/F");    
  myTree->Branch("hadFracJet",	 hadFracJet,   "hadFracJet[20]/F");   
					        	     
  myTree->Branch("nMuons",	&nMuons,	"nMuons/I");  
  myTree->Branch("nGoodMuons",	&nGoodMuons,	"nGoodMuons/I");  
  myTree->Branch("ZMassForVeto", &ZMassForVeto, "ZMassForVeto/F");  
  myTree->Branch("nZmumuCand",	&nZmumuCand,	"nZmumuCand/I");  
  myTree->Branch("bestMuEta",     bestMuEta,		"bestMuEta[2]/F");  
  myTree->Branch("bestMuPhi",     bestMuPhi,		"bestMuPhi[2]/F");  
  myTree->Branch("bestMuPt",      bestMuPt,		"bestMuPt[2]/F");  
  myTree->Branch("bestMuvtxX",    bestMuvtxX,	  "bestMuvtxX[2]/F");
  myTree->Branch("bestMuvtxY",    bestMuvtxY,	  "bestMuvtxY[2]/F");
  myTree->Branch("bestMuvtxZ",    bestMuvtxZ,	  "bestMuvtxZ[2]/F");
  myTree->Branch("bestMucharge",     bestMucharge,   "bestMucharge[2]/F"); 
  myTree->Branch("bestMuPx",	      bestMuPx,   "bestMuPx[2]/F");  
  myTree->Branch("bestMuPy",	      bestMuPy,   "bestMuPy[2]/F"); 
  myTree->Branch("muTrackDxy",       muTrackDxy,   "muTrackDxy[2]/F");     
  myTree->Branch("muTrackD0",        muTrackD0,   "muTrackD0[2]/F");	   
  myTree->Branch("muTrackDsz",       muTrackDsz,   "muTrackDsz[2]/F");     
  myTree->Branch("muTrackDz",        muTrackDz,   "muTrackDz[2]/F");	   
  myTree->Branch("muTrackDxyError",  muTrackDxyError,   "muTrackDxyError[2]/F");
  myTree->Branch("muTrackD0Error",   muTrackD0Error,   " muTrackD0Error[2]/F"); 
  myTree->Branch("muTrackDszError",  muTrackDszError,   "muTrackDszError[2]/F");
  myTree->Branch("muTrackDzError",   muTrackDzError,   "muTrackDzError[2]/F"); 
  myTree->Branch("sumPt03",	      sumPt03,   "sumPt03[2]/F");	   
  myTree->Branch("emEt03",	      emEt03,   "emEt03[2]/F");	   
  myTree->Branch("hadEt03",	      hadEt03,   "hadEt03[2]/F");	   
  myTree->Branch("hoEt03",	      hoEt03,   "hoEt03[2]/F");	   
  myTree->Branch("nTrk03",	      nTrk03,   "nTrk03[2]/F");	   
  myTree->Branch("nJets03",	      nJets03,   "nJets03[2]/F");	   
  myTree->Branch("sumPt05",	      sumPt05,   "sumPt05[2]/F");	   
  myTree->Branch("emEt05",	      emEt05,   "emEt05[2]/F");	   
  myTree->Branch("hadEt05",	      hadEt05,   "hadEt05[2]/F");	   
  myTree->Branch("hoEt05",	      hoEt05,   "hoEt05[2]/F");	   
  myTree->Branch("nTrk05",	      nTrk05,   "nTrk05[2]/F");	   
  myTree->Branch("nJets05",	      nJets05,   "nJets05[2]/F");	   
  myTree->Branch("EcalExpDepo",      EcalExpDepo,   "EcalExpDepo[2]/F");    
  myTree->Branch("HcalExpDepo",      HcalExpDepo,   "HcalExpDepo[2]/F");    
  myTree->Branch("HoExpDepo",        HoExpDepo,   "HoExpDepo[2]/F");	   
  myTree->Branch("emS9",	      emS9,   "emS9[2]/F");	   
  myTree->Branch("hadS9",	      hadS9,   "hadS9[2]/F");	   
  myTree->Branch("hoS9",	      hoS9,   "hoS9[2]/F");	   
  myTree->Branch("CaloComp",	      CaloComp,   "CaloComp[2]/F");	   
  myTree->Branch("SumMuonsPx",	      &SumMuonsPx,   "SumMuonPx/F");	   
  myTree->Branch("SumMuonsPy",	      &SumMuonsPy,   "SumMuonPy/F");	   
  
  myTree->Branch("Z0Mass",                 &Z0Mass,                 "Z0Mass/F");  
  myTree->Branch("Z0Pt",                 &Z0Pt,                 "Z0Pt/F");  
  myTree->Branch("Z0Eta",                 &Z0Eta,                 "Z0Eta/F");  
  myTree->Branch("Z0Phi",                 &Z0Phi,                 "Z0Phi/F");  
  
  

  myTree->Branch("met",                 &myMet,                 "met/F");  
  myTree->Branch("metx",                &myMetx,                "metx/F");  
  myTree->Branch("mety",                &myMety,                "mety/F");  

  myTree->Branch("transvMass",          &myTransvMass,          "transvMass/F");  

  myTree->Branch("helicity",            &myHelicity,            "helicity/F");  

  myTree->Branch("nPVTX",       &nPVTX,      "nPVTX/I");
  myTree->Branch("PVTXxPV",     PVTXxPV,     "PVTXxPV[10]/F");
  myTree->Branch("PVTXyPV",     PVTXyPV,     "PVTXyPV[10]/F");
  myTree->Branch("PVTXzPV",     PVTXzPV,     "PVTXzPV[10]/F");
  myTree->Branch("PVTXErrxPV",  PVTXErrxPV,  "PVTXErrxPV[10]/F");
  myTree->Branch("PVTXErryPV",  PVTXErryPV,  "PVTXErryPV[10]/F");
  myTree->Branch("PVTXErrzPV",  PVTXErrzPV,  "PVTXErrzPV[10]/F");
  myTree->Branch("SumPtPVTX",   SumPtPVTX,   "SumPtPVTX[10]/F"); 
  myTree->Branch("ndofPVTX",    ndofPVTX,    "ndofPVTX[10]/F");  
  myTree->Branch("chi2PVTX",    chi2PVTX,    "chi2PVTX[10]/F");  

}

RedVecbosMuTree::~RedVecbosMuTree() {delete myFile;}

void RedVecbosMuTree::addCSA07Infos() {
  myTree->Branch("CSA07weight",     &myWeight,     "CSA07weight/D");
  myTree->Branch("CSA07processId",  &myProcesId,   "CSA07processId/D");
  myTree->Branch("CSA07lumi",       &myLumi,       "CSA07lumi/F");
}

void RedVecbosMuTree::store() {myTree->Fill();}

//void RedVecbosMuTree::saveEfficiencies(std::vector<float>);


void RedVecbosMuTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedVecbosMuTree::fillAll(int nj, float eta1, float phi1, float pt1, float eta2, float phi2, float pt2, float helic, float imass, float tmass, float met) {
 }

void RedVecbosMuTree::fillAll(int nj, float eta1, float phi1, float pt1, float eta2, float phi2, float pt2, float helic, float imass, float tmass, float met, float metx, float mety) {
}

void RedVecbosMuTree::fillCSA07(double weight, double processId, float lumi) {
  myWeight = weight;
  myProcesId = processId;
  myLumi = lumi;
}

void RedVecbosMuTree::newEvent(){
  myNJets =-999;
  allNJets = -999;
  nMuons=-999;
  nGoodMuons=-999;
  nZmumuCand=-999;
  
  for(int index=0;index<2;++index){
    bestMuEta[index]=-999 ;  
    bestMuPhi[index]=-999 ;  
    bestMuPt[index]=-999 ;  
    bestMuvtxX[index]=-999 ;
    bestMuvtxY[index]=-999 ;
    bestMuvtxZ[index]=-999 ;
    bestMucharge[index]=-999.; 
    bestMuPx[index]=-999.;   
    bestMuPy[index]=-999.;  
    muTrackDxy[index]=-999.;         
    muTrackD0[index]=-999.;	       
    muTrackDsz[index]=-999.;         
    muTrackDz[index]=-999.;	       
    muTrackDxyError[index]=-999.;  
    muTrackD0Error[index]=-999.;   
    muTrackDszError[index]=-999.;  
    muTrackDzError[index]=-999.;   
    sumPt03[index]=-999.;    
    emEt03[index]=-999.;     
    hadEt03[index]=-999.;    
    hoEt03[index]=-999.;     
    nTrk03[index]=-999.;     
    nJets03[index]=-999.;    
    sumPt05[index]=-999.;    
    emEt05[index]=-999.;     
    hadEt05[index]=-999.;    
    hoEt05[index]=-999.;     
    nTrk05[index]=-999.;     
    nJets05[index]=-999.;    
    EcalExpDepo[index]=-999.;        
    HcalExpDepo[index]=-999.;        
    HoExpDepo[index]=-999.;	       
    emS9[index]=-999.;       
    hadS9[index]=-999.;      
    hoS9[index]=-999.;       
    CaloComp[index]=-999.;	       
  }     
SumMuonsPx =-999.;
SumMuonsPy =-999.; 
  
  for(int jetInd=0; jetInd<20 ; ++jetInd){
    energyJet[jetInd]=-999;
    etJet[jetInd]=-999; 
    etaJet[jetInd]=-999; 
    phiJet[jetInd]=-999; 
    vertexXJet[jetInd]=-999;
    vertexYJet[jetInd]=-999;
    vertexZJet[jetInd]=-999;
    emFracJet[jetInd]=-999;
    hadFracJet[jetInd]=-999;
  }

  ZMassForVeto=-999; 
 
  Z0Mass=-999;
  Z0Pt  =-999;
  Z0Eta =-999;
  Z0Phi =-999;

  weight=-999.;

  myHelicity =-999 ;
  myTransvMass=-999 ;
  myMet      =-999 ;
  myMetx     =-999 ;
  myMety     =-999 ;
  myWeight  =-999.;
  myProcesId=-999.;
  myLumi    =-999 ;
  nPVTX=-999 ;
  for(int vtxIndex=0; vtxIndex<10;++vtxIndex){
	PVTXxPV[vtxIndex]=-999 ;   
	PVTXyPV[vtxIndex]=-999 ;   
	PVTXzPV[vtxIndex]=-999 ;   
	PVTXErrxPV[vtxIndex]=-999 ;
	PVTXErryPV[vtxIndex]=-999 ;
	PVTXErrzPV[vtxIndex]=-999 ;
	SumPtPVTX[vtxIndex]=-999 ; 
	ndofPVTX[vtxIndex]=-999 ;
	chi2PVTX[vtxIndex]=-999 ;
  }


}

