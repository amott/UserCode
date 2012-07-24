// std includes
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>

using namespace std;

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

// VecbosApp includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "VecbosAnalysis.hh"
#define DUMMY 0

VecbosAnalysis::VecbosAnalysis(TTree *tree, string namefile) : Vecbos(tree) {
  // Set the analysis parameters to default values
  InitParameters(namefile);
}

VecbosAnalysis::~VecbosAnalysis(){}

void VecbosAnalysis::Loop() {
  
  // emanuele: FIXME
  float genWeight = 1.0;
  float genProcessId = 1;
  float genAlpgenID = 1;

  if(fChain == 0) return;

  // Here I put Ilaria's beginJob()
  this->beginJob();

  // Seems that Ilaria had a different way of booking histos.
  // Check some lines down.
  this->bookProcessIDHistoes();

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  
  // Here we begin the loop over events itself.
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    
    processIDInt = 0;
    Float_t   processIDIntOriginal;

    if(!useGlobalWeight){
      
      EventWeight = genWeight;
      EventWeightFloat = (float) EventWeight; 
      
      int NotAlpgen=soupType_.compare("ALPGEN"); //Returns false if equal
      
      if(NotAlpgen){
	///For PYTHIA SOUPS:
	float processID = genProcessId;
	processIDInt = (int) processID;
	processIDIntOriginal=processIDInt;
      } else{
	///FOR ALPGEN SOUPS:
	processIDIntOriginal = (int) genAlpgenID;    
	if ( processIDIntOriginal<1500) 			      processIDInt =0; //W+Jets
	if ((processIDIntOriginal>1999) & (processIDIntOriginal<2020))  processIDInt =1; //Z+Jets
	if ( processIDIntOriginal>2500) 			      processIDInt =2; //ttbar+Jets
      }
      
      double binx = (double) processIDIntOriginal;
      TAxis   * axis = Weight_vs_ProcessID->GetXaxis();
      Int_t binNumber = axis->FindBin(binx);
      Weight_vs_ProcessID->SetBinContent(binNumber,EventWeight);
    }

    // I substituted these lines for booking the histograms before the loop.
    //    if(AllEventsMET.find(processIDInt)==AllEventsMET.end())
    //      this->bookProcessIDHistoes();
    
    // Trigger still to be implemented.    
    //     iEvent.getByLabel(triggerResultsTag_,hltresults);
    //     if(hltresults.isValid()){
    //       trigNames.init((edm::TriggerResults&)*hltresults);
    //       numTriggers =   trigNames.size();
    //       unsigned int evento = (iEvent.id()).event() ;
    //       this->fillTriggerHist(TriggersBeforeSelection[processIDInt],TriggersBeforeSelection1D[processIDInt], evento);
    //     }else{
    //       std::cout << "NOTRIGGER" << std::endl;
    //     }
    
    this->PlotsBeforeSelection();

    // Ask Ilaria how to work with this.
    foundZCandidate[ZRecoFlagNumber]=false;
    if(checkOtherRecoTypes){
      for(unsigned int previousCandidateTypes = 0; previousCandidateTypes<ZRecoFlagNumber;
     	  ++previousCandidateTypes){
     	if(foundZCandidate[previousCandidateTypes]) continue;
      }
    }
    
    if(ZSelection){
      if(!this->doZed() ) continue;
      ++EventsAfterZedSelection;
    }
    
    if(METSelection){
      if (!this->doMET() ) continue;
      
      if(WSelection){
	if( !this->doW() ) continue;
      }
      
    }
    
    // More trigger stuff.
    //unsigned int evento = iEvent.id().event() ;
    //this->fillTriggerHist(TriggersAfterEWKSelection[processIDInt],TriggersAfterEWKSelection1D[processIDInt], evento );

    if(JetSelection){	
      if (!this->doRecoJets() ) continue;
      ProcessIDAfterJets->Fill(processIDIntOriginal, EventWeight); 
      ++EventsAfterJetSelection;
    }  
    
    // More trigger stuff.
    // this->fillTriggerHist(TriggersAfterJETSelection[processIDInt],TriggersAfterJETSelection1D[processIDInt], evento );
    
    this->fillHistoes();
    
    foundZCandidate[ZRecoFlagNumber]=true;
    TypeOfZCounter[ZRecoFlagNumber] += EventWeight;
    
  } // Closes the loop itself.

  // Here I put Ilaria's beginJob()
  this->endJob();

} // Closes VecbosAnalysis::Loop

/////////////////////////
// Book the histograms //
/////////////////////////

void VecbosAnalysis::bookProcessIDHistoes(){
  if(verbose) std::cout<<"In VecbosAnalysis::bookProcessIDHistoes"<<std::endl;
  
  fOutputFile ->cd();
  char histoid[50];
  
  
  sprintf(histoid,"TriggersBeforeSelection_Proc%d",processIDInt);	  
  TriggersBeforeSelection[processIDInt]  = new TH2D(histoid, "TriggersBeforeSelection",90,0.5,90.5,90,0.5,90.5);      
  sprintf(histoid,"TriggersBeforeSelection1D_Proc%d",processIDInt);	  
  TriggersBeforeSelection1D[processIDInt]  = new TH1D(histoid, "Triggers BeforeSelection",90,0.5,90.5);  
  
  sprintf(histoid,"TriggersAfterEWKSelection_Proc%d",processIDInt);	  
  TriggersAfterEWKSelection[processIDInt]  = new TH2D(histoid, "TriggersAfterEWKSelection",90,0.5,90.5,90,0.5,90.5);  
  sprintf(histoid,"TriggersAfterEWKSelection1D_Proc%d",processIDInt);	  
  TriggersAfterEWKSelection1D[processIDInt]  = new TH1D(histoid, "TriggersAfterEWKSelection",90,0.5,90.5);  
  
  sprintf(histoid,"TriggersAfterJETSelection_Proc%d",processIDInt);	  
  TriggersAfterJETSelection[processIDInt]  = new TH2D(histoid, "TriggersAfterJETSelection",90,0.5,90.5,90,0.5,90.5);  
  sprintf(histoid,"TriggersAfterJETSelection1D_Proc%d",processIDInt);	  
  TriggersAfterJETSelection1D[processIDInt]  = new TH1D(histoid, "TriggersAfterJETSelection",90,0.5,90.5);  
    
  sprintf(histoid,"TriggersInZMassRegion_Proc%d",processIDInt);	
  TriggersInZMassRegion[processIDInt]  = new TH2D(histoid, "TriggersInZMassRegion",90,0.5,90.5,90,0.5,90.5);  
  sprintf(histoid,"TriggersInZMassRegion1D_Proc%d",processIDInt);	
  TriggersInZMassRegion1D[processIDInt]  = new TH1D(histoid, "TriggersInZMassRegion", 90,0.5,90.5);  
  
  
  
  sprintf(histoid,"METNoCuts_Proc%d",processIDInt);	    
  AllEventsMET[processIDInt]  = new TH1D(histoid, "METNoCuts",200,0.,200.);
  
  sprintf(histoid,"AllMuonsIso1_Proc%d",processIDInt);	    
  AllMuonsIso1[processIDInt]  = new TH1D(histoid, "AllMuonsIso1",120,0.,100.);
  
  sprintf(histoid,"AllMuonsIso2_Proc%d",processIDInt);	    
  AllMuonsIso2[processIDInt]  = new TH1D(histoid, "AllMuonsIso2",120,0.,100.);
  
  sprintf(histoid,"METNoCutsvsIso1_Proc%d",processIDInt);	    
  AllEventsMET_vs_Iso1[processIDInt]  = new TH2D(histoid, "METNoCuts vs Iso1",200,0.,200.,120,0.,100.);
    
  sprintf(histoid,"METNoCutsvsIso2_Proc%d",processIDInt);	    
  AllEventsMET_vs_Iso2[processIDInt]  = new TH2D(histoid, "METNoCuts vs Iso2",200,0.,200.,120,0.,100.);
  
  
  sprintf(histoid,"WTransvMass_Proc%d",processIDInt);	    
  WTransvMass[processIDInt] = new TH1D(histoid, "WTransvMass",200,0.,200.);
  
  sprintf(histoid,"NumberOfWCandidates_Proc%d",processIDInt);	    
  NumberOfWCandidates[processIDInt] = new TH1D(histoid, "NumberOfWCandidates",10,0.5,10.5);
  
  
  sprintf(histoid,"TypeofZReco_Proc%d",processIDInt);	    
  TypeOfZCandidates[processIDInt] = new TH1D(histoid, "TypeofZReco",3,0.5,3.5);

  sprintf(histoid,"NumberOfZCandidates_Proc%d",processIDInt);	    
  NumberOfZCandidates[processIDInt] = new TH1D(histoid, "NumberOfZCandidates",10,0.5,10.5);
  
  sprintf(histoid,"ZCandidate_All_Proc%d",processIDInt);	    
  ZCandidatesMass[0][processIDInt]=  new TH1D(histoid, histoid,120,0.,180.);
   
  sprintf(histoid,"ZCandidate_BARREL-BARREL_Proc%d",processIDInt);	    
  ZCandidatesMass[1][processIDInt]=  new TH1D(histoid, histoid,120,0.,180.);
  
  sprintf(histoid,"ZCandidate_BARREL-ENDCAP_Proc%d",processIDInt);	    
  ZCandidatesMass[2][processIDInt]=  new TH1D(histoid, histoid,120,0.,180.);
   
  sprintf(histoid,"ZCandidate_ENDCAP-BARREL_Proc%d",processIDInt);	    
  ZCandidatesMass[3][processIDInt]=  new TH1D(histoid, histoid,120,0.,180.);
  
  sprintf(histoid,"ZCandidate_ENDCAP-ENDCAP_Proc%d",processIDInt);	    
  ZCandidatesMass[4][processIDInt]=  new TH1D(histoid, histoid,120,0.,180.);
  

  
  sprintf(histoid,"ZCandidatePt_Proc%d",processIDInt);	    
  ZCandidatesPt[processIDInt]  = new TH1D(histoid, "Z(mu mu) candidate Transverse Momentum",400,0.,400.);
  
  
  
  for(unsigned int muIndex=0; muIndex<3; ++muIndex){
    
    sprintf(histoid,"Muon_%d_FromZ_P_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_P_Proc%d",processIDInt);
    muP[muIndex][processIDInt] = new TH1D(histoid, histoid,80,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Pt_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Pt_Proc%d",processIDInt);
    muPt[muIndex][processIDInt] = new TH1D(histoid, histoid,80,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Px_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Px_Proc%d",processIDInt);
    muPx[muIndex][processIDInt] = new TH1D(histoid, histoid,80,-200.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Py_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Py_Proc%d",processIDInt);
    muPy[muIndex][processIDInt] = new TH1D(histoid, histoid,80,-200.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Pz_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Pz_Proc%d",processIDInt);
    muPz[muIndex][processIDInt] = new TH1D(histoid, histoid,80,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Eta_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Eta_Proc%d",processIDInt);
    muEta[muIndex][processIDInt] = new TH1D(histoid, histoid,30,-3.0,3.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Phi_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Phi_Proc%d",processIDInt);
    muPhi[muIndex][processIDInt] = new TH1D(histoid, histoid,30,-3.15,3.15);

    sprintf(histoid,"Muon_%d_FromZ_Eta_vs_Pt_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Eta_vs_Pt_Proc%d",processIDInt);
    muEtaPt[muIndex][processIDInt] = new TH2D(histoid, histoid,30,-3.15,3.15,50,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Isolation_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Isolation_Proc%d",processIDInt);
    muIso[muIndex][processIDInt] = new TH1D(histoid, histoid,100,0.0,50.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Eta_vs_ZMass_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Eta_vs_ZMass_Proc%d",processIDInt);
    ZCandidatesMassvsMuonEta[muIndex][processIDInt] = new TH2D(histoid, histoid,30,-3.15,3.15,90,0.,180.);
    
    /// Side Bands Plots
    sprintf(histoid,"Muon_%d_FromZ_P_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_P_SideBands_Proc%d",processIDInt);
    muP_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,80,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Pt_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Pt_SideBands_Proc%d",processIDInt);
    muPt_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,80,0.0,200.0);
	  
    sprintf(histoid,"Muon_%d_FromZ_Px_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Px_SideBands_Proc%d",processIDInt);
    muPx_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,80,-200.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Py_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Py_SideBands_Proc%d",processIDInt);
    muPy_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,80,-200.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Pz_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Pz_SideBands_Proc%d",processIDInt);
    muPz_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,80,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Eta_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Eta_SideBands_Proc%d",processIDInt);
    muEta_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,30,-3.0,3.0);
	 
    sprintf(histoid,"Muon_%d_FromZ_Phi_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Phi_SideBands_Proc%d",processIDInt);
    muPhi_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,30,-3.15,3.15);
    
    sprintf(histoid,"Muon_%d_FromZ_Eta_vs_Pt_SideBands_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Eta_vs_Pt_SideBands_Proc%d",processIDInt);
    muEtaPt_SB[muIndex][processIDInt] = new TH2D(histoid, histoid,30,-3.15,3.15,100,0.0,200.0);
    
    sprintf(histoid,"Muon_%d_FromZ_Isolation_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Isolation_SideBands_Proc%d",processIDInt);
    muIso_SB[muIndex][processIDInt] = new TH1D(histoid, histoid,100,0.0,100.0);

    sprintf(histoid,"Muon_%d_FromZ_Eta_vs_ZMass_Proc%d",muIndex,processIDInt);
    if(muIndex==2) sprintf(histoid,"Muons_FromZ_Eta_vs_ZMass_SideBands_Proc%d",processIDInt);
    ZCandidatesMassvsMuonEta_SB[muIndex][processIDInt] = new TH2D(histoid, histoid,30,-3.15,3.15,120,0.,180.);
  }
  
  sprintf(histoid,"Muons_FromZ_Eta_Correlation_Proc%d",processIDInt);
  muEtaCorrelation[processIDInt] = new TH1D(histoid, histoid,100,-6.0,6.0);
  
  sprintf(histoid,"Muons_FromZ_Phi_Correlation_Proc%d",processIDInt);
  muPhiCorrelation[processIDInt] = new TH1D(histoid, histoid,100,-7,7);
  
  
  sprintf(histoid,"Reco MET_Proc%d",processIDInt);	    
  MET[processIDInt]=new TH1D(histoid, "Reco METt" ,200,0.0,200.0);
  sprintf(histoid,"Reco METx_Proc%d",processIDInt);	    
  METx[processIDInt]=new TH1D(histoid, "Reco METxx",400,-200.0,200.0);
  sprintf(histoid,"Reco METy_Proc%d",processIDInt);	    
  METy[processIDInt]=new TH1D(histoid, "Reco METyy",400,-200.0,200.0); 
   
  sprintf(histoid,"ZGen MET_Proc%d",processIDInt);	    
  GenMET[processIDInt]=new TH1D(histoid, "Gen METt" ,200,0.0,200.0);
  sprintf(histoid,"Gen METx_Proc%d",processIDInt);	    
  GenMETx[processIDInt]=new TH1D(histoid, "Gen METxx",400,-200.0,200.0);
  sprintf(histoid,"Gen METy_Proc%d",processIDInt);	    
  GenMETy[processIDInt]=new TH1D(histoid, "Gen METyy",400,-200.0,200.0);  
  
  sprintf(histoid,"Gen NoNu MET_Proc%d",processIDInt);	    
  GenNoNuMET[processIDInt]=new TH1D(histoid, "Gen NoNu METt" ,200,0.0,200.0);
  sprintf(histoid,"Gen NoNu METx_Proc%d",processIDInt);	    
  GenNoNuMETx[processIDInt]=new TH1D(histoid, "Gen NoNu METxx",400,-200.0,200.0);
  sprintf(histoid,"Gen NoNu METy_Proc%d",processIDInt);	    
  GenNoNuMETy[processIDInt]=new TH1D(histoid, "Gen NoNu METyy",400,-200.0,200.0);  
  
  sprintf(histoid,"GenNuPt_Proc%d",processIDInt);	    
  NuPtGen[processIDInt]=new TH1D(histoid, "Nu Pt (GEN)",200,0.0,200.0);  
  sprintf(histoid,"NuPtResidual_Proc%d",processIDInt);	    
  NuPtResidual[processIDInt]=new TH1D(histoid, "Nu Pt REC-GEN",200,-100.0,100.0);
  
  sprintf(histoid,"HigPtJetsPT_Proc%d",processIDInt);	    
  JetPTAll[processIDInt]  = new TH1D(histoid, "JetPTAll",400,0.0,400.0);
  sprintf(histoid,"AHigPtJetsEta_Proc%d",processIDInt);	    
  JetEtaAll[processIDInt] = new TH1D(histoid, "JetEta",100,-3.0,3.0);
  sprintf(histoid,"HigPtJetsPhi_Proc%d",processIDInt);	    
  JetPhiAll[processIDInt] = new TH1D(histoid, "JetPhi",100,-3.2,3.2);
  sprintf(histoid,"HigPtJetsMultiplicity_Proc%d",processIDInt);	    
  JetMult[processIDInt]   = new TH1D(histoid, "JetMultiplicity",31,-0.5,30.5);
  sprintf(histoid,"HigPtJetsEta_vs_PT_Proc%d",processIDInt);	    
  JetEtaVSPTAll[processIDInt] = new TH2D(histoid, "JetEtaVsPt",100,-3.0,3.0,400,0.0,400.0);

  
}

/////////////////////////////
// Histogramming functions //
/////////////////////////////

void VecbosAnalysis::PlotsBeforeSelection() {
  if(verbose) std::cout<<"In VecbosAnalysis::PlotsBeforeSelection "<<std::endl;
  
  ProcessIDBeforeCuts->Fill(processIDIntOriginal, EventWeight); 

  // Gets DaughterIso1, DaughterIso2, and MET
  
  metXBeforecuts = pxMet[0];
  metYBeforecuts = pyMet[0];
  metBeforecuts = transverse(metXBeforecuts, metYBeforecuts);

  AllEventsMET[processIDInt]->Fill(metBeforecuts,EventWeight);
 
  // Loop over the muon isolations, in some way.
  // Still has to be implemented.
  {
    const double isolations = DUMMY;
    AllMuonsIso1[processIDInt]->Fill(isolations,EventWeight); 
    AllEventsMET_vs_Iso1[processIDInt]->Fill(metBeforecuts,isolations,EventWeight);
  }
 
  {
    const double isolations = DUMMY;
    AllMuonsIso2[processIDInt]->Fill(isolations,EventWeight); 
    AllEventsMET_vs_Iso2[processIDInt]->Fill(metBeforecuts,isolations,EventWeight);
  }
  
  treeBeforeCuts->Fill();
  
}

///////////////////////
// Working functions //
///////////////////////

/// Checks for W candidate.
bool VecbosAnalysis::doW(){
  if(verbose) std::cout<<"In VecbosAnalysis::doW "<<std::endl;
  TransWMass.clear();
  
  // Gets muons from W. I will use the collection of allMuons for now,
  // since this is what I understood from the .cfg files.
  // Also, gets MET.

  double metx= pxMet[0];
  double mety=pyMet[0];
  double metval=transverse(metx,mety);
   
  bool WCandFound=false;
  unsigned int countwcand=0;
  if(verbose) std::cout<<"Number of Candidates: "<< nMuon <<std::endl;
  
  for(unsigned int wIndex=0;wIndex<5;++wIndex) {
    WMuonPt[wIndex]=-100;
    WMuonPx[wIndex]=-100;
    WMuonPy[wIndex]=-100;
    WMuonEta[wIndex]=-100;
    WMuonPhi[wIndex]=-100;
    MuDauIsolation[wIndex]=-100;
  }
  
  // Loop over the muons;
  for(int muonCand = 0; muonCand != nMuon; ++muonCand) {
    //double isolationMuonFromW = (*DaughterIso[0])[muonRef];
    double isolationMuonFromW = DUMMY;
    //if(MuPt>LeptonLegsPtCut[0]){
    double MuPx=pxMuon[muonCand];
    double MuPy=pyMuon[muonCand];
    double MuPt=transverse(MuPx,MuPy);
    double tmass=transverse(MuPx+metx,MuPy+mety);
    TransWMass.push_back(  tmass );
    if(muonCand<5) {
      WMassTransverse[countwcand]=(float)tmass;
      WMuonPt[muonCand]=(float) MuPt;
      WMuonPx[muonCand]=(float) MuPx;
      WMuonPy[muonCand]=(float) MuPy;
      WMuonEta[muonCand]=(float) etaMuon[muonCand];
      WMuonPhi[muonCand]=(float) phiMuon[muonCand];
      MuDauIsolation[muonCand]=(float) isolationMuonFromW;
      if(verbose) std::cout<<"For Cand "<<countwcand<<"Pt ="<<MuPt <<std::endl;
    }
    BestWMuonPt=(float)MuPt;
    BestWMuonPx=(float)MuPx;
    BestWMuonPy=(float)MuPy;
    BestWMuonEta=(float)etaMuon[muonCand];
    BestWMuonPhi=(float)phiMuon[muonCand];    
    ++countwcand;
    WCandFound=true;
    //}
  }   
   
  if (countwcand>maxCandidatesCut)  WCandFound = false;
  
  return WCandFound;

}

/// Checks for Z candidate.
bool VecbosAnalysis::doZed(){
  if(verbose) std::cout<<"In VecbosAnalysis::doZed "<<std::endl;
   
  ZCandidates.clear();
   
  NumberoOfZedCand=0;   

  // Gets Z candidates. In the cfg files, this is the combination
  // of two muons (of collection allMuons), different charge, in
  // window 70 to 110. We can change this window in our settings file!
  
  // FIXME - add something to change the window.

  NumberoOfZedCand=nZ0ToMuMu;
  if(verbose)  std::cout<<"Found "<<NumberoOfZedCand<<" Z Candidate(s)!! "<<  std::endl;
  if(!NumberoOfZedCand) return false;

  double PdgZMass=91.3;
  double mass_previous=1000000.;
  int BestZ = -1;
  bool foundGoodZ=false;
  
  // Loop over Z candidates
  for(int zCand = 0; zCand != nZ0ToMuMu; ++zCand){
    int foundGoodLegs = 0;
    
    // Get the two daughters' indices.
    int daughters[2];
    
    double pt1 = transverse(pxMuon[d1IndexZ0ToMuMu[zCand]],pyMuon[d1IndexZ0ToMuMu[zCand]]);
    double pt2 = transverse(pxMuon[d2IndexZ0ToMuMu[zCand]],pyMuon[d2IndexZ0ToMuMu[zCand]]);
		
    // look for the most energetic muon
    if(pt1 >= pt2) {
      daughters[0] = d1IndexZ0ToMuMu[zCand];
      daughters[1] = d2IndexZ0ToMuMu[zCand];
    } else {
      daughters[0] = d2IndexZ0ToMuMu[zCand];
      daughters[1] = d1IndexZ0ToMuMu[zCand];
    }    
    // Loop over the daughters.
    for(int dauIndex =0; dauIndex<2;++dauIndex){
      double dau_px = pxMuon[daughters[dauIndex]];
      double dau_py = pyMuon[daughters[dauIndex]];
      double dau_pt = transverse(dau_px,dau_py);
      double dau_eta = etaMuon[daughters[dauIndex]];
      if(dau_pt< LeptonLegsPtCut[dauIndex])  break;  
      if(fabs(dau_eta)> LeptonLegsEtaCut[dauIndex]) break; 
      ++foundGoodLegs; 
    }
    
    // If we found a good Z candidate. 
    if(foundGoodLegs==2){
      
      foundGoodZ = true;
      if(UseBestCandidate){
	double ZCandmass = massZ0ToMuMu[zCand];
	if( (fabs(PdgZMass-ZCandmass))<(fabs(PdgZMass-mass_previous)) ) {
	  BestZ = zCand;
	  mass_previous=ZCandmass;
	}
      }else{
	ZCandidates.push_back(zCand);
      }
    }
  }
  
  if(UseBestCandidate) {
    if(foundGoodZ) ZCandidates.push_back(BestZ);
  }
  
  return foundGoodZ;
  
}

/// Checks for MET
bool VecbosAnalysis::doMET(){

  // In the framework, this function was simply getting the MET from the event.

  // Gets MET. In the cfg, this is simply ''met''.
  double met_px = pxMet[0];
  double met_py = pyMet[0];
  double met_pt = transverse(met_px,met_py);
  
  if(met_pt<metMinCut) return false;
  if(met_pt>metMaxCut) return false;

  // Gets GenMET. In the cfg, this is simply ''genMet''
  double genmet_px = pxGenMet[0];
  double genmet_py = pyGenMet[0];
  double genmet_pt = transverse(genmet_px, genmet_py);
  
  // Sorry, this collection is NOT in the Vecbos Ntuples (AFAIK).
  //   iEvent.getByLabel(_src_gen_met_nonu, genMetCollnoNU); 
  //   const GenMETCollection * genMETColnoNU = genMetCollnoNU.product();
  //   genMetnoNU = genMETColnoNU->front();
  
  return true;      
}

/// Checks for jets.
bool VecbosAnalysis::doRecoJets(){
  if(verbose) std::cout<<"In VecbosAnalysis::doRecoJets "<<std::endl;
  highPtJets.clear();
  bool foundJets = false;

  // Iterating over different jet algos.
  std::vector<std::string>::iterator jetAlgoName;
  for(jetAlgoName=src_jet.begin();jetAlgoName!=src_jet.end();++jetAlgoName){
   	
    // Gets the jets. How do we do if we have different collections of jets?
    // iEvent.getByLabel(*jetAlgoName, recoJets); 		
	
    if(nSisConeJet<jetMultiplicityCut) continue;
    unsigned int leadingJets(0);
    std::vector<int> selectedJets;
    
    // Loops over the jets, and selects the ones that pass the cuts.
    for(int itJet = 0;itJet!=nSisConeJet;++itJet){
      double ptJet_itJet = transverse(pxSisConeJet[itJet],pySisConeJet[itJet]);
      if(ptJet_itJet<jetPtCut)      continue;
      if(etaSisConeJet[itJet]<jetEtaMinCut) continue;
      if(etaSisConeJet[itJet]>jetEtaMaxCut) continue;
      selectedJets.push_back(itJet);
      ++leadingJets;
    }     

    if(leadingJets<jetMultiplicityCut) continue;
    if(jetExclusive) {if(leadingJets>jetMultiplicityCut) continue;}

    foundJets = true;
    // We save the selectedJets for a given algo.
    highPtJets[*jetAlgoName]=selectedJets;
  }	 

  return foundJets;

}

void VecbosAnalysis::fillHistoes(){
  if(verbose) std::cout<<"\n***\nFilling Histograms"<<std::endl;

  ////W////
  if(WSelection){
    WNumberofCandidates=TransWMass.size();
    NumberOfWCandidates[processIDInt]->Fill(WNumberofCandidates, EventWeight);
    for(std::vector<float>::iterator Wcand=TransWMass.begin();Wcand!=TransWMass.end();++Wcand){
      WTransvMass[processIDInt]->Fill(*Wcand,EventWeight);
    }
  }   

  ////Z////  
  if(ZSelection){
    ProcessIDAfterZed->Fill(processIDIntOriginal, EventWeight);
    NumberOfZCandidates[processIDInt]->Fill(ZCandidates.size(), EventWeight);
  
    std::vector<float> ZedCandPt;

    for(std::vector<int>::const_iterator it = ZCandidates.begin();it!=ZCandidates.end();++it){
      for(int mm=0;mm<2;++mm){
	ZMuonPt[mm]=-100;
	ZMuonPx[mm]=-100;
	ZMuonPy[mm]=-100;
	ZMuonEta[mm]=-100;
	ZMuonPhi[mm]=-100;
      }	
      double DiMuMass= massZ0ToMuMu[*it];
      double Zpx=pxZ0ToMuMu[*it];
      double Zpy=pyZ0ToMuMu[*it];
      double Zpt = transverse(Zpx,Zpy);
		
      ZCandidatesPt[processIDInt]->Fill(Zpt,EventWeight);
      ZedCandPt.push_back(Zpt);

      unsigned int muCount(0);
      double muonPhi[2];
      double muonEta[2];
      
      // What is she iterating over here? The daughters of the Z, I assume.
      // for(reco::Candidate::const_iterator muonFromZ = (*it)->begin(); muonFromZ< (*it)->end(); ++muonFromZ){
      // If that is the case, get the two daughters' indices.
      int daughters[2];
      daughters[0] = d1IndexZ0ToMuMu[*it];
      daughters[1] = d2IndexZ0ToMuMu[*it];
    
      for(int i = 0; i!= nDauZ0ToMuMu[*it]; ++i){
	if(muCount>1) {
	  std::cout<<"Too Many Muons from Z!!! "<<std::endl;
	  break;
	}
	
	int dauIndex = daughters[i];
	muonPhi[muCount]=phiMuon[dauIndex];
	muonEta[muCount]=etaMuon[dauIndex];
	double ptMuon = transverse(pxMuon[dauIndex],pyMuon[dauIndex]);
	// Know nothing about isolation.
	// double dauIsolation = (*DaughterIso[muCount])[muone];
	
	ZMuonPt[muCount]=ptMuon;
	ZMuonPx[muCount]=pxMuon[dauIndex];
	ZMuonPy[muCount]=pyMuon[dauIndex];
	ZMuonEta[muCount]=etaMuon[dauIndex];
	ZMuonPhi[muCount]=phiMuon[dauIndex];
	
	if((DiMuMass>ZMassRegion[0])&(DiMuMass<ZMassRegion[1])){
	  if(!muCount) ProcessIDInZedRegion->Fill(processIDIntOriginal, EventWeight);
	  
	  // Plots in the mass region.
	  // Know nothing about trigger.
	  //this->fillTriggerHist(TriggersInZMassRegion[processIDInt],TriggersInZMassRegion1D[processIDInt],0 );
				
	  //Individual Muon Plots
	  muP[muCount][processIDInt]->Fill(momentumMuon[dauIndex],EventWeight);
	  muPt[muCount][processIDInt]->Fill(ptMuon,EventWeight);
	  muPx[muCount][processIDInt]->Fill(pxMuon[dauIndex],EventWeight);
	  muPy[muCount][processIDInt]->Fill(pyMuon[dauIndex],EventWeight);
	  muPz[muCount][processIDInt]->Fill(pzMuon[dauIndex],EventWeight);
	  muEta[muCount][processIDInt]->Fill(etaMuon[dauIndex],EventWeight);
	  muPhi[muCount][processIDInt]->Fill(phiMuon[dauIndex],EventWeight);
	  muEtaPt[muCount][processIDInt]->Fill(etaMuon[dauIndex],ptMuon,EventWeight);
	  //muIso[muCount][processIDInt]->Fill(dauIsolation,EventWeight);
	  ZCandidatesMassvsMuonEta[muCount][processIDInt]->Fill(etaMuon[dauIndex],DiMuMass,EventWeight);
			
	  //Both Muons Plots
	  muP[2][processIDInt]->Fill(momentumMuon[dauIndex],EventWeight);
	  muPt[2][processIDInt]->Fill(ptMuon,EventWeight);
	  muPx[2][processIDInt]->Fill(pxMuon[dauIndex],EventWeight);
	  muPy[2][processIDInt]->Fill(pyMuon[dauIndex],EventWeight);
	  muPz[2][processIDInt]->Fill(pzMuon[dauIndex],EventWeight);
	  muEta[2][processIDInt]->Fill(etaMuon[dauIndex],EventWeight);
	  muPhi[2][processIDInt]->Fill(phiMuon[dauIndex],EventWeight);
	  muEtaPt[2][processIDInt]->Fill(etaMuon[dauIndex],ptMuon,EventWeight);
	  //muIso[2][processIDInt]->Fill(dauIsolation,EventWeight);
	  ZCandidatesMassvsMuonEta[2][processIDInt]->Fill(etaMuon[dauIndex],DiMuMass,EventWeight);
	}else{
			
	  // Plots in the sidebands.
	  //Individual Muon Plots
	  muP_SB[muCount][processIDInt]->Fill(momentumMuon[dauIndex],EventWeight);
	  muPt_SB[muCount][processIDInt]->Fill(ptMuon,EventWeight);
	  muPx_SB[muCount][processIDInt]->Fill(pxMuon[dauIndex],EventWeight);
	  muPy_SB[muCount][processIDInt]->Fill(pyMuon[dauIndex],EventWeight);
	  muPz_SB[muCount][processIDInt]->Fill(pzMuon[dauIndex],EventWeight);
	  muEta_SB[muCount][processIDInt]->Fill(etaMuon[dauIndex],EventWeight);
	  muPhi_SB[muCount][processIDInt]->Fill(phiMuon[dauIndex],EventWeight);
	  muEtaPt_SB[muCount][processIDInt]->Fill(etaMuon[dauIndex],ptMuon,EventWeight);
	  //muIso_SB[muCount][processIDInt]->Fill(dauIsolation,EventWeight);
	  ZCandidatesMassvsMuonEta_SB[muCount][processIDInt]->Fill(etaMuon[dauIndex],DiMuMass,EventWeight);
			
	  //Both Muons Plots
	  muP_SB[2][processIDInt]->Fill(momentumMuon[dauIndex],EventWeight);
	  muPt_SB[2][processIDInt]->Fill(ptMuon,EventWeight);
	  muPx_SB[2][processIDInt]->Fill(pxMuon[dauIndex],EventWeight);
	  muPy_SB[2][processIDInt]->Fill(pyMuon[dauIndex],EventWeight);
	  muPz_SB[2][processIDInt]->Fill(pzMuon[dauIndex],EventWeight);
	  muEta_SB[2][processIDInt]->Fill(etaMuon[dauIndex],EventWeight);
	  muPhi_SB[2][processIDInt]->Fill(phiMuon[dauIndex],EventWeight);
	  muEtaPt_SB[2][processIDInt]->Fill(etaMuon[dauIndex],ptMuon,EventWeight);
	  //muIso_SB[2][processIDInt]->Fill(dauIsolation,EventWeight);
	  ZCandidatesMassvsMuonEta_SB[2][processIDInt]->Fill(etaMuon[dauIndex],DiMuMass,EventWeight);
			
	}
	++muCount;
      }
		
      muEtaCorrelation[processIDInt]->Fill(muonEta[0]-muonEta[1]);
      muPhiCorrelation[processIDInt]->Fill(muonPhi[0]-muonPhi[1]);

      // Should be configurable.
      double etaLimit = 1.2;
      int region=0;
      if( (fabs(muonEta[0])<etaLimit) & (fabs(muonEta[1])<etaLimit) ) region =1;
      if( (fabs(muonEta[0])<etaLimit) & (fabs(muonEta[1])>etaLimit) ) region =2;
      if( (fabs(muonEta[0])>etaLimit) & (fabs(muonEta[1])<etaLimit) ) region =3;
      if( (fabs(muonEta[0])>etaLimit) & (fabs(muonEta[1])>etaLimit) ) region =4;

      ZCandidatesMass[0][processIDInt]->Fill(DiMuMass,EventWeight); // Always fill the all region
      ZCandidatesMass[region][processIDInt]->Fill(DiMuMass,EventWeight);

    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////                   MET
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(METSelection){
    if(verbose) std::cout<<"\n***\nFillingMET Histograms"<<std::endl;
    
    double ptMet = transverse(pxMet[0],pyMet[0]);
    
    MET[processIDInt]->Fill(ptMet,EventWeight);
    METx[processIDInt]->Fill(pxMet[0],EventWeight);
    METy[processIDInt]->Fill(pyMet[0],EventWeight);
    caloMET=ptMet;
    caloMETx=pxMet[0];
    caloMETy=pyMet[0];

    double ptGenMet = transverse(pxGenMet[0],pyGenMet[0]);

    GenMET[processIDInt]->Fill(ptGenMet,EventWeight);
    GenMETx[processIDInt]->Fill(pxGenMet[0],EventWeight);
    GenMETy[processIDInt]->Fill(pyGenMet[0],EventWeight);

    // I am not sure if the GenMet in the Ntuples has the neutrino or not. To discuss.
    //     GenNoNuMET[processIDInt]->Fill(genMetnoNU.pt(),EventWeight);
    //     GenNoNuMETx[processIDInt]->Fill(genMetnoNU.px(),EventWeight);
    //     GenNoNuMETy[processIDInt]->Fill(genMetnoNU.py(),EventWeight);
    
    //     NuPtGen[processIDInt]->Fill(genMet.pt()-genMetnoNU.pt(),EventWeight);
    //     NuPtResidual[processIDInt]->Fill(genMet.pt()-genMetnoNU.pt()-met.pt(),EventWeight);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////                   JETS
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(JetSelection){
    if(verbose) std::cout<<"\n***\nFillingJETS Histograms"<<std::endl;
    for(int jetind=0; jetind<15;++jetind){   
       JetsPt[jetind]=999;
       JetsEta[jetind]=999;
       JetsPhi[jetind]=999;
    }

    std::map<std::string,std::vector<int> >::iterator recoJetCollctnItr;
    for(recoJetCollctnItr=highPtJets.begin();recoJetCollctnItr!=highPtJets.end()
	  ;++recoJetCollctnItr){
      std::string jetAlgoName=recoJetCollctnItr->first;
	
      std::vector<int> selectedJets =recoJetCollctnItr->second;
      JetMultiplicity=selectedJets.size();
      JetMult[processIDInt]->Fill(JetMultiplicity,EventWeight);

      unsigned int jetCount(0);
      //float SumHiMultJetPx(0);float SumHiMultJetPy(0);
      std::vector<int>::const_iterator itJet;
      for(itJet=selectedJets.begin(); itJet!=selectedJets.end();++itJet){
	//if(verbose) std::cout<<(*itJet)->print();
	//JetPTAll[jetAlgoName]->Fill((*itJet)->pt());
	//JetEtaAll[jetAlgoName]->Fill((*itJet)->eta());
	//JetPhiAll[jetAlgoName]->Fill((*itJet)->phi());
	//JetEtaVSPTAll[jetAlgoName]->Fill((*itJet)->eta(),(*itJet)->phi());
	double ptJet = transverse(pxSisConeJet[*itJet],pySisConeJet[*itJet]);
	JetPTAll[processIDInt]->Fill(ptJet,EventWeight);
	JetEtaAll[processIDInt]->Fill(etaSisConeJet[*itJet],EventWeight);
	JetPhiAll[processIDInt]->Fill(phiSisConeJet[*itJet],EventWeight);
	JetEtaVSPTAll[processIDInt]->Fill(etaSisConeJet[*itJet],ptJet,EventWeight);
	if(jetCount<15){   
	  JetsPt[jetCount]=ptJet;
	  JetsEta[jetCount]=etaSisConeJet[*itJet];
	  JetsPhi[jetCount]=phiSisConeJet[*itJet];
	}
	//if(jetCount<jetMultiplicityCut){
	//SumHiMultJetPx+=(*itJet)->px();
	//SumHiMultJetPy+=(*itJet)->py();
	//}
	jetCount++;
	//float SumHiMultJetPt=sqrt(SumHiMultJetPx*SumHiMultJetPx + SumHiMultJetPy*SumHiMultJetPy);
	//HiMultJetSumPx[jetAlgoName]->Fill(SumHiMultJetPx);
	//HiMultJetSumPy[jetAlgoName]->Fill(SumHiMultJetPy);
	//HiMultJetSumPt[jetAlgoName]->Fill(SumHiMultJetPt);
	//for (std::vector<float>::iterator zCand=ZedCandPt.begin();zCand!=ZedCandPt.end();++zCand){
	//	ZSumPt_FromMuons_FromJets[jetAlgoName]->Fill((*zCand)-SumHiMultJetPt);
	//}
      }
    }
  }
  
  treeForFit->Fill();

}

///////////////////////
/// Setup functions ///
///////////////////////

//////////////////////
// InitParameters() //
//////////////////////
void VecbosAnalysis::InitParameters(string namefile) {

  useGlobalWeight=true;
  EventWeight=1.0;
  processIDInt=0;
  MassForFit=0.;
  processIDIntOriginal=0;
  numTriggers=0;

  soupType_ = "ALPGEN";
  outputHistoFileName= namefile;
  verbose=false; 
  
  WSelection=false;
  ZSelection=true;
  METSelection=false;
  JetSelection=true;
  // triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("TriggerResults");
  
  ///Initialization for Vec Boson
  UseBestCandidate=true;
  ZRecoFlagNumber=1;
  for(int ii=0;ii<3;++ii) foundZCandidate[ii]=false;
  checkOtherRecoTypes = false;
  ZMassRegion.push_back(30.);
  ZMassRegion.push_back(150.);
  // std::string Zsrc_string; (not used?)

  ///Initialization for MET
  metMinCut=10.0;
  metMaxCut=pow(10.,10.);
  
  ///Initialization for Jets
  jetAlgos=1;
  src_jet.clear();
  
  //for(unsigned int index=0;index<jetAlgos;++index){
  //  char algoName[20];
  //  sprintf(algoName,"src_jet_%d",index);
  //  src_jet.push_back(iConfig.getParameter<std::string>(algoName) ) ;
  //}

  src_jet.push_back("iterativeCone5");
  jetExclusive=false;
  jetMultiplicityCut=1;
  maxCandidatesCut=1;
  jetPtCut=30. ;
  jetEtaMinCut=0.;
  jetEtaMaxCut=3.;

  ///Initialization for Muons
  LeptonLegsPtCut.clear();
  LeptonLegsPtCut.push_back(15.);
  LeptonLegsPtCut.push_back(15.);
  LeptonLegsEtaCut.clear();
  LeptonLegsEtaCut.push_back(2.5);
  LeptonLegsEtaCut.push_back(2.5);
  
  TotalEventsProcessed=0;
  EventsAfterJetSelection=0;
  EventsAfterZedSelection=0;
}

/////////////////////////
// Ilaria's beginJob() //
/////////////////////////
void VecbosAnalysis::beginJob()
{
  if(verbose) std::cout<<"In VecbosAnalysis::beginJob "<<std::endl;
  
  for(int index =0; index<3; ++index) TypeOfZCounter[index]=0;
  
  fOutputFile   = new TFile( outputHistoFileName.c_str(), "RECREATE" ) ;

  //// TREE AFTER CUTS
  
  treeForFit= new TTree("ForFits", "Fit Variables");

  treeForFit->Branch("processIDIntOriginal", &processIDIntOriginal, "processIDIntOriginal");
  treeForFit->Branch("WNumberofCandidates", &WNumberofCandidates, "WNumberofCandidates");
  treeForFit->Branch("EventWeight", &EventWeightFloat, "EventWeight");
  treeForFit->Branch("DiMuMass", &MassForFit, "ZCandAllEvents");
  treeForFit->Branch("JetMultiplicity", &JetMultiplicity, "JetMultiplicity");
  treeForFit->Branch("HltBits", &HltBits, "HltBits[90]/F");
  
  treeForFit->Branch("WTransvMass", WMassTransverse, "WTransvMass");
  treeForFit->Branch("WMuonPt",  WMuonPt, "WMuonPt[5]/F");
  treeForFit->Branch("WMuonPx",  WMuonPx, "WMuonPx[5]/F");
  treeForFit->Branch("WMuonPy",  WMuonPy, "WMuonPy[5]/F");
  treeForFit->Branch("WMuonEta", WMuonEta, "WMuonEta[5]/F");
  treeForFit->Branch("WMuonPhi", WMuonPhi, "WMuonPhi[5]/F");
  treeForFit->Branch("WMuonIsolation", MuDauIsolation, "WMuonIsolation[5]/F");

  treeForFit->Branch("BestWMuonPt",  &BestWMuonPt, "BestWMuonPt");
  treeForFit->Branch("BestWMuonPx",  &BestWMuonPx, "BestWMuonPx");
  treeForFit->Branch("BestWMuonPy",  &BestWMuonPy, "BestWMuonPy");
  treeForFit->Branch("BestWMuonEta", &BestWMuonEta, "BestWMuonEta");
  treeForFit->Branch("BestWMuonPhi", &BestWMuonPhi, "BestWMuonPhi");
 
  treeForFit->Branch("caloMET", &caloMET, "caloMET");
  treeForFit->Branch("caloMETx", &caloMETx, "caloMETx");
  treeForFit->Branch("caloMETy", &caloMETy, "caloMETy");
  
  treeForFit->Branch("metMuCorrAftercuts", &metMuCorrAftercuts, "metMuCorrAftercuts");
  treeForFit->Branch("metXMuCorrAftercuts", &metXMuCorrAftercuts, "metMuXCorrAftercuts");
  treeForFit->Branch("metYMuCorrAftercuts", &metYMuCorrAftercuts, "metMuYCorrAftercuts");

  treeForFit->Branch("JetsPt",JetsPt, "JetsPt[15]/F");
  treeForFit->Branch("JetsEta", JetsEta, "JetsEta[15]/F");
  treeForFit->Branch("JetsPhi", JetsPhi, "JetsPhi[15]/F");

  treeForFit->Branch("ZMuonPt",  ZMuonPt, "ZMuonPt[2]/F");
  treeForFit->Branch("ZMuonPx",  ZMuonPx, "ZMuonPx[2]/F");
  treeForFit->Branch("ZMuonPy",  ZMuonPy, "ZMuonPy[2]/F");
  treeForFit->Branch("ZMuonEta", ZMuonEta, "ZMuonEta[2]/F");
  treeForFit->Branch("ZMuonPhi", ZMuonPhi, "ZMuonPhi[2]/F");
  //// TREE BEFORE CUTS
  
  treeBeforeCuts= new TTree("treeBeforeCuts", "Before Any Selection");
  
  treeBeforeCuts->Branch("metBeforecuts", &metBeforecuts, "metBeforecuts");
  treeBeforeCuts->Branch("metMuCorrBeforecuts", &metMuCorrBeforecuts, "metMuCorrBeforecuts");
  treeBeforeCuts->Branch("metXBeforecuts", &metXBeforecuts, "metXBeforecuts");
  treeBeforeCuts->Branch("metXMuCorrBeforecuts", &metXMuCorrBeforecuts, "metXMuCorrBeforecuts");
  treeBeforeCuts->Branch("metYBeforecuts", &metYBeforecuts, "metYBeforecuts");
  treeBeforeCuts->Branch("metYMuCorrBeforecuts", &metYMuCorrBeforecuts, "metYMuCorrBeforecuts");
  
  Weight_vs_ProcessID = new TH1D("Weight_vs_ProcessID", "Weight_vs_ProcessID",4001,-0.5,4000.5);

  ProcessIDBeforeCuts     = new TH1D("ProcessIDBeforeCuts", "ProcessIDBeforeCuts",4001,-0.5,4000.5);
  ProcessIDAfterJets      = new TH1D("ProcessIDAfterJets", "ProcessIDAfterJets Cuts",4001,-0.5,4000.5);
  ProcessIDAfterZed       = new TH1D("ProcessIDAfterZed", "ProcessIDAfterZed Cuts",4001,-0.5,4000.5);
  ProcessIDInZedRegion    = new TH1D("ProcessIDInZedRegion", "ProcessIDInZedRegion Cuts",4001,-0.5,4000.5);

  NumberOfEvents = new TH1D("NumberOfEvents", "NumberOfEvents",3,0.5,3.5);	
  
  NumberoOfZedCand=0;
}

///////////////////////
// Ilaria's endJob() //
///////////////////////
void VecbosAnalysis::endJob() {

  for(int bin=0;bin<3;++bin)  TypeOfZCandidates[processIDInt]->SetBinContent( bin+1, TypeOfZCounter[bin] );
  
  NumberOfEvents->SetBinContent( 1, TotalEventsProcessed );
  NumberOfEvents->SetBinContent( 2, EventsAfterJetSelection );
  NumberOfEvents->SetBinContent( 3, EventsAfterZedSelection );

  fOutputFile->Write() ;
  fOutputFile->Close() ;
}

/////////////////////
// AssignParameters//
/////////////////////
void VecbosAnalysis::AssignParameters(map<string, double> data) {
  // Should put code to assign parameters here.
}

///////////////
// Utilities //
///////////////

double VecbosAnalysis::transverse(double x, double y) {
  return sqrt(x*x+y*y);
}
