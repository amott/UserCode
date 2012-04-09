// -*- C++ -*-
//
// Package:    RAzrTO
// Class:      RAzrTO
// 
/**\class RAzrTO RAzrTO.cc TriggerTurnOn/RAzrTO/src/RAzrTO.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alex Mott,32 3-A13,+41227677936,
//         Created:  Tue Oct 18 11:14:14 CEST 2011
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <tuple>
//#include <boost/tuple/tuple.hpp>

// std root include files
#include <TH3D.h>
#include <TFile.h>
#include <TString.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class RAzrTO : public edm::EDAnalyzer {
   public:
      explicit RAzrTO(const edm::ParameterSet&);
      ~RAzrTO();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  
  void computeRMR(const edm::Event&);
  void passTrigger(const edm::Event&);
      // ----------member data ---------------------------
  edm::InputTag t_TrigResults;
  edm::InputTag t_JetCollection;
  edm::InputTag t_MetCollection;
  edm::InputTag t_ThirdObject;
  edm::InputTag t_VertexCollection;
  edm::InputTag t_NoiseFilter;
  std::string s_FileName;
  std::vector< std::string > vs_TrigNames;
  std::vector< int > vi_TrigDec;
  std::vector< float > vf_RCuts, vf_MRCuts, vf_ThirdObjectCuts;
  std::vector< std::string > vs_PreSelection;
  std::vector< std::int > vi_PreSelDec;
  TFile *outputFile;
  TH3D* th3d_total;
  TH3D* th3d_noisetotal;
  TH3D* th3d_noisepvtotal;
  std::vector< TH3D* > vth3d_pass; 

  double R,MR,objPT;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RAzrTO::RAzrTO(const edm::ParameterSet& iConfig):
  t_TrigResults(iConfig.getParameter<edm::InputTag>("TriggerResults")),
  t_JetCollection(iConfig.getParameter<edm::InputTag>("JetCollection")),
  t_MetCollection(iConfig.getParameter<edm::InputTag>("MetCollection")),
  t_ThirdObject(iConfig.getParameter<edm::InputTag>("ThirdObject")),
  t_VertexCollection(iConfig.getParameter<edm::InputTag>("VertexCollection")),
  t_NoiseFilter(iConfig.getParameter<edm::InputTag>("NoiseTag")),
  vs_TrigNames(iConfig.getParameter<std::vector< std::string > >("TrigNames")),
  vf_RCuts(iConfig.getParameter<std::vector< std::string > >("RCuts")),
  vf_MRCuts(iConfig.getParameter<std::vector< std::string > >("MRCuts")),
  vf_ThirdObjectCuts(iConfig.getParameter<std::vector< std::string > >("ThirdObjectCuts")),
  vs_PreSelection(iConfig.getParameter<std::vector< std::string > >("PreSelectionTriggers")),
  s_FileName(iConfig.getParameter<std::string>("OutputFileName")),
{
  outputFile = new TFile(s_FileName.c_str(),"RECREATE");
   //now do what ever initialization is needed  
}


RAzrTO::~RAzrTO()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  outputFile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
RAzrTO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get collections
   Handle< reco::CandidateCollection > h_Obj;
   Handle< reco::VertexCollection > h_Vtx;
   Handle< bool > h_NoiseFilter;

   // do the preselection

   iEvent.getObjectByLabel(t_TrigResults,h_Trg);
   if(t_ThirdObject != "none")
     iEvent.getObjectByLabel(t_ThirdObject,h_Obj);
   iEvent.getObjectByLabel(t_VertexCollection,h_Vtx);
   iEvent.getObjectByLabel(t_NoiseFilter,h_NoiseFilter);



   R=-1; MR=-1; objPT=0;
   this->computeRMR();
   if(R==-1 || MR==-1) return; // no valid hemispheres found!

   if(h_Obj.isValid()){
     objPT = h_Obj->front().pt();
   }
   th3d_total.Fill(R,MR,objPT);

   if(! (*h_NoiseFilter) ) return; // failed the hcal noise filter

   th3d_noisetotal.Fill(R,MR,objPT);
   
   if(!h_Vtx->size()>0) return; // no good vertex

   th3d_noisepvtotal.Fill(R,MR,objPT);

   std::vector< std::string >::iterator trigIter;
   

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
RAzrTO::beginJob()
{
  assert( vf_RCuts.size() == vf_MRCuts.size() );
  assert( vf_ThirdObjectCuts.size()==0 || (vf_ThirdObjectCuts.size() == vfMRSuts.size())  );

  th3d_total = new TH3D("total","",200,0,2,2000,0,2000,200,0,200); // R/MR/Obj  --> make bins much too fine to rebin later
  std::vector< std::string >::iterator vsit_TrigNames;
  for(vsit_TrigNames  = vs_TrigNames.begin(); 
      vsit_TrigNames != vs_TrigNames.end();
      vsit_TrigNames++){
    vth3d_pass.push_back( new TH3D(Form("%s_pass",vsit_TrigNames->c_str()),"",200,0,2,2000,0,2000,200,0,200) ); // R/MR/Obj  --> make bins much too fine to rebin later
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RAzrTO::endJob() 
{
  th3d_total->Write();
  for(int i=0;i<vth3d_pass.size();i++) vth3d_pass.at(i)->Write();
}

// ------------ method called when starting to processes a run  ------------
void 
RAzrTO::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RAzrTO::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RAzrTO::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RAzrTO::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RAzrTO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
RAzrTO::computeRMR(const edm::Event& iEvent, const edm::EventSetup& iSetup){
   Handle< reco::JetCollection > h_Jet;
   Handle< reco::MetCollection > h_Met;
   iEvent.getObjectByLabel(t_JetCollection,h_Jet);
   iEvent.getObjectByLabel(t_MetCollection,h_Met);
   
   reco::JetCollection passingJets;
   
   reco::JetCollection::const_iterator JetIter;
   for(JetIter = h_Jet->begin(); JetIter !=h_Jet->end(); JetIter++){
     if( fabs( JetIter->eta() ) > 3. ||
	 JetIter->pt() < 40. ) continue;
     passingJets.push_back(*JetIter);
   }

   const int n = passingJets.size();

   if(n<2) return;

   int N_comb(1); // compute the number of combinations of jets possible
   for(unsigned int i = 0; i < passingJets.size(); i++){
     N_comb *= 2;                
   }
  //Make the hemispheres
  XYZTLorentzVector j1R(0.1, 0., 0., 0.1);
  XYZTLorentzVector j2R(0.1, 0., 0., 0.1);
  double M_minR = 9999999999.0;
  int j_count;
  for (int i = 0; i < N_comb; i++) {       
    XYZTLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while (j_count > 0) {
      if (itemp/j_count == 1){
	j_temp1 += passingJets.at(count).p4();
      } else {
	j_temp2 += passingJets.at(count).p4();
      }
      itemp -= j_count * (itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.M2() + j_temp2.M2();
    if (M_temp < M_minR) {
      M_minR = M_temp;
      j1R = j_temp1; 
      j2R = j_temp2; 
    }
  }

  //j1R and j2R are now the hemisphere vectors

  //compute R and MR

  j1R.SetPtEtaPhiM(j1R.Pt(),j1R.Eta(),j1R.Phi(),0.0);
  j2R.SetPtEtaPhiM(j2R.Pt(),j2R.Eta(),j2R.Phi(),0.0);
  
  if(j1R.Pt() > j2R.Pt()){
    TLorentzVector temp = j1R;
    j1R = j2R;
    j2R = temp;
  }
  
  double A = j1R.P();
  double B = j2R.P();
  double az = j1R.Pz();
  double bz = j2R.Pz();
  TVector3 j1RT, j2RT;
  j1RT.SetXYZ(j1R.Px(),j1R.Py(),0.0);
  j2RT.SetXYZ(j2R.Px(),j2R.Py(),0.0);
  double ATBT = (j1RT+j2RT).Mag2();
  
  //set MR
  MR = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
		   (j2RT.Dot(j2RT)-j1RT.Dot(j1RT))*(j2RT.Dot(j2RT)-j1RT.Dot(j1RT))/(j1RT+j2RT).Mag2());
  
  double mybeta = (j2RT.Dot(j2RT)-j1RT.Dot(j1RT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  
  double mygamma = 1./sqrt(1.-mybeta*mybeta);
  
  //use gamma times MRstar
  MR *= mygamma;
  
  //now we can calculate MTR
  TVector3 met;
  met.SetPtEtaPhi((inputMet->front()).pt(),0.0,(inputMet->front()).phi());
  double MTR = sqrt(0.5*(met.Mag()*(j1R.Pt()+j2R.Pt()) - met.Dot(j1R.Vect()+j2R.Vect())));
  
  //set R
  R = float(MTR)/float(MR);
  
  
}

void
RAzrTO::passTrigger(const edm::Event& iEvent){

   // get trigger names
   Handle< edm::TriggerResults > h_Trg;
   const edm::TriggerNames & triggerNames = e.triggerNames(*triggerResults);


   if(vs_PreSelection.size() == 0) vi_PreSelDec.puh_back(1);
   else{
     std::vector< std::string >::iterator preIter;
     for(preIter = vs_PreSelection.begin(); preIter !=vs_PreSelection.end(); preIter++){
       h_Trg->accept(triggerNames.triggerIndex(*preIter));
     }
   }
}


//define this as a plug-in
DEFINE_FWK_MODULE(RAzrTO);
