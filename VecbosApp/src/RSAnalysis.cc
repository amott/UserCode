/// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/// ROOT includes
#include <TTree.h>
#include <TH3D.h>
#include <TVector3.h>
#include <TLorentzVector.h>

/// local includes
#include "VecbosBase.hh"
#include "RSAnalysis.hh"

RSAnalysis::RSAnalysis(TTree *tree) : Vecbos(tree) {
}

RSAnalysis::~RSAnalysis(){
  
}

void RSAnalysis::Loop(string outfile) {
  if(fChain == 0) return;

  vector<TH1D*> Hpt1;
  vector<TH1D*> Hpt2;

  vector<TH3D*> Hmass1;
  vector<TH3D*> Hmass2;
  vector<TH3D*> HmassRS;
  vector<TH3D*> HnJets;

  vector<TH3D*> HDeltamass;
  vector<TH3D*> HDeltapT;
  vector<TH3D*> HnTrk1;
  vector<TH3D*> HnTrk2;
  
  vector<TH2D*> Hmass;


  vector<string>labels;
  labels.push_back("_dr_0.1");
  labels.push_back("_dr_0.2");
  labels.push_back("_dr_0.3");
  labels.push_back("_dr_0.4");
  labels.push_back("_dr_0.5");
  labels.push_back("_dr_0.6");
  labels.push_back("_dr_0.7");
  labels.push_back("_dr_0.8");
  labels.push_back("_dr_0.9");
  labels.push_back("_dr_1.0");

  for(int i=0; i< int(labels.size()); i++) {
    Hpt1.push_back(new TH1D(string("Hpt1"+labels[i]).c_str(), string("Hpt1"+labels[i]).c_str(), 100, 0., 1500.));
    Hpt2.push_back(new TH1D(string("Hpt2"+labels[i]).c_str(), string("Hpt2"+labels[i]).c_str(), 100, 0., 1500.));

    Hmass1.push_back(new TH3D(string("Hmass1"+labels[i]).c_str(), string("Hmass1"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 100, 0., 500.));
    Hmass2.push_back(new TH3D(string("Hmass2"+labels[i]).c_str(), string("Hmass2"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 100, 0., 500.));
    HmassRS.push_back(new TH3D(string("HmassRS"+labels[i]).c_str(), string("HmassRS"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 200, 100., 2500.));
    HnJets.push_back(new TH3D(string("HnJets"+labels[i]).c_str(),  string("HnJets"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 20, 0., 20.));

    HDeltamass.push_back(new TH3D(string("HDeltamass"+labels[i]).c_str(), string("Deltamass"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 100, -300., 300.));    
    HDeltapT.push_back(new TH3D(string("HDeltapT"+labels[i]).c_str(), string("DeltapT"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 100, -300., 300.));    
    HnTrk1.push_back(new TH3D(string("HnTrk1"+labels[i]).c_str(), string("nTrk1"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 100, 0., 100.));    
    HnTrk2.push_back(new TH3D(string("HnTrk2"+labels[i]).c_str(), string("nTrk2"+labels[i]).c_str(), 100, 0., 1500., 100, 0., 1500., 100, 0., 100.));    
    
    Hmass.push_back(new TH2D(string("Hmass"+labels[i]).c_str(),  string("Hmass"+labels[i]).c_str(),  100, 0., 500., 100, 0., 500.));

  }

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
   
    vector<float> thresh;
    thresh.push_back(.9);
    thresh.push_back(1.1);
    thresh.push_back(1.4);
    thresh.push_back(1.4);
    thresh.push_back(1.2);
    thresh.push_back(1.8);
    thresh.push_back(.09);
    thresh.push_back(.45);
    thresh.push_back(.2);
    thresh.push_back(.45);
    
    double pt_cut = 30.0;
    CaloF = 0.06;

    if(nPV < 1) continue;

    iPV = -1;
    double maxpt = -1;
    for(int i = 0; i < nPV; i++){
      if(SumPtPV[i] > maxpt){
	maxpt = SumPtPV[i];
	iPV = i;
      }
    }

    //Create CaloTowers
    c_uncorr = CreateCaloTowers(thresh, 0.0, 0);
    
    // Make my jets
    
    vector<Jet> j_u;
    vector<Jet> j_u2;    
    vector<vector<Jet> > jets;
    
    for(int i=0; i<int(labels.size()); i++) {
      
      j_u.clear();
      j_u2.clear();
      double dR = 0.1*(i+1);
      
      j_u = SISCone(c_uncorr, dR, 0.0);
      j_uncorr_all = SortJet(j_u);
      
      for(int i = 0; i < j_uncorr_all.size();i++){
	j_u2.push_back(CorrectJet(j_uncorr_all[i], PVzPV[iPV]));
	if(j_uncorr_all[i].pt() > pt_cut && fabs(j_uncorr_all[i].eta()) < 3.0){
	  j_uncorr.push_back(j_uncorr_all[i]);
	}
      }
      
      jets.push_back(j_uncorr);
      j_uncorr.clear();
      j_uncorr_all.clear();
      
    }
    
    for(int i=0; i<int(labels.size()); i++) {
      
      if(int(jets[i].size()) < 2) continue; 

      // The jets are sorted. The pT and eta cuts are applied
      // let's make the plots with the first two
      
      double pt1 = jets[i][0].pt();
      double pt2 = jets[i][1].pt();
      double m1  = jets[i][0].mass();
      double m2  = jets[i][1].mass();
      double mRS = (jets[i][0].Sum(jets[i][1])).M();

      int ntk1 =0;
      int ntk2 =0;
      for(int iTk=0; iTk<nTrack; iTk++) {
	if(DeltaR(jets[i][0].eta(), jets[i][0].phi(), double(etaTrack[iTk]), double(phiTrack[iTk]))< (i+1)*0.1) ntk1++;
	if(DeltaR(jets[i][1].eta(), jets[i][1].phi(), double(etaTrack[iTk]), double(phiTrack[iTk]))< (i+1)*0.1) ntk2++;
      }

      Hpt1[i]->Fill(pt1);
      Hpt2[i]->Fill(pt2);
      Hmass1[i]->Fill(pt1,pt2,m1);
      Hmass2[i]->Fill(pt1,pt2,m2);
      Hmass[i]->Fill(m1,m2);
      HmassRS[i]->Fill(pt1,pt2,mRS);
      HnJets[i]->Fill(pt1,pt2,jets[i].size());

      HDeltamass[i]->Fill(pt1,pt2,m1-m2);
      HDeltapT[i]->Fill(pt1,pt2,pt1-pt2);
      HnTrk1[i]->Fill(pt1,pt2,ntk1);
      HnTrk2[i]->Fill(pt1,pt2,ntk2);

    }
  }
  
  TFile* outFile = new TFile(outfile.c_str(),"recreate");
  for(int i=0; i<int(labels.size()); i++) {
    Hpt1[i]->Write();
    Hpt2[i]->Write();  
    Hmass1[i]->Write();
    Hmass2[i]->Write();
    Hmass[i]->Write();
    HmassRS[i]->Write();
    HnJets[i]->Write();

    HDeltamass[i]->Write();
    HDeltapT[i]->Write();
    HnTrk1[i]->Write();
    HnTrk2[i]->Write();

  }
  outFile->Close();
  
}

