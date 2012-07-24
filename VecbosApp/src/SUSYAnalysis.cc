// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "SUSYAnalysis.hh"

SUSYAnalysis::SUSYAnalysis(TTree *tree) : Vecbos(tree) {}

SUSYAnalysis::~SUSYAnalysis(){}

vector<TH2D*> SUSYAnalysis::CreateHistos2D(string dirname){
  vector<TH2D*> h_2D;

  string name;

  name = dirname+"_ECHF_v_EEMF"; //0
  h_2D.push_back(new TH2D(name.c_str(),name.c_str(), 200, 0.0, 1.0, 200, 0.0, 2.5));
  name = dirname+"_dj1_v_dj2"; //1
  h_2D.push_back(new TH2D(name.c_str(),name.c_str(), 200, 0.0, TMath::Pi(), 200, 0.0, TMath::Pi()));
  return h_2D;
}

vector<TH1D*> SUSYAnalysis::CreateHistosD(string dirname){
  vector<TH1D*> h_d;

  string name;

  name = dirname+"_EEMF"; //0
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,1.0));
  name = dirname+"_ECHF"; //1
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,2.5));
  name = dirname+"_j1_EEMF"; //2
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,1.0));
  name = dirname+"_j1_ECHF"; //3
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,2.5));
  name = dirname+"_j1_Ntrack"; //4
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),20,0.,40));
  name = dirname+"_j1_et"; //5
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.,400));
  name = dirname+"_j2_EEMF"; //6
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,1.0));
  name = dirname+"_j2_ECHF"; //7
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,2.5));
  name = dirname+"_j2_Ntrack"; //8
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),20,0.,40));
  name = dirname+"_j2_et"; //9
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.,400));
  name = dirname+"_j3_EEMF"; //10
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,1.0));
  name = dirname+"_j3_ECHF"; //11
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),200,0.,2.5));
  name = dirname+"_j3_Ntrack"; //12
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),20,0.,40));
  name = dirname+"_j3_et"; //13
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.,400));

  name = dirname+"_j1_dphi"; //14
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.,TMath::Pi()));
  name = dirname+"_j2_dphi"; //15
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.,TMath::Pi()));
  name = dirname+"_j3_dphi"; //16
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.,TMath::Pi()));
  name = dirname+"_dMET_o_MET"; //17
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),1000,-5.0,5.0));
  name = dirname+"_dMETx"; //18
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),2000,-250.,250.));
  name = dirname+"_dMETy"; //19
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),2000,-250.,250.));
  name = dirname+"_dMETphi"; //20
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),500,0.0, TMath::Pi()));
  name = dirname+"_MET"; //21
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),2000,0.0, 2000.));
  
  name = dirname+"_EFF_num"; //22
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),12,0.0, 12.0));
  name = dirname+"_EFF_cum"; //23
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),12,0.0, 12.0));
  name = dirname+"_EFF_den"; //24
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),12,0.0, 12.0));
  name = dirname+"_dRj1"; //25
  h_d.push_back(new TH1D(name.c_str(),name.c_str(),2000,0.0, 5.0));
  return h_d;
  
}

vector<TProfile*> SUSYAnalysis::CreateHistosT(string dirname){
  vector<TProfile*> h_p;
  string name;

  name = dirname+"_j1_EMF_v_et"; //0
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,1.0));
  name = dirname+"_j1_CHF_v_et"; //1
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,4.5));
  name = dirname+"_j1_Ntrack_v_et"; //2
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  
  name = dirname+"_j2_EMF_v_et"; //3
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,1.0));
  name = dirname+"_j2_CHF_v_et"; //4
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,4.5));
  name = dirname+"_j2_Ntrack_v_et"; //5
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  name = dirname+"_jall_CHF_v_et_0"; //6
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  name = dirname+"_jall_CHF_v_et_1"; //7
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  name = dirname+"_jall_CHF_v_et_2"; //8
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  name = dirname+"_jall_CHF_v_et_3"; //9
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  name = dirname+"_jall_CHF_v_et_4"; //10
  h_p.push_back(new TProfile(name.c_str(),name.c_str(),50,0.0,400.,0.0,50.0));
  return h_p;


}

void SUSYAnalysis::FillHistos(vector<TH1D*> h_d, vector<TProfile*> h_p, vector<TH2D*> h_2D){
  double MET_phi = m_dum[0].phi();
  if((m_dum[0].met()+m_gen[0].met()) > 0.0){
    double dMET = (m_dum[0].met()-m_gen[0].met())/(m_dum[0].met()+m_gen[0].met());
    h_d[17]->Fill(dMET);
  }
  h_d[18]->Fill(m_dum[0].mex()-m_gen[0].mex());
  h_d[19]->Fill(m_dum[0].mey()-m_gen[0].mey());
  if(energyGenMet[0] > 0.0){
    h_d[20]->Fill(m_dum[0].phi()-m_gen[0].phi());
  }
  h_d[21]->Fill(m_dum[0].met());

  h_2D[0]->Fill(EEMF, ECHF);
  h_d[0]->Fill(EEMF);
  h_d[1]->Fill(ECHF);
  if(j_dum.size()){
    int Ntrack;
    double jCHF = JetCHF(j_dum[0], 0.75, Ntrack);
    h_d[2]->Fill(j_dum[0].EmFrac());
    h_d[3]->Fill(jCHF);
    h_d[4]->Fill(double(Ntrack));
    h_d[5]->Fill(j_dum[0].pt());
    h_p[0]->Fill(j_dum[0].pt(), j_dum[0].EmFrac());
    h_p[1]->Fill(j_dum[0].pt(), jCHF);
    h_p[2]->Fill(j_dum[0].pt(), double(Ntrack));
    h_d[14]->Fill(fabs(DeltaPhi(j_dum[0].phi(), MET_phi)));
    TVector3 vPV(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
    for(int i = 0; i < nTrack; i++){
      TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
      vT = vT-vPV;
      if(fabs(vT.z()) > 0.1) continue;
      if(fabs(vT.Mag()) > 1.0) continue;
      //if(abs(chargeTrack[i]) > 1) continue;
      double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
      if(pt < 1.2) continue;
      if(pt > 500.) continue;
      if(trackNormalizedChi2Track[i] > 20.0) continue;
      if(trackDxyPVTrack[i] > .6) continue;
      if(fabs(etaTrack[i]) > 2.4) continue;
      if(trackValidHitsTrack[i] < 5) continue;
      double dR = DeltaR(double(etaTrack[i]), double(phiTrack[i]), j_dum[0].eta(), j_dum[0].phi());
      h_d[25]->Fill(dR);
    }
  }
  if(j_dum.size() > 1){
    int Ntrack;
    double jCHF = JetCHF(j_dum[1], 0.75, Ntrack);
    h_p[3]->Fill(j_dum[1].pt(), j_dum[1].EmFrac());
    h_p[4]->Fill(j_dum[1].pt(), jCHF);
    h_p[5]->Fill(j_dum[1].pt(), double(Ntrack));
    h_d[6]->Fill(j_dum[1].EmFrac());
    h_d[7]->Fill(jCHF);
    h_d[8]->Fill(double(Ntrack));
    h_d[9]->Fill(j_dum[1].pt());
    h_d[15]->Fill(fabs(DeltaPhi(j_dum[1].phi(), MET_phi)));
    h_2D[1]->Fill(fabs(DeltaPhi(j_dum[0].phi(), MET_phi)),fabs(DeltaPhi(j_dum[1].phi(), MET_phi)));
  }
  if(j_dum.size() > 2){
    int Ntrack;
    double jCHF = JetCHF(j_dum[2], 0.75, Ntrack);
    h_d[10]->Fill(j_dum[2].EmFrac());
    h_d[11]->Fill(jCHF);
    h_d[12]->Fill(double(Ntrack));
    h_d[13]->Fill(j_dum[2].pt());
    h_d[16]->Fill(fabs(DeltaPhi(j_dum[2].phi(), MET_phi)));
  }
  for(int i = 0; i < j_dum.size(); i++){
    for(int j = 0; j < 5; j++){
      int Ntrack;
      double dR = 0.5+double(j)*0.25;
      double jCHF = JetCHF(j_dum[i], dR, Ntrack);
      h_p[6+j]->Fill(j_dum[i].pt(), double(Ntrack));
    }
  }

  bool EFF[12];
  for(int i = 0; i < 12; i++){
    EFF[i] = false;
  }
  if(m_dum[0].met() > 200.0) EFF[0] = true;
  if(j_dum.size() >= 3) EFF[1] = true;
  if(j_dum.size()){
    if(fabs(j_dum[0].eta()) <= 1.7) EFF[2] = true;
    EFF[5] = true;
  }
  if(EEMF >=0.174) EFF[3] = true;
  if(ECHF >= 0.1) EFF[4] = true;
  double Ht = 0.0;
  for(int i = 0; i < j_dum.size(); i++){
    if(fabs(DeltaPhi(MET_phi, j_dum[i].phi())) < 0.3) EFF[5] = false;
    if(i > 0 && i < 4) Ht += j_dum[i].pt();
  }
  Ht += m_dum[0].met();
  if(Ht > 500.0) EFF[11] = true;
  if(j_dum.size() > 1){
    if(fabs(DeltaPhi(MET_phi, j_dum[1].phi())) > TMath::Pi()*(20./180.)) EFF[6] = true;
    double dp1 = fabs(DeltaPhi(MET_phi, j_dum[0].phi()));
    double dp2 = fabs(DeltaPhi(MET_phi, j_dum[1].phi()));
    double R1 = sqrt(dp2*dp2 + (TMath::Pi()-dp1)*(TMath::Pi()-dp1));
    double R2 = sqrt(dp1*dp1 + (TMath::Pi()-dp2)*(TMath::Pi()-dp2));
    if(R1 > 0.5 && R2 > 0.5) EFF[7] = true;
    if(j_dum[0].EmFrac() < 0.9 && j_dum[1].EmFrac() < 0.9) EFF[9] = true;
    if(j_dum[0].pt() > 180. && j_dum[1].pt() < 110.) EFF[10] = true;
  }
  if(ISO == false) EFF[8] = true;
  bool cum = true;
  for(int i = 0; i < 12; i++){
    h_d[24]->Fill(double(i)+0.5);
    if(EFF[i]) h_d[22]->Fill(double(i)+0.5);
    if(EFF[i] && cum) h_d[23]->Fill(double(i)+0.5);
    if(EFF[i] == false) cum = false;
  }
  
  
    
  
}

void SUSYAnalysis::Loop() {
  if(fChain == 0) return;

  vector< vector<TH1D*> > H_d;
  vector< vector<TProfile*> > H_t;
  vector< vector<TH2D*> > H_2D;

 
  H_d.push_back(CreateHistosD("uncorr_IC"));
  H_t.push_back(CreateHistosT("uncorr_IC"));
  H_2D.push_back(CreateHistos2D("uncorr_IC"));

  H_d.push_back(CreateHistosD("corr_0_IC"));
  H_t.push_back(CreateHistosT("corr_0_IC"));
  H_2D.push_back(CreateHistos2D("corr_0_IC"));
  
  H_d.push_back(CreateHistosD("corr_1_IC"));
  H_t.push_back(CreateHistosT("corr_1_IC"));
  H_2D.push_back(CreateHistos2D("corr_1_IC"));

  H_d.push_back(CreateHistosD("corr_2_IC"));
  H_t.push_back(CreateHistosT("corr_2_IC"));
  H_2D.push_back(CreateHistos2D("corr_2_IC"));

  H_d.push_back(CreateHistosD("corr_3_IC"));
  H_t.push_back(CreateHistosT("corr_3_IC"));
  H_2D.push_back(CreateHistos2D("corr_3_IC"));

  H_d.push_back(CreateHistosD("uncorr_kT"));
  H_t.push_back(CreateHistosT("uncorr_kT"));
  H_2D.push_back(CreateHistos2D("uncorr_kT"));

  H_d.push_back(CreateHistosD("corr_0_kT"));
  H_t.push_back(CreateHistosT("corr_0_kT"));
  H_2D.push_back(CreateHistos2D("corr_0_kT"));
  
  H_d.push_back(CreateHistosD("corr_1_kT"));
  H_t.push_back(CreateHistosT("corr_1_kT"));
  H_2D.push_back(CreateHistos2D("corr_1_kT"));

  H_d.push_back(CreateHistosD("corr_2_kT"));
  H_t.push_back(CreateHistosT("corr_2_kT"));
  H_2D.push_back(CreateHistos2D("corr_2_kT"));

  H_d.push_back(CreateHistosD("corr_3_kT"));
  H_t.push_back(CreateHistosT("corr_3_kT"));
  H_2D.push_back(CreateHistos2D("corr_3_kT"));

  H_d.push_back(CreateHistosD("uncorr_SIS"));
  H_t.push_back(CreateHistosT("uncorr_SIS"));
  H_2D.push_back(CreateHistos2D("uncorr_SIS"));

  H_d.push_back(CreateHistosD("corr_0_SIS"));
  H_t.push_back(CreateHistosT("corr_0_SIS"));
  H_2D.push_back(CreateHistos2D("corr_0_SIS"));
  
  H_d.push_back(CreateHistosD("corr_1_SIS"));
  H_t.push_back(CreateHistosT("corr_1_SIS"));
  H_2D.push_back(CreateHistos2D("corr_1_SIS"));

  H_d.push_back(CreateHistosD("corr_2_SIS"));
  H_t.push_back(CreateHistosT("corr_2_SIS"));
  H_2D.push_back(CreateHistos2D("corr_2_SIS"));

  H_d.push_back(CreateHistosD("corr_3_SIS"));
  H_t.push_back(CreateHistosT("corr_3_SIS"));
  H_2D.push_back(CreateHistos2D("corr_3_SIS"));
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  nentries = 2000;
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    // check the alpgen ID
    
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

    ISO = false;

    maxpt = -1.0;
    int imax = -1;
    TVector3 vPV(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
    for(int i = 0; i < nTrack; i++){
      TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
      vT = vT-vPV;
      if(fabs(vT.z()) > 0.1) continue;
      if(fabs(vT.Mag()) > 1.0) continue;
      //if(abs(chargeTrack[i]) > 1) continue;
      double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
      if(pt < 1.2) continue;
      if(pt > 500.) continue;
      if(trackNormalizedChi2Track[i] > 20.0) continue;
      if(trackDxyPVTrack[i] > .6) continue;
      if(fabs(etaTrack[i]) > 2.4) continue;
      if(trackValidHitsTrack[i] < 5) continue;
      if(pt > maxpt){
	maxpt = pt;
	imax = i;
      }
    }
    if(maxpt < 15.0){
      ISO = false;
    } else {
      double sumpt = 0.0;
      for(int i = 0; i < nTrack; i++){
	if(i == imax) continue;
	TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
	vT = vT-vPV;
	if(fabs(vT.z()) > 0.1) continue;
	if(fabs(vT.Mag()) > 1.0) continue;
	//if(abs(chargeTrack[i]) > 1) continue;
	double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
	if(pt < 1.2) continue;
	if(pt > 500.) continue;
	if(trackNormalizedChi2Track[i] > 20.0) continue;
	if(trackDxyPVTrack[i] > .6) continue;
	if(fabs(etaTrack[i]) > 2.4) continue;
	if(trackValidHitsTrack[i] < 5) continue;
	double dR = DeltaR(etaTrack[i], phiTrack[i], etaTrack[imax], phiTrack[imax]);
	if(dR < 0.35){
	  sumpt += pt;
	}
      }
      
      if(sumpt/maxpt <= .1){
	ISO = true;
      }
    }


    //Create CaloTowers
    c_uncorr = CreateCaloTowers(thresh, 0.0, 0);
    c_uncorr_v = CreateCaloTowers(thresh, 0.0, 4);
    c_var = CreateCaloTowers(thresh, PVzPV[iPV], 3);
    c_fixed = CreateCaloTowers(thresh, PVzPV[iPV], 1);

    //cout << c_uncorr.size() << " " <<  c_uncorr_v.size() << " "; 
    //cout << c_fixed.size() << " " << c_var.size() << "here" << endl;

    //Create MET objects
    m_gen.push_back(MET(pxGenMet[0], pyGenMet[0], 0.0, 0.0));
    m_uncorr.push_back(CreateMET(c_uncorr));
    m_uncorr_v.push_back(CreateMET(c_uncorr_v));
    m_fixed.push_back(CreateMET(c_fixed));
    m_var.push_back(CreateMET(c_var));
    
    for(int ialg = 0; ialg < 3; ialg++){
      vector<Jet> j_u;
      vector<Jet> j_u2;
      
      j_u.clear();
      j_u2.clear();
      if(ialg == 0){
	j_u = CMSIterativeConeAlgorithm(c_uncorr, 0.5, 1.0);
      }
      if(ialg == 1){
	j_u = FastJetAlgorithm(c_uncorr, 0.6, 1.0);
      }
      if(ialg == 2){
	j_u = SISCone(c_uncorr, 0.5, 0.0);
      }
      j_uncorr_all = SortJet(j_u);
      for(int i = 0; i < j_uncorr_all.size();i++){
	j_u2.push_back(CorrectJet(j_uncorr_all[i], PVzPV[iPV]));
	if(j_uncorr_all[i].pt() > pt_cut && fabs(j_uncorr_all[i].eta()) < 3.0){
	  j_uncorr.push_back(j_uncorr_all[i]);
	}
      }
      j_corr_2_all = SortJet(j_u2);
      for(int i = 0; i < j_corr_2_all.size();i++){
	if(j_corr_2_all[i].pt() > pt_cut && fabs(j_corr_2_all[i].eta()) < 3.0){
	  j_corr_2.push_back(j_corr_2_all[i]);
	}
      }
      
      j_u.clear();
      j_u2.clear();
      
      if(ialg == 0){
	j_u = CMSIterativeConeAlgorithm(c_uncorr_v, 0.5, 1.0);
      }
      if(ialg == 1){
	j_u = FastJetAlgorithm(c_uncorr_v, 0.6, 1.0);
      }
      if(ialg == 2){
	j_u = SISCone(c_uncorr_v, 0.5, 0.0);
      }
      for(int i = 0; i < j_u.size(); i++){
	j_u2.push_back(CorrectJet(j_u[i], PVzPV[iPV]));
      }
      j_corr_3_all = SortJet(j_u2);
      for(int i = 0; i < j_corr_3_all.size(); i++){
	if(j_corr_3_all[i].pt() > pt_cut && fabs(j_corr_3_all[i].eta()) < 3.0){
	  j_corr_3.push_back(j_corr_3_all[i]);
	}
      }
      if(ialg == 0){
	j_corr_0_all = SortJet(CMSIterativeConeAlgorithm(c_fixed, 0.5, 1.0));
	j_corr_1_all = SortJet(CMSIterativeConeAlgorithm(c_var, 0.5, 1.0));
      }
      if(ialg == 1){
	j_corr_0_all = SortJet(FastJetAlgorithm(c_fixed, 0.6, 1.0));
	j_corr_1_all = SortJet(FastJetAlgorithm(c_var, 0.6, 1.0));
      }
      if(ialg == 2){
	j_corr_0_all = SortJet(SISCone(c_fixed, 0.5, 0.0));
	j_corr_1_all = SortJet(SISCone(c_var, 0.5, 0.0));
      }
      for(int i = 0; i < j_corr_0_all.size(); i++){
	if(j_corr_0_all[i].pt() > pt_cut && fabs(j_corr_0_all[i].eta()) < 3.0){
	  j_corr_0.push_back(j_corr_0_all[i]);
	}
      }
      for(int i = 0; i < j_corr_1_all.size(); i++){
	if(j_corr_1_all[i].pt() > pt_cut && fabs(j_corr_1_all[i].eta()) < 3.0){
	  j_corr_1.push_back(j_corr_1_all[i]);
	}
      }
      j_u.clear();
      j_u2.clear();
    
    
      j_dum = j_uncorr;
      j_dum_all = j_uncorr_all;
      EEMF = EventEMF(j_uncorr, 30.0, 3.0);
      ECHF = EventCHF(j_uncorr, 30.0, 1.7);
      m_dum = m_uncorr;
      FillHistos(H_d[0+ialg*5], H_t[0+ialg*5], H_2D[0+ialg*5]);
      
      
      j_dum  = j_corr_0;
      j_dum_all = j_corr_0_all;
      EEMF = EventEMF(j_corr_0, 30.0, 3.0);
      ECHF = EventCHF(j_corr_0, 30.0, 1.7);
      m_dum = m_fixed;
      FillHistos(H_d[1+ialg*5], H_t[1+ialg*5], H_2D[1+ialg*5]);
      
      
      j_dum  = j_corr_1;
      j_dum_all = j_corr_1_all;
      EEMF = EventEMF(j_corr_1, 30.0, 3.0);
      ECHF = EventCHF(j_corr_1, 30.0, 1.7);
      m_dum = m_var;
      FillHistos(H_d[2+ialg*5], H_t[2+ialg*5], H_2D[2+ialg*5]);
      
      j_dum  = j_corr_2;
      j_dum_all = j_corr_2_all;
      EEMF = EventEMF(j_corr_2, 30.0, 3.0);
      ECHF = EventCHF(j_corr_2, 30.0, 1.7);
      m_dum = m_uncorr;
      FillHistos(H_d[3+ialg*5], H_t[3+ialg*5], H_2D[3+ialg*5]);
      
      j_dum  = j_corr_3;
      j_dum_all = j_corr_3_all;
      EEMF = EventEMF(j_corr_3, 30.0, 3.0);
      ECHF = EventCHF(j_corr_3, 30.0, 1.7);
      m_dum = m_uncorr_v;
      FillHistos(H_d[4+ialg*5], H_t[4+ialg*5], H_2D[4+ialg*5]);

      j_uncorr.clear();
      j_corr_0.clear();
      j_corr_1.clear();
      j_corr_2.clear();
      j_corr_3.clear();
      j_uncorr_all.clear();
      j_corr_0_all.clear();
      j_corr_1_all.clear();
      j_corr_2_all.clear();
      j_corr_3_all.clear();
      
    }
    
    j_dum.clear();
    j_dum_all.clear();
    m_dum.clear();

    m_gen.clear();
    c_uncorr.clear();
    c_uncorr_v.clear();
    c_fixed.clear();
    c_var.clear();
    m_uncorr.clear();
    m_uncorr_v.clear();
    m_fixed.clear();
    m_var.clear();
  }
  
  TFile *file = new TFile("Histograms.root","RECREATE");
  WriteHistos(H_d[0],file,"uncorr_IC");
  WriteHistos(H_t[0],file,"uncorr_IC");
  WriteHistos(H_2D[0],file,"uncorr_IC");
  WriteHistos(H_d[1],file,"corr_0_IC");
  WriteHistos(H_t[1],file,"corr_0_IC");
  WriteHistos(H_2D[1],file,"corr_0_IC");
  WriteHistos(H_d[2],file,"corr_1_IC");
  WriteHistos(H_t[2],file,"corr_1_IC");
  WriteHistos(H_2D[2],file,"corr_1_IC");
  WriteHistos(H_d[3],file,"corr_2_IC");
  WriteHistos(H_t[3],file,"corr_2_IC");
  WriteHistos(H_2D[3],file,"corr_2_IC");
  WriteHistos(H_d[4],file,"corr_3_IC");
  WriteHistos(H_t[4],file,"corr_3_IC");
  WriteHistos(H_2D[4],file,"corr_3_IC");

  WriteHistos(H_d[5],file,"uncorr_kT");
  WriteHistos(H_t[5],file,"uncorr_kT");
  WriteHistos(H_2D[5],file,"uncorr_kT");
  WriteHistos(H_d[6],file,"corr_0_kT");
  WriteHistos(H_t[6],file,"corr_0_kT");
  WriteHistos(H_2D[6],file,"corr_0_kT");
  WriteHistos(H_d[7],file,"corr_1_kT");
  WriteHistos(H_t[7],file,"corr_1_kT");
  WriteHistos(H_2D[7],file,"corr_1_kT");
  WriteHistos(H_d[8],file,"corr_2_kT");
  WriteHistos(H_t[8],file,"corr_2_kT");
  WriteHistos(H_2D[8],file,"corr_2_kT");
  WriteHistos(H_d[9],file,"corr_3_kT");
  WriteHistos(H_t[9],file,"corr_3_kT");
  WriteHistos(H_2D[9],file,"corr_3_kT");

  WriteHistos(H_d[10],file,"uncorr_SIS");
  WriteHistos(H_t[10],file,"uncorr_SIS");
  WriteHistos(H_2D[10],file,"uncorr_SIS");
  WriteHistos(H_d[11],file,"corr_0_SIS");
  WriteHistos(H_t[11],file,"corr_0_SIS");
  WriteHistos(H_2D[11],file,"corr_0_SIS");
  WriteHistos(H_d[12],file,"corr_1_SIS");
  WriteHistos(H_t[12],file,"corr_1_SIS");
  WriteHistos(H_2D[12],file,"corr_1_SIS");
  WriteHistos(H_d[13],file,"corr_2_SIS");
  WriteHistos(H_t[13],file,"corr_2_SIS");
  WriteHistos(H_2D[13],file,"corr_2_SIS");
  WriteHistos(H_d[14],file,"corr_3_SIS");
  WriteHistos(H_t[14],file,"corr_3_SIS");
  WriteHistos(H_2D[14],file,"corr_3_SIS");

  // TFile *file = new TFile("Histograms.root","RECREATE");
//   //for(int i=0; i<41; i++) {
//   //  char name[32];
//   sprintf(name,"0_0_%i",0);//-20+i);
//   //  WriteHistos(Histos[i], file,name);
//   WriteHistos(Histos[0], file,name);
//   //} 
}
