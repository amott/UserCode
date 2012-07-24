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

#include "Vecbos.hh"
#include "EventShapeAnalysis.hh"

EventShapeAnalysis::EventShapeAnalysis(TTree *tree) : Vecbos(tree) {}

EventShapeAnalysis::~EventShapeAnalysis(){}

vector<TH2D*> EventShapeAnalysis::CreateHistos2D(string dirname){
  vector<TH2D*> h_2D;

  string name;

  string type[4];
  type[0] = "_tj";
  type[1] = "_cj";
  type[2] = "_ct";
  type[3] = "_aj";

  for(int i = 0; i < 4; i++){
    name = dirname+type[i]+"_Stheta_v_Ttheta"; //0
    h_2D.push_back(new TH2D(name.c_str(),name.c_str(), 200, 0.0, TMath::Pi()/2, 200, 0.0, TMath::Pi()/2));
    name = dirname+type[i]+"_Stheta_v_TMtheta"; //1
    h_2D.push_back(new TH2D(name.c_str(),name.c_str(), 200, 0.0, TMath::Pi()/2, 200, 0.0, TMath::Pi()/2));
    name = dirname+type[i]+"_Ttheta_v_TMtheta"; //2
    h_2D.push_back(new TH2D(name.c_str(),name.c_str(), 200, 0.0, TMath::Pi()/2, 200, 0.0, TMath::Pi()/2));
  }

  return h_2D;
}

vector<TH1D*> EventShapeAnalysis::CreateHistosD(string dirname){
  vector<TH1D*> h_d;

  string name;

  string type[4];
  type[0] = "_tj";
  type[1] = "_cj";
  type[2] = "_ct";
  type[3] = "_aj";

  for(int i = 0; i < 4; i++){
    name = dirname+type[i]+"_FW_0"; //0
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_FW_1"; //1
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_FW_2"; //2
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_FW_3"; //3
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_FW_4"; //4
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_FW_5"; //5
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_FW_6"; //6
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_TFW_1"; //7
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_TFW_2"; //8
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_TFW_3"; //9
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_TFW_4"; //10
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_TFW_5"; //11
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    name = dirname+type[i]+"_TFW_6"; //12
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_Aplan"; //13
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 0.5));
    
    name = dirname+type[i]+"_Spher"; //14
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_detTS"; //15
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_TranSpher"; //16
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_Thrust"; //17
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.5, 1.0));
    
    name = dirname+type[i]+"_TranThrust"; //18
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.5, 1.0));
    
    name = dirname+type[i]+"_Oblat"; //19
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_ThrustMaj"; //20
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_ThrustMin"; //21
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_TranThrustMin"; //22
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_Spher_theta"; //23
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, TMath::Pi()));
    
    name = dirname+type[i]+"_Thrust_theta"; //24
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, TMath::Pi()));
    
    name = dirname+type[i]+"_ThrustMaj_theta"; //25
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, TMath::Pi()));
    
    name = dirname+type[i]+"_ThrustMin_theta"; //26
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, TMath::Pi()));
    
    name = dirname+type[i]+"_SdotT"; //27
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_TSdotTT"; //28
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_SdotTM"; //29
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_TSdotMHT"; //30
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
    
    name = dirname+type[i]+"_TTdotMHT"; //31
    h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
  }
  
  string axes[7];
  axes[0] = "S";
  axes[1] = "TS";
  axes[2] = "T";
  axes[3] = "TT";
  axes[4] = "TM";
  axes[5] = "TMi";
  axes[6] = "MHT";

  string part[4];
  part[0] = "lead_mu";
  part[1] = "lead_j";
  part[2] = "lead_tr";
  part[3] = "Z";

  for(int i = 0; i < 2; i++){
    for(int iaxes = 0; iaxes < 7; iaxes++){
      for(int ipart = 0; ipart < 4; ipart++){
	name = dirname+type[2+i]+"_"+axes[iaxes]+"dot"+part[ipart];
	h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, 1.0));
	name = dirname+type[2+i]+"_"+axes[iaxes]+"phi"+part[ipart];
	h_d.push_back(new TH1D(name.c_str(),name.c_str(), 500, 0.0, TMath::Pi()));
      }
    }
  }


  return h_d;
}

vector<TProfile*> EventShapeAnalysis::CreateHistosT(string dirname){
  vector<TProfile*> h_p;

  return h_p;
}

void EventShapeAnalysis::FillHistos(vector<TH1D*> h_d, vector<TProfile*> h_p, vector<TH2D*> h_2D){

  TVector3 axis[7];
  
  float genWeight = 1.0;

  for(int i = 0; i < 4; i++){
    axis[0] = CT[i]->SphericityAxis();
    axis[1] = CT[i]->TSphericityAxis();
    axis[2] = CT[i]->ThrustAxis();
    axis[3] = CT[i]->TranThrustAxis();
    axis[4] = CT[i]->ThrustMajorAxis();
    axis[5] = CT[i]->ThrustMinorAxis();
    axis[6] = MHT[i].Vect();

    h_2D[0+i*3]->Fill(fabs(axis[0].Theta()-TMath::Pi()/2), fabs(axis[2].Theta()-TMath::Pi()/2),genWeight);
    h_2D[1+i*3]->Fill(fabs(axis[0].Theta()-TMath::Pi()/2), fabs(axis[4].Theta()-TMath::Pi()/2),genWeight);
    h_2D[2+i*3]->Fill(fabs(axis[2].Theta()-TMath::Pi()/2), fabs(axis[4].Theta()-TMath::Pi()/2),genWeight);

    for(int j = 0; j < 7; j++){
      if(j > 0){
	h_d[j+i*32]->Fill(CT[i]->FoxWolfram(JETS[i], j)/CT[i]->FoxWolfram(JETS[i], 0),genWeight);
	h_d[j+6+i*32]->Fill(CT[i]->TranFoxWolfram(JETS[i], j),genWeight);
      } else {
	h_d[j+i*32]->Fill(CT[i]->FoxWolfram(JETS[i], j),genWeight);
      }
    }
    h_d[13+i*32]->Fill(CT[i]->Aplanarity(),genWeight);
    h_d[14+i*32]->Fill(CT[i]->Sphericity(),genWeight);
    h_d[15+i*32]->Fill(CT[i]->detTS(),genWeight);
    h_d[16+i*32]->Fill(CT[i]->TranSphericity(),genWeight);
    h_d[17+i*32]->Fill(CT[i]->Thrust(),genWeight);
    h_d[18+i*32]->Fill(CT[i]->TranThrust(),genWeight);
    h_d[19+i*32]->Fill(CT[i]->Oblateness(),genWeight);
    h_d[20+i*32]->Fill(CT[i]->ThrustMajor(),genWeight);
    h_d[21+i*32]->Fill(CT[i]->ThrustMinor(),genWeight);
    h_d[22+i*32]->Fill(CT[i]->TranThrustMinor(),genWeight);
    h_d[23+i*32]->Fill(axis[0].Theta(),genWeight);
    h_d[24+i*32]->Fill(axis[2].Theta(),genWeight);
    h_d[25+i*32]->Fill(axis[4].Theta(),genWeight);
    h_d[26+i*32]->Fill(axis[5].Theta(),genWeight);
    
    h_d[27+i*32]->Fill(fabs(axis[0].Dot(axis[2])/(axis[0].Mag()*axis[2].Mag())),genWeight);
    h_d[28+i*32]->Fill(fabs(axis[1].Dot(axis[3])/(axis[1].Mag()*axis[3].Mag())),genWeight);
    h_d[29+i*32]->Fill(fabs(axis[0].Dot(axis[5])/(axis[0].Mag()*axis[5].Mag())),genWeight);
    h_d[30+i*32]->Fill(fabs(axis[1].Dot(axis[6])/(axis[1].Mag()*axis[6].Mag())),genWeight);
    h_d[31+i*32]->Fill(fabs(axis[3].Dot(axis[6])/(axis[3].Mag()*axis[6].Mag())),genWeight);
  }

  for(int i = 0; i < 2; i++){
    TVector3 axes[7];
    axes[0] = CT[i+2]->SphericityAxis();
    axes[1] = CT[i+2]->TSphericityAxis();
    axes[2] = CT[i+2]->ThrustAxis();
    axes[3] = CT[i+2]->TranThrustAxis();
    axes[4] = CT[i+2]->ThrustMajorAxis();
    axes[5] = CT[i+2]->ThrustMinorAxis();
    axes[6] = MHT[i+2].Vect();
    
    TVector3 part[4];
    part[0] = lead_mu.Vect();
    part[1] = lead_j.Vect();
    part[2] = lead_track.Vect();
    part[3] = Z.Vect();

    for(int iaxes = 0; iaxes < 7; iaxes++){
      for(int ipart = 0; ipart < 4; ipart++){
	h_d[128+i*56+iaxes*8+ipart*2]->Fill(fabs(axes[iaxes].Dot(part[ipart])/
						 (axes[iaxes].Mag()*part[ipart].Mag())),genWeight);
	if(iaxes == 6){
	  h_d[128+i*56+iaxes*8+ipart*2+1]->Fill(fabs(DeltaPhi(axes[iaxes].Phi(),part[ipart].Phi())),genWeight);
	} else {
	  double dphi = axes[iaxes].Phi();
	  if(fabs(DeltaPhi(axes[iaxes].Phi(),part[ipart].Phi())) > TMath::Pi()/2){
	    if(fabs(dphi+TMath::Pi()) > fabs(dphi-TMath::Pi())){
	      dphi -= TMath::Pi();
	    } else {
	      dphi += TMath::Pi();
	    }
	  }
	  h_d[128+i*56+iaxes*8+ipart*2+1]->Fill(fabs(DeltaPhi(dphi,part[ipart].Phi())),genWeight);
	}
      }
    }
  }
}
void EventShapeAnalysis::Loop(string outname) {
  outfilename = outname;
  if(fChain == 0) return;

  vector< vector<TH1D*> > H_d;
  vector< vector<TProfile*> > H_t;
  vector< vector<TH2D*> > H_2D;

  H_d.push_back(CreateHistosD("0calojet"));
  H_t.push_back(CreateHistosT("0calojet"));
  H_2D.push_back(CreateHistos2D("0calojet"));
  H_d.push_back(CreateHistosD("1calojet"));
  H_t.push_back(CreateHistosT("1calojet"));
  H_2D.push_back(CreateHistos2D("1calojet"));
  H_d.push_back(CreateHistosD("2calojet"));
  H_t.push_back(CreateHistosT("2calojet"));
  H_2D.push_back(CreateHistos2D("2calojet"));
  H_d.push_back(CreateHistosD("3calojet"));
  H_t.push_back(CreateHistosT("3calojet"));
  H_2D.push_back(CreateHistos2D("3calojet"));
  H_d.push_back(CreateHistosD("4calojet"));
  H_t.push_back(CreateHistosT("4calojet"));
  H_2D.push_back(CreateHistos2D("4calojet"));
  H_d.push_back(CreateHistosD("5calojet"));
  H_t.push_back(CreateHistosT("5calojet"));
  H_2D.push_back(CreateHistos2D("5calojet"));
  H_d.push_back(CreateHistosD("0trackjet"));
  H_t.push_back(CreateHistosT("0trackjet"));
  H_2D.push_back(CreateHistos2D("0trackjet"));
  H_d.push_back(CreateHistosD("1trackjet"));
  H_t.push_back(CreateHistosT("1trackjet"));
  H_2D.push_back(CreateHistos2D("1trackjet"));
  H_d.push_back(CreateHistosD("2trackjet"));
  H_t.push_back(CreateHistosT("2trackjet"));
  H_2D.push_back(CreateHistos2D("2trackjet"));
  H_d.push_back(CreateHistosD("3trackjet"));
  H_t.push_back(CreateHistosT("3trackjet"));
  H_2D.push_back(CreateHistos2D("3trackjet"));
  H_d.push_back(CreateHistosD("4trackjet"));
  H_t.push_back(CreateHistosT("4trackjet"));
  H_2D.push_back(CreateHistos2D("4trackjet"));
  H_d.push_back(CreateHistosD("5trackjet"));
  H_t.push_back(CreateHistosT("5trackjet"));
  H_2D.push_back(CreateHistos2D("5trackjet"));
  H_d.push_back(CreateHistosD("Wselec_0calojet"));
  H_t.push_back(CreateHistosT("Wselec_0calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_0calojet"));
  H_d.push_back(CreateHistosD("Wselec_1calojet"));
  H_t.push_back(CreateHistosT("Wselec_1calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_1calojet"));
  H_d.push_back(CreateHistosD("Wselec_2calojet"));
  H_t.push_back(CreateHistosT("Wselec_2calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_2calojet"));
  H_d.push_back(CreateHistosD("Wselec_3calojet"));
  H_t.push_back(CreateHistosT("Wselec_3calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_3calojet"));
  H_d.push_back(CreateHistosD("Wselec_4calojet"));
  H_t.push_back(CreateHistosT("Wselec_4calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_4calojet"));
  H_d.push_back(CreateHistosD("Wselec_5calojet"));
  H_t.push_back(CreateHistosT("Wselec_5calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_5calojet"));
  H_d.push_back(CreateHistosD("Wselec_0trackjet"));
  H_t.push_back(CreateHistosT("Wselec_0trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_0trackjet"));
  H_d.push_back(CreateHistosD("Wselec_1trackjet"));
  H_t.push_back(CreateHistosT("Wselec_1trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_1trackjet"));
  H_d.push_back(CreateHistosD("Wselec_2trackjet"));
  H_t.push_back(CreateHistosT("Wselec_2trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_2trackjet"));
  H_d.push_back(CreateHistosD("Wselec_3trackjet"));
  H_t.push_back(CreateHistosT("Wselec_3trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_3trackjet"));
  H_d.push_back(CreateHistosD("Wselec_4trackjet"));
  H_t.push_back(CreateHistosT("Wselec_4trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_4trackjet"));
  H_d.push_back(CreateHistosD("Wselec_5trackjet"));
  H_t.push_back(CreateHistosT("Wselec_5trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_5trackjet"));
  H_d.push_back(CreateHistosD("Wselec_boost_0calojet"));
  H_t.push_back(CreateHistosT("Wselec_boost_0calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_0calojet"));
  H_d.push_back(CreateHistosD("Wselec_boost_1calojet"));
  H_t.push_back(CreateHistosT("Wselec_boost_1calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_1calojet"));
  H_d.push_back(CreateHistosD("Wselec_boost_2calojet"));
  H_t.push_back(CreateHistosT("Wselec_boost_2calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_2calojet"));
  H_d.push_back(CreateHistosD("Wselec_boost_3calojet"));
  H_t.push_back(CreateHistosT("Wselec_boost_3calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_3calojet"));
  H_d.push_back(CreateHistosD("Wselec_boost_4calojet"));
  H_t.push_back(CreateHistosT("Wselec_boost_4calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_4calojet"));
  H_d.push_back(CreateHistosD("Wselec_boost_5calojet"));
  H_t.push_back(CreateHistosT("Wselec_boost_5calojet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_5calojet"));
  H_d.push_back(CreateHistosD("Wselec_boost_0trackjet"));
  H_t.push_back(CreateHistosT("Wselec_boost_0trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_0trackjet"));
  H_d.push_back(CreateHistosD("Wselec_boost_1trackjet"));
  H_t.push_back(CreateHistosT("Wselec_boost_1trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_1trackjet"));
  H_d.push_back(CreateHistosD("Wselec_boost_2trackjet"));
  H_t.push_back(CreateHistosT("Wselec_boost_2trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_2trackjet"));
  H_d.push_back(CreateHistosD("Wselec_boost_3trackjet"));
  H_t.push_back(CreateHistosT("Wselec_boost_3trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_3trackjet"));
  H_d.push_back(CreateHistosD("Wselec_boost_4trackjet"));
  H_t.push_back(CreateHistosT("Wselec_boost_4trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_4trackjet"));
  H_d.push_back(CreateHistosD("Wselec_boost_5trackjet"));
  H_t.push_back(CreateHistosT("Wselec_boost_5trackjet"));
  H_2D.push_back(CreateHistos2D("Wselec_boost_5trackjet"));
  H_d.push_back(CreateHistosD("Zselec_0calojet"));
  H_t.push_back(CreateHistosT("Zselec_0calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_0calojet"));
  H_d.push_back(CreateHistosD("Zselec_1calojet"));
  H_t.push_back(CreateHistosT("Zselec_1calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_1calojet"));
  H_d.push_back(CreateHistosD("Zselec_2calojet"));
  H_t.push_back(CreateHistosT("Zselec_2calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_2calojet"));
  H_d.push_back(CreateHistosD("Zselec_3calojet"));
  H_t.push_back(CreateHistosT("Zselec_3calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_3calojet"));
  H_d.push_back(CreateHistosD("Zselec_4calojet"));
  H_t.push_back(CreateHistosT("Zselec_4calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_4calojet"));
  H_d.push_back(CreateHistosD("Zselec_5calojet"));
  H_t.push_back(CreateHistosT("Zselec_5calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_5calojet"));
  H_d.push_back(CreateHistosD("Zselec_0trackjet"));
  H_t.push_back(CreateHistosT("Zselec_0trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_0trackjet"));
  H_d.push_back(CreateHistosD("Zselec_1trackjet"));
  H_t.push_back(CreateHistosT("Zselec_1trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_1trackjet"));
  H_d.push_back(CreateHistosD("Zselec_2trackjet"));
  H_t.push_back(CreateHistosT("Zselec_2trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_2trackjet"));
  H_d.push_back(CreateHistosD("Zselec_3trackjet"));
  H_t.push_back(CreateHistosT("Zselec_3trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_3trackjet"));
  H_d.push_back(CreateHistosD("Zselec_4trackjet"));
  H_t.push_back(CreateHistosT("Zselec_4trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_4trackjet"));
  H_d.push_back(CreateHistosD("Zselec_5trackjet"));
  H_t.push_back(CreateHistosT("Zselec_5trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_5trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_0calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_0calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_0calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_1calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_1calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_1calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_2calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_2calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_2calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_3calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_3calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_3calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_4calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_4calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_4calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_5calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_5calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_5calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_0trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_0trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_0trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_1trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_1trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_1trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_2trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_2trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_2trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_3trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_3trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_3trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_4trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_4trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_4trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_Z_5trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_Z_5trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_Z_5trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_0calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_0calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_0calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_1calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_1calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_1calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_2calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_2calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_2calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_3calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_3calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_3calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_4calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_4calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_4calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_5calojet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_5calojet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_5calojet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_0trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_0trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_0trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_1trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_1trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_1trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_2trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_2trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_2trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_3trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_3trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_3trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_4trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_4trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_4trackjet"));
  H_d.push_back(CreateHistosD("Zselec_boost_M_5trackjet"));
  H_t.push_back(CreateHistosT("Zselec_boost_M_5trackjet"));
  H_2D.push_back(CreateHistos2D("Zselec_boost_M_5trackjet"));
  
  

  for(int i = 0; i < 4; i++){
    CT[i] = new CoolTools();
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
    // check the alpgen ID
    
    // RecHit thresholds for CaloTowers
    vector<float> thresh;
    thresh.push_back(.9);       //HB
    thresh.push_back(1.1);      //etc...
    thresh.push_back(1.4);   
    thresh.push_back(1.4);
    thresh.push_back(1.2);
    thresh.push_back(1.8);
    thresh.push_back(.09);
    thresh.push_back(.45);
    thresh.push_back(.2);
    thresh.push_back(.45);

    CaloF = 0.06;
    double calo_pt_cut = 30;
    double track_pt_cut = 15;
    if(nPV < 1) continue;

    int iPV = -1;
    double maxpt = -1;
    for(int i = 0; i < nPV; i++){
      if(SumPtPV[i] > maxpt){
	maxpt = SumPtPV[i];
	iPV = i;
      }
    }

    x0 = PVxPV[iPV]; 
    y0 = PVyPV[iPV]; 
    z0 = PVzPV[iPV];
    
    // 3 := making CaloTowers with depth from 'CaloF' depth (%)
    c_calo = CreateCaloTowers(thresh, z0, 3);
   
    calo_jet.clear();
    calo_jet = SortJet(CMSIterativeConeAlgorithm(c_calo, 0.5, 1.0));

    N_calo_jet = 0;
    for(int i = 0; i < calo_jet.size(); i++){
      if(calo_jet[i].pt() > calo_pt_cut && fabs(calo_jet[i].eta()) <= 3.0)
	N_calo_jet++; 
    }

    if(calo_jet.size() == 0) continue;

    // Let's do EWK selection before clustering the track jet collection so we 
    // know the 'expected' muon content and can remove them

    int mu_track[2];
    mu_track[0] = -1;
    mu_track[1] = -1;

    // W selec
    W_pass = false;
    int w_mu_index = -1;
    int mu_count = 0;
    for(int imu = 0; imu < nMuon; imu++){
      double mu_pt = sqrt(pxMuon[imu]*pxMuon[imu]+pyMuon[imu]*pyMuon[imu]);
      if(mu_pt > 15 && fabs(etaMuon[imu]) < 2.4){
	mu_count++;
	w_mu_index = imu;
      }
    }
    if(mu_count == 1){
      if(sumPt05Muon[w_mu_index] < 10 &&
	 fabs(vertexZMuon[w_mu_index]-z0) < 1.0 &&
	 fabs(muTrackDxyMuon[w_mu_index]/muTrackDxyErrorMuon[w_mu_index]) < 3.0){
	// Got our W muon - find the track associated with it to remove from 
	// the track jet clustering

	W_pass = true;
	double dRmin_mu = 99999999999.;
	for(int i=0; i< nTrack; i++) {
	  double dR = DeltaR(etaMuon[w_mu_index],phiMuon[w_mu_index],etaTrack[i],phiTrack[i]);
	  if(dR<dRmin_mu) {
	    mu_track[0] = i;
	    dRmin_mu = dR;
	  }
	}

	if(dRmin_mu > 0.01) mu_track[0] = -1;
      }
    }

   

    // Now for the Z selection
    Z_pass = false;
    z_index = -1;
    int z_mu_index[2];
    double z_mu_pt[2];
    z_mu_pt[0] = 0.0;
    z_mu_pt[1] = 0.0;
    if(W_pass == false){
      if(nZ0ToMuMuVtx != 0){
	//Get best Z-candidate
	vector<double> masses;
	for(int iZ = 0; iZ < nZ0ToMuMuVtx; iZ++){
	  masses.push_back(massZ0ToMuMuVtx[iZ]);
	}
	double bestd = 9999999999.;
	for(int i=0; i<masses.size(); i++) {
	  double thisd = fabs(masses[i]-Ztruemass);
	  if(thisd < bestd) {
	    bestd = thisd;
	    z_index = i;
	  }
	}
	// Check that it's within mass window
	if(bestd < 15){
	  // Check muons for vertex compatibility and pt/eta requiremnt
	  z_mu_index[0] = d1IndexZ0ToMuMuVtx[z_index];
	  z_mu_index[1] = d2IndexZ0ToMuMuVtx[z_index];
	  if(sumPt05Muon[z_mu_index[0]] < 10 &&
	     fabs(vertexZMuon[z_mu_index[0]]-z0) < 1.0 &&
	     fabs(muTrackDxyMuon[z_mu_index[0]]/muTrackDxyErrorMuon[z_mu_index[0]]) < 3.0 &&
	     sumPt05Muon[z_mu_index[1]] < 10 &&
	     fabs(vertexZMuon[z_mu_index[1]]-z0) < 1.0 &&
	     fabs(muTrackDxyMuon[z_mu_index[1]]/muTrackDxyErrorMuon[z_mu_index[1]]) < 3.0){
	    z_mu_pt[0] = sqrt(pxMuon[z_mu_index[0]]*pxMuon[z_mu_index[0]]+
			      pyMuon[z_mu_index[0]]*pyMuon[z_mu_index[0]]);
	    z_mu_pt[1] = sqrt(pxMuon[z_mu_index[1]]*pxMuon[z_mu_index[1]]+
			      pyMuon[z_mu_index[1]]*pyMuon[z_mu_index[1]]);
	    if(z_mu_pt[0] > 15 && z_mu_pt[1] > 15 && 
	       fabs(etaMuon[z_mu_index[0]]) < 2.4 && fabs(etaMuon[z_mu_index[0]])){
	      Z_pass = true;
	      for(int i = 0; i < nMuon; i++){
		if(i == z_mu_index[0] || i == z_mu_index[1]) continue;
		double mu_pt = sqrt(pxMuon[i]*pxMuon[i]+
				    pyMuon[i]*pyMuon[i]);
		if(mu_pt > z_mu_pt[0] || mu_pt > z_mu_pt[1]){
		  Z_pass = false;
		  break;
		}
	      }
	    }
	  }
	}
      }
    }

    if(Z_pass){
      for(int j = 0; j < 2; j++){
	double dRmin_mu = 9999999;
	for(int i=0; i< nTrack; i++) {
	  double dR = DeltaR(etaMuon[z_mu_index[j]],phiMuon[z_mu_index[j]],
			     etaTrack[i],phiTrack[i]);
	  if(dR < dRmin_mu) {
	    mu_track[j] = i;
	    dRmin_mu = dR;
	  }
	}
	
	if(dRmin_mu > 0.01) mu_track[j] = -1;
      }
    }
    lead_mu_pt = 0.0;
    for(int i = 0; i < nMuon; i++){
      double mu_pt = sqrt(pxMuon[i]*pxMuon[i]+
			  pyMuon[i]*pyMuon[i]);
      if(mu_pt > lead_mu_pt){
	lead_mu_pt = mu_pt;
	lead_mu_index = i;
      }
    }

    // make track jet collection - remove muons from W or Z
    // also, get lead track
    double track_pt_max = -1.0;
    int i_track_max = -1;
    
    track_collection.clear();
    TVector3 vPV(x0, y0, z0);
    for(int i = 0; i < nTrack; i++){
      if(i == mu_track[0] || i == mu_track[1]) continue;
      TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
      vT = vT-vPV;
      if(fabs(vT.z()) > 0.1) continue;
      if(fabs(vT.Mag()) > 1.0) continue;
      double pt = sqrt(pxTrack[i]*pxTrack[i]+pyTrack[i]*pyTrack[i]);
      if(pt < 0.5) continue;
      if(pt > 500.0) continue;
      if(pt > track_pt_max){
	track_pt_max = pt;
	i_track_max = i;
      }
      if(trackNormalizedChi2Track[i] > 20.0) continue;
      if(fabs(trackDxyPVTrack[i]) > 0.6) continue;
      if(fabs(trackDxyPVTrack[i]/trackDxyErrorTrack[i]) > 6) continue;
      if(fabs(etaTrack[i]) > 2.4) continue;
      if(trackValidHitsTrack[i] < 5) continue;
      TVector3 v;
      v.SetPtEtaPhi(pt, etaTrack[i], phiTrack[i]);
      track_collection.push_back(CaloTower(pt*cosh(etaTrack[i]), 0.0, v, v, v));
    }
    //cluster track jets (this without the muons id'ed for the EWK sel.)
    track_jet.clear();
    track_jet = SortJet(FastJetAlgorithm(track_collection, 0.4, 1.0));

    N_track_jet = 0;
    for(int i = 0; i < track_jet.size(); i++){
      if(track_jet[i].pt() > track_pt_cut && fabs(track_jet[i].eta()) <= 3.0)
	N_track_jet++; 
    }

    if(i_track_max >= 0.0){
      lead_track.SetPtEtaPhiE(track_pt_max, etaTrack[i_track_max],
			      phiTrack[i_track_max], cosh(etaTrack[i_track_max])*etTrack[i_track_max]);
    } else {
      lead_track.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
    }

    if(Z_pass){
      Z.SetPxPyPzE(pxZ0ToMuMuVtx[z_index], pyZ0ToMuMuVtx[z_index], 
		   pzZ0ToMuMuVtx[z_index], energyZ0ToMuMuVtx[z_index]);
    } else {
      Z.SetPxPyPzE(0.0,0.0,0.0,0.0);
    }

    lead_j.SetPxPyPzE(calo_jet[0].px(),calo_jet[0].py(),
		      calo_jet[0].pz(),calo_jet[0].e());

    if(lead_mu_index >= 0){
      lead_mu.SetPxPyPzE(pxMuon[lead_mu_index],pyMuon[lead_mu_index],
			 pzMuon[lead_mu_index],energyMuon[lead_mu_index]);
    } else {
      lead_mu.SetPxPyPzE(0.0,0.0,0.0,0.0);
    }


    //Now let's get down to business......
  

  

    MHT[0] = CreateMET(track_jet).metVector();
    MHT[1] = CreateMET(calo_jet).metVector();
    MHT[2] = CreateMET(c_calo).metVector();
    MHT[3] = CreateMET(track_collection).metVector();

    JETS.clear();

    JETS.push_back(CT[0]->Get4Vectors(track_jet));
    JETS.push_back(CT[0]->Get4Vectors(calo_jet));
    JETS.push_back(CT[0]->Get4Vectors(CT[0]->CaloTowers2Jets(c_calo, 0)));
    JETS.push_back(CT[0]->Get4Vectors(CT[0]->CaloTowers2Jets(track_collection, 0)));

    for(int i = 0; i < 4; i++){
      if(JETS[i].size() == 0) continue;
      CT[i]->CalcSphericity(JETS[i]);
      CT[i]->CalcTSphericity(JETS[i]);
      CT[i]->CalcThrust(JETS[i]);
      CT[i]->CalcTranThrust(JETS[i]);//Now let's get down to business......
    } 

    for(int ijet = 0; ijet < 6; ijet++){
      if(N_calo_jet >= ijet){
	FillHistos(H_d[ijet], H_t[ijet], 
		   H_2D[ijet]);
      }
      if(N_track_jet >= ijet){
	FillHistos(H_d[ijet+6], H_t[ijet+6], 
		   H_2D[ijet+6]);
      }
    }

   

    if(W_pass || Z_pass){
      int i_index = 3;
      if(W_pass) i_index = 1;
      for(int ijet = 0; ijet < 6; ijet++){
	if(N_calo_jet >= ijet){
	  FillHistos(H_d[i_index*12+ijet], 
		     H_t[i_index*12+ijet], 
		     H_2D[i_index*12+ijet]);
	}
	if(N_track_jet >= ijet){
	  FillHistos(H_d[i_index*12+ijet+6], 
		     H_t[i_index*12+ijet+6], 
		     H_2D[i_index*12+ijet+6]);
	}
      }
    }
    // Use MET to boost jets
    if(W_pass){
      TLorentzVector W_from_MET;
      double E = sqrt(Wtruemass*Wtruemass+MHT[2].E()*MHT[2].E());
      W_from_MET.SetPxPyPzE(MHT[2].Px(), MHT[2].Py(), 0.0, E);
      TVector3 b = W_from_MET.BoostVector();
      JETS.clear();

      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(track_jet, b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(calo_jet, b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(CT[0]->CaloTowers2Jets(c_calo, 0), b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(CT[0]->CaloTowers2Jets(track_collection, 0), b)));

      for(int i = 0; i < 4; i++){
	if(JETS[i].size() == 0) continue;
	CT[i]->CalcSphericity(JETS[i]);
	CT[i]->CalcTSphericity(JETS[i]);
	CT[i]->CalcThrust(JETS[i]);
	CT[i]->CalcTranThrust(JETS[i]);
      }

      for(int ijet = 0; ijet < 6; ijet++){
	if(N_calo_jet >= ijet){
	  FillHistos(H_d[24+ijet], 
		     H_t[24+ijet], 
		     H_2D[24+ijet]);
	}
	if(N_track_jet >= ijet){
	  FillHistos(H_d[24+ijet+6], 
		     H_t[24+ijet+6], 
		     H_2D[24+ijet+6]);
	}
      }
      
    }
    if(Z_pass){
      TLorentzVector Z_from_MET;
      double E = sqrt(Ztruemass*Ztruemass+MHT[2].E()*MHT[2].E());
      Z_from_MET.SetPxPyPzE(MHT[2].Px(), MHT[2].Py(), 0.0, E);
      TVector3 b = Z_from_MET.BoostVector();
      JETS.clear();

      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(track_jet, b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(calo_jet, b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(CT[0]->CaloTowers2Jets(c_calo, 0), b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(CT[0]->CaloTowers2Jets(track_collection, 0), b)));

      for(int i = 0; i < 4; i++){
	if(JETS[i].size() == 0) continue;
	CT[i]->CalcSphericity(JETS[i]);
	CT[i]->CalcTSphericity(JETS[i]);
	CT[i]->CalcThrust(JETS[i]);
	CT[i]->CalcTranThrust(JETS[i]);
      }

      for(int ijet = 0; ijet < 6; ijet++){
	if(N_calo_jet >= ijet){
	  FillHistos(H_d[60+ijet], 
		     H_t[60+ijet], 
		     H_2D[60+ijet]);
	}
	if(N_track_jet >= ijet){
	  FillHistos(H_d[60+ijet+6], 
		     H_t[60+ijet+6], 
		     H_2D[60+ijet+6]);
	}
      }

      TLorentzVector Z_from_mu;
      E = sqrt(pxZ0ToMuMuVtx[z_index]*pxZ0ToMuMuVtx[z_index]+pyZ0ToMuMuVtx[z_index]*pyZ0ToMuMuVtx[z_index]+
	       pzZ0ToMuMuVtx[z_index]*pzZ0ToMuMuVtx[z_index]+Ztruemass*Ztruemass);
      Z_from_mu.SetPxPyPzE(pxZ0ToMuMuVtx[z_index], pyZ0ToMuMuVtx[z_index], pzZ0ToMuMuVtx[z_index], E);
      b = Z_from_MET.BoostVector();
      JETS.clear();

      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(track_jet, b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(calo_jet, b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(CT[0]->CaloTowers2Jets(c_calo, 0), b)));
      JETS.push_back(CT[0]->Get4Vectors(CT[0]->BoostJets(CT[0]->CaloTowers2Jets(track_collection, 0), b)));

      for(int i = 0; i < 4; i++){
	if(JETS[i].size() == 0) continue;
	CT[i]->CalcSphericity(JETS[i]);
	CT[i]->CalcTSphericity(JETS[i]);
	CT[i]->CalcThrust(JETS[i]);
	CT[i]->CalcTranThrust(JETS[i]);
      }

      for(int ijet = 0; ijet < 6; ijet++){
	if(N_calo_jet >= ijet){
	  FillHistos(H_d[48+ijet], 
		     H_t[48+ijet], 
		     H_2D[48+ijet]);
	}
	if(N_track_jet >= ijet){
	  FillHistos(H_d[48+ijet+6], 
		     H_t[48+ijet+6], 
		     H_2D[48+ijet+6]);
	}
      }
    }
    

    JETS.clear();
    
    c_calo.clear();
    track_collection.clear();
    c_all.clear();
    c_all_boost.clear();

    calo_jet.clear();
    track_jet.clear();
    boost_calo_jet.clear();
    boost_track_jet.clear();
    j_dum.clear();

  }
  int index_d = -1;
  int index_t = -1;
  int index_2D = -1;

    
  TFile *file = new TFile(outfilename.c_str(),"RECREATE");
  WriteHistos(H_d[++index_d],file,"0calojet");
  WriteHistos(H_t[++index_t],file,"0calojet");
  WriteHistos(H_2D[++index_2D],file,"0calojet");
  WriteHistos(H_d[++index_d],file,"1calojet");
  WriteHistos(H_t[++index_t],file,"1calojet");
  WriteHistos(H_2D[++index_2D],file,"1calojet");
  WriteHistos(H_d[++index_d],file,"2calojet");
  WriteHistos(H_t[++index_t],file,"2calojet");
  WriteHistos(H_2D[++index_2D],file,"2calojet");
  WriteHistos(H_d[++index_d],file,"3calojet");
  WriteHistos(H_t[++index_t],file,"3calojet");
  WriteHistos(H_2D[++index_2D],file,"3calojet");
  WriteHistos(H_d[++index_d],file,"4calojet");
  WriteHistos(H_t[++index_t],file,"4calojet");
  WriteHistos(H_2D[++index_2D],file,"4calojet");
  WriteHistos(H_d[++index_d],file,"5calojet");
  WriteHistos(H_t[++index_t],file,"5calojet");
  WriteHistos(H_2D[++index_2D],file,"5calojet");
  WriteHistos(H_d[++index_d],file,"0trackjet");
  WriteHistos(H_t[++index_t],file,"0trackjet");
  WriteHistos(H_2D[++index_2D],file,"0trackjet");
  WriteHistos(H_d[++index_d],file,"1trackjet");
  WriteHistos(H_t[++index_t],file,"1trackjet");
  WriteHistos(H_2D[++index_2D],file,"1trackjet");
  WriteHistos(H_d[++index_d],file,"2trackjet");
  WriteHistos(H_t[++index_t],file,"2trackjet");
  WriteHistos(H_2D[++index_2D],file,"2trackjet");
  WriteHistos(H_d[++index_d],file,"3trackjet");
  WriteHistos(H_t[++index_t],file,"3trackjet");
  WriteHistos(H_2D[++index_2D],file,"3trackjet");
  WriteHistos(H_d[++index_d],file,"4trackjet");
  WriteHistos(H_t[++index_t],file,"4trackjet");
  WriteHistos(H_2D[++index_2D],file,"4trackjet");
  WriteHistos(H_d[++index_d],file,"5trackjet");
  WriteHistos(H_t[++index_t],file,"5trackjet");
  WriteHistos(H_2D[++index_2D],file,"5trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_0calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_0calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_0calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_1calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_1calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_1calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_2calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_2calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_2calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_3calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_3calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_3calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_4calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_4calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_4calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_5calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_5calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_5calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_0trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_0trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_0trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_1trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_1trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_1trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_2trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_2trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_2trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_3trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_3trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_3trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_4trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_4trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_4trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_5trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_5trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_5trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_0calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_0calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_0calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_1calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_1calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_1calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_2calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_2calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_2calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_3calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_3calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_3calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_4calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_4calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_4calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_5calojet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_5calojet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_5calojet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_0trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_0trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_0trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_1trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_1trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_1trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_2trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_2trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_2trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_3trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_3trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_3trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_4trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_4trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_4trackjet");
  WriteHistos(H_d[++index_d],file,"Wselec_boost_5trackjet");
  WriteHistos(H_t[++index_t],file,"Wselec_boost_5trackjet");
  WriteHistos(H_2D[++index_2D],file,"Wselec_boost_5trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_0calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_0calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_0calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_1calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_1calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_1calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_2calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_2calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_2calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_3calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_3calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_3calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_4calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_4calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_4calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_5calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_5calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_5calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_0trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_0trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_0trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_1trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_1trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_1trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_2trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_2trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_2trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_3trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_3trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_3trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_4trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_4trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_4trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_5trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_5trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_5trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_0calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_0calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_0calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_1calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_1calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_1calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_2calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_2calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_2calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_3calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_3calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_3calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_4calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_4calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_4calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_5calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_5calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_5calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_0trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_0trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_0trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_1trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_1trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_1trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_2trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_2trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_2trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_3trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_3trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_3trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_4trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_4trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_4trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_Z_5trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_Z_5trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_Z_5trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_0calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_0calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_0calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_1calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_1calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_1calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_2calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_2calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_2calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_3calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_3calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_3calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_4calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_4calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_4calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_5calojet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_5calojet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_5calojet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_0trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_0trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_0trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_1trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_1trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_1trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_2trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_2trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_2trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_3trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_3trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_3trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_4trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_4trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_4trackjet");
  WriteHistos(H_d[++index_d],file,"Zselec_boost_M_5trackjet");
  WriteHistos(H_t[++index_t],file,"Zselec_boost_M_5trackjet");
  WriteHistos(H_2D[++index_2D],file,"Zselec_boost_M_5trackjet");
  
  
}
