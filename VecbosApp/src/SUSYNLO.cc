// std includes
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

#include <TH2D.h>
#include <TString.h>
#include <SUSYNLO.hh>

using namespace std;

SUSYNLO::SUSYNLO(TString filename, TString label, int iM0, double minM0, double maxM0, int iM12, double minM12, double maxM12) {
  _filename = filename;
  ng = new TH2D("ng"+label, "ng"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  ns = new TH2D("ns"+label, "ns"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  nn = new TH2D("nn"+label, "nn"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  ll = new TH2D("ll"+label, "ll"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  sb = new TH2D("sb"+label, "sb"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  ss = new TH2D("ss"+label, "ss"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  tb = new TH2D("tb"+label, "tb"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  bb = new TH2D("bb"+label, "bb"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  gg = new TH2D("gg"+label, "gg"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);
  sg = new TH2D("sg"+label, "sg"+label, iM0, minM0, maxM0, iM12, minM12, maxM12);  
}

SUSYNLO::~SUSYNLO(){
}

void SUSYNLO::SetXsec() {
  string line;
  float m0, m12, tanb, A0;
  float scale;
  char signMu[10];
  float myng = -99.;
  float myns = -99.;
  float mynn = -99.;
  float myll = -99.;
  float mysb = -99.;
  float myss = -99.;
  float mytb = -99.;
  float mybb = -99.;
  float mygg = -99.;
  float mysg = -99.;
  ifstream myfile (_filename.Data());
  if (myfile.is_open()) {
    while ( myfile.good()) {
      getline (myfile,line);
      if(line.find("Sub-processes") != string::npos) continue;      
      sscanf(line.c_str()," |(scale=%e) m0=%e, m1/2=%e, tanbeta=%e, A0=%e, sign(mu)=%s %e | %e | %e | %e | %e | %e | %e | %e | %e | %e |",
	     &scale, &m0, &m12, &tanb, &A0, signMu, &myng, &myns, &mynn, &myll, &mysb, &myss, &mytb, &mybb, &mygg, &mysg);
      ng->SetBinContent(ng->FindBin(m0,m12), double(myng));
      ns->SetBinContent(ns->FindBin(m0,m12), double(myns));
      nn->SetBinContent(nn->FindBin(m0,m12), double(mynn));
      ll->SetBinContent(ll->FindBin(m0,m12), double(myll));
      sb->SetBinContent(sb->FindBin(m0,m12), double(mysb));
      ss->SetBinContent(ss->FindBin(m0,m12), double(myss));
      tb->SetBinContent(tb->FindBin(m0,m12), double(mytb));
      bb->SetBinContent(bb->FindBin(m0,m12), double(mybb));
      gg->SetBinContent(gg->FindBin(m0,m12), double(mygg));
      sg->SetBinContent(sg->FindBin(m0,m12), double(mysg));
    }
    myfile.close();
  }
  else cout << "Error in SUSYNLO: Unable to open file " << _filename.Data(); 
}

void SUSYNLO::SetKfactor() {
  string line;
  int m0, m12, tanb, A0, signMu;
  float myng = -99.;
  float myns = -99.;
  float mynn = -99.;
  float myll = -99.;
  float mysb = -99.;
  float myss = -99.;
  float mytb = -99.;
  float mybb = -99.;
  float mygg = -99.;
  float mysg = -99.;
  ifstream myfile (_filename.Data());
  if (myfile.is_open()) {
    while ( myfile.good()) {
      getline (myfile,line);
      if(line.find("Sub-processes") != string::npos) continue;      
      sscanf(line.c_str()," | msugra_%i_%i_%i_%i_%i.slha | %e | %e | %e | %e | %e | %e | %e | %e | %e | %e |",
	     &m0, &m12, &tanb, &signMu, &A0, &myng, &myns, &mynn, &myll, &mysb, &myss, &mytb, &mybb, &mygg, &mysg);
      ng->SetBinContent(ng->FindBin(m0,m12), double(myng));
      ns->SetBinContent(ns->FindBin(m0,m12), double(myns));
      nn->SetBinContent(nn->FindBin(m0,m12), double(mynn));
      ll->SetBinContent(ll->FindBin(m0,m12), double(myll));
      sb->SetBinContent(sb->FindBin(m0,m12), double(mysb));
      ss->SetBinContent(ss->FindBin(m0,m12), double(myss));
      tb->SetBinContent(tb->FindBin(m0,m12), double(mytb));
      bb->SetBinContent(bb->FindBin(m0,m12), double(mybb));
      gg->SetBinContent(gg->FindBin(m0,m12), double(mygg));
      sg->SetBinContent(sg->FindBin(m0,m12), double(mysg));
    }
    myfile.close();
  }
  else cout << "Error in SUSYNLO: Unable to open file " << _filename.Data(); 
}

double SUSYNLO::GetVal(double m0, double m12, int iP1, int iP2) {
  double _iP1 = abs(iP1);
  double _iP2 = abs(iP2);
  double xsec = 0; 

  if((NeuCha(_iP1)  && Gluino(_iP2)) || (NeuCha(_iP2) && Gluino(_iP1)) ) xsec = ng->GetBinContent(ng->FindBin(m0, m12));
  else if((NeuCha(_iP1)  && Squark(_iP2)) || (NeuCha(_iP2) && Squark(_iP1)) ) xsec = ns->GetBinContent(ns->FindBin(m0, m12));
  else if((NeuCha(_iP1)  && NeuCha(_iP2)) || (NeuCha(_iP2) && NeuCha(_iP1)) ) xsec = nn->GetBinContent(nn->FindBin(m0, m12));
  else if((Slepton(_iP1)  && Slepton(_iP2)) || (Slepton(_iP2) && Slepton(_iP1)) ) xsec = ll->GetBinContent(ll->FindBin(m0, m12));
  else if(((Squark(_iP1)  && Squark(_iP2)) || (Squark(_iP2) && Squark(_iP1))) && iP1*iP2<0) xsec = sb->GetBinContent(sb->FindBin(m0, m12));
  else if(((Squark(_iP1)  && Squark(_iP2)) || (Squark(_iP2) && Squark(_iP1))) && iP1*iP2>0) xsec = ss->GetBinContent(ss->FindBin(m0, m12));
  else if((Stop(_iP1)  && Stop(_iP2)) || (Stop(_iP2) && Stop(_iP1)) ) xsec = tb->GetBinContent(tb->FindBin(m0, m12));
  else if((Sbottom(_iP1)  && Sbottom(_iP2)) || (Sbottom(_iP2) && Sbottom(_iP1)) ) xsec = bb->GetBinContent(bb->FindBin(m0, m12));
  else if((Gluino(_iP1)  && Gluino(_iP2)) || (Gluino(_iP2) && Gluino(_iP1)) ) xsec = gg->GetBinContent(gg->FindBin(m0, m12));
  else if((Squark(_iP1)  && Gluino(_iP2)) || (Squark(_iP2) && Gluino(_iP1)) ) xsec = sg->GetBinContent(sg->FindBin(m0, m12));

  return xsec;
}

bool SUSYNLO::NeuCha(int iP) {
  bool out = false;
  if(abs(iP) >= 1000022 && abs(iP) <= 1000037) out = true;
  return out;
}

bool SUSYNLO::Squark(int iP) {
  bool out = false;
  if(abs(iP) >= 1000001 && abs(iP) <= 1000004) out = true; 
  if(abs(iP) >= 2000001 && abs(iP) <= 2000004) out = true; 
  return out;
}

bool SUSYNLO::Slepton(int iP) {
  bool out = false;
  if(abs(iP) >= 1000011 && abs(iP) <= 1000016) out = true; 
  if(abs(iP) >= 2000011 && abs(iP) <= 20000164) out = true; 
  return out;
}

bool SUSYNLO::Gluino(int iP) {
  bool out = false;
  if(abs(iP) == 1000021) out = true; 
  return out;
}

bool SUSYNLO::Stop(int iP) {
  bool out = false;
  if(abs(iP) == 1000006 || 
     abs(iP) == 2000006) out = true; 
  return out;
}

bool SUSYNLO::Sbottom(int iP) {
  bool out = false;
  if(abs(iP) == 1000005 || 
     abs(iP) == 2000005) out = true; 
  return out;
}


