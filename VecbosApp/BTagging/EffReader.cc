
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdlib.h>

using namespace std;

struct EffSample{
  vector<double> E;
  string sample;
};

double mean(vector <double> v){
  double N=0;
  double m=0;
  for(unsigned int i=0; i < v.size(); i++){
    m+=v[i];
    N+=1;
  }
  if(N>0)return (m/N);
  else return 0;

}

vector<EffSample> Read(string _file=""){
  
  vector<EffSample> Effs;
  EffSample tmp;
  string Eff;
  ifstream inFile;
  char content[512];
  int Nsample=0;

  inFile.open(_file.c_str(), ifstream::in);
  if(inFile.good()) {
    while (!inFile.eof()) {
      inFile.getline(content,512);
      inFile >> Eff; 
      if(Eff.find("/")!=string::npos){
	if(Nsample>0)Effs.push_back(tmp);
	tmp.sample=Eff;
	tmp.E.clear();
	Nsample++;
      }else{
	if(Eff!="nan")tmp.E.push_back(atof(Eff.c_str()));
	else tmp.E.push_back(0);
      }
    }
  }
  inFile.close();
  return Effs;

}  

bool isMean(string sample=""){
  bool isMean=true;
  if(sample.find("TTbarJets")!=string::npos)isMean=false;
  return isMean;
}

string outfilename(string sample, string flavor){

  string sampleShort="Unknown";

  if(sample.find("TTbarJets")!=string::npos)sampleShort="TTbarJets";
  if(sample.find("WJets")!=string::npos)sampleShort="WJets";
  if(sample.find("DYJets")!=string::npos)sampleShort="DYJets";
  if(sample.find("QCD_Pt-30")!=string::npos)sampleShort="QCD30";
  if(sample.find("QCD_Pt-50")!=string::npos)sampleShort="QCD50";
  if(sample.find("QCD_Pt-80")!=string::npos)sampleShort="QCD80";
  if(sample.find("QCD_Pt-120")!=string::npos)sampleShort="QCD120";
  if(sample.find("QCD_Pt-170")!=string::npos)sampleShort="QCD170";
  if(sample.find("QCD_Pt-300")!=string::npos)sampleShort="QCD300";
  if(sample.find("QCD_Pt-470")!=string::npos)sampleShort="QCD470";
  if(sample.find("QCD_Pt-600")!=string::npos)sampleShort="QCD600";
  if(sample.find("QCD_Pt-800")!=string::npos)sampleShort="QCD800";
  if(sample.find("QCD_Pt-1000")!=string::npos)sampleShort="QCD1000";
  if(sample.find("QCD_Pt-1400")!=string::npos)sampleShort="QCD1400";
  if(sample.find("tSchannel")!=string::npos)sampleShort="tS";
  if(sample.find("tTchannel")!=string::npos)sampleShort="tT";
  if(sample.find("tbarSchannel")!=string::npos)sampleShort="tbarS";
  if(sample.find("tbarTchannel")!=string::npos)sampleShort="tbarT";
  if(sample.find("tWDR")!=string::npos)sampleShort="tWDR";
  if(sample.find("tWDS")!=string::npos)sampleShort="tWDS";
  if(sample.find("tbarWDR")!=string::npos)sampleShort="tbarWDR";
  if(sample.find("tbarWDS")!=string::npos)sampleShort="tbarWDS";

  string filename=sampleShort+"_Eff"+flavor+".txt";
  return filename;  
}

void Print(EffSample S, string flavor){
  
  ofstream outfile;
  string filename=outfilename(S.sample, flavor);
  outfile.open(filename.c_str());
  outfile<<endl;

  bool ismean=isMean(S.sample);

  if(ismean){
    double Mean=mean(S.E);
    for(unsigned int j=0; j< S.E.size(); j++){
      outfile<<Mean<<endl;
    }
  }else{
    for(unsigned int j=0; j< S.E.size(); j++){
      outfile<<S.E[j]<<endl;
    }
  }

  outfile.close();

}

int main(){

  vector<EffSample> esL=Read("EffLTot.txt");
  vector<EffSample> esB=Read("EffBTot.txt");
  vector<EffSample> esC=Read("EffCTot.txt");

  cout<<esL.size()<<esC.size()<<esB.size()<<endl;

  for(unsigned int k=0; k<esL.size(); k++){
    Print(esL[k], "L");
  }
  for(unsigned int k=0; k<esB.size(); k++){
    Print(esB[k], "B");
  }
  for(unsigned int k=0; k<esC.size(); k++){
    Print(esC[k], "C");
  }


  return 0;
}
