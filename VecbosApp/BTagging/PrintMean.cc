#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdlib.h>

using namespace std;

int main(){

  string _file="tbarWDS_EffBRaw.txt";

  string outfile="tbarWDS_EffB.txt";
  
  ifstream inFile;
  char content[512];
  string Eff;
  double Mean=0;
  double counter=0;
  inFile.open(_file.c_str(), ifstream::in);
  if(inFile.good()) {
    while (!inFile.eof()) {
      inFile.getline(content,512);
      inFile >> Eff;
      counter+=1;
      if(Eff!="nan"){
	Mean+=atof(Eff.c_str());
      }
      else Mean+=0;
    }
  }
  inFile.close();

  ofstream OutFile;
  cout<<endl<<counter<<endl;
  Mean=Mean/counter;
  OutFile.open(outfile.c_str());
  for(double i=0; i < counter; i+=1)OutFile<<Mean<<endl;
  OutFile.close();

  return 0;
}
