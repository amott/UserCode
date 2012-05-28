//-------------------------------------------------------
// Description:
//    Routine to run Vecbos selection
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini
//    CERN
//-------------------------------------------------------

// C++ includes
#include <iostream>

#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

// Vecbos includes
#include <CommonTools/include/TriggerMask.hh>
#include <include/Vecbos.hh>

#include <include/HggReducer.hh>
#include <include/HggMakePhotonTree.hh>

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char** argv) {
  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[2000];
  char outFileName[2000];
  char cfg[400];
  char json[400]="none";

  if ( argc < 4 ){
    cout << "Error at Input: please specify an input file including the list of input ROOT files" << endl; 
    cout << "Example:        ./HggApp list.txt output.root ConfigFile" << endl;
    cout << "Available options: " <<endl;
    cout << "-weight=w  weight of the MC" << endl;  
    cout << "-start=N start from event N in the chain" << endl; 
    cout << "-stop=N stop at event N in the chain" << endl; 
    cout << "-signal=N 0=W(prompt),1=Z(prompt),2=W(other),3=Z(other), 4= no mctruth" << endl; 
    cout << "--isData to run on good runs on data" << endl; 
    cout << "--PhotonTree to make the photon tree for the regression" << endl;
    cout << "-json=file path of the json file" << endl;
    return 1;
  }

  // rad running options
  strcpy(inputFileName,argv[1]);
  strcpy(outFileName,argv[2]);
  strcpy(cfg,argv[3]);

  TChain *theChain = new TChain("ntp1");
  char Buffer[2000];
  char MyRootFile[2000];  
  ifstream *inputFile = new ifstream(inputFileName);
  // get the tree with the conditions from the first file
  //  TTree *treeCond = new TTree();
  //  int nfiles=1;
  char tmpFileName[2000];
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,2000);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	if(string(MyRootFile).find("eos") != std::string::npos) {
	  theChain->Add("root:/"+TString(MyRootFile));
        } else if(string(MyRootFile).find("castor") != std::string::npos) {
	  theChain->Add("rfio:"+TString(MyRootFile));
	} else{
	  theChain->Add(TString(MyRootFile));	 
	}
        // theChain->Add("root://castorcms/"+TString(MyRootFile));
	//        theChain->Add(TString(MyRootFile));
	std::cout << "chaining " << MyRootFile << std::endl;
	//	if ( nfiles==1 ) {
	//	  TFile *firstfile = TFile::Open("root://castorcms/"+TString(MyRootFile));
	//	  treeCond = (TTree*)firstfile->Get("Conditions");
	//	}
	//        nfiles++;
      }
  }
 

  //  theChain->MakeClass("thisiswhyitcrashed");

  inputFile->close();
  //delete inputFile;
  // get additional input options
  int signal = 0;
  int start = 0;
  int stop  = -1;
  bool isData = false;
  float lumi = -999.;
  float xsec = -999.;
  float weight = 1.;
  bool photonTree = false;
  for (int i=1;i<argc;i++){
    if (strncmp(argv[i],"-start",6)==0) sscanf(argv[i],"-start=%i",&start);
    if (strncmp(argv[i],"-stop",5)==0)  sscanf(argv[i],"-stop=%i",&stop);
    if (strncmp(argv[i],"-signal",7)==0)  sscanf(argv[i],"-signal=%i",&signal);
    if (strncmp(argv[i],"-weight",7)==0)  sscanf(argv[i],"-weight=%f",&weight);
    if (strncmp(argv[i],"--isData",8)==0)  isData = true;
    if (strncmp(argv[i],"-lumi",5)==0)  sscanf(argv[i],"-lumi=%f",&lumi);
    if (strncmp(argv[i],"-xsec",5)==0)  sscanf(argv[i],"-xsec=%f",&xsec);
    if (strncmp(argv[i],"-json",5)==0)  sscanf(argv[i],"-json=%s",&json);
    if (strncmp(argv[i],"--PhotonTree",12)==0)  photonTree = true;
  }
  
  cout << "Running on: " << stop << " Entries" << endl;

  if(!photonTree){
    HggReducer vecbos(theChain, string(json), isData, isData);
    cout << "Setting Config File:" << endl;
    cout << string(cfg) << endl;
    vecbos.setConfig(cfg);
    vecbos.Loop(string(outFileName), start, stop);
  }else{
    HggMakePhotonTree vecbos(theChain,string(outFileName),isData);
    vecbos.Loop(start,stop);
  }
  //  system("rm thisiswhyitcrashed*");

  //exit(0);
  return 0;
}
