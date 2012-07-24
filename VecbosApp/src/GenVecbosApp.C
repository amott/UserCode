//-------------------------------------------------------
// Description:
//    Routine to run Vecbos selection
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini
//    CERN
//-------------------------------------------------------

// W+jets GenApplication = 1 
// Z+jets GenApplication = 2
#define GenApplication 2

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

#if GenApplication == 1
#include <include/GenWjets.hh>
#endif

#if GenApplication == 2
#include <include/GenZjets.hh>
#endif

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  if ( argc < 5 ){
    cout << "Error at Input: please specify a list of input files, an output file a ref lumi (in pb-1) and the process cross section (in pb)" << endl; 
    cout << "Example:        ./GenVecbosApp list.txt output.root 100 1" << endl;
    return 1;
  }
  
  char inputFileName[150];
  char outFileName[150];
  float lumi, xsec;
  
  strcpy(inputFileName,argv[1]);
  strcpy(outFileName,argv[2]);
  sscanf(argv[3],"%f",&lumi);
  sscanf(argv[4],"%f",&xsec);

  std::string line;
  std::vector<std::string> theListOfInputFiles;
  ifstream myFileWithListOfInputFiles (inputFileName);
  if (myFileWithListOfInputFiles.is_open()) {
    while (! myFileWithListOfInputFiles.eof() )
      {
	getline (myFileWithListOfInputFiles,line);
	theListOfInputFiles.push_back(line);
      }
    myFileWithListOfInputFiles.close();
  }
  else {
    std::cout << "Cannot open " << inputFileName << std::endl; 
    return 1;
  }
  
  TChain *theChain = new TChain("t3");
  for(std::vector<std::string>::const_iterator i = theListOfInputFiles.begin(); i != theListOfInputFiles.end(); ++i) {
    theChain->Add(i->c_str());
    std::cout << "File added to TChain: " << *i << std::endl;
  }

#if GenApplication == 1
  GenWjets vecbos(theChain, xsec, lumi);
  vecbos.SetEtaMax(3.);
  vecbos.SetPtWMin(0.);
#endif
#if GenApplication == 2
  GenZjets vecbos(theChain, xsec, lumi);
  vecbos.SetEtaMax(3.);
  vecbos.SetPtZMin(0.);
#endif
  vecbos.Loop(string(outFileName));  
  
  return 0;
}
