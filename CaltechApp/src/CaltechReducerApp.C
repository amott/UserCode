//-------------------------------------------------------
// Description:
//    Routine to run the Caltech Reducer
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//    Maurizio Pierini
//    CERN
//    Alex Mott
//    Caltech
//-------------------------------------------------------

// C++ includes
#include <iostream>

#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

// Vecbos includes
#include <CommonTools/include/TriggerMask.hh>
#include <include/Vecbos.hh>
#include <ArgParser.hh>

#include <include/CaltechReducer.hh>

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char** argv) {
  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[2000];
  bool isList=true;
  char outFileName[2000];
  char cfg[400];
  char json[400]="none";
  bool verbose_flag=0;
  bool isData=false;
  float xsec=-999;
  int start=0;
  int stop=-1;

  ArgParser parser(argc,argv);
  parser.addLongOption("isData",ArgParser::noArg,"Specify data or MC");
  parser.addLongOption("json",ArgParser::reqArg,"Specify the JSON to use");
  parser.addLongOption("config",ArgParser::reqArg, "Specify the configuration file");
  parser.addLongOption("start",ArgParser::reqArg, "first event to process");
  parser.addLongOption("stop",ArgParser::reqArg, "last event to process");
  parser.addShortOption('f',ArgParser::noArg,"Specify whether the input file is a single file"); //use this if we only specify one file, not a list
  parser.addArgument("inputFile",ArgParser::required,"input file, expects a list unless -f is set");
  parser.addArgument("outputFile",ArgParser::required,"output file");

  std::string error;
  int retCode = parser.process(error);

  if(retCode != 0){
    cout << "ERROR Parsing option: " << error << "   return code: " << retCode << endl;
    parser.printOptions("CaltechReducerApp");
    return 0;
  }
  isData = parser.longFlagPres("isData");
  strncpy(json,parser.getLongFlag("json").c_str(),400);
  strncpy(cfg,parser.getLongFlag("config").c_str(),400);
  isList = !parser.shortFlagPres('f');

  if(parser.longFlagPres("start")) start = atoi(parser.getLongFlag("start").c_str());
  if(parser.longFlagPres("stop")) stop = atoi(parser.getLongFlag("stop").c_str());

  strncpy(inputFileName,parser.getArgument("inputFile").c_str(),2000);
  strncpy(outFileName,parser.getArgument("outputFile").c_str(),2000);
  
  cout << "Running with Arguments: " << endl
       << "isData: " << isData << endl
       << "json: " << json << endl
       << "isList: " << isList << endl
       << "config: " << cfg << endl
       << "inputFile: " << inputFileName << endl
       << "outputFile: " << outFileName << endl;

  TChain *theChain = new TChain("ntp1");
  char Buffer[2000];
  char MyRootFile[2000];  
  if(isList){
    ifstream *inputFile = new ifstream(inputFileName);
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
    inputFile->close();
    
  }else{
    if(string(inputFileName).find("eos") != std::string::npos) {
      theChain->Add("root:/"+TString(inputFileName));
    } else if(string(inputFileName).find("castor") != std::string::npos) {
      theChain->Add("rfio:"+TString(inputFileName));
    } else{
      theChain->Add(TString(inputFileName));	 
    }
    std::cout << "chaining " << inputFileName << std::endl;
  }

    //delete inputFile;
  // get additional input options

  CaltechReducer vecbos(theChain, string(json), isData, isData);
  cout << "Setting Config File:" << endl;
  cout << string(cfg) << endl;
  vecbos.setConfig(cfg);
  vecbos.Loop(string(outFileName), start,stop);
  return 0;
}
