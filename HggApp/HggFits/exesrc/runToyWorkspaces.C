#include "ArgParser.hh"
#include "MakeSpinToy.h"
#include "ReadConfig.hh"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  ArgParser a(argc,argv);
  a.addArgument("WorkspaceFile",ArgParser::required,"path to the workspace file");
  a.addArgument("outputFile",ArgParser::required,"path to the output file");

  a.addArgument("N",ArgParser::required,"Number of Toys");
  a.addArgument("TargetLumi",ArgParser::required,"target luminosity");
  a.addArgument("InputLumi",ArgParser::required,"input luminosity");

  a.addLongOption("useR9",ArgParser::noArg,"use r9 categories (default: off)");


  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  string wsFile = a.getArgument("WorkspaceFile");
  string outFile = a.getArgument("outputFile");

  int N = atoi(a.getArgument("N").c_str());
  float tlumi = atof(a.getArgument("TargetLumi").c_str());
  float nlumi = atof(a.getArgument("InputLumi").c_str());
  
  bool useR9 = a.longFlagPres("useR9");
  MakeSpinToy mst(wsFile);

  mst.setTargetLumi(tlumi);
  mst.setNominalLumi(nlumi);

  mst.setUseR9(useR9);
  mst.runN(N);
  mst.save(outFile);



  return 0;

}
