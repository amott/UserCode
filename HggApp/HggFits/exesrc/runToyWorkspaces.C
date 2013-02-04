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

  a.addLongOption("SaveToyWorkspaces",ArgParser::noArg,"save the generated toy workspaces (for debugging!) (default: off)");
  a.addLongOption("RunData",ArgParser::noArg,"determine the S value for data");


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
  
  bool saveWS = a.longFlagPres("SaveToyWorkspaces");
  bool runData = a.longFlagPres("RunData");
  MakeSpinToy mst(wsFile);

  mst.setTargetLumi(tlumi);
  mst.setNominalLumi(nlumi);

  mst.setSaveWorkspaces(saveWS);
  mst.setDoData(runData);
  mst.runN(N);
  mst.save(outFile);



  return 0;

}
