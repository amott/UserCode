#include "ArgParser.hh"
#include "MakeSpinToyWorkspace.h"
#include "ReadConfig.hh"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  ArgParser a(argc,argv);
  a.addArgument("WorkspaceFile",ArgParser::required,"path to the workspace file");
  a.addArgument("outputFile",ArgParser::required,"path to the output file");

  a.addArgument("N",ArgParser::required,"Number of Toys");
  a.addArgument("tag",ArgParser::required,"");
  a.addArgument("lumi",ArgParser::required,"effective luminosity");
  //a.addLongOption("SaveWorkspaces",Arg);

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  string wsFile = a.getArgument("WorkspaceFile");
  string outFile = a.getArgument("outputFile");

  int N = atoi(a.getArgument("N").c_str());
  string tag = a.getArgument("tag");
  float lumi = atof(a.getArgument("lumi").c_str());
  

  MakeSpinToyWorkspace mstw(wsFile);

  mstw.setup(tag,lumi);

  mstw.generateN(N);
  mstw.save(outFile);



  return 0;

}
