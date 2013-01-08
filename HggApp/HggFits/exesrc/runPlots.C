#include "ArgParser.hh"
#include "MakeSpinPlots.h"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("Lumi",ArgParser::required,"Luminosity");
  a.addArgument("output path",ArgParser::required,"output directory");
  a.addArgument("output tag",ArgParser::required,"tag for the output plots");
  a.addLongOption("useR9",ArgParser::noArg,"use r9 categories (default: off)");

  string ret;
  if(a.process(ret) != 0){
    a.printOptions(argv[0]);
    return 0;
  }

  string inputWS = a.getArgument("InputWorkspace");
  bool useR9 = a.longFlagPres("useR9");
  float lumi = atof(a.getArgument("Lumi").c_str());
  string bp  = a.getArgument("output path");
  string tag = a.getArgument("output tag");

  MakeSpinPlots msp(inputWS,tag);

  msp.addMCType("Hgg125");
  msp.addMCType("RSG125");
  msp.setUseR9(useR9);
  msp.setLumi(lumi);
  msp.setBasePath(bp);

  msp.runAll("Hgg125");
}
