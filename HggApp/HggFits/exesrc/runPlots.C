#include "ArgParser.hh"
#include "MakeSpinPlots.h"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("Lumi",ArgParser::required,"Luminosity");
  a.addArgument("OutputPath",ArgParser::required,"output directory");
  a.addArgument("OutputTag",ArgParser::required,"tag for the output plots");
  a.addLongOption("PrintOnly",ArgParser::noArg,"Only print the yields, don't make plots");

  string ret;
  if(a.process(ret) != 0){
    a.printOptions(argv[0]);
    return 0;
  }

  string inputWS = a.getArgument("InputWorkspace");
  float lumi = atof(a.getArgument("Lumi").c_str());
  string bp  = a.getArgument("OutputPath");
  string tag = a.getArgument("OutputTag");
  bool pOnly = a.longFlagPres("PrintOnly");

  MakeSpinPlots msp(inputWS,tag);

  msp.setLumi(lumi);
  msp.setBasePath(bp);

  if(!pOnly)  msp.runAll("Hgg125");

  
  std::cout << "Event Yields for H-->gg MC" <<std::endl <<std::endl;
  msp.printYields("Hgg125");

  std::cout << "Event Yields for RS Graviton-->gg MC" <<std::endl <<std::endl;
  msp.printYields("RSG125");
}
