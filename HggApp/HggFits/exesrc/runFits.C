#include "ArgParser.hh"
#include "MakeSpinFits.h"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("OutputWorkspace",ArgParser::required,"output workspace");
  a.addLongOption("useR9",ArgParser::noArg,"use r9 categories (default: off)");
  a.addLongOption("CombinedFit",ArgParser::noArg,"use a combined fit to the categories (default: off)");
  a.addLongOption("BkgFit",ArgParser::reqArg,"Background Fit Type [poly,exp] (default: exp)");

  string ret;
  if(a.process(ret) != 0){
    a.printOptions(argv[0]);
    return 0;
  }

  string inputWS = a.getArgument("InputWorkspace");
  string outputWS = a.getArgument("OutputWorkspace");
  bool useR9 = a.longFlagPres("useR9");
  bool useCombFit = a.longFlagPres("CombinedFit");
  string fit = "exp";
  if(a.longFlagPres("BkgFit")){
    fit = a.getLongFlag("BkgFit");    
  }

  
  MakeSpinFits msf(inputWS,outputWS);

  if(fit.compare("poly")==0) msf.setBkgFit(MakeSpinFits::kPoly);
  else msf.setBkgFit(MakeSpinFits::kExp);

  msf.setAddSWeight(true);
  msf.addMCLabel("Hgg125");
  msf.addMCLabel("RSG125");
  msf.setUseR9(useR9);
  msf.setUseCombinedFit(useCombFit);
  msf.run();
  msf.save();

}
