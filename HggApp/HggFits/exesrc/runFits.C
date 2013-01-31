#include "ArgParser.hh"
#include "MakeSpinFits.h"

#include <iostream>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("OutputWorkspace",ArgParser::required,"output workspace");
  a.addLongOption("BkgFit",ArgParser::reqArg,"Background Fit Type [poly,exp] (default: exp)");

  string ret;
  if(a.process(ret) != 0){
    a.printOptions(argv[0]);
    return 0;
  }

  string inputWS = a.getArgument("InputWorkspace");
  string outputWS = a.getArgument("OutputWorkspace");
  string fit = "exp";
  if(a.longFlagPres("BkgFit")){
    fit = a.getLongFlag("BkgFit");    
  }

  
  MakeSpinFits msf(inputWS,outputWS);

  if(fit.compare("poly")==0) msf.setBkgFit(MakeSpinFits::kPoly);
  else msf.setBkgFit(MakeSpinFits::kExp);

  msf.setAddSWeight(true);
  msf.run();
  msf.save();

}
