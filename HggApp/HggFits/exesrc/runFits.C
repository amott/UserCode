#include "ArgParser.hh"
#include "MakeSpinFits.h"
#include "ReadConfig.hh"

#include <iostream>
#include <string>

int main(int argc, char** argv){
  using namespace std;
  
  ArgParser a(argc,argv);
  a.addArgument("InputWorkspace",ArgParser::required,"input workspace");
  a.addArgument("OutputWorkspace",ArgParser::required,"output workspace");
  a.addLongOption("BkgFit",ArgParser::reqArg,"Background Fit Type [poly,exp] (default: exp)");
  a.addLongOption("MCSamples",ArgParser::reqArg,"Specify the MC sample to process (comma separated) or none (default: all)");

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

  if(a.longFlagPres("MCSamples")) // specify the samples to test
    {
      vector<string> samples = ReadConfig::tokenizeString( a.getLongFlag("MCSamples"), ",");
      if(samples.at(0).compare("none")==0) samples.clear();
      if(msf.specifySamples(samples) != 0) 
	{//error in a sample specification
	  return -1;
	}
    }


  msf.setAddSWeight(true);
  msf.run();
  msf.save();

}
