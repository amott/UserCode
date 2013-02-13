#ifndef MixSpinDatasets_h
#define MixSpinDatasets_h

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"

#include "TString.h"

#include <vector>
#include "MakeSpinFits.h"

class MixSpinDatasets{
public:
  MixSpinDatasets(RooWorkspace *w);

  void mix(const char* mc1, const char* mc2, float f1,TString outputName="");
protected:
  RooWorkspace *ws;

  std::vector<TString> catNames;

  void internalMix(const char* mc1, const char* mc2, float f1,TString outputName,TString cat);
};

#endif
