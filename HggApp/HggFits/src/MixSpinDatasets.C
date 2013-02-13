#include "MixSpinDatasets.h"

MixSpinDatasets::MixSpinDatasets(RooWorkspace *w):ws(w)
{
  MakeSpinFits::getLabels("evtcat",&catNames,ws);
}

void MixSpinDatasets::mix(const char* mc1, const char* mc2, 
			  float f1, TString outputName){
  if(outputName=="") outputName = Form("%s_%0.2f_%s_$0.2f",mc1,f1,mc2,1-f1);

  internalMix(mc1,mc2,f1,outputName,"Combined");
  
}

void MixSpinDatasets::internalMix(const char* mc1, const char* mc2, 
				  float f1,TString outputName,TString cat){

  RooDataSet *mc1_d = (RooDataSet*)ws->data(Form("%s_%s",mc1,cat.Data()));
  RooDataSet *mc2_d = (RooDataSet*)ws->data(Form("%s_%s",mc2,cat.Data()));

  RooDataSet *mixed = new RooDataSet(outputName+"_"+cat,"",*(mc1_d->get(0)),"evtWeight");

  Long64_t iEntry=-1;
  const RooArgSet *set;
  while( (set=mc1_d->get(++iEntry)) ){ //loop over first DS
    float w = ((RooRealVar*)set->find("evtWeight"))->getVal();
    float we = ((RooRealVar*)set->find("evtWeight"))->getError();
    mixed->add( *set,w*f1,we*f1);
  }

  iEntry=-1;
  while( (set=mc2_d->get(++iEntry)) ){ //loop over first DS
    float w = ((RooRealVar*)set->find("evtWeight"))->getVal();
    float we = ((RooRealVar*)set->find("evtWeight"))->getError();
    mixed->add( *set,w*(1-f1),we*(1-f1));
  }
  ws->import( *mixed);
  delete mixed;
}
