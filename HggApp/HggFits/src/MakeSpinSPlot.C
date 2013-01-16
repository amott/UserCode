#include "MakeSpinSPlot.h"


MakeSpinSPlot::MakeSpinSPlot(RooAbsData *data){
  __dataSet = data;
  __nSpec=0;
  __covMatrix=0;
}

void MakeSpinSPlot::addSpecies(TString name, RooAbsPdf* pdf, double expYield){
  __pdfs.push_back(pdf);
  __speciesNames.push_back(name);
  __expectedYields.push_back(expYield);
  __nSpec++;
}


void MakeSpinSPlot::computeCovMatrix(){

  __covMatrix = new TMatrixD(__nSpec,__nSpec);

  Long64_t iEntry=-1;
  while( __dataSet->get(++iEntry) ){
    
    for(int iRow = 0; iRow<__nSpec;iRow++){
      for(int iCol = 0; iCol<__nSpec;iCol++){
	double den = computeDenom();
	double num = __pdfs.at(iRow)->getVal()*__pdfs.at(iCol)->getVal();
	(*__covMatrix)[iRow][iCol] += num/den;
      }
    }
  }//while

  __covMatrix->Invert();
}

double MakeSpinSPlot::computeDenom(){ //compute the denominator for the covariance matrix

  std::vector<RooAbsPdf*>::iterator pdfIt = __pdfs.begin();
  std::vector<double>::iterator expIt = __expectedYields.begin();

  double denom=0;
  for(; pdfIt != __pdfs.end(); pdfIt++, expIt++){
    denom+= (*expIt) * ((*pdfIt)->getVal());
  }
  return denom;
}

void MakeSpinSPlot::computeSWeight(){

  RooArgSet swVars;
  std::vector<TString>::const_iterator specIt = __speciesNames.begin();
  for(; specIt != __speciesNames.end(); specIt++){
    swVars.add( *(new RooRealVar( Form("%s_sw",specIt->Data()),"",-1e+6,1e+6)) );
  }
  
  __sWeightDataSet = new RooDataSet("SWeightDataSet","",swVars);
  
  Long64_t iEntry=-1;

  while( __dataSet->get(++iEntry) ){

    specIt = __speciesNames.begin();
    for(; specIt != __speciesNames.end(); specIt++){ // loop over the species
      int iRow = specIt - __speciesNames.begin();

      double num=0;
      std::vector<RooAbsPdf*>::iterator pdfIt = __pdfs.begin();
      for(; pdfIt != __pdfs.end(); pdfIt++){
	int iCol = pdfIt - __pdfs.begin();
	
	num += (*__covMatrix)[iRow][iCol] * ((*pdfIt)->getVal());
      }

      double denom = computeDenom();
      ((RooRealVar*)swVars.find( Form( "%s_sw", specIt->Data()) ))->setVal(num/denom);
    }//end loop over species

    __sWeightDataSet->add(swVars);
  }  //end while loop

}

void MakeSpinSPlot::calculate(){
  if(__nSpec<0) return;
  computeCovMatrix();
  computeSWeight();
}
