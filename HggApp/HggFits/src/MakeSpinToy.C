#include "MakeSpinToy.h"
#include <iostream>

MakeSpinToy::MakeSpinToy(TString fileName,TString wsName):
  nCat(MakeSpinWorkspace::nCat){
  TFile *f = new TFile(fileName);
  ws = (RooWorkspace*)f->Get(wsName);
  mass = ws->var("mass");
  cosT = ws->var("cosT");
  cosT->setRange(-0.8,1);
  cosT->setBins(9);
  S = new RooRealVar("S","",0,-1e6,1e6);

  GenMinusFit = new RooRealVar("NgenMinusNfit","",0,-1e6,1e6);

  S_TruthHgg = new RooDataSet("S_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_TruthRSG = new RooDataSet("S_TruthRSG","",RooArgSet(*S,*GenMinusFit));

  S_splot_TruthHgg = new RooDataSet("S_splot_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_splot_TruthRSG = new RooDataSet("S_splot_TruthRSG","",RooArgSet(*S,*GenMinusFit));

  S_2D_TruthHgg = new RooDataSet("S_2D_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_2D_TruthRSG = new RooDataSet("S_2D_TruthRSG","",RooArgSet(*S,*GenMinusFit));

  S_tot_TruthHgg = new RooDataSet("S_tot_TruthHgg","",RooArgSet(*S,*GenMinusFit));
  S_tot_TruthRSG = new RooDataSet("S_tot_TruthRSG","",RooArgSet(*S,*GenMinusFit));

  setSaveWorkspaces(false);

  mcLabels[0]="Hgg125";
  mcLabels[1]="RSG125";
}

void MakeSpinToy::generateToyWorkspace(RooWorkspace &toyws,genType gen){
  RooCategory *cat = new RooCategory("evtcat","evtcat");
  int i=0;
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++,i++){
    std::cout << *it << "   " << i << std::endl;
    cat->defineType( *it,i);
  }
  RooDataSet* dataComb = new RooDataSet("Data_Combined","",RooArgSet(*mass,*cosT,*cat) );

  TRandom3 rng(0);
  float Nsig = rng.Poisson( targetLumi/12.*607 );
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++,i++){
    generateToyWorkspace(toyws,it->Data(),gen,Nsig);
    RooDataSet* d = (RooDataSet*)toyws.data( Form("Data_%s",it->Data()) );
    RooDataSet* tmp = new RooDataSet("DataCat+"+*it,"",RooArgSet(*mass,*cosT),RooFit::Index(*cat),RooFit::Import(*it,*d) );
    dataComb->append(*tmp);
    toyws.import(*(ws->var(Form("Hgg125_FIT_%s_sigmaEff",it->Data()) )));
    toyws.import(*(ws->var(Form("RSG125_FIT_%s_sigmaEff",it->Data()) )));
  }
  toyws.import(*cat);
  toyws.import(*dataComb);
  toyws.import(*(ws->data("Hgg125_Combined")));
  toyws.import(*(ws->data("RSG125_Combined")));
  toyws.import(*(ws->var("Hgg125_FIT_Combined_sigmaEff")));
  toyws.import(*(ws->var("RSG125_FIT_Combined_sigmaEff")));
}

void MakeSpinToy::generateToyWorkspace(RooWorkspace &toyws,const char* cat,genType gen,float nSigTot){
  //get the generation PDFs

  //mass
  RooAbsPdf* hggMassPdf = ws->pdf( Form("Hgg125_FIT_%s",cat) );
  RooAbsPdf* rsgMassPdf = ws->pdf( Form("RSG125_FIT_%s",cat) );

  ws->var(Form("Hgg125_FIT_%s_mean",cat))->setVal(125.);
  ws->var(Form("RSG125_FIT_%s_mean",cat))->setVal(125.);

  RooAbsPdf* bkgMassPdf = ws->pdf( Form("Data_BKGFIT_%s_bkgModel",cat) );

  //cosT
  RooDataSet *tmp  = (RooDataSet*)(ws->data(Form("Data_%s",cat))->reduce("(mass>110 && mass<120) || (mass>130 && mass<140)"));
  RooDataSet *tmpH = (RooDataSet*)ws->data(Form("Hgg125_%s",cat));
  RooDataSet *tmpR = (RooDataSet*)ws->data(Form("RSG125_%s",cat));

  RooKeysPdf* bkgGenPdf = new RooKeysPdf("bkgGenPdf","",*cosT,*tmp);
  RooKeysPdf* hggGenPdf = new RooKeysPdf("hggGenPdf","",*cosT,*tmpH);
  RooKeysPdf* rsgGenPdf = new RooKeysPdf("rsgGenPdf","",*cosT,*tmpR);
  /*
  tmp->SetName(Form("OriginalData_%s",cat));
  RooDataHist tmpHist("bkgGenHist","",*cosT,*tmp);
  RooHistPdf* bkgGenPdf = new RooHistPdf("bkgGenPdf","",*cosT,tmpHist);
  
  std::cout << "BKG GEN PDF:  " << bkgGenPdf <<std::endl;
  if(!bkgGenPdf) return;
  RooHistPdf* hggGenPdf = (RooHistPdf*)ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",cat));
  RooHistPdf* rsgGenPdf = (RooHistPdf*)ws->pdf(Form("RSG125_FIT_%s_cosTpdf",cat));
  */
  std::cout <<hggGenPdf <<std::endl;


  //get the expected number of events in this categoy
  //Data_Hgg125_FULAALFIT_EB_0_Nbkg
  int Nbkg = ((RooFormulaVar*)ws->obj(Form("Data_%s_FULLFIT_%s_Nbkg",mcLabels[gen].Data(),cat)))->getVal()*targetLumi/nominalLumi;
  int Nsig = nSigTot * ws->var(Form("Data_%s_FULLFIT_%s_fsig",mcLabels[gen].Data(),cat))->getVal()*targetLumi/nominalLumi;

  TRandom3 rng(0);
  int thisNbkg = (int)rng.Poisson(Nbkg);
  int thisNsig = Nsig;//(int)rng.Poisson(Nsig); //number of events to generate

  std::cout << "\nGenerating  " << thisNbkg << " background and " << thisNsig << " signal events\n" <<std::endl;

  RooRealVar nSigGen(Form("N_gen_sig_%s",cat),"",thisNsig);
  RooRealVar nBkgGen(Form("N_gen_bkg_%s",cat),"",thisNbkg);

  RooDataSet *M = bkgMassPdf->generate(*mass,thisNbkg);
  RooDataSet *C = bkgGenPdf->generate(*cosT,thisNbkg,RooFit::AutoBinned(kFALSE));


  RooDataSet *mcM, *mcC;
  switch(gen){
  case Hgg125:
    std::cout << "GENERATING HIGGS SPIN DISTRIBUTION" <<std::endl;
    mcM = hggMassPdf->generate(*mass,thisNsig);
    mcC = hggGenPdf->generate(*cosT,thisNsig,RooFit::AutoBinned(kFALSE));
    break;
    
  case RSG125:
    std::cout << "GENERATING RS GRAVITON SPIN DISTRIBUTION" <<std::endl;
    mcM = hggMassPdf->generate(*mass,thisNsig); // still generate the peak according to the higgs line shape
    mcC = rsgGenPdf->generate(*cosT,thisNsig,RooFit::AutoBinned(kFALSE));
    break;
    
  default:
    return;
  }
  std::cout << "M size: " << M->sumEntries() << std::endl << "C size: " << C->sumEntries() <<std::endl;
  std::cout << "M size: " << mcM->sumEntries() << std::endl << "C size: " << mcC->sumEntries() <<std::endl;
  
  M->Print();
  C->Print();
  mcM->Print();
  mcC->Print();
  
  M->merge(C);
  mcM->merge(mcC);
  M->append(*mcM);

  mcM->SetName(Form("ToySigData_%s",cat));
  M->SetName(Form("Data_%s",cat));
  
  toyws.import(*M);
  toyws.import(*mcM);
  toyws.import(*hggMassPdf);
  toyws.import(*rsgMassPdf);
  toyws.import(nSigGen);
  toyws.import(nBkgGen);

  delete bkgGenPdf;
}
 
double MakeSpinToy::computeLL(RooAbsPdf* pdf, RooAbsData* data,RooRealVar* var, int rebin){
  //TH1F dataHist("dataHist","",nBins,-1,1);

  //data->fillHistogram(&dataHist,*var);

  TH1F* dataHist = getHistogram(data,"dataHist",rebin);
  TH1F* pdfHist = (TH1F*)pdf->createHistogram("cosT",nBins/rebin);

  std::cout << dataHist->Integral() << std::endl;
  pdfHist->Scale(dataHist->Integral()/pdfHist->Integral());


  std::cout << "Bins:" <<std::endl;
  for(int i=0;i<cosT->getBins();i++){
    std::cout << "\t" << i+1 << ":    " << dataHist->GetBinContent(i+1) << " +- " << dataHist->GetBinError(i+1) << "          pdf: " 
	      << pdfHist->GetBinContent(i+1) << std::endl;
  }

  double p = dataHist->Chi2Test(pdfHist,"WW");

  std::cout << "p: " << p << std::endl;

  delete pdfHist;
  delete dataHist;
  if(p==0) return 0;
  return TMath::Log(p);

  /*
  RooDataHist dataHist("tmp","",*var,*data);
  return pdf->createNLL(dataHist)->getVal();
  pdf->fitTo(*data,RooFit::Extended(),RooFit::Strategy(0),RooFit::SumW2Error(kFALSE));
  RooFitResult* res = pdf->fitTo(*data,RooFit::Extended(),RooFit::Save(),RooFit::Strategy(2),RooFit::SumW2Error(kFALSE));
  std::cout << res->minNll() << std::endl;
  return res->minNll();
  */
}

TH1F* MakeSpinToy::getHistogram(RooAbsData* data, TString title, int rebin){
  TH1F* out = new TH1F(title,"",cosT->getBins()/rebin,-1,1);

  int i=0;
  while(data->get(i)){
    out->SetBinContent(i+1,data->weight());
    out->SetBinError(i+1,data->weightError());
    i++;
  }
  return out;
}

TTree* MakeSpinToy::makeForCombineTool(TString treeName, RooAbsData* hggData, RooAbsData* rsgData){
  TTree * out = new TTree(treeName,"");
  float q;
  float gMf;
  int type;

  out->Branch("q",&q);
  out->Branch("NgenMinusNfit",&gMf);
  out->Branch("type",&type,"type/I");

  Long64_t iEntry=-1;
  type=1; //SM HIggs
  const RooArgSet *set;

  while( (set=hggData->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    gMf = ((RooRealVar*)set->find("NgenMinusNfit"))->getVal();
    out->Fill();
  }
  iEntry=-1;
  type=-1;

  while( (set=rsgData->get(++iEntry)) ){ 
    q = ((RooRealVar*)set->find("S"))->getVal();
    gMf = ((RooRealVar*)set->find("NgenMinusNfit"))->getVal();
    out->Fill();
  }
  return out;
}

float MakeSpinToy::getExpEvents(float lumi, TString cat,TString mcType,RooWorkspace *toyws){
  //double tot    = ws->data(Form("%s_Combined",mcType.Data()))->sumEntries();
  //double thisN  = ws->data(Form("%s_%s",mcType.Data(),cat.Data()))->sumEntries();

  double f = ws->var(Form("Data_%s_FULLFIT_%s_fsig",mcType.Data(),cat.Data()))->getVal();
  
  if(toyws){
    toyws->import(*ws->var(Form("%s_EB_totalEvents",mcType.Data())));
    toyws->import(*ws->var(Form("%s_EE_totalEvents",mcType.Data())));
  }

  return f*lumi/12*607; //607 events in 12/fb @ 8 TeV
}


double* MakeSpinToy::run1(genType gen, int& N){
  TString mcType = "Hgg125";
  if(gen==RSG125) mcType="RSG125";

  //we currently have 3 outputs
  N=4;
  double * out = new double[N];

  RooWorkspace *toyws = new RooWorkspace(Form("toyws_%s",mcType.Data()),"toyws");
  generateToyWorkspace(*toyws,gen);


  MakeSpinFits fits("","");
  fits.setWorkspace(toyws);
  
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++){
    toyws->import(*ws->data(Form("Hgg125_%s",it->Data())));
    toyws->import(*ws->data(Form("RSG125_%s",it->Data())));

    toyws->import(*ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",it->Data())));
    toyws->import(*ws->pdf(Form("RSG125_FIT_%s_cosTpdf",it->Data())));

    fits.MakeSignalFitForFit(*it,"Hgg125");
    fits.MakeSignalFitForFit(*it,"RSG125");
    fits.MakeBackgroundOnlyFit(*it);    
  }
  fits.MakeCombinedSignalTest("Hgg125");
  fits.MakeCombinedSignalTest("RSG125");
  fits.AddCombinedBkgOnlySWeight("Hgg125");
  fits.AddCombinedBkgOnlySWeight("RSG125");

  fits.Make2DCombinedSignalTest("Hgg125","Hgg125");
  fits.Make2DCombinedSignalTest("Hgg125","RSG125");
  toyws->Print();

  double hggll=0,rsgll=0;
  double shggll=0,srsgll=0;
  double nFit=0,nGen=0;
  for( std::vector<TString>::const_iterator it = catLabels.begin();
       it != catLabels.end(); it++){
    //fits.MakeBackgroundFit("Hgg125",*it,125,2,false);

    fits.getSimpleBkgSubtraction("Hgg125",*it);
    fits.getSimpleBkgSubtraction("RSG125",*it);
    RooDataHist *thisBkgHgg = (RooDataHist*)toyws->data( Form("Data_Hgg125_%s_bkgSub_cosT",it->Data()));
    RooDataHist *thisBkgRSG = (RooDataHist*)toyws->data( Form("Data_RSG125_%s_bkgSub_cosT",it->Data()));
    RooDataSet *thisSdataHgg = (RooDataSet*)toyws->data( Form("Data_%s_Hgg125_sigWeight",it->Data()));
    RooDataSet *thisSdataRSG = (RooDataSet*)toyws->data( Form("Data_%s_RSG125_sigWeight",it->Data()));
    RooDataHist *thisShistHgg = new RooDataHist("tmpHgg","",*cosT,*thisSdataHgg);
    RooDataHist *thisShistRSG = new RooDataHist("tmpRSG","",*cosT,*thisSdataRSG);

    hggPdf = (RooHistPdf*)ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",it->Data()));
    rsgPdf = (RooHistPdf*)ws->pdf(Form("RSG125_FIT_%s_cosTpdf",it->Data()));


    if(gen==Hgg125) nFit+= thisSdataHgg->sumEntries();
    else nFit+= thisSdataRSG->sumEntries();
    nGen+=toyws->var( Form("N_gen_sig_%s",it->Data()) )->getVal();

    if(saveWorkspaces){
      toyws->import(*hggPdf);
      toyws->import(*rsgPdf);
    }
    
    std::cout << "hgg" <<std::endl;
    hggll += computeLL(hggPdf,thisBkgHgg,cosT,3);
    std::cout << "rsg" <<std::endl;
    //rsgll += computeLL(rsgPdf,thisBkgRSG,cosT,3);
    rsgll += computeLL(rsgPdf,thisBkgHgg,cosT,3);
    
    std::cout << "hgg splot" <<std::endl;
    shggll += computeLL(hggPdf,thisShistHgg,cosT,3);
    std::cout << "rsg splot" <<std::endl;
    //srsgll += computeLL(rsgPdf,thisShistRSG,cosT,3);
    srsgll += computeLL(rsgPdf,thisShistHgg,cosT,3);

    std::cout << "Bkg-sub:" << std::endl;
    std::cout << "Hgg: " << hggll << std::endl;
    std::cout << "RSG: " << rsgll << std::endl;

    std::cout << "splot:" << std::endl;
    std::cout << "Hgg: " << shggll << std::endl;
    std::cout << "RSG: " << srsgll << std::endl;

    delete thisShistHgg;
    delete thisShistRSG;
  }
  fits.getSimpleTotalBkgSubtraction("Hgg125");
  fits.getSimpleTotalBkgSubtraction("RSG125");
  RooDataHist *thisBkgHgg = (RooDataHist*)toyws->data(Form("Data_Hgg125_Combined_bkgSub_cosT"));
  RooDataHist *thisBkgRSG = (RooDataHist*)toyws->data(Form("Data_RSG125_Combined_bkgSub_cosT"));
  hggPdf = (RooHistPdf*)ws->pdf("Hgg125_FIT_cosTpdf");
  rsgPdf = (RooHistPdf*)ws->pdf("RSG125_FIT_cosTpdf");

  double totHggLL = computeLL(hggPdf,thisBkgHgg,cosT);
  //double totRSGLL = computeLL(rsgPdf,thisBkgRSG,cosT);
  double totRSGLL = computeLL(rsgPdf,thisBkgHgg,cosT);

  //use 2D fit to test compatibility
  RooRealVar *yield1D   = toyws->var("Data_Hgg125_FULLFIT_Nsig");
  RooRealVar *yieldHgg  = toyws->var("Data_m_Hgg125_c_Hgg125_FULL2DFIT_Nsig");
  RooRealVar *yieldRSG  = toyws->var("Data_m_Hgg125_c_RSG125_FULL2DFIT_Nsig"); // yield fitting hgg mass and rsg cosT distributions
  
  double errHgg = TMath::Sqrt(pow(yield1D->getError(),2) + pow(yieldHgg->getError(),2));
  double errRSG = TMath::Sqrt(pow(yield1D->getError(),2) + pow(yieldRSG->getError(),2));

  double hgg2Dll  = TMath::Log(TMath::Prob( pow(yield1D->getVal()-yieldHgg->getVal(),2)/pow(errHgg,2),1));
  double rsg2Dll  = TMath::Log(TMath::Prob( pow(yield1D->getVal()-yieldRSG->getVal(),2)/pow(errRSG,2),1));

  GenMinusFit->setVal( nGen-nFit);
  std::cout << "Gen - Fit: " << nGen - nFit <<std::endl;
  if(saveWorkspaces) toyWSs.AddLast(toyws);
  else delete toyws;

  //setup outputs:
  out[0] = 2*(rsgll-hggll);
  out[1] = 2*(totRSGLL-totHggLL);
  out[2] = 2*(srsgll-shggll);
  out[3] = 2*(rsg2Dll-hgg2Dll);

  return out;
}

void MakeSpinToy::runN(int N){
  RooMsgService::instance().setSilentMode(true);
  int Nout;
  for(int i=0;i<N;i++){
    double * hggOut = run1(Hgg125,Nout);
    double * rsgOut = run1(RSG125,Nout);

    S->setVal(hggOut[0]);
    S_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[1]);
    S_tot_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[2]);
    S_splot_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(hggOut[3]);
    S_2D_TruthHgg->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(rsgOut[0]);
    S_TruthRSG->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(rsgOut[1]);
    S_tot_TruthRSG->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(rsgOut[2]);
    S_splot_TruthRSG->add(RooArgSet(*S,*GenMinusFit));   

    S->setVal(rsgOut[3]);
    S_2D_TruthRSG->add(RooArgSet(*S,*GenMinusFit));   

    delete hggOut;
    delete rsgOut;
  }
}

void MakeSpinToy::save(TString outputFile){
  TFile *f = new TFile(outputFile,"RECREATE");
  makeForCombineTool("q",S_TruthHgg,S_TruthRSG)->Write();
  makeForCombineTool("qtot",S_tot_TruthHgg,S_tot_TruthRSG)->Write();
  makeForCombineTool("qsplot",S_splot_TruthHgg,S_splot_TruthRSG)->Write();
  makeForCombineTool("q2D",S_2D_TruthHgg,S_2D_TruthRSG)->Write();
  if(toyWSs.GetEntries()>0) toyWSs.Write();
  f->Close();
}

void MakeSpinToy::setUseR9(bool b){
  useR9 = b;
  catLabels.clear();
  if(useR9){
    for(int i=0;i<nCat;i++){
      catLabels.push_back( Form("EB_%d",i) );
      catLabels.push_back( Form("EE_%d",i) );
    }
  }else{
    for(int i=0;i<nCat;i++){ for(int j=0;j<nCat;j++){ 
	catLabels.push_back( Form("EB_%d_%d",i,j) );
	catLabels.push_back( Form("EE_%d_%d",i,j) );
      } }
  }
}
