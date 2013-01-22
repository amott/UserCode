#include "MakeSpinPlots.h"
#include "MakeSpinWorkspace.h"
#include "subtract.cc"


MakeSpinPlots::MakeSpinPlots(TString inputFileName, TString outTag): nCat(MakeSpinWorkspace::nCat){
  inputFile = new TFile(inputFileName);
  ws = (RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace");

  basePath = "figs/";
  outputTag=outTag;

  setStyle();

  if(ws->pdf("Data_Hgg125_FIT")==0) combinedFit=false;
  else combinedFit=true;
}

MakeSpinPlots::~MakeSpinPlots(){
  inputFile->Close();
}

void MakeSpinPlots::addMCType(TString s){
  mcNames.push_back(s);
}

void MakeSpinPlots::runAll(TString mcName){
  if(useR9){
    for(int i=0;i<nCat;i++){
      for(int eb=0;eb<2;eb++){
	TString tag = Form("%s_%d",(eb==0?"EB":"EE"),i);
	runAll(tag,mcName);
      }
    }
  }else{
      for(int i=0;i<nCat;i++){
	for(int j=0;j<nCat;j++){
	  for(int eb=0;eb<2;eb++){
	    TString tag = Form("%s_%d_%d",(eb==0?"EB":"EE"),i,j);
	    runAll(tag,mcName);
	  }
	}
      }
  }
  DrawSpinSubTotBackground(mcName,false);
  DrawSpinSubTotBackground(mcName,true);
}

void MakeSpinPlots::runAll(TString tag, TString mcName){
  getFitValues(tag,mcName);
  DrawBlindFit(tag,mcName);
  //if(!combinedFit) DrawBlindFit(tag,mcName);
  DrawFit(tag,mcName);
  if(!combinedFit) PlotSignalFits(tag,mcName);
  DrawSpinSubBackground(tag,mcName,false);
  DrawSpinSubBackground(tag,mcName,true);
}

void MakeSpinPlots::getFitValues(TString tag,TString mcName){
  std::cout << tag << std::endl;
  //Data_RSG125_FULLFIT_Nsig
  double sig  = ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))->getVal();
  double sige = ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))->getError();
  double fsig = ws->var(Form("Data_%s_FULLFIT_%s_fsig",mcName.Data(),tag.Data()))->getVal();
  std::cout << sig*fsig << " +- " << sige*fsig <<std::endl;
  double bkg  = ws->var(Form("Data_%s_FULLFIT_Nbkg",mcName.Data()))->getVal();
  double bkge = ws->var(Form("Data_%s_FULLFIT_Nbkg",mcName.Data()))->getError();
  double fbkg = ws->var(Form("Data_%s_FULLFIT_%s_fbkg",mcName.Data(),tag.Data()))->getVal();
  std::cout << bkg*fbkg << " +- " << bkge*fbkg <<std::endl;
  double f    = ws->var(Form("%s_FIT_%s_f",mcName.Data(),tag.Data()))->getVal();
  double fe   = ws->var(Form("%s_FIT_%s_f",mcName.Data(),tag.Data()))->getError();
  std::cout << f << " +- " << fe <<std::endl;
  double s1   = ws->var(Form("%s_FIT_%s_sigma1",mcName.Data(),tag.Data()))->getVal();
  double s2   = ws->var(Form("%s_FIT_%s_sigma2",mcName.Data(),tag.Data()))->getVal();

  double s1e  = ws->var(Form("%s_FIT_%s_sigma1",mcName.Data(),tag.Data()))->getError();
  double s2e  = ws->var(Form("%s_FIT_%s_sigma2",mcName.Data(),tag.Data()))->getError();

  std::cout << s1 << " +- " << s1e <<std::endl;  
  std::cout << s2 << " +- " << s2e <<std::endl;  
  double mean = ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data()))->getVal();
  double meane= ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data()))->getError();

  std::cout << mean << " +- " << meane <<std::endl;  

  double sigEff = s1*f + s2*(1-f); // s2 + f*(s1-s2)
  double s1ms2e = TMath::Sqrt(s1e*s1e+s2e*s2e);
  double e2 = f*(s1-s2)*TMath::Sqrt(fe*fe/f/f+s1ms2e*s1ms2e/(s1-s2)/(s1-s2));
  double sigEffE = TMath::Sqrt(s2e*s2e+e2*e2);
  std::cout << sigEff << " +- " << sigEffE <<std::endl;  
  //Data_BKGFIT_EB_0_bkgModel
  RooAbsPdf * bkgPdf = ws->pdf(Form("Data_BKGFIT_%s_bkgModel",tag.Data()));
  std::cout << bkgPdf <<std::endl;
  RooRealVar range("range","",mean-sigEff,mean+sigEff);
  RooRealVar all("all","",110,180);
  double BkgInRange  = bkgPdf->createIntegral(range)->getVal()/bkgPdf->createIntegral(all)->getVal()*bkg*fbkg;
  double BkgInRangeE = bkgPdf->createIntegral(range)->getVal()/bkgPdf->createIntegral(all)->getVal()*bkge*fbkg;
  

  tPair lbl(mcName,tag);

  nSignal[lbl]      = dPair(sig*fsig,sige*fsig);
  nBackground[lbl]  = dPair(bkg*fbkg,bkge*fbkg);
  fitMean[lbl]      = dPair(mean,meane);
  fitSigEff[lbl]    = dPair(sigEff,sigEffE);
  fitBkg1Sigma[lbl] = dPair(BkgInRange,BkgInRangeE);

}


void MakeSpinPlots::DrawBlindFit(TString tag, TString mcName){
  TCanvas *cv = new TCanvas(Form("%s_%s",mcName.Data(),tag.Data()));

  
  RooRealVar* mass = ws->var("mass");
  RooPlot* frame  = mass->frame(110,170,40);
  double Nb = ws->var(Form("Data_BKGFIT_%s_Nbkg",tag.Data()))->getVal();
  cout << Nb << endl;
  RooDataSet *blind = (RooDataSet*)ws->data(Form("Data_%s",tag.Data()))->reduce("(mass>110 && mass<119) || (mass>135.5 && mass<170)");
  blind->plotOn(frame);

  tPair lbl(mcName,tag);
  double nBkg = ws->data(Form("Data_%s",tag.Data()))->sumEntries();

  ws->pdf( Form("Data_BKGFIT_%s_bkgModel",tag.Data()) )->plotOn(frame,RooFit::Range("all"),//RooFit::Normalization(nBkg/blind->sumEntries()),
									       RooFit::LineColor(kRed));

  //TLatex *prelim = new TLatex(250,x->GetXmax()-40.,"CMS Preliminary");
  TLatex *prelim = new TLatex(0.12,0.96,"CMS Preliminary");
  TLatex *lum = new TLatex(0.7,0.96,Form("#sqrt{s}=8 TeV  L = %0.1f fb^{-1}",lumi));
  prelim->SetNDC();
  lum->SetNDC();
  prelim->SetTextSize(0.045);
  prelim->SetTextColor(kBlack);
  lum->SetTextSize(0.045);
  lum->SetTextColor(kBlack);

  TLatex *owner = new TLatex(0.6,0.88,"Alex Mott (Nov. 13, 2012)");
  owner->SetNDC();
  owner->SetTextSize(0.045);
  owner->SetTextColor(kBlack);

  TLatex *Nbkg = new TLatex(0.7,0.8,Form("N_{bkg}= %0.1f #pm %0.1f",nBackground[lbl].first,nBackground[lbl].second));
  Nbkg->SetNDC();
  Nbkg->SetTextSize(0.045);

  TLatex *sig = new TLatex(0.7,0.72,Form("#sigma_{eff} = %0.1f #pm %0.2f",fitSigEff[lbl].first,fitSigEff[lbl].second));
  sig->SetNDC();
  sig->SetTextSize(0.045);

  TLatex *expBkg = new TLatex(0.7,0.64,Form("B @ 125 = %0.1f",fitBkg1Sigma[lbl].first));
  expBkg->SetNDC();
  expBkg->SetTextSize(0.045);


  frame->addObject(prelim);
  frame->addObject(lum);
  //frame->addObject(owner);
  frame->addObject(Nbkg);
  frame->addObject(sig);
  frame->addObject(expBkg);
  frame->Draw();
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s_BLIND.png",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s_BLIND.pdf",outputTag.Data(),mcName.Data(),tag.Data()) );
  delete cv;
}

void MakeSpinPlots::DrawFit(TString tag, TString mcName){
  TCanvas *cv = new TCanvas(Form("%s_%s",mcName.Data(),tag.Data()));
  
  RooRealVar* mass = ws->var("mass");
  RooPlot* frame  = mass->frame(110,170,40);

  tPair lbl(mcName,tag);

  double Ns = ((RooFormulaVar*)ws->obj( Form("Data_%s_FULLFIT_%s_Nsig",mcName.Data(),tag.Data()) ))->getVal();
  double Nb = ((RooFormulaVar*)ws->obj( Form("Data_%s_FULLFIT_%s_Nbkg",mcName.Data(),tag.Data()) ))->getVal();

  double Nblind = ws->data(Form("Data_%s",tag.Data()))->reduce("(mass>100 && mass<119) || (mass>135.5 && mass<170)")->sumEntries();
  double Ntot   = ws->data(Form("Data_%s",tag.Data()))->sumEntries();

  RooFitResult* fitres = (RooFitResult*)ws->obj(Form("Data_%s_FULLFIT_fitResult",mcName.Data())); 
  std::cout << fitres << std::endl;
  ws->data(Form("Data_%s",tag.Data()))->plotOn(frame,RooFit::LineColor(kWhite),RooFit::MarkerColor(kWhite));
  //Data_Hgg125_FULLFIT_EB_0
  ws->pdf(Form("Data_%s_FULLFIT_%s",mcName.Data(),tag.Data()))->plotOn(frame, RooFit::FillColor(kGreen),RooFit::VisualizeError(*fitres,2.0));
  ws->pdf(Form("Data_%s_FULLFIT_%s",mcName.Data(),tag.Data()))->plotOn(frame, RooFit::FillColor(kYellow),RooFit::VisualizeError(*fitres,1.0));
  ws->pdf(Form("Data_%s_FULLFIT_%s",mcName.Data(),tag.Data()))->plotOn(frame, RooFit::LineColor(kRed));
  std::cout << "1" << std::endl;
  ws->pdf(Form("Data_BKGFIT_%s_bkgModel",tag.Data()))->plotOn(frame, RooFit::Normalization(Nb/(Nb+Ns)*Nblind/Ntot),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
  std::cout << "2" << std::endl;

  ws->data(Form("Data_%s",tag.Data()))->plotOn(frame);
  frame->Draw();

  //TLatex *prelim = new TLatex(250,x->GetXmax()-40.,"CMS Preliminary");
  TLatex *prelim = new TLatex(0.12,0.96,"CMS Preliminary");
  TLatex *lum = new TLatex(0.7,0.96,Form("#sqrt{s}=8 TeV  L = %0.1f fb^{-1}",lumi));
  prelim->SetNDC();
  lum->SetNDC();
  prelim->SetTextSize(0.045);
  prelim->SetTextColor(kBlack);
  lum->SetTextSize(0.045);
  lum->SetTextColor(kBlack);

  TLatex *owner = new TLatex(0.6,0.88,"Alex Mott (Nov. 13, 2012)");
  owner->SetNDC();
  owner->SetTextSize(0.045);
  owner->SetTextColor(kBlack);

  TLatex *mu = new TLatex(0.7,0.8,Form("#mu = %0.1f #pm %0.2f", fitMean[lbl].first,fitMean[lbl].second));
  mu->SetNDC();
  mu->SetTextSize(0.045);

  TLatex *sig = new TLatex(0.7,0.72,Form("#sigma_{eff} = %0.1f #pm %0.2f", fitSigEff[lbl].first,fitSigEff[lbl].second));
  sig->SetNDC();
  sig->SetTextSize(0.045);

  TLatex *Nsig = new TLatex(0.7,0.64,Form("N_{sig}= %0.1f #pm %0.1f",nSignal[lbl].first,nSignal[lbl].second));
  Nsig->SetNDC();
  Nsig->SetTextSize(0.045);


  frame->addObject(prelim);
  frame->addObject(lum);
  frame->addObject(owner);
  frame->addObject(mu);
  frame->addObject(sig);
  frame->addObject(Nsig);
  frame->Draw();
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s.png",outputTag.Data(),mcName.Data(),tag.Data()) );
  cv->SaveAs( basePath+Form("/mgg-%s-%s-%s.pdf",outputTag.Data(),mcName.Data(),tag.Data()) );
  delete cv;
}
/*
void printYields(TString wsFile,TString mcName){
  TFile *f = TFile::Open(wsFile);
  RooWorkspace *ws = (RooWorkspace*)f->Get("cms_hgg_spin_workspace");

  cout << "\\begin{tabular}{c|c|c||c|c|c|c|c}" << endl;
  cout << "Region & $\\gamma_1$ & $\\gamma_2$ & $N_{bkg}$ & $N_{sig}$ & mean & $\\sigma_{eff}$ & Signal Efficiency \\\\" <<endl;
  cout << "\\hline" <<endl;
  float N=0,Esq=0;
  double selEB=0;

  double totEB  = ws->var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws->var("Hgg125_EE_totalEvents")->getVal();

  for(int i=0;i<nCat;i++){
    for(int j=0;j<nCat;j++){
      double thisEB  = ws->data(Form("%s_EB_%d_%d",mcName.Data(),i,j))->sumEntries();
      selEB+=thisEB;
      double v[8];
      getSignalFitVals(ws,mcName,Form("EB_%d_%d",i,j),v);
      N+=v[0];
      Esq+=v[1]*v[1];
      cout  << "EB & " << (i==0?"pass":"fail") << " & " << (j==0?"pass":"fail") << " & "
	    << v[2]  << " $\\pm$ " << v[3]  << " & "
	    << v[0]  << " $\\pm$ " << v[1]  << " & "
	    << v[4] << " $\\pm$ " << v[5] << " & "
	    << v[6] <<  " $\\pm$ " << v[7] << " & "
	    << thisEB/(totEB+totEE) << " \\\\ " << endl;
	
    }
  }
  cout << "\\hline" << endl;
  double selEE = 0;
  for(int i=0;i<nCat;i++){
    for(int j=0;j<nCat;j++){
      double thisEE  = ws->data(Form("%s_EE_%d_%d",mcName.Data(),i,j))->sumEntries();
      selEE+=thisEE;
      double v[8];
      getSignalFitVals(ws,mcName,Form("EE_%d_%d",i,j),v);
      N+=v[0];
      Esq+=v[1]*v[1];
      cout  << "EE & " << (i==0?"pass":"fail") << " & " << (j==0?"pass":"fail") << " & "
	    << v[2]  << " $\\pm$ " << v[3]  << " & "
	    << v[0]  << " $\\pm$ " << v[1]  << " & "
	    << v[4] << " $\\pm$ " << v[5] << " & "
	    << v[6] <<  " $\\pm$ " << v[7] << " & " 
	    << thisEE/(totEB+totEE) <<" \\\\\ " << endl;
	
    }
  }

  cout << "\\end{tabular}" << endl;

  cout << "Total Nsig:  " << N << " +- " << TMath::Sqrt(Esq) << endl;
  cout << "Total Eff: \n\tEB: " << selEB/totEB << "\n\tEE: " << selEE/totEE << "\n\tTot: " << (selEB+selEE)/(totEB+totEE) << endl;

}
*/

void MakeSpinPlots::DrawSpinBackground(TString tag, TString mcName,bool signal){
  double totEB  = ws->var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws->var("Hgg125_EE_totalEvents")->getVal();

  TCanvas cv;
  double thisN  = ws->data(Form("%s_%s",mcName.Data(),tag.Data()))->sumEntries();
  float norm = 607*lumi/12.*thisN/(totEB+totEE);
  if(signal) norm = ws->var(Form("Data_%s_FIT_%s_Nsig",mcName.Data(),tag.Data()))->getVal();
  RooPlot *frame = ws->var("cosT")->frame(-1,1,3);

  RooDataSet* bkgWeight = (RooDataSet*)ws->data(Form("%s_bkgWeight_%s",mcName.Data(),tag.Data()));
  RooDataSet* tmp = (RooDataSet*)ws->data(Form("Data_%s",tag.Data()))->reduce("(mass>115 && mass<120) || (mass>130 && mass<135)");
  tmp->plotOn(frame,RooFit::Rescale(norm/tmp->sumEntries()));

  ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",tag.Data()))->plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(norm/tmp->sumEntries()));
  ws->pdf(Form("RSG125_FIT_%s_cosTpdf",tag.Data()))->plotOn(frame,RooFit::LineColor(kGreen),RooFit::Normalization(norm/tmp->sumEntries()));

  bkgWeight->plotOn(frame,RooFit::Rescale(norm/bkgWeight->sumEntries()),RooFit::MarkerColor(kBlue) );  
  if(signal){
    ws->data(Form("%s_sigWeight_%s",mcName.Data(),tag.Data()))->plotOn(frame,RooFit::MarkerStyle(4));
  }
  
  frame->SetMaximum(frame->GetMaximum()*(signal?0.8:0.4)*norm/tmp->sumEntries());
  frame->SetMinimum(-1*frame->GetMaximum());
  TLegend l(0.6,0.2,0.95,0.45);
  l.SetFillColor(0);
  l.SetBorderSize(0);
  l.SetHeader(tag);
  l.AddEntry(frame->getObject(0),"Data m#in [115,120]#cup[130,135]","p");
  l.AddEntry(frame->getObject(1),"SM Higgs","l");
  l.AddEntry(frame->getObject(2),"RS Graviton","l");
  l.AddEntry(frame->getObject(3),"background weighted Data","p");
  if(signal) l.AddEntry(frame->getObject(4),"signal weighted Data","p");
  
  frame->Draw();
  l.Draw("SAME");
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_%s%s_%s_%s.png",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_%s%s_%s_%s.pdf",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
}

void MakeSpinPlots::DrawSpinSubBackground(TString tag, TString mcName,bool signal){
  double totEB  = ws->var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws->var("Hgg125_EE_totalEvents")->getVal();

  TCanvas cv;
  double thisN  = ws->data(Form("%s_%s",mcName.Data(),tag.Data()))->sumEntries();
  float norm = 607*lumi/12.*thisN/(totEB+totEE);
  tPair lbl(mcName,tag);


  if(signal) norm = nSignal[lbl].first;   //((RooFormulaVar*)ws->obj(Form("Data_%s_FULLFIT_%s_Nsig",mcName.Data(),tag.Data())) )->getVal();
  RooPlot *frame = ws->var("cosT")->frame(-1,1,10);

  RooDataSet* tmp = (RooDataSet*)ws->data(Form("Data_%s",tag.Data()))->reduce("(mass>115 && mass<120) || (mass>130 && mass<135)");
  tmp->plotOn(frame,RooFit::Rescale(norm/tmp->sumEntries()));

  ws->pdf(Form("Hgg125_FIT_%s_cosTpdf",tag.Data()))->plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(norm/tmp->sumEntries()));
  ws->pdf(Form("RSG125_FIT_%s_cosTpdf",tag.Data()))->plotOn(frame,RooFit::LineColor(kGreen),RooFit::Normalization(norm/tmp->sumEntries()));
  if(signal){
    RooDataHist *h = (RooDataHist*)ws->data( Form("Data_%s_%s_bkgSub_cosT",mcName.Data(),tag.Data()) );
    h->plotOn(frame,RooFit::MarkerStyle(4));
    std::cout << "Nsig: " << h->sumEntries() << std::endl;
  }
  
  frame->SetMaximum(frame->GetMaximum()*(signal?10.:1.2)*norm/tmp->sumEntries());
  frame->SetMinimum(-1*frame->GetMaximum());
  TLegend l(0.6,0.2,0.95,0.45);
  l.SetFillColor(0);
  l.SetBorderSize(0);
  l.SetHeader(tag);
  l.AddEntry(frame->getObject(0),"Data m#in [115,120]#cup[130,135]","p");
  l.AddEntry(frame->getObject(1),"SM Higgs","l");
  l.AddEntry(frame->getObject(2),"RS Graviton","l");
  if(signal) l.AddEntry(frame->getObject(3),"bkg-subtracted Data","p");
  
  frame->Draw();
  l.Draw("SAME");
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s_%s.png",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s_%s.pdf",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data(),tag.Data()) );
}

void MakeSpinPlots::DrawSpinSubTotBackground(TString mcName,bool signal){
  double totEB  = ws->var("Hgg125_EB_totalEvents")->getVal();
  double totEE  = ws->var("Hgg125_EE_totalEvents")->getVal();

  TCanvas cv;
  double thisN  = ws->data(Form("%s_Combined",mcName.Data()))->sumEntries();
  float norm = 607*lumi/12.*thisN/(totEB+totEE);


  if(signal) norm = ws->var(Form("Data_%s_FULLFIT_Nsig",mcName.Data()))->getVal();
  RooPlot *frame = ws->var("cosT")->frame(-1,1,10);

  RooDataSet* tmp = (RooDataSet*)ws->data(Form("Data_Combined"))->reduce("(mass>115 && mass<120) || (mass>130 && mass<135)");
  tmp->plotOn(frame,RooFit::Rescale(norm/tmp->sumEntries()));

  ws->pdf(Form("Hgg125_FIT_cosTpdf"))->plotOn(frame,RooFit::LineColor(kRed),RooFit::Normalization(norm/tmp->sumEntries()));
  ws->pdf(Form("RSG125_FIT_cosTpdf"))->plotOn(frame,RooFit::LineColor(kGreen),RooFit::Normalization(norm/tmp->sumEntries()));
  if(signal){
    RooDataHist *h = (RooDataHist*)ws->data( Form("Data_%s_Combined_bkgSub_cosT",mcName.Data()) );
    h->plotOn(frame,RooFit::MarkerStyle(4));
    std::cout << "Nsig: " << h->sumEntries() << std::endl;
  }
  
  frame->SetMaximum(frame->GetMaximum()*(signal?2.:1.2)*norm/tmp->sumEntries());
  frame->SetMinimum(-1*frame->GetMaximum());
  TLegend l(0.6,0.2,0.95,0.45);
  l.SetFillColor(0);
  l.SetBorderSize(0);
  l.SetHeader("Combined");
  l.AddEntry(frame->getObject(0),"Data m#in [115,120]#cup[130,135]","p");
  l.AddEntry(frame->getObject(1),"SM Higgs","l");
  l.AddEntry(frame->getObject(2),"RS Graviton","l");
  if(signal) l.AddEntry(frame->getObject(3),"bkg-subtracted Data","p");
  
  frame->Draw();
  l.Draw("SAME");
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s.png",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data()) );
  cv.SaveAs( basePath+Form("/cosThetaPlots/CosThetaDist_SimpleSub_%s%s_%s.pdf",outputTag.Data(),(signal ? "":"_BLIND"),mcName.Data()) );
}


void MakeSpinPlots::PlotSignalFits(TString tag, TString mcName){
  TCanvas cv;
  float mean = ws->var(Form("%s_FIT_%s_mean",mcName.Data(),tag.Data()))->getVal();
  RooPlot *frame = ws->var("mass")->frame(mean-10,mean+10,40);
  ws->data(Form("%s_%s",mcName.Data(),tag.Data()))->plotOn(frame); //data
  RooFitResult *res = (RooFitResult*)ws->obj(Form("%s_FIT_%s_fitResult",mcName.Data(),tag.Data()));
  RooAbsPdf * pdf = ws->pdf(Form("%s_FIT_%s",mcName.Data(),tag.Data())); //signal model
  std::cout << pdf << "\t" << res << std::endl;
  pdf->plotOn(frame,RooFit::FillColor(kGreen),RooFit::VisualizeError(*res,2.0));
  pdf->plotOn(frame,RooFit::FillColor(kYellow),RooFit::VisualizeError(*res,1.0));
  pdf->plotOn(frame,RooFit::LineColor(kRed));
  ws->data(Form("%s_%s",mcName.Data(),tag.Data()))->plotOn(frame); //data
  
  tPair lbl(mcName,tag);

  TLatex *prelim = new TLatex(0.57,0.9,"CMS Preliminary Simulation");
  TLatex *sigL  = new TLatex(0.67,0.6,Form("#sigma_{eff} = %0.1f #pm %0.2f GeV",fitSigEff[lbl].first,fitSigEff[lbl].second));
  prelim->SetNDC();
  sigL->SetNDC();
  prelim->SetTextSize(0.05);
  sigL->SetTextSize(0.05);
  
  frame->addObject(prelim);
  frame->addObject(sigL);
  frame->Draw();
  cv.SaveAs(basePath+Form("/signalModels/sig_%s_%s.png",outputTag.Data(),tag.Data()));
  cv.SaveAs(basePath+Form("/signalModels/sig_%s_%s.pdf",outputTag.Data(),tag.Data()));

}

void MakeSpinPlots::setStyle(){

  vecbosStyle = new TStyle("vecbosStyle","Style for P-TDR");

  // For the canvas:
  vecbosStyle->SetCanvasBorderMode(0);
  vecbosStyle->SetCanvasColor(kWhite);
  vecbosStyle->SetCanvasDefH(600); //Height of canvas
  vecbosStyle->SetCanvasDefW(900); //Width of canvas
  vecbosStyle->SetCanvasDefX(0);   //POsition on screen
  vecbosStyle->SetCanvasDefY(0);

  // For the Pad:
  vecbosStyle->SetPadBorderMode(0);
  // vecbosStyle->SetPadBorderSize(Width_t size = 1);
  vecbosStyle->SetPadColor(kWhite);
  vecbosStyle->SetPadGridX(true);
  vecbosStyle->SetPadGridY(true);
  vecbosStyle->SetGridColor(0);
  vecbosStyle->SetGridStyle(3);
  vecbosStyle->SetGridWidth(1);

  // For the frame:
  vecbosStyle->SetFrameBorderMode(0);
  vecbosStyle->SetFrameBorderSize(1);
  vecbosStyle->SetFrameFillColor(0);
  vecbosStyle->SetFrameFillStyle(0);
  vecbosStyle->SetFrameLineColor(1);
  vecbosStyle->SetFrameLineStyle(1);
  vecbosStyle->SetFrameLineWidth(1);

  // set the paper & margin sizes
  vecbosStyle->SetPaperSize(20,26);
  vecbosStyle->SetPadTopMargin(0.05);
  vecbosStyle->SetPadRightMargin(0.05);
  vecbosStyle->SetPadBottomMargin(0.16);
  vecbosStyle->SetPadLeftMargin(0.12);

  // use large Times-Roman fonts
  vecbosStyle->SetTitleFont(132,"xyz");  // set the all 3 axes title font
  vecbosStyle->SetTitleFont(132," ");    // set the pad title font
  vecbosStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  vecbosStyle->SetTitleSize(0.06," ");   // set the pad title size
  vecbosStyle->SetLabelFont(132,"xyz");
  vecbosStyle->SetLabelSize(0.05,"xyz");
  vecbosStyle->SetLabelColor(1,"xyz");
  vecbosStyle->SetTextFont(132);
  vecbosStyle->SetTextSize(0.08);
  vecbosStyle->SetStatFont(132);

  vecbosStyle->SetTitleOffset(0.9,"Y");
  // use bold lines and markers
  vecbosStyle->SetMarkerStyle(8);
  vecbosStyle->SetMarkerSize(1.2);
  vecbosStyle->SetHistLineWidth(1.85);
  vecbosStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //..Get rid of X error bars
  //vecbosStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  vecbosStyle->SetOptTitle(0);
  vecbosStyle->SetOptStat(0);
  vecbosStyle->SetOptFit(11111111);

  // put tick marks on top and RHS of plots
  vecbosStyle->SetPadTickX(1);
  vecbosStyle->SetPadTickY(1);

  // set a decent palette
  vecbosStyle->SetPalette(1);

  vecbosStyle->cd();

  gROOT->SetStyle("vecbosStyle");
  gROOT->ForceStyle();

}
