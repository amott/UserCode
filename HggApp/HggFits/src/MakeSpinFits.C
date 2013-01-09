#include "MakeSpinFits.h"
using namespace std;
MakeSpinFits::MakeSpinFits(TString inputFileName, TString outputFileName):
  nCat(MakeSpinWorkspace::nCat),
  addSWeight(true),
  useR9(false)
{
  if(inputFileName != ""){
    inputFile = new TFile(inputFileName);
    ws = ((RooWorkspace*)inputFile->Get("cms_hgg_spin_workspace"));
  }
  if(outputFileName != ""){
    outputFile = new TFile(outputFileName,"RECREATE");
  }

  setUseR9(false);
}

MakeSpinFits::~MakeSpinFits(){

}

void MakeSpinFits::setUseR9(bool b){
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

void MakeSpinFits::MakeSignalFit(TString tag, TString mcName){
  if(ws==0) return;
  TString inputTag = Form("%s_%s",mcName.Data(),tag.Data());
  TString outputTag = Form("%s_FIT_%s",mcName.Data(),tag.Data());

  RooRealVar mass = *(ws->var("mass"));
  RooRealVar cosT = *(ws->var("cosT"));
  cosT.setBins(10);

  RooRealVar mean(Form("%s_mean",outputTag.Data()),Form("%s_mean",outputTag.Data()),125,100,180);
  RooRealVar sig1(Form("%s_sigma1",outputTag.Data()),Form("%s_sigma1",outputTag.Data()),1,0.2,10);
  RooRealVar sig2(Form("%s_sigma2",outputTag.Data()),Form("%s_sigma2",outputTag.Data()),1,0.2,10);
  RooRealVar f(Form("%s_f",outputTag.Data()),Form("%s_f",outputTag.Data()),0.1,0,1);
  RooGaussian g1(Form("%s_g1",outputTag.Data()),Form("%s_g1",outputTag.Data()),mass,mean,sig1);
  RooGaussian g2(Form("%s_g2",outputTag.Data()),Form("%s_g2",outputTag.Data()),mass,mean,sig2);
  RooAddPdf SignalModel(outputTag.Data(),"Signal Model",RooArgList(g1,g2),f);

  RooRealVar meanCB(Form("%s_meanCB",outputTag.Data()),Form("%s_meanCB",outputTag.Data()),125,100,180);
  RooRealVar sig1CB(Form("%s_sigma1CB",outputTag.Data()),Form("%s_sigma1CB",outputTag.Data()),1,0.2,10);
  RooRealVar sig2CB(Form("%s_sigma2CB",outputTag.Data()),Form("%s_sigma2CB",outputTag.Data()),1,0.2,10);
  RooRealVar alpha1CB(Form("%s_alpha1CB",outputTag.Data()),Form("%s_alpha1CB",outputTag.Data()),1,0.,1000);
  RooRealVar alpha2CB(Form("%s_alpha2CB",outputTag.Data()),Form("%s_alpha2CB",outputTag.Data()),1,0.,1000);
  RooRealVar n1CB(Form("%s_n1CB",outputTag.Data()),Form("%s_n1CB",outputTag.Data()),1,0.,1000);
  RooRealVar n2CB(Form("%s_n2CB",outputTag.Data()),Form("%s_n2CB",outputTag.Data()),1,0.,1000);
  RooRealVar fCB(Form("%s_fCB",outputTag.Data()),Form("%s_fCB",outputTag.Data()),0.1,0,1);
  
  RooCBShape cb1(Form("%s_cb1",outputTag.Data()),Form("%s_cb1",outputTag.Data()),mass,meanCB,sig1CB,alpha1CB,n1CB);
  RooCBShape cb2(Form("%s_cb2",outputTag.Data()),Form("%s_cb2",outputTag.Data()),mass,meanCB,sig2CB,alpha2CB,n2CB);

  RooAddPdf SignalModelCB(outputTag.Data(),"Signal Model",RooArgList(cb1,cb2),fCB);
  

  RooDataSet *ds = (RooDataSet*)ws->data(inputTag);

  //RooKeysPdf cosTkde(Form("%s_cosTpdf",outputTag.Data()),"KDE for cos(theta) dist",cosT,*ds);
  RooDataHist hist(Form("%s_cosThist",outputTag.Data()),"Data Hist for cos(theta)",RooArgSet(cosT),*ds);
  RooHistPdf cosTkde(Form("%s_cosTpdf",outputTag.Data()),"Hist PDF for cos(theta)",RooArgSet(cosT),hist);

  SignalModel.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(0));
  RooFitResult *res = SignalModel.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2));
  SignalModelCB.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(0));
  RooFitResult *resCB = SignalModelCB.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2));
  std::cout << res <<std::endl;
  res->SetName(Form("%s_FitResult",outputTag.Data()));
  resCB->SetName(Form("%s_FitResultCB",outputTag.Data()));

  cosTkde.fitTo(*ds);
  ws->import(SignalModel);
  ws->import(SignalModelCB);
  ws->import(cosTkde);
  ws->import(*res);
}

void MakeSpinFits::AddSWeight(TString mcName,TString catTag,TString inputTag){
  if(ws==0) return;
  if(inputTag=="") inputTag=catTag;
  RooRealVar mean = *(ws->var(Form("%s_FIT_%s_mean",mcName.Data(),inputTag.Data())));
  mean.setConstant(kTRUE);

  RooDataSet *ds  = (RooDataSet*)ws->data(Form("Data_%s",inputTag.Data()));
  RooAbsPdf  *pdf = ws->pdf(Form("Data_%s_FIT_%s_fitModel",mcName.Data(),inputTag.Data()));
  RooRealVar *Ns = ws->var(Form("Data_%s_FIT_%s_Nsig",mcName.Data(),inputTag.Data()));
  RooRealVar *Nb = ws->var(Form("Data_%s_FIT_%s_Nbkg",mcName.Data(),inputTag.Data()));

  RooStats::SPlot* sData = new RooStats::SPlot(Form("%s_sData_%s",mcName.Data(),catTag.Data()),"",
					       *ds,pdf,RooArgList(*Ns,*Nb));


  RooArgList sweights = sData->GetSWeightVars();

  RooArgSet dsVars;
  dsVars.add(*(RooRealVar*)(sweights.at(0))); // add the sweights
  dsVars.add(*(RooRealVar*)(sweights.at(1)));
  
  //add all variables from the original dataset
  const RooArgSet *orig = ds->get(0);
  RooLinkedListIter it = orig->iterator();
  while(it.Next()){
    std::cout << (*it)->GetName() << std::endl;
    if(ws->var( (*it)->GetName()) ==0 ) continue;
    dsVars.add(*(ws->var( (*it)->GetName()) ));
  }
  //make weighted datasets
  RooDataSet sig(Form("%s_sigWeight_%s",mcName.Data(),catTag.Data()),"",ds,dsVars,0,sweights.at(0)->GetName());
  RooDataSet bkg(Form("%s_bkgWeight_%s",mcName.Data(),catTag.Data()),"",ds,dsVars,0,sweights.at(1)->GetName());
  

  ws->import(*ds,RooFit::Rename(Form("%s_sData_%s",mcName.Data(),catTag.Data())));
  ws->import(sig);
  ws->import(bkg);

}


void MakeSpinFits::MakeBackgroundFitCosTBin(TString mcName,TString catTag,float minCosT,float maxCosT){
  if(ws==0) return;
  RooRealVar *mass = ws->var("mass");

  const char *ctTag = Form("%0.1f_cosT_%0.1f",minCosT,maxCosT);

  RooRealVar* mean = ws->var(Form("Data_%s_FIT_%s_mean",mcName.Data(), catTag.Data()));
  
  RooDataSet *dR = (RooDataSet*)ws->data(Form("Data_%s",catTag.Data()))->reduce( Form("%f <= cosT && %f > cosT",minCosT,maxCosT) );
  ws->import(*dR,RooFit::Rename(Form("Data_%s_%s",ctTag,catTag.Data())));


  MakeBackgroundFit( mcName,Form("%s_%s",ctTag,catTag.Data()), mean->getVal(),mean->getError()/3.,true,catTag);
}

void MakeSpinFits::MakeCombinedBackgroundFit(TString mcName,float initMass,float range,TString inputTag){
  //the mean and relative signal fraction are constrained between the categories

  RooRealVar mass = *(ws->var("mass"));
  RooCategory *cat = ((RooCategory*)ws->obj("evtcat"));
  cout << "1" <<endl;

  std::vector<TString>::const_iterator catIt = catLabels.begin();

  RooArgSet allConstraints;

  //common variables
  RooRealVar mean(Form("Data_%s_COMBFIT_mean",mcName.Data()),"mean",initMass,initMass-range,initMass+range);
  RooRealVar Nsig(Form("Data_%s_COMBFIT_Nsig",mcName.Data()),"N signal Events",100,0,100000);

  RooSimultaneous combinedFit(Form("Data_%s_COMBFIT",mcName.Data()),"",*cat);
  cout << "2" <<endl;
  
  for(; catIt != catLabels.end(); catIt++){
    cout << *catIt <<endl;
    if(inputTag=="") inputTag=*catIt;
    RooRealVar sig1Base = *(ws->var(Form("%s_FIT_%s_sigma1",mcName.Data(),inputTag.Data())));
    RooRealVar sig2Base = *(ws->var(Form("%s_FIT_%s_sigma2",mcName.Data(),inputTag.Data())));
    RooRealVar fBase    = *(ws->var(Form("%s_FIT_%s_f",mcName.Data(),inputTag.Data())));
    cout << "\tA" <<endl;
  
    //fit constraints
    RooRealVar *sig1Width = new RooRealVar("sig1Wdith","",sig1Base.getError()/8);
    RooRealVar *sig2Width = new RooRealVar("sig2Wdith","",sig2Base.getError()/8);
    RooRealVar *fWidth    = new RooRealVar("fWdith","",fBase.getError()/8);
    
    RooRealVar *sig1 = new RooRealVar(Form("Data_%s_COMBFIT_%s_sigma1",mcName.Data(), catIt->Data()),"sig1",sig1Base.getVal(),sig1Base.getVal()-4*sig1Base.getError(),sig1Base.getVal()+4*sig1Base.getError());
    RooRealVar *sig2 = new RooRealVar(Form("Data_%s_COMBFIT_%s_sigma2",mcName.Data(), catIt->Data()),"sig2",sig2Base.getVal(),sig2Base.getVal()-4*sig2Base.getError(),sig2Base.getVal()+4*sig2Base.getError());
    RooRealVar *f    = new RooRealVar(Form("Data_%s_COMBFIT_%s_f",mcName.Data(), catIt->Data()),"f",fBase.getVal(),fBase.getVal()-4*fBase.getError(),fBase.getVal()+4*fBase.getError());
  
    RooGaussian *sig1Const = new RooGaussian(Form("sig1Constraint_%s",catIt->Data()),"",*sig1,sig1Base,*sig1Width);
    RooGaussian *sig2Const = new RooGaussian(Form("sig2Constraint_%s",catIt->Data()),"",*sig2,sig2Base,*sig2Width);
    RooGaussian *fConst    = new RooGaussian(Form("fConstraint_%s",catIt->Data()),"",*f,fBase,*fWidth);
    cout << "\tB" <<endl;
    allConstraints.add(*sig1Const);
    allConstraints.add(*sig2Const);
    allConstraints.add(*fConst);
    cout << "\tC" <<endl;
    //sig1.setConstant(kTRUE);
    //sig2.setConstant(kTRUE);
    //f.setConstant(kTRUE);

    //build the signal model
    RooGaussian *g1 = new RooGaussian(Form("Data_%s_COMBFIT_%s_g1",mcName.Data(),catIt->Data()),"g1",mass,mean,*sig1);
    RooGaussian *g2 = new RooGaussian(Form("Data_%s_COMBFIT_%s_g2",mcName.Data(),catIt->Data()),"g2",mass,mean,*sig2);
    RooAddPdf *SignalModel = new RooAddPdf(Form("Data_%s_COMBFIT_%s_signalModel",mcName.Data(),catIt->Data()),
					   "Signal Model",RooArgList(*g1,*g2),*f);
    cout << "\tD" <<endl;
    //build the background model
    RooRealVar *alpha1 = new RooRealVar(Form("Data_%s_COMBFIT_%s_alpha1",mcName.Data(),catIt->Data()),"alpha1",-0.1,-1.,0.);
    RooRealVar *alpha2 = new RooRealVar(Form("Data_%s_COMBFIT_%s_alpha2",mcName.Data(),catIt->Data()),"alpha2",-0.1,-1.,0.);
    RooRealVar *f_bkg  = new RooRealVar( Form("Data_%s_COMBFIT_%s_fbkg",mcName.Data(),catIt->Data()),"f_bkg",0.1,0,1);
    RooExponential *exp1 = new RooExponential(Form("Data_%s_COMBFIT_%s_exp1",mcName.Data(),catIt->Data()),"exp1",mass,*alpha1);
    RooExponential *exp2 = new RooExponential(Form("Data_%s_COMBFIT_%s_exp2",mcName.Data(),catIt->Data()),"exp2",mass,*alpha2);
    
    RooAddPdf *BkgModel = new RooAddPdf(Form("Data_%s_COMBFIT_%s_bkgModel",mcName.Data(),catIt->Data()),"Background Model",
					RooArgList(*exp1,*exp2),*f_bkg);
    cout << "\tE" <<endl;
    RooRealVar *Nbkg = new RooRealVar(Form("Data_%s_COMBFIT_%s_Nbkg",mcName.Data(),catIt->Data()),"N background Events",100,0,1e+09);
    //fraction of the signal in this category
    double totalEvents  = ws->var("Hgg125_EB_totalEvents")->getVal() + ws->var("Hgg125_EE_totalEvents")->getVal();
    double thisCat =  ws->data(Form("Hgg125_%s",catIt->Data()))->sumEntries();
    double thisFrac = thisCat/totalEvents;
    double thisFracE = thisFrac * TMath::Sqrt(1/thisCat+1/totalEvents);
    RooRealVar *fSig = new RooRealVar(Form("Data_%s_COMBFIT_%s_fSig",mcName.Data(),catIt->Data()),"",thisFrac,thisFrac-2*thisFracE,thisFrac+2*thisFracE);
    cout << "\tF" <<endl;
    RooFormulaVar *NsigThisCat = new RooFormulaVar(Form("Data_%s_COMBFIT_%s_Nsig",mcName.Data(),catIt->Data()),"","@0*@1",RooArgSet(*fSig,Nsig) );

    //fit model
    RooAddPdf *FitModel = new RooAddPdf(Form("Data_%s_COMBFIT_%s_fitModel",mcName.Data(),catIt->Data()),"Fit Model",
					RooArgList(*SignalModel,*BkgModel),RooArgList(*NsigThisCat,*Nbkg));

    combinedFit.addPdf(*FitModel,*catIt);
  }//category loop

  //data 
  RooDataSet *data = (RooDataSet*)ws->data("Data_Combined");
  
  combinedFit.fitTo(*data,RooFit::ExternalConstraints(allConstraints),RooFit::Strategy(0),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  //combinedFit.fitTo(*data);
  RooFitResult *res=combinedFit.fitTo(*data,RooFit::ExternalConstraints(allConstraints),RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_%s_COMBFIT_fitResult",mcName.Data()) );

  ws->import(*res);
  ws->import(combinedFit);
}

void MakeSpinFits::MakeBackgroundFit(TString mcName, TString catTag,float initMass,float range,bool gausPen,TString inputTag){
  RooRealVar mass = *(ws->var("mass"));
  cout << "1" <<endl;
  if(inputTag=="") inputTag = catTag;
  //build the signal model
  float min = initMass - (gausPen ? 3 : 1)*range; //allow 3 sigma up down
  float max = initMass + (gausPen ? 3 : 1)*range; //allow 3 sigma up down
  RooRealVar mean(Form("Data_%s_FIT_%s_mean",mcName.Data(), catTag.Data()),"mean",initMass,min,max);
  cout << "5" <<endl;
  //won't use this unless gaussian penalty is on
  RooRealVar initialMass(Form("Data_%s_FIT_%s_initialMass",mcName.Data(), catTag.Data()),"Initial mass val",
			 initMass);
  cout << "3" <<endl;
  RooRealVar penaltyMass(Form("Data_%s_FIT_%s_penaltyMass",mcName.Data(), catTag.Data()),"Mass Penalty",
			 range);
  cout << "4" <<endl;
  RooGaussian meanConst("meanConstraint","",mean,initialMass,penaltyMass);
  cout << "5   " << Form("%s_FIT_%s_sigma1",mcName.Data(),inputTag.Data()) << endl;
  RooRealVar sig1Base = *(ws->var(Form("%s_FIT_%s_sigma1",mcName.Data(),inputTag.Data())));
  cout << "6" <<endl;
  RooRealVar sig2Base = *(ws->var(Form("%s_FIT_%s_sigma2",mcName.Data(),inputTag.Data())));
  cout << "7" <<endl;
  RooRealVar fBase    = *(ws->var(Form("%s_FIT_%s_f",mcName.Data(),inputTag.Data())));


  RooRealVar sig1Width("sig1Wdith","",sig1Base.getError()/8);
  RooRealVar sig2Width("sig2Wdith","",sig2Base.getError()/8);
  RooRealVar fWidth("fWdith","",fBase.getError()/8);

  RooRealVar sig1(Form("Data_%s_FIT_%s_sigma1",mcName.Data(), catTag.Data()),"sig1",sig1Base.getVal(),sig1Base.getVal()-4*sig1Base.getError(),sig1Base.getVal()+4*sig1Base.getError());
  RooRealVar sig2(Form("Data_%s_FIT_%s_sigma2",mcName.Data(), catTag.Data()),"sig2",sig2Base.getVal(),sig2Base.getVal()-4*sig2Base.getError(),sig2Base.getVal()+4*sig2Base.getError());
  RooRealVar f(Form("Data_%s_FIT_%s_f",mcName.Data(), catTag.Data()),"f",fBase.getVal(),fBase.getVal()-4*fBase.getError(),fBase.getVal()+4*fBase.getError());
  
  RooGaussian sig1Const("sig1Constraint","",sig1,sig1Base,sig1Width);
  RooGaussian sig2Const("sig2Constraint","",sig2,sig2Base,sig2Width);
  RooGaussian fConst("fConstraint","",f,fBase,fWidth);
  //sig1.setConstant(kTRUE);
  //sig2.setConstant(kTRUE);
  //f.setConstant(kTRUE);
  RooGaussian g1(Form("Data_%s_FIT_%s_g1",mcName.Data(),catTag.Data()),"g1",mass,mean,sig1);
  RooGaussian g2(Form("Data_%s_FIT_%s_g2",mcName.Data(),catTag.Data()),"g2",mass,mean,sig2);
  RooAddPdf SignalModel(Form("Data_%s_FIT_%s_signalModel",mcName.Data(),catTag.Data()),
			"Signal Model",RooArgList(g1,g2),f);

  RooAbsData *ds = ws->data(Form("Data_%s",catTag.Data()));

  //background model
  cout << "8  " << ds  <<endl;
  
  RooRealVar alpha1(Form("Data_%s_FIT_%s_alpha1",mcName.Data(),catTag.Data()),"alpha1",-0.1,-1.,0.);
  RooRealVar alpha2(Form("Data_%s_FIT_%s_alpha2",mcName.Data(),catTag.Data()),"alpha2",-0.1,-1.,0.);
  RooRealVar f_bkg( Form("Data_%s_FIT_%s_fbkg",mcName.Data(),catTag.Data()),"f_bkg",0.1,0,1);
  RooExponential exp1(Form("Data_%s_FIT_%s_exp1",mcName.Data(),catTag.Data()),"exp1",mass,alpha1);
  RooExponential exp2(Form("Data_%s_FIT_%s_exp2",mcName.Data(),catTag.Data()),"exp2",mass,alpha2);
    
  RooAddPdf BkgModel(Form("Data_%s_FIT_%s_bkgModel",mcName.Data(),catTag.Data()),"Background Model",
		     RooArgList(exp1,exp2),f_bkg);
  
  /*
  RooRealVar pC(Form("Data_%s_FIT_%s_pC",mcName.Data(),catTag.Data()),"pC",1,1,1);
  RooRealVar p0(Form("Data_%s_FIT_%s_p0",mcName.Data(),catTag.Data()),"p0",0,-1e6,1e6);
  RooRealVar p1(Form("Data_%s_FIT_%s_p1",mcName.Data(),catTag.Data()),"p1",0,-1e6,1e6);
  RooRealVar p2(Form("Data_%s_FIT_%s_p2",mcName.Data(),catTag.Data()),"p2",0,-1e6,1e6);
  RooRealVar p3(Form("Data_%s_FIT_%s_p3",mcName.Data(),catTag.Data()),"p3",0,-1e6,1e6);

  RooFormulaVar pCmod(Form("Data_%s_FIT_%s_pCmod",mcName.Data(),catTag.Data()),"","@0*@0",pC);
  RooFormulaVar p0mod(Form("Data_%s_FIT_%s_p0mod",mcName.Data(),catTag.Data()),"","@0*@0",p0);
  RooFormulaVar p1mod(Form("Data_%s_FIT_%s_p1mod",mcName.Data(),catTag.Data()),"","@0*@0",p1);
  RooFormulaVar p2mod(Form("Data_%s_FIT_%s_p2mod",mcName.Data(),catTag.Data()),"","@0*@0",p2);
  RooFormulaVar p3mod(Form("Data_%s_FIT_%s_p3mod",mcName.Data(),catTag.Data()),"","@0*@0",p3);

  

  RooArgList *args;
  if(ds->sumEntries() < 300) args = new RooArgList(pCmod,p0mod,p1mod,p2mod);
  else args = new RooArgList(pCmod,p0mod,p1mod,p2mod,p3mod);

  RooBernstein BkgModel(Form("Data_%s_FIT_%s_bkgModel",mcName.Data(),catTag.Data()),"Background Model",mass,*args);
  */
  RooRealVar Nsig(Form("Data_%s_FIT_%s_Nsig",mcName.Data(),catTag.Data()),"N signal Events",100,0,100000);
  RooRealVar Nbkg(Form("Data_%s_FIT_%s_Nbkg",mcName.Data(),catTag.Data()),"N background Events",100,0,1e+09);

  //fit model
  RooAddPdf FitModel(Form("Data_%s_FIT_%s_fitModel",mcName.Data(),catTag.Data()),"Fit Model",
		     RooArgList(SignalModel,BkgModel),RooArgList(Nsig,Nbkg));

  //do the fit

  RooArgSet allConstraints(sig1Const,sig2Const,fConst);
  if(gausPen) allConstraints.add(meanConst);
  


  FitModel.fitTo(*ds,RooFit::ExternalConstraints(allConstraints),RooFit::Strategy(0),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  RooFitResult *res=FitModel.fitTo(*ds,RooFit::ExternalConstraints(allConstraints),RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::Minos(kFALSE),RooFit::Extended(kTRUE));
  res->SetName(Form("Data_%s_FIT_%s_fitResult",mcName.Data(),catTag.Data()) );

  ws->import(FitModel);
  ws->import(*res);

  if(addSWeight) AddSWeight(mcName,catTag,inputTag);
}

void MakeSpinFits::MakeBackgroundOnlyFit(TString mcName,TString catTag){
  if(ws==0) return;
  RooRealVar mass = *(ws->var("mass"));
  cout << "1" <<endl;

  RooAbsData *ds = ws->data(Form("Data_%s",catTag.Data()));

  //background model
  
  RooRealVar alpha1(Form("Data_%s_BKGFIT_%s_alpha1",mcName.Data(),catTag.Data()),"alpha1",-0.1,-1.,0.);
  RooRealVar alpha2(Form("Data_%s_BKGFIT_%s_alpha2",mcName.Data(),catTag.Data()),"alpha2",-0.1,-1.,0.);
  RooRealVar f_bkg( Form("Data_%s_BKGFIT_%s_f",mcName.Data(),catTag.Data()),"f_bkg",0.1,0,1);
  RooExponential exp1(Form("Data_%s_BKGFIT_%s_exp1",mcName.Data(),catTag.Data()),"exp1",mass,alpha1);
  RooExponential exp2(Form("Data_%s_BKGFIT_%s_exp2",mcName.Data(),catTag.Data()),"exp2",mass,alpha2);
    
  RooAddPdf BkgModel(Form("Data_%s_BKGFIT_%s_bkgModel",mcName.Data(),catTag.Data()),"Background Model",
		     RooArgList(exp1,exp2),f_bkg);
  
  /*
  RooRealVar pC(Form("Data_%s_FIT_%s_pC",mcName.Data(),catTag.Data()),"pC",1,1,1);
  RooRealVar p0(Form("Data_%s_FIT_%s_p0",mcName.Data(),catTag.Data()),"p0",0,-1e6,1e6);
  RooRealVar p1(Form("Data_%s_FIT_%s_p1",mcName.Data(),catTag.Data()),"p1",0,-1e6,1e6);
  RooRealVar p2(Form("Data_%s_FIT_%s_p2",mcName.Data(),catTag.Data()),"p2",0,-1e6,1e6);
  RooRealVar p3(Form("Data_%s_FIT_%s_p3",mcName.Data(),catTag.Data()),"p3",0,-1e6,1e6);

  RooFormulaVar pCmod(Form("Data_%s_FIT_%s_pCmod",mcName.Data(),catTag.Data()),"","@0*@0",pC);
  RooFormulaVar p0mod(Form("Data_%s_FIT_%s_p0mod",mcName.Data(),catTag.Data()),"","@0*@0",p0);
  RooFormulaVar p1mod(Form("Data_%s_FIT_%s_p1mod",mcName.Data(),catTag.Data()),"","@0*@0",p1);
  RooFormulaVar p2mod(Form("Data_%s_FIT_%s_p2mod",mcName.Data(),catTag.Data()),"","@0*@0",p2);
  RooFormulaVar p3mod(Form("Data_%s_FIT_%s_p3mod",mcName.Data(),catTag.Data()),"","@0*@0",p3);

  

  RooArgList *args;
  if(ds->sumEntries() < 300) args = new RooArgList(pCmod,p0mod,p1mod,p2mod);
  else args = new RooArgList(pCmod,p0mod,p1mod,p2mod,p3mod);

  RooBernstein BkgModel(Form("Data_%s_FIT_%s_bkgModel",mcName.Data(),catTag.Data()),"Background Model",mass,*args);
  */
  RooRealVar Nbkg(Form("Data_%s_BKGFIT_%s_Nbkg",mcName.Data(),catTag.Data()),"N background Events",ds->sumEntries(),0,1e9);
  
  mass.setRange("fitLow", 110,120);
  mass.setRange("fitHigh",135,170);

  BkgModel.fitTo(*ds,RooFit::Strategy(0),RooFit::Minos(kFALSE),RooFit::Range("fitlow,fithigh"));
  RooFitResult *res=BkgModel.fitTo(*ds,RooFit::Save(kTRUE),RooFit::Strategy(2),RooFit::Minos(kFALSE),RooFit::Range("fitlow,fithigh"));
  res->SetName(Form("Data_%s_BKGFIT_%s_fitResult",mcName.Data(),catTag.Data()) );

  ws->import(BkgModel);
  ws->import(Nbkg);
  ws->import(*res);

}

void MakeSpinFits::MakeAllBackgroundFits(TString cat, TString mcName){
  if(ws==0) return;
  if(cat != "EB_0_0" && cat !="EB_0") 
    MakeBackgroundFit(mcName,cat,mean0,meanE0,true);      
  
  MakeBackgroundFitCosTBin(mcName,cat,-1,-0.3);
  MakeBackgroundFitCosTBin(mcName,cat,-0.3,0.3);
  MakeBackgroundFitCosTBin(mcName,cat,0.3,1);
  //make background only fits
  MakeBackgroundOnlyFit(mcName,cat);
}


void MakeSpinFits::run(){
  if(ws==0) return;

  std::vector<TString>::const_iterator mcIt = mcLabel.begin();
  for(; mcIt != mcLabel.end(); mcIt++){
    std::vector<TString>::const_iterator catIt = catLabels.begin();
    for(; catIt != catLabels.end(); catIt++){
      MakeSignalFit(*catIt,*mcIt);

      if(catIt == catLabels.begin()){
	MakeBackgroundFit(*mcIt,*catIt,124.5,2,false); // do the initial fit
	mean0  = ws->var(Form("Data_%s_FIT_EB_0_mean",mcIt->Data()))->getVal();
	meanE0 = ws->var(Form("Data_%s_FIT_EB_0_mean",mcIt->Data()))->getError();    
      }
      MakeAllBackgroundFits(*catIt,*mcIt);
    }
    MakeCombinedBackgroundFit(*mcIt,125,2);
  }

}

void MakeSpinFits::save(){
  if(ws==0 || outputFile==0) return;
  outputFile->cd();
  ws->Write();
  outputFile->Close();
}

