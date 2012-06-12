//modified code from Y. Yang
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooHistFunc.h"
#include "RooMoment.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooProdPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooBernstein.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include <fstream>

using namespace RooFit;
void setTCanvasNicev1(TCanvas *can0){
  
  
  can0->SetFillColor(0);
  can0->SetBorderMode(0);
  can0->SetBorderSize(2);
  can0->SetTickx(1);
  can0->SetTicky(1);
  can0->SetLeftMargin(0.15);
  can0->SetRightMargin(0.05);
  can0->SetTopMargin(0.05);
  can0->SetBottomMargin(0.13);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
}

int getCatMVA(float mva,float mjj,
	      const int nCat,
	      const float* min,const float* max,const float *mjjmin,const float *mjjmax){
  for(int iCat=0;iCat<nCat;iCat++){
    if(mva>min[iCat] && mva<=max[iCat]
       && mjj > mjjmin[iCat] && mjj <=mjjmax[iCat]) return iCat;
  }
  return -1;
}

int getCatPFCiC( float r9_1,  float r9_2,  float eta_1,  float eta_2,  float mjj){
  if(mjj > 500.) return 0;
  if(mjj > 250.) return 1;
  return ( r9_1<0.94 || r9_2<0.94) + 2*( fabs(eta_1) > 1.48 || fabs(eta_2) > 1.48 );
}		

void makeMinimalTrees(string inputFiles,float mMin=100.,float mMax=180,
		      TString outputFile="rootFiles/smallTree.root", bool doPFCiC=false,bool applyTrigger=true){
  
  const int nCat = 6;
  const float CatMin[nCat] = {0.05,0.05,0.89,0.74,0.55,0.05};
  const float CatMax[nCat] = {9999,9999,9999,0.89,0.74,0.55}; // predefined cuts for the diPhotonMVA categories
  const float MjjMin[nCat] = {500.,250.,-999,-999,-999,-999};
  const float MjjMax[nCat] = {9999,500.,250.,250.,250.,250.};
  
  TChain *fChain = new TChain("HggOutput");
  
  ifstream st(inputFiles.c_str(),fstream::in);
  string file;
  while(st.good()){
    getline(st,file);
    if(file.empty()) continue;
    if(file.find("/castor/cern.ch")!=string::npos) file.insert(0,"rfio://");
    fChain->AddFile(file.c_str());
  }
  Float_t         mPair;
  Float_t         diPhotonMVA;
  Int_t           trigger;
  Float_t           Mjj;
  Int_t            nPho;
  Float_t          pho_r9[2];
  Float_t          pho_px[2];
  Float_t          pho_py[2];
  Float_t          pho_pz[2];
  Float_t          pho_E[2];
  fChain->SetBranchAddress("mPair", &mPair);
  fChain->SetBranchAddress("diPhotonMVA", &diPhotonMVA);
  fChain->SetBranchAddress("trigger",&trigger);
  fChain->SetBranchAddress("Mjj",&Mjj);
  
  fChain->SetBranchAddress("Photon.r9",pho_r9);
  fChain->SetBranchAddress("Photon.p4.fP.fX",pho_px);
  fChain->SetBranchAddress("Photon.p4.fP.fY",pho_py);
  fChain->SetBranchAddress("Photon.p4.fP.fZ",pho_pz);
  fChain->SetBranchAddress("Photon.p4.fE",pho_E);

  TTree *outTree = new TTree("HggOutputReduced","");
  Float_t         mPairOut;
  Int_t           catOut;
  outTree->Branch("mPair", &mPairOut);
  outTree->Branch("cat",&catOut,"cat/I");

  Long64_t ientry = -1;
  while(fChain->GetEntry(++ientry)){
    cout << pho_E[0] << endl;
    TLorentzVector p1(1.,0.,0.,1.);//(pho_px[0],pho_py[0],pho_pz[0],pho_E[0]);
    TLorentzVector p2(1.,0.,0.,1.);//(pho_px[1],pho_py[1],pho_pz[1],pho_E[1]);
    int category;
    if(doPFCiC) category = getCatPFCiC(pho_r9[0],pho_r9[1],p1.Eta(),p2.Eta(),Mjj);
    else category = getCatMVA(diPhotonMVA,Mjj,nCat,CatMin,CatMax,MjjMin,MjjMax);
    if(category == -1 || mPair < mMin || mPair > mMax) continue;
    mPairOut = mPair;
    catOut   = category;
    outTree->Fill();
  }
  TFile *f = new TFile(outputFile,"RECREATE");
  outTree->Write();
  f->Close();
}

void makeBkgWorkSpace(TChain *fChain, float lumi, float mMin=100.,float mMax=180,
		      TString outputFile="interpolated/BkgWorkSpace.root", bool doPFCiC=false,bool applyTrigger=true){

  const int nCat = 6;
  const float CatMin[nCat] = {0.05,0.05,0.89,0.74,0.55,0.05};
  const float CatMax[nCat] = {9999,9999,9999,0.89,0.74,0.55}; // predefined cuts for the diPhotonMVA categories
  const float MjjMin[nCat] = {500.,250.,-999,-999,-999,-999};
  const float MjjMax[nCat] = {9999,500.,250.,250.,250.,250.};
  //gStyle->SetErrorX(0); 
  //gStyle->SetOptStat(0);


  Float_t         mPair;
  Int_t           cat;
  fChain->SetBranchAddress("mPair", &mPair);
  fChain->SetBranchAddress("cat",&cat);
  
  RooRealVar *rv_mass = new RooRealVar("mass","mass",100,mMin,mMax);
  //rv_mass->setRange(100,150);
  rv_mass->setRange(mMin,mMax);
  rv_mass->setBins(mMax-mMin);
  rv_mass->SetTitle("m_{#gamma#gamma}");
  rv_mass->setUnit("GeV/c^{2}");
  rv_mass->setRange("fitrange",mMin,mMax);
    
   RooRealVar *IntLumi = new RooRealVar("IntLumi","",lumi,0,10E3);
   IntLumi->setConstant();
   
   RooDataSet*  rds_data_mass[nCat+1]; 
   for(int j=0;j<nCat+1;j++){ //last one is the combination
     string rname = string( Form("data_mass_cat%d",j) );
     rds_data_mass[j] = new RooDataSet(rname.c_str(),rname.c_str(),*rv_mass);
   }
  
  int Nevents = fChain->GetEntries();
  
  vector<float> allmass; 

  vector<TString> catnames;
  for(int i=0;i<nCat;i++){
    catnames.push_back(Form("cat_%d",i));
  }

  //loop over events
  int nSelected[nCat] = {0,0,0,0,0,0};
  for(int iEvent=0; iEvent< Nevents; iEvent++){
    fChain->GetEntry(iEvent);

    if(iEvent%500==0){
      cout << "Processing Event " << iEvent << "  Selected: " << endl;
      for(int i=0;i<nCat;i++) cout << "\t cat" << i  <<  ": " << nSelected[i];
      cout << endl;
    }

    if(cat==-1) continue; // not a good diPhoton event
    nSelected[cat]++;
    if(mPair >=mMin && mPair <= mMax){
      rv_mass->setVal(mPair);
      rds_data_mass[cat]->add(*rv_mass);
      rds_data_mass[nCat]->add(*rv_mass);
      allmass.push_back(mPair);
    }
  }  
  RooWorkspace *w = new RooWorkspace("cms_hgg_workspace","wbkg") ;

  std::vector<RooAbsData*> datav;
  std::vector<RooAbsPdf*> pdfv;
  
  std::vector<RooRealVar*> coeffv;
  std::vector<RooRealVar*> normv;
  
  RooRealVar *p1first=0;
  RooRealVar *p2first=0;
  RooRealVar *nfirst=0;
  
  RooCategory finalcat("finalcat","finalcat") ;  
  RooSimultaneous fullbkgpdf("fullbkgpdf","fullbkgpdf",finalcat);
  RooDataSet datacomb("datacomb","",RooArgList(*rv_mass,finalcat)) ;
  RooDataSet *datacombcat = new RooDataSet("data_combcat","",RooArgList(*rv_mass)) ;
  
  
  
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    
    TString catname = catnames.at(icat);
    finalcat.defineType(catname);
    
    RooDataSet *indata = rds_data_mass[icat];
    RooDataSet *data = new RooDataSet(TString("data_")+catname,"",*rv_mass,Import(*indata));
    RooDataSet *datacat = new RooDataSet(TString("datacat")+catname,"",*rv_mass,Index(finalcat),Import(catname,*data)) ;
    datacomb.append(*datacat);
    
    
    RooRealVar *p1 = new RooRealVar(TString::Format("CMS_hgg_%s_p1",catname.Data()),"",0.1,-10,10);
    RooRealVar *p2 = new RooRealVar(TString::Format("CMS_hgg_%s_p2",catname.Data()),"",0.1,-10,10);
    RooRealVar *p3 = new RooRealVar(TString::Format("CMS_hgg_%s_p3",catname.Data()),"",0.1,-10,10);
    RooRealVar *p4 = new RooRealVar(TString::Format("CMS_hgg_%s_p4",catname.Data()),"",0.1,-10,10);
    RooRealVar *p5 = new RooRealVar(TString::Format("CMS_hgg_%s_p5",catname.Data()),"",0.1,-10,10);
    
//     p1->setVal(1);
//     p2->setVal(1);
//     p3->setVal(1);
//     p4->setVal(1);
//     p5->setVal(1);
    
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString("p1mod"+catname),"","@0*@0",*p1);
    RooFormulaVar *p2mod = new RooFormulaVar(TString("p2mod"+catname),"","@0*@0",*p2);
    RooFormulaVar *p3mod = new RooFormulaVar(TString("p3mod"+catname),"","@0*@0",*p3);
    RooFormulaVar *p4mod = new RooFormulaVar(TString("p4mod"+catname),"","@0*@0",*p4);
    RooFormulaVar *p5mod = new RooFormulaVar(TString("p5mod"+catname),"","@0*@0",*p5);

  
//     coeffv.push_back(p1);
//     coeffv.push_back(p2);
//     coeffv.push_back(p3);
//     coeffv.push_back(p4);
//     coeffv.push_back(p5);
    
    
    RooRealVar *alpha1 = new RooRealVar(TString("alpha1"+catname),"",-0.1,-5.0,5.0);
    RooRealVar *alpha2 = new RooRealVar(TString("alpha2"+catname),"",-0.1,-5.0,5.0);
    //alpha2->removeRange();
    
    
    RooRealVar *nbkg = new RooRealVar(TString::Format("CMS_hgg_%s_bkgshape_norm",catname.Data()),"",100.0,0.0,10e3);
    nbkg->setVal(indata->numEntries());
    normv.push_back(nbkg);
    
    RooBernstein *bkgshape = new RooBernstein(TString::Format("CMS_hgg_%s_bkgshape",catname.Data()),"",*rv_mass,RooArgList(RooConst(1.0),*p1mod,*p2mod,*p3mod,*p4mod,*p5mod));
    RooExtendPdf *bkgpdf = new RooExtendPdf(TString("bkgpdf")+catname,"",*bkgshape,*nbkg);
    
    fullbkgpdf.addPdf(*bkgpdf,catname);
    datav.push_back(data);
    pdfv.push_back(bkgpdf);
    
    RooDataHist *databinned = new RooDataHist(TString("databinned_")+catname,"",*rv_mass,*data);
    

    //    TH1D* dataH = (TH1D*)databinned->createHistogram(TString("dataH_")+catname,data);
    //TH1D* fitH  = (TH1D*)bkgpdf->createHistogram(TString("fitH_")+catname,nbkg);
    //TFile *f = TFile::Open("out.root");
    //fitH->Write();
    //dataH->Write();
    //f->Close();


    w->import(*data);
    w->import(*databinned);
    
  }
  catnames.push_back("combcat");    
  ///last one is combined
  datav.push_back(rds_data_mass[nCat]);  
  
  w->import(datacomb);  

  //gStyle->SetOptStat(1110);
  
  ///here is to get all p1 and p2 for 8 categories by simulateous fit 

  //first fit with strategy(0) fast , and strategy(2) to improve robustness
  fullbkgpdf.fitTo(datacomb,Strategy(0),Minos(kFALSE),Save(kTRUE));
  RooFitResult *fullbkgfitres = fullbkgpdf.fitTo(datacomb,Strategy(2),Minos(kFALSE),Save(kTRUE));
  
  
  int n_alldata = 0; 
  for(int j=0;j <nCat; j++){
    cout<< "rds_data_mass[j]->numEntries()" << rds_data_mass[j]->numEntries()<<endl; 
    n_alldata += rds_data_mass[j]->numEntries(); 
  }
  //return; 
  
  //return; 
  
  
  //const UInt_t nbins = 0.5*(mMax - mMin);
  const UInt_t nbins = 1*(mMax - mMin);
  
   
  RooRealVar *p1test = new RooRealVar(TString::Format("CMS_hgg_p1test"),"",0.1,-10,10);
  RooRealVar *p2test = new RooRealVar(TString::Format("CMS_hgg_p2test"),"",0.1,-10,10);
  RooRealVar *p3test= new RooRealVar(TString::Format("CMS_hgg_p3test"),"",0.1,-10,10);
  RooRealVar *p4test = new RooRealVar(TString::Format("CMS_hgg_p4test"),"",0.1,-10,10);
  RooRealVar *p5test = new RooRealVar(TString::Format("CMS_hgg_p5test"),"",0.1,-10,10);
  
  RooFormulaVar *p1testmod = new RooFormulaVar(TString("p1testmod"),"","@0*@0",*p1test);
  RooFormulaVar *p2testmod = new RooFormulaVar(TString("p2testmod"),"","@0*@0",*p2test);
  RooFormulaVar *p3testmod = new RooFormulaVar(TString("p3testmod"),"","@0*@0",*p3test);
  RooFormulaVar *p4testmod = new RooFormulaVar(TString("p4testmod"),"","@0*@0",*p4test);
  RooFormulaVar *p5testmod = new RooFormulaVar(TString("p5testmod"),"","@0*@0",*p5test);
  
  p1test->setVal(1);
  p2test->setVal(1);
  p3test->setVal(1);
  p4test->setVal(1);
  p5test->setVal(1);
  
  
  RooBernstein *bkgshape = new RooBernstein("bkgshape","",*rv_mass,RooArgList(RooConst(1.0),*p1testmod,*p2testmod,*p3testmod,*p4testmod,*p5testmod));
  RooRealVar *nbkgindv = new RooRealVar(TString("CMS_hgg_indv_bkgshape_norm"),"",500.0,0.0,10e5);
  RooExtendPdf *bkgpdfindv = new RooExtendPdf(TString("bkgpdfall"),"",*bkgshape,*nbkgindv);
  
  
  
  //float ymax[10] = {600,800,600,800,2000};
  float ymax[10] = {20,40,40,100,200};
  
  string ptname[2] = {"p^{#gamma#gamma}_{T} > 40 GeV/c","p^{#gamma#gamma}_{T} < 40 GeV/c"};
  string mvaname[6] = { "M_{jj} > 500 GeV  MVA > 0.05", "500 > M_{jj} > 250 GeV  MVA > 0.05", "MVA > 0.89", 
			"0.89 > MVA > 0.74" , "0.74 > MVA > 0.55","0.55 > MVA > 0.05" };
  string cicname[6] = { "M_{jj} > 500 GeV", "500 > M_{jj} > 250 GeV", "Barrel High R9", "Barrel Low R9",
			"Endcap High R9", "Endcap Low R9"};
			
  
  
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {
    
    TString catname = catnames.at(icat);
        
    p1test->setVal(0.5);
    p2test->setVal(0.5);
    p3test->setVal(0.5);
    p4test->setVal(0.5);
    p5test->setVal(0.5);
    nbkgindv->setVal(datav.at(icat)->numEntries());
    
    //RooFitResult *fitres = bkgshape->fitTo(*datav.at(icat),Strategy(2),Minos(kFALSE),Save(kTRUE));
    ///extended fit on individual category
    RooFitResult *fitres =  bkgpdfindv->fitTo(*datav.at(icat),Strategy(2),Minos(kFALSE),Save(kTRUE));
    
    
    TCanvas *cbkg = new TCanvas;
    setTCanvasNicev1(cbkg);
    
    RooPlot *plot = rv_mass->frame(Bins(nbins),Range("fitrange"));
    datav.at(icat)->plotOn(plot,LineColor(kWhite),MarkerColor(kWhite));
    
    if( icat < nCat ){
      
      //here from simulatenous fit (extended, normalization is floating in each category) 
      pdfv.at(icat)->plotOn(plot,FillColor(kGreen),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fullbkgfitres,2.0)); 
      pdfv.at(icat)->plotOn(plot,FillColor(kYellow),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fullbkgfitres,1.0));
      pdfv.at(icat)->plotOn(plot,RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
      
      ///from this individual extended fit
      //bkgpdfindv->plotOn(plot,FillColor(kGreen),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fitres,2.0));
      //bkgpdfindv->plotOn(plot,FillColor(kYellow),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fitres,1.0));
      //bkgpdfindv->plotOn(plot,RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));
      
      
    }else{
      bkgpdfindv->plotOn(plot,FillColor(kGreen),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fitres,2.0)); /// ,Normalization(1.0));    
      bkgpdfindv->plotOn(plot,FillColor(kYellow),Range("fitrange"),NormRange("fitrange"),VisualizeError(*fitres,1.0)); ////,Normalization(1.0));
      bkgpdfindv->plotOn(plot,RooFit::LineColor(kRed),Range("fitrange"),NormRange("fitrange"));

    }
    
    datav.at(icat)->plotOn(plot);
    plot->SetTitle("");      
    plot->SetMinimum(0.0);
    plot->SetMaximum(ymax[icat]);
    
    
    
    plot->Draw();    
    TLegend *legmc = new TLegend(0.63,0.7,0.97,0.93);  
    legmc->AddEntry(plot->getObject(4),"Data","LPE");
    legmc->AddEntry(plot->getObject(3),"Bkg Model","L");
    legmc->AddEntry(plot->getObject(2),"#pm1 #sigma","F");  
    legmc->AddEntry(plot->getObject(1),"#pm2 #sigma","F");  
    
    //legmc->AddEntry(plot->getObject(5),"Chebychev Fit","L");  
    
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();   
    
    
    if( icat < nCat ){
      string textname = (doPFCiC ? cicname[icat] : mvaname[icat]);
      TLatex *tex = new TLatex(0.2,0.8,textname.c_str());
      tex->SetNDC();
      tex->SetTextSize(0.04);
      tex->Draw();
    }else{
      TLatex *   tex = new TLatex(0.2,0.8,"Combined categories");
      tex->SetNDC();
      tex->SetTextSize(0.04);
      tex->Draw();
    }
    
    TString texname = TString(Form("L = %4.3f fb^{-1}",lumi));
    TLatex *tex = new TLatex(0.2,0.9, texname);
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->Draw();
    

    cbkg->SaveAs(TString("figs/databkg") + catname + TString(".eps"));
    cbkg->SaveAs(TString("figs/databkg") + catname + TString(".gif"));
    cbkg->SaveAs(TString("figs/databkg") + catname + TString(".C"));
    ///normv.at(icat)->setRange(0.0,normv.at(icat)->getVal()+8.0*sqrt(normv.at(icat)->getVal()));
    //normv.at(icat)->setConstant();
    //return; 
    
    
  }
  

  //  return; 
  

//   for (UInt_t i=0; i<coeffv.size(); ++i) {
//     coeffv.at(i)->setRange(0.0,20.0);
//   }
  
  

  TFile *fnew = new TFile(outputFile,"recreate");
  
  
  w->import(datacomb);  
  w->import(fullbkgpdf);
  w->import(*fullbkgfitres);
  w->Print();
  w->Write();
  

  printf("IntLumi = %5f\n",IntLumi->getVal());
  printf("ndata:\n");
  for (UInt_t icat=0; icat<catnames.size(); ++icat) {    
    //printf("%u: ndata = %i\n",icat,datav.at(icat)->numEntries());
    printf("%i ",datav.at(icat)->numEntries());
    
  }   
  printf("\n");

  fnew->Write();
  fnew->Close();
  
  
  

}




void makeBkgWorkSpace(string inputFiles, float lumi, float mMin=100.,float mMax=180,
		      TString outputFile="interpolated/BkgWorkSpace.root", bool doPFCiC=false,bool applyTrigger=true){
  TChain *fChain = new TChain("HggOutputReduced");
  
  ifstream st(inputFiles.c_str(),fstream::in);
  string file;
  while(st.good()){
    getline(st,file);
    if(file.empty()) continue;
    if(file.find("/castor/cern.ch")!=string::npos) file.insert(0,"rfio://");
    fChain->AddFile(file.c_str());
  }

  makeBkgWorkSpace(fChain,lumi,mMin,mMax,outputFile,doPFCiC,applyTrigger);
}
