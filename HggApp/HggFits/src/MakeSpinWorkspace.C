#include "MakeSpinWorkspace.h"

using namespace std;
#include "selectionMaps.C"

MakeSpinWorkspace::MakeSpinWorkspace(TString outputFileName):
  mixer(0),
  fileKFactor(0),
  fileRescaleFactor(0)
{
  //setup the defaults

  requireCiC=true;
  tightPt=false;
  selectionMap=0;

  runLow=0;
  runHigh= 9999999;

  useR9=false;

  setMassRange(100.,180.);
  setPtCuts(32.,24.);

  EfficiencyCorrectionFile="";
  filenameKFactor="";
  filenameRescaleFactor="";

  //open output file and create the new workspace
  outputFile = TFile::Open(outputFileName,"RECREATE");
  ws = new RooWorkspace();
  ws->SetName("cms_hgg_spin_workspace");  
  labels = new RooCategory("labels","labels");
  
}

MakeSpinWorkspace::~MakeSpinWorkspace(){
  if(mixer) delete mixer;
}


TH2F* MakeSpinWorkspace::getSelectionMap(int map,bool isData){
  // get selection map
  // needs to be cleaned up, but useful for testing atm
  switch(map){
  case 0:
    return getSelectionMap0();
  case 1:
    return getSelectionMap1();
  case 2:
    return getSelectionMap2();
  case 3:
    return getSelectionMap3();
  case 4:
    return getSelectionMap4(isData);
  case 5:
    return getSelectionMap5(isData);
  case 6:
    return getSelectionMap6(isData);
  case 7:
    return getSelectionMap7(isData);
  }
  return 0;
}

int MakeSpinWorkspace::passSelection(TH2F* map,float sigEoE,float etaSC, float pt){
  if(sigEoE < map->GetBinContent(map->FindFixBin(fabs(etaSC),pt))) return 0;
  return 1;
}

int MakeSpinWorkspace::passSelection(float r9){
  if(r9 > 0.94) return 0;
  return 1;
}
bool MakeSpinWorkspace::getBaselineSelection(HggOutputReader2* h,int maxI,int minI,float mass){
  if(h->nPhoton < 2
     || h->diPhotonMVA<-1
     || mass < mMin
     || mass > mMax
     || h->Photon_pt[minI] < pt2Min
     || h->Photon_pt[maxI] < pt1Min) return false;
  return true;
}

void MakeSpinWorkspace::AddToWorkspace(TString inputFile,TString tag, bool isData){
  //adds the data/MC to the workspace with the given tag
  cout << "saveDataSet" <<endl;
  TFile *f = TFile::Open(inputFile);
  TTree *tree = (TTree*)f->Get("HggOutput");

  //define a reader class for the input tree
  HggOutputReader2 h(tree);
  
  std::cout << "Get Selection Map" << std::endl;
  TH2F* map = getSelectionMap(selectionMap,isData);
  std::cout << map->GetBinContent(1,1) << std::endl;

  TFile *efficiencyCorrection=0;
  if(EfficiencyCorrectionFile!=""){ //Use MC derived correction weights for the photons
    //as a function of eta/R9/phi
    efficiencyCorrection = new TFile(EfficiencyCorrectionFile);
  }

  if(filenameKFactor!=""){
    if(tag=="Hgg125") fileKFactor = new TFile(filenameKFactor);
    else fileKFactor=0;
  }

  if(filenameRescaleFactor!=""){
    fileRescaleFactor = new TFile(filenameRescaleFactor);
  }

  //define the variables for the workspace
  RooRealVar* mass   = new RooRealVar("mass",  "Mass [GeV]", 100., 180.);
  RooRealVar* cosT   = new RooRealVar("cosT",  "cos(theta)", -1, 1);
  RooRealVar* sige1  = new RooRealVar("sigEoE1","#sigma_{E}/E Lead Photon",0,0.1);
  RooRealVar* sige2  = new RooRealVar("sigEoE2","#sigma_{E}/E SubLead Photon",0,0.1);
  RooRealVar* evtW   = new RooRealVar("evtWeight","Event Weight",1,0,1e6);

  RooRealVar* eta1 = new RooRealVar("eta1","#eta Lead Photon",0,-3,3);
  RooRealVar* eta2 = new RooRealVar("eta2","#eta SubLead Photon",0,-3,3);

  RooRealVar* etaSC1 = new RooRealVar("etaSC1","SC #eta Lead Photon",0,-3,3);
  RooRealVar* etaSC2 = new RooRealVar("etaSC2","SC #eta SubLead Photon",0,-3,3);

  RooRealVar* phi1 = new RooRealVar("phi1","#phi Lead Photon",0,0,6.3);
  RooRealVar* phi2 = new RooRealVar("phi2","#phi SubLead Photon",0,0,6.3);

  RooRealVar* pt1 = new RooRealVar("pt1","p_{T} Lead Photon",0,0,2e3);
  RooRealVar* pt2 = new RooRealVar("pt2","p_{T} SubLead Photon",0,0,2e3);

  RooRealVar* r91 = new RooRealVar("r91","R_{9} Lead Photon",0,0,1.3);
  RooRealVar* r92 = new RooRealVar("r92","R_{9} SubLead Photon",0,0,1.3);

  RooRealVar* idMVA1 = new RooRealVar("idMVA1","ID MVA Lead Photon",0,-1,1.);
  RooRealVar* idMVA2 = new RooRealVar("idMVA2","ID MVA SubLead Photon",0,-1,1.);

  RooRealVar* diPhotonMVA = new RooRealVar("diPhotonMVA","DiPhoton MVA",0,-1,1);

  RooCategory *cat = new RooCategory("evtcat","evtcat");

  //create the RooArgSet to initialize the datasets
  RooArgSet set;
  set.add(*mass);
  set.add(*cosT);
  set.add(*evtW);
  /*
  set.add(*sige1);
  set.add(*sige2);
  set.add(*eta1); set.add(*eta2);
  set.add(*etaSC1); set.add(*etaSC2);
  set.add(*phi1); set.add(*phi2);
  set.add(*pt1); set.add(*pt2);
  set.add(*r91); set.add(*r92);
  set.add(*idMVA1); set.add(*idMVA2);
  set.add(*diPhotonMVA);
  */
  RooRealVar *totEB = new RooRealVar(Form("%s_EB_totalEvents",tag.Data()),"",0,0,1e9);
  RooRealVar *totEE = new RooRealVar(Form("%s_EE_totalEvents",tag.Data()),"",0,0,1e9);

  std::map<std::pair<int,int>, RooDataSet*> dataMapEB, dataMapEE;
  std::map<std::pair<int,int>, RooDataSet*> *datamap;

  if(useR9){
    //if we are using the CiC categories, we only have nCat categories per detector region, so we don't use the 
    //second index of the pair and define fewer categories
    for(int j=0;j<nCat;j++){
      cat->defineType( Form("EB_%d",j),2*j );
      cat->defineType( Form("EE_%d",j),2*j+1 );
      dataMapEB[std::pair<int,int>(j,0)] = new RooDataSet(Form("%s_EB_%d",tag.Data(),j),"",set);
      dataMapEE[std::pair<int,int>(j,0)] = new RooDataSet(Form("%s_EE_%d",tag.Data(),j),"",set);
    }
  }else{
    //using the sigmaE/E categories, we have nCat^2 categories
    for(int i=0;i<nCat;i++){
      for(int j=0;j<nCat;j++){
	cat->defineType( Form("EB_%d_%d",i,j),4*i+2*j );
	cat->defineType( Form("EE_%d_%d",i,j), 4*i+2*j+1 );
	dataMapEB[std::pair<int,int>(i,j)] = new RooDataSet(Form("%s_EB_%d_%d",tag.Data(),i,j),"",set);
	dataMapEE[std::pair<int,int>(i,j)] = new RooDataSet(Form("%s_EE_%d_%d",tag.Data(),i,j),"",set);
      }
    }
  }

  Long64_t iEntry=-1;
  cout << "Making DataSet" << endl;
  double nEB=0,nEE=0; // store the total number of events before any cuts or selection

  //begin main loop over tree
  while(h.GetEntry(++iEntry)){

    if( !(iEntry%10000) ) cout << "Processing " << iEntry << "\t\t" << h.evtWeight << endl;
    //determine  the leading and trailing photons
    int maxI = (h.Photon_pt[1] > h.Photon_pt[0] ? 1:0);
    int minI = (h.Photon_pt[1] > h.Photon_pt[0] ? 0:1);
    if(isData && (h.runNumber < runLow || h.runNumber > runHigh) ) continue;

    if(fabs(h.Photon_etaSC[1]) < 1.48 && fabs(h.Photon_etaSC[0]) < 1.48) nEB+=h.evtWeight;
    else nEE+=h.evtWeight;

    float m = (useUncorrMass ? h.mPairNoCorr : h.mPair);

    // apply selections
    if(!getBaselineSelection(&h,maxI,minI,m)) continue;
    if(tightPt && (h.Photon_pt[maxI]/m < 1./3. || h.Photon_pt[minI]/m < 1./4.) ) continue;


    if(requireCiC){
      if(h.Photon_passPFCiC[1]==false || h.Photon_passPFCiC[0]==false) continue;
    }

    float se1 = h.Photon_EError[maxI]/h.Photon_E[maxI];
    float se2 = h.Photon_EError[minI]/h.Photon_E[minI];

    int p1,p2;
    //determine the photon categories
    if(useR9){
      p1 = passSelection(h.Photon_r9[maxI]);
      p2 = passSelection(h.Photon_r9[minI]);
    }else{
      p1 = passSelection(map,se1,h.Photon_etaSC[maxI],h.Photon_pt[maxI]);
      p2 = passSelection(map,se2,h.Photon_etaSC[minI],h.Photon_pt[minI]);
    }
    if(p1 >= nCat || p2 >= nCat) continue; //we can veto photons here
    datamap = ((fabs(h.Photon_etaSC[maxI]) < 1.48 && fabs(h.Photon_etaSC[minI]) < 1.48) ?
	       &dataMapEB : & dataMapEE); //choses the correct dataset to add the event to
      
    //set all the variables
    mass->setVal(m);
    cosT->setVal(h.cosThetaLead);
    sige1->setVal(se1);
    sige2->setVal(se2);

    eta1->setVal(h.Photon_eta[maxI]);
    eta2->setVal(h.Photon_eta[minI]);
    
    etaSC1->setVal(h.Photon_etaSC[maxI]);
    etaSC2->setVal(h.Photon_etaSC[minI]);
    
    phi1->setVal(h.Photon_phi[maxI]);
    phi2->setVal(h.Photon_phi[minI]);
    
    pt1->setVal(h.Photon_pt[maxI]);
    pt2->setVal(h.Photon_pt[minI]);
    
    r91->setVal(h.Photon_r9[maxI]);
    r92->setVal(h.Photon_r9[minI]);
    
    idMVA1->setVal(h.Photon_idMVA[maxI]);
    idMVA2->setVal(h.Photon_idMVA[minI]);
    
    diPhotonMVA->setVal(h.diPhotonMVA);
    //efficiencyCorrection=0
    float pho1EffWeight = getEffWeight(efficiencyCorrection,eta1->getVal(),pt1->getVal(),phi1->getVal(),r91->getVal());
    float pho2EffWeight = getEffWeight(efficiencyCorrection,eta2->getVal(),pt2->getVal(),phi2->getVal(),r92->getVal());

    //set the event weight
    if(!isData) evtW->setVal(h.evtWeight*pho1EffWeight*pho2EffWeight*getEfficiency(h,125));
    else evtW->setVal(1*pho1EffWeight*pho2EffWeight);
    if( !(iEntry%1000) ) cout <<  "\t\t\t" << evtW->getVal() << endl;


    //if(evtW->getVal() > 7.4) cout << iEntry << ": " << mass->getVal() << " : " << cosT->getVal() << " : " << evtW->getVal() << endl;
    //cout << pho1EffWeight << "\t" <<pho2EffWeight << "\t" << evtW->getVal()<<endl;


    if(useR9){
      int cat = (p1==0 && p2==0 ? 0 : 1);
      (*datamap)[std::pair<int,int>(cat,0)]->add(set);
    }else{
      (*datamap)[std::pair<int,int>(p1,p2)]->add(set);
    }
  }
  cout << "Processed " << iEntry << " Entries" <<endl;

  totEB->setVal(nEB);
  totEE->setVal(nEE);

  //build a combined DataSet using the cat RooCategory
  RooArgSet setCat(set);
  setCat.add(*cat);
  RooDataSet* dataComb = new RooDataSet(tag+"_Combined","",setCat);
  
  std::map<std::pair<int,int>, RooDataSet*>::iterator dIt;
  for(dIt = dataMapEB.begin();dIt!=dataMapEB.end();dIt++){
    //loop over the barrel datasets and add the individual data to the combined dataset
    //dIt->setWeightVar(*evtW);
    TString cattag;
    if(useR9) cattag = Form("EB_%d", (dIt->first).first );
    else cattag = Form("EB_%d_%d", (dIt->first).first, (dIt->first).second );
    std::cout << cattag <<std::endl;
    RooDataSet *tmp = new RooDataSet("DataCat_"+cattag,"",set,RooFit::Index(*cat),RooFit::Import(cattag,*(dIt->second)) );
    dataComb->append(*tmp);
    //ws->import(*(dIt->second));
  }
  for(dIt = dataMapEE.begin();dIt!=dataMapEE.end();dIt++){
    //loop over the endcap datasets and add the individual data to the combined dataset
    //dIt->setWeightVar(*evtW);
    TString cattag;
    if(useR9) cattag = Form("EE_%d", (dIt->first).first );
    else cattag = Form("EE_%d_%d", (dIt->first).first, (dIt->first).second );
    std::cout << cattag <<std::endl;
    RooDataSet *tmp = new RooDataSet("DataCat_"+cattag,"",set,RooFit::Index(*cat),RooFit::Import(cattag,*(dIt->second)) );
    dataComb->append(*tmp);    
    //ws->import(*(dIt->second));
  }

  RooDataSet *dataComb_w = new RooDataSet(tag+"_Combined","",dataComb,setCat,0,"evtWeight");
  //import everything
  ws->import(*dataComb_w);
  ws->import(*totEB);
  ws->import(*totEE);
  ws->import(*evtW);
  cout << "Done" <<endl;

  if(efficiencyCorrection) delete efficiencyCorrection;
  if(fileKFactor) fileKFactor->Close();
  if(fileRescaleFactor) fileRescaleFactor->Close();
}

float MakeSpinWorkspace::getEffWeight(TFile *effFile, float eta, float pt, float phi, float r9){

  if(effFile==0) return 1;

  TH3F* effMap=0;
  TList *keyList = effFile->GetListOfKeys();
  const char* prefix = (requireCiC ? "cic" : "pre");
  for(int i=0;i<keyList->GetSize(); i++){
    TObjArray *o = TString(keyList->At(i)->GetName()).Tokenize("_");
    //the keys have names of the form <effType>_minPt_maxPt
    if( strcmp(prefix,o->At(0)->GetName()) != 0) continue;
    if( atof(o->At(1)->GetName()) > pt) continue; //if the min pt is above the photon pt
    if( atof(o->At(2)->GetName()) <= pt) continue; // if the max pt is below the photon pt
    effMap = (TH3F*)effFile->Get(keyList->At(i)->GetName());
    break;
  }
  if(effMap==0) return 1;
  if(effMap->GetBinContent(effMap->FindFixBin(eta,phi,r9))==0){ // try to interpolate the above and below in phi
    int ix,iy,iz;
    int index = effMap->FindFixBin(eta,phi,r9); // global coordinate
    effMap->GetBinXYZ(index,ix,iy,iz);
    int indexYup   = (iy == effMap->GetNbinsY() ? 1 : iy+1); //wrap around if we are at the top edge of the map
    int indexYdown = (iy == 1 ? effMap->GetNbinsY() : iy-1); //wrap around if we are at the bottom edge of the map
    
    float up = effMap->GetBinContent(ix,indexYup,iz);
    float down = effMap->GetBinContent(ix,indexYdown,iz);
    if(up==0){
      if(down!=0) return 1./down;
      else return 1;
    }
    if(down==0) return 1./up;
    return 2./(up+down); //return 1/avg
  }


  return 1./effMap->GetBinContent(effMap->FindFixBin(eta,phi,r9));
}

void MakeSpinWorkspace::addFile(TString fName,TString l,bool is){
  fileName.push_back(fName);
  label.push_back(l);
  isData.push_back(is);
  if(!is) labels->defineType(l,labels->numBins(""));
}


void MakeSpinWorkspace::MakeWorkspace(){
  //run AddToWorkspace(...) on all the input files
  std::vector<TString>::const_iterator fnIt,lIt;
  std::vector<bool>::const_iterator idIt;

  fnIt = fileName.begin();
  lIt  = label.begin();
  idIt = isData.begin();
  for(; fnIt != fileName.end(); fnIt++, lIt++, idIt++){
    AddToWorkspace(*fnIt,*lIt,*idIt);
  }
  ws->import(*labels);

  if(mixer) mixer->mixAll();
  
  outputFile->cd();
  ws->Write();
  outputFile->Close();
}

float MakeSpinWorkspace::getEfficiency(HggOutputReader2 &h, int massPoint){
  float KFac = getKFactor(h,massPoint);
  float rescale = getRescaleFactor(h);

  return KFac*rescale;
}

float MakeSpinWorkspace::getKFactor(HggOutputReader2 &h, int massPoint){
  if(fileKFactor==0) return 1;
  TString name = Form("kfact%d_0",massPoint);

  TH1D* hist = (TH1D*)fileKFactor->Get(name);

  float kf=1;
  assert(hist!=0);
  if(hist){
    kf = hist->GetBinContent(hist->FindFixBin(h.genHiggsPt));
  }
  delete hist;
  return kf;
}

float MakeSpinWorkspace::getRescaleFactor(HggOutputReader2 &h){
  if(fileRescaleFactor==0) return 1;

  float EFF=1;

  int cicCat = 2*(fabs(h.Photon_etaSC[0])>1.48 || fabs(h.Photon_etaSC[1])>1.48)
    + (h.Photon_r9[0]<0.94 || h.Photon_r9[1]<0.94);

  TLorentzVector p4_1; p4_1.SetPtEtaPhiM(h.Photon_pt[0],h.Photon_eta[0],h.Photon_phi[0],0);
  TLorentzVector p4_2; p4_2.SetPtEtaPhiM(h.Photon_pt[1],h.Photon_eta[1],h.Photon_phi[1],0);

  float pt_sys = (p4_1+p4_1).Pt();

  //L1-HLT factor
  TGraphAsymmErrors* e = (TGraphAsymmErrors*)fileRescaleFactor->Get( Form("effL1HLT_cat%d",cicCat) );
  EFF*=getEffFromTGraph(e,pt_sys);

  //vertex efficiencies
  e = (TGraphAsymmErrors*)fileRescaleFactor->Get( Form("ratioVertex_cat%d_pass",cicCat) );
  EFF*=getEffFromTGraph(e,pt_sys);
  e = (TGraphAsymmErrors*)fileRescaleFactor->Get( Form("ratioVertex_cat%d_fail",cicCat) );
  EFF*=getEffFromTGraph(e,pt_sys);

  //per-photon efficiencies
  TString labels[4] = {"EBHighR9","EBLowR9","EEHighR9","EELowR9"};
  for(int i=0;i<2;i++){
    int index = 2*(fabs(h.Photon_etaSC[i])>1.48)+(h.Photon_r9[i] > 0.94);
    e = (TGraphAsymmErrors*)fileRescaleFactor->Get( TString("ratioTP_")+labels[index] );
    EFF*=getEffFromTGraph(e,h.Photon_pt[i]);
    e = (TGraphAsymmErrors*)fileRescaleFactor->Get( TString("ratioR9_")+labels[index] );
    EFF*=getEffFromTGraph(e,h.Photon_pt[i]);
  }
  return EFF;
}


float MakeSpinWorkspace::getEffFromTGraph(TGraphAsymmErrors* e,float pt){
  int nBins = e-> GetN();
  int selectedBin=0;
  for(selectedBin=0;selectedBin<nBins;selectedBin++){
    Double_t x,y;
    e->GetPoint(selectedBin,x,y);
    if(pt < x) break; // found the first bin above the pT
  }
  if(selectedBin==nBins) selectedBin--;                       

  if(selectedBin==0){
    double x,y;
    e->GetPoint(selectedBin,x,y);
    return y;
  }else{
    double xLow,yLow;
    e->GetPoint(selectedBin-1,xLow,yLow);
    double xHigh,yHigh;
    e->GetPoint(selectedBin,xHigh,yHigh);
    float fracLow = (pt-xLow)/(xHigh-xLow); // the weight to give to the lower bin
    float weight = yLow*fracLow+yHigh*(1-fracLow);
    return weight;
  }
  return 1;
}
