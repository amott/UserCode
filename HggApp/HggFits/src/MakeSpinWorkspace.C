#include "MakeSpinWorkspace.h"

using namespace std;
#include "selectionMaps.C"

MakeSpinWorkspace::MakeSpinWorkspace(TString outputFileName){
  requireCiC=true;
  tightPt=false;
  selectionMap=0;

  runLow=0;
  runHigh= 9999999;

  useR9=false;

  outputFile = TFile::Open(outputFileName,"RECREATE");
  ws = new RooWorkspace();
  ws->SetName("cms_hgg_spin_workspace");  
  
}

MakeSpinWorkspace::~MakeSpinWorkspace(){
  if(outputFile){
    outputFile->cd();
    ws->Write();
    outputFile->Close();
  }
}

void MakeSpinWorkspace::runOnAll( void (*f)(TString,TString), TString mcName){
  for(int i=0;i<nCat;i++){
    for(int j=0;j<nCat;j++){
      f(Form("EB_%d_%d",i,j),mcName);
      f(Form("EE_%d_%d",i,j),mcName);
    }
  }
}

void MakeSpinWorkspace::runOnAllR9( void (*f)(TString,TString), TString mcName){
  for(int i=0;i<nCat;i++){
      f(Form("EB_%d",i),mcName);
      f(Form("EE_%d",i),mcName);
  }
}


TH2F* MakeSpinWorkspace::getSelectionMap(int map,bool isData){
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
bool MakeSpinWorkspace::getBaselineSelection(HggOutputReader2* h,int maxI,int minI){
  if(h->nPhoton < 2
     || h->diPhotonMVA<-1
     || h->mPair < 100
     || h->mPair > 180
     || h->Photon_pt[minI] < 24
     || h->Photon_pt[maxI] < 32) return false;
  return true;
}

void MakeSpinWorkspace::AddToWorkspace(TString inputFile,TString tag, bool isData){
  //gROOT->ProcessLine(".L /home/amott/HggApp/spin/LeptonSPlots/src/HggOutputReader2.C+");
  cout << "saveDataSet" <<endl;
  TFile *f = TFile::Open(inputFile);
  TTree *tree = (TTree*)f->Get("HggOutput");
  
  HggOutputReader2 h(tree);
  
  std::cout << "Get Selection Map" << std::endl;
  TH2F* map = getSelectionMap(selectionMap,isData);
  std::cout << map->GetBinContent(1,1) << std::endl;

  float **offset=0,**scale=0;
  // fit variables
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

  RooArgSet set;
  set.add(*mass);
  set.add(*cosT);
  set.add(*sige1);
  set.add(*sige2);
  set.add(*evtW);
  set.add(*eta1); set.add(*eta2);
  set.add(*etaSC1); set.add(*etaSC2);
  set.add(*phi1); set.add(*phi2);
  set.add(*pt1); set.add(*pt2);
  set.add(*r91); set.add(*r92);
  set.add(*idMVA1); set.add(*idMVA2);
  set.add(*diPhotonMVA);

  RooRealVar *totEB = new RooRealVar(Form("%s_EB_totalEvents",tag.Data()),"",0,0,1e9);
  RooRealVar *totEE = new RooRealVar(Form("%s_EE_totalEvents",tag.Data()),"",0,0,1e9);

  std::map<std::pair<int,int>, RooDataSet*> dataMapEB, dataMapEE;
  std::map<std::pair<int,int>, RooDataSet*> *datamap;

  if(useR9){
    for(int j=0;j<nCat;j++){
      cat->defineType( Form("EB_%d",j),2*j );
      cat->defineType( Form("EE_%d",j),2*j+1 );
      dataMapEB[std::pair<int,int>(j,0)] = new RooDataSet(Form("%s_EB_%d",tag.Data(),j),"",set,"evtWeight");
      dataMapEE[std::pair<int,int>(j,0)] = new RooDataSet(Form("%s_EE_%d",tag.Data(),j),"",set,"evtWeight");
    }
  }else{
    for(int i=0;i<nCat;i++){
      for(int j=0;j<nCat;j++){
	cat->defineType( Form("EB_%d_%d",i,j),4*i+2*j );
	cat->defineType( Form("EE_%d_%d",i,j), 4*i+2*j+1 );
	dataMapEB[std::pair<int,int>(i,j)] = new RooDataSet(Form("%s_EB_%d_%d",tag.Data(),i,j),"",set,"evtWeight");
	dataMapEE[std::pair<int,int>(i,j)] = new RooDataSet(Form("%s_EE_%d_%d",tag.Data(),i,j),"",set,"evtWeight");
      }
    }
  }

  Long64_t iEntry=-1;
  cout << "Making DataSet" << endl;
  Long64_t nEB=0,nEE=0;
  while(h.GetEntry(++iEntry)){
    if( !(iEntry%10000) ) cout << "Processing " << iEntry << "\t\t" << h.evtWeight << endl;
    int maxI = (h.Photon_pt[1] > h.Photon_pt[0] ? 1:0);
    int minI = (h.Photon_pt[1] > h.Photon_pt[0] ? 0:1);
    if(!getBaselineSelection(&h,maxI,minI)) continue;
    if(tightPt && (h.Photon_pt[maxI]/h.mPair < 1./3. || h.Photon_pt[minI]/h.mPair < 1./4.) ) continue;

    if(isData && (h.runNumber < runLow || h.runNumber > runHigh) ) continue;

    if(fabs(h.Photon_etaSC[1]) < 1.48 && fabs(h.Photon_etaSC[0]) < 1.48) nEB+=h.evtWeight;
    else nEE+=h.evtWeight;

    if(requireCiC){
      if(h.Photon_passPFCiC[1]==false || h.Photon_passPFCiC[0]==false) continue;
    }

    float se1 = h.Photon_EError[maxI]/h.Photon_E[maxI];
    float se2 = h.Photon_EError[minI]/h.Photon_E[minI];

    int p1,p2;
    if(useR9){
      p1 = passSelection(h.Photon_r9[maxI]);
      p2 = passSelection(h.Photon_r9[minI]);
    }else{
      p1 = passSelection(map,se1,h.Photon_etaSC[maxI],h.Photon_pt[maxI]);
      p2 = passSelection(map,se2,h.Photon_etaSC[minI],h.Photon_pt[minI]);
    }
    if(p1 >= nCat || p2 >= nCat) continue;
    datamap = ((fabs(h.Photon_etaSC[maxI]) < 1.48 && fabs(h.Photon_etaSC[minI]) < 1.48) ?
	       &dataMapEB : & dataMapEE);
      
    mass->setVal(h.mPair);
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

    if(!isData) evtW->setVal(h.evtWeight);
    else evtW->setVal(1);
    if( !(iEntry%10000) ) cout <<  "\t\t\t" << evtW->getVal() << endl;

    int cat = (p1==0 && p2==0 ? 0 : 1);
    (*datamap)[std::pair<int,int>(cat,0)]->add(set);
  }
  cout << "Processed " << iEntry << " Entries" <<endl;

  totEB->setVal(nEB);
  totEE->setVal(nEE);

  RooArgSet setCat(set);
  setCat.add(*cat);
  RooDataSet* dataComb = new RooDataSet(tag+"_Combined","",setCat );
  
  std::map<std::pair<int,int>, RooDataSet*>::iterator dIt;
  for(dIt = dataMapEB.begin();dIt!=dataMapEB.end();dIt++){
    TString cattag;
    if(useR9) cattag = Form("EB_%d", (dIt->first).first );
    else cattag = Form("EB_%d_%d", (dIt->first).first, (dIt->first).second );
    std::cout << cattag <<std::endl;
    RooDataSet *tmp = new RooDataSet("DataCat_"+cattag,"",set,RooFit::Index(*cat),RooFit::Import(cattag,*(dIt->second)) );
    dataComb->append(*tmp);
    ws->import(*(dIt->second));
  }
  for(dIt = dataMapEE.begin();dIt!=dataMapEE.end();dIt++){
    TString cattag;
    if(useR9) cattag = Form("EE_%d", (dIt->first).first );
    else cattag = Form("EE_%d_%d", (dIt->first).first, (dIt->first).second );
    std::cout << cattag <<std::endl;
    RooDataSet *tmp = new RooDataSet("DataCat_"+cattag,"",set,RooFit::Index(*cat),RooFit::Import(cattag,*(dIt->second)) );
    dataComb->append(*tmp);    
    ws->import(*(dIt->second));
  }
  ws->import(*dataComb);
  ws->import(*totEB);
  ws->import(*totEE);
  cout << "Done" <<endl;
  if(offset) delete offset;
  if(scale)  delete scale;
}

void MakeSpinWorkspace::addFile(TString fName,TString l,bool is){
  fileName.push_back(fName);
  label.push_back(l);
  isData.push_back(is);
}


void MakeSpinWorkspace::MakeWorkspace(){

  std::vector<TString>::const_iterator fnIt,lIt;
  std::vector<bool>::const_iterator idIt;

  fnIt = fileName.begin();
  lIt  = label.begin();
  idIt = isData.begin();
  for(; fnIt != fileName.end(); fnIt++, lIt++, idIt++){
    AddToWorkspace(*fnIt,*lIt,*idIt);
  }
  
  outputFile->cd();
  ws->Write();
  outputFile->Close();
}
