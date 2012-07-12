// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "ReadConfig.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "HggVertexingNew.hh"

#define debugVertexing 0

HggVertexing::HggVertexing(VecbosBase *b):
  base(b)
{
  useConversion = true;
  isInit=false;
}

 
HggVertexing::~HggVertexing() {}

void HggVertexing::init(){
  //read the config file to setup TMVA
  cout << "Doing Vertexing Initialization" << endl;
  ReadConfig cfg;
  int errorCode = cfg.read(configFilePath);
  cout << "ReadConfig exited with code: " << errorCode << endl;
  if(errorCode) return;
  cfg.printAll();
  perVtxMvaWeights = cfg.getParameter("perVtxMvaWeights");
  perVtxMvaMethod  = cfg.getParameter("perVtxMvaMethod");
  perEvtMvaWeights = cfg.getParameter("perEvtMvaWeights");
  perEvtMvaMethod  = cfg.getParameter("perEvtMvaMethod");

  perVtxReader = new TMVA::Reader( "!Color:!Silent" );
  perEvtReader = new TMVA::Reader( "!Color:!Silent" );

  PerVtxVars = new float[5];
  PerEvtVars = new float[8];

  perVtxReader->AddVariable( "ptbal", &PerVtxVars[0] );
  perVtxReader->AddVariable( "ptasym", &PerVtxVars[1] );
  perVtxReader->AddVariable( "logsumpt2", &PerVtxVars[2] );
  perVtxReader->AddVariable( "limPullToConv", &PerVtxVars[3] );
  perVtxReader->AddVariable( "nConv", &PerVtxVars[4] );
  perVtxReader->BookMVA( perVtxMvaMethod, perVtxMvaWeights);

  perEvtReader->AddVariable( "diphoPt0", &PerEvtVars[0] );
  perEvtReader->AddVariable( "nVert", &PerEvtVars[1] );
  perEvtReader->AddVariable( "MVA0", &PerEvtVars[2] );
  perEvtReader->AddVariable( "MVA1", &PerEvtVars[3] );
  perEvtReader->AddVariable( "dZ1", &PerEvtVars[4] );
  perEvtReader->AddVariable( "MVA2", &PerEvtVars[5] );
  perEvtReader->AddVariable( "dZ2", &PerEvtVars[6] );
  perEvtReader->AddVariable( "nConv", &PerEvtVars[7] );
  perEvtReader->BookMVA(perEvtMvaMethod, perEvtMvaWeights);
  
  isInit=true;
}

float HggVertexing::evalPerVtxMVA(float ptbal, float ptasym, float logsumpt2, float limPullToConv, float nConv){

  PerVtxVars[0] = ptbal;
  PerVtxVars[1] = ptasym;
  PerVtxVars[2] = logsumpt2;
  PerVtxVars[3] = limPullToConv;
  PerVtxVars[4] = nConv;

  return perVtxReader->EvaluateMVA(perVtxMvaMethod);
}


std::vector<std::pair<int,float> > HggVertexing::evalPerVtxMVA(VecbosPho* pho1, VecbosPho* pho2){
  if(debugVertexing) cout << "evalPerVtxMVA" << endl;
  double *ptbal          = new double[base->nPV];
  double *ptasym         = new double[base->nPV];
  double *sumpt2         = new double[base->nPV];
  double *limPullToConv  = new double[base->nPV];

  //lists of the 4 vectors from each vertex
  std::vector<TLorentzVector> pho1_fromVtx;
  std::vector<TLorentzVector> pho2_fromVtx;
  //loop over the vertices
  for(int iVtx=0; iVtx<base->nPV; iVtx++){
    TVector3 thisVtxPos(base->PVxPV[iVtx],base->PVyPV[iVtx],base->PVzPV[iVtx]);

    TLorentzVector pho1_p4 = pho1->p4FromVtx(thisVtxPos,pho1->finalEnergy);
    TLorentzVector pho2_p4 = pho2->p4FromVtx(thisVtxPos,pho2->finalEnergy);

    pho1_fromVtx.push_back(pho1_p4);
    pho2_fromVtx.push_back(pho2_p4);
    
    //clear
    ptbal      [iVtx] = 0.0;
    ptasym     [iVtx] = 0.0;
    sumpt2     [iVtx] = 0.0;
    limPullToConv [iVtx] = -1.0; 

  }//for(int iVtx=0; iVtx<base->nPV; iVtx++)


  std::vector<TLorentzVector> trackMomentum;
  for(int iVtx=0; iVtx<base->nPV; iVtx++) trackMomentum.push_back(TLorentzVector(0.,0.,0.,0.));

  for(int iTrk=0; iTrk<base->nTrack; iTrk++){
    int iVtx = base->vtxIndexTrack[iTrk];
    if(iVtx<0 || iVtx >= base->nPV) continue;
    TLorentzVector thisMomentum; //use M=0
    thisMomentum.SetPxPyPzE(base->pxTrack[iTrk],base->pyTrack[iTrk],base->pzTrack[iTrk],
			    TMath::Sqrt(TMath::Power(base->pxTrack[iTrk],2)+
					TMath::Power(base->pyTrack[iTrk],2)+
					TMath::Power(base->pzTrack[iTrk],2)));
    if(thisMomentum.Pt() <0.001) continue;
    sumpt2[iVtx]+= thisMomentum.Pt()*thisMomentum.Pt();

    //if(debugVertexing) cout << "This Track " << iTrk << "  Pt: " << thisMomentum.Pt() << endl;


    if(thisMomentum.DeltaR(pho1_fromVtx[iVtx]) < 0.05 ||
       thisMomentum.DeltaR(pho2_fromVtx[iVtx]) < 0.05 ) continue;

    trackMomentum[iVtx] += thisMomentum;
    
  }//for(int iTrk=0; iTrk<base->nTrack; iTrk++)

  std::pair<float,float> convZ = getZConv(pho1,pho2); //z,dz
  int maxVertices = ( (pho1_fromVtx[0] + pho2_fromVtx[0]).Pt() > 30 ? 3 : 5);
  double minDz = 999;


  std::vector<std::pair<int,float> > vertexMVAs;
  for(int iVtx=0; iVtx<base->nPV; iVtx++){ //loop over the vertices AGAIN to compute
    //this is inelegant, but far quicker, rather than looping over the track collection for
    //every vertex

    TLorentzVector thisTrackMom = trackMomentum[iVtx];
    TLorentzVector thisHiggsMom = pho1_fromVtx[iVtx] + pho2_fromVtx[iVtx];

    ptbal[iVtx] =  thisTrackMom.Px()*thisHiggsMom.Px() + thisTrackMom.Py()*thisHiggsMom.Py();
    ptbal[iVtx] = -ptbal[iVtx]/thisHiggsMom.Pt();

    ptasym[iVtx] = (thisTrackMom.Pt() - thisHiggsMom.Pt())/(thisTrackMom.Pt() + thisHiggsMom.Pt());

    //if there is a conversion, compute the pull
    if(pho1->conversion.index != -1 || pho2->conversion.index != -1){
      limPullToConv[iVtx] = TMath::Abs(base->PVzPV[iVtx]-convZ.first)/convZ.second;
    }
    
    double mva = evalPerVtxMVA(ptbal[iVtx],ptasym[iVtx],log(sumpt2[iVtx]),
			       limPullToConv[iVtx],getNConv(pho1,pho2));
    vertexMVAs.push_back(std::pair<int,float>(iVtx,mva));
  }

  //hmmm std::sort doesn't appear to work
  //std::sort(vertexMVAs.begin(),vertexMVAs.end(),sort_pred()); //sort by MVA value

  delete ptbal;
  delete ptasym;
  delete sumpt2;
  delete limPullToConv;

  if(debugVertexing) cout << "done with perVtxMVA" << endl;
  return vertexMVAs;
}


float HggVertexing::evalPerEvtMVA(VecbosPho* pho1,VecbosPho* pho2, 
				  std::vector<std::pair<int,float> > *perVertexRank){
  if(debugVertexing) cout << "PerEvtMVA" << endl;
  int bestVtx = perVertexRank->at(0).first;
  if(bestVtx <0 || bestVtx >=base->nPV) return -1e6;
  TVector3 bestVtxPos(base->PVxPV[bestVtx],base->PVyPV[bestVtx],base->PVzPV[bestVtx]);
    
  TLorentzVector pho1_p4 = pho1->p4FromVtx(bestVtxPos,pho1->finalEnergy);
  TLorentzVector pho2_p4 = pho2->p4FromVtx(bestVtxPos,pho2->finalEnergy);
  TLorentzVector higgs_p4 = pho1_p4+pho2_p4;

  if(perVertexRank->size()==0) return -1e6;
  while(perVertexRank->size()<3) perVertexRank->push_back(std::pair<int,float>(-1,-1e6));

  PerEvtVars[0] = higgs_p4.Pt();
  PerEvtVars[1] = base->nPV;
  PerEvtVars[2] = perVertexRank->at(0).second;
  PerEvtVars[3] = perVertexRank->at(1).second;
  PerEvtVars[4] = (perVertexRank->at(1).first!=-1 ? base->PVzPV[perVertexRank->at(1).first] - bestVtxPos.Z() : 0);
  PerEvtVars[5] = perVertexRank->at(2).second;
  PerEvtVars[6] = (perVertexRank->at(2).first!=-1 ? base->PVzPV[perVertexRank->at(2).first] - bestVtxPos.Z() : 0);
  PerEvtVars[7] = getNConv(pho1,pho2);

  if(debugVertexing) cout << "PerEvtMVA --Done" << endl;
  
  return perEvtReader->EvaluateMVA( perEvtMvaMethod );
}



//method to do the vertexing
vector<pair<int,float> > HggVertexing::vertex_tmva(VecbosPho *pho1, VecbosPho *pho2,float& evtMVA){
  if(debugVertexing) cout << "vertexTMVA" << endl;
  
  //if(base->nPV == 1) return pair<int,float>(0,1); 
  if(base->nPV == 0) return std::vector<pair<int,float> >();

  std::vector<std::pair<int,float> > vertexRanks = evalPerVtxMVA(pho1,pho2);
  evtMVA = evalPerEvtMVA(pho1,pho2,&vertexRanks);
  if(debugVertexing) cout << "Done --VertexTMVA" << endl;

  return vertexRanks;
}

std::pair<float,float> HggVertexing::getZConv(VecbosPho* pho1,VecbosPho* pho2){

  std::pair<float,float> pho1_z = getZConv(pho1);
  std::pair<float,float> pho2_z = getZConv(pho2);

  float z1 = pho1_z.first;
  float dz1 = pho1_z.second;
  float z2 = pho2_z.second;
  float dz2 = pho2_z.second;

  float zconv=0;
  float dzconv=0;
  if(pho1->conversion.index!=-1){
    zconv = z1; dzconv = dz1;
  }
  if(pho2->conversion.index!=-1){
    zconv = z2; dzconv = dz2;
  }
  if(pho1->conversion.index!=-1 && pho2->conversion.index!=-1){
    float zconv = (z1/dz1/dz1 + z2/dz2/dz2)/(1./dz1/dz1 + 1./dz2/dz2 );
    float dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
  }
  return std::pair<float,float>(zconv,dzconv);
}

//copied from MIT framework
 std::pair<float,float> HggVertexing::getZConv(VecbosPho* pho){
   const double dzpxb = 0.016;
   const double dztib = 0.331;
   const double dztob = 1.564;
   const double dzpxf = 0.082;
   const double dztid = 0.321;
   const double dztec = 0.815;

   const double dzpxbsingle = 0.036;
   const double dztibsingle = 0.456;
   const double dztobsingle = 0.362;
   const double dzpxfsingle = 0.130;
   const double dztidsingle = 0.465;
   const double dztecsingle = 1.018;

   double zconv  = -99;
   double dzconv = -99;

   VecbosConversion c = pho->conversion;
   if(c.index == -1) return std::pair<float,float>(0.,0.);


  TVector3 beamSpot(base->beamSpotX,base->beamSpotY,base->beamSpotZ);
  double zconvtrk = convCorrectedDz(&c,beamSpot) + beamSpot.Z();
  double zconvsc  = Z0EcalVtxCiC(&c,beamSpot,pho->SC.CaloPos);

   if( c.vtxNTracks == 2){
     if(pho->SC.CaloPos.Eta() < 1.5){
       double rho = c.CaloPos.Mag();
       if(rho<15)           { dzconv = dzpxb; zconv = zconvtrk; }
       else if( rho < 60. ) { dzconv = dztib; zconv = zconvsc; }
       else                 { dzconv = dztob; zconv = zconvsc; }
     }else{ 
       double z = c.CaloPos.Z();
       if     ( TMath::Abs(z) < 50. ) { dzconv = dzpxf; zconv = zconvtrk; }
       else if( TMath::Abs(z) < 100.) { dzconv = dztid; zconv = zconvtrk; }
       else                           { dzconv = dztec; zconv = zconvsc; }
     }
   }//if( c.vtxNTrack == 2)
   else if( c.vtxNTracks == 1){
     if(pho->SC.CaloPos.Eta() < 1.5){
       double rho = c.CaloPos.Mag();
       if(rho<15)           { dzconv = dzpxbsingle; zconv = zconvsc; }
       else if( rho < 60. ) { dzconv = dztibsingle; zconv = zconvsc; }
       else                 { dzconv = dztobsingle; zconv = zconvsc; }
     }else{
       double z = c.CaloPos.Z();
       if     ( TMath::Abs(z) < 50. ) { dzconv = dzpxfsingle; zconv = zconvsc; }
       else if( TMath::Abs(z) < 100.) { dzconv = dztidsingle; zconv = zconvsc; }
       else                           { dzconv = dztecsingle; zconv = zconvsc; }
     }     
   }

   return std::pair<float,float>(zconv,dzconv);

}

 float HggVertexing::convCorrectedDz(VecbosConversion* c,TVector3 basePos){
   TVector3 momPerp(c->pRefittedPair.Px(),c->pRefittedPair.Py(),0);
   TVector3 posPerp(c->CaloPos.X()-basePos.X(),c->CaloPos.Y()-basePos.Y(),0);

   return c->CaloPos.Z() - basePos.Z() - posPerp.Dot(momPerp)/momPerp.Pt() * (c->pRefittedPair.Pt()/momPerp.Pt());
 }

float HggVertexing::Z0EcalVtxCiC(VecbosConversion* c, TVector3 basePos,TVector3 caloPos){
   TVector3 dirscvtx = caloPos - c->CaloPos;
   TVector3 momPerp(c->pRefittedPair.Px(),c->pRefittedPair.Py(),0);
   TVector3 posPerp(c->CaloPos.X()-basePos.X(),c->CaloPos.Y()-basePos.Y(),0);

   return c->CaloPos.Z() - posPerp.Mag() * (dirscvtx.Z()/dirscvtx.Mag());
   
 }
