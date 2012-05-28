// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

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
#include "HggVertexing.hh"

HggVertexing::HggVertexing(VecbosBase *b):
  base(b),
  vAna(vtxAlgoParams)
{
  useConversion = false;
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

  rankmethod.clear();
  varNameUsed.clear();
  
  rankmethod.push_back(1); ///from largest to smallest ( from signal like to background like) 
  rankmethod.push_back(1);
  rankmethod.push_back(1);
  
  varNameUsed.push_back("ptbal");
  varNameUsed.push_back("ptasym");
  varNameUsed.push_back("logsumpt2");

  indvertexSelected_allpairpresel  = new vector<short>;
  photontrkisoselvtxdr03 = new vector< vector<float> >;

  //configure vertex analyzer
  vAna.setupWithDefaultOptions(perVtxMvaWeights, perEvtMvaWeights, rankVariables, perVtxReader, perVtxMvaMethod, perEvtReader, perEvtMvaMethod);
  isInit=true;
}


//method to do the vertexing
pair<int,float> HggVertexing::vertex_tmva(VecbosPho *pho1, VecbosPho *pho2){
  
  if(base->nPV == 1) return pair<int,float>(0,1); 
  if(base->nPV == 0) return pair<int,float>(-1,1); 

  vAna.clear();
  if(useConversion){ 
    //cout << "\tMatching Photons to Conversions ... " << flush;
    pho1->matchConversion(base,false); //try to match photons to conversions
    pho2->matchConversion(base,false); //use deta dphi matching (false flag)
    //cout << "Done" << endl;
  }
  //general info
  TVector3 beamspot(base->beamSpotX,base->beamSpotY,base->beamSpotZ);
  
  //fill photon info adapter classes
  //  cout << "\tFilling Photon Info ... " << flush;
  PhotonInfo photonInfo1 = PhotonInfo(pho1->index,
				      pho1->CaloPos,
				      beamspot,
				      pho1->conversion.vtx,
				      pho1->conversion.pRefittedPair,
				      pho1->energy,
				      (int)pho1->isBarrel(),
				      pho1->conversion.vtxNTracks,
				      pho1->conversion.vtxIsValid,
				      pho1->conversion.vtxChi2Prob,
				      pho1->conversion.eOverP
				      );//fill me
  
  PhotonInfo photonInfo2 = PhotonInfo(pho2->index,
				      pho2->CaloPos,
				      beamspot,
				      pho2->conversion.vtx,
				      pho2->conversion.pRefittedPair,
				      pho2->energy,
				      (int)pho2->isBarrel(),
				      pho2->conversion.vtxNTracks,
				      pho2->conversion.vtxIsValid,
				      pho2->conversion.vtxChi2Prob,
				      pho2->conversion.eOverP
				      );//fill me
    
  /*
  cout << pho1->index << endl
       << pho1->energy << endl
       << pho1->isBarrel() << endl
       << pho1->conv()->vtxNTracks << endl
       << pho1->conv()->vtxIsValid << endl
       << pho1->conv()->vtxChi2Prob << endl
       << pho1->conv()->eOverP << endl << endl;
  cout << pho2->index << endl
       << pho2->energy << endl
       << pho2->isBarrel() << endl
       << pho2->conv()->vtxNTracks << endl
       << pho2->conv()->vtxIsValid << endl
       << pho2->conv()->vtxChi2Prob << endl
       << pho2->conv()->eOverP << endl;
  */     
    
  
  //cout << "Done" << endl;
  //It looks like we have to explicitly fill these arrays, gah!
  const int maxPV = 40;
  int nPV;
  float PVxPV[maxPV],PVyPV[maxPV],PVzPV[maxPV];
  nPV = base->nPV;
  if(nPV > maxPV) nPV = maxPV; //safety!
  for(int i=0;i<nPV;i++){
    PVxPV[i] = base->PVxPV[i];
    PVyPV[i] = base->PVyPV[i];
    PVzPV[i] = base->PVzPV[i];
  }

  const int maxTracks = 1000;
  int nTracks;
  float pxTrack[maxTracks],pyTrack[maxTracks],pzTrack[maxTracks];
  float ptErrorTrack[maxTracks],vtxWeightTrack[maxTracks];
  float d0Track[maxTracks],d0ErrorTrack[maxTracks],dzTrack[maxTracks],dzErrorTrack[maxTracks];
  int vtxIndexTrack[maxTracks];

  nTracks = base->nTrack;
  if(nTracks> maxTracks) return pair<int,float>(0,0); //safety!
  if(nTracks> maxTracks) nTracks = maxTracks; //safety!
  for(int i=0;i<nTracks;i++){
    pxTrack[i] = base->pxTrack[i]; pyTrack[i] = base->pyTrack[i]; pzTrack[i] = base->pzTrack[i]; 
    ptErrorTrack[i] = base->ptErrorTrack[i]; vtxIndexTrack[i] = base->vtxIndexTrack[i]; vtxWeightTrack[i] = base->vtxWeightTrack[i]; 
    d0Track[i] = base->d0Track[i]; d0ErrorTrack[i] = base->d0ErrorTrack[i]; dzTrack[i] = base->dzTrack[i]; dzErrorTrack[i] = base->dzErrorTrack[i]; 
  }
  bool highPurityTrack[nTracks];
  //get High Purity bit from quality flag
  //cout << "\tGetting Track Purity ... " << flush;
  for(int iTrack=0; iTrack<nTracks; iTrack++){
    //purity flag
    int thisQuality = base->qualityMaskTrack[iTrack];
    const int highPurityFlag = 3;
    highPurityTrack[iTrack] = (bool)( thisQuality & 1 << highPurityFlag); // applies the mask 1000 to the quality mask 
  }
  
  TupleVertexInfo vertexInfoAdapter = TupleVertexInfo( nPV,
						       &(PVxPV[0]),
						       &(PVyPV[0]),
						       &(PVzPV[0]),
						       nTracks,
						       &(pxTrack[0]),
						       &(pyTrack[0]),
						       &(pzTrack[0]),
						       &(ptErrorTrack[0]),
						       &(vtxIndexTrack[0]),
						       &(vtxWeightTrack[0]),
						       &(d0Track[0]), //add this
						       &(d0ErrorTrack[0]), //add this
						       &(dzTrack[0]), //add thiss
						       &(dzErrorTrack[0]), //add this
						       &highPurityTrack[0]
						       );//fill me
  
  
  //cout<<"tmva vertex ranking....\n";
  //cout << "\tAnalyzing ... " << flush;
  vAna.analyze(vertexInfoAdapter,photonInfo1,photonInfo2);     
  //cout << "Done" << endl;
  /// rank product vertex selection. Including pre-selection based on conversions information.
  vector<int> rankprod = vAna.rankprod(rankVariables);
  
  /// MVA vertex selection
  vector<int> vtx_ranked_tmva = vAna.rank(*perVtxReader,perVtxMvaMethod);

  return pair<int,float>(vtx_ranked_tmva[0],vAna.mva(vtx_ranked_tmva[0]));  
}


