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
}

 
HggVertexing::~HggVertexing() {}

void HggVertexing::init(){
  //read the config file to setup TMVA
  ReadConfig cfg(configFilePath);
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

}


//method to do the vertexing
int HggVertexing::vertex_tmva(int ipho1, int ipho2){
  
  if(base->nPV == 1) return 0; 
  if(base->nPV == 0) return -1; 
  
  VecbosPho pho1(base,ipho1);
  VecbosPho pho2(base,ipho2);

  if(useConversion){ 
    pho1.matchConversion(base,false); //try to match photons to conversions
    pho2.matchConversion(base,false); //use deta dphi matching (false flag)
  }
  //general info
  TVector3 beamspot(base->beamSpotX,base->beamSpotY,base->beamSpotZ);
  
  //fill photon info adapter classes
  PhotonInfo photonInfo1 = PhotonInfo(ipho1,
				      pho1.CaloPos,
				      beamspot,
				      pho1.conversion->vtx,
				      pho1.conversion->pRefittedPair,
				      pho1.energy,
				      (int)pho1.isBarrel(),
				      pho1.conversion->vtxNTracks,
				      pho1.conversion->vtxIsValid,
				      pho1.conversion->vtxChi2Prob,
				      pho1.conversion->eOverP
				      );//fill me
  
  PhotonInfo photonInfo2 = PhotonInfo(ipho2,
				      pho2.CaloPos,
				      beamspot,
				      pho2.conversion->vtx,
				      pho2.conversion->pRefittedPair,
				      pho2.energy,
				      (int)pho2.isBarrel(),
				      pho2.conversion->vtxNTracks,
				      pho2.conversion->vtxIsValid,
				      pho2.conversion->vtxChi2Prob,
				      pho2.conversion->eOverP
				      );//fill me
  
    
    
  
  bool highPurityTrack[base->nTrack];
  //get High Purity bit from quality flag
  for(int iTrack=0; iTrack<base->nTrack; iTrack++){
    //purity flag
    int thisQuality = base->qualityMaskTrack[iTrack];
    const int highPurityFlag = 3;
    highPurityTrack[iTrack] = (bool)( thisQuality & 1 << highPurityFlag); // applies the mask 1000 to the quality mask 
  }

  TupleVertexInfo vertexInfoAdapter = TupleVertexInfo( base->nPV,
						       &(base->PVxPV[0]),
						       &(base->PVyPV[0]),
						       &(base->PVzPV[0]),
						       base->nTrack,
						       &(base->pxTrack[0]),
						       &(base->pyTrack[0]),
						       &(base->pzTrack[0]),
						       &(base->ptErrorTrack[0]),
						       &(base->vtxIndexTrack[0]),
						       &(base->vtxWeightTrack[0]),
						       &(base->d0Track[0]), //add this
						       &(base->d0ErrorTrack[0]), //add this
						       &(base->dzTrack[0]), //add thiss
						       &(base->dzErrorTrack[0]), //add this
						       &highPurityTrack[0]
						       );//fill me
  
  
  //cout<<"tmva vertex ranking....\n";
  vAna.analyze(vertexInfoAdapter,photonInfo1,photonInfo2);     
 
  /// rank product vertex selection. Including pre-selection based on conversions information.
  vector<int> rankprod = vAna.rankprod(rankVariables);
  
  /// MVA vertex selection
  vector<int> vtx_ranked_tmva = vAna.rank(*perVtxReader,perVtxMvaMethod);

  return vtx_ranked_tmva[0];  
}


