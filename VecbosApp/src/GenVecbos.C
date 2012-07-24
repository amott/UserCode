#define GenVecbos_cxx
#include "GenVecbos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <Jet.hh>

void GenVecbos::Loop() {}

vector<Jet> GenVecbos::SISCone(vector<TLorentzVector> InputCollection, double Rparam, double thePtMin){
  fastjet::SISConePlugin* mPlugin;

  fastjet::SISConePlugin::SplitMergeScale scale = fastjet::SISConePlugin::SM_pttilde;
  mPlugin = new fastjet::SISConePlugin(Rparam, 0.75, 0, thePtMin, false, scale);

  std::vector<fastjet::PseudoJet> input_vectors;

  int index_=0;
  for(int i = 0; i < InputCollection.size(); i++){
    
    double px = InputCollection[i].Px();
    double py = InputCollection[i].Py();
    double pz = InputCollection[i].Pz();
    double E = InputCollection[i].E();
    fastjet::PseudoJet PsJet(px,py,pz,E);
    PsJet.set_user_index(index_);
    input_vectors.push_back(PsJet);
    index_++;
  }

  vector<Jet> output;
  if(index_ == 0) return output;

  fastjet::JetDefinition jetDefinition(mPlugin);

  fastjet::ClusterSequence clusterSequence(input_vectors, jetDefinition);

  vector<fastjet::PseudoJet> inclusiveJets = clusterSequence.inclusive_jets(1.0);
  
  
  for(std::vector<fastjet::PseudoJet>::const_iterator itJet=inclusiveJets.begin();
      itJet!=inclusiveJets.end();itJet++){
    
    double px = (*itJet).px();
    double py = (*itJet).py();
    double pz = (*itJet).pz();
    double E = (*itJet).E();
    TLorentzVector J(px,py,pz,E);
    
    output.push_back(Jet(J,0.0,0.0));
  }
 

  delete mPlugin;
  return output;
}
