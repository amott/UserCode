/// std includes
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/// ROOT includes
#include <TTree.h>
#include <TH1D.h>
#include <TVector3.h>
#include <TLorentzVector.h>

/// local includes
#include "RSZZAnalysis.hh"
#include "Plane.hh"
#include "MaxFlow.hh"
#include "Comparator.hh"

// Declare pointer to data as global (not elegant but TMinuit needs this).
vector<TVector2>* TracksPtr;
// HACK
double weight = 1.0;

RSZZAnalysis::RSZZAnalysis(TTree *tree) : Vecbos(tree) {
  InitParameters();
  caloThresholds.reserve(10);
  caloThresholds = DefaultCaloThresholds();
}

RSZZAnalysis::~RSZZAnalysis(){
  
}

void RSZZAnalysis::Loop(string outfile) {
  
  if(fChain == 0) return;

  // Prepare the pdgId vectors.
  vector<int> GIds;
  GIds.push_back(5000039);
  vector<int> ZIds;
  ZIds.push_back(23);
  vector<int> excludedForJets;
  excludedForJets.reserve(6);
  excludedForJets.push_back(12);
  excludedForJets.push_back(-12);
  excludedForJets.push_back(14);
  excludedForJets.push_back(-14);
  excludedForJets.push_back(16);
  excludedForJets.push_back(-16);
  if(verbose) cout << "pdgId ready" << endl;

//   // Create sets of histograms.
//   histosJet1 = CreateJetHistos("jet1");
//   histosJet2 = CreateJetHistos("jet2");
//   histos2DJet1 = Create2DHistos("jet1");
//   histos2DJet2 = Create2DHistos("jet2");
//   histosJet1_AC = CreateJetHistos("jet1_AC");
//   histosJet2_AC = CreateJetHistos("jet2_AC");
//   histos2DJet1_AC = Create2DHistos("jet1_AC");
//   histos2DJet2_AC = Create2DHistos("jet2_AC");
//   if(verbose) cout << "Sets of histograms ready" << endl;
    
//   // Individual histograms.
//   flowJet1 = new TH1D("flowJet1", "flowJet1", 100, 0., 1.);
//   antiFlowJet1 = new TH1D("antiFlowJet1", "antiFlowJet1", 100, 0., 1.);
//   aFOverFJet1 = new TH1D("aFOverFJet1", "aFOverFJet1", 100, 0., 1.);
//   flowJet1_AC = new TH1D("flowJet1_AC", "flowJet1_AC", 100, 0., 1.);
//   antiFlowJet1_AC = new TH1D("antiFlowJet1_AC", "antiFlowJet1_AC", 100, 0., 1.);
//   aFOverFJet1_AC = new TH1D("aFOverFJet1_AC", "aFOverFJet1_AC", 100, 0., 1.);
//   flowJet2 = new TH1D("flowJet2", "flowJet2", 100, 0., 1.);
//   antiFlowJet2 = new TH1D("antiFlowJet2", "antiFlowJet2", 100, 0., 1.);
//   aFOverFJet2 = new TH1D("aFOverFJet2", "aFOverFJet2", 100, 0., 1.);
//   flowJet2_AC = new TH1D("flowJet2_AC", "flowJet2_AC", 100, 0., 1.);
//   antiFlowJet2_AC = new TH1D("antiFlowJet2_AC", "antiFlowJet2_AC", 100, 0., 1.);
//   aFOverFJet2_AC = new TH1D("aFOverFJet2_AC", "aFOverFJet2_AC", 100, 0., 1.);
  
  H_jm1Xjm2 = new TH2D("jm1Xjm2", "jm1Xjm2", 1000, 0., 1000., 1000, 0., 1000.);
  H_j1_manytracks_et = new TH1D("j1_manytracks_et", "j1_manytracks_et", 5000, 0., 5000.);
  H_j1_manytracks_eta = new TH1D("j1_manytracks_eta", "j1_manytracks_eta", 100, -5., 5.);
  H_j2_manytracks_et = new TH1D("j2_manytracks_et", "j2_manytracks_et", 5000, 0., 5000.);
  H_j2_manytracks_eta = new TH1D("j2_manytracks_eta", "j2_manytracks_eta", 100, -5., 5.);
    
  
  if(verbose) cout << "Histograms ready" << endl;

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  Long64_t nentries = fChain->GetEntries();
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries ;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(verbose) {
      cout << ">>>>>> Processing event # " << jentry << endl;
    }
    if (jentry%1000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    
    // Substituted this part for a VecBos function that returns the "default"
    // calorimetric thresholds. In case the defaults are not wanted, one could just
    // uncomment the code below and fill in new numbers. See the constructor.
    //     caloThresholds.clear();
    //     caloThresholds.push_back(.9);
    //     caloThresholds.push_back(1.1);
    //     caloThresholds.push_back(1.4);
    //     caloThresholds.push_back(1.4);
    //     caloThresholds.push_back(1.2);
    //     caloThresholds.push_back(1.8);
    //     caloThresholds.push_back(.09);
    //     caloThresholds.push_back(.45);
    //     caloThresholds.push_back(.2);
    //     caloThresholds.push_back(.45);

    // I don't want to spend too much time running this, so I use the
    // already computed jets (IC5) to make a preliminary cut.
    if(nSisConeJet < minNumJets) continue;

    // First cut - we want at least ONE primary vertex in the event.
    if(nPV < 1) continue;

    // Get the hardest primary vertex.
    int iPV = GetHardestPV(); 
   
//     // MC information.
//     vector<TLorentzVector> theGraviton = ParticlesFromMcWithId(3, GIds);
//     vector<TLorentzVector> theZs = ParticlesFromMcWithId(3, ZIds);
//     vector<TLorentzVector> particlesForJets = ParticlesFromMcWithNotId(1, excludedForJets);
//     //vector<TLorentzVector> chargedParticles = ParticlesFromMcCharged();

//     // Make the class aware of these objects.
//     graviton = &theGraviton.front();
//     if(theZs.at(0).Pt() > theZs.at(1).Pt()) {
//       Z1 = &theZs.at(0);
//       Z2 = &theZs.at(1);
//     }
//     else {
//       Z1 = &theZs.at(1);
//       Z2 = &theZs.at(0);
//     }
//     //pCharged = &chargedParticles;
//     pCharged = 0;

    if(verbose) cout << "MC info ready" << endl;

    if(verbose) cout << "Preparing CaloTowers... " << nCaloTowers << endl;

    // Two vectors of CaloTowers, with different settings.
    //vector<CaloTower> defaultCaloTowers = CreateCaloTowers(caloThresholds, 0.0, 0); // CMSSW default.
    vector<CaloTower> goodCaloTowers = CreateCaloTowers(caloThresholds, PVzPV[iPV], 3); // Good settings from Chris.
    
    if(verbose) cout << "CaloTowers ready: " << goodCaloTowers.size() << endl;

    // Make jets, cut on them, and sort them by Et.
    vector<Jet> goodJets = SortJetByEt(FilterJetsByEt(SISCone(goodCaloTowers, jetRadius, 0.0), jetEtCut));
    // First cut: minNumJets;
    if(goodJets.size() < minNumJets) continue;

    if(verbose) cout << "Jets ready" << endl;
    
    // The tracks in the event, with the agreed pT cut.
    vector<TLorentzVector> goodTracks = Tracks(trackPtCut);

    if(verbose) cout << "Tracks ready" << endl;

    // Make the class aware of these objects.
    pGoodJets   = &goodJets;
    pGoodTracks = &goodTracks;

    // Create the jets with associations.
    Jet jet1 = goodJets.at(0);
    Jet jet2 = goodJets.at(1);
    vector<TLorentzVector> tracks1 = CloseInEtaPhi(goodTracks, jet1.Get4Vector(), jetRadius);
    vector<TLorentzVector> tracks2 = CloseInEtaPhi(goodTracks, jet2.Get4Vector(), jetRadius);
    JetWithAssociation jetTracks1 = make_pair(jet1, tracks1);
    JetWithAssociation jetTracks2 = make_pair(jet2, tracks2);
          
    if(verbose) cout << "Associations ready" << endl;

    H_jm1Xjm2->Fill(jet1.mass(), jet2.mass());
    if(tracks1.size() > 120) {
      H_j1_manytracks_et->Fill(jet1.et());
      H_j1_manytracks_eta->Fill(jet1.eta());
    }
    if(tracks2.size() > 120) {
      H_j2_manytracks_et->Fill(jet2.et());
      H_j2_manytracks_eta->Fill(jet2.eta());
    }
  




























































    // For each jet: get the plane of the jet, project the tracks
    // inside of the jet in that plane, and calculate the thrust
    // axis on that plane.
    
    //########
    // Jet 1 #
    //########

//     // First jet: get the plane.
//     TVector3 theJetVector3 = jet1.Get4Vector().Vect();
//     Plane p1(theJetVector3);

//     // Project the tracks.
//     vector<TVector2> projTracks;
//     vector<TLorentzVector> pairedTracks = jetTracks1.second;
//     for(int i=0; i!=pairedTracks.size(); i++) {
//       TVector2 thisTrack = p1.projection2d(pairedTracks.at(i).Vect());
//       projTracks.push_back(thisTrack);
//     }
//     sort(projTracks.begin(), projTracks.end(), Mod2Comparator<TVector2>);
    
//     // Connect these tracks with the global pointer...
//     TracksPtr = &projTracks;
//     // ... and pray.

//     if(verbose) cout << "Projections ready" << endl;    
//     if(verbose) cout << "Num of tracks = " << projTracks.size() << endl;        
    
//     // At least two tracks in the jet, or minimization will fail badly...
//     if(projTracks.size() < 2) continue;
    
//     // Make the flow.
    
//     // The fitting procedure
//     TFitter minimizer(1);
//     // MAKE IT QUIET!!
//     {
//       double p1 = -1;
//       minimizer.ExecuteCommand("SET PRINTOUT",&p1,1);
//     }
//     // Tell the minimizer about the function to be minimized
//     // Look inside "MaxFlow.hh" to see the function.
//     minimizer.SetFCN(minuitFunction);
//     // Define the parameters
//     //   arg1 - parameter number
//     //   arg2 - parameter name
//     //   arg3 - first guess at parameter value
//     //   arg4 - estimated distance to minimum
//     //   arg5, arg6 - ignore for now
//     minimizer.SetParameter(0,"phi",projTracks.at(0).Phi(),1.5,0,0);
//     minimizer.ExecuteCommand("SIMPLEX",0,0);
//     minimizer.ExecuteCommand("MIGRAD",0,0);
    
//     if(verbose) cout << "Flow ready 1" << endl;
    
//     // Get the flow.
    
//     double bestPhi = minimizer.GetParameter(0);
//     TVector2 theFlowVector; theFlowVector.SetMagPhi(1.0,bestPhi);
//     TVector2 theAntiFlowVector (theFlowVector.Y(), -theFlowVector.X());

//     if(verbose) {
//       cout << "Flow and AntiFlow vector" << endl;
//       theFlowVector.Print();
//       theAntiFlowVector.Print();
//     }
    
//     double scalarFlow_Jet1 = 0.0;
//     double scalarAntiFlow_Jet1 = 0.0;
//     double denominator_Jet1 = 0.0;
//     for(int i=0; i!=projTracks.size(); i++) {
//       scalarFlow_Jet1 += fabs(projTracks.at(i)*theFlowVector);
//       scalarAntiFlow_Jet1 += fabs(projTracks.at(i)*theAntiFlowVector);
//       denominator_Jet1 += projTracks.at(i).Mod();
//     }
    
//     scalarFlow_Jet1 = (scalarFlow_Jet1/denominator_Jet1);
//     scalarAntiFlow_Jet1 = (scalarAntiFlow_Jet1/denominator_Jet1);

//     if(verbose) cout << "Flow ready 2" << endl;
    
//     if(verbose) {
//       cout << "Flow and AntiFlow" << endl;
//       cout << scalarFlow_Jet1 << endl;
//       cout << scalarAntiFlow_Jet1 << endl;
//     }
    
//     //########
//     // Jet 2 #
//     //########

//     // Second jet: get the plane.
//     theJetVector3 = jet2.Get4Vector().Vect();
//     Plane p2(theJetVector3);

//     // Project the tracks.
//     projTracks.clear();
//     pairedTracks = jetTracks2.second;
//     for(int i=0; i!=pairedTracks.size(); i++) {
//       TVector2 thisTrack = p2.projection2d(pairedTracks.at(i).Vect());
//       projTracks.push_back(thisTrack);
//     }
//     sort(projTracks.begin(), projTracks.end(), Mod2Comparator<TVector2>);
    
//     // Connect these tracks with the global pointer...
//     TracksPtr = &projTracks;
//     // ... and pray.

//     if(verbose) cout << "Projections ready" << endl;    
//     if(verbose) cout << "Num of tracks = " << projTracks.size() << endl;        
    
//     // At least two tracks in the jet, or minimization will fail badly...
//     if(projTracks.size() < 2) continue;
    
//     // Make the flow.
    
//     // Reusing the same minimizer.
//     minimizer.SetParameter(0,"phi",projTracks.at(0).Phi(),1.5,0,0);
//     minimizer.ExecuteCommand("SIMPLEX",0,0);
//     minimizer.ExecuteCommand("MIGRAD",0,0);
    
//     if(verbose) cout << "Flow ready 1" << endl;
    
//     // Get the flow.
    
//     bestPhi = minimizer.GetParameter(0);
//     theFlowVector.SetMagPhi(1.0,bestPhi);
//     theAntiFlowVector.Set(theFlowVector.Y(), -theFlowVector.X());

//     if(verbose) {
//       cout << "Flow and AntiFlow vector" << endl;
//       theFlowVector.Print();
//       theAntiFlowVector.Print();
//     }
    
//     double scalarFlow_Jet2 = 0.0;
//     double scalarAntiFlow_Jet2 = 0.0;
//     double denominator_Jet2 = 0.0;
//     for(int i=0; i!=projTracks.size(); i++) {
//       scalarFlow_Jet2 += fabs(projTracks.at(i)*theFlowVector);
//       scalarAntiFlow_Jet2 += fabs(projTracks.at(i)*theAntiFlowVector);
//       denominator_Jet2 += projTracks.at(i).Mod();
//     }
    
//     scalarFlow_Jet2 = (scalarFlow_Jet2/denominator_Jet2);
//     scalarAntiFlow_Jet2 = (scalarAntiFlow_Jet2/denominator_Jet2);

//     if(verbose) cout << "Flow ready 2" << endl;
     
//     if(verbose) {
//       cout << "Flow and AntiFlow" << endl;
//       cout << scalarFlow_Jet2 << endl;
//       cout << scalarAntiFlow_Jet2 << endl;
//     }
    
    // ############
    // Before cuts.
    // ############

//     // Fill the jet histos.
//     flowJet1->Fill(scalarFlow_Jet1,weight);
//     antiFlowJet1->Fill(scalarAntiFlow_Jet1,weight);
//     aFOverFJet1->Fill(scalarAntiFlow_Jet1/scalarFlow_Jet1,weight);
//     flowJet2->Fill(scalarFlow_Jet2,weight);
//     antiFlowJet2->Fill(scalarAntiFlow_Jet2,weight);
//     aFOverFJet2->Fill(scalarAntiFlow_Jet2/scalarFlow_Jet2,weight);
        
//     // Fill the corresponding basic histograms.
//     if(verbose) cout << "Filling 1D histos." << endl;
//     FillJetHistos(histosJet1, jetTracks1);
//     FillJetHistos(histosJet2, jetTracks2);

//     if(verbose) cout << "Filling 2D histos." << endl;
//     Fill2DHistos(histos2DJet1, jetTracks1);
//     Fill2DHistos(histos2DJet2, jetTracks2);
    
    // ###########
    // After cuts.
    // ###########
    
//     if(jet1.mass() < j1MinMass) continue;
//     if(jet1.mass() > j1MaxMass) continue;
//     if(jet2.mass() < j2MinMass) continue;
//     if(jet2.mass() < j1MaxMass) continue;
//     if(jet1.et()   < j1EtCut) continue;
//     if(jet2.et()   < j2EtCut) continue;

//     // Fill the jet histos.
//     flowJet1_AC->Fill(scalarFlow_Jet1,weight);
//     antiFlowJet1_AC->Fill(scalarAntiFlow_Jet1,weight);
//     aFOverFJet1_AC->Fill(scalarAntiFlow_Jet1/scalarFlow_Jet1,weight);
//     flowJet2_AC->Fill(scalarFlow_Jet2,weight);
//     antiFlowJet2_AC->Fill(scalarAntiFlow_Jet2,weight);
//     aFOverFJet2_AC->Fill(scalarAntiFlow_Jet2/scalarFlow_Jet2,weight);
        
//     // Fill the corresponding basic histograms.
//     if(verbose) cout << "Filling 1D histos." << endl;
//     FillJetHistos(histosJet1_AC, jetTracks1);
//     FillJetHistos(histosJet2_AC, jetTracks2);

//     if(verbose) cout << "Filling 2D histos." << endl;
//     Fill2DHistos(histos2DJet1_AC, jetTracks1);
//     Fill2DHistos(histos2DJet2_AC, jetTracks2);
    
  } //Closes the loop itself.
  
  // Open the TFile.
  TFile* theFile = new TFile(outfile.c_str(), "RECREATE");

  // Write the histograms to file.
//   flowJet1->Write();
//   antiFlowJet1->Write();
//   aFOverFJet1->Write();
//   flowJet2->Write();
//   antiFlowJet2->Write();
//   aFOverFJet2_AC->Write();
//   flowJet1_AC->Write();
//   antiFlowJet1_AC->Write();
//   aFOverFJet1_AC->Write();
//   flowJet2_AC->Write();
//   antiFlowJet2_AC->Write();
//   aFOverFJet2_AC->Write();
//   WriteHistos(histosJet1,theFile,"jet1_1D");
//   WriteHistos(histosJet2,theFile,"jet2_1D");
//   WriteHistos(histos2DJet1,theFile,"jet1_2D");
//   WriteHistos(histos2DJet2,theFile,"jet2_2D");
//   WriteHistos(histosJet1_AC,theFile,"jet1_1D_AC");
//   WriteHistos(histosJet2_AC,theFile,"jet2_1D_AC");
//   WriteHistos(histos2DJet1_AC,theFile,"jet1_2D_AC");
//   WriteHistos(histos2DJet2_AC,theFile,"jet2_2D_AC");

  H_jm1Xjm2->Write();
  H_j1_manytracks_et->Write();
  H_j1_manytracks_eta->Write();
  H_j2_manytracks_et->Write();
  H_j2_manytracks_eta->Write();

  theFile->Close();
  delete theFile;

} //Closes RSZZAnalysis::Loop

vector<TH1D*> RSZZAnalysis::CreateJetHistos(string dirname) {
 
  vector<TH1D*> theHistos;
  char name[256];
 
  sprintf(name, "%s_et", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 5000, 0.0, 5000.0));
  sprintf(name, "%s_eta", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 100, -5.0, 5.0));
  sprintf(name, "%s_mass", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 1000, 0.0, 1000.0));
  sprintf(name, "%s_tracks", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 1000, -0.5, 999.5));
  // Total: 4 histos.

  return theHistos;
}

vector<TH1D*> RSZZAnalysis::CreateBasicHistos(string dirname) {
 
  vector<TH1D*> theHistos;
  char name[256];
 
  sprintf(name, "%s_et", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 5000, 0.0, 5000.0));
  sprintf(name, "%s_eta", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 100, -5.0, 5.0));
  sprintf(name, "%s_mass", dirname.c_str());
  theHistos.push_back(new TH1D(name, name, 1000, 0.0, 1000.0));
  // Total: 3 histos.

  return theHistos;
}

vector<TH2D*> RSZZAnalysis::Create2DHistos(string dirname) {
  
  vector<TH2D*> theHistos;
  char name[256];
  
  sprintf(name, "%s_et_mass", dirname.c_str());
  theHistos.push_back(new TH2D(name, name, 5000, 0., 5000., 1000, 0., 1000.));
  sprintf(name, "%s_mass_numtracks", dirname.c_str());
  theHistos.push_back(new TH2D(name, name, 1000, 0., 1000., 200, -0.5, 199.5));
  sprintf(name, "%s_numtracks_et", dirname.c_str());
  theHistos.push_back(new TH2D(name, name, 200, -0.5, 199.5, 5000, 0., 5000.));
  // Total: 3 histos.

  return theHistos;
}

void RSZZAnalysis::FillJetHistos(vector<TH1D*>& theHistos, JetWithAssociation& theJet) {
  int ih = 0;
  int numTracks = 0;

  // Et, eta, mass, tracks.
  numTracks = theJet.second.size();
  theHistos[ih++]->Fill(theJet.first.et(),weight);
  theHistos[ih++]->Fill(theJet.first.eta(),weight);
  theHistos[ih++]->Fill(theJet.first.mass(),weight);
  theHistos[ih++]->Fill(numTracks,weight);
}

void RSZZAnalysis::FillBasicHistos(vector<TH1D*>& theHistos, TLorentzVector* theCandidate) {
  int ih = 0;
  
  // pt, eta, mass.
  theHistos[ih++]->Fill(theCandidate->Pt(),weight);
  theHistos[ih++]->Fill(theCandidate->Eta(),weight);
  theHistos[ih++]->Fill(theCandidate->M(),weight);
}

void RSZZAnalysis::Fill2DHistos(vector<TH2D*>& theHistos, JetWithAssociation& theJet) {
  int ih = 0;
  
  Jet j = theJet.first;
  
  double et = j.et();
  double mass = j.mass();
  int numtracks  = theJet.second.size();
  
  theHistos[ih++]->Fill(et, mass,weight);
  theHistos[ih++]->Fill(mass, numtracks,weight);
  theHistos[ih++]->Fill(numtracks, mass,weight);
}

int RSZZAnalysis::GetHardestPV(double threshold) {
  int iPV = -1;
  double maxpt = -1;
  for(int i = 0; i < nPV; i++) {
    if(SumPtPV[i] > maxpt && SumPtPV[i] > threshold){
      maxpt = SumPtPV[i];
      iPV = i;
    }
  }
  return iPV;
}

int RSZZAnalysis::NumTracks(TLorentzVector v, double r) {
  int numTracks = 0;
  for(int i=0; i!=nTrack; i++)
    if(DeltaR(double(etaTrack[i]), double(phiTrack[i]), v.Eta(), v.Phi()) < r) {
      double ptTrack = sqrt(pxTrack[i]*pxTrack[i] +
			    pyTrack[i]*pyTrack[i]);
      if(ptTrack > trackPtCut)
	numTracks++;
    }
  return numTracks;
}

int RSZZAnalysis::NumCharged(TLorentzVector v, double r) {
  int numCharged = 0;
  int theSize = pCharged->size();
  for(int i=0; i!=theSize; i++)
    if(DeltaR((*pCharged)[i].Eta(), (*pCharged)[i].Phi(), v.Eta(), v.Phi()) < r)
      numCharged++;
  return numCharged;
}

vector<Jet> RSZZAnalysis::FilterJetsByEt(vector<Jet> theJets, double jetEtCut) {
  
  vector<Jet> filteredJets;
  
  for(int i=0; i!=theJets.size(); ++i)
    if(theJets.at(i).et() > jetEtCut)
      filteredJets.push_back(theJets.at(i));

  return filteredJets;
  
}

void RSZZAnalysis::AssignParameters(map<string, string> sdata) {
  process = string (sdata["process"]);
}

void RSZZAnalysis::AssignParameters(map<string, double> ndata) {
  // Assign parameters ONLY if they are in the map.
  // This way, we keep the ones defined in InitParameters();
  
  map<string, double>::const_iterator theEnd = ndata.end();
  
  if(ndata.find("barrellimit") != theEnd)
    barrellimit  = ndata["barrellimit"];
  if(ndata.find("endcaplimit") != theEnd)
    endcaplimit  = ndata["endcaplimit"];
  if(ndata.find("minNumJets") != theEnd)
    minNumJets   = int(ndata["minNumJets"]);
  if(ndata.find("jetEtCut") != theEnd)
    jetEtCut     = ndata["jetEtCut"];
  if(ndata.find("j1EtCut") != theEnd)
    j1EtCut      = ndata["j1EtCut"];
  if(ndata.find("j2EtCut") != theEnd)
    j2EtCut      = ndata["j2EtCut"];
  if(ndata.find("j1MinMass") != theEnd)
    j1MinMass    = ndata["j1MinMass"];
  if(ndata.find("j2MinMass") != theEnd)
    j2MinMass    = ndata["j2MinMass"];
  if(ndata.find("j1MaxMass") != theEnd)
    j1MaxMass    = ndata["j1MaxMass"];
  if(ndata.find("j2MaxMass") != theEnd)
    j2MaxMass    = ndata["j2MaxMass"];
  if(ndata.find("trackPtCut") != theEnd)
    trackPtCut   = ndata["trackPtCut"];
  if(ndata.find("jetTracksCut") != theEnd)
    jetTracksCut = int(ndata["jetTracksCut"]);
  if(ndata.find("jetRadius") != theEnd)
    jetRadius    = ndata["jetRadius"];
  if(ndata.find("verbose") != theEnd)
    verbose      = ndata["verbose"];
}

void RSZZAnalysis::InitParameters() {
  process      = "RSZZ";
  barrellimit  = 1.3;
  endcaplimit  = 3.0;
  minNumJets   = 2;
  jetEtCut     = 30.0;
  j1EtCut      = 30.0;
  j2EtCut      = 30.0;
  j1MinMass    = 0.;
  j2MinMass    = 0.;  
  j1MaxMass    = 9999.;
  j2MaxMass    = 9999.;  
  trackPtCut   = 0.8;
  jetTracksCut = 0;
  jetRadius    = 0.7;
  verbose      = false;
}
