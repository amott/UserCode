//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun  4 12:20:46 2010 by ROOT version 5.22/00d
// from TTree t3/Reconst ntuple
// found on file: /data1/tomei/GenVecbos/z1j_0ptz100.root
//////////////////////////////////////////////////////////

#ifndef GenVecbos_h
#define GenVecbos_h


// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

// FASTJET includes
#include "FASTJET/include/fastjet/PseudoJet.hh"
#include "FASTJET/include/fastjet/ClusterSequence.hh"
#include "FASTJET/include/fastjet/ClusterSequenceActiveArea.hh"
#include "FASTJET/include/fastjet/SISConePlugin.hh"

#include <vector>
#include <Jet.hh>

using namespace std;

class GenVecbos {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           npar;
   Float_t         px[10000];   //[npar]
   Float_t         py[10000];   //[npar]
   Float_t         pz[10000];   //[npar]
   Float_t         E[10000];   //[npar]
   Float_t         m[10000];   //[npar]
   Int_t           id[10000];   //[npar]
   Int_t           index[10000];   //[npar]
   Int_t           status[10000];   //[npar]
   Int_t           jmo1[10000];   //[npar]
   Int_t           jmo2[10000];   //[npar]
   Int_t           jda1[10000];   //[npar]
   Int_t           jda2[10000];   //[npar]

   // List of branches
   TBranch        *b_npar;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_m;   //!
   TBranch        *b_id;   //!
   TBranch        *b_index;   //!
   TBranch        *b_status;   //!
   TBranch        *b_jmo1;   //!
   TBranch        *b_jmo2;   //!
   TBranch        *b_jda1;   //!
   TBranch        *b_jda2;   //!

   GenVecbos(TTree *tree=0);
   virtual ~GenVecbos();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

 protected:
   vector<Jet> SISCone(vector<TLorentzVector> InputCollection, double Rparam, double thePtMin);
};

#endif

#ifdef GenVecbos_cxx
GenVecbos::GenVecbos(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data1/tomei/GenVecbos/z1j_0ptz100.root");
      if (!f) {
         f = new TFile("/data1/tomei/GenVecbos/z1j_0ptz100.root");
      }
      tree = (TTree*)gDirectory->Get("t3");

   }
   Init(tree);
}

GenVecbos::~GenVecbos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GenVecbos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GenVecbos::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GenVecbos::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npar", &npar, &b_npar);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("index", index, &b_index);
   fChain->SetBranchAddress("status", status, &b_status);
   fChain->SetBranchAddress("jmo1", jmo1, &b_jmo1);
   fChain->SetBranchAddress("jmo2", jmo2, &b_jmo2);
   fChain->SetBranchAddress("jda1", jda1, &b_jda1);
   fChain->SetBranchAddress("jda2", jda2, &b_jda2);
   Notify();
}

Bool_t GenVecbos::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GenVecbos::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GenVecbos::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GenVecbos_cxx
