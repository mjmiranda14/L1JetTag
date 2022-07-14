//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  1 14:10:02 2022 by ROOT version 6.16/00
// from TTree objects/objects
// found on file: pfTuple.root
//////////////////////////////////////////////////////////

#ifndef L1JetTagPlotter_h
#define L1JetTagPlotter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.
#include <iostream>
#include <vector>
#include "TLorentzVector.h"
#include <utility> //for pair 

using namespace std;

class L1JetTagPlotter {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::vector<pair<TLorentzVector,int> > *emcalo;
   std::vector<pair<TLorentzVector,int> > *egcalo;
   std::vector<pair<TLorentzVector,int> > *calo;
   std::vector<pair<TLorentzVector,int> > *pf;
   std::vector<pair<TLorentzVector,int> > *pup;
   std::vector<pair<TLorentzVector,int> > *gen;
   std::vector<pair<TLorentzVector,int> > *l1jet;
   std::vector<pair<TLorentzVector,int> > *recojet;
   std::vector<float>   *pf_vx;
   std::vector<float>   *pf_vy;
   std::vector<float>   *pf_vz;
   std::vector<float>   *pup_vx;
   std::vector<float>   *pup_vy;
   std::vector<float>   *pup_vz;
   std::vector<float>   *vz;

   // List of branches
   TBranch        *b_emcalo;   //!
   TBranch        *b_egcalo;   //!
   TBranch        *b_calo;   //!
   TBranch        *b_pf;   //!
   TBranch        *b_pup;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_l1jet;   //!
   TBranch        *b_recojet;   //!
   TBranch        *b_pf_vx;   //!
   TBranch        *b_pf_vy;   //!
   TBranch        *b_pf_vz;   //!
   TBranch        *b_pup_vx;   //!
   TBranch        *b_pup_vy;   //!
   TBranch        *b_pup_vz;   //!
   TBranch        *b_vz;   //!

   TH1F*          h_pf_pt;
   TH1F*          h_pf_vz;


   L1JetTagPlotter(TTree *tree=0);
   virtual ~L1JetTagPlotter();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString outFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef L1JetTagPlotter_cxx
L1JetTagPlotter::L1JetTagPlotter(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pfTuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pfTuple.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("pfTuple.root:/ntuple0");
      dir->GetObject("objects",tree);

   }
   Init(tree);
}

L1JetTagPlotter::~L1JetTagPlotter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t L1JetTagPlotter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t L1JetTagPlotter::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void L1JetTagPlotter::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   emcalo = 0;
   egcalo = 0;
   calo = 0;
   pf = 0;
   pup = 0;
   gen = 0;
   l1jet = 0;
   recojet = 0;
   pf_vx = 0;
   pf_vy = 0;
   pf_vz = 0;
   pup_vx = 0;
   pup_vy = 0;
   pup_vz = 0;
   vz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("emcalo", &emcalo, &b_emcalo);
   fChain->SetBranchAddress("egcalo", &egcalo, &b_egcalo);
   fChain->SetBranchAddress("calo", &calo, &b_calo);
   fChain->SetBranchAddress("pf", &pf, &b_pf);
   fChain->SetBranchAddress("pup", &pup, &b_pup);
   fChain->SetBranchAddress("gen", &gen, &b_gen);
   fChain->SetBranchAddress("l1jet", &l1jet, &b_l1jet);
   fChain->SetBranchAddress("recojet", &recojet, &b_recojet);
   fChain->SetBranchAddress("pf_vx", &pf_vx, &b_pf_vx);
   fChain->SetBranchAddress("pf_vy", &pf_vy, &b_pf_vy);
   fChain->SetBranchAddress("pf_vz", &pf_vz, &b_pf_vz);
   fChain->SetBranchAddress("pup_vx", &pup_vx, &b_pup_vx);
   fChain->SetBranchAddress("pup_vy", &pup_vy, &b_pup_vy);
   fChain->SetBranchAddress("pup_vz", &pup_vz, &b_pup_vz);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   Notify();
}

Bool_t L1JetTagPlotter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void L1JetTagPlotter::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t L1JetTagPlotter::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef L1JetTagPlotter_cxx
