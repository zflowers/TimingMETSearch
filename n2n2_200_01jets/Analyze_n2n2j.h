//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 25 10:44:12 2018 by ROOT version 6.10/08
// from TTree n2n2j/Output
// found on file: n2n2j_output.root
//////////////////////////////////////////////////////////

#ifndef Analyze_n2n2j_h
#define Analyze_n2n2j_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TVector3.h"
#include "vector"
#include "TLorentzVector.h"

class Analyze_n2n2j {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   TVector3        *MET;
   Float_t         HT;
   Int_t           HT_Size;
   vector<double>  *Electron_Delphes_PT;
   vector<double>  *Electron_Delphes_Eta;
   vector<double>  *Electron_Delphes_Phi;
   vector<double>  *Electron_Delphes_E;
   vector<double>  *Electron_Delphes_Charge;
   Int_t           Electron_Delphes_Size;
   vector<double>  *Muon_Delphes_PT;
   vector<double>  *Muon_Delphes_Eta;
   vector<double>  *Muon_Delphes_Phi;
   vector<double>  *Muon_Delphes_E;
   vector<double>  *Muon_Delphes_Charge;
   Int_t           Muon_Delphes_Size;
   vector<double>  *Jet_Delphes_PT;
   vector<double>  *Jet_Delphes_Eta;
   vector<double>  *Jet_Delphes_Phi;
   vector<double>  *Jet_Delphes_E;
   vector<double>  *Jet_Delphes_Mass;
   Int_t           Jet_Delphes_Size;
   TLorentzVector  *n2A_MC;
   TLorentzVector  *n2B_MC;
   TLorentzVector  *n1A_MC;
   TLorentzVector  *n1B_MC;
   vector<double>  *Electron_MC_PT;
   vector<double>  *Electron_MC_Eta;
   vector<double>  *Electron_MC_Phi;
   vector<double>  *Electron_MC_E;
   Int_t           Electron_MC_Size;
   vector<double>  *Muon_MC_PT;
   vector<double>  *Muon_MC_Eta;
   vector<double>  *Muon_MC_Phi;
   vector<double>  *Muon_MC_E;
   Int_t           Muon_MC_Size;
   vector<double>  *Jet_MC_PT;
   vector<double>  *Jet_MC_Eta;
   vector<double>  *Jet_MC_Phi;
   vector<double>  *Jet_MC_E;
   vector<double>  *Jet_MC_Mass;
   Int_t           Jet_MC_Size;

   // List of branches
   TBranch        *b_MET;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_HT_Size;   //!
   TBranch        *b_Electron_Delphes_PT;   //!
   TBranch        *b_Electron_Delphes_Eta;   //!
   TBranch        *b_Electron_Delphes_Phi;   //!
   TBranch        *b_Electron_Delphes_E;   //!
   TBranch        *b_Electron_Delphes_Charge;   //!
   TBranch        *b_Electron_Delphes_Size;   //!
   TBranch        *b_Muon_Delphes_PT;   //!
   TBranch        *b_Muon_Delphes_Eta;   //!
   TBranch        *b_Muon_Delphes_Phi;   //!
   TBranch        *b_Muon_Delphes_E;   //!
   TBranch        *b_Muon_Delphes_Charge;   //!
   TBranch        *b_Muon_Delphes_Size;   //!
   TBranch        *b_Jet_Delphes_PT;   //!
   TBranch        *b_Jet_Delphes_Eta;   //!
   TBranch        *b_Jet_Delphes_Phi;   //!
   TBranch        *b_Jet_Delphes_E;   //!
   TBranch        *b_Jet_Delphes_Mass;   //!
   TBranch        *b_Jet_Delphes_Size;   //!
   TBranch        *b_n2A_MC;   //!
   TBranch        *b_n2B_MC;   //!
   TBranch        *b_n1A_MC;   //!
   TBranch        *b_n1B_MC;   //!
   TBranch        *b_Electron_MC_PT;   //!
   TBranch        *b_Electron_MC_Eta;   //!
   TBranch        *b_Electron_MC_Phi;   //!
   TBranch        *b_Electron_MC_E;   //!
   TBranch        *b_Electron_MC_Size;   //!
   TBranch        *b_Muon_MC_PT;   //!
   TBranch        *b_Muon_MC_Eta;   //!
   TBranch        *b_Muon_MC_Phi;   //!
   TBranch        *b_Muon_MC_E;   //!
   TBranch        *b_Muon_MC_Size;   //!
   TBranch        *b_Jet_MC_PT;   //!
   TBranch        *b_Jet_MC_Eta;   //!
   TBranch        *b_Jet_MC_Phi;   //!
   TBranch        *b_Jet_MC_E;   //!
   TBranch        *b_Jet_MC_Mass;   //!
   TBranch        *b_Jet_MC_Size;   //!

   Analyze_n2n2j(TTree *tree=0);
   virtual ~Analyze_n2n2j();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyze_n2n2j_cxx
Analyze_n2n2j::Analyze_n2n2j(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("n2n2j_output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("n2n2j_output.root");
      }
      f->GetObject("n2n2j",tree);

   }
   Init(tree);
}

Analyze_n2n2j::~Analyze_n2n2j()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyze_n2n2j::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyze_n2n2j::LoadTree(Long64_t entry)
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

void Analyze_n2n2j::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   MET = 0;
   Electron_Delphes_PT = 0;
   Electron_Delphes_Eta = 0;
   Electron_Delphes_Phi = 0;
   Electron_Delphes_E = 0;
   Electron_Delphes_Charge = 0;
   Muon_Delphes_PT = 0;
   Muon_Delphes_Eta = 0;
   Muon_Delphes_Phi = 0;
   Muon_Delphes_E = 0;
   Muon_Delphes_Charge = 0;
   Jet_Delphes_PT = 0;
   Jet_Delphes_Eta = 0;
   Jet_Delphes_Phi = 0;
   Jet_Delphes_E = 0;
   Jet_Delphes_Mass = 0;
   n2A_MC = 0;
   n2B_MC = 0;
   n1A_MC = 0;
   n1B_MC = 0;
   Electron_MC_PT = 0;
   Electron_MC_Eta = 0;
   Electron_MC_Phi = 0;
   Electron_MC_E = 0;
   Muon_MC_PT = 0;
   Muon_MC_Eta = 0;
   Muon_MC_Phi = 0;
   Muon_MC_E = 0;
   Jet_MC_PT = 0;
   Jet_MC_Eta = 0;
   Jet_MC_Phi = 0;
   Jet_MC_E = 0;
   Jet_MC_Mass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("HT_Size", &HT_Size, &b_HT_Size);
   fChain->SetBranchAddress("Electron_Delphes_PT", &Electron_Delphes_PT, &b_Electron_Delphes_PT);
   fChain->SetBranchAddress("Electron_Delphes_Eta", &Electron_Delphes_Eta, &b_Electron_Delphes_Eta);
   fChain->SetBranchAddress("Electron_Delphes_Phi", &Electron_Delphes_Phi, &b_Electron_Delphes_Phi);
   fChain->SetBranchAddress("Electron_Delphes_E", &Electron_Delphes_E, &b_Electron_Delphes_E);
   fChain->SetBranchAddress("Electron_Delphes_Charge", &Electron_Delphes_Charge, &b_Electron_Delphes_Charge);
   fChain->SetBranchAddress("Electron_Delphes_Size", &Electron_Delphes_Size, &b_Electron_Delphes_Size);
   fChain->SetBranchAddress("Muon_Delphes_PT", &Muon_Delphes_PT, &b_Muon_Delphes_PT);
   fChain->SetBranchAddress("Muon_Delphes_Eta", &Muon_Delphes_Eta, &b_Muon_Delphes_Eta);
   fChain->SetBranchAddress("Muon_Delphes_Phi", &Muon_Delphes_Phi, &b_Muon_Delphes_Phi);
   fChain->SetBranchAddress("Muon_Delphes_E", &Muon_Delphes_E, &b_Muon_Delphes_E);
   fChain->SetBranchAddress("Muon_Delphes_Charge", &Muon_Delphes_Charge, &b_Muon_Delphes_Charge);
   fChain->SetBranchAddress("Muon_Delphes_Size", &Muon_Delphes_Size, &b_Muon_Delphes_Size);
   fChain->SetBranchAddress("Jet_Delphes_PT", &Jet_Delphes_PT, &b_Jet_Delphes_PT);
   fChain->SetBranchAddress("Jet_Delphes_Eta", &Jet_Delphes_Eta, &b_Jet_Delphes_Eta);
   fChain->SetBranchAddress("Jet_Delphes_Phi", &Jet_Delphes_Phi, &b_Jet_Delphes_Phi);
   fChain->SetBranchAddress("Jet_Delphes_E", &Jet_Delphes_E, &b_Jet_Delphes_E);
   fChain->SetBranchAddress("Jet_Delphes_Mass", &Jet_Delphes_Mass, &b_Jet_Delphes_Mass);
   fChain->SetBranchAddress("Jet_Delphes_Size", &Jet_Delphes_Size, &b_Jet_Delphes_Size);
   fChain->SetBranchAddress("n2A_MC", &n2A_MC, &b_n2A_MC);
   fChain->SetBranchAddress("n2B_MC", &n2B_MC, &b_n2B_MC);
   fChain->SetBranchAddress("n1A_MC", &n1A_MC, &b_n1A_MC);
   fChain->SetBranchAddress("n1B_MC", &n1B_MC, &b_n1B_MC);
   fChain->SetBranchAddress("Electron_MC_PT", &Electron_MC_PT, &b_Electron_MC_PT);
   fChain->SetBranchAddress("Electron_MC_Eta", &Electron_MC_Eta, &b_Electron_MC_Eta);
   fChain->SetBranchAddress("Electron_MC_Phi", &Electron_MC_Phi, &b_Electron_MC_Phi);
   fChain->SetBranchAddress("Electron_MC_E", &Electron_MC_E, &b_Electron_MC_E);
   fChain->SetBranchAddress("Electron_MC_Size", &Electron_MC_Size, &b_Electron_MC_Size);
   fChain->SetBranchAddress("Muon_MC_PT", &Muon_MC_PT, &b_Muon_MC_PT);
   fChain->SetBranchAddress("Muon_MC_Eta", &Muon_MC_Eta, &b_Muon_MC_Eta);
   fChain->SetBranchAddress("Muon_MC_Phi", &Muon_MC_Phi, &b_Muon_MC_Phi);
   fChain->SetBranchAddress("Muon_MC_E", &Muon_MC_E, &b_Muon_MC_E);
   fChain->SetBranchAddress("Muon_MC_Size", &Muon_MC_Size, &b_Muon_MC_Size);
   fChain->SetBranchAddress("Jet_MC_PT", &Jet_MC_PT, &b_Jet_MC_PT);
   fChain->SetBranchAddress("Jet_MC_Eta", &Jet_MC_Eta, &b_Jet_MC_Eta);
   fChain->SetBranchAddress("Jet_MC_Phi", &Jet_MC_Phi, &b_Jet_MC_Phi);
   fChain->SetBranchAddress("Jet_MC_E", &Jet_MC_E, &b_Jet_MC_E);
   fChain->SetBranchAddress("Jet_MC_Mass", &Jet_MC_Mass, &b_Jet_MC_Mass);
   fChain->SetBranchAddress("Jet_MC_Size", &Jet_MC_Size, &b_Jet_MC_Size);
   Notify();
}

Bool_t Analyze_n2n2j::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyze_n2n2j::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyze_n2n2j::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyze_n2n2j_cxx
