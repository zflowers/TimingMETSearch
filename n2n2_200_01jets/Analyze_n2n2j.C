#define Analyze_n2n2j_cxx
#include "Analyze_n2n2j.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analyze_n2n2j::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Analyze_n2n2j.C
//      root> Analyze_n2n2j t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

    TCanvas *c1 = new TCanvas("c1","Canvas for PT",750,500);
    c1->cd();
    
    TH2F *hist_pythia_Pt_vs_delphes_Pt_n2n2j = new TH2F("hist_pythia_Pt_vs_delphes_Pt_n2n2j","Pythia Jets Pt vs Delphes Jets Pt", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
    
    double pythia_JET_PT=0.0;
    double delphes_JET_PT=0.0;
    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       for(int i=0; i<Jet_MC_Size; i++)
       {
           pythia_JET_PT+=Jet_MC_PT->at(i);
       }
       for(int i=0; i<HT_Size; i++)
       {
           delphes_JET_PT=HT;
       }
       hist_pythia_Pt_vs_delphes_Pt_n2n2j->Fill(pythia_JET_PT,delphes_JET_PT);
       pythia_JET_PT=0.0;
       delphes_JET_PT=0.0;
   }
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->SetStats(false);
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->GetXaxis()->SetTitle("Pythia (MC) Scalar Sum of PT of Jet(s) [GeV]");
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->GetYaxis()->SetTitle("Delphes (RECO) Scalar HT [GeV]");
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->Draw("COLZ");
    c1->SaveAs("JetPT.png");
}
