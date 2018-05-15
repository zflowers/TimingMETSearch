/* ---------------------------------------------------------------------------- *
 Author: Zach Flowers
 Description:  Analyze the tree output after collecting relevant branches
    from pythia and delphes
 Date: April 23rd 2018
 * ---------------------------------------------------------------------------- */

#define Analyze_n2n2j_cxx
#include "Analyze_n2n2j.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analyze_n2n2j::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

    //create canvas for later
    
    TCanvas *c1 = new TCanvas("c1","Canvas for PT",750,500);
    
    TCanvas *c2 = new TCanvas("c2","Canvas for PT Difference",750,500);
    
    //create histograms to look at delphes smearing of jet PT
    
    TH2F *hist_pythia_Pt_vs_delphes_Pt_n2n2j = new TH2F("hist_pythia_Pt_vs_delphes_Pt_n2n2j","Pythia Jets Pt vs Delphes Jets Pt", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
    
    TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,0.0,700.0);
    
    //dummy variables to fill histograms
    
    double pythia_JET_PT=0.0;
    double delphes_JET_PT=0.0;
    double PT_Difference=0.0;
    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //event loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       for(int i=0; i<Jet_MC_Size; i++) //loop through every MC jet in the tree
       {
           pythia_JET_PT+=Jet_MC_PT->at(i); //scalar sum the PT of the MC jets in an event
       }
       for(int i=0; i<HT_Size; i++) //Fill dummy variable with Delphes Jet PT (Scalar_HT)
       {
           delphes_JET_PT=HT;
       }
       PT_Difference=(pythia_JET_PT-delphes_JET_PT); //compute absolute value of the difference of the PT from the jets in pythia and delphes
       //fill histograms
       hist_pythia_Pt_vs_delphes_Pt_n2n2j->Fill(pythia_JET_PT,delphes_JET_PT);
       hist_PT_Difference->Fill(abs(PT_Difference));
       //reset variables for the next event
       pythia_JET_PT=0.0;
       delphes_JET_PT=0.0;
   }
    //Draw and save histograms
    c1->cd();
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->SetStats(false);
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->GetXaxis()->SetTitle("Pythia (MC) Scalar Sum of PT of Jet(s) [GeV]");
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->GetYaxis()->SetTitle("Delphes (RECO) Scalar HT [GeV]");
    hist_pythia_Pt_vs_delphes_Pt_n2n2j->Draw("COLZ");
    c1->SaveAs("JetPT.png");
    c2->cd();
    //hist_PT_Difference->SetStats(false);
    hist_PT_Difference->GetXaxis()->SetTitle("Difference in Summed PT of Jets from Pythia and Delphes [GeV]");
    hist_PT_Difference->GetYaxis()->SetTitle("Number of Events");
    hist_PT_Difference->SetLineColor(kBlack);
    hist_PT_Difference->Draw();
    c2->SaveAs("PT_Difference.png");
}
