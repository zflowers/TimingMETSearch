#define Analyze_n2n2j_cxx
#include "Analyze_n2n2j.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analyze_n2n2j::Loop()
{
    double MET_Parallel_Mag=0.0;
    double MET_Perpendicular_Mag=0.0;
    TVector3 MET_Perpendicular;
    TVector3 BeamSpot;
    BeamSpot.SetXYZ(0.0,0.0,1.0);
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

    TCanvas *c1 = new TCanvas("c1","Canvas for PT",750,500);
    TCanvas *c2 = new TCanvas("c2","Canvas for PT Difference",750,500);
    TCanvas *c3 = new TCanvas("c3","Canvas for MET Parallel",750,500);
    TCanvas *c4 = new TCanvas("c4","Canvas for MET Perpendicular",750,500);
    
    TH2F *hist_pythia_Pt_vs_delphes_Pt_n2n2j = new TH2F("hist_pythia_Pt_vs_delphes_Pt_n2n2j","Pythia Jets Pt vs Delphes Jets Pt", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
    
    TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,-700.0,200.0);
    //TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,0.0,700.0);
    
    TH1D *hist_MET_Parallel = new TH1D("hist_MET_Parallel","MET_Parallel",100,-200.0,200.0);
    TH1D *hist_MET_Perpendicular = new TH1D("hist_MET_Perpendicular","MET_Perpendicular",100,-200.0,200.0);
    
    double pythia_JET_PT=0.0;
    double delphes_JET_PT=0.0;
    double PT_Difference=0.0;
    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       //Break up MET into parallel and perpendicular components
       TVector3 PT_true;
       TVector3 n1A_vect=n1A_MC->Vect();
       TVector3 n1B_vect=n1B_MC->Vect();
       PT_true=n1A_vect+n1B_vect;
       TVector3 MET_RECO=*MET;
       PT_true.SetZ(0.0);
       MET_RECO.SetZ(0.0);

       MET_Perpendicular_Mag=PT_true.Cross(BeamSpot).Unit().Dot(MET_RECO);
       MET_Parallel_Mag=(MET_RECO.Dot(PT_true))/PT_true.Mag();
       
       hist_MET_Parallel->Fill(MET_Parallel_Mag-PT_true.Mag());
       hist_MET_Perpendicular->Fill(MET_Perpendicular_Mag);
       
       //reset Magnitudes
       MET_Parallel_Mag=0.0;
       MET_Perpendicular_Mag=0.0;
       
       
       for(int i=0; i<Jet_MC_Size; i++)
       {
           pythia_JET_PT+=Jet_MC_PT->at(i);
       }
       for(int i=0; i<HT_Size; i++)
       {
           delphes_JET_PT=HT;
       }
       PT_Difference=(pythia_JET_PT-delphes_JET_PT);
       hist_pythia_Pt_vs_delphes_Pt_n2n2j->Fill(pythia_JET_PT,delphes_JET_PT);
       hist_PT_Difference->Fill(PT_Difference);
       //hist_PT_Difference->Fill(abs(PT_Difference));
       pythia_JET_PT=0.0;
       delphes_JET_PT=0.0;
   }
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
    c2->SaveAs("PT_Pythia_Minus_Delphes.png");
    //c2->SaveAs("abs_PT_Difference.png");
    c3->cd();
    //hist_PT_Difference->SetStats(false);
    hist_MET_Parallel->GetXaxis()->SetTitle("MET Parallel [GeV]");
    hist_MET_Parallel->GetYaxis()->SetTitle("Number of Events");
    hist_MET_Parallel->SetLineColor(kBlack);
    hist_MET_Parallel->Draw();
    c3->SaveAs("MET_Parallel.png");
    c4->cd();
    //hist_PT_Difference->SetStats(false);
    hist_MET_Perpendicular->GetXaxis()->SetTitle("MET Perpendicular [GeV]");
    hist_MET_Perpendicular->GetYaxis()->SetTitle("Number of Events");
    hist_MET_Perpendicular->SetLineColor(kBlack);
    hist_MET_Perpendicular->Draw();
    c4->SaveAs("MET_Perpendicular.png");
}
