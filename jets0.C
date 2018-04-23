#define jets0_cxx
#include "jets0.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void jets0::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
    
    TFile* output;
    output = new TFile("../output.root", "UPDATE");
    //output = new TFile("../output.root", "RECREATE");
    
    int PTxaxis=1200.0;
    int PTyaxis=130.0;
    int PTbins=50;
    int Etaxaxis=10.0;
    int Etayaxis=130.0;
    int Etabins=50;
    int Rapidityxaxis=4.0;
    int Rapidityyaxis=130.0;
    int Rapiditybins=50;
    int Massxaxis=2000.0;
    int Massyaxis=130.0;
    int Massbins=50;
    
    TCanvas *c1 = new TCanvas("canvas_PTn2n2","Canvas for PTn2n2 Histogram",750,500);
    TCanvas *c2 = new TCanvas("canvas_Etan2n2","Canvas for Etan2n2 Histogram",750,500);
    TCanvas *c3 = new TCanvas("canvas_Rapidityn2n2","Canvas for Rapidityn2n2 Histogram",750,500);
    TCanvas *c4 = new TCanvas("canvas_Massn2n2","Canvas for Massn2n2 Histogram",750,500);
    
    TH1F *hist_PT_n2n2 = new TH1F("hist_PT_n2n2", "PT of n2n2", PTbins, 0.0, PTxaxis);
    TH1F *hist_Eta_n2n2 = new TH1F("hist_Eta_n2n2", "Eta of n2n2", Etabins, -Etaxaxis, Etaxaxis);
    TH1F *hist_Rapidity_n2n2 = new TH1F("hist_Rapidity_n2n2", "Rapidity of n2n2", Rapiditybins, -Rapidityxaxis, Rapidityxaxis);
    TH1F *hist_Mass_n2n2 = new TH1F("hist_Mass_n2n2", "Mass of n2n2", Massbins, 400.0, Massxaxis);
    
    //hist_PT_n2n2->GetYaxis()->SetRangeUser(0.0,PTyaxis);
    //hist_Eta_n2n2->GetYaxis()->SetRangeUser(0.0,Etayaxis);
    //hist_Rapidity_n2n2->GetYaxis()->SetRangeUser(0.0,Rapidityyaxis);
    //hist_Mass_n2n2->GetYaxis()->SetRangeUser(0.0,Massyaxis);
    
    int dummycount=0;
    //cout << nentries << endl;
    TLorentzVector n2A;
    TLorentzVector n2B;
    TLorentzVector n2;
    bool n2check=true;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        for(int i=0; i<Event_Nparticles[0]; i++)
        {
            if(Particle_PID[i]==1000023 && Particle_Status[i]==2 && n2check==true)
            {
                n2check=false;
                n2A.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
            }
            if(Particle_PID[i]==1000023 && Particle_Status[i]==2 && n2check==false)
            {
                n2B.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
            }
        }
        n2check=true;
        n2=n2A+n2B;
        hist_PT_n2n2->Fill(n2.Pt());
        hist_Eta_n2n2->Fill(n2.Eta());
        hist_Rapidity_n2n2->Fill(n2.Rapidity());
        hist_Mass_n2n2->Fill(n2.M());
    }
    //cout << dummycount << endl;
    
    hist_PT_n2n2->Scale(1.0/hist_PT_n2n2->Integral("width"));
    hist_PT_n2n2->SetLineColor(kRed);
    hist_Eta_n2n2->Scale(1.0/hist_Eta_n2n2->Integral("width"));
    hist_Eta_n2n2->SetLineColor(kRed);
    hist_Rapidity_n2n2->Scale(1.0/hist_Rapidity_n2n2->Integral("width"));
    hist_Rapidity_n2n2->SetLineColor(kRed);
    hist_Mass_n2n2->Scale(1.0/hist_Mass_n2n2->Integral("width"));
    hist_Mass_n2n2->SetLineColor(kRed);
    output->Write();
    output->Close();
}
