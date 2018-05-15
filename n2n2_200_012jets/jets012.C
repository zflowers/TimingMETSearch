/* ----------------------------------------------------------------------------
 Author: Zach Flowers
 Description:  Adds the 012 jet histograms to the output root file
    Typically ran before 01 and 0 jet contribution files
 Date: April 23rd 2018
 * ---------------------------------------------------------------------------- */

#define jets012_cxx
#include "jets012.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void jets012::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
    
    TFile* output;
    //output = new TFile("../output.root", "UPDATE");
    output = new TFile("../output.root", "RECREATE"); //create the ROOT File for storing histograms
    
    //initial histogram bins and ranges
    
    double PTxaxis=1200.0;
    double PTyaxis=130.0;
    int PTbins=50;
    double Etaxaxis=10.0;
    double Etayaxis=130.0;
    int Etabins=50;
    double Rapidityxaxis=4.0;
    double Rapidityyaxis=130.0;
    int Rapiditybins=50;
    double Massxaxis=2000.0;
    double Massyaxis=130.0;
    int Massbins=50;
    
    //canvas for histograms
    
    TCanvas *c1 = new TCanvas("canvas_PTn2n2jj","Canvas for PTn2n2jj Histogram",750,500);
    TCanvas *c2 = new TCanvas("canvas_Etan2n2jj","Canvas for Etan2n2jj Histogram",750,500);
    TCanvas *c3 = new TCanvas("canvas_Rapidityn2n2jj","Canvas for Rapidityn2n2jj Histogram",750,500);
    TCanvas *c4 = new TCanvas("canvas_Massn2n2jj","Canvas for Massn2n2jj Histogram",750,500);
    
    //histograms
    
    TH1F *hist_PT_n2n2jj = new TH1F("hist_PT_n2n2jj", "PT of n2n2jj", PTbins, 0.0, PTxaxis);
    TH1F *hist_Eta_n2n2jj = new TH1F("hist_Eta_n2n2jj", "Eta of n2n2jj", Etabins, -Etaxaxis, Etaxaxis);
    TH1F *hist_Rapidity_n2n2jj = new TH1F("hist_Rapidity_n2n2jj", "Rapidity of n2n2jj", Rapiditybins, -Rapidityxaxis, Rapidityxaxis);
    TH1F *hist_Mass_n2n2jj = new TH1F("hist_Mass_n2n2jj", "Mass of n2n2jj", Massbins, 400.0, Massxaxis);

    //Uncomment below to set the range on the y axis to some specified value
    
    //hist_PT_n2n2jj->GetYaxis()->SetRangeUser(0.0,PTyaxis);
    //hist_Eta_n2n2jj->GetYaxis()->SetRangeUser(0.0,Etayaxis);
    //hist_Rapidity_n2n2jj->GetYaxis()->SetRangeUser(0.0,Rapidityyaxis);
    //hist_Mass_n2n2jj->GetYaxis()->SetRangeUser(0.0,Massyaxis);
    
    TLorentzVector n2A; //TLV for the first neutralino in the event
    TLorentzVector n2B; //TLV for the second neutralino in the event
    TLorentzVector n2; //Dummy TLV for summing the TLV of the neutralino TLVs
                        // Useful for getting vectorial sums of particles
    bool n2check=true; //boolean for checking if the first neutralino in an event has been identified
    
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //event loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       for(int i=0; i<Event_Nparticles[0]; i++) //particle loop
       {
           if(Particle_PID[i]==1000023 && Particle_Status[i]==2 && n2check==true) //if block for finding the first neutralino in an event and setting a TLV for it
           {
               n2check=false;
               n2A.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
           }
           if(Particle_PID[i]==1000023 && Particle_Status[i]==2 && n2check==false) //if block for finding the second neutralino in an event and setting a TLV for it
           {
               n2B.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
           }
       }
       n2check=true; //reset the boolean checker
       n2=n2A+n2B; //sum the neutralinos in the event
       //fill histograms
       hist_PT_n2n2jj->Fill(n2.Pt());
       hist_Eta_n2n2jj->Fill(n2.Eta());
       hist_Rapidity_n2n2jj->Fill(n2.Rapidity());
       hist_Mass_n2n2jj->Fill(n2.M());
    }
    //normalize the histograms and set color scheme
    hist_PT_n2n2jj->Scale(1.0/hist_PT_n2n2jj->Integral("width"));
    hist_PT_n2n2jj->SetLineColor(kBlack);
    hist_Eta_n2n2jj->Scale(1.0/hist_Eta_n2n2jj->Integral("width"));
    hist_Eta_n2n2jj->SetLineColor(kBlack);
    hist_Rapidity_n2n2jj->Scale(1.0/hist_Rapidity_n2n2jj->Integral("width"));
    hist_Rapidity_n2n2jj->SetLineColor(kBlack);
    hist_Mass_n2n2jj->Scale(1.0/hist_Mass_n2n2jj->Integral("width"));
    hist_Mass_n2n2jj->SetLineColor(kBlack);
    //Write and close TFile
    output->Write();
    output->Close();
}
