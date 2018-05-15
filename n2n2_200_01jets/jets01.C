/* ---------------------------------------------------------------------------- *
 Author: Zach Flowers
 Description:  Adds the 01 jet histograms to the output root file
 Run after 012 jets contribution file
 Contains new algorithm for tagging other particles in the event
 Date: April 23rd 2018
 * ---------------------------------------------------------------------------- */

#define jets01_cxx
#include "jets01.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>
#include <TTree.h>
#include <TFile.h>

void jets01::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
    
    TFile* output;
    output = new TFile("../output.root", "UPDATE"); //open the output root file to store histograms
    //output = new TFile("../output.root", "RECREATE");
    
    
    
    TLorentzVector zeroTLV; //make a TLorentzVector (TLV) and set all entries to zero
    zeroTLV.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
    //create TLVs for each particle in the event
    TLorentzVector n2A_MC; //TLV for the first neutralino in the event
    TLorentzVector n2B_MC; //TLV for the second neutralino in the event
    TLorentzVector n1A_MC;
    TLorentzVector n1B_MC;
    TLorentzVector lA_MC;
    TLorentzVector lB_MC;
    TLorentzVector lC_MC;
    TLorentzVector lD_MC;
    TLorentzVector n2; //Dummy TLV for summing the TLV of the neutralino TLVs
                        // Useful for getting vectorial sums of particles
    
    //set all of the branches to zero
    n2A_MC=zeroTLV;
    n2B_MC=zeroTLV;
    n1A_MC=zeroTLV;
    n1B_MC=zeroTLV;
    lA_MC=zeroTLV;
    lB_MC=zeroTLV;
    lC_MC=zeroTLV;
    lD_MC=zeroTLV;
    
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
    
    TCanvas *c1 = new TCanvas("canvas_PTn2n2j","Canvas for PTn2n2j Histogram",750,500);
    TCanvas *c2 = new TCanvas("canvas_Etan2n2j","Canvas for Etan2n2j Histogram",750,500);
    TCanvas *c3 = new TCanvas("canvas_Rapidityn2n2j","Canvas for Rapidityn2n2j Histogram",750,500);
    TCanvas *c4 = new TCanvas("canvas_Massn2n2j","Canvas for Massn2n2j Histogram",750,500);
    TCanvas *c5 = new TCanvas("canvas_PtvsEta","Canvas for Pt vs Eta Histogram",750,500);
    TCanvas *c6 = new TCanvas("canvas_PtvsRapidity","Canvas for Pt vs Rapidity Histogram",750,500);
    TCanvas *c7 = new TCanvas("canvas_PtvsMass","Canvas for Pt vs Mass Histogram",750,500);
    
    //histograms
    
    TH1F *hist_PT_n2n2j = new TH1F("hist_PT_n2n2j", "PT of n2n2j", PTbins, 0.0, PTxaxis);
    TH1F *hist_Eta_n2n2j = new TH1F("hist_Eta_n2n2j", "Eta of n2n2j", Etabins, -Etaxaxis, Etaxaxis);
    TH1F *hist_Rapidity_n2n2j = new TH1F("hist_Rapidity_n2n2j", "Rapidity of n2n2j", Rapiditybins, -Rapidityxaxis, Rapidityxaxis);
    TH1F *hist_Mass_n2n2j = new TH1F("hist_Mass_n2n2j", "Mass of n2n2j", Massbins, 400.0, Massxaxis);
    TH2F *hist_PtvsEta_n2n2j = new TH2F("hist_PtvsEta_n2n2j","Pt vs Eta of n2n2j", PTbins, 0.0, PTxaxis, Etabins, -Etaxaxis, Etaxaxis);
    TH2F *hist_PtvsRapidity_n2n2j = new TH2F("hist_PtvsRapidity_n2n2j","Pt vs Rapidity of n2n2j", PTbins, 0.0, PTxaxis, Rapiditybins, -Rapidityxaxis, Rapidityxaxis);
    TH2F *hist_PtvsMass_n2n2j = new TH2F("hist_PtvsMass_n2n2j","Pt vs Mass of n2n2j", PTbins, 0.0, PTxaxis, Massbins, 400.0, Massxaxis);
    
    //set title of axis on histograms
    
    hist_PtvsEta_n2n2j->GetYaxis()->SetTitle("Eta");
    hist_PtvsEta_n2n2j->GetXaxis()->SetTitle("Pt [GeV]");
    hist_PtvsRapidity_n2n2j->GetYaxis()->SetTitle("Rapidity");
    hist_PtvsRapidity_n2n2j->GetXaxis()->SetTitle("Pt [GeV]");
    hist_PtvsMass_n2n2j->GetYaxis()->SetTitle("Mass [GeV]");
    hist_PtvsMass_n2n2j->GetXaxis()->SetTitle("Pt [GeV]");
    
    //Uncomment below to set the range on the y axis to some specified value
    
    //hist_PT_n2n2j->GetYaxis()->SetRangeUser(0.0,PTyaxis);
    //hist_Eta_n2n2j->GetYaxis()->SetRangeUser(0.0,Etayaxis);
    //hist_Rapidity_n2n2j->GetYaxis()->SetRangeUser(0.0,Rapidityyaxis);
    //hist_Mass_n2n2j->GetYaxis()->SetRangeUser(0.0,Massyaxis);
    
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) { //event loop
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        
        for(int i=0; i<Event_Nparticles[0]; i++) //particle loop
        {
            if(Particle_PID[i]==1000023 && Particle_Status[i]==2)
            {
                if(n2A_MC==zeroTLV) //use this to see if the first neutralino has been identified. if not then set a TLV for it
                {
                    n2A_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
                else //else then we are past the first neutralino and will fill the next
                {
                    n2B_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
            }
            if(Particle_PID[i]==1000022 && Particle_Status[i]==2)
            {
                if(n1A_MC==zeroTLV)
                {
                    n1A_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
                else
                {
                    n1B_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
            }
            if(Particle_PID[Particle_Mother1[i]]==23 && (abs(Particle_PID[i])==11 || abs(Particle_PID[i])==13)) //similar filling process to the neutralinos
            {
                if(lA_MC==zeroTLV)
                {
                    lA_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
                else if(lB_MC==zeroTLV)
                {
                    lB_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
                else if(lC_MC==zeroTLV)
                {
                    lC_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
                else
                {
                    lD_MC.SetPtEtaPhiE(Particle_PT[i],Particle_Eta[i],Particle_Phi[i],Particle_E[i]);
                }
            }
        }
        
        //check to make sure TLVs are unique
        if(n2A_MC==n2B_MC)
        {
            cout << "FAIL N2" << endl;
        }
        if(n1A_MC==n1B_MC)
        {
            cout << "FAIL N1" << endl;
        }
        if(lA_MC==lB_MC || lA_MC==lC_MC || lA_MC==lD_MC || lB_MC==lC_MC || lB_MC==lD_MC || lC_MC==lD_MC)
        {
            cout << "FAIL l" << endl;
        }
        
        n2=n2A+n2B; //sum the n2 and fill histograms
        hist_PT_n2n2j->Fill(n2.Pt());
        hist_Eta_n2n2j->Fill(n2.Eta());
        hist_Rapidity_n2n2j->Fill(n2.Rapidity());
        hist_Mass_n2n2j->Fill(n2.M());
        hist_PtvsEta_n2n2j->Fill(n2.Pt(), n2.Eta());
        hist_PtvsRapidity_n2n2j->Fill(n2.Pt(), n2.Rapidity());
        hist_PtvsMass_n2n2j->Fill(n2.Pt(), n2.M());
        
        //reset TLVs
        n2A_MC=zeroTLV;
        n2B_MC=zeroTLV;
        n1A_MC=zeroTLV;
        n1B_MC=zeroTLV;
        lA_MC=zeroTLV;
        lB_MC=zeroTLV;
        lC_MC=zeroTLV;
        lD_MC=zeroTLV;
    }
    //scale histograms and set colors
    hist_PT_n2n2j->Scale(1.0/hist_PT_n2n2j->Integral("width"));
    hist_PT_n2n2j->SetLineColor(kBlue);
    hist_Eta_n2n2j->Scale(1.0/hist_Eta_n2n2j->Integral("width"));
    hist_Eta_n2n2j->SetLineColor(kBlue);
    hist_Rapidity_n2n2j->Scale(1.0/hist_Rapidity_n2n2j->Integral("width"));
    hist_Rapidity_n2n2j->SetLineColor(kBlue);
    hist_Mass_n2n2j->Scale(1.0/hist_Mass_n2n2j->Integral("width"));
    hist_Mass_n2n2j->SetLineColor(kBlue);
    output->Write(); //write tree
    //for 2D histograms
    c5->cd();
    hist_PtvsEta_n2n2j->Draw("COLZ");
    c5->SaveAs("PtvsEta_n2n2j.pdf");
    c6->cd();
    hist_PtvsRapidity_n2n2j->Draw("COLZ");
    c6->SaveAs("PtvsRapidity_n2n2j.pdf");
    c7->cd();
    hist_PtvsMass_n2n2j->Draw("COLZ");
    c7->SaveAs("PtvsMass_n2n2j.pdf");
    output->Close(); //close TFile
}
