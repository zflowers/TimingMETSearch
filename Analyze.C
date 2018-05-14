/* ----------------------------------------------------------------------------
Author: Zach Flowers
Description:  Combines the output histograms from the jet files
Date: April 23rd 2018
* ---------------------------------------------------------------------------- */

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include "TLegend.h"

void Analyze()
{
    gROOT->Clear();
    gStyle->SetOptStat(0);
    bool high=true;
    //bool high=false;
    
    TFile *f1 = new TFile("output.root", "READ"); //open the root file with the histograms
    //Make a Canvas and pad for the PT histograms
    TCanvas *c1 = new TCanvas("Merged_PT_hist", "Canvas for Merged PT Hist", 800, 800);
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
    pad1->SetBottomMargin(0.07);
    //pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
    
    //Get the histograms for merging from the root file
    TH1F* h1_PT = (TH1F*)f1->Get("hist_PT_n2n2");
    TH1F* h2_PT = (TH1F*)f1->Get("hist_PT_n2n2j");
    TH1F* h3_PT = (TH1F*)f1->Get("hist_PT_n2n2jj");
    h1_PT->SetDirectory(0);
    h2_PT->SetDirectory(0);
    h3_PT->SetDirectory(0);
    h1_PT->SetTitle("Pt of Each Process");
    h2_PT->SetTitle("");
    h3_PT->SetTitle("");
    h1_PT->GetXaxis()->SetTitle("Pt [GeV]");
    h1_PT->GetYaxis()->SetTitle("Number of Events");
    if(high==true)
    {
        h1_PT->GetXaxis()->SetRangeUser(420.0,1200.0);
        h2_PT->GetXaxis()->SetRangeUser(420.0,1200.0);
        h3_PT->GetXaxis()->SetRangeUser(420.0,1200.0);
    }
    h1_PT->SetStats(false);
    h2_PT->SetStats(false);
    h3_PT->SetStats(false);
    
    //Make histograms for the ratio plots
    
    TH1F* h1_res_PT = (TH1F*)h1_PT->Clone("h1_res_PT");
    TH1F* h2_res_PT = (TH1F*)h2_PT->Clone("h2_res_PT");
    TH1F* h3_res_PT = (TH1F*)h3_PT->Clone("h3_res_PT");
    
    h1_res_PT->SetDirectory(0);
    h2_res_PT->SetDirectory(0);
    h3_res_PT->SetDirectory(0);
    
    //Divide by the histogram with 1 jet
    h1_res_PT->Divide(h2_PT);
    h2_res_PT->Divide(h2_PT);
    h3_res_PT->Divide(h2_PT);
    
    //Turn off the stats box
    h1_res_PT->SetStats(false);
    h2_res_PT->SetStats(false);
    h3_res_PT->SetStats(false);
    
    //Create a legend for the combined histograms
    TLegend* leg_PT = new TLegend(0.7,0.7,.9,.9, "");
    leg_PT->SetTextFont(42);
    leg_PT->SetTextSize(0.04);
    leg_PT->SetFillColor(kWhite);
    leg_PT->SetTextColor(kBlack);
    leg_PT->AddEntry(h1_PT,"n2n2","L");
    leg_PT->AddEntry(h2_PT,"n2n2j","L");
    leg_PT->AddEntry(h3_PT,"n2n2jj","L");
    h1_PT->Draw();
    h2_PT->Draw("SAME");
    h3_PT->Draw("SAME");
    leg_PT->Draw("SAME");
    
    //Make the pad for the residuals
    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.05);
    pad2->Update();

    if(high==true)
    {
        h1_res_PT->GetYaxis()->SetRangeUser(0.0,8.0);
        h2_res_PT->GetYaxis()->SetRangeUser(0.0,8.0);
        h3_res_PT->GetYaxis()->SetRangeUser(0.0,8.0);
    }
    //Draw the residuals
    h1_res_PT->Draw();
    h2_res_PT->Draw("SAME");
    h3_res_PT->Draw("SAME");
    pad2->Update();
    
    //Save the outputs
    if(high==true)
    {
        c1->SaveAs("PTResidual_Truncated_highPT.pdf");
    }
    else
    {
        c1->SaveAs("PTResidual.pdf");
    }
    //Same setup as before but now on a log scale
    TCanvas *c2 = new TCanvas("Merged_PT_histLOG", "Canvas for Merged PT HistLOG", 800, 800);
    TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1.0);
    pad3->SetBottomMargin(0.07);
    pad3->SetLogy(); //Change to Log scale
    pad3->Draw();
    pad3->cd();
    
    h1_PT->Draw();
    h2_PT->Draw("SAME");
    h3_PT->Draw("SAME");
    leg_PT->Draw("SAME");
    
    c2->cd();
    TPad *pad4 = new TPad("pad4", "pad4", 0, 0.05, 1, 0.3);
    pad4->Draw();
    pad4->cd();
    pad4->SetTopMargin(0.05);
    pad4->Update();
    
    h1_res_PT->Draw();
    h2_res_PT->Draw("SAME");
    h3_res_PT->Draw("SAME");
    pad4->Update();
    if(high==true)
    {
        c2->SaveAs("PTResidual_Truncated_highPTLOG.pdf");
    }
    else
    {
        c2->SaveAs("PTResidual_TruncatedLOG.pdf");
    }
        
    TCanvas *c3 = new TCanvas("Merged_Eta_hist", "Canvas for Merged Eta Hist", 800, 800);
    
    TPad *pad5 = new TPad("pad5","pad5",0,0.3,1,1.0);
    pad5->SetBottomMargin(0.07);
    //pad5->SetLogy();
    pad5->Draw();
    pad5->cd();
    
    TH1F* h1_Eta = (TH1F*)f1->Get("hist_Eta_n2n2");
    TH1F* h2_Eta = (TH1F*)f1->Get("hist_Eta_n2n2j");
    TH1F* h3_Eta = (TH1F*)f1->Get("hist_Eta_n2n2jj");
    h1_Eta->SetDirectory(0);
    h2_Eta->SetDirectory(0);
    h3_Eta->SetDirectory(0);
    h1_Eta->SetTitle("Eta of Each Process");
    h2_Eta->SetTitle("");
    h3_Eta->SetTitle("");
    h1_Eta->GetXaxis()->SetTitle("Eta");
    h1_Eta->GetYaxis()->SetTitle("Number of Events");
    h1_Eta->SetStats(false);
    h2_Eta->SetStats(false);
    h3_Eta->SetStats(false);
    
    TH1F* h1_res_Eta = (TH1F*)h1_Eta->Clone("h1_res_Eta");
    TH1F* h2_res_Eta = (TH1F*)h2_Eta->Clone("h2_res_Eta");
    TH1F* h3_res_Eta = (TH1F*)h3_Eta->Clone("h3_res_Eta");
    
    h1_res_Eta->SetDirectory(0);
    h2_res_Eta->SetDirectory(0);
    h3_res_Eta->SetDirectory(0);
    
    h1_res_Eta->Divide(h2_Eta);
    h2_res_Eta->Divide(h2_Eta);
    h3_res_Eta->Divide(h2_Eta);
    
    h1_res_Eta->SetStats(false);
    h2_res_Eta->SetStats(false);
    h3_res_Eta->SetStats(false);
    
    TLegend* leg_Eta = new TLegend(0.7,0.7,.9,.9, "");
    leg_Eta->SetTextFont(42);
    leg_Eta->SetTextSize(0.04);
    leg_Eta->SetFillColor(kWhite);
    leg_Eta->SetTextColor(kBlack);
    leg_Eta->AddEntry(h1_Eta,"n2n2","L");
    leg_Eta->AddEntry(h2_Eta,"n2n2j","L");
    leg_Eta->AddEntry(h3_Eta,"n2n2jj","L");
    h1_Eta->Draw();
    h2_Eta->Draw("SAME");
    h3_Eta->Draw("SAME");
    leg_Eta->Draw("SAME");
    
    c3->cd();
    TPad *pad6 = new TPad("pad6", "pad6", 0, 0.05, 1, 0.3);
    pad6->Draw();
    pad6->cd();
    pad6->SetTopMargin(0.05);
    pad6->Update();
    
    h1_res_Eta->Draw();
    h2_res_Eta->Draw("SAME");
    h3_res_Eta->Draw("SAME");
    pad6->Update();
    c3->SaveAs("EtaResidual.pdf");
    
    TCanvas *c4 = new TCanvas("Merged_Eta_histLOG", "Canvas for Merged Eta HistLOG", 800, 800);
    TPad *pad7 = new TPad("pad7","pad7",0,0.3,1,1.0);
    pad7->SetBottomMargin(0.07);
    pad7->SetLogy();
    pad7->Draw();
    pad7->cd();
    
    h1_Eta->Draw();
    h2_Eta->Draw("SAME");
    h3_Eta->Draw("SAME");
    leg_Eta->Draw("SAME");
    
    c4->cd();
    TPad *pad8 = new TPad("pad8", "pad8", 0, 0.05, 1, 0.3);
    pad8->Draw();
    pad8->cd();
    pad8->SetTopMargin(0.05);
    pad8->Update();
    
    h1_res_Eta->Draw();
    h2_res_Eta->Draw("SAME");
    h3_res_Eta->Draw("SAME");
    pad8->Update();
    c4->SaveAs("EtaResidualLOG.pdf");
    
    TCanvas *c5 = new TCanvas("Merged_Rapidity_hist", "Canvas for Merged Rapidity Hist", 800, 800);
    
    TPad *pad9 = new TPad("pad9","pad9",0,0.3,1,1.0);
    pad9->SetBottomMargin(0.07);
    //pad9->SetLogy();
    pad9->Draw();
    pad9->cd();
    
    TH1F* h1_Rapidity = (TH1F*)f1->Get("hist_Rapidity_n2n2");
    TH1F* h2_Rapidity = (TH1F*)f1->Get("hist_Rapidity_n2n2j");
    TH1F* h3_Rapidity = (TH1F*)f1->Get("hist_Rapidity_n2n2jj");
    h1_Rapidity->SetDirectory(0);
    h2_Rapidity->SetDirectory(0);
    h3_Rapidity->SetDirectory(0);
    h1_Rapidity->SetTitle("Rapidity of Each Process");
    h2_Rapidity->SetTitle("");
    h3_Rapidity->SetTitle("");
    h1_Rapidity->GetXaxis()->SetTitle("Rapidity");
    h1_Rapidity->GetYaxis()->SetTitle("Number of Events");
    h1_Rapidity->SetStats(false);
    h2_Rapidity->SetStats(false);
    h3_Rapidity->SetStats(false);
    
    TH1F* h1_res_Rapidity = (TH1F*)h1_Rapidity->Clone("h1_res_Rapidity");
    TH1F* h2_res_Rapidity = (TH1F*)h2_Rapidity->Clone("h2_res_Rapidity");
    TH1F* h3_res_Rapidity = (TH1F*)h3_Rapidity->Clone("h3_res_Rapidity");
    
    h1_res_Rapidity->SetDirectory(0);
    h2_res_Rapidity->SetDirectory(0);
    h3_res_Rapidity->SetDirectory(0);
    
    h1_res_Rapidity->Divide(h2_Rapidity);
    h2_res_Rapidity->Divide(h2_Rapidity);
    h3_res_Rapidity->Divide(h2_Rapidity);
    
    h1_res_Rapidity->SetStats(false);
    h2_res_Rapidity->SetStats(false);
    h3_res_Rapidity->SetStats(false);
    
    TLegend* leg_Rapidity = new TLegend(0.7,0.7,.9,.9, "");
    leg_Rapidity->SetTextFont(42);
    leg_Rapidity->SetTextSize(0.04);
    leg_Rapidity->SetFillColor(kWhite);
    leg_Rapidity->SetTextColor(kBlack);
    leg_Rapidity->AddEntry(h1_Rapidity,"n2n2","L");
    leg_Rapidity->AddEntry(h2_Rapidity,"n2n2j","L");
    leg_Rapidity->AddEntry(h3_Rapidity,"n2n2jj","L");
    h1_Rapidity->Draw();
    h2_Rapidity->Draw("SAME");
    h3_Rapidity->Draw("SAME");
    leg_Rapidity->Draw("SAME");
    
    c5->cd();
    TPad *pad10 = new TPad("pad10", "pad10", 0, 0.05, 1, 0.3);
    pad10->Draw();
    pad10->cd();
    pad10->SetTopMargin(0.05);
    pad10->Update();
    
    h1_res_Rapidity->Draw();
    h2_res_Rapidity->Draw("SAME");
    h3_res_Rapidity->Draw("SAME");
    pad10->Update();
    c5->SaveAs("RapidityResidual.pdf");
    
    TCanvas *c6 = new TCanvas("Merged_Rapidity_histLOG", "Canvas for Merged Rapidity HistLOG", 800, 800);
    
    TPad *pad11 = new TPad("pad11","pad11",0,0.3,1,1.0);
    pad11->SetBottomMargin(0.07);
    pad11->SetLogy();
    pad11->Draw();
    pad11->cd();
    
    h1_Rapidity->Draw();
    h2_Rapidity->Draw("SAME");
    h3_Rapidity->Draw("SAME");
    leg_Rapidity->Draw("SAME");
    
    c6->cd();
    TPad *pad12 = new TPad("pad12", "pad12", 0, 0.05, 1, 0.3);
    pad12->Draw();
    pad12->cd();
    pad12->SetTopMargin(0.05);
    pad12->Update();
    
    h1_res_Rapidity->Draw();
    h2_res_Rapidity->Draw("SAME");
    h3_res_Rapidity->Draw("SAME");
    pad12->Update();
    c6->SaveAs("RapidityResidualLOG.pdf");
    
    TCanvas *c7 = new TCanvas("Merged_Mass_hist", "Canvas for Merged Mass Hist", 800, 800);
    
    TPad *pad13 = new TPad("pad13","pad13",0,0.3,1,1.0);
    pad13->SetBottomMargin(0.07);
    //pad13->SetLogy();
    pad13->Draw();
    pad13->cd();
    
    TH1F* h1_Mass = (TH1F*)f1->Get("hist_Mass_n2n2");
    TH1F* h2_Mass = (TH1F*)f1->Get("hist_Mass_n2n2j");
    TH1F* h3_Mass = (TH1F*)f1->Get("hist_Mass_n2n2jj");
    h1_Mass->SetDirectory(0);
    h2_Mass->SetDirectory(0);
    h3_Mass->SetDirectory(0);
    h1_Mass->SetTitle("Invariant Mass of Each Process");
    h2_Mass->SetTitle("");
    h3_Mass->SetTitle("");
    h1_Mass->GetXaxis()->SetTitle("Mass [GeV]");
    h1_Mass->GetYaxis()->SetTitle("Number of Events");
    h1_Mass->SetStats(false);
    h2_Mass->SetStats(false);
    h3_Mass->SetStats(false);
    
    TH1F* h1_res_Mass = (TH1F*)h1_Mass->Clone("h1_res_Mass");
    TH1F* h2_res_Mass = (TH1F*)h2_Mass->Clone("h2_res_Mass");
    TH1F* h3_res_Mass = (TH1F*)h3_Mass->Clone("h3_res_Mass");
    
    h1_res_Mass->SetDirectory(0);
    h2_res_Mass->SetDirectory(0);
    h3_res_Mass->SetDirectory(0);
    
    h1_res_Mass->Divide(h2_Mass);
    h2_res_Mass->Divide(h2_Mass);
    h3_res_Mass->Divide(h2_Mass);
    
    h1_res_Mass->SetStats(false);
    h2_res_Mass->SetStats(false);
    h3_res_Mass->SetStats(false);
    
    TLegend* leg_Mass = new TLegend(0.7,0.7,.9,.9, "");
    leg_Mass->SetTextFont(42);
    leg_Mass->SetTextSize(0.04);
    leg_Mass->SetFillColor(kWhite);
    leg_Mass->SetTextColor(kBlack);
    leg_Mass->AddEntry(h1_Mass,"n2n2","L");
    leg_Mass->AddEntry(h2_Mass,"n2n2j","L");
    leg_Mass->AddEntry(h3_Mass,"n2n2jj","L");
    h1_Mass->Draw();
    h2_Mass->Draw("SAME");
    h3_Mass->Draw("SAME");
    leg_Mass->Draw("SAME");
    
    c7->cd();
    TPad *pad14 = new TPad("pad14", "pad14", 0, 0.05, 1, 0.3);
    pad14->Draw();
    pad14->cd();
    pad14->SetTopMargin(0.05);
    pad14->Update();
    
    h1_res_Mass->Draw();
    h2_res_Mass->Draw("SAME");
    h3_res_Mass->Draw("SAME");
    pad14->Update();
    c7->SaveAs("MassResidual.pdf");
    
    TCanvas *c8 = new TCanvas("Merged_Mass_histLOG", "Canvas for Merged Mass HistLOG", 800, 800);
    
    TPad *pad15 = new TPad("pad15","pad15",0,0.3,1,1.0);
    pad15->SetBottomMargin(0.07);
    pad15->SetLogy();
    pad15->Draw();
    pad15->cd();
    
    h1_Mass->Draw();
    h2_Mass->Draw("SAME");
    h3_Mass->Draw("SAME");
    leg_Mass->Draw("SAME");
    
    c8->cd();
    TPad *pad16 = new TPad("pad16", "pad16", 0, 0.05, 1, 0.3);
    pad16->Draw();
    pad16->cd();
    pad16->SetTopMargin(0.05);
    pad16->Update();
    
    h1_res_Mass->Draw();
    h2_res_Mass->Draw("SAME");
    h3_res_Mass->Draw("SAME");
    pad16->Update();
    c8->SaveAs("MassResidualLOG.pdf");
    
    f1->Close();
}
