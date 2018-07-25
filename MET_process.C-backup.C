//01 jet
#define MET_process_cxx
#include "MET_process.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../Detector.hh"
#include "Bonus.h"

void MET_process::Loop()
{
    TFile* Plots = new TFile("PlotsMET.root","UPDATE");
    CreatePalette();
    TLorentzVector zeroTLV;
    zeroTLV.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
    TLorentzVector n2A_MC, n2B_MC, n1A_MC, n1B_MC, eA_MC, eB_MC, eC_MC, eD_MC, mA_MC, mB_MC, mC_MC, mD_MC;
    TLorentzVector PV_RECO, SV_A_RECO, SV_B_RECO;
    TLorentzVector eA_RECO, eB_RECO, eC_RECO, eD_RECO, mA_RECO, mB_RECO, mC_RECO, mD_RECO;
    
    TCanvas *c1 = new TCanvas("METCanvasforSumPTandPT_n2n2","MET Canvas for SumPT and PT_n2n2",750,500);
    TCanvas *c2 = new TCanvas("METCanvasforSumPTdividedbyPT_n2n2","MET Canvas for SumPT divided by PT_n2n2",750,500);
    TCanvas *c3 = new TCanvas("METCanvasforSumPTdividedbyPT_n2n2vsPT_n2n2","MET Canvas for SumPT divided by PT_n2n2 vs PT_n2n2",750,500);
    TCanvas *c4 = new TCanvas("METCanvasforPT_n2n2","MET Canvas for PT_n2n2",750,500);
    TCanvas *c5 = new TCanvas("METCanvasforMETResolution_SumPT","MET Canvas for MET Resolution SumPT",750,500);
    TCanvas *c6 = new TCanvas("METCanvasforMETResolution_PT_n2n2","MET Canvas for MET Resolution PT_n2n2",750,500);
    
    
    TH2F *hist_SumPT_vs_PT_n2n2 = new TH2F("hist_SumPT_vs_PT_n2n2","SumPT vs PT of n2n2", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
    //TH2F *hist_SumPT_vs_PT_n2n2 = new TH2F("hist_SumPT_vs_PT_n2n2","SumPT vs PT of n2n2", 100, 0.0, 500.0, 100, 0.0, 1200.0);
    TH1D *hist_SumPT_over_PT_n2n2 = new TH1D("hist_SumPT_over_PT_n2n2","SumPT divided by PT of n2n2",100,0.0,10.0);
    TH2F *hist_SumPT_over_PT_n2n2_vs_PT_n2n2 = new TH2F("hist_SumPT_over_PT_n2n2_vs_PT_n2n2","SumPT divided by PT_n2n2 vs PT of n2n2", 100, 0.0, 700.0, 100, 0.0, 20.0);
    TH1D *hist_PT_n2n2 = new TH1D("hist_PT_n2n2","PT of chi2 chi2 system",100,0.0,700);
    TH1F *hist_True_SumPT = new TH1F("hist_True_SumPT","True SumPT",100,0.0,1500.0);
    TH1F *hist_Test_SumPT = new TH1F("hist_Test_SumPT","Test SumPT",100,0.0,1500.0);
    
    double sum_PT=0.0;
    zeroTLV.SetPtEtaPhiE(0.0,0.0,0.0,0.0);
    gRandom = new TRandom3();
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    //double Decay_Time=1.0; //need to decide on units
    double MET_Parallel_Mag=0.0;
    double MET_Perpendicular_Mag=0.0;
    TVector3 MET_Perpendicular;
    TVector3 zhat;
    zhat.SetXYZ(0.0,0.0,1.0);
    
    if (fChain == 0) return;
    
    std::vector<TH1D*> vect_hist_MET_Parallel_SumPT;
    std::vector<TH1D*> vect_hist_MET_Perpendicular_SumPT;
    std::vector<TH1D*> vect_hist_MET_Parallel_PT_n2n2;
    std::vector<TH1D*> vect_hist_MET_Perpendicular_PT_n2n2;
    std::vector<double> vect_n2n2_Pt;
    std::vector<double> vect_SumPt;
    std::vector<TCanvas*> vect_canvas_MET_Parallel_SumPT;
    std::vector<TCanvas*> vect_canvas_MET_Perpendicular_SumPT;
    std::vector<TCanvas*> vect_canvas_MET_Parallel_PT_n2n2;
    std::vector<TCanvas*> vect_canvas_MET_Perpendicular_PT_n2n2;
    std::vector<double> vect_MET_Parallel;
    std::vector<double> vect_MET_Perpendicular;
    std::vector<TH1D*> vect_hist_sum_PT_over_PT_n2n2;
    std::vector<TCanvas*> vect_canvas_sum_PT_over_PT_n2n2;
    std::vector<double> vect_sum_PT_over_PT_n2n2;
    std::vector<TH1D*> vect_hist_sum_PT_over_scalar_PT_n2n2;
    std::vector<TCanvas*> vect_canvas_sum_PT_over_scalar_PT_n2n2;
    std::vector<double> vect_sum_PT_over_scalar_PT_n2n2;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    double arr_sum_PT_over_PT_n2n2[nentries];
    double arr_PT_n2n2[nentries];
    double arr_sum_PT[nentries];
    
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        
        sum_PT = 0.0;
        MET_Parallel_Mag=0.0;
        MET_Perpendicular_Mag=0.0;
        
        for(int i=0; i<Jet_Delphes_Size; i++)
        {
            sum_PT+=Jet_Delphes_PT->at(i);
        }
         
        n2A_MC.SetPtEtaPhiE(n2_MC_PT->at(0), n2_MC_Eta->at(0), n2_MC_Phi->at(0), n2_MC_E->at(0));
        n2B_MC.SetPtEtaPhiE(n2_MC_PT->at(1), n2_MC_Eta->at(1), n2_MC_Phi->at(1), n2_MC_E->at(1));
        
        n1A_MC.SetPtEtaPhiE(n1_MC_PT->at(0), n1_MC_Eta->at(0), n1_MC_Phi->at(0), n1_MC_E->at(0));
        n1B_MC.SetPtEtaPhiE(n1_MC_PT->at(1), n1_MC_Eta->at(1), n1_MC_Phi->at(1), n1_MC_E->at(1));

        Detector test;
        //TVector3 MET = test.Smear_MET(sum_PT);
        
        TLorentzVector PV_RECO, SV_A_RECO, SV_B_RECO;
        TLorentzVector eA_RECO, eB_RECO, eC_RECO, eD_RECO, mA_RECO, mB_RECO, mC_RECO, mD_RECO;
        
        TLorentzVector n2n2 = n2A_MC + n2B_MC;
        
        TVector3 PT_true = (n1A_MC + n1B_MC).Vect();
        PT_true.SetZ(0.0);
        
        TVector3 MET_Delphes;
        MET_Delphes.SetPtEtaPhi(MET_PT,MET_Eta,MET_Phi);
        MET_Delphes.SetZ(0.0);
        
        MET_Perpendicular_Mag=PT_true.Cross(zhat).Unit().Dot(MET_Delphes);
        MET_Parallel_Mag=(MET_Delphes.Dot(PT_true))/PT_true.Mag();
        
        vect_MET_Parallel.push_back(MET_Parallel_Mag-PT_true.Mag());
        vect_MET_Perpendicular.push_back(MET_Perpendicular_Mag);
        
        //double n2n2_Pt = n2n2.Pt();
        vect_n2n2_Pt.push_back(n2n2.Pt());
        vect_SumPt.push_back(sum_PT);
        vect_sum_PT_over_PT_n2n2.push_back(sum_PT/n2n2.Pt());
        vect_sum_PT_over_scalar_PT_n2n2.push_back(sum_PT/(n2A_MC.Pt()+n2B_MC.Pt()));
        
        hist_SumPT_vs_PT_n2n2->Fill(sum_PT,n2n2.Pt());
        hist_SumPT_over_PT_n2n2->Fill(sum_PT/n2n2.Pt());
        hist_PT_n2n2->Fill(n2n2.Pt());
        arr_sum_PT_over_PT_n2n2[jentry]=sum_PT/n2n2.Pt();
        arr_PT_n2n2[jentry]=n2n2.Pt();
        arr_sum_PT[jentry]=sum_PT;
        hist_SumPT_over_PT_n2n2_vs_PT_n2n2->Fill(n2n2.Pt(), sum_PT/n2n2.Pt());
        //hist_SumPT_over_PT_n2n2_vs_PT_n2n2->Fill(sum_PT, sum_PT/n2n2.Pt());
    }
    cout << "Looped Through Events" << endl;
    //bool Save=false;
    bool Save=true;
    //vector<double> vBinEdge_SumPT;
    vector<double> vBinEdge_PT_n2n2;
    
    //Bins for SumPT
    /*
    int Ninterval_SumPT = 17;
    vBinEdge_SumPT.push_back(0.);
    vBinEdge_SumPT.push_back(400.);
    vBinEdge_SumPT.push_back(450.);
    vBinEdge_SumPT.push_back(500.);
    vBinEdge_SumPT.push_back(525.);
    vBinEdge_SumPT.push_back(550.);
    vBinEdge_SumPT.push_back(575.);
    vBinEdge_SumPT.push_back(600.);
    vBinEdge_SumPT.push_back(625.);
    vBinEdge_SumPT.push_back(650.);
    vBinEdge_SumPT.push_back(675.);
    vBinEdge_SumPT.push_back(700.);
    vBinEdge_SumPT.push_back(750.);
    vBinEdge_SumPT.push_back(800.);
    vBinEdge_SumPT.push_back(900.);
    vBinEdge_SumPT.push_back(1000.);
    vBinEdge_SumPT.push_back(1200.);
     */
    
    //Bins for PT_n2n2
    int Ninterval_PT_n2n2 = 16;
    vBinEdge_PT_n2n2.push_back(0.);
    vBinEdge_PT_n2n2.push_back(25.);
    vBinEdge_PT_n2n2.push_back(50.);
    vBinEdge_PT_n2n2.push_back(75.);
    vBinEdge_PT_n2n2.push_back(100.);
    vBinEdge_PT_n2n2.push_back(125.);
    vBinEdge_PT_n2n2.push_back(150.);
    vBinEdge_PT_n2n2.push_back(175.);
    vBinEdge_PT_n2n2.push_back(200.);
    vBinEdge_PT_n2n2.push_back(225.);
    vBinEdge_PT_n2n2.push_back(250.);
    vBinEdge_PT_n2n2.push_back(275.);
    vBinEdge_PT_n2n2.push_back(300.);
    vBinEdge_PT_n2n2.push_back(350.);
    vBinEdge_PT_n2n2.push_back(400.);
    vBinEdge_PT_n2n2.push_back(475.);
    vBinEdge_PT_n2n2.push_back(600.);
    
    
    /*vector<double> vBinSum_SumPT;
    vector<double> vBinSum2_SumPT;
    vector<double> vBinN_SumPT;
    */
    vector<double> vBinSum_PT_n2n2;
    vector<double> vBinSum2_PT_n2n2;
    vector<double> vBinN_PT_n2n2;
    
    /*for(int i=0; i<Ninterval_SumPT; i++)
    {
        vBinSum_SumPT.push_back(0.);
        vBinSum2_SumPT.push_back(0.);
        vBinN_SumPT.push_back(0.);
        
        //SumPT
        char* canvasname_Parallel_SumPT = new char[100];
        std::sprintf(canvasname_Parallel_SumPT,"MET_canvas_MET_Parallel_SumPT_%d",i);
        string canvasname_Parallel_SumPT_str = "MET_canvas_MET_Parallel_SumPT"+std::to_string(i);
        TCanvas* dummy_canvas_Parallel_SumPT = new TCanvas((canvasname_Parallel_SumPT_str).c_str(),(canvasname_Parallel_SumPT_str).c_str(),750,500);
        vect_canvas_MET_Parallel_SumPT.push_back(dummy_canvas_Parallel_SumPT);
        
        char* canvasname_Perpendicular_SumPT = new char[100];
        std::sprintf(canvasname_Perpendicular_SumPT,"MET_canvas_MET_Perpendicular_SumPT_%d",i);
        string canvasname_Perpendicular_SumPT_str = "MET_canvas_MET_Perpendicular_SumPT"+std::to_string(i);
        TCanvas* dummy_canvas_Perpendicular_SumPT = new TCanvas((canvasname_Perpendicular_SumPT_str).c_str(),(canvasname_Perpendicular_SumPT_str).c_str(),750,500);
        vect_canvas_MET_Perpendicular_SumPT.push_back(dummy_canvas_Perpendicular_SumPT);
        
        char* histname_Parallel_SumPT = new char[100];
        std::sprintf(histname_Parallel_SumPT,"MET_hist_MET_Parallel_SumPT_%d",i);
        string histname_Parallel_SumPT_str = "MET_hist_MET_Parallel_SumPT"+std::to_string(i);
        
        TH1D* dummy_hist_Parallel_SumPT = new TH1D((histname_Parallel_SumPT_str).c_str(),(histname_Parallel_SumPT_str).c_str(),200,-200.0,200.0);
        
        vect_hist_MET_Parallel_SumPT.push_back(dummy_hist_Parallel_SumPT);
        
        char* histname_Perpendicular_SumPT = new char[100];
        std::sprintf(histname_Perpendicular_SumPT,"MET_hist_MET_Perpendicular_SumPT_%d",i);
        string histname_Perpendicular_SumPT_str = "MET_hist_MET_Perpendicular_SumPT"+std::to_string(i);
        
        TH1D* dummy_hist_Perpendicular_SumPT = new TH1D((histname_Perpendicular_SumPT_str).c_str(),(histname_Perpendicular_SumPT_str).c_str(),200,-200.0,200.0);
        
        vect_hist_MET_Perpendicular_SumPT.push_back(dummy_hist_Perpendicular_SumPT);
    }
     */
    for(int i=0; i<Ninterval_PT_n2n2; i++)
    {
        vBinSum_PT_n2n2.push_back(0.);
        vBinSum2_PT_n2n2.push_back(0.);
        vBinN_PT_n2n2.push_back(0.);
        //PT_n2n2
        char* canvasname_Parallel_PT_n2n2 = new char[100];
        std::sprintf(canvasname_Parallel_PT_n2n2,"MET_canvas_MET_Parallel_PT_n2n2_%d",i);
        string canvasname_Parallel_PT_n2n2_str = "MET_canvas_MET_Parallel_PT_n2n2"+std::to_string(i);
        TCanvas* dummy_canvas_Parallel_PT_n2n2 = new TCanvas((canvasname_Parallel_PT_n2n2_str).c_str(),(canvasname_Parallel_PT_n2n2_str).c_str(),750,500);
        vect_canvas_MET_Parallel_PT_n2n2.push_back(dummy_canvas_Parallel_PT_n2n2);
        
        char* canvasname_Perpendicular_PT_n2n2 = new char[100];
        std::sprintf(canvasname_Perpendicular_PT_n2n2,"MET_canvas_MET_Perpendicular_PT_n2n2_%d",i);
        string canvasname_Perpendicular_PT_n2n2_str = "MET_canvas_MET_Perpendicular_PT_n2n2"+std::to_string(i);
        TCanvas* dummy_canvas_Perpendicular_PT_n2n2 = new TCanvas((canvasname_Perpendicular_PT_n2n2_str).c_str(),(canvasname_Perpendicular_PT_n2n2_str).c_str(),750,500);
        vect_canvas_MET_Perpendicular_PT_n2n2.push_back(dummy_canvas_Perpendicular_PT_n2n2);
        
        char* histname_Parallel_PT_n2n2 = new char[100];
        std::sprintf(histname_Parallel_PT_n2n2,"MET_hist_MET_Parallel_PT_n2n2_%d",i);
        string histname_Parallel_PT_n2n2_str = "MET_hist_MET_Parallel_PT_n2n2"+std::to_string(i);
        
        TH1D* dummy_hist_Parallel_PT_n2n2 = new TH1D((histname_Parallel_PT_n2n2_str).c_str(),(histname_Parallel_PT_n2n2_str).c_str(),200,-200.0,200.0);
        
        vect_hist_MET_Parallel_PT_n2n2.push_back(dummy_hist_Parallel_PT_n2n2);
        
        char* histname_Perpendicular_PT_n2n2 = new char[100];
        std::sprintf(histname_Perpendicular_PT_n2n2,"MET_hist_MET_Perpendicular_PT_n2n2_%d",i);
        string histname_Perpendicular_PT_n2n2_str = "MET_hist_MET_Perpendicular_PT_n2n2"+std::to_string(i);
        
        TH1D* dummy_hist_Perpendicular_PT_n2n2 = new TH1D((histname_Perpendicular_PT_n2n2_str).c_str(),(histname_Perpendicular_PT_n2n2_str).c_str(),200,-200.0,200.0);
        
        vect_hist_MET_Perpendicular_PT_n2n2.push_back(dummy_hist_Perpendicular_PT_n2n2);
    }
    
    
    //slice by SumPT
    /*
    for(int k=0; k<nentries; k++)
    {
        for(int i = Ninterval_SumPT-1; i >= 0; i--){
            if(vect_SumPt[k] > vBinEdge_SumPT[i]){
                vBinSum_SumPT[i] += vect_SumPt[k];
                vBinSum2_SumPT[i] += vect_SumPt[k]*vect_SumPt[k];
                vBinN_SumPT[i]++;
                vect_hist_MET_Parallel_SumPT[i]->Fill(vect_MET_Parallel[k]);
                vect_hist_MET_Perpendicular_SumPT[i]->Fill(vect_MET_Perpendicular[k]);
                break;
            }
        }
    }
    */
    //slice by PT_n2n2
    for(int k=0; k<nentries; k++)
    {
        for(int i = Ninterval_PT_n2n2-1; i >= 0; i--){
            if(vect_n2n2_Pt[k] > vBinEdge_PT_n2n2[i]){
                vBinSum_PT_n2n2[i] += vect_n2n2_Pt[k];
                vBinSum2_PT_n2n2[i] += vect_n2n2_Pt[k]*vect_n2n2_Pt[k];
                vBinN_PT_n2n2[i]++;
                vect_hist_MET_Parallel_PT_n2n2[i]->Fill(vect_MET_Parallel[k]);
                vect_hist_MET_Perpendicular_PT_n2n2[i]->Fill(vect_MET_Perpendicular[k]);
                break;
            }
        }
    }
    
    c1->cd();
    hist_SumPT_vs_PT_n2n2->GetXaxis()->SetTitle("sum_PT [GeV]");
    hist_SumPT_vs_PT_n2n2->GetYaxis()->SetTitle("PT of chi2 chi2 system [GeV]");
    hist_SumPT_vs_PT_n2n2->Draw("COLZ");
    if(Save)
        c1->SaveAs("MET_SumPT_vs_n2n2PT.png");
    c1->Write();
    
    c2->cd();
    hist_SumPT_over_PT_n2n2->GetXaxis()->SetTitle("sum_PT divided by PT of chi2 chi2 system");
    hist_SumPT_over_PT_n2n2->GetYaxis()->SetTitle("Number of Entries");
    hist_SumPT_over_PT_n2n2->SetLineColor(kBlack);
    hist_SumPT_over_PT_n2n2->Draw("");
    if(Save)
        c2->SaveAs("MET_SumPT_over_n2n2PT.png");
    c2->Write();
    
    c3->cd();
    hist_SumPT_over_PT_n2n2_vs_PT_n2n2->GetXaxis()->SetTitle("PT of chi2 chi2 system [GeV]");
    //hist_SumPT_over_PT_n2n2_vs_PT_n2n2->GetXaxis()->SetTitle("Sum_PT [GeV]");
    hist_SumPT_over_PT_n2n2_vs_PT_n2n2->GetYaxis()->SetTitle("sum_PT divided by PT of chi2 chi2 system");
    hist_SumPT_over_PT_n2n2_vs_PT_n2n2->Draw("COLZ");
    if(Save)
        c3->SaveAs("MET_SumPT_over_n2n2PT_vs_n2n2PT.png");
    c3->Write();
    
    c4->cd();
    hist_PT_n2n2->GetXaxis()->SetTitle("PT of chi2 chi2 system [GeV]");
    hist_PT_n2n2->GetYaxis()->SetTitle("Number of Entries");
    hist_PT_n2n2->Draw();
    if(Save)
        c4->SaveAs("MET_PT_n2n2.png");
    c4->Write();
    
    /*
    double sigmaMETperp_SumPT[Ninterval_SumPT];
    double sigmaMETpar_SumPT[Ninterval_SumPT];
    double sigmaMETperp_err_SumPT[Ninterval_SumPT];
    double sigmaMETpar_err_SumPT[Ninterval_SumPT];
    */
    double sigmaMETperp_PT_n2n2[Ninterval_PT_n2n2];
    double sigmaMETpar_PT_n2n2[Ninterval_PT_n2n2];
    double sigmaMETperp_err_PT_n2n2[Ninterval_PT_n2n2];
    double sigmaMETpar_err_PT_n2n2[Ninterval_PT_n2n2];
    /*
    double sumPT[Ninterval_SumPT];
    double sumPT_err[Ninterval_SumPT];
    */
    double PT_n2n2[Ninterval_PT_n2n2];
    double PT_n2n2_err[Ninterval_PT_n2n2];
    /*
    for(int i=0; i < Ninterval_SumPT; i++)
    {
        sumPT[i] = vBinSum_SumPT[i]/vBinN_SumPT[i];
        sumPT_err[i] = sqrt( vBinSum2_SumPT[i]/vBinN_SumPT[i] - vBinSum_SumPT[i]*vBinSum_SumPT[i]/vBinN_SumPT[i]/vBinN_SumPT[i])/sqrt(vBinN_SumPT[i]);
        
        //SumPT
        
        vect_canvas_MET_Parallel_SumPT.at(i)->cd();
        vect_hist_MET_Parallel_SumPT.at(i)->Draw();
        cout << endl;
        vect_hist_MET_Parallel_SumPT.at(i)->Fit("gaus","Q");
        vect_hist_MET_Parallel_SumPT.at(i)->GetXaxis()->SetTitle("Magnitude of the Parallel component of MET [GeV]");
        vect_hist_MET_Parallel_SumPT.at(i)->GetYaxis()->SetTitle("Number of Entries");
        string MET_Parallel_Name_SumPT(vect_hist_MET_Parallel_SumPT.at(i)->GetName());
        MET_Parallel_Name_SumPT=MET_Parallel_Name_SumPT+"_MET.png";
        if(Save)
            vect_canvas_MET_Parallel_SumPT.at(i)->SaveAs(MET_Parallel_Name_SumPT.c_str());
        vect_canvas_MET_Parallel_SumPT.at(i)->Update();
        vect_canvas_MET_Parallel_SumPT.at(i)->Write();
        TF1* func = vect_hist_MET_Parallel_SumPT[i]->GetFunction("gaus");
        
        sigmaMETpar_SumPT[i] = func->GetParameter(2)/sumPT[i];
        sigmaMETpar_err_SumPT[i] = func->GetParError(2)/sumPT[i];
        
        vect_canvas_MET_Perpendicular_SumPT.at(i)->cd();
        vect_hist_MET_Perpendicular_SumPT.at(i)->Draw();
        cout << endl;
        vect_hist_MET_Perpendicular_SumPT.at(i)->Fit("gaus","Q");
        vect_hist_MET_Perpendicular_SumPT.at(i)->GetXaxis()->SetTitle("Magnitude of the Perpendicular component of MET [GeV]");
        vect_hist_MET_Perpendicular_SumPT.at(i)->GetYaxis()->SetTitle("Number of Entries");
        string MET_Perpendicular_Name_SumPT(vect_hist_MET_Perpendicular_SumPT.at(i)->GetName());
        MET_Perpendicular_Name_SumPT=MET_Perpendicular_Name_SumPT+"_MET.png";
        if(Save)
            vect_canvas_MET_Perpendicular_SumPT.at(i)->SaveAs(MET_Perpendicular_Name_SumPT.c_str());
        vect_canvas_MET_Perpendicular_SumPT.at(i)->Update();
        vect_canvas_MET_Perpendicular_SumPT.at(i)->Write();
        func = (TF1*) vect_hist_MET_Perpendicular_SumPT[i]->GetFunction("gaus");
        sigmaMETperp_SumPT[i] = func->GetParameter(2)/sumPT[i];
        sigmaMETperp_err_SumPT[i] = func->GetParError(2)/sumPT[i];
    }
     */
    for(int i=0; i<Ninterval_PT_n2n2; i++)
    {
        //PT_n2n2
        PT_n2n2[i] = vBinSum_PT_n2n2[i]/vBinN_PT_n2n2[i];
        PT_n2n2_err[i] = sqrt( vBinSum2_PT_n2n2[i]/vBinN_PT_n2n2[i] - vBinSum_PT_n2n2[i]*vBinSum_PT_n2n2[i]/vBinN_PT_n2n2[i]/vBinN_PT_n2n2[i])/sqrt(vBinN_PT_n2n2[i]);
        vect_canvas_MET_Parallel_PT_n2n2.at(i)->cd();
        vect_hist_MET_Parallel_PT_n2n2.at(i)->Draw();
        cout << endl;
        TF1* func_gaus_PT_n2n2_par_MET = new TF1("func_gaus_PT_n2n2_par_MET","gaus");
        //func_gaus_PT_n2n2_par_MET->FixParameter(1,0.0);
        vect_hist_MET_Parallel_PT_n2n2.at(i)->Fit(func_gaus_PT_n2n2_par_MET,"");
        vect_hist_MET_Parallel_PT_n2n2.at(i)->GetXaxis()->SetTitle("Magnitude of the Parallel component of MET [GeV]");
        vect_hist_MET_Parallel_PT_n2n2.at(i)->GetYaxis()->SetTitle("Number of Entries");
        string MET_Parallel_Name_PT_n2n2(vect_hist_MET_Parallel_PT_n2n2.at(i)->GetName());
        MET_Parallel_Name_PT_n2n2=MET_Parallel_Name_PT_n2n2+"_MET.png";
        if(Save)
            vect_canvas_MET_Parallel_PT_n2n2.at(i)->SaveAs(MET_Parallel_Name_PT_n2n2.c_str());
        vect_canvas_MET_Parallel_PT_n2n2.at(i)->Update();
        vect_canvas_MET_Parallel_PT_n2n2.at(i)->Write();
        
        //TF1* func = vect_hist_MET_Parallel_PT_n2n2[i]->GetFunction("gaus");
        
        sigmaMETpar_PT_n2n2[i] = func_gaus_PT_n2n2_par_MET->GetParameter(2)/PT_n2n2[i];
        sigmaMETpar_err_PT_n2n2[i] = func_gaus_PT_n2n2_par_MET->GetParError(2)/PT_n2n2[i];
        
        vect_canvas_MET_Perpendicular_PT_n2n2.at(i)->cd();
        vect_hist_MET_Perpendicular_PT_n2n2.at(i)->Draw();
        cout << endl;
        TF1* func_gaus_PT_n2n2_perp_MET = new TF1("func_gaus_PT_n2n2_perp_MET","gaus");
        //func_gaus_PT_n2n2_perp_MET->FixParameter(1,0.0);
        vect_hist_MET_Perpendicular_PT_n2n2.at(i)->Fit(func_gaus_PT_n2n2_perp_MET,"");
        vect_hist_MET_Perpendicular_PT_n2n2.at(i)->GetXaxis()->SetTitle("Magnitude of the Perpendicular component of MET [GeV]");
        vect_hist_MET_Perpendicular_PT_n2n2.at(i)->GetYaxis()->SetTitle("Number of Entries");
        string MET_Perpendicular_Name_PT_n2n2(vect_hist_MET_Perpendicular_PT_n2n2.at(i)->GetName());
        MET_Perpendicular_Name_PT_n2n2=MET_Perpendicular_Name_PT_n2n2+"_MET.png";
        if(Save)
            vect_canvas_MET_Perpendicular_PT_n2n2.at(i)->SaveAs(MET_Perpendicular_Name_PT_n2n2.c_str());
        vect_canvas_MET_Perpendicular_PT_n2n2.at(i)->Update();
        vect_canvas_MET_Perpendicular_PT_n2n2.at(i)->Write();
        //func = (TF1*) vect_hist_MET_Perpendicular_PT_n2n2[i]->GetFunction("gaus");
        
        sigmaMETperp_PT_n2n2[i] = func_gaus_PT_n2n2_perp_MET->GetParameter(2)/PT_n2n2[i];
        sigmaMETperp_err_PT_n2n2[i] = func_gaus_PT_n2n2_perp_MET->GetParError(2)/PT_n2n2[i];
    }
    
    //TGraphErrors *gr_perp_SumPT = new TGraphErrors(Ninterval_SumPT, sumPT, sigmaMETperp_SumPT, sumPT_err, sigmaMETperp_err_SumPT);
    
    //TGraphErrors *gr_par_SumPT = new TGraphErrors(Ninterval_SumPT, sumPT, sigmaMETpar_SumPT, sumPT_err, sigmaMETpar_err_SumPT);
    
    TGraphErrors *gr_perp_PT_n2n2 = new TGraphErrors(Ninterval_PT_n2n2, PT_n2n2, sigmaMETperp_PT_n2n2, PT_n2n2_err, sigmaMETperp_err_PT_n2n2);
    
    TGraphErrors *gr_par_PT_n2n2 = new TGraphErrors(Ninterval_PT_n2n2, PT_n2n2, sigmaMETpar_PT_n2n2, PT_n2n2_err, sigmaMETpar_err_PT_n2n2);
    
    /*
    Plot_Res_vs_SumPT(gr_perp_SumPT,gr_par_SumPT,c5);
    if(Save)
        c5->SaveAs("Jet_Resolution_vs_SumPT.png");
    c5->Write();
    */
    Plot_Res_vs_PT_n2n2(gr_perp_PT_n2n2,gr_par_PT_n2n2,c6);
    if(Save)
        c6->SaveAs("Jet_Resolution_vs_PT_n2n2.png");
    c6->Write();
    
    Plots->Close();
    cout << "Finished!" << endl;
}

