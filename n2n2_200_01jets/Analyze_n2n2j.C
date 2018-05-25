#define Analyze_n2n2j_cxx
#include "Analyze_n2n2j.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <string>
#include "DetectorSmear.h"
#include "Bonus.h"

    void Analyze_n2n2j::Loop()
    {
        gRandom = new TRandom3();
        //gStyle->SetOptStat(000002000);
        gStyle->SetOptStat(1);
        gStyle->SetOptFit(1);
        //double Decay_Time=1.0; //need to decide on units
        double MET_Parallel_Mag=0.0;
        double MET_Perpendicular_Mag=0.0;
        TVector3 MET_Perpendicular;
        TVector3 zhat;
        zhat.SetXYZ(0.0,0.0,1.0);
       if (fChain == 0) return;

       Long64_t nentries = fChain->GetEntriesFast();

        TCanvas *c1 = new TCanvas("c1","Canvas for PT",750,500);
        TCanvas *c2 = new TCanvas("c2","Canvas for PT Difference",750,500);
        TCanvas *c3 = new TCanvas("c3","Canvas for MET Parallel",750,500);
        TCanvas *c4 = new TCanvas("c4","Canvas for MET Perpendicular",750,500);
        TCanvas *c5 = new TCanvas("c5","Canvas for Sum PT (MC)",750,500);
        TCanvas *c7 = new TCanvas("c7","Canvas for MET RECO",750,500);
        
        TH2F *hist_pythia_Pt_vs_delphes_Pt_n2n2j = new TH2F("hist_pythia_Pt_vs_delphes_Pt_n2n2j","Pythia Jets Pt vs Delphes Jets Pt", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
        
        TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,-700.0,200.0);
        //TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,0.0,700.0);
        
        TH1D *hist_MET_Parallel = new TH1D("hist_MET_Parallel","MET_Parallel",100,-100.0,100.0);
        TH1D *hist_MET_Perpendicular = new TH1D("hist_MET_Perpendicular","MET_Perpendicular",100,-100.0,100.0);
        TH1D *hist_PT_true = new TH1D("hist_PT_true","Histogram of Sum PT MC",100,0.0,1200.0);
        
        TH1D *hist_MET_RECO = new TH1D("hist_MET_RECO","Smeared PT True (MET)",100,0.0,800.0);
        TH1D *hist_MET_Delphes = new TH1D("hist_MET_Delphes","RECO MET",100,0.0,800.0);
        
        std::vector<TH1D*> vect_hist_MET_Parallel;
        std::vector<TH1D*> vect_hist_MET_Perpendicular;
        std::vector<double> vect_SumPt;
        std::vector<TCanvas*> vect_canvas_MET_Parallel;
        std::vector<TCanvas*> vect_canvas_MET_Perpendicular;
        std::vector<double> vect_MET_Parallel;
        std::vector<double> vect_MET_Perpendicular;
        std::vector<TLorentzVector> vect_Electron_MC;
        std::vector<TLorentzVector> vect_Muon_MC;
        std::vector<TLorentzVector> vect_Electron_RECO;
        std::vector<TLorentzVector> vect_Muon_RECO;
        std::vector<TLorentzVector> vect_SV_MC;
        std::vector<TLorentzVector> vect_SV_RECO;
        
        double arr_Sum_Sum_Pt[17];
        double arr_Total_Sum_Pt[17];
        double arr_Ave_Sum_Pt[17];
        for(int i=0; i<17; i++)
        {
            arr_Sum_Sum_Pt[i]=0.0;
            arr_Total_Sum_Pt[i]=0.0;
        }
        
        double pythia_JET_PT=0.0;
        double delphes_JET_PT=0.0;
        double PT_Difference=0.0;
        
       Long64_t nbytes = 0, nb = 0;
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
          Long64_t ientry = LoadTree(jentry);
          if (ientry < 0) break;
          nb = fChain->GetEntry(jentry);   nbytes += nb;
          // if (Cut(ientry) < 0) continue;
           
           //Setup Muon and Electron TLVs (and smear them)
           for(int i=0; i<Electron_MC_Size; i++)
           {
               TLorentzVector Electron_TLV_dummy;
               Electron_TLV_dummy.SetPtEtaPhiE(Electron_MC_PT->at(i), Electron_MC_Eta->at(i), Electron_MC_Phi->at(i),Electron_MC_E->at(i));
               vect_Electron_MC.push_back(Electron_TLV_dummy);
               vect_Electron_RECO.push_back(Smear_Electron(Electron_TLV_dummy));
           }
           for(int i=0; i<Muon_MC_Size; i++)
           {
               TLorentzVector Muon_TLV_dummy;
               Muon_TLV_dummy.SetPtEtaPhiE(Muon_MC_PT->at(i),Muon_MC_Eta->at(i), Muon_MC_Phi->at(i), Muon_MC_E->at(i));
               vect_Muon_MC.push_back(Muon_TLV_dummy);
               vect_Muon_RECO.push_back(Smear_Muon(Muon_TLV_dummy));
           }
           
           //Setup Primary Vertex TLV (PV)
           TLorentzVector PV_MC;
           //Smear PV
           TLorentzVector PV_RECO = Smear_PV(PV_MC);
           
           
           //setup and smear Secondary Vertices (SV)
           TLorentzVector dummy_SV = Create_SV(*n2A_MC);
           vect_SV_MC.push_back(dummy_SV);
           vect_SV_RECO.push_back(Smear_SV(dummy_SV));
           dummy_SV = Create_SV(*n2B_MC);
           vect_SV_MC.push_back(dummy_SV);
           vect_SV_RECO.push_back(Smear_SV(dummy_SV));
           
           //Break up MET into parallel and perpendicular components
           TVector3 PT_true = PT_true_TV3(n1A_MC,n1B_MC);
           TVector3 MET_Delphes=*MET;
           MET_Delphes.SetZ(0.0);

           MET_Perpendicular_Mag=PT_true.Cross(zhat).Unit().Dot(MET_Delphes);
           MET_Parallel_Mag=(MET_Delphes.Dot(PT_true))/PT_true.Mag();
           
           hist_MET_Parallel->Fill(MET_Parallel_Mag-PT_true.Mag());
           hist_MET_Perpendicular->Fill(MET_Perpendicular_Mag);
           vect_MET_Parallel.push_back(MET_Parallel_Mag-PT_true.Mag());
           vect_MET_Perpendicular.push_back(MET_Perpendicular_Mag);
           
           for(int i=0; i<Jet_MC_Size; i++)
           {
               pythia_JET_PT+=Jet_MC_PT->at(i);
           }
           
           hist_PT_true->Fill(pythia_JET_PT);
           vect_SumPt.push_back(pythia_JET_PT);
           TVector3 MET_RECO = Smear_MET(PT_true, pythia_JET_PT);
           hist_MET_RECO->Fill(MET_RECO.Mag());
           hist_MET_Delphes->Fill(MET_Delphes.Mag());
           
           for(int i=0; i<HT_Size; i++)
           {
               delphes_JET_PT=HT;
           }
           
           PT_Difference=(pythia_JET_PT-delphes_JET_PT);
           hist_pythia_Pt_vs_delphes_Pt_n2n2j->Fill(pythia_JET_PT,delphes_JET_PT);
           hist_PT_Difference->Fill(PT_Difference);
           //hist_PT_Difference->Fill(abs(PT_Difference));
           
           //reset fill values
           MET_Parallel_Mag=0.0;
           MET_Perpendicular_Mag=0.0;
           
           pythia_JET_PT=0.0;
           delphes_JET_PT=0.0;
           
       }
        
        int Ninterval = 21;
        vector<double> vBinEdge;
        vBinEdge.push_back(0.);
        vBinEdge.push_back(25.);
        vBinEdge.push_back(50.);
        vBinEdge.push_back(75.);
        vBinEdge.push_back(100.);
        vBinEdge.push_back(125.);
        vBinEdge.push_back(150.);
        vBinEdge.push_back(175.);
        vBinEdge.push_back(200.);
        vBinEdge.push_back(225.);
        vBinEdge.push_back(250.);
        vBinEdge.push_back(275.);
        vBinEdge.push_back(300.);
        vBinEdge.push_back(325.);
        vBinEdge.push_back(350.);
        vBinEdge.push_back(375.);
        vBinEdge.push_back(400.);
        vBinEdge.push_back(450.);
        vBinEdge.push_back(500.);
        vBinEdge.push_back(600.);
        vBinEdge.push_back(800.);
        
        vector<double> vBinSum;
        vector<double> vBinSum2;
        vector<double> vBinN;
        
        for(int i=0; i<Ninterval; i++)
        {
            vBinSum.push_back(0.);
            vBinSum2.push_back(0.);
            vBinN.push_back(0.);
            
            char* canvasname_Parallel = new char[100];
            std::sprintf(canvasname_Parallel,"canvas_MET_Parallel_%d",i);
            string canvasname_Parallel_str = "canvas_MET_Parallel"+std::to_string(i);
            TCanvas* dummy_canvas_Parallel = new TCanvas((canvasname_Parallel_str).c_str(),(canvasname_Parallel_str).c_str(),750,500);
            vect_canvas_MET_Parallel.push_back(dummy_canvas_Parallel);
            
            char* canvasname_Perpendicular = new char[100];
            std::sprintf(canvasname_Perpendicular,"canvas_MET_Perpendicular_%d",i);
            string canvasname_Perpendicular_str = "canvas_MET_Perpendicular"+std::to_string(i);
            TCanvas* dummy_canvas_Perpendicular = new TCanvas((canvasname_Perpendicular_str).c_str(),(canvasname_Perpendicular_str).c_str(),750,500);
            vect_canvas_MET_Perpendicular.push_back(dummy_canvas_Perpendicular);
            
            char* histname_Parallel = new char[100];
            std::sprintf(histname_Parallel,"hist_MET_Parallel_%d",i);
            string histname_Parallel_str = "hist_MET_Parallel"+std::to_string(i);
            
            TH1D* dummy_hist_Parallel = new TH1D((histname_Parallel_str).c_str(),(histname_Parallel_str).c_str(),200,-200.0,200.0);
            
            vect_hist_MET_Parallel.push_back(dummy_hist_Parallel);
            
            char* histname_Perpendicular = new char[100];
            std::sprintf(histname_Perpendicular,"hist_MET_Perpendicular_%d",i);
            string histname_Perpendicular_str = "hist_MET_Perpendicular"+std::to_string(i);
            
            TH1D* dummy_hist_Perpendicular = new TH1D((histname_Perpendicular_str).c_str(),(histname_Perpendicular_str).c_str(),200,-200.0,200.0);
            
            vect_hist_MET_Perpendicular.push_back(dummy_hist_Perpendicular);
        }

        for(int k=0; k<vect_SumPt.size(); k++)
        {
            
            for(int i = Ninterval-1; i >= 0; i--){
                if(vect_SumPt[k] > vBinEdge[i]){
                    vBinSum[i] += vect_SumPt[k];
                    vBinSum2[i] += vect_SumPt[k]*vect_SumPt[k];
                    vBinN[i]++;
                    vect_hist_MET_Parallel[i]->Fill(vect_MET_Parallel[k]);
                    vect_hist_MET_Perpendicular[i]->Fill(vect_MET_Perpendicular[k]);
                    break;
                }
            }
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
        hist_MET_Parallel->GetXaxis()->SetTitle("MET Parallel [GeV]");
        hist_MET_Parallel->GetYaxis()->SetTitle("Number of Events");
        hist_MET_Parallel->SetLineColor(kBlack);
        hist_MET_Parallel->Draw();
        c3->SaveAs("MET_Parallel.png");
        c4->cd();
        hist_MET_Perpendicular->GetXaxis()->SetTitle("MET Perpendicular [GeV]");
        hist_MET_Perpendicular->GetYaxis()->SetTitle("Number of Events");
        hist_MET_Perpendicular->SetLineColor(kBlack);
        hist_MET_Perpendicular->Draw();
        c4->SaveAs("MET_Perpendicular.png");
        c5->cd();
        hist_PT_true->GetXaxis()->SetTitle("Sum PT (MC) [GeV]");
        hist_PT_true->GetYaxis()->SetTitle("Number of Events");
        hist_PT_true->SetLineColor(kBlack);
        hist_PT_true->Draw();
        c5->SaveAs("Sum_PT_MC.png");
        
        c7->cd();
        hist_MET_RECO->GetXaxis()->SetTitle("MET Magnitude [GeV]");
        hist_MET_RECO->GetYaxis()->SetTitle("Number of Events");
        hist_MET_RECO->SetLineColor(kBlack);
        hist_MET_Delphes->SetLineColor(kRed);
        hist_MET_RECO->Draw();
        /*TPaveStats* stats_delphes=(TPaveStats*)(hist_MET_Delphes->FindObject("stats"));
        stats_delphes->SetY1NDC(0.5);
        stats_delphes->SetY2NDC(0.7);
        stats_delphes->SetTextColor(kRed);*/
        hist_MET_Delphes->Draw("SAMES");
        c7->SaveAs("MET_RECO_Delphes.png");
        
        int Parallel_entries=0;
        int Perpendicular_entries=0;
        double sigmaMETperp[Ninterval];
        double sigmaMETpar[Ninterval];
        double sigmaMETperp_err[Ninterval];
        double sigmaMETpar_err[Ninterval];
        double sumPT[Ninterval];
        double sumPT_err[Ninterval];
        
        for(int i=0; i < Ninterval; i++)
        {
            sumPT[i] = vBinSum[i]/vBinN[i];
            sumPT_err[i] = sqrt( vBinSum2[i]/vBinN[i] - vBinSum[i]*vBinSum[i]/vBinN[i]/vBinN[i])/sqrt(vBinN[i]);
            
            vect_canvas_MET_Parallel.at(i)->cd();
            vect_hist_MET_Parallel.at(i)->Draw();
            cout << endl;
            //cout << "Ave Sum Pt of " << i << " " << arr_Ave_Sum_Pt[i] << endl;
            //cout << "Fit Result of " << vect_hist_MET_Parallel.at(i)->GetName() << endl;
            vect_hist_MET_Parallel.at(i)->Fit("gaus","Q");
            vect_hist_MET_Parallel.at(i)->GetXaxis()->SetName("Magnitude of the Parallel component of MET [GeV]");
            vect_hist_MET_Parallel.at(i)->GetYaxis()->SetName("Number of Entries");
            string MET_Parallel_Name(vect_hist_MET_Parallel.at(i)->GetName());
            MET_Parallel_Name=MET_Parallel_Name+".png";
            vect_canvas_MET_Parallel.at(i)->SaveAs(MET_Parallel_Name.c_str());
            
            TF1* func = vect_hist_MET_Parallel[i]->GetFunction("gaus");
            sigmaMETpar[i] = func->GetParameter(2)/sumPT[i];
            sigmaMETpar_err[i] = func->GetParError(2)/sumPT[i];
            
            vect_canvas_MET_Perpendicular.at(i)->cd();
            vect_hist_MET_Perpendicular.at(i)->Draw();
            cout << endl;
            //cout << "Fit Result of " << vect_hist_MET_Perpendicular.at(i)->GetName() << endl;
            vect_hist_MET_Perpendicular.at(i)->Fit("gaus","Q");
            vect_hist_MET_Perpendicular.at(i)->GetXaxis()->SetName("Magnitude of the Perpendicular component of MET [GeV]");
            vect_hist_MET_Perpendicular.at(i)->GetYaxis()->SetName("Number of Entries");
            string MET_Perpendicular_Name(vect_hist_MET_Perpendicular.at(i)->GetName());
            MET_Perpendicular_Name=MET_Perpendicular_Name+".png";
            vect_canvas_MET_Perpendicular.at(i)->SaveAs(MET_Perpendicular_Name.c_str());
            
            func = (TF1*) vect_hist_MET_Perpendicular[i]->GetFunction("gaus");
            sigmaMETperp[i] = func->GetParameter(2)/sumPT[i];
            sigmaMETperp_err[i] = func->GetParError(2)/sumPT[i];
        
        }
        
        TGraphErrors *gr_perp = new TGraphErrors(Ninterval, sumPT, sigmaMETperp, sumPT_err, sigmaMETperp_err);
        
        TGraphErrors *gr_par = new TGraphErrors(Ninterval, sumPT, sigmaMETpar, sumPT_err, sigmaMETpar_err);
        
        //Uncomment to draw TGraphs
        Plot_Res_vs_SumPT(gr_perp,gr_par);
    }
