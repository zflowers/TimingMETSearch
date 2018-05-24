    #define Analyze_n2n2j_cxx
    #include "Analyze_n2n2j.h"
    #include <TH2.h>
    #include <TStyle.h>
    #include <TCanvas.h>
    #include <stdio.h>
    #include <string>

// smeared MET given true MET and sumPT
TVector3 Smear_MET(const TVector3& PT_true, double sum_PT);

//take two tgrapherrors, plot, them and fit them, save output to a png
void Plot_Res_vs_SumPT(TGraphErrors* gr_perp, TGraphErrors* gr_par);

Double_t Fit_sqrt(Double_t *x,Double_t *par)
{
    double con0=par[0];
    double con1=par[1];
    double con2=par[2];
    double fitnum=0.0;
    /*if(fitnum < 1e-12 || fitnum == nan)
     {
     cout << "FIT FAIL" << endl;
     return 0.;
     }*/
    fitnum=sqrt((con0*con0)/(x[0]*x[0])+(con1*con1)/x[0]+con2*con2);
    return fitnum;
}

TGraphErrors* TGE_TF1(TGraphErrors *plot, TF1 *fit_func)
{
    int N = plot->GetN();
    double xnew[N];
    double ynew[N];
    double xnew_err[N];
    double ynew_err[N];
    for(int i=0; i < N; i++)
    {
        double x,y,ex,ey;
        plot->GetPoint(i, x, y);
        xnew[i] = x;
        ynew[i] = fit_func->Eval(x)-y;
    }
    TGraphErrors* res_plot = new TGraphErrors(N,xnew,ynew);
    for(int i=0; i < N; i++)
    {
        res_plot->SetPointError(i,plot->GetErrorX(i),plot->GetErrorY(i));
    }
    res_plot->SetTitle("");
    return res_plot;
}

    void Analyze_n2n2j::Loop()
    {
        gRandom = new TRandom3();
        //gStyle->SetOptStat(000002000);
        //gStyle->SetOptStat(1);
        gStyle->SetOptFit(1);
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
        TCanvas *c5 = new TCanvas("c5","Canvas for Sum PT (MC)",750,500);
        TCanvas *c7 = new TCanvas("c7","Canvas for MET RECO",750,500);
        
        TH2F *hist_pythia_Pt_vs_delphes_Pt_n2n2j = new TH2F("hist_pythia_Pt_vs_delphes_Pt_n2n2j","Pythia Jets Pt vs Delphes Jets Pt", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
        
        TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,-700.0,200.0);
        //TH1D *hist_PT_Difference = new TH1D("hist_PT_Difference","Difference between Pythia and Delphes PT",100,0.0,700.0);
        
        TH1D *hist_MET_Parallel = new TH1D("hist_MET_Parallel","MET_Parallel",100,-100.0,100.0);
        TH1D *hist_MET_Perpendicular = new TH1D("hist_MET_Perpendicular","MET_Perpendicular",100,-100.0,100.0);
        TH1D *hist_PT_true = new TH1D("hist_PT_true","Histogram of Sum PT MC",100,0.0,1200.0);
        
        TH1D *hist_MET_RECO = new TH1D("hist_MET_RECO","Smeared PT True (MET)",100,0.0,800.0);
        
        std::vector<TH1D*> vect_hist_MET_Parallel;
        std::vector<TH1D*> vect_hist_MET_Perpendicular;
        std::vector<double> vect_SumPt;
        std::vector<TCanvas*> vect_canvas_MET_Parallel;
        std::vector<TCanvas*> vect_canvas_MET_Perpendicular;
        std::vector<double> vect_MET_Parallel;
        std::vector<double> vect_MET_Perpendicular;
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
           
           //Break up MET into parallel and perpendicular components
           TVector3 PT_true;
           TVector3 n1A_vect=n1A_MC->Vect();
           TVector3 n1B_vect=n1B_MC->Vect();
           PT_true=n1A_vect+n1B_vect;
           TVector3 MET_Delphes=*MET;
           PT_true.SetZ(0.0);
           MET_Delphes.SetZ(0.0);

           MET_Perpendicular_Mag=PT_true.Cross(BeamSpot).Unit().Dot(MET_Delphes);
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
        hist_MET_RECO->GetXaxis()->SetTitle("MET RECO [GeV]");
        hist_MET_RECO->GetYaxis()->SetTitle("Number of Events");
        hist_MET_RECO->SetLineColor(kBlack);
        hist_MET_RECO->Draw();
        c7->SaveAs("MET_RECO.png");
        
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
        
        //push this into one function
        Plot_Res_vs_SumPT(gr_perp,gr_par);
        
    }


void Plot_Res_vs_SumPT(TGraphErrors* gr_perp, TGraphErrors* gr_par)
{
    TCanvas *c6 = new TCanvas("c6","Canvas for TGE Perp",750,500);

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
    pad1->SetBottomMargin(0.95);
    //pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    gr_perp->SetMarkerColor(kRed);
    gr_perp->SetLineColor(kRed);
    gr_par->SetMarkerColor(kBlue);
    gr_par->SetLineColor(kBlue);


    TMultiGraph* mg = new TMultiGraph();


    double con0_perp=0.2;
    double con1_perp=0.02;
    double con2_perp=0.1;
    double con0_par=0.2;
    double con1_par=0.05;
    double con2_par=0.1;

    cout << endl << endl << "Fitting Perpendicular TGE" << endl;

    TF1 *func_perp = new TF1("func_perp", Fit_sqrt,0.0, 1200.0, 3);
    func_perp->SetParameter(0,con0_perp);
    func_perp->SetParameter(1,con1_perp);
    func_perp->SetParameter(2,con2_perp);
    func_perp->SetParName(0,"con0_perp");
    func_perp->SetParName(1,"con1_perp");
    func_perp->SetParName(2,"con2_perp");
    func_perp->SetParLimits(0,0.0,10.0);
    func_perp->SetParLimits(1,0.0,1.0);
    func_perp->SetParLimits(2,0.0,1.0);
    func_perp->SetNpx(10000000);
    func_perp->SetLineColor(kYellow);
    gr_perp->Fit(func_perp,"EM");
    mg->Add(gr_perp);

    cout << "Fitting Parallel TGE" << endl;

    TF1 *func_par = new TF1("func_par", Fit_sqrt,0.0, 1200.0, 3);
    func_par->SetParameter(0,con0_par);
    func_par->SetParameter(1,con1_par);
    func_par->SetParameter(2,con2_par);
    func_par->SetParName(0,"con0_par");
    func_par->SetParName(1,"con1_par");
    func_par->SetParName(2,"con2_par");
    func_par->SetParLimits(0,0.0,10.0);
    func_par->SetParLimits(1,0.0,1.0);
    func_par->SetParLimits(2,0.0,1.0);
    func_par->SetNpx(10000000);
    func_par->SetLineColor(kGreen);
    gr_par->Fit(func_par,"EM");
    mg->Add(gr_par);

    pad1->Update();

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("SumPt [GeV]");
    mg->GetYaxis()->SetTitle("Resolution/SumPt [GeV]");
    mg->GetXaxis()->SetTitleSize(0.047);
    mg->GetYaxis()->SetTitleSize(0.047);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05);
    pad1->Update();

    TPaveStats* stats_par=(TPaveStats*)(gr_par->GetListOfFunctions()->FindObject("stats"));
    stats_par->SetY1NDC(0.5);
    stats_par->SetY2NDC(0.7);
    stats_par->SetTextColor(kBlue);
    TPaveStats* stats_perp=(TPaveStats*)(gr_perp->GetListOfFunctions()->FindObject("stats"));
    stats_perp->SetTextColor(kRed);

    pad1->Modified();

    c6->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.05);
    pad2->Update();
    c6->Update();

    TGraphErrors* residual_perp=TGE_TF1(gr_perp,func_perp);
    residual_perp->SetMarkerColor(kRed);
    residual_perp->SetLineColor(kRed);

    TGraphErrors* residual_par=TGE_TF1(gr_par,func_par);
    residual_par->SetMarkerColor(kBlue);
    residual_par->SetLineColor(kBlue);
    TMultiGraph* mg_res = new TMultiGraph();
    mg_res->Add(residual_perp);
    mg_res->Add(residual_par);
    mg_res->Draw("AP");
    mg_res->GetXaxis()->SetLabelSize(0.12);
    mg_res->GetYaxis()->SetLabelSize(0.12);

    TLine* line = new TLine(0.0,0.0,1048.0,0.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw();
    pad2->Update();
    c6->Update();

    c6->SaveAs("Resolution_vs_SumPT.png");
}



TVector3 Smear_MET(const TVector3& PT_true, double sum_PT)
{
    TVector3 zhat;
    zhat.SetXYZ(0.0,0.0,1.0);
  
    double sigma_perp = sum_PT*sqrt((4.20687*4.20687)/(sum_PT*sum_PT)+(0.548129*0.548129)/sum_PT+0.0116634*0.0116634);
    double sigma_par  = sum_PT*sqrt((4.02779*4.02779)/(sum_PT*sum_PT)+(0.512707*0.512707)/sum_PT+0.0431817*0.0431817);
    
    TVector3 parhat  = PT_true.Unit();
    TVector3 perphat = PT_true.Cross(zhat).Unit();
    
    TVector3 MET = PT_true + gRandom->Gaus(0.,sigma_par)*parhat + gRandom->Gaus(0.,sigma_perp)*perphat;
    
    return MET;
    
}
