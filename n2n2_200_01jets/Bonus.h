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
    
    cout << endl << "Fitting Parallel TGE" << endl;
    
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
