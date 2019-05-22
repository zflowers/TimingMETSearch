#include <TGraphErrors.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMinuit.h>

double One_Sigma_Interval(std::vector<double> Sigma_Var) //Pass a vector of analytic calculations of some mass
{
    std::sort(Sigma_Var.begin(),Sigma_Var.end()); //sort the vector in ascending order
    std::vector<double> Interval; //For keeping the widths of the intervals
    int N = Sigma_Var.size();
    int N68 = N*0.68;
    for(int i = 0; i < (N - N68); i++) //loop over the given vector and stop at the last entry
    {
        Interval.push_back(Sigma_Var.at(N68+i)-Sigma_Var.at(i)); //push back the width of an interval
    }
    std::sort(Interval.begin(),Interval.end());
    return Interval.at(0); //return the smallest interval
}

void CreatePalette()
{
    //From Prof. Chris Rogan
    //Useful ROOT5 palette
    Int_t MyPalette[28];
    Int_t alpha=1;
    
    Double_t Stops[5] = { 0.00, 0.50, 0.70, 0.82, 1.00 };
    Double_t Red[5]   = { 0.00, 0.00, 0.74, 1.00, 1.00 };
    Double_t Green[5] = { 0.00, 0.61, 0.82, 0.70, 1.00 };
    Double_t Blue[5]  = { 0.31, 0.73, 0.08, 0.00, 1.00 };
    
    Int_t test = TColor::CreateGradientColorTable(5, Stops, Red, Green, Blue, 28);
    for(int i=0; i<28; i++)
    {
        MyPalette[i] = test+i;
    }
    gStyle->SetPalette(28,MyPalette);
}

Double_t Expo_sqrt(Double_t *x,Double_t *par)
{
    double c0=par[0];
    double c1=par[1];
    double fitnum=0.0;
    fitnum = c0*exp(-c1*sqrt(x[0]));
    return fitnum;
}

Double_t Fit_sqrt_three_term(Double_t *x,Double_t *par)
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

Double_t Fit_sqrt_two_term(Double_t *x,Double_t *par)
{
    double con0=par[0];
    double con1=par[1];
    double fitnum=0.0;
    /*if(fitnum < 1e-12 || fitnum == nan)
     {
     cout << "FIT FAIL" << endl;
     return 0.;
     }*/
    fitnum=sqrt((con0*con0)/(x[0]*x[0])+con1*con1);
    return fitnum;
}

Double_t Fit_Voigt(Double_t *x, Double_t *par)
{
    double norm = par[0];
    double mean = par[1];
    double sigma = par[2];
    double lg = par[3];
    
    if ((sigma < 0 || lg < 0) || (sigma==0 && lg==0))
    {
        //cout << "CHECK FIT" << endl;
        //return 10.0;
        //sigma+=1.0e-06;
    }
    
    return norm*TMath::Voigt(x[0]-mean,sigma,lg);
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

void Plot_Res_vs_SumPT(TGraphErrors* gr_perp, TGraphErrors* gr_par, TCanvas* canvas)
{
    canvas->cd();
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
    
    double con0_perp=5.0;
    double con1_perp=0.02;
    double con2_perp=0.01;
    double con0_par=5.0;
    double con1_par=0.05;
    double con2_par=0.01;
    
    cout << endl << endl << "Fitting Perpendicular TGE SumPT" << endl;
    
    TF1 *func_perp = new TF1("func_perp_SumPT", Fit_sqrt_three_term,gr_perp->GetXaxis()->GetXmin(), gr_perp->GetXaxis()->GetXmax(), 3);
    func_perp->SetParameter(0,con0_perp);
    func_perp->SetParameter(1,con1_perp);
    func_perp->SetParameter(2,con2_perp);
    func_perp->SetParName(0,"con0_perp");
    func_perp->SetParName(1,"con1_perp");
    func_perp->SetParName(2,"con2_perp");
    //func_perp->SetParLimits(0,0.0,100.0);
    //func_perp->SetParLimits(1,0.0,100.0);
    //func_perp->SetParLimits(2,0.0,100.0);
    func_perp->SetLineColor(kYellow);
    gr_perp->Fit(func_perp,"EM+");
    func_perp->SetNpx(10000000);
    mg->Add(gr_perp);
    
    cout << endl << "Fitting Parallel TGE SumPT" << endl;
    
    TF1 *func_par = new TF1("func_par_SumPT", Fit_sqrt_three_term,gr_perp->GetXaxis()->GetXmin(), gr_perp->GetXaxis()->GetXmax(), 3);
    func_par->SetParameter(0,con0_par);
    func_par->SetParameter(1,con1_par);
    func_par->SetParameter(2,con2_par);
    func_par->SetParName(0,"con0_par");
    func_par->SetParName(1,"con1_par");
    func_par->SetParName(2,"con2_par");
    //func_par->SetParLimits(0,0.0,100.0);
    //func_par->SetParLimits(1,0.0,100.0);
    //func_par->SetParLimits(2,0.0,100.0);
    func_par->SetLineColor(kGreen);
    gr_par->Fit(func_par,"EM+");
    func_par->SetNpx(10000000);
    mg->Add(gr_par);
    
    pad1->Update();
    
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("Pt of chi2 chi2 system [GeV]");
    mg->GetYaxis()->SetTitle("Resolution/Pt of chi2 chi2 system [GeV]");
    //mg->GetXaxis()->SetTitle("SumPt [GeV]");
    //mg->GetYaxis()->SetTitle("Resolution/SumPt [GeV]");
    mg->GetXaxis()->SetTitleSize(0.047);
    mg->GetYaxis()->SetTitleSize(0.047);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05);
    pad1->Update();

    TLegend* leg = new TLegend(0.78,0.3,0.98,.5, "");
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetFillColor(kWhite);
    leg->SetTextColor(kBlack);
    leg->AddEntry(func_par,"Par Fit","L");
    leg->AddEntry(func_perp,"Perp Fit","L");
    leg->Draw("SAMES");
    pad1->Update();
    
    func_perp->Draw("SAMES");
    func_par->Draw("SAMES");
    pad1->Update();
    canvas->Update();
    
    TPaveStats* stats_par=(TPaveStats*)(gr_par->GetListOfFunctions()->FindObject("stats"));
    if(stats_par == NULL)
    {
        cout << "Error: Stats Boxes Turned Off" << endl;
    }
    stats_par->SetY1NDC(0.7);
    stats_par->SetY2NDC(0.9);
    stats_par->SetTextColor(kBlue);
    
    TPaveStats* stats_perp=(TPaveStats*)(gr_perp->GetListOfFunctions()->FindObject("stats"));
    stats_perp->SetY1NDC(0.5);
    stats_perp->SetY2NDC(0.7);
    stats_perp->SetTextColor(kRed);

    pad1->Modified();
    canvas->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->Draw();
    pad2->cd();
    pad2->SetTopMargin(0.05);
    pad2->Update();
    canvas->Update();
    
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
    
    TLine* line = new TLine(residual_par->GetXaxis()->GetXmin(),0.0,residual_par->GetXaxis()->GetXmax(),0.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw();
    pad2->Update();
    canvas->Update();
}

void Plot_Res_vs_PT_n2n2(TGraphErrors* gr_perp, TGraphErrors* gr_par, TCanvas* canvas)
{
    canvas->cd();
    TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1.0);
    pad3->SetBottomMargin(0.95);
    //pad3->SetLogy();
    pad3->Draw();
    pad3->cd();
    
    gr_perp->SetMarkerColor(kRed);
    gr_perp->SetLineColor(kRed);
    gr_par->SetMarkerColor(kBlue);
    gr_par->SetLineColor(kBlue);
    
    
    TMultiGraph* mg = new TMultiGraph();
    
    double con0_perp=5.0;
    double con1_perp=0.02;
    double con2_perp=0.01;
    double con0_par=5.0;
    double con1_par=0.05;
    double con2_par=0.01;
    
    cout << endl << endl << "Fitting Perpendicular TGE PT_n2n2" << endl;
    
    TF1 *func_perp = new TF1("func_perp_PT_n2n2", Fit_sqrt_three_term,gr_perp->GetXaxis()->GetXmin(), gr_perp->GetXaxis()->GetXmax(), 3);
    func_perp->SetParameter(0,con0_perp);
    func_perp->SetParameter(1,con1_perp);
    func_perp->SetParameter(2,con2_perp);
    func_perp->SetParName(0,"con0_perp");
    func_perp->SetParName(1,"con1_perp");
    func_perp->SetParName(2,"con2_perp");
    //func_perp->SetParLimits(0,0.0,100.0);
    //func_perp->SetParLimits(1,0.0,100.0);
    //func_perp->SetParLimits(2,0.0,100.0);
    func_perp->SetLineColor(kYellow);
    gr_perp->Fit(func_perp,"EM+");
    func_perp->SetNpx(10000000);
    mg->Add(gr_perp);
    
    cout << endl << "Fitting Parallel TGE PT_n2n2" << endl;
    
    TF1 *func_par = new TF1("func_par_PT_n2n2", Fit_sqrt_three_term,gr_perp->GetXaxis()->GetXmin(), gr_perp->GetXaxis()->GetXmax(), 3);
    func_par->SetParameter(0,con0_par);
    func_par->SetParameter(1,con1_par);
    func_par->SetParameter(2,con2_par);
    func_par->SetParName(0,"con0_par");
    func_par->SetParName(1,"con1_par");
    func_par->SetParName(2,"con2_par");
    //func_par->SetParLimits(0,0.0,100.0);
    //func_par->SetParLimits(1,0.0,100.0);
    //func_par->SetParLimits(2,0.0,100.0);
    func_par->SetLineColor(kGreen);
    gr_par->Fit(func_par,"EM+");
    func_par->SetNpx(10000000);
    mg->Add(gr_par);
    
    pad3->Update();
    
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("Pt of chi2 chi2 system [GeV]");
    mg->GetYaxis()->SetTitle("Resolution/Pt of chi2 chi2 system [GeV]");
    mg->GetXaxis()->SetTitleSize(0.047);
    mg->GetYaxis()->SetTitleSize(0.047);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05);
    pad3->Update();
    
    TLegend* leg = new TLegend(0.78,0.3,0.98,.5, "");
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetFillColor(kWhite);
    leg->SetTextColor(kBlack);
    leg->AddEntry(func_par,"Par Fit","L");
    leg->AddEntry(func_perp,"Perp Fit","L");
    leg->Draw("SAMES");
    pad3->Update();
    
    func_perp->Draw("SAMES");
    func_par->Draw("SAMES");
    pad3->Update();
    canvas->Update();
    
    TPaveStats* stats_par=(TPaveStats*)(gr_par->GetListOfFunctions()->FindObject("stats"));
    //TF1* func = (TF1*)(gr_errs->GetListOfFunctions()->FindObject("nameoffunc"));
    if(stats_par == NULL)
    {
        cout << "Error: Stats Boxes Turned Off" << endl;
    }
    stats_par->SetY1NDC(0.7);
    stats_par->SetY2NDC(0.9);
    stats_par->SetTextColor(kBlue);
    
    TPaveStats* stats_perp=(TPaveStats*)(gr_perp->GetListOfFunctions()->FindObject("stats"));
    stats_perp->SetY1NDC(0.5);
    stats_perp->SetY2NDC(0.7);
    stats_perp->SetTextColor(kRed);
    
    pad3->Modified();
    canvas->cd();
    TPad *pad4 = new TPad("pad4", "pad4", 0, 0.05, 1, 0.3);
    pad4->Draw();
    pad4->cd();
    pad4->SetTopMargin(0.05);
    pad4->Update();
    canvas->Update();
    
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
    
    TLine* line = new TLine(residual_par->GetXaxis()->GetXmin(),0.0,residual_par->GetXaxis()->GetXmax(),0.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw();
    pad4->Update();
    canvas->Update();
}

Double_t Fit_TGE_func(Double_t *x,Double_t *par)
{
    double c0=par[0];
    double c1=par[1];
    double c2=par[2];
    double fitnum=0.0;
    /*if(fitnum < 1e-12 || fitnum == nan)
     {
     cout << "FIT FAIL" << endl;
     return 0.;
     }*/
    fitnum=c0*x[0]*x[0]+c1*x[0]+c2;
    return fitnum;
}

TH1F* res_TH1_TF1(TH1* hist, TF1* fit_func)
{
    int Bins=hist->GetNbinsX();
    TH1F* res_hist = new TH1F("res_hist","",Bins,hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    for(int i=0; i < Bins; i++)
    {
        //res_hist->SetBinContent(i+2, (hist->GetBinError(i+2) > 0. ? 1.0/((fit_func->Eval(hist->GetBinCenter(i+2)))-hist->GetBinContent(i+2))/hist->GetBinError(i+2) : 0.) );
        /*double val=(hist->GetBinError(i+2) > 0. ? ((fit_func->Eval(hist->GetBinCenter(i+2)))-hist->GetBinContent(i+2))/hist->GetBinError(i+2) : 0.);
         if(val != 0.0)
         {
         res_hist->SetBinContent(i+2, 1.0/val);
         }
         else
         {
         res_hist->SetBinContent(i+2, 0.0);
         }*/
        //res_hist->SetBinContent(i+1, (fit_func->Eval(hist->GetBinCenter(i+1)))-hist->GetBinCenter(i+1));
        res_hist->SetBinContent(i+1, 1.0-((fit_func->Eval(hist->GetBinCenter(i+1)))-hist->GetBinContent(i+1))/hist->GetBinError(i+1)); //subtract
        //res_hist->SetBinContent(i+1, ((fit_func->Eval(hist->GetBinCenter(i+1)))-hist->GetBinContent(i+1))/hist->GetBinError(i+1));
        //res_hist->SetBinContent(i+1, ((fit_func->Eval(hist->GetBinCenter(i+1)))/hist->GetBinContent(i+1)));///hist->GetBinError(i+1)); //divide
    }
    res_hist->SetBinContent(1,0.0);
    return res_hist;
}

void Fit_TGE(TGraphErrors* gr, TCanvas *canvas)
{
    cout << "Fitting " << gr->GetName() << endl;
    canvas->cd();
    
    TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1.0);
    pad3->SetBottomMargin(0.95);
    //pad3->SetLogy();
    pad3->Draw();
    pad3->cd();

    double c0=-0.0004;
    double c1=0.001;
    double c2=2.0;

    TF1 *func = new TF1("func", Fit_TGE_func,gr->GetXaxis()->GetXmin(), gr->GetXaxis()->GetXmax(), 3);
    func->SetParameter(0,c0);
    func->SetParameter(1,c1);
    func->SetParameter(2,c2);
    func->SetParName(0,"c0");
    func->SetParName(1,"c1");
    func->SetParName(2,"c2");
    //func->SetParLimits(0,0.0,10.0);
    //func->SetParLimits(1,0.0,100000.0);
    //func->SetParLimits(2,0.0,10.0);
 
 
    gr->Fit(func,"EMR");
    func->SetNpx(10000000);
    pad3->Update();
    //pad3->Modified();
    
    gr->Draw("AP");
    gr->SetTitle("");
    gr->GetYaxis()->SetTitle("Value of P1");
    gr->GetXaxis()->SetTitle("PT of n2n2");
    gr->GetXaxis()->SetTitleSize(0.047);
    gr->GetYaxis()->SetTitleSize(0.047);
    gr->GetXaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetLabelSize(0.05);
    pad3->Update();
    
    
    canvas->cd();
    TPad *pad4 = new TPad("pad4", "pad4", 0, 0.05, 1, 0.3);
    pad4->Draw();
    pad4->cd();
    pad4->SetTopMargin(0.05);
    pad4->Update();
    canvas->Update();
    
    TGraphErrors* residual=TGE_TF1(gr,func);
    residual->Draw("AP");
    residual->GetXaxis()->SetLabelSize(0.12);
    residual->GetYaxis()->SetLabelSize(0.12);
    
    TLine* line = new TLine(gr->GetXaxis()->GetXmin(),0.0,gr->GetXaxis()->GetXmax(),0.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw();
    pad4->Update();
    canvas->Update();
    
    canvas->SaveAs("FitParametersGraph.png");
}

void Draw_Two_Hists(TH1* hist1, TH1* hist2, TCanvas* canvas)
{
    canvas->cd();
    hist1->Draw();
    hist2->Draw("SAMES");
    canvas->Update();
    
    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);
    
    TPaveStats *stats1 = (TPaveStats*)hist1->FindObject("stats");
    string stats1_name = hist1->GetName();
    stats1_name+="_stats";
    stats1->SetName(stats1_name.c_str());
    stats1->SetLineColor(kRed);
    stats1->SetTextColor(kRed);
    
    TPaveStats *stats2 = (TPaveStats*)hist2->FindObject("stats");
    string stats2_name = hist2->GetName();
    stats2_name+="_stats";
    stats2->SetName(stats2_name.c_str());
    stats2->SetY1NDC(.7);
    stats2->SetY2NDC(.5);
    stats2->SetLineColor(kBlue);
    stats2->SetTextColor(kBlue);
    canvas->Update();
}

void Draw_Hists(vector<TH1F*> hists, TCanvas* canvas)
{
    canvas->cd();
    int size = hists.size();
    for(int i=0; i<size; i++)
    {
        TRandom3 rndm_color;
        rndm_color.SetSeed(0);
        Double_t dbl_color_index = rndm_color.Uniform(300.0,1000.0);
        int color_index = (int)dbl_color_index;
        hists.at(i)->Draw("SAMES");
        canvas->Update();
        hists.at(i)->SetLineColor(color_index);
        
        TPaveStats *stats = (TPaveStats*)hists.at(i)->FindObject("stats");
        string stats_name = hists.at(i)->GetName();
        stats_name+="_stats";
        stats->SetName(stats_name.c_str());
        stats->SetTextColor(color_index);
        stats->SetLineColor(color_index);
        canvas->Update();
    }
}

void Fit_Hist_Two_Funcs(TH1* hist, TF1* func1, TF1* func2, TCanvas* canvas)
{
    hist->SetLineColor(kBlack);
    TH1* hist_clone = (TH1*)(hist->Clone());
    func1->SetLineColor(kRed);
    func2->SetLineColor(kBlue);
    
    hist->Fit(func1,"EMR");
    TString status_func1 = gMinuit->fCstatu;
    func1->Draw("SAMES");
    canvas->Update();
    
    TPaveStats *stats1 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    stats1->SetTextColor(kRed);
    stats1->SetLineColor(kRed);
    stats1->SetY1NDC(0.9);
    stats1->SetY2NDC(0.7);
    stats1->Draw("SAMES");
    canvas->Update();
    
    hist_clone->Fit(func2,"EMR","SAMES");
    canvas->Update();
    
    TPaveStats *stats2 = (TPaveStats*)hist_clone->GetListOfFunctions()->FindObject("stats");
    stats2->SetTextColor(kBlue);
    stats2->SetLineColor(kBlue);
    stats2->SetY1NDC(0.7);
    stats2->SetY2NDC(0.5);
    canvas->Update();
    
    TString status_func2 = gMinuit->fCstatu;
    
    TString status_func1_leg = "Func1 Status = ";
    status_func1_leg+=status_func1;
    TString status_func2_leg = "Func2 Status = ";
    status_func2_leg+=status_func2;
    
    TLegend* leg1 = new TLegend(0.55,0.4,0.98,0.5,"");
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(kWhite);
    leg1->SetTextColor(kRed);
    leg1->AddEntry((TObject*)0, status_func1_leg, "");
    leg1->Draw("SAMES");
    
    canvas->Update();
    
    TLegend* leg2 = new TLegend(0.55,0.3,0.98,0.4,"");
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04);
    leg2->SetFillColor(kWhite);
    leg2->SetTextColor(kBlue);
    leg2->AddEntry((TObject*)0, status_func2_leg, "");
    leg2->Draw("SAMES");
    
    canvas->Update();
}

