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
#include <TLatex.h>
#include <TFile.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMinuit.h>
#include <numeric>
#include <algorithm>

double Get_Mode(TH1D& hist)
{
    return hist.GetXaxis()->GetBinCenter(hist.GetMaximumBin());
}

Double_t Fit_Mode(Double_t *x,Double_t *par)
{
    return par[0]*(1./x[0])+par[1];
}
/*
double Hist_Mode(TH1D hist)
{
    return hist.GetXaxis()->GetBinCenter(hist.GetMaximumBin());
}
*/
double Hist_Mode(TH1D hist)
{
    TGraph* width_bins = new TGraph(0);
    width_bins->SetTitle("Mode");
    int base = 2.;
    for(int i = 0; i < 11; i++)
    {
        if(TMath::Power(base,i) == hist.GetNbinsX())
        {
            width_bins->Set(i);
            break;
        }
    }
    if(width_bins->GetN() == 0)
    {
        cout << "Improper number of Histogram Bins ";
        cout << "Must be a power of: " << base << " ";
        cout << "Returing zero now" << endl;
        return 0.;
    }
    double constant = 0.;
    for(int j = width_bins->GetN(); j >= 0; j--)
    {
        width_bins->SetPoint(j,TMath::Power(base,j),Get_Mode(hist));
        if(j == width_bins->GetN()/2) { constant = Get_Mode(hist); }
        if(hist.GetNbinsX() == 1) break;
        hist.Rebin(base);
    }
    string name = hist.GetName();
    name+="Mode";
    TCanvas* canv = new TCanvas(name.c_str(),"",750,500);
    canv->cd();
    width_bins->Draw("AP");
    TF1* fit_func = new TF1(hist.GetName(),Fit_Mode,width_bins->GetXaxis()->GetXmin(),width_bins->GetXaxis()->GetXmax(),2);
    fit_func->SetParameter(0,Get_Mode(hist));
    fit_func->SetParameter(1,constant);
    width_bins->Fit(fit_func,"QEMR");
    double mode = fit_func->GetParameter(1);
    //canv->SaveAs((name+".pdf").c_str());
    delete fit_func;
    delete width_bins;
    delete canv;
    return mode;
}

double Get_FWHM(TH1D& hist)
{
    int bin1 = hist.FindFirstBinAbove(hist.GetMaximum()/2);
    int bin2 = hist.FindLastBinAbove(hist.GetMaximum()/2);
    return (hist.GetXaxis()->GetBinUpEdge(bin2) - hist.GetXaxis()->GetBinLowEdge(bin1));
}
/*
double Hist_FWHM(TH1D& hist)
{
    int bin1 = hist.FindFirstBinAbove(hist.GetMaximum()/2);
    int bin2 = hist.FindLastBinAbove(hist.GetMaximum()/2);
    return (hist.GetXaxis()->GetBinUpEdge(bin2) - hist.GetXaxis()->GetBinLowEdge(bin1));
}
*/
Double_t Fit_FWHM(Double_t *x,Double_t *par)
{
    return par[0]*(1./x[0])+par[1];
}

double Hist_FWHM(TH1D hist)
{
    TGraph* width_bins = new TGraph(0);
    width_bins->SetTitle("FWHM");
    int base = 2.;
    for(int i = 0; i < 11; i++)
    {
        if(TMath::Power(base,i) == hist.GetNbinsX())
        {
            width_bins->Set(i);
            break;
        }
    }
    if(width_bins->GetN() == 0)
    {
        cout << "Improper number of Histogram Bins ";
        cout << "Must be a power of: " << base << " ";
        cout << "Returing zero now" << endl;
        return 0.;
    }
    double constant = 0.;
    for(int j = width_bins->GetN(); j >= 0; j--)
    {
        width_bins->SetPoint(j,TMath::Power(base,j),Get_FWHM(hist));
        if(j == width_bins->GetN()/2) { constant = Get_FWHM(hist); }
        if(hist.GetNbinsX() == 1) break;
        hist.Rebin(base);
    }
    string name = hist.GetName();
    name+="FWHM";
    TCanvas* canv = new TCanvas(name.c_str(),"",750,500);
    canv->cd();
    width_bins->Draw("AP");
    TF1* fit_func = new TF1(hist.GetName(),Fit_FWHM,width_bins->GetXaxis()->GetXmin(),width_bins->GetXaxis()->GetXmax(),2);
    fit_func->SetParameter(0,Get_FWHM(hist));
    fit_func->SetParameter(1,constant);
    width_bins->Fit(fit_func,"QEMR");
    double fwhm = fit_func->GetParameter(1);
    //canv->SaveAs((name+".pdf").c_str());
    delete fit_func;
    delete width_bins;
    delete canv;
    return fwhm;
}

double Hist_68_Interval(const TH1F& hist_user)
{
    TH1F hist = hist_user;
    hist.Scale(1./hist.GetEntries());
    int mode_Bin = hist.GetBin(hist.GetMaximumBin());
    int left_bin = mode_Bin-1;
    int right_bin = mode_Bin+1;
    double left_content = hist.GetBinContent(left_bin);
    double right_content = hist.GetBinContent(right_bin);
    double prob = hist.GetBinContent(mode_Bin);
    if(prob > 0.68) return hist.GetBinCenter(mode_Bin)/2.;
    else if((prob + right_content) > 0.68) return (hist.GetBinCenter(right_bin)-hist.GetBinCenter(mode_Bin))/2.;
    else if((prob + left_content) > 0.68) return (hist.GetBinCenter(mode_Bin)-hist.GetBinCenter(left_bin))/2.;
    else if((prob + left_content + right_content) > 0.68) return (hist.GetBinCenter(right_bin)-hist.GetBinCenter(left_bin))/2.;
    else
    {
        prob = prob + left_content + right_content;
        int left_it = 1;
        int right_it = 1;
        while(prob < 0.68)
        {
            if((prob + hist.GetBinContent(left_bin-left_it)) > (prob + hist.GetBinContent(right_bin+right_it)))
            {
                prob+=hist.GetBinContent(left_bin-left_it);
                left_it++;
            }
            else
            {
                prob+=hist.GetBinContent(right_bin+right_it);
                right_it++;
            }
        }
        left_bin-=left_it;
        right_bin+=right_it;
        return (hist.GetBinCenter(right_bin)-hist.GetBinCenter(left_bin))/2.;
    }
}

vector<TH1*> list_histos(const char *fname)
{
    std::vector<TH1*> vect_hist;
    std::vector<string> vect_names;
    TKey *key;
    TKey *key2;
    bool skip = true;
    TFile *f = TFile::Open(fname, "READ");
    if(!f || f->IsZombie())
    {
        cout << "Unable to open " << fname << " for reading..." << endl;
        return vect_hist;
    }
    TDirectoryFile* folder = nullptr;
    f->GetObject("Plots",folder);
    TIter next((TList *)folder->GetListOfKeys());
    while((key = (TKey *)next()))
    {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if(cl->InheritsFrom("TCanvas"))
        {
            TCanvas* c = (TCanvas*)key->ReadObj();
            TIter next2(c->GetListOfPrimitives());
            while((key2 = (TKey *)next2())){
                if(key2->InheritsFrom("TH1"))
                {
                    if(skip)
                    {
                        skip = false;
                        continue;
                    }
                    TH1* h = (TH1*)key2;
                    vect_hist.push_back(h);
                }
            }
        }
    }
    return vect_hist;
}

double Vector_Mean(const std::vector<double>& vect) //returns mean of a vector
{
    return accumulate(vect.begin(),vect.end(),0.0)/vect.size();
}

double Vector_Mode(std::vector<double>& vect)
{
    std::sort(vect.begin(),vect.end());
    double binning = 0.5;
    int N = vect.size();
    double mode = 0.;
    int NWidth = N*binning;
    int left_edge = 0;
    int left = left_edge;
    int right_edge = NWidth;
    int right = right_edge;
    double width = fabs((vect.at(right)-vect.at(left)));
    while(right < int(vect.size()))
    {
        if(fabs(((vect.at(right)-vect.at(left)))) < width)
        {
            width = fabs(((vect.at(right)-vect.at(left))));
            left_edge = left;
            right_edge = right;
        }
        left++;
        right++;
    }
    mode = ((vect.at(left_edge)+vect.at(right_edge))/2.);
    return mode;
}

double One_Sigma_Interval(std::vector<double> vect) //Pass a vector of analytic calculations of some mass
{
    std::sort(vect.begin(),vect.end()); //sort the vector in ascending order
    int N = vect.size();
    int NWidth = N*0.68;
    int left_edge = 0;
    int left = left_edge;
    int right_edge = NWidth;
    int right = right_edge;
    double width = fabs((vect.at(right)-vect.at(left)));
    while(right < int(vect.size()))
    {
        if(fabs(((vect.at(right)-vect.at(left)))) < width)
        {
            width = fabs(((vect.at(right)-vect.at(left))));
            left_edge = left;
            right_edge = right;
        }
        left++;
        right++;
    }
    //return ((vect.at(right_edge)-vect.at(left_edge)));
    return ((vect.at(right_edge)-vect.at(left_edge))/2.);
}

void CreatePalette()
{
    //From Prof. Chris Rogan
    //Useful ROOT5 palette
    //Call with CreatePalette()
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

void setMyStyle()
{
    CreatePalette();
    TStyle* myStyle = new TStyle("myStyle","myStyle");
    
    // For the canvas:
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetCanvasColor(kWhite);
    myStyle->SetCanvasDefX(0);
    myStyle->SetCanvasDefY(0);
    
    // For the Legend
    myStyle->SetLegendBorderSize(0);
    myStyle->SetLegendFillColor(0);
    myStyle->SetLegendFont(42);
    myStyle->SetLegendTextSize(0.04);
    
    // For the Pad:
    myStyle->SetPadBorderMode(0);
    myStyle->SetPadColor(kWhite);
    myStyle->SetGridColor(0);
    myStyle->SetGridStyle(3);
    myStyle->SetGridWidth(1);
    
    // For the frame:
    myStyle->SetFrameBorderMode(0);
    myStyle->SetFrameBorderSize(1);
    myStyle->SetFrameFillColor(0);
    myStyle->SetFrameFillStyle(0);
    myStyle->SetFrameLineColor(1);
    myStyle->SetFrameLineStyle(1);
    myStyle->SetFrameLineWidth(1);
    
    // set the paper & margin sizes
    myStyle->SetPaperSize(20,26);
    myStyle->SetPadTopMargin(0.12);
    myStyle->SetPadRightMargin(0.12);
    myStyle->SetPadBottomMargin(0.18);
    myStyle->SetPadLeftMargin(0.12);
    
    myStyle->SetTitleX(0.5);
    myStyle->SetTitleAlign(23);
    myStyle->SetTitleY(0.985);
    myStyle->SetTitleFont(42);
    myStyle->SetTitleColor(1);
    myStyle->SetTitleTextColor(1);
    myStyle->SetTitleFillColor(0);
    myStyle->SetTitleFontSize(0.05);
    myStyle->SetTitleBorderSize(0);
    
    // use large Helvetica fonts
    myStyle->SetTitleFont(42,"xyz");  // set the all 3 axes title font
    myStyle->SetTitleFont(42," ");    // set the pad title font
    myStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
    myStyle->SetTitleSize(0.06," ");   // set the pad title size
    myStyle->SetLabelFont(42,"xyz");
    myStyle->SetLabelSize(0.05,"xyz");
    myStyle->SetLabelColor(1,"xyz");
    myStyle->SetTextFont(42);
    myStyle->SetTextSize(0.08);
    myStyle->SetStatFont(42);
    
    // use bold lines and markers
    myStyle->SetLineWidth(2);
    myStyle->SetMarkerStyle(8);
    myStyle->SetHistLineWidth(2);
    myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    
    //..Get rid of X error bars
    myStyle->SetErrorX(0.001);
    
    // do not display any of the standard histogram decorations
    myStyle->SetOptTitle(0);
    myStyle->SetOptStat(0);
    myStyle->SetOptFit(11111111);
    
    // put tick marks on top and RHS of plots
    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);
    
    Int_t MyPalette[28];
    Int_t alpha=1;
    
    Double_t Stops[5] = { 0.00, 0.50, 0.70, 0.82, 1.00 };
    Double_t Red[5]   = { 0.00, 0.00, 0.74, 1.00, 1.00 };
    Double_t Green[5] = { 0.00, 0.61, 0.82, 0.70, 1.00 };
    Double_t Blue[5]  = { 0.31, 0.73, 0.08, 0.00, 1.00 };
    
    Int_t zpalette = TColor::CreateGradientColorTable(5, Stops, Red, Green, Blue, 28);
    for(int i=0; i<28; i++)
    {
        MyPalette[i] = zpalette+i;
    }
    myStyle->SetPalette(28,MyPalette);
    
    myStyle->cd();
}

void Draw_Graphs(TFile& fout, vector<TGraph*>& vect_graph, const vector<string>& leg_text, const string& YaxisText, const string& XaxisText, const string& plotName, bool type)
{
    setMyStyle();
    gStyle->SetOptTitle(0);
    if(!(fout.IsOpen())) {cout << "Output File Not Open..." << endl;}
    TLegend* leg = new TLegend(0.18,0.65,0.31,0.86);
    vector<TLegendEntry*> vect_leg_entry;
    {for(int i = 0; i<int(leg_text.size()); i++) { TLegendEntry* leg_entry = leg->AddEntry(vect_graph.at(i),leg_text.at(i).c_str(),"L"); }}
    leg->SetTextSize(0.04);
    TCanvas* canvas_graph = new TCanvas(("canvas_graph"+plotName).c_str(),"canvas_graph",750,500);
    canvas_graph->SetGridx();
    canvas_graph->SetLogx();
    canvas_graph->SetGridy();
    canvas_graph->cd();
    TMultiGraph* mg = new TMultiGraph();
    string title = vect_graph.at(0)->GetTitle();
    title += (";"+XaxisText+";"+YaxisText);
    mg->SetTitle(title.c_str());
    for(int i = 0; i < int(vect_graph.size()); i++) { mg->Add(vect_graph.at(i)); }
    vect_graph[0]->SetLineStyle(1);
    vect_graph[0]->SetMarkerStyle(22);
    vect_graph[0]->SetLineColor(kAzure-1);
    vect_graph[0]->SetMarkerColor(kAzure-1);
    if(vect_graph.size() > 1){
        vect_graph[1]->SetLineStyle(1);
        vect_graph[1]->SetMarkerStyle(22);
        vect_graph[1]->SetLineColor(kOrange+1);
        vect_graph[1]->SetMarkerColor(kOrange+1);}
    if(vect_graph.size() > 2){
        vect_graph[2]->SetLineStyle(1);
        vect_graph[2]->SetMarkerStyle(22);
        vect_graph[2]->SetLineColor(kGreen+2);
        vect_graph[2]->SetMarkerColor(kGreen+2);}
    if(vect_graph.size() == 4){
        vect_graph[0]->SetLineStyle(1);
        vect_graph[0]->SetMarkerStyle(22);
        vect_graph[0]->SetLineColor(kRed);
        vect_graph[0]->SetMarkerColor(kRed);
        vect_graph[1]->SetLineStyle(1);
        vect_graph[1]->SetMarkerStyle(22);
        vect_graph[1]->SetLineColor(kGreen);
        vect_graph[1]->SetMarkerColor(kGreen);
        vect_graph[2]->SetLineStyle(9);
        vect_graph[2]->SetMarkerStyle(32);
        vect_graph[2]->SetLineColor(kBlue);
        vect_graph[2]->SetMarkerColor(kBlue);
        vect_graph[3]->SetLineStyle(9);
        vect_graph[3]->SetMarkerStyle(32);
        vect_graph[3]->SetLineColor(kMagenta);
        vect_graph[3]->SetMarkerColor(kMagenta);}
    if(vect_graph.size() == 5)
    {
        vect_graph[0]->SetLineStyle(1);
        vect_graph[0]->SetMarkerStyle(22);
        vect_graph[0]->SetLineColor(kRed);
        vect_graph[0]->SetMarkerColor(kRed);
        vect_graph[1]->SetLineStyle(1);
        vect_graph[1]->SetMarkerStyle(22);
        vect_graph[1]->SetLineColor(kGreen-4);
        vect_graph[1]->SetMarkerColor(kGreen-4);
        vect_graph[2]->SetLineStyle(1);
        vect_graph[2]->SetMarkerStyle(22);
        vect_graph[2]->SetLineColor(kCyan);
        vect_graph[2]->SetMarkerColor(kCyan);
        vect_graph[3]->SetLineStyle(9);
        vect_graph[3]->SetMarkerStyle(32);
        vect_graph[3]->SetLineColor(kMagenta);
        vect_graph[3]->SetMarkerColor(kMagenta);
        vect_graph[4]->SetLineStyle(9);
        vect_graph[4]->SetMarkerStyle(32);
        vect_graph[4]->SetLineColor(kYellow-4);
        vect_graph[4]->SetMarkerColor(kYellow-4);
    }
    if(vect_graph.size() == 6){
        vect_graph[0]->SetLineStyle(1);
        vect_graph[0]->SetMarkerStyle(22);
        vect_graph[0]->SetLineColor(kBlue);
        vect_graph[0]->SetMarkerColor(kBlue);
        vect_graph[1]->SetLineStyle(1);
        vect_graph[1]->SetMarkerStyle(22);
        vect_graph[1]->SetLineColor(kRed);
        vect_graph[1]->SetMarkerColor(kRed);
        vect_graph[2]->SetLineStyle(1);
        vect_graph[2]->SetMarkerStyle(22);
        vect_graph[2]->SetLineColor(kGreen+2);
        vect_graph[2]->SetMarkerColor(kGreen+2);
        vect_graph[3]->SetLineStyle(9);
        vect_graph[3]->SetMarkerStyle(32);
        vect_graph[3]->SetLineColor(kBlue);
        vect_graph[3]->SetMarkerColor(kBlue);
        vect_graph[4]->SetLineStyle(9);
        vect_graph[4]->SetMarkerStyle(32);
        vect_graph[4]->SetLineColor(kRed);
        vect_graph[4]->SetMarkerColor(kRed);
        vect_graph[5]->SetLineStyle(9);
        vect_graph[5]->SetMarkerStyle(32);
        vect_graph[5]->SetLineColor(kGreen+2);
        vect_graph[5]->SetMarkerColor(kGreen+2);}
    if(type) { mg->Draw("AP"); }
    else { mg->Draw("AL"); }
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetYaxis()->SetTitleOffset(1.05);
    mg->GetXaxis()->SetTitleOffset(1.3);
    mg->GetYaxis()->SetTitleSize(.04);
    mg->GetXaxis()->SetTitleSize(.04);
    mg->GetYaxis()->SetLabelSize(.04);
    mg->GetXaxis()->SetLabelSize(.04);
    leg->Draw("SAMES");
    TLatex l(0.6,0.943,"#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");
    l.SetNDC();
    l.SetTextSize(0.045);
    l.SetTextFont(42);
    l.SetTextColor(kBlack);
    l.DrawLatex(0.62,0.95,"#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");
    l.SetTextSize(0.04);
    l.SetTextFont(42);
    l.SetTextColor(kBlack);
    l.DrawLatex(0.02,0.95,"#bf{#it{RestFrames}} Event Generation");
    //invert colors
    
    mg->GetXaxis()->SetAxisColor(kWhite);
    mg->GetYaxis()->SetAxisColor(kWhite);
    mg->GetXaxis()->SetTitleColor(kWhite);
    mg->GetYaxis()->SetTitleColor(kWhite);
    mg->GetXaxis()->SetLabelColor(kWhite);
    mg->GetYaxis()->SetLabelColor(kWhite);
    canvas_graph->SetFillColor(kBlack);
    leg->SetTextColor(kWhite);
    leg->SetFillColor(kBlack);
    
    canvas_graph->SaveAs((plotName+".pdf").c_str());
    canvas_graph->Write();
}

double evaluateZbi(double Nsig, double Nbkg, double sys)
{
    double Nobs = Nsig+Nbkg;
    double tau = 1./Nbkg/(sys*sys/10000.);
    double aux = Nbkg*tau;
    double Pvalue = TMath::BetaIncomplete(1./(1.+tau),Nobs,aux+1.);
    return sqrt(2.)*TMath::ErfcInverse(Pvalue*2.);
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

TH1F* Histogram_1D(string name, double Nbins, double Xmin, double Xmax, string Xname)
{
    TH1F* hist = new TH1F(("hist_"+name).c_str(),"",Nbins,Xmin,Xmax);
    hist->SetLineColor(kBlack);
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleOffset(1.06);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitle(Xname.c_str());
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleOffset(1.12);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetTitle("Number of Entries");
    return hist;
}

TCanvas* Make_CMS_Canvas(string PlotTitle)
{
    TCanvas* can = new TCanvas(("can_"+PlotTitle).c_str(),PlotTitle.c_str(),700,600);
    can->SetLeftMargin(0.15);
    can->SetRightMargin(0.18);
    can->SetBottomMargin(0.15);
    can->SetGridx();
    can->SetGridy();
    can->Draw();
    can->cd();
    TLatex latex_CMS;
    latex_CMS.SetTextFont(42);
    latex_CMS.SetNDC();
    latex_CMS.SetTextSize(0.035);
    latex_CMS.SetTextFont(42);
    latex_CMS.DrawLatex(0.41,0.943,PlotTitle.c_str());
    latex_CMS.SetTextSize(0.04);
    latex_CMS.SetTextFont(42);
    latex_CMS.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");
    return can;
}

/*
 void Loop_1DHist(TTree fChain, TH1F*& hist, string branchname, int Category)
 {
 if (fChain == 0) return;
 Long64_t nentries = fChain->GetEntriesFast();
 TBranch* branch = fChain->GetBranch(branchname.c_str());
 vector<double> *v_branch_var = 0;
 double branch_var = 0;
 if(Category < 0)
 {
 cout << "Category is less than zero... Assume filling with non-vector variable" << endl;
 fChain->SetBranchAddress(branchname.c_str(), &branch_var, &branch);
 }
 else
 {
 cout << "Category is greater than zero... Assume filling with vector variable" << endl;
 fChain->SetBranchAddress(branchname.c_str(), &v_branch_var, &branch);
 }
 cout << "Filling Histogram with Branch: " << branch->GetName() << endl;
 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
 Long64_t ientry = LoadTree(jentry);
 if (ientry < 0) break;
 // nb = fChain->GetEntry(jentry);   nbytes += nb;
 // if (Cut(ientry) < 0) continue;
 b_weight->GetEntry(ientry);
 branch->GetEntry(ientry);
 // if(jentry > 1000) break; // for quick checks
 if(Category < 0)
 {
 hist->Fill(branch_var, weight);
 }
 else
 {
 hist->Fill(v_branch_var->at(Category), weight);
 }
 }
 }
 */

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
    for(int i = 0; i < int(hists.size()); i++){
        int icolor0 = 7003 + (i%8)*10;
        int icolor1 = 7000 + (i%8)*10;
        hists[i]->SetFillColor(icolor1);
        hists[i]->SetFillStyle(3002);
        hists[i]->SetLineColor(icolor0);
        hists[i]->SetLineWidth(3);
        hists[i]->SetMarkerColor(icolor0);
        hists[i]->SetMarkerSize(0);
        hists[i]->Draw("same");
    }
}

void Draw_Hists(vector<TH1D*> hists, TCanvas* canvas)
{
    canvas->cd();
    for(int i = 0; i < int(hists.size()); i++){
        int icolor0 = 7003 + (i%8)*10;
        int icolor1 = 7000 + (i%8)*10;
        hists[i]->SetFillColor(icolor1);
        hists[i]->SetFillStyle(3002);
        hists[i]->SetLineColor(icolor0);
        hists[i]->SetLineWidth(3);
        hists[i]->SetMarkerColor(icolor0);
        hists[i]->SetMarkerSize(0);
        hists[i]->Draw("same");
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

