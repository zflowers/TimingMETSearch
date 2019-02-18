#include "Detector.hh"
#include "Physics.hh"
#include "Resolution.hh"
#include <TSystem.h>

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

void Fit(std::vector<TH1*> vect_hist)
{
 TCanvas* canvas_Chi2_ctau = new TCanvas("canvas_Chi2_ctau","canvas_Chi2_ctau",750,500);
 TGraph* Chi2_ctau = new TGraph(vect_hist.size());
 Chi2_ctau->SetTitle("#chi2/NDF from fit of Pull Distribution of E_{Za}^{#tilde{#chi}_{2a}^{0}} vs c#tau");
 Chi2_ctau->SetMarkerStyle(22);
 Chi2_ctau->SetMarkerColor(kBlue);
 Chi2_ctau->GetXaxis()->SetTitle("c#tau");
 Chi2_ctau->GetYaxis()->SetTitle("#chi2/NDF");
 canvas_Chi2_ctau->SetLogy();
 canvas_Chi2_ctau->SetGrid();
 TCanvas* canvas_mean_ctau = new TCanvas("canvas_mean_ctau","canvas_mean_ctau",750,500);
 TGraph* mean_ctau = new TGraph(vect_hist.size());
 mean_ctau->SetMarkerStyle(22);
 mean_ctau->SetMarkerColor(kBlue);
 mean_ctau->GetXaxis()->SetTitle("c#tau");
 mean_ctau->GetYaxis()->SetTitle("#mu");
 canvas_mean_ctau->SetGrid();
 mean_ctau->SetTitle("#mu from fit of Pull Distribution of E_{Za}^{#tilde{#chi}_{2a}^{0}} vs c#tau");
 TCanvas* canvas_sigma_ctau = new TCanvas("canvas_sigma_ctau","canvas_sigma_ctau",750,500);
 TGraph* sigma_ctau = new TGraph(vect_hist.size());
 sigma_ctau->SetTitle("#sigma from fit of Pull Distribution of E_{Za}^{#tilde{#chi}_{2a}^{0}} vs c#tau");
 sigma_ctau->SetMarkerStyle(22);
 sigma_ctau->SetMarkerColor(kBlue);
 sigma_ctau->GetXaxis()->SetTitle("c#tau");
 sigma_ctau->GetYaxis()->SetTitle("#sigma");
 canvas_sigma_ctau->SetGrid();
 for(int i = vect_hist.size()-1; i >= 0; i--)
 {
  string name = vect_hist.at(i)->GetName();
  TCanvas* dummy = new TCanvas(name.c_str(),name.c_str(),750,500);
     //cout << name << endl;
  dummy->cd();
  vect_hist.at(i)->Draw();
  vect_hist.at(i)->Fit("gaus","QEM");
  TF1* fit_func = vect_hist.at(i)->GetFunction("gaus");
  string output = "Fit_Plots_EZ/"+name+".pdf";
  dummy->SaveAs(output.c_str());
  Chi2_ctau->SetPoint(i,vect_hist.size()-i,fit_func->GetChisquare()/fit_func->GetNDF());
  mean_ctau->SetPoint(i,vect_hist.size()-i,fit_func->GetParameter(1));
  sigma_ctau->SetPoint(i,vect_hist.size()-i,fit_func->GetParameter(2));
 }
 canvas_Chi2_ctau->cd();
 Chi2_ctau->Draw("AP");
 canvas_Chi2_ctau->SaveAs("Fit_Plots_EZ/Chi2_ctau.pdf");
 canvas_mean_ctau->cd();
 mean_ctau->Draw("AP");
 canvas_mean_ctau->SaveAs("Fit_Plots_EZ/mean_ctau.pdf");
 canvas_sigma_ctau->cd();
 sigma_ctau->Draw("AP");
 canvas_sigma_ctau->SaveAs("Fit_Plots_EZ/sigma_ctau.pdf");
}

void Fitter()
{
 gStyle->SetOptFit(1111);
 string fname = "output_Vertex_LLP_Detector_X2X2_to_ZallXZbllX_EZ.root";
 std::vector<TH1*> vect_hist = list_histos(fname.c_str());
 Fit(vect_hist);
}
