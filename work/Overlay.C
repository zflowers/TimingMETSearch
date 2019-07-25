#include "Bonus.h"
#include "RestFrames/RestFrames.hh"

using namespace RestFrames;

vector<TH1D*> list_histos(string fname, vector<string> dir_names, string hist_name)
{   
    std::vector<TH1D*> vect_hist;
    TKey *key; 
    TFile *f = TFile::Open(fname.c_str(), "READ");
    if(!f || f->IsZombie())
    {   
        cout << "Unable to open " << fname << " for reading..." << endl;
        return vect_hist;
    }
    for(int i = 0; i < int(dir_names.size()); i++)
    {
     TH1D* h = nullptr;
     f->GetObject(("Plots"+dir_names[i]+"/"+"hist/"+hist_name+"_Plots"+dir_names[i]).c_str(),h);
     vect_hist.push_back(h);
    }
    return vect_hist;
}

void Get_MultiHist(TCanvas*& canv, vector<TH1D*>& hists, vector<string> labels)
{
 canv->SetLeftMargin(0.18);
 canv->SetRightMargin(0.12);
 canv->SetBottomMargin(0.15);
 canv->SetGridx();
 canv->SetGridy();
 canv->SetGridx();
 canv->SetGridy();
 canv->cd();
 double ymin = hists[0]->GetYaxis()->GetXmin();
 double ymax = hists[0]->GetYaxis()->GetXmax();
 for(int i = 0; i < int(hists.size()); i++)
 {
  if(hists[i]->GetYaxis()->GetXmax() > ymax) {ymax = hists[i]->GetYaxis()->GetXmax();}
 }
 /*hists[0]->SetLineColor(kBlack);
 if(hists.size() > 1) {hists[1]->SetLineColor(kRed);}
 if(hists.size() > 2) {hists[2]->SetLineColor(kBlue);}
 if(hists.size() > 3) {hists[3]->SetLineColor(kGreen+2);}
 if(hists.size() > 4) {hists[4]->SetLineColor(kMagenta);}
 if(hists.size() > 5) {hists[5]->SetLineColor(kOrange-2);}*/
 for(int k = 0; k < int(hists.size()); k++) 
 {
  //hists[k]->GetXaxis()->SetTitle("M_{#tilde{#chi}_{2a}^{0}} [GeV]");
  //hists[k]->GetYaxis()->SetTitle("Fraction of Events in Each Bin");
  //hists[k]->SetFillStyle(0);
  hists[k]->Draw("SAMES HIST"); 
 }
 TLegend* leg = new TLegend(0.65,0.6,0.92,0.91);
 {for(int i = 0; i<int(labels.size()); i++) { TLegendEntry* leg_entry = leg->AddEntry(hists.at(i),labels.at(i).c_str(),"L"); }}
 leg->Draw("SAMES");
}

void Overlay()
{
 SetStyle();
 vector<string> directories{"", "_Cos"};

 string inFile = "output_ctau_X2X2_to_ZallXZbllX.root";
 vector<TH1D*> mLLP_20cm = list_histos(inFile, directories, "MXa2_ctau_0");
 vector<TH1D*> mLLP_10cm = list_histos(inFile, directories, "MXa2_ctau_1");
 vector<TH1D*> mLLP_5cm = list_histos(inFile, directories, "MXa2_ctau_2");

 TFile* outFile = new TFile(("out_"+inFile).c_str(),"RECREATE");

 string can_name_mLLP_20cm = "canv_mLLP_20cm";
 TCanvas* canv_mLLP_20cm = new TCanvas(can_name_mLLP_20cm.c_str(),"",750,600);
 Get_MultiHist(canv_mLLP_20cm, mLLP_20cm, directories);
 canv_mLLP_20cm->Write();

 string can_name_mLLP_10cm = "canv_mLLP_10cm";
 TCanvas* canv_mLLP_10cm = new TCanvas(can_name_mLLP_10cm.c_str(),"",750,600);
 Get_MultiHist(canv_mLLP_10cm, mLLP_10cm, directories);
 canv_mLLP_10cm->Write();

 string can_name_mLLP_5cm = "canv_mLLP_5cm";
 TCanvas* canv_mLLP_5cm = new TCanvas(can_name_mLLP_5cm.c_str(),"",750,600);
 Get_MultiHist(canv_mLLP_5cm, mLLP_5cm, directories);
 canv_mLLP_5cm->Write();
}
