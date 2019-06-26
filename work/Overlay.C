#include "Bonus.h"
//#include "RestFrames/RestFrames.hh"
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
     f->GetObject(("Plots_"+dir_names[i]+"/"+"hist/"+hist_name+"_Plots_"+dir_names[i]).c_str(),h);
     cout << h->GetName() << endl;
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
 hists[0]->SetLineColor(kBlack);
 if(hists.size() > 1) {hists[1]->SetLineColor(kRed);}
 if(hists.size() > 2) {hists[2]->SetLineColor(kBlue);}
 if(hists.size() > 3) {hists[3]->SetLineColor(kGreen+2);}
 if(hists.size() > 4) {hists[4]->SetLineColor(kMagenta);}
 if(hists.size() > 5) {hists[5]->SetLineColor(kOrange-2);}
 for(int k = 0; k < int(hists.size()); k++) 
 {
  hists[k]->SetFillStyle(0); 
  hists[k]->Draw("SAMES HIST"); 
 }
 TLegend* leg = new TLegend(0.65,0.6,0.92,0.91);
 {for(int i = 0; i<int(labels.size()); i++) { TLegendEntry* leg_entry = leg->AddEntry(hists.at(i),labels.at(i).c_str(),"L"); }}
 leg->Draw("SAMES");
}

void Overlay()
{
 setMyStyle();
 vector<string> directories{"NoCut", "Cut"};

 string inFile = "output_LSP_testing.root";
 vector<TH1D*> mLLP_25cm = list_histos(inFile, directories, "MXa2_ctau_0");

 TFile* outFile = new TFile(("out_"+inFile).c_str(),"RECREATE");

 string can_name_mLLP_25cm = "canv_mLLP_25cm";
 TCanvas* canv_mLLP_25cm = new TCanvas(can_name_mLLP_25cm.c_str(),"",750,600);
 Get_MultiHist(canv_mLLP_25cm, mLLP_25cm, directories);
 canv_mLLP_25cm->Write();
}
