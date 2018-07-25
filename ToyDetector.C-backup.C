#define ToyDetector_cxx
#include "ToyDetector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../Detector.hh"
#include "Bonus.h"

void ToyDetector::Loop()
{
    //TFile* Plots = new TFile("PlotsDetector.root","UPDATE");
    TFile* Plots = new TFile("PlotsDetector.root","RECREATE");
    CreatePalette();
    TLorentzVector n2A_MC, n2B_MC, n1A_MC, n1B_MC, eA_MC, eB_MC, eC_MC, eD_MC, mA_MC, mB_MC, mC_MC, mD_MC;
    TLorentzVector PV_RECO, SV_A_RECO, SV_B_RECO;
    TLorentzVector eA_RECO, eB_RECO, eC_RECO, eD_RECO, mA_RECO, mB_RECO, mC_RECO, mD_RECO;
    int dummycount=0;
    
    double true_SumPT=0.0;
    gRandom = new TRandom3();
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    
   if (fChain == 0) return;
    
    TCanvas* c1 = new TCanvas("CanvasforPUPPIMET","Canvas for PUPPI MET",750,500);
    TH1F *hist_True_MET_PUPPI = new TH1F("hist_True_MET_PUPPI","True MET PUPPI",100,0.0,1400.0);
    TH1F *hist_Test_MET_PUPPI = new TH1F("hist_Test_MET_PUPPI","Test MET PUPPI",100,0.0,1400.0);
    
    TCanvas* c2 = new TCanvas("CanvasforMET","Canvas for MET",750,500);
    TH1F *hist_True_MET = new TH1F("hist_True_MET","True MET",100,0.0,800.0);
    TH1F *hist_Test_MET = new TH1F("hist_Test_MET","Test MET",100,0.0,800.0);
    
    TCanvas* c3 = new TCanvas("Canvasfor2DPTvsEta","Canvas for 2D PT vs Eta",750,500);
    TH2F *hist_PTvsEta = new TH2F("hist_PTvsEta","PT vs Eta",100,-10.0,10.0,100,0.0,1000.0);
    
    TCanvas* c4 = new TCanvas("Canvasfor2DPTvsMass","Canvas for 2D PT vs Mass",750,500);
    TH2F *hist_PTvsMass = new TH2F("hist_PTvsMass","PT vs Mass",100,350.0,1500.0,100,0.0,800.0);
    
    TCanvas* c5 = new TCanvas("CanvasforElectron","Canvas for Electron",750,500);
    TH1F* hist_Test_Electron = new TH1F("hist_Test_Electron","Test Electron",100,0.0,800.0);
    TH1F* hist_True_Electron = new TH1F("hist_True_Electron","True Electron",100,0.0,800.0);
    
    TCanvas* c6 = new TCanvas("CanvasforMuon","Canvas for Muon",750,500);
    TH1F* hist_Test_Muon = new TH1F("hist_Test_Muon","Test Muon",100,0.0,800.0);
    TH1F* hist_True_Muon = new TH1F("hist_True_Muon","True Muon",100,0.0,800.0);
    
    TCanvas* c7 = new TCanvas("CanvasforDeltaRElectron","Canvas for DeltaR Electron",750,500);
    TH1F* hist_DeltaR_Electron = new TH1F("hist_DeltaR_Electron","Hist DeltaR Electron",100,0.0,0.015);
    
    TCanvas* c8 = new TCanvas("CanvasforDeltaPtElectron","Canvas for DeltaPt Electron",750,500);
    TH1F* hist_DeltaPt_Electron = new TH1F("hist_DeltaPt_Electron","Hist DeltaPt Electron",100,-25.0,25.0);
    
    TCanvas* c9 = new TCanvas("CanvasforDeltaRMuon","Canvas for DeltaR Muon",750,500);
    TH1F* hist_DeltaR_Muon = new TH1F("hist_DeltaR_Muon","Hist DeltaR Muon",100,0.0,0.015);
    
    TCanvas* c10 = new TCanvas("CanvasforDeltaPtMuon","Canvas for DeltaPt Muon",750,500);
    TH1F* hist_DeltaPt_Muon = new TH1F("hist_DeltaPt_Muon","Hist DeltaPt Muon",100,-25.0,25.0);
    
    TCanvas* c11 = new TCanvas("Canvasfor2DElectronDeltaRDeltaPt","Canvas for 2D Electron DeltaR and DeltaPt",750,500);
    TH2F* hist_DeltaR_DeltaPt_Electron = new TH2F("hist_DeltaR_DeltaPtElectron","Hist DeltaR and DeltaPt Electron",100,0.0,0.015,100,-25.0,25.0);

    TCanvas* c12 = new TCanvas("Canvasfor2DMuonDeltaRDeltaPt","Canvas for 2D Muon DeltaR and DeltaPt",750,500);
    TH2F* hist_DeltaR_DeltaPt_Muon = new TH2F("hist_DeltaR_DeltaPtMuon","Hist DeltaR and DeltaPt Muon",100,0.0,0.015,100,-25.0,25.0);
    
    TCanvas* c13 = new TCanvas("CanvasforDeltaRElectronMC","Canvas for DeltaR Electron MC",750,500);
    TH1F* hist_DeltaR_ElectronMC = new TH1F("hist_DeltaR_ElectronMC","Hist DeltaR Electron MC",100,0.0,0.015);
    
    TCanvas* c14 = new TCanvas("CanvasforDeltaPtElectronMC","Canvas for DeltaPt Electron MC",750,500);
    TH1F* hist_DeltaPt_ElectronMC = new TH1F("hist_DeltaPt_ElectronMC","Hist DeltaPt Electron MC",100,-25.0,25.0);
    
    TCanvas* c15 = new TCanvas("CanvasforDeltaRMuonMC","Canvas for DeltaR Muon MC",750,500);
    TH1F* hist_DeltaR_MuonMC = new TH1F("hist_DeltaR_MuonMC","Hist DeltaR Muon MC",100,0.0,0.015);
    
    TCanvas* c16 = new TCanvas("CanvasforDeltaPtMuonMC","Canvas for DeltaPt Muon MC",750,500);
    TH1F* hist_DeltaPt_MuonMC = new TH1F("hist_DeltaPt_MuonMC","Hist DeltaPt Muon MC",100,-25.0,25.0);
    
    TCanvas* c17 = new TCanvas("Canvasfor2DElectronMCDeltaRDeltaPt","Canvas for 2D Electron DeltaR and DeltaPt MC",750,500);
    TH2F* hist_DeltaR_DeltaPt_ElectronMC = new TH2F("hist_DeltaR_DeltaPtElectronMC","Hist DeltaR and DeltaPt MC",100,0.0,0.015,100,-25.0,25.0);
    
    TCanvas* c18 = new TCanvas("Canvasfor2DMuonMCDeltaRDeltaPt","Canvas for 2D Muon DeltaR and DeltaPt MC",750,500);
    TH2F* hist_DeltaR_DeltaPt_MuonMC = new TH2F("hist_DeltaR_DeltaPtMuonMC","Hist DeltaR and DeltaPt MC",100,0.0,0.015,100,-25.0,25.0);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       n2A_MC.SetPtEtaPhiE(n2_MC_PT->at(0), n2_MC_Eta->at(0), n2_MC_Phi->at(0), n2_MC_E->at(0));
       n2B_MC.SetPtEtaPhiE(n2_MC_PT->at(1), n2_MC_Eta->at(1), n2_MC_Phi->at(1), n2_MC_E->at(1));
       
       n1A_MC.SetPtEtaPhiE(n1_MC_PT->at(0), n1_MC_Eta->at(0), n1_MC_Phi->at(0), n1_MC_E->at(0));
       n1B_MC.SetPtEtaPhiE(n1_MC_PT->at(1), n1_MC_Eta->at(1), n1_MC_Phi->at(1), n1_MC_E->at(1));
       
       Detector test_Lepton;
       Detector test_MET;
       Detector test_MET_PUPPI;
       TLorentzVector PV_RECO, SV_A_RECO, SV_B_RECO;
       TLorentzVector eA_RECO, eB_RECO, eC_RECO, eD_RECO, mA_RECO, mB_RECO, mC_RECO, mD_RECO;

       test_MET.Set_con0_par(59.8722);
       test_MET.Set_con1_par(1.01476);
       test_MET.Set_con2_par(0.0407631);
       test_MET.Set_con0_perp(59.9967);
       test_MET.Set_con1_perp(0.904125);
       test_MET.Set_con2_perp(0.0226999);
       
       test_MET_PUPPI.Set_con0_par(14.7467);
       test_MET_PUPPI.Set_con1_par(1.68788);
       test_MET_PUPPI.Set_con2_par(-2.95569e-07);
       test_MET_PUPPI.Set_con0_perp(15.4667);
       test_MET_PUPPI.Set_con1_perp(1.41597);
       test_MET_PUPPI.Set_con2_perp(-6.37947e-07);
       
       TLorentzVector n2 = n2A_MC + n2B_MC;
       TVector3 inv = (n1A_MC + n1B_MC).Vect();
       inv.SetZ(0.0);
       
       TVector3 MET_Delphes;
       MET_Delphes.SetPtEtaPhi(MET_PT,MET_Eta,MET_Phi);
       MET_Delphes.SetZ(0.0);
       hist_True_MET->Fill(MET_Delphes.Mag());
       
       TVector3 MET_Delphes_PUPPI;
       MET_Delphes_PUPPI.SetPtEtaPhi(PUPPI_MET_PT,PUPPI_MET_Eta,PUPPI_MET_Phi);
       MET_Delphes_PUPPI.SetZ(0.0);
       hist_True_MET_PUPPI->Fill(MET_Delphes_PUPPI.Mag());

       TVector3 MET_RECO = test_MET.Smear_MET(n2, inv);
       hist_Test_MET->Fill(MET_RECO.Mag());
       
       TVector3 MET_RECO_PUPPI = test_MET_PUPPI.Smear_MET(n2, inv);
       hist_Test_MET_PUPPI->Fill(MET_RECO_PUPPI.Mag());
       
       hist_PTvsEta->Fill(n2.Eta(),n2.Pt());
       hist_PTvsMass->Fill(n2.M(),n2.Pt());
       
       //Run on leptons coming from a Z
       
       std::vector<TLorentzVector> vect_Electron_Delphes;
       for(int k = 0; k<ElectronSignal_MC_Size; k++)
       {
           if(ElectronSignal_MC_Size == 0) break;
           TLorentzVector ElectronSignal_MC;
           ElectronSignal_MC.SetPtEtaPhiE(ElectronSignal_MC_PT->at(k),ElectronSignal_MC_Eta->at(k),ElectronSignal_MC_Phi->at(k),ElectronSignal_MC_E->at(k));
           for(int i = 0; i<Electron_Delphes_Size; i++)
           {
               TLorentzVector Electron_Test;
               Electron_Test.SetPtEtaPhiE(Electron_Delphes_PT->at(i),Electron_Delphes_Eta->at(i),Electron_Delphes_Phi->at(i),Electron_Delphes_E->at(i));
               if(ElectronSignal_MC.DeltaR(Electron_Test) < 0.01)
               {
                   //cout << Electron_Delphes_E->at(i) << endl;
                   vect_Electron_Delphes.push_back(Electron_Test);
               }
           }
       }
       if(vect_Electron_Delphes.size() == ElectronSignal_MC_Size)
       {
           for(int i = 0; i<ElectronSignal_MC_Size; i++)
           {
               TLorentzVector ElectronSignal_MC;
               ElectronSignal_MC.SetPtEtaPhiE(ElectronSignal_MC_PT->at(i),ElectronSignal_MC_Eta->at(i),ElectronSignal_MC_Phi->at(i),ElectronSignal_MC_E->at(i));
               TLorentzVector Smeared = test_Lepton.Smear_Electron(ElectronSignal_MC);
               hist_Test_Electron->Fill(Smeared.Pt());
               hist_True_Electron->Fill(vect_Electron_Delphes.at(i).Pt());
               hist_DeltaR_Electron->Fill(vect_Electron_Delphes.at(i).DeltaR(Smeared));
               hist_DeltaPt_Electron->Fill(vect_Electron_Delphes.at(i).Pt()-Smeared.Pt());
               hist_DeltaR_DeltaPt_Electron->Fill(vect_Electron_Delphes.at(i).DeltaR(Smeared), vect_Electron_Delphes.at(i).Pt()-Smeared.Pt());
               hist_DeltaR_ElectronMC->Fill(vect_Electron_Delphes.at(i).DeltaR(ElectronSignal_MC));
               hist_DeltaPt_ElectronMC->Fill(vect_Electron_Delphes.at(i).Pt()-ElectronSignal_MC.Pt());
               hist_DeltaR_DeltaPt_ElectronMC->Fill(vect_Electron_Delphes.at(i).DeltaR(ElectronSignal_MC), vect_Electron_Delphes.at(i).Pt()-ElectronSignal_MC.Pt());
           }
       }
       
       std::vector<TLorentzVector> vect_Muon_Delphes;
       for(int k = 0; k<MuonSignal_MC_Size; k++)
       {
           if(MuonSignal_MC_Size == 0) break;
           TLorentzVector MuonSignal_MC;
           MuonSignal_MC.SetPtEtaPhiE(MuonSignal_MC_PT->at(k),MuonSignal_MC_Eta->at(k),MuonSignal_MC_Phi->at(k),MuonSignal_MC_E->at(k));
           for(int i = 0; i<MuonLoose_Delphes_Size; i++)
           {
               TLorentzVector Muon_Test;
               Muon_Test.SetPtEtaPhiE(MuonLoose_Delphes_PT->at(i),MuonLoose_Delphes_Eta->at(i),MuonLoose_Delphes_Phi->at(i),MuonLoose_Delphes_E->at(i));
               if(MuonSignal_MC.DeltaR(Muon_Test) < 0.01)
               {
                   vect_Muon_Delphes.push_back(Muon_Test);
               }
           }
       }
       if(vect_Muon_Delphes.size() == MuonSignal_MC_Size)
       {
           for(int i = 0; i<MuonSignal_MC_Size; i++)
           {
               TLorentzVector MuonSignal_MC;
               MuonSignal_MC.SetPtEtaPhiE(MuonSignal_MC_PT->at(i),MuonSignal_MC_Eta->at(i),MuonSignal_MC_Phi->at(i),MuonSignal_MC_E->at(i));
               TLorentzVector Smeared = test_Lepton.Smear_Muon(MuonSignal_MC);
               hist_Test_Muon->Fill(Smeared.Pt());
               hist_True_Muon->Fill(vect_Muon_Delphes.at(i).Pt());
               hist_DeltaR_Muon->Fill(vect_Muon_Delphes.at(i).DeltaR(Smeared));
               hist_DeltaPt_Muon->Fill(vect_Muon_Delphes.at(i).Pt()-Smeared.Pt());
               hist_DeltaR_DeltaPt_Muon->Fill(vect_Muon_Delphes.at(i).DeltaR(Smeared), vect_Muon_Delphes.at(i).Pt()-Smeared.Pt());
               hist_DeltaR_MuonMC->Fill(vect_Muon_Delphes.at(i).DeltaR(MuonSignal_MC));
               hist_DeltaPt_MuonMC->Fill(vect_Muon_Delphes.at(i).Pt()-MuonSignal_MC.Pt());
               hist_DeltaR_DeltaPt_MuonMC->Fill(vect_Muon_Delphes.at(i).DeltaR(MuonSignal_MC), vect_Muon_Delphes.at(i).Pt()-MuonSignal_MC.Pt());
           }
       }
       
       //Run on all leptons
       /*
       std::vector<TLorentzVector> vect_Electron_Delphes;
       for(int k = 0; k<Electron_MC_Size; k++)
       {
           if(Electron_MC_Size == 0) break;
           TLorentzVector Electron_MC;
           Electron_MC.SetPtEtaPhiE(Electron_MC_PT->at(k),Electron_MC_Eta->at(k),Electron_MC_Phi->at(k),Electron_MC_E->at(k));
           for(int i = 0; i<Electron_Delphes_Size; i++)
           {
               TLorentzVector Electron_Test;
               Electron_Test.SetPtEtaPhiE(Electron_Delphes_PT->at(i),Electron_Delphes_Eta->at(i),Electron_Delphes_Phi->at(i),Electron_Delphes_E->at(i));
               if(Electron_MC.DeltaR(Electron_Test) < 0.01)
               {
                   //cout << Electron_Delphes_E->at(i) << endl;
                   vect_Electron_Delphes.push_back(Electron_Test);
               }
           }
       }
       if(vect_Electron_Delphes.size() == Electron_MC_Size)
       {
           for(int i = 0; i<Electron_MC_Size; i++)
           {
               TLorentzVector Electron_MC;
               Electron_MC.SetPtEtaPhiE(Electron_MC_PT->at(i),Electron_MC_Eta->at(i),Electron_MC_Phi->at(i),Electron_MC_E->at(i));
               TLorentzVector Smeared = test_Lepton.Smear_Electron(Electron_MC);
               hist_Test_Electron->Fill(Smeared.Pt());
               hist_True_Electron->Fill(vect_Electron_Delphes.at(i).Pt());
           }
       }
       
       std::vector<TLorentzVector> vect_Muon_Delphes;
       for(int k = 0; k<Muon_MC_Size; k++)
       {
           if(Muon_MC_Size == 0) break;
           TLorentzVector Muon_MC;
           Muon_MC.SetPtEtaPhiE(Muon_MC_PT->at(k),Muon_MC_Eta->at(k),Muon_MC_Phi->at(k),Muon_MC_E->at(k));
           for(int i = 0; i<MuonLoose_Delphes_Size; i++)
           {
               TLorentzVector Muon_Test;
               Muon_Test.SetPtEtaPhiE(MuonLoose_Delphes_PT->at(i),MuonLoose_Delphes_Eta->at(i),MuonLoose_Delphes_Phi->at(i),MuonLoose_Delphes_E->at(i));
               if(Muon_MC.DeltaR(Muon_Test) < 0.01)
               {
                   vect_Muon_Delphes.push_back(Muon_Test);
               }
           }
       }
       if(vect_Muon_Delphes.size() == Muon_MC_Size)
       {
           for(int i = 0; i<Muon_MC_Size; i++)
           {
               TLorentzVector Muon_MC;
               Muon_MC.SetPtEtaPhiE(Muon_MC_PT->at(i),Muon_MC_Eta->at(i),Muon_MC_Phi->at(i),Muon_MC_E->at(i));
               TLorentzVector Smeared = test_Lepton.Smear_Muon(Muon_MC);
               hist_Test_Muon->Fill(Smeared.Pt());
               hist_True_Muon->Fill(vect_Muon_Delphes.at(i).Pt());
           }
       }
        */
   }
    cout << endl << dummycount << endl;
    //cout << endl << dummycount << "/" << nentries << endl;
    
    Draw_Two_Hists(hist_Test_MET_PUPPI,hist_True_MET_PUPPI,c1);
    hist_Test_MET_PUPPI->GetXaxis()->SetTitle("PUPPI MET [GeV]");
    hist_Test_MET_PUPPI->GetYaxis()->SetTitle("Number of Entries");
    c1->SaveAs("MET_PUPPI_Test.png");
    c1->Update();
    c1->Write();
    
    Draw_Two_Hists(hist_Test_MET,hist_True_MET,c2);
    hist_Test_MET->GetXaxis()->SetTitle("MET [GeV]");
    hist_Test_MET->GetYaxis()->SetTitle("Number of Entries");
    c2->SaveAs("MET_Test.png");
    c2->Update();
    c2->Write();
    
    c3->cd();
    hist_PTvsEta->GetXaxis()->SetTitle("Eta");
    hist_PTvsEta->GetYaxis()->SetTitle("PT [GeV]");
    hist_PTvsEta->Draw("COLZ");
    c3->SaveAs("PTvsEta.png");
    c3->Update();
    //c3->Write();
    
    c4->cd();
    hist_PTvsEta->GetXaxis()->SetTitle("Mass [GeV]");
    hist_PTvsEta->GetXaxis()->SetTitle("PT [GeV]");
    hist_PTvsMass->Draw("COLZ");
    c4->SaveAs("PTvsMass.png");
    c4->Update();
    //c4->Write();
    
    Draw_Two_Hists(hist_Test_Electron,hist_True_Electron,c5);
    hist_Test_Electron->GetXaxis()->SetTitle("Electron PT [GeV]");
    hist_Test_Electron->GetYaxis()->SetTitle("Number of Entries");
    c5->SaveAs("Electron_Test.png");
    c5->Update();
    c5->Write();
    
    Draw_Two_Hists(hist_Test_Muon,hist_True_Muon,c6);
    hist_Test_Muon->GetXaxis()->SetTitle("Muon PT [GeV]");
    hist_Test_Muon->GetYaxis()->SetTitle("Number of Entries");
    c6->SaveAs("Muon_Test.png");
    c6->Update();
    c6->Write();
    
    c7->cd();
    hist_DeltaR_Electron->GetXaxis()->SetTitle("DeltaR Between Toy Smeared Electrons and Delphes Smeared Electrons");
    hist_DeltaR_Electron->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaR_Electron->Draw();
    c7->Update();
    c7->Write();
    
    c8->cd();
    hist_DeltaPt_Electron->GetXaxis()->SetTitle("DeltaPt Between Toy Smeared Electrons and Delphes Smeared Electrons [GeV]");
    hist_DeltaPt_Electron->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaPt_Electron->Draw();
    c8->Update();
    c8->Write();
    
    c9->cd();
    hist_DeltaR_Muon->GetXaxis()->SetTitle("DeltaR Between Toy Smeared Muons and Delphes Smeared Muons");
    hist_DeltaR_Muon->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaR_Muon->Draw();
    c9->Update();
    c9->Write();
    
    c10->cd();
    hist_DeltaPt_Muon->GetXaxis()->SetTitle("DeltaPt Between Toy Smeared Muons and Delphes Smeared Muons [GeV]");
    hist_DeltaPt_Muon->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaPt_Muon->Draw();
    c10->Update();
    c10->Write();

    c11->cd();
    hist_DeltaR_DeltaPt_Electron->GetXaxis()->SetTitle("DeltaR Electron");
    hist_DeltaR_DeltaPt_Electron->GetYaxis()->SetTitle("DeltaPt Electron");
    hist_DeltaR_DeltaPt_Electron->Draw("COLZ");
    c11->Update();
    c11->Write();
    
    c12->cd();
    hist_DeltaR_DeltaPt_Muon->GetXaxis()->SetTitle("DeltaR Muon");
    hist_DeltaR_DeltaPt_Muon->GetYaxis()->SetTitle("DeltaPt Muon");
    hist_DeltaR_DeltaPt_Muon->Draw("COLZ");
    c12->Update();
    c12->Write();
    
    c13->cd();
    hist_DeltaR_ElectronMC->GetXaxis()->SetTitle("DeltaR ElectronMC");
    hist_DeltaR_ElectronMC->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaR_ElectronMC->Draw();
    c13->Update();
    c13->Write();
    
    c14->cd();
    hist_DeltaPt_ElectronMC->GetXaxis()->SetTitle("DeltaPt ElectronMC");
    hist_DeltaPt_ElectronMC->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaPt_ElectronMC->Draw();
    c14->Update();
    c14->Write();
    
    c15->cd();
    hist_DeltaR_MuonMC->GetXaxis()->SetTitle("DeltaR MuonMC");
    hist_DeltaR_MuonMC->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaR_MuonMC->Draw();
    c15->Update();
    c15->Write();
    
    c16->cd();
    hist_DeltaPt_MuonMC->GetXaxis()->SetTitle("DeltaPt MuonMC");
    hist_DeltaPt_MuonMC->GetYaxis()->SetTitle("Number of Entries");
    hist_DeltaPt_MuonMC->Draw();
    c16->Update();
    c16->Write();
    
    c17->cd();
    hist_DeltaR_DeltaPt_ElectronMC->GetXaxis()->SetTitle("DeltaR ElectronMC");
    hist_DeltaR_DeltaPt_ElectronMC->GetYaxis()->SetTitle("DeltaPt ElectronMC");
    hist_DeltaR_DeltaPt_ElectronMC->Draw("COLZ");
    c17->Update();
    c17->Write();
    
    c18->cd();
    hist_DeltaR_DeltaPt_MuonMC->GetXaxis()->SetTitle("DeltaR MuonMC");
    hist_DeltaR_DeltaPt_MuonMC->GetYaxis()->SetTitle("DeltaPt MuonMC");
    hist_DeltaR_DeltaPt_MuonMC->Draw("COLZ");
    c18->Update();
    c18->Write();
    
    Plots->Close();
}
