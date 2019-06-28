/////////////////////////////////////////////////////////////////////////
//   RestFrames: particle physics event analysis library
//   --------------------------------------------------------------------
//   Copyright (c) 2014-2016, Christopher Rogan
/////////////////////////////////////////////////////////////////////////
///
///  \file   example_X2X2_to_ZllXHggX.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2015 Nov
//
//   This file is part of RestFrames.
//
//   RestFrames is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   RestFrames is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with RestFrames. If not, see <http://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////////////////////////

#include "RestFrames/RestFrames.hh"

#include "Detector.hh"
#include "Physics.hh"
#include "Resolution.hh"
#include <TGraph.h>
#include <TSystem.h>
#include "Bonus.h"

using namespace RestFrames;

void Draw_Graphs(TFile& fout, vector<TGraph*>& vect_graph, const vector<string>& leg_text, const string& YaxisText, const string& XaxisText, const string& plotName);

void Efficiency(std::string output_name =
			      "Efficiency.root"){
    setMyStyle();
    Long64_t start = gSystem->Now();
    Long64_t end = 0.;
    
    //setting masses and widths
    double mX2 = 1000.0;
    double mX1 = 500.0;
    double mZ = 91.19;
    double wZ = 2.50;
    
    vector<double> ctau;
    
    //ctau.push_back(25.);
    ctau.push_back(10.);
    //ctau.push_back(5.);
    
    int Nctau = ctau.size();
    vector<double> sigmaT;
    vector<double> sigmaMET;
    
    for(double i = 10.; i <= 350.; i+=10.)
    {
        sigmaMET.push_back(i);
        sigmaT.push_back(i);
    }
    
    //sigmaT.push_back(30.);
    int NsigmaT = sigmaT.size();
    int NsigmaMET = sigmaMET.size();
    
    //Number of events
    int Ngen = 10000;
    
    vector<vector<int>> vect_timing_acp(Nctau, vector<int>(NsigmaT,0));
    vector<vector<int>> vect_timing_decayangle_acp(Nctau, vector<int>(NsigmaT,0));
    vector<vector<int>> vect_met_acp(Nctau, vector<int>(NsigmaMET,0));
    
    vector<TGraph*> vect_graph_Efficiency_SigmaT;
    vector<TGraph*> vect_graph_Efficiency_SigmaT_DecayAngle;
    vector<TGraph*> vect_graph_Efficiency_SigmaMET;
    for(int j = 0; j < Nctau; j++)
    {
        TGraph* graph_Efficiency_SigmaT = new TGraph(NsigmaT);
        graph_Efficiency_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        vect_graph_Efficiency_SigmaT.push_back(graph_Efficiency_SigmaT);
        TGraph* graph_Efficiency_SigmaT_DecayAngle = new TGraph(NsigmaT);
        graph_Efficiency_SigmaT_DecayAngle->SetName(("graph_sigmaT_DecayAngle"+std::to_string(j)).c_str());
        vect_graph_Efficiency_SigmaT_DecayAngle.push_back(graph_Efficiency_SigmaT_DecayAngle);
        TGraph* graph_Efficiency_SigmaMET = new TGraph(NsigmaMET);
        graph_Efficiency_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        vect_graph_Efficiency_SigmaMET.push_back(graph_Efficiency_SigmaMET);
    }
    
    g_Log << LogInfo << "Initializing generator frames and tree..." << LogEnd;
    
    //Setting up the tree: x2x2->X1X1 + ZZ->4l
    ppLabGenFrame LAB_Gen("LAB_Gen","LAB");
    DecayGenFrame     X2X2_Gen("X2X2_Gen","#tilde{#chi}^{ 0}_{2} #tilde{#chi}^{ 0}_{2}");
    DecayGenFrame     X2a_Gen("X2a_Gen","#tilde{#chi}^{ 0}_{2 a}");
    DecayGenFrame     X2b_Gen("X2b_Gen","#tilde{#chi}^{ 0}_{2 b}");
    ResonanceGenFrame Za_Gen("Za_Gen","Z_{a}");
    ResonanceGenFrame Zb_Gen("Zb_Gen","Z_{b}");
    VisibleGenFrame   L1a_Gen("L1a_Gen","#it{l}_{1a}");
    VisibleGenFrame   L2a_Gen("L2a_Gen","#it{l}_{2a}");
    VisibleGenFrame   L1b_Gen("L1b_Gen","#it{l}_{1b}");
    VisibleGenFrame   L2b_Gen("L2b_Gen","#it{l}_{2b}");
    InvisibleGenFrame X1a_Gen("X1a_Gen","#tilde{#chi}^{ 0}_{1 a}");
    InvisibleGenFrame X1b_Gen("X1b_Gen","#tilde{#chi}^{ 0}_{1 b}");
    
    LabRecoFrame LAB_Reco("LAB_Reco","LAB");
    DecayRecoFrame     X2X2_Reco("X2X2_Reco","#tilde{#chi}^{ 0}_{2} #tilde{#chi}^{ 0}_{2}");
    DecayRecoFrame     X2a_Reco("X2a_Reco","#tilde{#chi}^{ 0}_{2 a}");
    DecayRecoFrame     X2b_Reco("X2b_Reco","#tilde{#chi}^{ 0}_{2 b}");
    DecayRecoFrame Za_Reco("Za_Reco","Z_{a}");
    DecayRecoFrame Zb_Reco("Zb_Reco","Z_{b}");
    VisibleRecoFrame   L1a_Reco("L1a_Reco","#it{l}_{1a}");
    VisibleRecoFrame   L2a_Reco("L2a_Reco","#it{l}_{2a}");
    VisibleRecoFrame   L1b_Reco("L1b_Reco","#it{l}_{1b}");
    VisibleRecoFrame   L2b_Reco("L2b_Reco","#it{l}_{2b}");
    VisibleRecoFrame X1a_Reco("X1a_Reco","#tilde{#chi}^{ 0}_{1 a}");
    VisibleRecoFrame X1b_Reco("X1b_Reco","#tilde{#chi}^{ 0}_{1 b}");
    
    //
    LAB_Gen.SetChildFrame(X2X2_Gen);
    X2X2_Gen.AddChildFrame(X2a_Gen);
    X2X2_Gen.AddChildFrame(X2b_Gen);
    X2a_Gen.AddChildFrame(Za_Gen);
    X2a_Gen.AddChildFrame(X1a_Gen);
    X2b_Gen.AddChildFrame(Zb_Gen);
    X2b_Gen.AddChildFrame(X1b_Gen);
    Za_Gen.AddChildFrame(L1a_Gen);
    Za_Gen.AddChildFrame(L2a_Gen);
    Zb_Gen.AddChildFrame(L1b_Gen);
    Zb_Gen.AddChildFrame(L2b_Gen);
    
    LAB_Reco.SetChildFrame(X2X2_Reco);
    X2X2_Reco.AddChildFrame(X2a_Reco);
    X2X2_Reco.AddChildFrame(X2b_Reco);
    X2a_Reco.AddChildFrame(Za_Reco);
    X2a_Reco.AddChildFrame(X1a_Reco);
    X2b_Reco.AddChildFrame(Zb_Reco);
    X2b_Reco.AddChildFrame(X1b_Reco);
    Za_Reco.AddChildFrame(L1a_Reco);
    Za_Reco.AddChildFrame(L2a_Reco);
    Zb_Reco.AddChildFrame(L1b_Reco);
    Zb_Reco.AddChildFrame(L2b_Reco);
    
    if(LAB_Gen.InitializeTree())
        g_Log << LogInfo << "...Successfully initialized generator tree" << LogEnd;
    else
        g_Log << LogError << "...Failed initializing generator tree" << LogEnd;
    if(LAB_Reco.InitializeTree())
        g_Log << LogInfo << "...Successfully initialized reconstructed tree" << LogEnd;
    else
        g_Log << LogError << "...Failed initializing reconstructed tree" << LogEnd;
    
  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
    
  
  X2X2_Gen.SetVariableMass(); //Reset this down below from a 3D histogram
  X2a_Gen.SetMass(mX2);       X2b_Gen.SetMass(mX2);
  X1a_Gen.SetMass(mX1);       X1b_Gen.SetMass(mX1);
  Za_Gen.SetMass(mZ);
  Za_Gen.SetWidth(wZ);
  Zb_Gen.SetMass(mZ);
  Zb_Gen.SetWidth(wZ);
    
    
    //For now we remove any cuts
  L1a_Gen.SetPtCut(10.);        L1a_Gen.SetEtaCut(2.5);
  L2a_Gen.SetPtCut(10.);        L2a_Gen.SetEtaCut(2.5);
  L1b_Gen.SetPtCut(10.);        L1b_Gen.SetEtaCut(2.5);
  L2b_Gen.SetPtCut(10.);        L2b_Gen.SetEtaCut(2.5);
  
    if(LAB_Reco.InitializeAnalysis())
        g_Log << LogInfo << "...Successfully initialized reconstructed analysis" << LogEnd;
    else
        g_Log << LogError << "...Failed initializing reconstructed analysis" << LogEnd;
  if(LAB_Gen.InitializeAnalysis())
    g_Log << LogInfo << "...Successfully initialized generator analysis" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator analysis" << LogEnd;
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  TreePlot* treePlot = new TreePlot("TreePlot","TreePlot");
 
  treePlot->SetTree(LAB_Gen);
  //treePlot->Draw("GenTree", "Generator Tree", true);
    
    //since there is a correlation between MET and the PT/Eta of the CM frame
    //from 200-1000 GeV (in 100 GeV steps) the correlation depending on the X2 mass
    TFile* input = new TFile("PTEta.root");
    
    string PTEta_histname = "hist_PTvsEtavsMass_";
    int hist_mX2 = mX2;
    PTEta_histname += std::to_string(hist_mX2);
    
    TH3* hist = (TH3*)input->Get(PTEta_histname.c_str());
    Physics physics;
    physics.SetEtaPtMCM(*hist);
    input->Close();
    delete input;
    
    //build the detector
    Detector PUPPI_Detector;
    //Refer to the formula for MET smearing (see Detector.hh)
    PUPPI_Detector.Set_con0_par(14.7467);
    PUPPI_Detector.Set_con1_par(1.68788);
    PUPPI_Detector.Set_con2_par(-2.95569e-07);
    PUPPI_Detector.Set_con0_perp(15.4667);
    PUPPI_Detector.Set_con1_perp(1.41597);
    PUPPI_Detector.Set_con2_perp(-6.37947e-07);
    PUPPI_Detector.Set_sigmaPV(20.0/10000.0); //Primary Vertex resolution
    PUPPI_Detector.Set_sigmaSV(65.0/10000.0); //secondary Vertex resolution
    Vertex PV(0.0,0.0,0.0,0.0); //nominal PV at 0
    Resolution test_Resolution(PUPPI_Detector); //setting up the Resolution "calculator"
    double LAB_Pt;
    double LAB_eta;
    double LAB_M;
    
    //relative uncertainty on the distance between PV and SV
    double sigmaDistance = sqrt((PUPPI_Detector.Get_sigmaPV()*PUPPI_Detector.Get_sigmaPV()+PUPPI_Detector.Get_sigmaSV()*PUPPI_Detector.Get_sigmaSV()));
    
    for(int igen = 0; igen < Ngen; igen++){
        if(igen%((std::max(Ngen,10))/10) == 0)
            g_Log << LogInfo << "Generating event " << igen << " of " << Ngen << LogEnd;
        
        // generate event
        LAB_Gen.ClearEvent();                           // clear the gen tree
        //set momentum based upon the mass
        physics.GetEtaPtMCM(LAB_eta,LAB_Pt,LAB_M);
        //Fix the momentum by hand
        
        LAB_Gen.SetTransverseMomentum(LAB_Pt);
        LAB_Gen.SetLongitudinalMomentum(LAB_Pt*TMath::SinH(LAB_eta));
        //X2X2_Gen.SetMass(LAB_M); //Uncomment for "Official" Plots
        //LAB_Gen.InitializeAnalysis(); //Uncomment for "Official" Plots
        
        LAB_Gen.AnalyzeEvent();                         // generate a new event
        TLorentzVector Pa = X2a_Gen.GetFourVector();
        TLorentzVector Pb = X2b_Gen.GetFourVector();
        TLorentzVector sys = Pa+Pb;
        TLorentzVector Ia  = X1a_Gen.GetFourVector();
        TLorentzVector Ib  = X1b_Gen.GetFourVector();
        TLorentzVector I   = Ia+Ib;
        
        for(int m = 0; m < Nctau; m++)
        {
            LAB_Gen.InitializeAnalysis(); //Comment for "Official" Plots
            
            double ToFa = physics.Get_ToF(ctau[m], Pa);
            double ToFb = physics.Get_ToF(ctau[m], Pb);
            
            if(ToFa <= 0.0)
            {
                igen--;
                continue;
            }
            if(ToFb <= 0.0)
            {
                igen--;
                continue;
            }
            Vertex SVa = physics.Get_SV(ToFa,Pa);
            Vertex SVb = physics.Get_SV(ToFb,Pb);
            Vertex Smeared_PV = PUPPI_Detector.Smear_PV(PV);
            Vertex Smeared_SVa = PUPPI_Detector.Smear_SV(SVa);
            Vertex Smeared_SVb = PUPPI_Detector.Smear_SV(SVb);
            TLorentzVector L1a_RECO = PUPPI_Detector.Smear_Muon(L1a_Gen.GetFourVector());
            TLorentzVector L1b_RECO = PUPPI_Detector.Smear_Muon(L1b_Gen.GetFourVector());
            TLorentzVector L2a_RECO = PUPPI_Detector.Smear_Muon(L2a_Gen.GetFourVector());
            TLorentzVector L2b_RECO = PUPPI_Detector.Smear_Muon(L2b_Gen.GetFourVector());
            TLorentzVector L1a_Gent = L1a_Gen.GetFourVector();
            L1a_Gent.SetZ(0.0);
            TLorentzVector L2a_Gent = L2a_Gen.GetFourVector();
            L2a_Gent.SetZ(0.0);
            TLorentzVector Ia_Gent = Ia;
            Ia_Gent.SetZ(0.0);
            TLorentzVector L1a_RECOt = PUPPI_Detector.Smear_Muon(L1a_Gent);
            TLorentzVector L2a_RECOt = PUPPI_Detector.Smear_Muon(L2a_Gent);
            TLorentzVector L1b_Gent = L1b_Gen.GetFourVector();
            L1b_Gent.SetZ(0.0);
            TLorentzVector L2b_Gent = L2b_Gen.GetFourVector();
            L2b_Gent.SetZ(0.0);
            TLorentzVector Ib_Gent = Ib;
            Ib_Gent.SetZ(0.0);
            TLorentzVector L1b_RECOt = PUPPI_Detector.Smear_Muon(L1b_Gent);
            TLorentzVector L2b_RECOt = PUPPI_Detector.Smear_Muon(L2b_Gent);
            TLorentzVector Va = L1a_RECOt+L2a_RECOt;
            TLorentzVector Vb = L1b_RECOt+L2b_RECOt;
            
             for(int k = 0; k < NsigmaT; k++){
                 PUPPI_Detector.Set_sigmaT((sigmaT[k]/1000.)/sqrt(2.));
                 
                 double Smeared_ToFa = PUPPI_Detector.Smear_ToF(ToFa);
                 double Smeared_ToFb = PUPPI_Detector.Smear_ToF(ToFb);
                 TVector3 Smeared_vBetaa = PUPPI_Detector.Smear_Beta(Smeared_PV,Smeared_SVa);
                 TVector3 Smeared_vBetab = PUPPI_Detector.Smear_Beta(Smeared_PV,Smeared_SVb);
                 if(Smeared_vBetaa.Mag() >= 1.) { continue; }
                 if(Smeared_vBetab.Mag() >= 1.) { continue; }
                 
                 double Da = 30.*Smeared_ToFa*Smeared_vBetaa.Mag();
                 double Db = 30.*Smeared_ToFb*Smeared_vBetab.Mag();
                 
                 //require significant displacement in space and time
                 if(fabs(Smeared_ToFa) < 2.*PUPPI_Detector.Get_sigmaT() || fabs(Smeared_ToFb) < 2.*PUPPI_Detector.Get_sigmaT() || fabs(Da) < 2.*sigmaDistance || fabs(Db) < 2.*sigmaDistance) { continue; }
                 
                 if(sigmaT[k] < 31. && sigmaT[k] > 29.)
                 {
                     for(int h = 0; h < NsigmaMET; h++)
                     {
                         PUPPI_Detector.Set_Sigma_Par(sys, sigmaMET[h]);
                         PUPPI_Detector.Set_Sigma_Perp(sys, sigmaMET[h]);
                         
                         //The smearing begins
                         TVector3 I_Vect = I.Vect();
                         TVector3 MET_RECO_PUPPI = PUPPI_Detector.Smear_MET(I_Vect);
                         MET_RECO_PUPPI.SetZ(0.0);
                         double MXa2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
                         double MXb2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);
                         
                         //Angle Analysis
                         TLorentzVector PX2a;
                         PX2a.SetPxPyPzE(0.0,0.0,0.0,MXa2);
                         TLorentzVector PX2b;
                         PX2b.SetPxPyPzE(0.0,0.0,0.0,MXb2);
                         PX2a.Boost(Smeared_vBetaa);
                         PX2b.Boost(Smeared_vBetab);
                         
                         //RECO Tree
                         LAB_Reco.ClearEvent();
                         L1a_Reco.SetLabFrameFourVector(L1a_RECO);
                         L1b_Reco.SetLabFrameFourVector(L1b_RECO);
                         L2a_Reco.SetLabFrameFourVector(L2a_RECO);
                         L2b_Reco.SetLabFrameFourVector(L2b_RECO);
                         X1a_Reco.SetLabFrameFourVector(PX2a-L1a_RECO-L2a_RECO);
                         X1b_Reco.SetLabFrameFourVector(PX2b-L1b_RECO-L2b_RECO);
                         
                         LAB_Reco.AnalyzeEvent();
                         
                         double CosX2a = X2a_Reco.GetCosDecayAngle();
                         double CosX2b = X2b_Reco.GetCosDecayAngle();
                         
                         if(fabs(CosX2a) > 0.9 || fabs(CosX2b) > 0.9){ continue; }
                         vect_met_acp.at(m).at(h)++;
                     }
                 }
                 PUPPI_Detector.Set_Sigma_Par(sys, 1.);
                 PUPPI_Detector.Set_Sigma_Perp(sys, 1.);
                 
                 //The smearing begins
                 TVector3 I_Vect = I.Vect();
                 TVector3 MET_RECO_PUPPI = PUPPI_Detector.Smear_MET(I_Vect);
                 MET_RECO_PUPPI.SetZ(0.0);
                 double MXa2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
                 double MXb2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);
                 
                 //Angle Analysis
                 TLorentzVector PX2a;
                 PX2a.SetPxPyPzE(0.0,0.0,0.0,MXa2);
                 TLorentzVector PX2b;
                 PX2b.SetPxPyPzE(0.0,0.0,0.0,MXb2);
                 PX2a.Boost(Smeared_vBetaa);
                 PX2b.Boost(Smeared_vBetab);
                 
                 //RECO Tree
                 LAB_Reco.ClearEvent();
                 L1a_Reco.SetLabFrameFourVector(L1a_RECO);
                 L1b_Reco.SetLabFrameFourVector(L1b_RECO);
                 L2a_Reco.SetLabFrameFourVector(L2a_RECO);
                 L2b_Reco.SetLabFrameFourVector(L2b_RECO);
                 X1a_Reco.SetLabFrameFourVector(PX2a-L1a_RECO-L2a_RECO);
                 X1b_Reco.SetLabFrameFourVector(PX2b-L1b_RECO-L2b_RECO);
                 
                 LAB_Reco.AnalyzeEvent();
                 
                 double CosX2a = X2a_Reco.GetCosDecayAngle();
                 double CosX2b = X2b_Reco.GetCosDecayAngle();
                 
                 if(fabs(CosX2a) > 0.9 || fabs(CosX2b) > 0.9){ continue; }
                 vect_timing_acp.at(m).at(k)++;
             }
        }
    }
    end = gSystem->Now();
    g_Log << LogInfo << "Time to Generate " << Ngen*Nctau*NsigmaT << " Events: " << (end-start)/1000.0 << " seconds" << LogEnd;
    for(int j = 0; j < Nctau; j++) { for(int l = 0; l < NsigmaT; l++){ vect_graph_Efficiency_SigmaT.at(j)->SetPoint(l,sigmaT[l],100.*(vect_timing_acp.at(j).at(l))/Ngen); }
        for(int l = 0; l < NsigmaMET; l++){ vect_graph_Efficiency_SigmaMET.at(j)->SetPoint(l,sigmaMET[l],100.*(vect_met_acp.at(j).at(l))/Ngen); } }
    TFile fout(output_name.c_str(),"RECREATE");
    vector<string> leg_text_Efficiency_SigmaT;
    for(int j = 0; j < Nctau; j++){leg_text_Efficiency_SigmaT.push_back("c#tau "+std::to_string(int(ctau.at(j))));}
    Draw_Graphs(fout, vect_graph_Efficiency_SigmaT, leg_text_Efficiency_SigmaT, "Reconstruction Efficiency [%]", "#sigma_{t} [ps]", "Efficiency_Timing");
    vector<string> leg_text_Efficiency_SigmaMET;
    for(int j = 0; j < Nctau; j++){leg_text_Efficiency_SigmaMET.push_back("c#tau "+std::to_string(int(ctau.at(j))));}
    Draw_Graphs(fout, vect_graph_Efficiency_SigmaMET, leg_text_Efficiency_SigmaMET, "Reconstruction Efficiency [%]", "#sigma_{MET} [%]", "Efficiency_MET");
  fout.Close();
  g_Log << LogInfo << "Finished" << LogEnd;
  g_Log << LogInfo << "Time to Process " << Ngen*Nctau << " Events: " << (Long64_t(gSystem->Now())-end)/1000.0 << " seconds" << LogEnd;
}

void Draw_Graphs(TFile& fout, vector<TGraph*>& vect_graph, const vector<string>& leg_text, const string& YaxisText, const string& XaxisText, const string& plotName)
{
    if(!(fout.IsOpen())) {cout << "Output File Not Open..." << endl;}
    TLegend* leg = new TLegend(0.15,0.6,0.403,0.91);
    vector<TLegendEntry*> vect_leg_entry;
    {for(int i = 0; i<int(leg_text.size()); i++) { TLegendEntry* leg_entry = leg->AddEntry(vect_graph.at(i),leg_text.at(i).c_str(),"P"); }}
    
    TCanvas* canvas_graph = new TCanvas(("canvas_graph"+plotName).c_str(),"canvas_graph",750,500);
    canvas_graph->SetGridx();
    canvas_graph->SetGridy();
    canvas_graph->cd();
    
    TMultiGraph* mg = new TMultiGraph();
    for(int i = 0; i < int(vect_graph.size()); i++) { mg->Add(vect_graph.at(i)); }
    vect_graph[0]->SetMarkerStyle(22);
    vect_graph[0]->SetMarkerColor(kBlue);
    if(vect_graph.size() > 1){
    vect_graph[1]->SetMarkerStyle(22);
    vect_graph[1]->SetMarkerColor(kRed);}
    if(vect_graph.size() > 2){
    vect_graph[2]->SetMarkerStyle(22);
    vect_graph[2]->SetMarkerColor(kGreen+2);}
    if(vect_graph.size() == 4){
        vect_graph[3]->SetMarkerStyle(22);
        vect_graph[3]->SetMarkerColor(kMagenta);}
    if(vect_graph.size() == 6){
        vect_graph[3]->SetMarkerStyle(32);
        vect_graph[3]->SetMarkerColor(kBlue);
        vect_graph[4]->SetMarkerStyle(32);
        vect_graph[4]->SetMarkerColor(kRed);
        vect_graph[5]->SetMarkerStyle(32);
        vect_graph[5]->SetMarkerColor(kGreen+2);}
    mg->Draw("AP");
    mg->GetYaxis()->SetTitleOffset(1.05);
    mg->GetYaxis()->SetTitleSize(.04);
    mg->GetXaxis()->SetTitleSize(.04);
    mg->GetYaxis()->SetLabelSize(.04);
    mg->GetXaxis()->SetLabelSize(.04);
    mg->GetYaxis()->SetTitle(YaxisText.c_str());
    mg->GetXaxis()->SetTitle(XaxisText.c_str());
    leg->Draw("SAMES");
    canvas_graph->SaveAs((plotName+".pdf").c_str());
    canvas_graph->Write();
}
