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

#include "../work/Detector.hh"
#include "../work/Physics.hh"
#include "../work/Resolution.hh"
#include "../work/Bonus.h"
#include <TGraph.h>
#include <TSystem.h>

using namespace RestFrames;

void Efficiency_BKG_ZZ(std::string output_name =
			      "Efficiency_BKG_ZZ.root"){
    setMyStyle();
    Long64_t start = gSystem->Now();
    Long64_t end = 0.;
    
    //setting masses and widths
    double mX2 = 400.0;
    double mX1 = 200.0;
    double mZ = 91.19;
    double wZ = 2.50;
    
    vector<double> displacement;
    
    displacement.push_back(0.);
    displacement.push_back(1.);
    displacement.push_back(2.);
    
    int Ndisplacement = displacement.size();
    vector<double> sigmaT;
    
    sigmaT.push_back(0.);
    sigmaT.push_back(10.);
    sigmaT.push_back(30.);
    sigmaT.push_back(50.);
    sigmaT.push_back(100.);
    sigmaT.push_back(300.);
    
    int NsigmaT = sigmaT.size();
    
    //Number of events
    int Ngen = 10000;
    
    bool points = false;
    
    vector<vector<int>> vect_timing_displacement(Ndisplacement, vector<int>(NsigmaT,0));
    vector<vector<int>> vect_timing_displacement_decayangle(Ndisplacement, vector<int>(NsigmaT,0));
    
    vector<TGraph*> vect_graph_timing_displacement;
    vector<TGraph*> vect_graph_timing_displacement_decayangle;
    
    for(int j = 0; j < Ndisplacement; j++)
    {
        TGraph* graph_Efficiency_displacement = new TGraph(NsigmaT);
        graph_Efficiency_displacement->SetName(("graph_displacement"+std::to_string(j)).c_str());
        vect_graph_timing_displacement.push_back(graph_Efficiency_displacement);
        TGraph* graph_Efficiency_displacement_decayangle = new TGraph(NsigmaT);
        graph_Efficiency_displacement_decayangle->SetName(("graph_displacement_decayangle"+std::to_string(j)).c_str());
        vect_graph_timing_displacement_decayangle.push_back(graph_Efficiency_displacement_decayangle);
    }
    
    g_Log << LogInfo << "Initializing generator frames and tree..." << LogEnd;
    
    //Setting up the tree: x2x2->X1X1 + ZZ->4l
    ppLabGenFrame     LAB_Gen("LAB_Gen","LAB");
    DecayGenFrame     ZZ_Gen("ZZ_Gen","Z_{a} Z_{b}");
    ResonanceGenFrame Za_Gen("Za_Gen","Z_{a}");
    ResonanceGenFrame Zb_Gen("Zb_Gen","Z_{b}");
    VisibleGenFrame   L1a_Gen("L1a_Gen","#it{l}_{1a}");
    VisibleGenFrame   L2a_Gen("L2a_Gen","#it{l}_{2a}");
    VisibleGenFrame   L1b_Gen("L1b_Gen","#it{l}_{1b}");
    VisibleGenFrame   L2b_Gen("L2b_Gen","#it{l}_{2b}");
    
    LabRecoFrame       LAB_Reco("LAB_Reco","LAB");
    DecayRecoFrame     ZZ_Reco("ZZ_Reco","Z_{a} Z_{b}");
    DecayRecoFrame     Za_Reco("Za_Reco","Z_{a}");
    DecayRecoFrame     Zb_Reco("Zb_Reco","Z_{b}");
    VisibleRecoFrame   L1a_Reco("L1a_Reco","#it{l}_{1a}");
    VisibleRecoFrame   L2a_Reco("L2a_Reco","#it{l}_{2a}");
    VisibleRecoFrame   L1b_Reco("L1b_Reco","#it{l}_{1b}");
    VisibleRecoFrame   L2b_Reco("L2b_Reco","#it{l}_{2b}");
    
    //
    LAB_Gen.SetChildFrame(ZZ_Gen);
    ZZ_Gen.AddChildFrame(Za_Gen);
    ZZ_Gen.AddChildFrame(Zb_Gen);
    Za_Gen.AddChildFrame(L1a_Gen);
    Za_Gen.AddChildFrame(L2a_Gen);
    Zb_Gen.AddChildFrame(L1b_Gen);
    Zb_Gen.AddChildFrame(L2b_Gen);
    
    LAB_Reco.SetChildFrame(ZZ_Reco);
    ZZ_Reco.AddChildFrame(Za_Reco);
    ZZ_Reco.AddChildFrame(Zb_Reco);
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
    
    ZZ_Gen.SetVariableMass(); //Reset this down below from a 3D histogram
    Za_Gen.SetMass(mZ);
    Za_Gen.SetWidth(wZ);
    Zb_Gen.SetMass(mZ);
    Zb_Gen.SetWidth(wZ);
    
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
    TFile* input = new TFile("PTEta_BKG_ZZ.root");
    
    string PTEta_histname = "hist_PTvsEtavsMass_BKG";
    
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
    Vertex SV_true(0.0,0.0,0.0,0.0); //nominal PV at 0
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
        TLorentzVector Pa  = Za_Gen.GetFourVector();
        TLorentzVector Pb  = Zb_Gen.GetFourVector();
        TLorentzVector sys = Pa+Pb;
        TLorentzVector Ia  = L1a_Gen.GetFourVector()+L2a_Gen.GetFourVector();
        TLorentzVector Ib  = L1b_Gen.GetFourVector()+L2b_Gen.GetFourVector();
        TLorentzVector I   = Ia+Ib;
        I.SetZ(0.0);
        
        for(int m = 0; m < Ndisplacement; m++)
        {
            LAB_Gen.InitializeAnalysis(); //Comment for "Official" Plots
            
            TLorentzVector L1a_RECO = PUPPI_Detector.Smear_Muon(L1a_Gen.GetFourVector());
            TLorentzVector L1b_RECO = PUPPI_Detector.Smear_Muon(L1b_Gen.GetFourVector());
            TLorentzVector L2a_RECO = PUPPI_Detector.Smear_Muon(L2a_Gen.GetFourVector());
            TLorentzVector L2b_RECO = PUPPI_Detector.Smear_Muon(L2b_Gen.GetFourVector());
            TLorentzVector Va = L1a_RECO+L2a_RECO;
            TLorentzVector Vb = L1b_RECO+L2b_RECO;
            
             for(int k = 0; k < NsigmaT; k++){
                 bool displaced = false;
                 PUPPI_Detector.Set_sigmaT((sigmaT[k]/1000.)/sqrt(2.));
                 Vertex Smeared_PV = PUPPI_Detector.Smear_PV(PV);
                 Vertex Smeared_SVa = PUPPI_Detector.Smear_SV(SV_true);
                 Vertex Smeared_SVb = PUPPI_Detector.Smear_SV(SV_true);
                 TVector3 Smeared_vBetaa = PUPPI_Detector.Get_Beta(Smeared_PV,Smeared_SVa);
                 TVector3 Smeared_vBetab = PUPPI_Detector.Get_Beta(Smeared_PV,Smeared_SVb);
                 
                 double Da = 30.*Smeared_SVa.GetTPos()*Smeared_vBetaa.Mag();
                 double Db = 30.*Smeared_SVa.GetTPos()*Smeared_vBetab.Mag();
                 
                 //require significant displacement in space and time
                 if(fabs(Smeared_SVa.GetTPos()) < displacement.at(m)*PUPPI_Detector.Get_sigmaT() || fabs(Smeared_SVb.GetTPos()) < displacement.at(m)*PUPPI_Detector.Get_sigmaT() || fabs(Da) < displacement.at(m)*sigmaDistance || fabs(Db) < displacement.at(m)*sigmaDistance || Smeared_vBetaa.Mag() >= 1. || Smeared_vBetab.Mag() >= 1.) {vect_timing_displacement.at(m).at(k)++; displaced = true;}
                 
                 PUPPI_Detector.Set_Sigma_Par(sys, 1.);
                 PUPPI_Detector.Set_Sigma_Perp(sys, 1.);
                 
                 //The smearing begins
                 TVector3 MET = L1a_RECO.Vect() + L1b_RECO.Vect() + L2a_RECO.Vect() + L2b_RECO.Vect();
                 
                 double MXa2 = test_Resolution.Mass_Parents2(MET,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
                 double MXb2 = test_Resolution.Mass_Parents2(MET,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);
                 
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
                 LAB_Reco.AnalyzeEvent();
                 
                 double CosX2a = Za_Reco.GetCosDecayAngle();
                 double CosX2b = Zb_Reco.GetCosDecayAngle();
                 if((fabs(CosX2a) > 0.8 || fabs(CosX2b) > 0.8) || displaced){ vect_timing_displacement_decayangle.at(m).at(k)++; }
             }
        }
    }
    end = gSystem->Now();
    g_Log << LogInfo << "Time to Generate " << Ngen << " Events: " << (end-start)/1000.0 << " seconds" << LogEnd;
    for(int j = 0; j < Ndisplacement; j++) {
    for(int l = 0; l < NsigmaT; l++) {
        vect_graph_timing_displacement.at(j)->SetPoint(l,sigmaT[l],100.*(Ngen-vect_timing_displacement.at(j).at(l))/Ngen);
        vect_graph_timing_displacement_decayangle.at(j)->SetPoint(l,sigmaT[l],100.*(Ngen-vect_timing_displacement_decayangle.at(j).at(l))/Ngen);
    }
    }
    TFile fout(output_name.c_str(),"RECREATE");
    vector<string> leg_text_displacement;
    for(int j = 0; j < Ndisplacement; j++){leg_text_displacement.push_back(std::to_string(int(displacement.at(j)))+"#sigma_{Displacement}");}
    Draw_Graphs(fout, vect_graph_timing_displacement, leg_text_displacement, "Reconstruction Efficiency #sigma [%]", "#sigma_{t} [ps]", "Efficiency_Sigma", points);
    Draw_Graphs(fout, vect_graph_timing_displacement_decayangle, leg_text_displacement, "Reconstruction Efficiency #sigma [%] & |Cos(#theta)| < 0.8 [%]", "#sigma_{t} [ps]", "Efficiency_Sigma_DecayAngle", points);
    vector<TGraph*> vect_graph_displacement;
    for(int k = 0; k < int(vect_graph_timing_displacement.size()); k++) {vect_graph_displacement.push_back(vect_graph_timing_displacement.at(k));}
    for(int k = 0; k < int(vect_graph_timing_displacement_decayangle.size()); k++) {vect_graph_displacement.push_back(vect_graph_timing_displacement_decayangle.at(k));}
    vector<string> leg_text_all;
    for(int j = 0; j < int(vect_graph_timing_displacement.size()); j++){leg_text_all.push_back(std::to_string(int(displacement.at(j)))+"#sigma_{Displacement}");}
    for(int j = 0; j < int(vect_graph_timing_displacement_decayangle.size()); j++){leg_text_all.push_back(std::to_string(int(displacement.at(j)))+"#sigma_{Displacement} & |Cos(#theta_{#tilde{#chi}^{ 0}_{2a}})| < 0.8");}
    Draw_Graphs(fout, vect_graph_displacement, leg_text_all, "Reconstruction Efficiency [%]", "#sigma_{t} [ps]", "Efficiency", points);
  fout.Close();
  g_Log << LogInfo << "Finished" << LogEnd;
  g_Log << LogInfo << "Time to Process " << Ngen*Ndisplacement << " Events: " << (Long64_t(gSystem->Now())-end)/1000.0 << " seconds" << LogEnd;
}
