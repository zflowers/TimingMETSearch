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

void LSP_testing(std::string output_name =
			      "output_LSP_testing.root"){
    setMyStyle();
    Long64_t start = gSystem->Now();
    Long64_t end = 0.;
    
    //setting masses and widths
    double mX2 = 1000.0;
    double mX1 = 100.0;
    double mZ = 91.19;
    double wZ = 2.50;
    
    vector<double> ctau;
    
    ctau.push_back(25.);
    ctau.push_back(10.);
    ctau.push_back(5.);
    //ctau.push_back(1.);
    
    int Nctau = ctau.size();
    vector<double> sigmaT;
    vector<double> sigmaMET;
    
    for(double i = 10.; i <= 350.; i+=100.)
    {
        sigmaMET.push_back(i);
        sigmaT.push_back(i);
    }
    
    //sigmaT.push_back(30.);
    int NsigmaT = sigmaT.size();
    int NsigmaMET = sigmaMET.size();
    bool ctau_flag = true; //set to false to turn off anything related to looping over ctau
    bool timing_flag = false; //set to false to turn off anything related to looping over sigmat
    bool MET_flag = false;
    bool LSP_flag = false; //set to false to turn off regenerating the event for non-real LSP reco masses
    bool LSP_Calc_flag = false; //leave this alone
    
    //Number of events
    int Ngen = 100000;
    
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
  
    // Declare observables for histogram booking
    HistPlot* histPlot_Cut = new HistPlot("Plots_Cut",
                                      std::string("#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}") +
                                      "#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}"+
                                      "Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");
    
    HistPlot* histPlot_NoCut = new HistPlot("Plots_NoCut",
                                      std::string("#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}") +
                                      "#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}"+
                                      "Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");
    
    histPlot_Cut->SetRebin(1);
    histPlot_NoCut->SetRebin(1);
    
    RFList<const HistPlotCategory> cat_list_ctau_cut;
    char smassX2_cut[200];
    string sctau_cut = "c#tau = ";
    for(int m = 0; m < Nctau; m++){
        char snamectau[200], scatctau[50];
        sprintf(scatctau, "ctau_%d", m);
        sprintf(snamectau, "%.1f cm", ctau[m]);
        cat_list_ctau_cut += histPlot_Cut->GetNewCategory(scatctau, sctau_cut+std::string(snamectau));
    }
    
    RFList<const HistPlotCategory> cat_list_ctau_nocut;
    char smassX2_nocut[200];
    string sctau_nocut = "c#tau = ";
    for(int m = 0; m < Nctau; m++){
        char snamectau[200], scatctau[50];
        sprintf(scatctau, "ctau_%d", m);
        sprintf(snamectau, "%.1f cm", ctau[m]);
        cat_list_ctau_nocut += histPlot_NoCut->GetNewCategory(scatctau, sctau_nocut+std::string(snamectau));
    }
    
    //setting up all the variables that could be potentially plotted
    const HistPlotVar& MET = histPlot_Cut->GetNewVar("MET","MET",0.,1000.0,"");
    const HistPlotVar& ToFa = histPlot_Cut->GetNewVar("ToFa","ToF(#tilde{#chi}_{2}^{0})",0.0,30.0,"[ns]");
    const HistPlotVar& ToFb = histPlot_Cut->GetNewVar("ToFb","ToF(#tilde{#chi}_{2}^{0})",0.0,30.0,"[ns]");
    const HistPlotVar& Da = histPlot_Cut->GetNewVar("Da", "D(#tilde{#chi}_{2}^{0})", 0.0, 0.05, "[cm]");
    const HistPlotVar& Db = histPlot_Cut->GetNewVar("Db", "D(#tilde{#chi}_{2}^{0})", -15.0, 15.0, "[cm]");
    const HistPlotVar& EZa_Cut = histPlot_Cut->GetNewVar("EZa", "E_{Za}^{#tilde{#chi}_{2a}^{0}}", 0., 800., "[GeV]");
    const HistPlotVar& EZa_NoCut = histPlot_Cut->GetNewVar("EZa", "E_{Za}^{#tilde{#chi}_{2a}^{0}}", 0., 800., "[GeV]");
    const HistPlotVar& Pull_Mass_Invisible = histPlot_Cut->GetNewVar("Pull_Mass_Inv","Pull of M(#tilde{#chi}_{1a}^{0})",-5.0,5.0,"");
    const HistPlotVar& MIa = histPlot_Cut->GetNewVar("MIa", "M(#tilde{#chi}_{1a}^{0})", 0., 1000., "[GeV]");
    const HistPlotVar& MIa2 = histPlot_Cut->GetNewVar("MIa2", "M(#tilde{#chi}_{1a}^{0})", -50., 300., "[GeV]");
    const HistPlotVar& MIb2 = histPlot_Cut->GetNewVar("MIb2", "M(#tilde{#chi}_{1b}^{0})", -50., 300., "[GeV]");
    const HistPlotVar& MIa3 = histPlot_Cut->GetNewVar("MIa2", "M(#tilde{#chi}_{1a}^{0}) miniRJR", -700., 600., "[GeV]");
    const HistPlotVar& MX2X2 = histPlot_Cut->GetNewVar("MX2X2", "M(#tilde{#chi}_{2}^{0})(#tilde{#chi}_{2}^{0})", 0., 3000., "[GeV]");
    const HistPlotVar& MV = histPlot_Cut->GetNewVar("MV", "M(#it{la})", 0., 150., "[GeV]");
    const HistPlotVar& MXa2_Cut = histPlot_Cut->GetNewVar("MXa2", "M(#tilde{#chi}_{2a}^{0})", 0., 1200., "[GeV]");
    const HistPlotVar& MXa2_NoCut = histPlot_NoCut->GetNewVar("MXa2", "M(#tilde{#chi}_{2a}^{0})", 0., 1200., "[GeV]");
    
    
    //histPlot_Cut->AddPlot(MIa2, cat_list_ctau_cut);
    //histPlot_Cut->AddPlot(MIa3, cat_list_ctau_cut);
    histPlot_Cut->AddPlot(MXa2_Cut, cat_list_ctau_cut);
    histPlot_Cut->AddPlot(EZa_Cut, cat_list_ctau_cut);
    histPlot_NoCut->AddPlot(MXa2_NoCut, cat_list_ctau_nocut);
    histPlot_NoCut->AddPlot(EZa_NoCut, cat_list_ctau_nocut);
    /*
    histPlot_Cut->AddPlot(MIa2, MV, cat_list_ctau_cut);
    histPlot_Cut->AddPlot(MIa2, MET, cat_list_ctau_cut);
    histPlot_Cut->AddPlot(MIa2, MXa2, cat_list_ctau_cut);
    histPlot_Cut->AddPlot(MIa2, EZa, cat_list_ctau_cut);
    histPlot_Cut->AddPlot(MIa2, MIb2, cat_list_ctau_cut);
    */
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
    PUPPI_Detector.Set_sigmaT((sigmaT[0]/1000.)/sqrt(2.)); //timing resolution
    //PUPPI_Detector.Set_sigmaT(0.); //timing resolution
    PUPPI_Detector.Set_sigmaPV(20.0/10000.0); //Primary Vertex resolution
    PUPPI_Detector.Set_sigmaSV(65.0/10000.0); //secondary Vertex resolution
    Vertex PV(0.0,0.0,0.0,0.0); //nominal PV at 0
    Resolution test_Resolution(PUPPI_Detector); //setting up the Resolution "calculator"
    double LAB_Pt;
    double LAB_eta;
    double LAB_M;
    
    TVector3 Zhat(0.,0.,1.);
    TVector3 Ref_Vect(1.,1.,1.);
    
    double Epa;
    double Epb;
    double log_10=log(10.);
    
    //relative uncertainty on the distance between PV and SV
    double sigmaDistance = sqrt((PUPPI_Detector.Get_sigmaPV()*PUPPI_Detector.Get_sigmaPV()+PUPPI_Detector.Get_sigmaSV()*PUPPI_Detector.Get_sigmaSV()));
    
    //for checking generator efficiency
    int gen_events = 0;
    int acp_events = 0;
    
    
    if(ctau_flag){
  for(int m = 0; m < Nctau; m++){
    g_Log << LogInfo << "Generating events for ";
    g_Log << "mX2 = " << mX2 << ", ";
    g_Log << "mX1 = " << mX1 << ", ";
    g_Log << "ctau = " << ctau[m] << LogEnd;
    
    LAB_Gen.InitializeAnalysis(); //Comment for "Official" Plots
      
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
        gen_events++;
        MX2X2 = X2X2_Gen.GetMass();
        TLorentzVector Pa = X2a_Gen.GetFourVector();
        TLorentzVector Pb = X2b_Gen.GetFourVector();
        TLorentzVector sys = Pa+Pb;
        TLorentzVector Ia  = X1a_Gen.GetFourVector();
        TLorentzVector Ib  = X1b_Gen.GetFourVector();
        TLorentzVector I   = Ia+Ib;
        
        PUPPI_Detector.Set_Sigma_Par(sys, 1.);
        PUPPI_Detector.Set_Sigma_Perp(sys, 1.);
        //PUPPI_Detector.Set_Sigma_Par(sys, .25);
        //PUPPI_Detector.Set_Sigma_Perp(sys, .25);
        
        //The smearing begins
        TVector3 I_Vect = I.Vect();
        TVector3 MET_RECO_PUPPI = I_Vect;// PUPPI_Detector.Smear_MET(I_Vect);
        MET_RECO_PUPPI.SetZ(0.0);
        
        TLorentzVector L1a_RECO = L1a_Gen.GetFourVector();//PUPPI_Detector.Smear_Muon(L1a_Gen.GetFourVector());
        TLorentzVector L1b_RECO = L1b_Gen.GetFourVector();//PUPPI_Detector.Smear_Muon(L1b_Gen.GetFourVector());
        TLorentzVector L2a_RECO = L2a_Gen.GetFourVector();//PUPPI_Detector.Smear_Muon(L2a_Gen.GetFourVector());
        TLorentzVector L2b_RECO = L2b_Gen.GetFourVector();//PUPPI_Detector.Smear_Muon(L2b_Gen.GetFourVector());
        
        TVector3 vBetaaGen = Pa.BoostVector();
        TVector3 vBetabGen = Pb.BoostVector();
        
        ToFa = physics.Get_ToF(ctau[m], Pa);
        ToFb = physics.Get_ToF(ctau[m], Pb);
        
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
        
        double Smeared_ToFa = PUPPI_Detector.Smear_ToF(ToFa);
        double Smeared_ToFb = PUPPI_Detector.Smear_ToF(ToFb);
        Vertex SVa = physics.Get_SV(ToFa,Pa);
        Vertex SVb = physics.Get_SV(ToFb,Pb);
        Vertex Smeared_PV = PV;//PUPPI_Detector.Smear_PV(PV);
        Vertex Smeared_SVa = SVa;//PUPPI_Detector.Smear_SV(SVa);
        Vertex Smeared_SVb = SVb;//PUPPI_Detector.Smear_SV(SVb);
        TVector3 Smeared_vBetaa = PUPPI_Detector.Smear_Beta(Smeared_PV,Smeared_SVa);
        TVector3 Smeared_vBetab = PUPPI_Detector.Smear_Beta(Smeared_PV,Smeared_SVb);
        TLorentzVector L1a_Gent = L1a_Gen.GetFourVector();
        L1a_Gent.SetZ(0.0);
        TLorentzVector L2a_Gent = L2a_Gen.GetFourVector();
        L2a_Gent.SetZ(0.0);
        TLorentzVector Ia_Gent = Ia;
        Ia_Gent.SetZ(0.0);
        TLorentzVector L1a_RECOt = L1a_Gent;//PUPPI_Detector.Smear_Muon(L1a_Gent);
        TLorentzVector L2a_RECOt = L2a_Gent;//PUPPI_Detector.Smear_Muon(L2a_Gent);
        TLorentzVector L1b_Gent = L1b_Gen.GetFourVector();
        L1b_Gent.SetZ(0.0);
        TLorentzVector L2b_Gent = L2b_Gen.GetFourVector();
        L2b_Gent.SetZ(0.0);
        TLorentzVector Ib_Gent = Ib;
        Ib_Gent.SetZ(0.0);
        TLorentzVector L1b_RECOt = L1b_Gent;//PUPPI_Detector.Smear_Muon(L1b_Gent);
        TLorentzVector L2b_RECOt = L2b_Gent;//PUPPI_Detector.Smear_Muon(L2b_Gent);
        
        if(Smeared_vBetaa.Mag() >= 1.)
        {
            igen--;
            continue;
        }
        if(Smeared_vBetab.Mag() >= 1.)
        {
            igen--;
            continue;
        }
        
        double Da = 30.*Smeared_ToFa*Smeared_vBetaa.Mag();
        double Db = 30.*Smeared_ToFb*Smeared_vBetab.Mag();
        
        if(fabs(Smeared_ToFa) < 2.*PUPPI_Detector.Get_sigmaT() || fabs(Smeared_ToFb) < 2.*PUPPI_Detector.Get_sigmaT() || fabs(Da) < 2.*sigmaDistance || fabs(Db) < 2.*sigmaDistance) //require significant displacement in space and time
        {
            igen--;
            continue;
        }
        
        TLorentzVector vZa = L1a_RECO + L2a_RECO;
        vZa.Boost(-Smeared_vBetaa);
        TLorentzVector vZb = L1b_RECO + L2b_RECO;
        vZb.Boost(-Smeared_vBetab);
        
        EZa_Cut = vZa.E();
        EZa_NoCut = vZa.E();
        TLorentzVector vZaGen = L1a_Gen.GetFourVector() + L2a_Gen.GetFourVector();
        vZaGen.Boost(-vBetaaGen);
        
        double EZb = vZb.E();
        
        TLorentzVector Va = L1a_RECOt+L2a_RECOt;
        TLorentzVector Vb = L1b_RECOt+L2b_RECOt;
        
        MXa2_Cut = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
        MXa2_NoCut = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
        double MXb2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);
        
        //Two LLPs
        MV = (L1a_RECO+L2a_RECO).M();
        double MVb = (L1b_RECO+L2b_RECO).M();
        MIa2 = test_Resolution.Mass_Invisible2(MXa2_Cut,EZa_Cut,MV);
        MIb2 = test_Resolution.Mass_Invisible2(MXb2,EZb,MVb);
        if(MIa2 <= 0. && LSP_flag)
        {
            igen--;
            continue;
        }
        
        //Angle Analysis
        TLorentzVector PX2a;
        PX2a.SetPxPyPzE(0.0,0.0,0.0,MXa2_Cut);
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
        MIa3 = X1a_Reco.GetFourVector().M();
        
        double CosX2a = X2a_Reco.GetCosDecayAngle();
        double CosX2b = X2b_Reco.GetCosDecayAngle();
        
        if(fabs(CosX2a) > 0.9 || fabs(CosX2b) > 0.9)
        {
            igen--;
            continue;
        }
        if(MIa2 < 1.e-6 && MIb2 < 1.e-6) {histPlot_Cut->Fill(cat_list_ctau_cut[m]); acp_events++;}
        else {histPlot_NoCut->Fill(cat_list_ctau_nocut[m]);}
        //histPlot_Cut->Fill(cat_list_ctau_cut[m]);
        //acp_events++;
    }
    }
    //LAB_Gen.PrintGeneratorEfficiency();
  }
    g_Log << LogInfo << "Generated a Total of " << gen_events << " Events " << LogEnd;
    g_Log << LogInfo << acp_events << " passed selection requirements " << LogEnd;
    g_Log << LogInfo << "Efficiency: " << 100.0*acp_events/gen_events << "%" << LogEnd;
    end = gSystem->Now();
    g_Log << LogInfo << "Time to Generate " << Ngen*Nctau << " Events: " << (end-start)/1000.0 << " seconds" << LogEnd;
    g_Log << LogInfo << "Processing " << Ngen*Nctau << " Events" << LogEnd;
    histPlot_Cut->Draw();
    histPlot_NoCut->Draw();
    TFile fout(output_name.c_str(),"RECREATE");
  fout.Close();
  histPlot_Cut->WriteOutput(output_name);
  histPlot_Cut->WriteHist(output_name);
  histPlot_NoCut->WriteOutput(output_name);
  histPlot_NoCut->WriteHist(output_name);
  treePlot->WriteOutput(output_name);
  g_Log << LogInfo << "Finished" << LogEnd;
  g_Log << LogInfo << "Time to Process " << Ngen*Nctau << " Events: " << (Long64_t(gSystem->Now())-end)/1000.0 << " seconds" << LogEnd;
}
