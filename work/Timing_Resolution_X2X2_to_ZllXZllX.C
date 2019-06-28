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
#include <TSystem.h>
#include "Bonus.h"

using namespace RestFrames;

void Timing_Resolution_X2X2_to_ZllXZllX(std::string output_name =
			      "output_Timing_Resolution_X2X2_to_ZallXZbllX.root"){

    Long64_t start = gSystem->Now();

    //setting masses and widths
    double mX2 = 700.0;
    double mX1 = 500.0;
    double mZ = 91.19;
    double wZ = 2.50;
    double ctau = 50.;

    vector<double> sigmaT;

    sigmaT.push_back(30.);
    //sigmaT.push_back(50.);
    //sigmaT.push_back(300.);
    
    //Number of events
    int Ngen = 100000;
    
    int NsigmaT = sigmaT.size();

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
    HistPlot* histPlot = new HistPlot("Plots",
                                      std::string("#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}") +
                                      "#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}"+
                                      "Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");

    histPlot->SetRebin(1);

    RFList<const HistPlotCategory> cat_list;
    string ssigmaT = "#sigma_{T} = ";
    for(int m = 0; m < NsigmaT; m++){
        char sname[200], scat[50];
        sprintf(scat, "sigmaT_%d", m);
        sprintf(sname, "%.1f ps", sigmaT[m]);
        cat_list += histPlot->GetNewCategory(scat, ssigmaT+std::string(sname));
    }

    //setting up all the variables that could be potentially plotted
    const HistPlotVar& Pull_Mass_Parent = histPlot->GetNewVar("Pull_Mass_Parent","Pull of M(#tilde{#chi}_{2}^{0})",-5.0,5.0,"");
    const HistPlotVar& Pull_E_Z = histPlot->GetNewVar("Pull_E_Z","Pull of #delta E_{Za}^{#tilde{#chi}_{2a}^{0}}",-5.0,5.0,"");
    const HistPlotVar& Pull_E_L = histPlot->GetNewVar("Pull_E_L","Pull of #delta E_{L1a}",-5.0,5.0,"");
    const HistPlotVar& Pull_D = histPlot->GetNewVar("Pull_D","Pull of Distance",-5.0,5.0,"");
    const HistPlotVar& Pull_T = histPlot->GetNewVar("Pull_T","Pull of Time",-5.0,5.0,"");
    const HistPlotVar& Pull_Beta_Mag = histPlot->GetNewVar("Pull_Beta_Mag","Pull of |#vec{#Beta}|",-5.0,5.0,"");
    const HistPlotVar& Pull_MET = histPlot->GetNewVar("Pull MET","Pull of MET",-5.0,5.0,"");
    const HistPlotVar& ToFaL = histPlot->GetNewVar("ToFaL", "log_{10} ToF(#tilde{#chi}_{2}^{0})", -3.5, 2.5, "[ns]");
    const HistPlotVar& ToFbL = histPlot->GetNewVar("ToFbL", "log_{10} ToF(#tilde{#chi}_{2}^{0})", -3.5, 2.5, "[ns]");
    const HistPlotVar& DaL = histPlot->GetNewVar("DaL", "log_{10} D(#tilde{#chi}_{2}^{0})", -2.5, 4.5, "[cm]");
    const HistPlotVar& DbL = histPlot->GetNewVar("DbL", "log_{10} D(#tilde{#chi}_{2}^{0})", -2.5, 4.5, "[cm]");
    const HistPlotVar& ToFa = histPlot->GetNewVar("ToFa","ToF(#tilde{#chi}_{2}^{0})",0.0,30.0,"[ns]");
    const HistPlotVar& ToFb = histPlot->GetNewVar("ToFb","ToF(#tilde{#chi}_{2}^{0})",0.0,30.0,"[ns]");
    const HistPlotVar& Da = histPlot->GetNewVar("Da", "D(#tilde{#chi}_{2}^{0})", 0.0, 0.05, "[cm]");
    const HistPlotVar& Db = histPlot->GetNewVar("Db", "D(#tilde{#chi}_{2}^{0})", -15.0, 15.0, "[cm]");
    const HistPlotVar& betaa = histPlot->GetNewVar("betaa", "#beta(#tilde{#chi}_{2}^{0})", 0., 1.);
    const HistPlotVar& betab = histPlot->GetNewVar("betab", "#beta(#tilde{#chi}_{2}^{0})", 0., 1.);
    const HistPlotVar& EZa = histPlot->GetNewVar("EZa", "E_{Za}^{#tilde{#chi}_{2a}^{0}}", 0., 800., "[GeV]");
    const HistPlotVar& Pull_Mass_Invisible = histPlot->GetNewVar("Pull_Mass_Inv","Pull of M(#tilde{#chi}_{1a}^{0})",-5.0,5.0,"");
    const HistPlotVar& MIa = histPlot->GetNewVar("MIa", "M(#tilde{#chi}_{1a}^{0})", 0., 1000., "[GeV]");
    const HistPlotVar& MIa2 = histPlot->GetNewVar("MIa2", "M(#tilde{#chi}_{1a}^{0})", 0., 1000., "[GeV]");
    const HistPlotVar& Delta_MIa2 = histPlot->GetNewVar("Delta_MIa2", "#DeltaM(#tilde{#chi}_{1a}^{0})", -1000.0, 1000.0, "");
    const HistPlotVar& Pull_MIa2 = histPlot->GetNewVar("Pull_MIa2", "Pull of M(#tilde{#chi}_{1a}^{0})", -5.0, 5.0, "");
    const HistPlotVar& Pull_MXa2 = histPlot->GetNewVar("Pull_MXa2", "Pull of M(#tilde{#chi}_{2a}^{0})", -5.0, 5.0, "");
    const HistPlotVar& Pull_MXb2 = histPlot->GetNewVar("Pull_MXb2", "Pull of M(#tilde{#chi}_{2b}^{0})", -5.0, 5.0, "");
    const HistPlotVar& Pull_Vis = histPlot->GetNewVar("Pull_Vis", "Pull of Vis", -5.0, 5.0, "");
    const HistPlotVar& MX2X2 = histPlot->GetNewVar("MX2X2", "M(#tilde{#chi}_{2}^{0})(#tilde{#chi}_{2}^{0})", 0., 3000., "[GeV]");
    const HistPlotVar& MXa2 = histPlot->GetNewVar("MXa2", "M(#tilde{#chi}_{2a}^{0})", 0., 2000., "[GeV]");
    const HistPlotVar& MXb2 = histPlot->GetNewVar("MXb2", "M(#tilde{#chi}_{2b}^{0})", 0., 2000., "[GeV]");
    const HistPlotVar& MXa = histPlot->GetNewVar("MXa", "M(#tilde{#chi}_{2}^{0})", 0., 2000., "[GeV]");
    const HistPlotVar& CosX2a = histPlot->GetNewVar("CosX2a", "Cos_X2a", -1, 1., "");
    const HistPlotVar& CosX2b = histPlot->GetNewVar("CosX2b", "Cos_X2b", -1, 1., "");
    const HistPlotVar& CosX2a_Gen = histPlot->GetNewVar("CosX2a_Gen", "Cos_X2a_Gen", -1, 1., "");
    const HistPlotVar& CosX2b_Gen = histPlot->GetNewVar("CosX2b_Gen", "Cos_X2b_Gen", -1, 1., "");
    const HistPlotVar& X2X2_Frame_X2a = histPlot->GetNewVar("X2X2_Frame_X2a","X2X2_Frame_X2a",-0.1,6.5,"");
    const HistPlotVar& X2X2_Frame_X2b = histPlot->GetNewVar("X2X2_Frame_X2b","X2X2_Frame_X2b",-0.1,6.5,"");
    const HistPlotVar& X2a_Frame_X2b = histPlot->GetNewVar("X2a_Frame_X2b","X2a_Frame_X2b",-0.1,6.5,"");
    const HistPlotVar& X2X2_Frame_X2a_Gen = histPlot->GetNewVar("X2X2_Frame_X2a_Gen","X2X2_Frame_X2a_Gen",-0.1,6.5,"");
    const HistPlotVar& X2X2_Frame_X2b_Gen = histPlot->GetNewVar("X2X2_Frame_X2b_Gen","X2X2_Frame_X2b_Gen",-0.1,6.5,"");
    const HistPlotVar& X2a_Frame_X2b_Gen = histPlot->GetNewVar("X2a_Frame_X2b_Gen","X2a_Frame_X2b_Gen",-0.1,6.5,"");
    const HistPlotVar& tReco_tTrue = histPlot->GetNewVar("tReco_tTrue","#theta^{R}_{X^{0}_{2a}} - #theta^{G}_{X^{0}_{2a}}",-1.,1.,"");
    //Various Pulls of mass
    const HistPlotVar& Pull_MXa2L = histPlot->GetNewVar("Pull_MXa2L", "Pull of M(#tilde{#chi}_{2a}^{0})L", -5.0, 5.0, ""); //turn off lepton
    const HistPlotVar& Pull_MXa2B = histPlot->GetNewVar("Pull_MXa2B", "Pull of M(#tilde{#chi}_{2a}^{0})B", -5.0, 5.0, ""); //turn off beta
    const HistPlotVar& Pull_MXa2D = histPlot->GetNewVar("Pull_MXa2D", "Pull of M(#tilde{#chi}_{2a}^{0})D", -5.0, 5.0, ""); //turn off MET Direction

    //comment in/out whatever plots are interesting
    //histPlot->AddPlot(Pull_Mass_Parent, cat_list); //need ~1 TeV MX2 & ~300 MX1, to look good
    //histPlot->AddPlot(Pull_MET, cat_list);
    //histPlot->AddPlot(ToFaL, cat_list);
    //histPlot->AddPlot(DaL, cat_list);
    //histPlot->AddPlot(ToFbL, cat_list);
    //histPlot->AddPlot(DbL, cat_list);
    //histPlot->AddPlot(ToFa, cat_list);
    //histPlot->AddPlot(Da, cat_list);
    //histPlot->AddPlot(ToFb, cat_list);
    //histPlot->AddPlot(Db, cat_list);
    //histPlot->AddPlot(betaa, cat_list);
    //histPlot->AddPlot(EZa, cat_list);
    //histPlot->AddPlot(Pull_E_Z, cat_list);
    //histPlot->AddPlot(Pull_E_L, cat_list);
    //histPlot->AddPlot(Pull_D, cat_list);
    //histPlot->AddPlot(Pull_T, cat_list);
    //histPlot->AddPlot(Pull_Beta_Mag, cat_list);
    //histPlot->AddPlot(Pull_Mass_Invisible, cat_list); //need ~1 TeV MX2 & ~500 MX1, to look ok
    //histPlot->AddPlot(MIa, cat_list);
    //histPlot->AddPlot(MXa, cat_list);
    histPlot->AddPlot(MXa2, cat_list);
    histPlot->AddPlot(MIa2, cat_list);
    //histPlot->AddPlot(Delta_MIa2, cat_list);
    histPlot->AddPlot(Pull_MIa2, cat_list);
    //histPlot->AddPlot(MXa2, MXb2, cat_list);
    histPlot->AddPlot(Pull_MXa2, cat_list);
    //histPlot->AddPlot(Pull_MXa2, Pull_MXb2, cat_list);
    //histPlot->AddPlot(Pull_MXa2L, cat_list);
    //histPlot->AddPlot(Pull_MXa2B, cat_list);
    //histPlot->AddPlot(Pull_MXa2D, cat_list);
    //histPlot->AddPlot(Pull_Vis, cat_list);
    //histPlot->AddPlot(CosX2a, cat_list);
    //histPlot->AddPlot(CosX2b, cat_list);
    //histPlot->AddPlot(Pull_Beta_Mag, CosX2a, cat_list);
    //histPlot->AddPlot(Pull_MXa2, CosX2a, cat_list);
    //histPlot->AddPlot(CosX2a_Gen, cat_list);
    //histPlot->AddPlot(CosX2b_Gen, cat_list);
    //histPlot->AddPlot(CosX2a_Gen, CosX2a, cat_list);
    //histPlot->AddPlot(CosX2b_Gen, CosX2b, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2a, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2b, cat_list);
    //histPlot->AddPlot(X2a_Frame_X2b, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2a, Pull_MXa2, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2a, Pull_MXb2, cat_list);
    //histPlot->AddPlot(X2a_Frame_X2b, Pull_MXa2, cat_list);
    //histPlot->AddPlot(X2a_Frame_X2b, Pull_MXb2, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2b, Pull_MXb2, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2a, CosX2a, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2a, X2X2_Frame_X2b, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2a_Gen, cat_list);
    //histPlot->AddPlot(X2X2_Frame_X2b_Gen, cat_list);
    //histPlot->AddPlot(X2a_Frame_X2b_Gen, cat_list);
    //histPlot->AddPlot(MX2X2, cat_list);
    //histPlot->AddPlot(tReco_tTrue, cat_list);

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
    double LAB_Pt;
    double LAB_eta;
    double LAB_M;

    TVector3 Zhat(0.,0.,1.);

    double log_10=log(10.);

    //relative uncertainty on the distance between PV and SV
    double sigmaDistance = sqrt((PUPPI_Detector.Get_sigmaPV()*PUPPI_Detector.Get_sigmaPV()+PUPPI_Detector.Get_sigmaSV()*PUPPI_Detector.Get_sigmaSV()));

    //for checking generator efficiency
    int gen_events = 0;
    int acp_events = 0;

    for(int m = 0; m < NsigmaT; m++){
    //if(flag && m > 0) continue;
    g_Log << LogInfo << "Generating events for ";
    g_Log << "mX2 = " << mX2 << ", ";
    g_Log << "mX1 = " << mX1 << ", ";
    g_Log << "sigmaT = " << sigmaT[m] << LogEnd;

      PUPPI_Detector.Set_sigmaT((sigmaT[m]/1000.0)/sqrt(2.));
        Resolution test_Resolution(PUPPI_Detector); //setting up the Resolution "calculator"
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
        I.SetZ(0.0);

        //The smearing begins
        double MET_Mag_Resolution = PUPPI_Detector.Get_Sigma_Par(sys);
        double MET_Dir_Resolution = PUPPI_Detector.Get_Sigma_Perp(sys);

        TVector3 I_Vect = I.Vect();
        TVector3 MET_RECO_PUPPI = PUPPI_Detector.Smear_MET(I_Vect);
        MET_RECO_PUPPI.SetZ(0.0);
        TVector3 Ia_RECO = PUPPI_Detector.Smear_MET(Ia.Vect());
        Ia_RECO.SetZ(0.0);
        TVector3 Ib_RECO = PUPPI_Detector.Smear_MET(Ib.Vect());
        Ib_RECO.SetZ(0.0);

        TLorentzVector L1a_RECO = PUPPI_Detector.Smear_Muon(L1a_Gen.GetFourVector());
        TLorentzVector L1b_RECO = PUPPI_Detector.Smear_Muon(L1b_Gen.GetFourVector());
        TLorentzVector L2a_RECO = PUPPI_Detector.Smear_Muon(L2a_Gen.GetFourVector());
        TLorentzVector L2b_RECO = PUPPI_Detector.Smear_Muon(L2b_Gen.GetFourVector());

        TVector3 vBetaaGen = Pa.BoostVector();
        TVector3 vBetabGen = Pb.BoostVector();

        ToFa = physics.Get_ToF(ctau, Pa);
        ToFb = physics.Get_ToF(ctau, Pb);

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
        Vertex Smeared_PV = PUPPI_Detector.Smear_PV(PV);
        Vertex Smeared_SVa = PUPPI_Detector.Smear_SV(SVa);
        Vertex Smeared_SVb = PUPPI_Detector.Smear_SV(SVb);
        TVector3 Smeared_vBetaa = PUPPI_Detector.Get_Beta(Smeared_PV,Smeared_SVa);
        TVector3 Smeared_vBetab = PUPPI_Detector.Get_Beta(Smeared_PV,Smeared_SVb);
        TLorentzVector L1a_Gent = L1a_Gen.GetFourVector();
        L1a_Gent.SetZ(0.0);
        TLorentzVector L2a_Gent = L2a_Gen.GetFourVector();
        L2a_Gent.SetZ(0.0);
        TLorentzVector L1a_RECOt = PUPPI_Detector.Smear_Muon(L1a_Gent);
        TLorentzVector L2a_RECOt = PUPPI_Detector.Smear_Muon(L2a_Gent);
        TLorentzVector L1b_Gent = L1b_Gen.GetFourVector();
        L1b_Gent.SetZ(0.0);
        TLorentzVector L2b_Gent = L2b_Gen.GetFourVector();
        L2b_Gent.SetZ(0.0);
        TLorentzVector L1b_RECOt = PUPPI_Detector.Smear_Muon(L1b_Gent);
        TLorentzVector L2b_RECOt = PUPPI_Detector.Smear_Muon(L2b_Gent);

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

        Da = 30.*ToFa*vBetaaGen.Mag();;
        Db = 30.*ToFb*vBetabGen.Mag();;

        /*
         if(betaa < VelocityUncertainty || betab < VelocityUncertainty) //require significant displacement
         {
         igen--;
         continue;
         }
         */



        TVector3 vPta = (L1a_RECO+L2a_RECO).Vect()+Ia_RECO;
        TVector3 vPtaGen = (L1a_Gen.GetFourVector()+L2a_Gen.GetFourVector()+Ia).Vect();
        vPtaGen.SetZ(0.0);
        TVector3 vPtb = (L1b_RECO+L2b_RECO).Vect()+Ib_RECO;
        TVector3 vPt  = (L1a_RECO+L2a_RECO+L1b_RECO+L2b_RECO).Vect()+MET_RECO_PUPPI;
        vPta.SetZ(0.);
        vPtb.SetZ(0.);
        vPt.SetZ(0.);
        TVector3 vBetaaT_Gen = vBetaaGen;
        vBetaaT_Gen.SetZ(0.0);

        double Sigma_Angle = test_Resolution.GetAngleError(vBetaaGen,vPtaGen,Smeared_vBetaa,vPta);
        double Sigma_E_Lepton1 = PUPPI_Detector.GetMuonResolution(L1a_Gen.GetFourVector());
        double Sigma_E_Lepton2 = PUPPI_Detector.GetMuonResolution(L2a_Gen.GetFourVector());
        double Sigma_Angle1 = test_Resolution.GetAngleError(vBetaaGen,L1a_Gen.GetFourVector().Vect(),Smeared_vBetaa,L1a_RECO.Vect());
        double Sigma_Angle2 = test_Resolution.GetAngleError(vBetaaGen,L2a_Gen.GetFourVector().Vect(),Smeared_vBetaa,L2a_RECO.Vect());

        double Sigma_Beta_Mag = sqrt((1.0/(Smeared_ToFa*Smeared_ToFa))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetaa.Mag()*Smeared_vBetaa.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));

        double Sigma_Beta_Magb = sqrt((1.0/(Smeared_ToFb*Smeared_ToFb))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetab.Mag()*Smeared_vBetab.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));


        ToFaL = log(ToFa)/log_10;
        DaL   = log(Da)/log_10;
        ToFbL = log(ToFb)/log_10;
        DbL   = log(Db)/log_10;

        TVector3 Da_Smear(Smeared_PV.GetXPos()-Smeared_SVa.GetXPos(),Smeared_PV.GetYPos()-Smeared_SVa.GetYPos(),Smeared_PV.GetZPos()-Smeared_SVa.GetZPos());
        TVector3 Da_True(PV.GetXPos()-SVa.GetXPos(),PV.GetYPos()-SVa.GetYPos(),PV.GetZPos()-SVa.GetZPos());

        Pull_D = (Da_Smear.Mag()-Da_True.Mag())/sigmaDistance;
        Pull_T = (ToFa-Smeared_ToFa)/PUPPI_Detector.Get_sigmaT();
        Pull_T = (SVa.GetTPos()-Smeared_SVa.GetTPos())/PUPPI_Detector.Get_sigmaT();
        //check lepton smearing in transverse plane
        Pull_E_L = (L1a_Gent.Pt() - L1a_RECOt.Pt())/(PUPPI_Detector.GetMuonResolution(L1a_Gent)*L1a_RECOt.Pt()); //check lepton smearing in the transverse plane

        TLorentzVector vZa = L1a_RECO + L2a_RECO;
        vZa.Boost(-Smeared_vBetaa);
        TLorentzVector vZb = L1b_RECO + L2b_RECO;
        vZb.Boost(-Smeared_vBetab);

        EZa = vZa.E();
        TLorentzVector vZaGen = L1a_Gen.GetFourVector() + L2a_Gen.GetFourVector();
        vZaGen.Boost(-vBetaaGen);


        Pull_Beta_Mag = (vBetaaGen.Mag() - Smeared_vBetaa.Mag())/Sigma_Beta_Mag;

        Pull_MET = (MET_RECO_PUPPI.Mag()-I.Vect().Mag())/MET_Mag_Resolution;

        Pull_E_Z = (EZa - vZaGen.E())/test_Resolution.Energy_Z_Parent_Resolution(vBetaaGen,L1a_RECO,L2a_RECO,L1a_RECO.E()*PUPPI_Detector.GetMuonResolution(L1a_Gen.GetFourVector()),L2a_RECO.E()*PUPPI_Detector.GetMuonResolution(L2a_Gen.GetFourVector()),Sigma_Angle1,Sigma_Angle2,Sigma_Beta_Mag); //see Resolution.hh for calculation

        double MP_Resolution = test_Resolution.Mass_Parent_Resolution(Smeared_vBetaa,Ia_RECO,L1a_RECOt,L2a_RECOt,MET_Mag_Resolution,L1a_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L1a_Gent),L2a_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L2a_Gent),Sigma_Beta_Mag);

        double MXa_Gen = test_Resolution.Mass_Parent(vPtaGen,vBetaaGen);
        MXa = test_Resolution.Mass_Parent(vPta,Smeared_vBetaa);
        Pull_Mass_Parent = (MXa-MXa_Gen)/MP_Resolution;
        //one LLP

        double Mass_Invisible_Resolution = test_Resolution.Mass_Invisible_Resolution(Smeared_vBetaa,Ia_RECO,L1a_RECO,L2a_RECO,MET_Mag_Resolution,Sigma_Beta_Mag);

        double MI_Gen = sqrt(MXa_Gen*MXa_Gen-2.*MXa_Gen*vZaGen.E()+((L1a_Gen.GetFourVector()+L2a_Gen.GetFourVector()).M2()));
        MIa = sqrt(MXa*MXa-2.*MXa*EZa+((L1a_RECO+L2a_RECO).M2()));
        Pull_Mass_Invisible = (MI_Gen-MIa)/Mass_Invisible_Resolution;

        
        
        TLorentzVector Va = L1a_RECO+L2a_RECO;
        TLorentzVector Vb = L1b_RECO+L2b_RECO;
        
        double Sigma_Vis = test_Resolution.Visible_Resolution(L1a_RECOt.Vect(),L2a_RECOt.Vect(),L1b_RECOt.Vect(),L2b_RECOt.Vect(),L1a_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L1a_Gent),L2a_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L2a_Gent),L1b_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L1b_Gent),L2b_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L2b_Gent));

        Pull_Vis = ((L1a_RECOt.Vect()+L2a_RECOt.Vect()+L1b_RECOt.Vect()+L2b_RECOt.Vect()).Mag()-(L1a_Gent.Vect()+L2a_Gent.Vect()+L1b_Gent.Vect()+L2b_Gent.Vect()).Mag())/Sigma_Vis;

        double MPa_Gen = test_Resolution.Mass_Parents2(I_Vect,L1a_Gent.Vect()+L2a_Gent.Vect()+L1b_Gent.Vect()+L2b_Gent.Vect(),vBetaaGen,vBetabGen);
        double MPb_Gen = test_Resolution.Mass_Parents2(I_Vect,L1a_Gent.Vect()+L2a_Gent.Vect()+L1b_Gent.Vect()+L2b_Gent.Vect(),vBetabGen,vBetaaGen);
        MXa2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
        MXb2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);

        double MXa2_ResolutionL = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,MET_Mag_Resolution,MET_Dir_Resolution,0.);
        double MXa2_ResolutionB = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,0.,MET_Mag_Resolution,MET_Dir_Resolution,0.);
        double MXa2_ResolutionD = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,0.,MET_Mag_Resolution,0.,0.);

        Pull_MXa2L = (MPa_Gen-MXa2)/MXa2_ResolutionL;
        Pull_MXa2B = (MPa_Gen-MXa2)/MXa2_ResolutionB;
        Pull_MXa2D = (MPa_Gen-MXa2)/MXa2_ResolutionD;

        double MXb2_Resolution = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa,Sigma_Beta_Magb,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);
        double MXa2_Resolution = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);

        Pull_MXa2 = (MPa_Gen-MXa2)/MXa2_Resolution;
        Pull_MXb2 = (MPb_Gen-MXb2)/MXb2_Resolution;

        //Two LLPs
        double MIa2_Gen = test_Resolution.Mass_Invisible2(MPa_Gen,vZaGen.E(),(L1a_Gen.GetFourVector()+L2a_Gen.GetFourVector()).M());
        MIa2 = test_Resolution.Mass_Invisible2(MXa2,EZa,(L1a_RECO+L2a_RECO).M());
        double MIa2_Res = test_Resolution.Mass_Invisible_Resolution2(MET_RECO_PUPPI,Va,Vb,Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,MET_Mag_Resolution,MET_Dir_Resolution);
        Delta_MIa2 = MIa2-MIa2_Gen;
        Pull_MIa2 = Delta_MIa2/MIa2_Gen;

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

        CosX2a = X2a_Reco.GetCosDecayAngle();
        CosX2b = X2b_Reco.GetCosDecayAngle();
        CosX2a_Gen = X2a_Gen.GetCosDecayAngle();
        CosX2b_Gen = X2b_Gen.GetCosDecayAngle();
        X2X2_Frame_X2a = X2X2_Reco.GetDeltaPhiDecayPlanes(X2a_Reco);
        X2X2_Frame_X2b = X2X2_Reco.GetDeltaPhiDecayPlanes(X2b_Reco);
        X2a_Frame_X2b = X2a_Reco.GetDeltaPhiDecayPlanes(X2b_Reco);
        tReco_tTrue = TMath::ACos(CosX2a) - TMath::ACos(CosX2a_Gen);

        X2X2_Frame_X2a_Gen = X2X2_Gen.GetDeltaPhiDecayPlanes(X2a_Gen);
        X2X2_Frame_X2b_Gen = X2X2_Gen.GetDeltaPhiDecayPlanes(X2b_Gen);
        X2a_Frame_X2b_Gen = X2a_Gen.GetDeltaPhiDecayPlanes(X2b_Gen);
        if(MIa2 > 0.)
        {
           histPlot->Fill(cat_list[m]);
           acp_events++;
        }
    }
    LAB_Gen.PrintGeneratorEfficiency();
  }

    histPlot->Draw();


    TFile fout(output_name.c_str(),"RECREATE");

    fout.Close();
    histPlot->WriteOutput(output_name);
    histPlot->WriteHist(output_name);
    treePlot->WriteOutput(output_name);

  g_Log << LogInfo << "Finished" << LogEnd;
    g_Log << LogInfo << "Generated a Total of " << gen_events << " Events " << LogEnd;
    g_Log << LogInfo << acp_events << " passed selection requirements " << LogEnd;
    g_Log << LogInfo << "Efficiency: " << 100.0*acp_events/gen_events << "%" << LogEnd;
    Long64_t end = gSystem->Now();
    g_Log << LogInfo << "Time to Generate " << Ngen*NsigmaT << " Events: " << (end-start)/1000.0 << " seconds" << LogEnd;
}
