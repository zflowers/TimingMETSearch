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

void SVz_ZZ_to_4l(std::string output_name =
			      "output_SV_ZZ_to_4l.root"){

    Long64_t start = gSystem->Now();
    Long64_t end = 0.;
    
    //setting masses and widths
    double mZ = 91.19;
    double wZ = 2.50;
    
    vector<double> SV; //resolution on SVs in microns
    
    //SV.push_back(30.);
    SV.push_back(65.);
    //SV.push_back(100.);
    
    int NSV = SV.size();
    
    //Number of events
    int Ngen = 1000000;
    int Entries = Ngen;
    double displacement_cut = 3.;
    
    g_Log << LogInfo << "Initializing generator frames and tree..." << LogEnd;
    
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
  treePlot->Draw("GenTree", "Generator Tree", true);
  
    // Declare observables for histogram booking
    HistPlot* histPlot = new HistPlot("Plots",
                                      std::string("Za(#it{l}#it{l}) Zb(#it{l}#it{l})"));
    
    histPlot->SetRebin(1);
    
    RFList<const HistPlotCategory> cat_list_SV;
    char smassX2[200];
    string sSV = "#delta_{SV} = ";
    for(int m = 0; m < NSV; m++){
        
        char snameSV[200], scatSV[50];
        sprintf(scatSV, "#delta_{SV}_%d", m);
        sprintf(snameSV, "%.1f #mu m", SV[m]);
        
        cat_list_SV += histPlot->GetNewCategory(scatSV, sSV+std::string(snameSV));
    }
    
    const HistPlotVar& EZa = histPlot->GetNewVar("EZ", "E_{Z}^{#tilde{#chi}_{2}^{0}}", 0., 500., "[GeV]");
    const HistPlotVar& MXa2 = histPlot->GetNewVar("MXa2", "M(#tilde{#chi}_{2a}^{0})", 0., 2000., "[GeV]");
    const HistPlotVar& MIa2 = histPlot->GetNewVar("MIa2", "M(#tilde{#chi}_{1a}^{0})", 0., 2000., "[GeV]");
    const HistPlotVar& MXb2 = histPlot->GetNewVar("MXb2", "M(#tilde{#chi}_{2b}^{0})", 0., 2000., "[GeV]");
    
    histPlot->AddPlot(EZa, cat_list_SV);
    histPlot->AddPlot(MXa2, cat_list_SV);
    histPlot->AddPlot(MIa2, cat_list_SV);
    
    //get the kinematics of the ZZ system
    //use 200 GeV for now
    TFile* input = new TFile("PTEta_BKG.root");
    TH3* hist = (TH3*)input->Get("hist_PTvsEtavsMass_BKG");
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
    PUPPI_Detector.Set_sigmaT(0.03/sqrt(2.)); //timing resolution
    PUPPI_Detector.Set_sigmaPV(20.0/10000.0); //Primary Vertex resolution
    //PUPPI_Detector.Set_sigmaSV(65.0/10000.0); //secondary Vertex resolution
    Vertex PV(0.0,0.0,0.0,0.0); //nominal PV at 0
    Vertex SV_true(0.0,0.0,0.0,0.0); //nominal PV at 0
    Resolution test_Resolution(PUPPI_Detector); //setting up the Resolution "calculator"
    double LAB_Pt;
    double LAB_eta;
    double LAB_M;
    
    //relative uncertainty on the distance between PV and SV
    //double sigmaDistance = sqrt((PUPPI_Detector.Get_sigmaPV()*PUPPI_Detector.Get_sigmaPV()+PUPPI_Detector.Get_sigmaSV()*PUPPI_Detector.Get_sigmaSV()));
    
    //for checking generator efficiency
    int gen_events = 0;
    int acp_events = 0;
    
  for(int m = 0; m < NSV; m++){
    g_Log << LogInfo << "Generating BKG events for ";
    g_Log << "SV = " << SV[m] << LogEnd;
    
    LAB_Gen.InitializeAnalysis(); //Comment for "Official" Plots
    PUPPI_Detector.Set_sigmaSV(SV[m]/10000.0); //secondary Vertex resolution
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
        gen_events++;
        TLorentzVector Pa  = Za_Gen.GetFourVector();
        TLorentzVector Pb  = Zb_Gen.GetFourVector();
        TLorentzVector sys = Pa+Pb;
        TLorentzVector Ia  = L1a_Gen.GetFourVector()+L2a_Gen.GetFourVector();
        TLorentzVector Ib  = L1b_Gen.GetFourVector()+L2b_Gen.GetFourVector();
        TLorentzVector I   = Ia+Ib;
        I.SetZ(0.0);
        
        double MET_Mag_Resolution = PUPPI_Detector.Get_Sigma_Par(sys);
        double MET_Dir_Resolution = PUPPI_Detector.Get_Sigma_Perp(sys);
        
        TLorentzVector L1a_RECO = PUPPI_Detector.Smear_Muon(L1a_Gen.GetFourVector());
        TLorentzVector L1b_RECO = PUPPI_Detector.Smear_Muon(L1b_Gen.GetFourVector());
        TLorentzVector L2a_RECO = PUPPI_Detector.Smear_Muon(L2a_Gen.GetFourVector());
        TLorentzVector L2b_RECO = PUPPI_Detector.Smear_Muon(L2b_Gen.GetFourVector());
        
        //TVector3 MET = PUPPI_Detector.Smear_MET(I.Vect());
        TVector3 MET = L1a_RECO.Vect() + L1b_RECO.Vect() + L2a_RECO.Vect() + L2b_RECO.Vect();
        
        Vertex Smeared_PV = PUPPI_Detector.Smear_PV(PV);
        Vertex Smeared_SVa = PUPPI_Detector.Smear_SV(SV_true);
        Vertex Smeared_SVb = PUPPI_Detector.Smear_SV(SV_true);
        TVector3 Smeared_vBetaa = PUPPI_Detector.Get_Beta(Smeared_PV,Smeared_SVa);
        TVector3 Smeared_vBetab = PUPPI_Detector.Get_Beta(Smeared_PV,Smeared_SVb);
        
        double Da = 30.*Smeared_SVa.GetTPos()*Smeared_vBetaa.Mag();
        double Db = 30.*Smeared_SVa.GetTPos()*Smeared_vBetab.Mag();
        
        //if(fabs(Smeared_SVa.GetTPos()) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Smeared_SVb.GetTPos()) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Da) < displacement_cut*sigmaDistance || fabs(Db) < displacement_cut*sigmaDistance || Smeared_vBetaa.Mag() >= 1. || Smeared_vBetab.Mag() >= 1.) { igen--; continue;}
        if(Smeared_vBetaa.Mag() >= 1. || Smeared_vBetab.Mag() >= 1.) { igen--; continue;}
        
        double Sigma_Beta_Mag = sqrt((1.0/(Smeared_SVa.GetTPos()*Smeared_SVa.GetTPos()))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetaa.Mag()*Smeared_vBetaa.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));
        
        TLorentzVector vZa = L1a_RECO + L2a_RECO;
        vZa.Boost(-Smeared_vBetaa);
        EZa = vZa.E();
        
        TLorentzVector Va = L1a_RECO+L2a_RECO;
        TLorentzVector Vb = L1b_RECO+L2b_RECO;
        
        MXa2 = test_Resolution.Mass_Parents2(MET,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
        MXb2 = test_Resolution.Mass_Parents2(MET,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);
        
        MIa2 = test_Resolution.Mass_Invisible2(MET,Va,Vb,Smeared_vBetaa,Smeared_vBetab);
        if(MIa2 <= 0.) {igen--; continue;}
        //Reco Analysis
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
        
        histPlot->Fill(cat_list_SV[m]);
        acp_events++;
    }
    LAB_Gen.PrintGeneratorEfficiency();
  }
    g_Log << LogInfo << "Generated a Total of " << gen_events << " Events " << LogEnd;
    g_Log << LogInfo << acp_events << " passed selection requirements " << LogEnd;
    g_Log << LogInfo << "Efficiency: " << 100.0*acp_events/gen_events << "%" << LogEnd;
    end = gSystem->Now();
    g_Log << LogInfo << "Time to Generate " << Ngen*NSV << " Events: " << (end-start)/1000.0 << " seconds" << LogEnd;
    g_Log << LogInfo << "Processing " << Ngen*NSV << " Events" << LogEnd;
    histPlot->Draw(false);
    
    TFile fout(output_name.c_str(),"RECREATE");
  fout.Close();
  histPlot->WriteOutput(output_name);
  histPlot->WriteHist(output_name);
  treePlot->WriteOutput(output_name);
  g_Log << LogInfo << "Finished" << LogEnd;
  g_Log << LogInfo << "Time to Process " << Ngen*NSV << " Events: " << (Long64_t(gSystem->Now())-end)/1000.0 << " seconds" << LogEnd;
}
