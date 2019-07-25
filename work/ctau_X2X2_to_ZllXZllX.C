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

void ctau_X2X2_to_ZllXZllX(std::string output_name =
			      "output_ctau_X2X2_to_ZallXZbllX.root"){

    Long64_t start = gSystem->Now();
    Long64_t end = 0.;
    gStyle->SetOptTitle(0);
    
    //setting masses and widths
    double mX2 = 400.0;
    double mX1 = 350.0;
    double mZ = 91.19;
    double wZ = 2.50;
    
    vector<double> ctau;
    
    ctau.push_back(20.);
    ctau.push_back(10.);
    ctau.push_back(5.);
    
    int Nctau = ctau.size();
    vector<double> sigmaT;
    vector<double> sigmaMET;
    /*
    for(double i = 0.; i <= 330.; i+=30.)
    {
        sigmaMET.push_back(i);
        sigmaT.push_back(i);
    }
    */
    sigmaT.push_back(0.);
    sigmaT.push_back(10.);
    sigmaT.push_back(30.);
    sigmaT.push_back(50.);
    sigmaT.push_back(300.);
    
    sigmaMET.push_back(25.);
    sigmaMET.push_back(75.);
    sigmaMET.push_back(100.);
    sigmaMET.push_back(150.);
    sigmaMET.push_back(300.);
    
    //sigmaT.push_back(30.);
    int NsigmaT = sigmaT.size();
    int NsigmaMET = sigmaMET.size();
    bool timing_flag = false; //set to false to turn off anything related to looping over sigmaT
    bool MET_flag = false; //set to false to turn off anything related to looping over sigmaMET
    bool points = false;
    bool decayangle = false;
    
    //Number of events
    int Ngen = 1000000;
    int Entries = 100000;
    //int Entries = 100000; //for comparing turning effects off
    double displacement_cut = 3.;
    
    int bins_MX2 = 1024;
    double xmin_MX2 = 0.;
    double xmax_MX2 = 1500.;
    int bins_MX1 = 1024;
    double xmin_MX1 = 0.;
    vector<double> xmax_MX1 = { 1200., 1200., 1200. };
    double ana_factor = 500.;
    
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_MET_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Timing_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Measured_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Measured_SigmaT_cos08(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Measured_SigmaT_cos09(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_MET_Measured_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Timing_Measured_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_MET_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_Timing_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_Measured_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_MET_Measured_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_Timing_Measured_SigmaT(Nctau, vector<TH1D*>(NsigmaT));
    
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_MET_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Timing_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Measured_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_MET_Measured_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX2_Timing_Measured_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_MET_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_Timing_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_Measured_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_MET_Measured_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    vector<vector<TH1D*>> vect_hist_Sigma_MX1_Timing_Measured_SigmaMET(Nctau, vector<TH1D*>(NsigmaMET));
    
    vector<TGraph*> vect_graph_Sigma_MX2_SigmaT;
    vector<TGraph*> vect_graph_Sigma_MX2_MET_SigmaT;
    vector<TGraph*> vect_graph_Sigma_MX2_Timing_SigmaT;
    vector<TGraph*> vect_graph_Sigma_MX2_SigmaT_Measured;
    vector<TGraph*> vect_graph_Sigma_MX2_SigmaT_Measured_cos08;
    vector<TGraph*> vect_graph_Sigma_MX2_SigmaT_Measured_cos09;
    vector<TGraph*> vect_graph_Sigma_MX2_MET_SigmaT_Measured;
    vector<TGraph*> vect_graph_Sigma_MX2_Timing_SigmaT_Measured;
    
    vector<TGraph*> vect_graph_Sigma_MX1_SigmaT;
    vector<TGraph*> vect_graph_Sigma_MX1_MET_SigmaT;
    vector<TGraph*> vect_graph_Sigma_MX1_Timing_SigmaT;
    vector<TGraph*> vect_graph_Sigma_MX1_SigmaT_Measured;
    vector<TGraph*> vect_graph_Sigma_MX1_MET_SigmaT_Measured;
    vector<TGraph*> vect_graph_Sigma_MX1_Timing_SigmaT_Measured;
    
    vector<TGraph*> vect_graph_Sigma_MX2_SigmaMET;
    vector<TGraph*> vect_graph_Sigma_MX2_MET_SigmaMET;
    vector<TGraph*> vect_graph_Sigma_MX2_Timing_SigmaMET;
    vector<TGraph*> vect_graph_Sigma_MX2_SigmaMET_Measured;
    vector<TGraph*> vect_graph_Sigma_MX2_MET_SigmaMET_Measured;
    vector<TGraph*> vect_graph_Sigma_MX2_Timing_SigmaMET_Measured;
    
    vector<TGraph*> vect_graph_Sigma_MX1_SigmaMET;
    vector<TGraph*> vect_graph_Sigma_MX1_MET_SigmaMET;
    vector<TGraph*> vect_graph_Sigma_MX1_Timing_SigmaMET;
    vector<TGraph*> vect_graph_Sigma_MX1_SigmaMET_Measured;
    vector<TGraph*> vect_graph_Sigma_MX1_MET_SigmaMET_Measured;
    vector<TGraph*> vect_graph_Sigma_MX1_Timing_SigmaMET_Measured;
    
    for(int j = 0; j < Nctau; j++)
    {
        if(timing_flag){
        TGraph* graph_Sigma_MX2_SigmaT = new TGraph(NsigmaT);
        graph_Sigma_MX2_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_MET_SigmaT = new TGraph(NsigmaT);
        graph_Sigma_MX2_MET_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_Timing_SigmaT = new TGraph(NsigmaT);
        graph_Sigma_MX2_Timing_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_SigmaT_Measured = new TGraph(NsigmaT);
        graph_Sigma_MX2_SigmaT_Measured->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_SigmaT_Measured_cos08 = new TGraph(NsigmaT);
        graph_Sigma_MX2_SigmaT_Measured_cos08->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_SigmaT_Measured_cos09 = new TGraph(NsigmaT);
        graph_Sigma_MX2_SigmaT_Measured_cos09->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_MET_SigmaT_Measured = new TGraph(NsigmaT);
        graph_Sigma_MX2_MET_SigmaT_Measured->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_Timing_SigmaT_Measured = new TGraph(NsigmaT);
        graph_Sigma_MX2_Timing_SigmaT_Measured->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        
        TGraph* graph_Sigma_MX1_SigmaT = new TGraph(NsigmaT);
        graph_Sigma_MX1_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_MET_SigmaT = new TGraph(NsigmaT);
        graph_Sigma_MX1_MET_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_Timing_SigmaT = new TGraph(NsigmaT);
        graph_Sigma_MX1_Timing_SigmaT->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_SigmaT_Measured = new TGraph(NsigmaT);
        graph_Sigma_MX1_SigmaT_Measured->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_MET_SigmaT_Measured = new TGraph(NsigmaT);
        graph_Sigma_MX1_MET_SigmaT_Measured->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_Timing_SigmaT_Measured = new TGraph(NsigmaT);
        graph_Sigma_MX1_Timing_SigmaT_Measured->SetName(("graph_sigmaT"+std::to_string(j)).c_str());
        
        vect_graph_Sigma_MX2_SigmaT.push_back(graph_Sigma_MX2_SigmaT);
        vect_graph_Sigma_MX2_MET_SigmaT.push_back(graph_Sigma_MX2_MET_SigmaT);
        vect_graph_Sigma_MX2_Timing_SigmaT.push_back(graph_Sigma_MX2_Timing_SigmaT);
        vect_graph_Sigma_MX2_SigmaT_Measured.push_back(graph_Sigma_MX2_SigmaT_Measured);
        vect_graph_Sigma_MX2_SigmaT_Measured_cos08.push_back(graph_Sigma_MX2_SigmaT_Measured_cos08);
        vect_graph_Sigma_MX2_SigmaT_Measured_cos09.push_back(graph_Sigma_MX2_SigmaT_Measured_cos09);
        vect_graph_Sigma_MX2_MET_SigmaT_Measured.push_back(graph_Sigma_MX2_MET_SigmaT_Measured);
        vect_graph_Sigma_MX2_Timing_SigmaT_Measured.push_back(graph_Sigma_MX2_Timing_SigmaT_Measured);
        
        vect_graph_Sigma_MX1_SigmaT.push_back(graph_Sigma_MX1_SigmaT);
        vect_graph_Sigma_MX1_MET_SigmaT.push_back(graph_Sigma_MX1_MET_SigmaT);
        vect_graph_Sigma_MX1_Timing_SigmaT.push_back(graph_Sigma_MX1_Timing_SigmaT);
        vect_graph_Sigma_MX1_SigmaT_Measured.push_back(graph_Sigma_MX1_SigmaT_Measured);
        vect_graph_Sigma_MX1_MET_SigmaT_Measured.push_back(graph_Sigma_MX1_MET_SigmaT_Measured);
        vect_graph_Sigma_MX1_Timing_SigmaT_Measured.push_back(graph_Sigma_MX1_Timing_SigmaT_Measured);
        for(int k = 0; k < NsigmaT; k++){
            TH1D* hist_Sigma_MX2_SigmaT = new TH1D(("hist_Sigma_MX2_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2/ana_factor);
            vect_hist_Sigma_MX2_SigmaT.at(j).at(k) = hist_Sigma_MX2_SigmaT;
            TH1D* hist_Sigma_MX2_MET_SigmaT = new TH1D(("hist_Sigma_MX2_MET_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2/ana_factor);
            vect_hist_Sigma_MX2_MET_SigmaT.at(j).at(k) = hist_Sigma_MX2_MET_SigmaT;
            TH1D* hist_Sigma_MX2_Timing_SigmaT = new TH1D(("hist_Sigma_MX2_Timing_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2/ana_factor);
            vect_hist_Sigma_MX2_Timing_SigmaT.at(j).at(k) = hist_Sigma_MX2_Timing_SigmaT;
            TH1D* hist_Sigma_MX2_Measured_SigmaT = new TH1D(("hist_Sigma_MX2_Measured_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_Measured_SigmaT.at(j).at(k) = hist_Sigma_MX2_Measured_SigmaT;
            TH1D* hist_Sigma_MX2_Measured_SigmaT_cos08 = new TH1D(("hist_Sigma_MX2_Measured_SigmaT_cos08"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_Measured_SigmaT_cos08.at(j).at(k) = hist_Sigma_MX2_Measured_SigmaT_cos08;
            TH1D* hist_Sigma_MX2_Measured_SigmaT_cos09 = new TH1D(("hist_Sigma_MX2_Measured_SigmaT_cos09"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_Measured_SigmaT_cos09.at(j).at(k) = hist_Sigma_MX2_Measured_SigmaT_cos09;
            TH1D* hist_Sigma_MX2_MET_Measured_SigmaT = new TH1D(("hist_Sigma_MX2_MET_Measured_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_MET_Measured_SigmaT.at(j).at(k) = hist_Sigma_MX2_MET_Measured_SigmaT;
            TH1D* hist_Sigma_MX2_Timing_Measured_SigmaT = new TH1D(("hist_Sigma_MX2_Timing_Measured_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_Timing_Measured_SigmaT.at(j).at(k) = hist_Sigma_MX2_Timing_Measured_SigmaT;
            TH1D* hist_Sigma_MX1_SigmaT = new TH1D(("hist_Sigma_MX1_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]/ana_factor);
            vect_hist_Sigma_MX1_SigmaT.at(j).at(k) = hist_Sigma_MX1_SigmaT;
            TH1D* hist_Sigma_MX1_MET_SigmaT = new TH1D(("hist_Sigma_MX1_MET_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]/ana_factor);
            vect_hist_Sigma_MX1_MET_SigmaT.at(j).at(k) = hist_Sigma_MX1_MET_SigmaT;
            TH1D* hist_Sigma_MX1_Timing_SigmaT = new TH1D(("hist_Sigma_MX1_Timing_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]/ana_factor);
            vect_hist_Sigma_MX1_Timing_SigmaT.at(j).at(k) = hist_Sigma_MX1_Timing_SigmaT;
            TH1D* hist_Sigma_MX1_Measured_SigmaT = new TH1D(("hist_Sigma_MX1_Measured_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]);
            vect_hist_Sigma_MX1_Measured_SigmaT.at(j).at(k) = hist_Sigma_MX1_Measured_SigmaT;
            TH1D* hist_Sigma_MX1_MET_Measured_SigmaT = new TH1D(("hist_Sigma_MX1_MET_Measured_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]);
            vect_hist_Sigma_MX1_MET_Measured_SigmaT.at(j).at(k) = hist_Sigma_MX1_MET_Measured_SigmaT;
            TH1D* hist_Sigma_MX1_Timing_Measured_SigmaT = new TH1D(("hist_Sigma_MX1_Timing_Measured_SigmaT"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaT[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]);
            vect_hist_Sigma_MX1_Timing_Measured_SigmaT.at(j).at(k) = hist_Sigma_MX1_Timing_Measured_SigmaT;
        }
        }if(MET_flag){
        TGraph* graph_Sigma_MX2_SigmaMET = new TGraph(NsigmaMET);
        graph_Sigma_MX2_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_MET_SigmaMET = new TGraph(NsigmaMET);
        graph_Sigma_MX2_MET_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_Timing_SigmaMET = new TGraph(NsigmaMET);
        graph_Sigma_MX2_Timing_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_SigmaMET_Measured = new TGraph(NsigmaMET);
        graph_Sigma_MX2_SigmaMET_Measured->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_MET_SigmaMET_Measured = new TGraph(NsigmaMET);
        graph_Sigma_MX2_MET_SigmaMET_Measured->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX2_Timing_SigmaMET_Measured = new TGraph(NsigmaMET);
        graph_Sigma_MX2_Timing_SigmaMET_Measured->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        
        TGraph* graph_Sigma_MX1_SigmaMET = new TGraph(NsigmaMET);
        graph_Sigma_MX1_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_MET_SigmaMET = new TGraph(NsigmaMET);
        graph_Sigma_MX1_MET_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_Timing_SigmaMET = new TGraph(NsigmaMET);
        graph_Sigma_MX1_Timing_SigmaMET->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_SigmaMET_Measured = new TGraph(NsigmaMET);
        graph_Sigma_MX1_SigmaMET_Measured->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_MET_SigmaMET_Measured = new TGraph(NsigmaMET);
        graph_Sigma_MX1_MET_SigmaMET_Measured->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        TGraph* graph_Sigma_MX1_Timing_SigmaMET_Measured = new TGraph(NsigmaMET);
        graph_Sigma_MX1_Timing_SigmaMET_Measured->SetName(("graph_sigmaMET"+std::to_string(j)).c_str());
        
        vect_graph_Sigma_MX2_SigmaMET.push_back(graph_Sigma_MX2_SigmaMET);
        vect_graph_Sigma_MX2_MET_SigmaMET.push_back(graph_Sigma_MX2_MET_SigmaMET);
        vect_graph_Sigma_MX2_Timing_SigmaMET.push_back(graph_Sigma_MX2_Timing_SigmaMET);
        vect_graph_Sigma_MX2_SigmaMET_Measured.push_back(graph_Sigma_MX2_SigmaMET_Measured);
        vect_graph_Sigma_MX2_MET_SigmaMET_Measured.push_back(graph_Sigma_MX2_MET_SigmaMET_Measured);
        vect_graph_Sigma_MX2_Timing_SigmaMET_Measured.push_back(graph_Sigma_MX2_Timing_SigmaMET_Measured);
        
        vect_graph_Sigma_MX1_SigmaMET.push_back(graph_Sigma_MX1_SigmaMET);
        vect_graph_Sigma_MX1_MET_SigmaMET.push_back(graph_Sigma_MX1_MET_SigmaMET);
        vect_graph_Sigma_MX1_Timing_SigmaMET.push_back(graph_Sigma_MX1_Timing_SigmaMET);
        vect_graph_Sigma_MX1_SigmaMET_Measured.push_back(graph_Sigma_MX1_SigmaMET_Measured);
        vect_graph_Sigma_MX1_MET_SigmaMET_Measured.push_back(graph_Sigma_MX1_MET_SigmaMET_Measured);
        vect_graph_Sigma_MX1_Timing_SigmaMET_Measured.push_back(graph_Sigma_MX1_Timing_SigmaMET_Measured);
        for(int k = 0; k < NsigmaMET; k++){
            TH1D* hist_Sigma_MX2_SigmaMET = new TH1D(("hist_Sigma_MX2_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2/ana_factor);
            vect_hist_Sigma_MX2_SigmaMET.at(j).at(k) = hist_Sigma_MX2_SigmaMET;
            TH1D* hist_Sigma_MX2_MET_SigmaMET = new TH1D(("hist_Sigma_MX2_MET_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2/ana_factor);
            vect_hist_Sigma_MX2_MET_SigmaMET.at(j).at(k) = hist_Sigma_MX2_MET_SigmaMET;
            TH1D* hist_Sigma_MX2_Timing_SigmaMET = new TH1D(("hist_Sigma_MX2_Timing_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2/ana_factor);
            vect_hist_Sigma_MX2_Timing_SigmaMET.at(j).at(k) = hist_Sigma_MX2_Timing_SigmaMET;
            TH1D* hist_Sigma_MX2_Measured_SigmaMET = new TH1D(("hist_Sigma_MX2_Measured_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_Measured_SigmaMET.at(j).at(k) = hist_Sigma_MX2_Measured_SigmaMET;
            TH1D* hist_Sigma_MX2_MET_Measured_SigmaMET = new TH1D(("hist_Sigma_MX2_MET_Measured_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_MET_Measured_SigmaMET.at(j).at(k) = hist_Sigma_MX2_MET_Measured_SigmaMET;
            TH1D* hist_Sigma_MX2_Timing_Measured_SigmaMET = new TH1D(("hist_Sigma_MX2_Timing_Measured_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX2,xmin_MX2,xmax_MX2);
            vect_hist_Sigma_MX2_Timing_Measured_SigmaMET.at(j).at(k) = hist_Sigma_MX2_Timing_Measured_SigmaMET;
            TH1D* hist_Sigma_MX1_SigmaMET = new TH1D(("hist_Sigma_MX1_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]/ana_factor);
            vect_hist_Sigma_MX1_SigmaMET.at(j).at(k) = hist_Sigma_MX1_SigmaMET;
            TH1D* hist_Sigma_MX1_MET_SigmaMET = new TH1D(("hist_Sigma_MX1_MET_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]/ana_factor);
            vect_hist_Sigma_MX1_MET_SigmaMET.at(j).at(k) = hist_Sigma_MX1_MET_SigmaMET;
            TH1D* hist_Sigma_MX1_Timing_SigmaMET = new TH1D(("hist_Sigma_MX1_Timing_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]/ana_factor);
            vect_hist_Sigma_MX1_Timing_SigmaMET.at(j).at(k) = hist_Sigma_MX1_Timing_SigmaMET;
            TH1D* hist_Sigma_MX1_Measured_SigmaMET = new TH1D(("hist_Sigma_MX1_Measured_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]);
            vect_hist_Sigma_MX1_Measured_SigmaMET.at(j).at(k) = hist_Sigma_MX1_Measured_SigmaMET;
            TH1D* hist_Sigma_MX1_MET_Measured_SigmaMET = new TH1D(("hist_Sigma_MX1_MET_Measured_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]);
            vect_hist_Sigma_MX1_MET_Measured_SigmaMET.at(j).at(k) = hist_Sigma_MX1_MET_Measured_SigmaMET;
            TH1D* hist_Sigma_MX1_Timing_Measured_SigmaMET = new TH1D(("hist_Sigma_MX1_Timing_Measured_SigmaMET"+std::to_string(int(ctau[j]))+std::to_string(int(sigmaMET[k]))).c_str(),"",bins_MX1,xmin_MX1,xmax_MX1[j]);
            vect_hist_Sigma_MX1_Timing_Measured_SigmaMET.at(j).at(k) = hist_Sigma_MX1_Timing_Measured_SigmaMET;
        }
        }
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
  
    // Declare observables for histogram booking
    HistPlot* histPlot = new HistPlot("Plots",
                                      std::string("#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}") +
                                      "#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}"+
                                      "Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");
    
    HistPlot* histPlot_Cos = new HistPlot("Plots_Cos",
                                            std::string("#tilde{#chi}_{2}^{ 0} #tilde{#chi}_{2}^{ 0}") +
                                            "#rightarrow Za(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}"+
                                            "Zb(#it{l}#it{l}) #tilde{#chi}_{1}^{ 0}");
    
    histPlot->SetRebin(1);
    histPlot_Cos->SetRebin(1);
    
    RFList<const HistPlotCategory> cat_list_ctau;
    RFList<const HistPlotCategory> cat_list_ctau_cos;
    char smassX2[200];
    string sctau = "c#tau = ";
    for(int m = 0; m < Nctau; m++){
        
        char snamectau[200], scatctau[50];
        sprintf(scatctau, "ctau_%d", m);
        sprintf(snamectau, "%.0f cm", ctau[m]);
        
        cat_list_ctau += histPlot->GetNewCategory(scatctau, sctau+std::string(snamectau));
        cat_list_ctau_cos += histPlot_Cos->GetNewCategory(scatctau, sctau+std::string(snamectau));
    }
    
    //setting up all the variables that could be potentially plotted
    const HistPlotVar& MIa2 = histPlot->GetNewVar("MIa2", "M(#tilde{#chi}_{1a}^{0})", 0., 600., "[GeV]");
    const HistPlotVar& MIb2 = histPlot->GetNewVar("MIb2", "M(#tilde{#chi}_{1b}^{0})", 0., 600., "[GeV]");
    const HistPlotVar& MXa2 = histPlot->GetNewVar("MXa2", "M(#tilde{#chi}_{2a}^{0})", 0., 800., "[GeV]");
    const HistPlotVar& MXa2_Cos = histPlot_Cos->GetNewVar("MXa2", "M(#tilde{#chi}_{2a}^{0})", 0., 800., "[GeV]");
    const HistPlotVar& MXb2 = histPlot->GetNewVar("MXb2", "M(#tilde{#chi}_{2b}^{0})", 0., 800., "[GeV]");
    const HistPlotVar& MIa = histPlot->GetNewVar("MIa", "M(#tilde{#chi}_{1}^{0})", 200., 500., "[GeV]");
    const HistPlotVar& MXa = histPlot->GetNewVar("MXa", "M(#tilde{#chi}_{2}^{0})", 0., 800., "[GeV]");
    const HistPlotVar& EZa = histPlot->GetNewVar("EZ", "E_{Z}^{#tilde{#chi}_{2}^{0}}", 50., 250., "[GeV]");
    const HistPlotVar& CosX2a_Plot = histPlot->GetNewVar("CosX2a", "Cos(#theta_{#tilde{#chi}_{2a}^{0}})", -1., 1., "");
    
    histPlot->AddPlot(MIa, cat_list_ctau);
    histPlot->AddPlot(MXa, cat_list_ctau);
    histPlot->AddPlot(MXa2, cat_list_ctau);
    histPlot->AddPlot(MIa2, cat_list_ctau);
    histPlot->AddPlot(MXa2, CosX2a_Plot, cat_list_ctau);
    histPlot->AddPlot(CosX2a_Plot, cat_list_ctau);
    histPlot->AddPlot(MXa2, MXb2, cat_list_ctau);
    histPlot->AddPlot(MIa2, MIb2, cat_list_ctau);
    
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
    PUPPI_Detector.Set_sigmaT(0.03/sqrt(2.)); //timing resolution
    PUPPI_Detector.Set_sigmaPV(20.0/10000.0); //Primary Vertex resolution
    PUPPI_Detector.Set_sigmaSV(65.0/10000.0); //secondary Vertex resolution
    Vertex PV(0.0,0.0,0.0,0.0); //nominal PV at 0
    Resolution test_Resolution(PUPPI_Detector); //setting up the Resolution "calculator"
    double LAB_Pt;
    double LAB_eta;
    double LAB_M;
    
    //relative uncertainty on the distance between PV and SV
    double sigmaDistance = sqrt((PUPPI_Detector.Get_sigmaPV()*PUPPI_Detector.Get_sigmaPV()+PUPPI_Detector.Get_sigmaSV()*PUPPI_Detector.Get_sigmaSV()));
    
    //for checking generator efficiency
    int gen_events = 0;
    int acp_events = 0;
    
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
        double MX2X2 = X2X2_Gen.GetMass();
        TLorentzVector Pa = X2a_Gen.GetFourVector();
        TLorentzVector Pb = X2b_Gen.GetFourVector();
        TLorentzVector sys = Pa+Pb;
        TLorentzVector Ia  = X1a_Gen.GetFourVector();
        TLorentzVector Ib  = X1b_Gen.GetFourVector();
        TLorentzVector I   = Ia+Ib;
        I.SetZ(0.0);
        decayangle = false;
        double MET_Mag_Resolution = PUPPI_Detector.Get_Sigma_Par(sys);
        double MET_Dir_Resolution = PUPPI_Detector.Get_Sigma_Perp(sys);
        
        //The smearing begins
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
        
        double ToFa = physics.Get_ToF(ctau[m], Pa);
        double ToFb = physics.Get_ToF(ctau[m], Pb);
        
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
        
        double Da = 30.*ToFa*vBetaaGen.Mag();
        double Db = 30.*ToFb*vBetabGen.Mag();
        
        if(fabs(Smeared_ToFa) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Smeared_ToFb) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Da) < displacement_cut*sigmaDistance || fabs(Db) < displacement_cut*sigmaDistance || Smeared_vBetaa.Mag() >= 1. || Smeared_vBetab.Mag() >= 1.) { igen--; continue;}
        
        //double Sigma_Beta_Mag = sqrt((1.0/(Smeared_ToFa*Smeared_ToFa))*(6.*sigmaDistance*sigmaDistance+2.*Smeared_vBetaa.Mag()*Smeared_vBetaa.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));
        //Factor of 6 since for the distance there is 6 dof: 1 for each spatial dimension and we have two vertices
        double Sigma_Beta_Mag = sqrt((1.0/(Smeared_ToFa*Smeared_ToFa))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetaa.Mag()*Smeared_vBetaa.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));
        
        TLorentzVector vZa = L1a_RECO + L2a_RECO;
        vZa.Boost(-Smeared_vBetaa);
        EZa = vZa.E();
        
        MXa = test_Resolution.Mass_Parent((L1a_RECO+L2a_RECO).Vect()+Ia_RECO,Smeared_vBetaa);
        MIa = test_Resolution.Mass_Invisible(Ia_RECO,L1a_RECO+L2a_RECO,Smeared_vBetaa);
        
        TLorentzVector Va = L1a_RECO+L2a_RECO;
        TLorentzVector Vb = L1b_RECO+L2b_RECO;
        
        MXa2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
        MXa2_Cos = MXa2;
        MXb2 = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetab,Smeared_vBetaa);
        
        double Sigma_Vis = test_Resolution.Visible_Resolution(L1a_RECOt.Vect(),L2a_RECOt.Vect(),L1b_RECOt.Vect(),L2b_RECOt.Vect(),L1a_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L1a_Gent),L2a_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L2a_Gent),L1b_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L1b_Gent),L2b_RECOt.Pt()*PUPPI_Detector.GetMuonResolution(L2b_Gent));
        
        double MXa2_Resolution = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);
        
        MIa2 = test_Resolution.Mass_Invisible2(MET_RECO_PUPPI,Va,Vb,Smeared_vBetaa,Smeared_vBetab);
        MIb2 = test_Resolution.Mass_Invisible2(MET_RECO_PUPPI,Vb,Va,Smeared_vBetab,Smeared_vBetaa);
        
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
        
        CosX2a_Plot = X2a_Reco.GetCosDecayAngle();
        double CosX2b = X2b_Reco.GetCosDecayAngle();
        if((fabs(CosX2a_Plot) < 0.8 || fabs(CosX2b) < 0.8))
        {
            decayangle = true;
        }
        
        if(timing_flag){
            PUPPI_Detector.Set_sigmaT(0.);
            Vertex Smeared_PV_No_Time = PUPPI_Detector.Smear_PV(PV);
            Vertex Smeared_SVa_No_Time = PUPPI_Detector.Smear_SV(SVa);
            Vertex Smeared_SVb_No_Time = PUPPI_Detector.Smear_SV(SVb);
            TVector3 Smeared_vBetaa_No_Time = PUPPI_Detector.Get_Beta(Smeared_PV_No_Time,Smeared_SVa_No_Time);
            TVector3 Smeared_vBetab_No_Time = PUPPI_Detector.Get_Beta(Smeared_PV_No_Time,Smeared_SVb_No_Time);
            
            double Da_No_Time = 30.*ToFa*Smeared_vBetaa_No_Time.Mag();
            double Db_No_Time = 30.*ToFb*Smeared_vBetab_No_Time.Mag();
            if(fabs(ToFa) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(ToFb) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Da_No_Time) < displacement_cut*sigmaDistance || fabs(Db_No_Time) < displacement_cut*sigmaDistance || Smeared_vBetaa_No_Time.Mag() >= 1. || Smeared_vBetab_No_Time.Mag() >= 1.)
            {
                igen--;
                continue;
            }
            
            double Sigma_Beta_Mag_No_Time = sqrt((1.0/(ToFa*ToFa))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetaa_No_Time.Mag()*Smeared_vBetaa_No_Time.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));
            
            double MXa2_Timing_Calc = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa_No_Time,Smeared_vBetab_No_Time);
            TLorentzVector vZa_Timing_Calc = L1a_RECO + L2a_RECO;
            vZa_Timing_Calc.Boost(-Smeared_vBetaa_No_Time);
            double EZa_Timing_Calc = vZa_Timing_Calc.E();
            double Mass_Vis = (L1a_RECO+L2a_RECO).M(); //Same for MET and Timing
            
            double MXa1_Calc[NsigmaT];
            double MXa1_MET_Calc[NsigmaT];
            double MXa1_Timing_Calc[NsigmaT];
            double MXa1_Res[NsigmaT];
            double MXa1_Res_MET[NsigmaT];
            double MXa1_Res_Timing[NsigmaT];
            
            for(int i = 0; i < NsigmaT; i++)
            {
                PUPPI_Detector.Set_sigmaT((sigmaT[i]/1000.0)/sqrt(2.));
                double Smeared_ToFa_Time = PUPPI_Detector.Smear_ToF(ToFa);
                double Smeared_ToFb_Time = PUPPI_Detector.Smear_ToF(ToFb);
                Vertex SVa_Time = physics.Get_SV(Smeared_ToFa_Time,Pa);
                Vertex SVb_Time = physics.Get_SV(Smeared_ToFb_Time,Pb);
                Vertex Smeared_PV_Time = PUPPI_Detector.Smear_PV(PV);
                Vertex Smeared_SVa_Time = PUPPI_Detector.Smear_SV(SVa);
                Vertex Smeared_SVb_Time = PUPPI_Detector.Smear_SV(SVb);
                TVector3 Smeared_vBetaa_Time = PUPPI_Detector.Get_Beta(Smeared_PV_Time,Smeared_SVa_Time);
                TVector3 Smeared_vBetab_Time = PUPPI_Detector.Get_Beta(Smeared_PV_Time,Smeared_SVb_Time);
                
                double Da_Time = 30.*Smeared_ToFa_Time*Smeared_vBetaa_Time.Mag();
                double Db_Time = 30.*Smeared_ToFb_Time*Smeared_vBetab_Time.Mag();
                if(fabs(Smeared_ToFa_Time) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Smeared_ToFb_Time) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Da_Time) < displacement_cut*sigmaDistance || fabs(Db_Time) < displacement_cut*sigmaDistance || Smeared_vBetaa_Time.Mag() >= 1. || Smeared_vBetab_Time.Mag() >= 1.)
                {
                    //i--;
                    continue;
                }
                
                double Sigma_Beta_Mag_Time = sqrt((1.0/(PUPPI_Detector.Smear_ToF(ToFa)*PUPPI_Detector.Smear_ToF(ToFa)))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetaa_Time.Mag()*Smeared_vBetaa_Time.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));
                
                //Begin Calculations:
                double MXa2_Calc = test_Resolution.Mass_Parents2(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa_Time,Smeared_vBetab_Time);
                double MXa2_MET_Calc = test_Resolution.Mass_Parents2(I_Vect,Va.Vect()+Vb.Vect(),Smeared_vBetaa_Time,Smeared_vBetab_Time);
                
                TLorentzVector vZa_Calc = L1a_RECO + L2a_RECO;
                vZa_Calc.Boost(-Smeared_vBetaa_Time);
                
                double EZa_Calc = vZa_Calc.E(); //Same for MET
                
                MXa1_Calc[i] = test_Resolution.Mass_Invisible2(MXa2_Calc, EZa_Calc, Mass_Vis);
                MXa1_MET_Calc[i] = test_Resolution.Mass_Invisible2(MXa2_MET_Calc, EZa_Calc, Mass_Vis);
                MXa1_Timing_Calc[i] = test_Resolution.Mass_Invisible2(MXa2_Timing_Calc, EZa_Timing_Calc, Mass_Vis);
                
                double MXa2_Res = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa_Time,Smeared_vBetab_Time,Sigma_Beta_Mag_Time,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);
                MXa1_Res[i] = test_Resolution.Mass_Invisible_Resolution2(MET_RECO_PUPPI,Va,Vb,Smeared_vBetaa_Time,Smeared_vBetab_Time,Sigma_Beta_Mag_Time,MET_Mag_Resolution,MET_Dir_Resolution);
                double MXa2_Res_MET = test_Resolution.Mass_Parents2_Resolution(I_Vect,Va.Vect()+Vb.Vect(),Smeared_vBetaa_Time,Smeared_vBetab_Time,Sigma_Beta_Mag_Time,0.,0.,Sigma_Vis);
                MXa1_Res_MET[i] = test_Resolution.Mass_Invisible_Resolution2(I_Vect,Va,Vb,Smeared_vBetaa_Time,Smeared_vBetab_Time,Sigma_Beta_Mag_Time,0.,0.);
                double MXa2_Res_Timing = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI,Va.Vect()+Vb.Vect(),Smeared_vBetaa_No_Time,Smeared_vBetab_No_Time,Sigma_Beta_Mag_No_Time,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);
                MXa1_Res_Timing[i] = test_Resolution.Mass_Invisible_Resolution2(MET_RECO_PUPPI,Va,Vb,Smeared_vBetaa_No_Time,Smeared_vBetab_No_Time,Sigma_Beta_Mag_No_Time,MET_Mag_Resolution,MET_Dir_Resolution);
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
                //if((fabs(CosX2a) > 0.8 || fabs(CosX2b) > 0.8))
                if((fabs(CosX2a) < 0.8 && fabs(CosX2b) < 0.8))// && vect_hist_Sigma_MX2_Measured_SigmaT_cos08.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX2_Measured_SigmaT_cos08.at(m).at(i)->Fill(MXa2_Calc);
                }
                if((fabs(CosX2a) < 0.9 && fabs(CosX2b) < 0.9))// && vect_hist_Sigma_MX2_Measured_SigmaT_cos09.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX2_Measured_SigmaT_cos09.at(m).at(i)->Fill(MXa2_Calc);
                }
                //Fill Vectors
                if(vect_hist_Sigma_MX2_SigmaT.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX2_SigmaT.at(m).at(i)->Fill(MXa2_Res/MXa2_Calc);
                    vect_hist_Sigma_MX2_MET_SigmaT.at(m).at(i)->Fill(MXa2_Res_MET/MXa2_MET_Calc);
                    vect_hist_Sigma_MX2_Timing_SigmaT.at(m).at(i)->Fill(MXa2_Res_Timing/MXa2_Timing_Calc);
                    vect_hist_Sigma_MX2_Measured_SigmaT.at(m).at(i)->Fill(MXa2_Calc);
                    vect_hist_Sigma_MX2_MET_Measured_SigmaT.at(m).at(i)->Fill(MXa2_MET_Calc);
                    vect_hist_Sigma_MX2_Timing_Measured_SigmaT.at(m).at(i)->Fill(MXa2_Timing_Calc);
                }
                if((MXa1_Calc[i] > 0.001) && vect_hist_Sigma_MX1_SigmaT.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX1_SigmaT.at(m).at(i)->Fill(MXa1_Res[i]/MXa1_Calc[i]);
                    vect_hist_Sigma_MX1_Measured_SigmaT.at(m).at(i)->Fill(MXa1_Calc[i]);
                }
                if((MXa1_MET_Calc[i] > 0.001) && vect_hist_Sigma_MX1_MET_SigmaT.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX1_MET_SigmaT.at(m).at(i)->Fill(MXa1_Res_MET[i]/MXa1_MET_Calc[i]);
                    vect_hist_Sigma_MX1_MET_Measured_SigmaT.at(m).at(i)->Fill(MXa1_MET_Calc[i]);
                }
                if((MXa1_Timing_Calc[i] > 0.001) && vect_hist_Sigma_MX1_Timing_SigmaT.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX1_Timing_SigmaT.at(m).at(i)->Fill(MXa1_Res_Timing[i]/MXa1_Timing_Calc[i]);
                    vect_hist_Sigma_MX1_Timing_Measured_SigmaT.at(m).at(i)->Fill(MXa1_Timing_Calc[i]);
                }
            }
            PUPPI_Detector.Set_sigmaT((sigmaT[0]/1000.)/sqrt(2.));
        }
        
        if(MET_flag){
            PUPPI_Detector.Set_sigmaT(0.);
            Vertex SVa_No_Time = physics.Get_SV(ToFa,Pa);
            Vertex SVb_No_Time = physics.Get_SV(ToFb,Pb);
            Vertex Smeared_PV_No_Time = PUPPI_Detector.Smear_PV(PV);
            Vertex Smeared_SVa_No_Time = PUPPI_Detector.Smear_SV(SVa);
            Vertex Smeared_SVb_No_Time = PUPPI_Detector.Smear_SV(SVb);
            TVector3 Smeared_vBetaa_No_Time = PUPPI_Detector.Get_Beta(Smeared_PV_No_Time,Smeared_SVa_No_Time);
            TVector3 Smeared_vBetab_No_Time = PUPPI_Detector.Get_Beta(Smeared_PV_No_Time,Smeared_SVb_No_Time);
            
            double Da_No_Time = 30.*ToFa*Smeared_vBetaa_No_Time.Mag();
            double Db_No_Time = 30.*ToFb*Smeared_vBetab_No_Time.Mag();
            if(fabs(ToFa) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(ToFb) < displacement_cut*PUPPI_Detector.Get_sigmaT() || fabs(Da_No_Time) < displacement_cut*sigmaDistance || fabs(Db_No_Time) < displacement_cut*sigmaDistance || Smeared_vBetaa_No_Time.Mag() >= 1. || Smeared_vBetab_No_Time.Mag() >= 1.)
            {
                igen--;
                continue;
            }
            
            double Sigma_Beta_Mag_No_Time = sqrt((1.0/(ToFa*ToFa))*(sigmaDistance*sigmaDistance+2.*Smeared_vBetaa_No_Time.Mag()*Smeared_vBetaa_No_Time.Mag()*PUPPI_Detector.Get_sigmaT()*PUPPI_Detector.Get_sigmaT()));
            PUPPI_Detector.Set_sigmaT((sigmaT[0]/1000.)/sqrt(2.));
            double MXa2_MET_Calc = test_Resolution.Mass_Parents2(I_Vect,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
            
            double MXa1_Calc[NsigmaMET];
            double MXa1_MET_Calc[NsigmaMET];
            double MXa1_Timing_Calc[NsigmaMET];
            double MXa1_Res[NsigmaMET];
            double MXa1_Res_MET[NsigmaMET];
            double MXa1_Res_Timing[NsigmaMET];
            PUPPI_Detector.Set_sigmaT((30./1000.)/sqrt(2.));
            
            for(int i = 0; i < NsigmaMET; i++)
            {
                PUPPI_Detector.Set_Sigma_Perp(sys,sigmaMET[i]/100.);
                PUPPI_Detector.Set_Sigma_Par(sys,sigmaMET[i]/100.);
                TVector3 MET_RECO_PUPPI_LOOP = PUPPI_Detector.Smear_MET(I_Vect);
                MET_RECO_PUPPI_LOOP.SetZ(0.0);
                MET_Mag_Resolution = PUPPI_Detector.Get_Sigma_Par();
                MET_Dir_Resolution = PUPPI_Detector.Get_Sigma_Perp();
                
                //Begin Calculations:
                double MXa2_Calc = test_Resolution.Mass_Parents2(MET_RECO_PUPPI_LOOP,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab);
                double MXa2_Timing_Calc = test_Resolution.Mass_Parents2(MET_RECO_PUPPI_LOOP,Va.Vect()+Vb.Vect(),Smeared_vBetaa_No_Time,Smeared_vBetab_No_Time);
                
                TLorentzVector vZa_Calc = L1a_RECO + L2a_RECO;
                vZa_Calc.Boost(-Smeared_vBetaa);
                TLorentzVector vZa_Timing_Calc = L1a_RECO + L2a_RECO;
                vZa_Timing_Calc.Boost(-Smeared_vBetaa_No_Time);
                
                double EZa_Calc = vZa_Calc.E(); //Same for MET
                double EZa_Timing_Calc = vZa_Timing_Calc.E();
                double Mass_Vis = (L1a_RECO+L2a_RECO).M(); //Same for MET and Timing
                
                MXa1_Calc[i] = test_Resolution.Mass_Invisible2(MXa2_Calc, EZa_Calc, Mass_Vis);
                MXa1_MET_Calc[i] = test_Resolution.Mass_Invisible2(MXa2_MET_Calc, EZa_Calc, Mass_Vis);
                MXa1_Timing_Calc[i] = test_Resolution.Mass_Invisible2(MXa2_Timing_Calc, EZa_Timing_Calc, Mass_Vis);
                
                double MXa2_Res = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI_LOOP,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);
                MXa1_Res[i] = test_Resolution.Mass_Invisible_Resolution2(MET_RECO_PUPPI_LOOP,Va,Vb,Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,MET_Mag_Resolution,MET_Dir_Resolution);
                double MXa2_Res_MET = test_Resolution.Mass_Parents2_Resolution(I_Vect,Va.Vect()+Vb.Vect(),Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,0.,0.,Sigma_Vis);
                MXa1_Res_MET[i] = test_Resolution.Mass_Invisible_Resolution2(I_Vect,Va,Vb,Smeared_vBetaa,Smeared_vBetab,Sigma_Beta_Mag,0.,0.);
                double MXa2_Res_Timing = test_Resolution.Mass_Parents2_Resolution(MET_RECO_PUPPI_LOOP,Va.Vect()+Vb.Vect(),Smeared_vBetaa_No_Time,Smeared_vBetab_No_Time,Sigma_Beta_Mag_No_Time,MET_Mag_Resolution,MET_Dir_Resolution,Sigma_Vis);
                MXa1_Res_Timing[i] = test_Resolution.Mass_Invisible_Resolution2(MET_RECO_PUPPI_LOOP,Va,Vb,Smeared_vBetaa_No_Time,Smeared_vBetab_No_Time,Sigma_Beta_Mag_No_Time,MET_Mag_Resolution,MET_Dir_Resolution);
                /*
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
                if((fabs(CosX2a) > 0.8 || fabs(CosX2b) > 0.8))
                {
                    //i--;
                    continue;
                }
                }*/
                //Fill Vectors
                if(vect_hist_Sigma_MX2_SigmaMET.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX2_SigmaMET.at(m).at(i)->Fill(MXa2_Res/MXa2_Calc);
                    vect_hist_Sigma_MX2_MET_SigmaMET.at(m).at(i)->Fill(MXa2_Res_MET/MXa2_MET_Calc);
                    vect_hist_Sigma_MX2_Timing_SigmaMET.at(m).at(i)->Fill(MXa2_Res_Timing/MXa2_Timing_Calc);
                    vect_hist_Sigma_MX2_Measured_SigmaMET.at(m).at(i)->Fill(MXa2_Calc);
                    vect_hist_Sigma_MX2_MET_Measured_SigmaMET.at(m).at(i)->Fill(MXa2_MET_Calc);
                    vect_hist_Sigma_MX2_Timing_Measured_SigmaMET.at(m).at(i)->Fill(MXa2_Timing_Calc);
                }
                if((MXa1_Calc[i] > 0.001) && vect_hist_Sigma_MX1_SigmaMET.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX1_SigmaMET.at(m).at(i)->Fill(MXa1_Res[i]/MXa1_Calc[i]);
                    vect_hist_Sigma_MX1_Measured_SigmaMET.at(m).at(i)->Fill(MXa1_Calc[i]);
                }
                if((MXa1_MET_Calc[i] > 0.001) && vect_hist_Sigma_MX1_MET_SigmaMET.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX1_MET_SigmaMET.at(m).at(i)->Fill(MXa1_Res_MET[i]/MXa1_MET_Calc[i]);
                    vect_hist_Sigma_MX1_MET_Measured_SigmaMET.at(m).at(i)->Fill(MXa1_MET_Calc[i]);
                }
                if((MXa1_Timing_Calc[i] > 0.001) && vect_hist_Sigma_MX1_Timing_SigmaMET.at(m).at(i)->GetEntries() < Entries)
                {
                    vect_hist_Sigma_MX1_Timing_SigmaMET.at(m).at(i)->Fill(MXa1_Res_Timing[i]/MXa1_Timing_Calc[i]);
                    vect_hist_Sigma_MX1_Timing_Measured_SigmaMET.at(m).at(i)->Fill(MXa1_Timing_Calc[i]);
                }
            }
            MET_Mag_Resolution = PUPPI_Detector.Get_Sigma_Par(sys);
            MET_Dir_Resolution = PUPPI_Detector.Get_Sigma_Perp(sys);
        }
        //if((MIa2 > 0. && MIb2 > 0.) && MIa > 0.){
        if(decayangle){
        histPlot_Cos->Fill(cat_list_ctau_cos[m]);
        histPlot->Fill(cat_list_ctau[m]);
        acp_events++;
        }
        else{
        histPlot->Fill(cat_list_ctau[m]);
        acp_events++;
        }
        //}
    }
    LAB_Gen.PrintGeneratorEfficiency();
  }
    g_Log << LogInfo << "Generated a Total of " << gen_events << " Events " << LogEnd;
    g_Log << LogInfo << acp_events << " passed selection requirements " << LogEnd;
    g_Log << LogInfo << "Efficiency: " << 100.0*acp_events/gen_events << "%" << LogEnd;
    end = gSystem->Now();
    g_Log << LogInfo << "Time to Generate " << Ngen*Nctau << " Events: " << (end-start)/1000.0 << " seconds" << LogEnd;
    g_Log << LogInfo << "Processing " << Ngen*Nctau << " Events" << LogEnd;
    histPlot->Draw(true);
    //histPlot_Cos->Draw(true);
    TFile fout(output_name.c_str(),"RECREATE");
    /*if(MET_flag || timing_flag){
    TCanvas* canv = new TCanvas("canv","",750,500);
    Draw_Hists(vect_hist_Sigma_MX1_Measured_SigmaT.at(2),canv);
    TCanvas* canv2 = new TCanvas("canv2","",750,500);
    Draw_Hists(vect_hist_Sigma_MX1_Measured_SigmaT.at(1),canv2);
    TCanvas* canv1 = new TCanvas("canv1","",750,500);
    Draw_Hists(vect_hist_Sigma_MX1_Measured_SigmaT.at(0),canv1);
    }*/
    if(timing_flag){
    for(int i = 0; i<Nctau; i++)
    {
        g_Log << LogInfo << "Computing Points For timing and ctau = " << ctau[i] << LogEnd;
    for(int j = 0; j<NsigmaT; j++)
    {
        vect_graph_Sigma_MX2_SigmaT.at(i)->SetPoint(j,sigmaT[j],Hist_Mode(*vect_hist_Sigma_MX2_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX2_MET_SigmaT.at(i)->SetPoint(j,sigmaT[j],Hist_Mode(*vect_hist_Sigma_MX2_MET_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX2_Timing_SigmaT.at(i)->SetPoint(j,sigmaT[j],Hist_Mode(*vect_hist_Sigma_MX2_Timing_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX2_SigmaT_Measured.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX2_Measured_SigmaT.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_Measured_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX2_SigmaT_Measured_cos08.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX2_Measured_SigmaT_cos08.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_Measured_SigmaT_cos08.at(i).at(j)));
        vect_graph_Sigma_MX2_SigmaT_Measured_cos09.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX2_Measured_SigmaT_cos09.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_Measured_SigmaT_cos09.at(i).at(j)));
        vect_graph_Sigma_MX2_MET_SigmaT_Measured.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX2_MET_Measured_SigmaT.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_MET_Measured_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX2_Timing_SigmaT_Measured.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX2_Timing_Measured_SigmaT.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_Timing_Measured_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX1_SigmaT.at(i)->SetPoint(j,sigmaT[j],Hist_Mode(*vect_hist_Sigma_MX1_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX1_MET_SigmaT.at(i)->SetPoint(j,sigmaT[j],Hist_Mode(*vect_hist_Sigma_MX1_MET_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX1_Timing_SigmaT.at(i)->SetPoint(j,sigmaT[j],Hist_Mode(*vect_hist_Sigma_MX1_Timing_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX1_SigmaT_Measured.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX1_Measured_SigmaT.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX1_Measured_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX1_MET_SigmaT_Measured.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX1_MET_Measured_SigmaT.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX1_MET_Measured_SigmaT.at(i).at(j)));
        vect_graph_Sigma_MX1_Timing_SigmaT_Measured.at(i)->SetPoint(j,sigmaT[j],Hist_FWHM(*vect_hist_Sigma_MX1_Timing_Measured_SigmaT.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX1_Timing_Measured_SigmaT.at(i).at(j)));
        
        delete vect_hist_Sigma_MX2_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX2_MET_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX2_Timing_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX2_Measured_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX2_Measured_SigmaT_cos08.at(i).at(j);
        delete vect_hist_Sigma_MX2_Measured_SigmaT_cos09.at(i).at(j);
        delete vect_hist_Sigma_MX2_MET_Measured_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX2_Timing_Measured_SigmaT.at(i).at(j);
        
        delete vect_hist_Sigma_MX1_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX1_MET_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX1_Timing_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX1_Measured_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX1_MET_Measured_SigmaT.at(i).at(j);
        delete vect_hist_Sigma_MX1_Timing_Measured_SigmaT.at(i).at(j);
    }
    }
    vector<string> leg_text_Sigma_MX2_SigmaT;
    for(int j = 0; j < Nctau; j++){leg_text_Sigma_MX2_SigmaT.push_back("c#tau "+std::to_string(int(ctau.at(j))));}
    gStyle->SetOptTitle(1);
    vect_graph_Sigma_MX2_SigmaT.at(0)->SetTitle("Analytical: Everything On");
    vect_graph_Sigma_MX2_MET_SigmaT.at(0)->SetTitle("Analytical: MET Off");
    vect_graph_Sigma_MX2_Timing_SigmaT.at(0)->SetTitle("Analytical: Timing Off");
    vect_graph_Sigma_MX2_SigmaT_Measured.at(0)->SetTitle("Measured: Everything On");
    vect_graph_Sigma_MX2_MET_SigmaT_Measured.at(0)->SetTitle("Measured: MET Off");
    vect_graph_Sigma_MX2_Timing_SigmaT_Measured.at(0)->SetTitle("Measured: Timing Off");
    vect_graph_Sigma_MX1_SigmaT.at(0)->SetTitle("Analytical: Everything On");
    vect_graph_Sigma_MX1_MET_SigmaT.at(0)->SetTitle("Analytical: MET Off");
    vect_graph_Sigma_MX1_Timing_SigmaT.at(0)->SetTitle("Analytical: Timing Off");
    vect_graph_Sigma_MX1_SigmaT_Measured.at(0)->SetTitle("Measured: Everything On");
    vect_graph_Sigma_MX1_MET_SigmaT_Measured.at(0)->SetTitle("Measured: MET Off");
    vect_graph_Sigma_MX1_Timing_SigmaT_Measured.at(0)->SetTitle("Measured: Timing Off");
    
    //Draw_Graphs(fout, vect_graph_Sigma_MX2_SigmaT, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LLP_Analytical_ctau_timing_Both", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX2_MET_SigmaT, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LLP_Analytical_ctau_timing_MET_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX2_Timing_SigmaT, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LLP_Analytical_ctau_timing_Timing_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX2_SigmaT_Measured, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LLP_Measured_ctau_timing_Both", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX2_MET_SigmaT_Measured, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LLP_Measured_ctau_timing_MET_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX2_Timing_SigmaT_Measured, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LLP_Measured_ctau_timing_Timing_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX1_SigmaT, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LSP_Analytical_ctau_timing_Both", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX1_MET_SigmaT, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LSP_Analytical_ctau_timing_MET_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX1_Timing_SigmaT, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LSP_Analytical_ctau_timing_Timing_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX1_SigmaT_Measured, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LSP_Measured_ctau_timing_Both", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX1_MET_SigmaT_Measured, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LSP_Measured_ctau_timing_MET_Off", points);
    //Draw_Graphs(fout, vect_graph_Sigma_MX1_Timing_SigmaT_Measured, leg_text_Sigma_MX2_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "Res_LSP_Measured_ctau_timing_Timing_Off", points);
       /*
        vector<TGraph*> vect_graph_ctau_SigmaT_SigmaMX2;
        //vect_graph_ctau_SigmaT_SigmaMX2.push_back(vect_graph_Sigma_MX2_SigmaT[0]);
        //vect_graph_ctau_SigmaT_SigmaMX2.push_back(vect_graph_Sigma_MX2_MET_SigmaT[0]);
        //vect_graph_ctau_SigmaT_SigmaMX2.push_back(vect_graph_Sigma_MX2_Timing_SigmaT[0]);
        vect_graph_ctau_SigmaT_SigmaMX2.push_back(vect_graph_Sigma_MX2_SigmaT_Measured[2]);
        vect_graph_ctau_SigmaT_SigmaMX2.push_back(vect_graph_Sigma_MX2_MET_SigmaT_Measured[2]);
        vect_graph_ctau_SigmaT_SigmaMX2.push_back(vect_graph_Sigma_MX2_Timing_SigmaT_Measured[2]);
        vector<string> leg_text_ctau_SigmaT;
        //leg_text_ctau_SigmaT.push_back("Analytical");
        leg_text_ctau_SigmaT.push_back("Measured");
        leg_text_ctau_SigmaT.push_back("#sigma_{MET} Off");
        leg_text_ctau_SigmaT.push_back("#sigma_{t} Off");
        //leg_text_ctau_SigmaT.push_back("Measured");
        Draw_Graphs(fout, vect_graph_ctau_SigmaT_SigmaMX2, leg_text_ctau_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "MLLP_Timing_ctau", points);
        */
        vector<TGraph*> vect_graph_ctau_SigmaT_SigmaMX1;
        vect_graph_ctau_SigmaT_SigmaMX1.push_back(vect_graph_Sigma_MX1_SigmaT[0]);
        vect_graph_ctau_SigmaT_SigmaMX1.push_back(vect_graph_Sigma_MX1_MET_SigmaT[0]);
        vect_graph_ctau_SigmaT_SigmaMX1.push_back(vect_graph_Sigma_MX1_Timing_SigmaT[0]);
        vect_graph_ctau_SigmaT_SigmaMX1.push_back(vect_graph_Sigma_MX1_SigmaT_Measured[0]);
        vect_graph_ctau_SigmaT_SigmaMX1.push_back(vect_graph_Sigma_MX1_MET_SigmaT_Measured[0]);
        vect_graph_ctau_SigmaT_SigmaMX1.push_back(vect_graph_Sigma_MX1_Timing_SigmaT_Measured[0]);
        //Draw_Graphs(fout, vect_graph_ctau_SigmaT_SigmaMX1, leg_text_ctau_SigmaT, "FWHM M [%]", "#sigma_{t} [ps]", "MLSP_Timing_ctau", points);
        //For cos decay angle
        vector<TGraph*> vect_graph_ctau_cos0;
        vect_graph_ctau_cos0.push_back(vect_graph_Sigma_MX2_SigmaT_Measured[0]);
        //vect_graph_ctau_cos0.push_back(vect_graph_Sigma_MX2_SigmaT_Measured_cos09[0]);
        vect_graph_ctau_cos0.push_back(vect_graph_Sigma_MX2_SigmaT_Measured_cos08[0]);
        vect_graph_ctau_cos0.push_back(vect_graph_Sigma_MX2_MET_SigmaT_Measured[0]);
        vect_graph_ctau_cos0.push_back(vect_graph_Sigma_MX2_Timing_SigmaT_Measured[0]);
        vector<string> leg_text_ctau_cos0;
        //leg_text_ctau_cos0.push_back("No |Cos(#theta_{#tilde{#chi}_{2a}^{0}})| Cut");
        leg_text_ctau_cos0.push_back("Expected");
        //leg_text_ctau_cos0.push_back("|Cos(#theta_{#tilde{#chi}_{2a}^{0}})| < 0.9");
        leg_text_ctau_cos0.push_back("|Cos(#theta_{#tilde{#chi}_{2a}^{0}})| < 0.8");
        leg_text_ctau_cos0.push_back("#sigma_{MET} Off");
        leg_text_ctau_cos0.push_back("#sigma_{t} Off");
        //Draw_Graphs(fout, vect_graph_ctau_cos0, leg_text_ctau_cos0, "FWHM M [%]", "#sigma_{t} [ps]", "MLLP_Timing_Cos0", points);
        vector<TGraph*> vect_graph_ctau_cos2;
        vect_graph_ctau_cos2.push_back(vect_graph_Sigma_MX2_SigmaT_Measured[2]);
        //vect_graph_ctau_cos2.push_back(vect_graph_Sigma_MX2_SigmaT_Measured_cos09[2]);
        vect_graph_ctau_cos2.push_back(vect_graph_Sigma_MX2_SigmaT_Measured_cos08[2]);
        vect_graph_ctau_cos2.push_back(vect_graph_Sigma_MX2_MET_SigmaT_Measured[2]);
        vect_graph_ctau_cos2.push_back(vect_graph_Sigma_MX2_Timing_SigmaT_Measured[2]);
        vector<string> leg_text_ctau_cos2;
        //leg_text_ctau_cos2.push_back("No |Cos(#theta_{#tilde{#chi}_{2a}^{0}})| Cut");
        leg_text_ctau_cos2.push_back("Expected");
        //leg_text_ctau_cos2.push_back("|Cos(#theta_{#tilde{#chi}_{2a}^{0}})| < 0.9");
        leg_text_ctau_cos2.push_back("|Cos(#theta_{#tilde{#chi}_{2a}^{0}})| < 0.8");
        leg_text_ctau_cos2.push_back("#sigma_{MET} Off");
        leg_text_ctau_cos2.push_back("#sigma_{t} Off");
        Draw_Graphs(fout, vect_graph_ctau_cos2, leg_text_ctau_cos2, "FWHM M [%]", "#sigma_{t} [ps]", "MLLP_Timing_Cos2", points);
    }
    if(MET_flag){
    for(int i = 0; i<Nctau; i++)
    {
        g_Log << LogInfo << "Computing Points For MET and ctau = " << ctau[i] << LogEnd;
        for(int j = 0; j<NsigmaMET; j++)
        {
            vect_graph_Sigma_MX2_SigmaMET.at(i)->SetPoint(j,sigmaMET[j],Hist_Mode(*vect_hist_Sigma_MX2_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX2_MET_SigmaMET.at(i)->SetPoint(j,sigmaMET[j],Hist_Mode(*vect_hist_Sigma_MX2_MET_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX2_Timing_SigmaMET.at(i)->SetPoint(j,sigmaMET[j],Hist_Mode(*vect_hist_Sigma_MX2_Timing_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX2_SigmaMET_Measured.at(i)->SetPoint(j,sigmaMET[j],Hist_FWHM(*vect_hist_Sigma_MX2_Measured_SigmaMET.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_Measured_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX2_MET_SigmaMET_Measured.at(i)->SetPoint(j,sigmaMET[j],Hist_FWHM(*vect_hist_Sigma_MX2_MET_Measured_SigmaMET.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_MET_Measured_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX2_Timing_SigmaMET_Measured.at(i)->SetPoint(j,sigmaMET[j],Hist_FWHM(*vect_hist_Sigma_MX2_Timing_Measured_SigmaMET.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX2_Timing_Measured_SigmaMET.at(i).at(j)));
            
            vect_graph_Sigma_MX1_SigmaMET.at(i)->SetPoint(j,sigmaMET[j],Hist_Mode(*vect_hist_Sigma_MX1_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX1_MET_SigmaMET.at(i)->SetPoint(j,sigmaMET[j],Hist_Mode(*vect_hist_Sigma_MX1_MET_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX1_Timing_SigmaMET.at(i)->SetPoint(j,sigmaMET[j],Hist_Mode(*vect_hist_Sigma_MX1_Timing_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX1_SigmaMET_Measured.at(i)->SetPoint(j,sigmaMET[j],Hist_FWHM(*vect_hist_Sigma_MX1_Measured_SigmaMET.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX1_Measured_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX1_MET_SigmaMET_Measured.at(i)->SetPoint(j,sigmaMET[j],Hist_FWHM(*vect_hist_Sigma_MX1_MET_Measured_SigmaMET.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX1_MET_Measured_SigmaMET.at(i).at(j)));
            vect_graph_Sigma_MX1_Timing_SigmaMET_Measured.at(i)->SetPoint(j,sigmaMET[j],Hist_FWHM(*vect_hist_Sigma_MX1_Timing_Measured_SigmaMET.at(i).at(j))/Hist_Mode(*vect_hist_Sigma_MX1_Timing_Measured_SigmaMET.at(i).at(j)));
            
            delete vect_hist_Sigma_MX2_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX2_MET_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX2_Timing_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX2_Measured_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX2_MET_Measured_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX2_Timing_Measured_SigmaMET.at(i).at(j);
                
            delete vect_hist_Sigma_MX1_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX1_MET_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX1_Timing_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX1_Measured_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX1_MET_Measured_SigmaMET.at(i).at(j);
            delete vect_hist_Sigma_MX1_Timing_Measured_SigmaMET.at(i).at(j);
        }
        }
        vector<string> leg_text_Sigma_MX2_SigmaMET;
        for(int j = 0; j < Nctau; j++){leg_text_Sigma_MX2_SigmaMET.push_back("c#tau "+std::to_string(int(ctau.at(j))));}
        gStyle->SetOptTitle(1);
        vect_graph_Sigma_MX2_SigmaMET.at(0)->SetTitle("Analytical: Everything On");
        vect_graph_Sigma_MX2_MET_SigmaMET.at(0)->SetTitle("Analytical: MET Off");
        vect_graph_Sigma_MX2_Timing_SigmaMET.at(0)->SetTitle("Analytical: Timing Off");
        vect_graph_Sigma_MX2_SigmaMET_Measured.at(0)->SetTitle("Measured: Everything On");
        vect_graph_Sigma_MX2_MET_SigmaMET_Measured.at(0)->SetTitle("Measured: MET Off");
        vect_graph_Sigma_MX2_Timing_SigmaMET_Measured.at(0)->SetTitle("Measured: Timing Off");
        vect_graph_Sigma_MX1_SigmaMET.at(0)->SetTitle("Analytical: Everything On");
        vect_graph_Sigma_MX1_MET_SigmaMET.at(0)->SetTitle("Analytical: MET Off");
        vect_graph_Sigma_MX1_Timing_SigmaMET.at(0)->SetTitle("Analytical: Timing Off");
        vect_graph_Sigma_MX1_SigmaMET_Measured.at(0)->SetTitle("Measured: Everything On");
        vect_graph_Sigma_MX1_MET_SigmaMET_Measured.at(0)->SetTitle("Measured: MET Off");
        vect_graph_Sigma_MX1_Timing_SigmaMET_Measured.at(0)->SetTitle("Measured: Timing Off");
        
        //Draw_Graphs(fout, vect_graph_Sigma_MX2_SigmaMET, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LLP_Analytical_ctau_met_Both", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX2_MET_SigmaMET, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LLP_Analytical_ctau_met_MET_Off", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX2_Timing_SigmaMET, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LLP_Analytical_ctau_met_Timing_Off", points);
        Draw_Graphs(fout, vect_graph_Sigma_MX2_SigmaMET_Measured, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LLP_Measured_ctau_met_Both", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX2_MET_SigmaMET_Measured, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LLP_Measured_ctau_met_MET_Off", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX2_Timing_SigmaMET_Measured, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LLP_Measured_ctau_met_Timing_Off", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX1_SigmaMET, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LSP_Analytical_ctau_met_Both", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX1_MET_SigmaMET, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LSP_Analytical_ctau_met_MET_Off", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX1_Timing_SigmaMET, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LSP_Analytical_ctau_met_Timing_Off", points);
        Draw_Graphs(fout, vect_graph_Sigma_MX1_SigmaMET_Measured, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LSP_Measured_ctau_met_Both", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX1_MET_SigmaMET_Measured, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LSP_Measured_ctau_met_MET_Off", points);
        //Draw_Graphs(fout, vect_graph_Sigma_MX1_Timing_SigmaMET_Measured, leg_text_Sigma_MX2_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "Res_LSP_Measured_ctau_met_Timing_Off", points);
        vector<TGraph*> vect_graph_ctau_SigmaMET_SigmaMX2;
        vect_graph_ctau_SigmaMET_SigmaMX2.push_back(vect_graph_Sigma_MX2_SigmaMET[0]);
        vect_graph_ctau_SigmaMET_SigmaMX2.push_back(vect_graph_Sigma_MX2_MET_SigmaMET[0]);
        vect_graph_ctau_SigmaMET_SigmaMX2.push_back(vect_graph_Sigma_MX2_Timing_SigmaMET[0]);
        vect_graph_ctau_SigmaMET_SigmaMX2.push_back(vect_graph_Sigma_MX2_SigmaMET_Measured[0]);
        vect_graph_ctau_SigmaMET_SigmaMX2.push_back(vect_graph_Sigma_MX2_MET_SigmaMET_Measured[0]);
        vect_graph_ctau_SigmaMET_SigmaMX2.push_back(vect_graph_Sigma_MX2_Timing_SigmaMET_Measured[0]);
        vector<string> leg_text_ctau_SigmaMET;
        leg_text_ctau_SigmaMET.push_back("Analytical");
        leg_text_ctau_SigmaMET.push_back("#sigma_{MET} Off");
        leg_text_ctau_SigmaMET.push_back("#sigma_{t} Off");
        leg_text_ctau_SigmaMET.push_back("Measured");
        //Draw_Graphs(fout, vect_graph_ctau_SigmaMET_SigmaMX2, leg_text_ctau_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "MLLP_MET_ctau", points);
        vector<TGraph*> vect_graph_ctau_SigmaMET_SigmaMX1;
        vect_graph_ctau_SigmaMET_SigmaMX1.push_back(vect_graph_Sigma_MX1_SigmaMET[0]);
        vect_graph_ctau_SigmaMET_SigmaMX1.push_back(vect_graph_Sigma_MX1_MET_SigmaMET[0]);
        vect_graph_ctau_SigmaMET_SigmaMX1.push_back(vect_graph_Sigma_MX1_Timing_SigmaMET[0]);
        vect_graph_ctau_SigmaMET_SigmaMX1.push_back(vect_graph_Sigma_MX1_SigmaMET_Measured[0]);
        vect_graph_ctau_SigmaMET_SigmaMX1.push_back(vect_graph_Sigma_MX1_MET_SigmaMET_Measured[0]);
        vect_graph_ctau_SigmaMET_SigmaMX1.push_back(vect_graph_Sigma_MX1_Timing_SigmaMET_Measured[0]);
        //Draw_Graphs(fout, vect_graph_ctau_SigmaMET_SigmaMX1, leg_text_ctau_SigmaMET, "FWHM M [%]", "#sigma_{MET} [%]", "MLSP_MET_ctau", points);
    }
  fout.Close();
  histPlot->WriteOutput(output_name);
  histPlot->WriteHist(output_name);
  //histPlot_Cos->WriteOutput(output_name);
  //histPlot_Cos->WriteHist(output_name);
  treePlot->WriteOutput(output_name);
  g_Log << LogInfo << "Finished" << LogEnd;
  g_Log << LogInfo << "Time to Process " << Ngen*Nctau << " Events: " << (Long64_t(gSystem->Now())-end)/1000.0 << " seconds" << LogEnd;
}
