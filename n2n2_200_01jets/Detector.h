#ifndef DETECTOR_H
#define DETECTOR_H
#include "Analyze_n2n2j.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <string>
#include "Bonus.h"

class Detector
{
private:
    TLorentzVector PV_MC;
    TLorentzVector PV_RECO;
    TLorentzVector SV_MC_A;
    TLorentzVector SV_MC_B;
    TLorentzVector SV_RECO_A;
    TLorentzVector SV_RECO_B;
    TVector3 PT_true;
    TVector3 MET;
    //TLorentzVector BeamSpot; // ???
    double Proper_Decay_Time;
    double Time_Resolution;

public:
    
    Detector();
    
};
#endif
