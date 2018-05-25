#include "Detector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <string>

double Decay_Time=1.0; //need to decide units
double Time_Resolution = 30.0;

TVector3 Detector::PT_true_TV3(TLorentzVector* n1A_MC, TLorentzVector* n1B_MC)
{
    TVector3 PT_true;
    TVector3 n1A_vect=n1A_MC->Vect();
    TVector3 n1B_vect=n1B_MC->Vect();
    PT_true=n1A_vect+n1B_vect;
    PT_true.SetZ(0.0);
    return PT_true;
}

// smear MET given true MET and sumPT
TVector3 Detector::Smear_MET(const TVector3& PT_true, double sum_PT)
{
    TVector3 zhat;
    zhat.SetXYZ(0.0,0.0,1.0);
    
    double sigma_perp = sum_PT*sqrt((4.20687*4.20687)/(sum_PT*sum_PT)+(0.548129*0.548129)/sum_PT+0.0116634*0.0116634);
    double sigma_par  = sum_PT*sqrt((4.02779*4.02779)/(sum_PT*sum_PT)+(0.512707*0.512707)/sum_PT+0.0431817*0.0431817);
    
    TVector3 parhat  = PT_true.Unit();
    TVector3 perphat = PT_true.Cross(zhat).Unit();
    
    TVector3 MET = PT_true + gRandom->Gaus(0.,sigma_par)*parhat + gRandom->Gaus(0.,sigma_perp)*perphat;
    
    return MET;
}

// smear Electron Given TLV of Electron
TLorentzVector Detector::Smear_Electron(const TLorentzVector& Electron_MC_TLV)
{
    TLorentzVector Electron_Smear;
    return Electron_Smear;
}

// smear Muon Given TLV of Muon
TLorentzVector Detector::Smear_Muon(const TLorentzVector& Muon_MC_TLV)
{
    TLorentzVector Muon_Smear;
    return Muon_Smear;
}

TLorentzVector Detector::Create_PV(const TLorentzVector& n2A_MC, const TLorentzVector& n2B_MC) //create a Primary Vertex
{
    TLorentzVector PV_MC;
    return PV_MC;
}

//smear Primary Vertex
TLorentzVector Detector::Smear_PV(const TLorentzVector& PV_MC)//, double Time_Resolution)
{
    //Time_Resolution = 100.0;
    TLorentzVector PV_Smear;
    return PV_Smear;
}

//create Secondary Vertex
TLorentzVector Detector::Create_SV(const TLorentzVector& n2_MC)//, double Decay_Time)
{
    //Decay_Time=1202.012;
    TLorentzVector SV_MC;
    return SV_MC;
}

//smear Secondary Vertex
TLorentzVector Detector::Smear_SV(const TLorentzVector& SV_MC)//, double Time_Resolution)
{
    TLorentzVector SV_Smear;
    return SV_Smear;
}
