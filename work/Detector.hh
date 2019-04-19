#ifndef DETECTOR_H
#define DETECTOR_H
#include <TLorentzVector.h>
#include <TRandom.h>
#include "Vertex.hh"

class Detector
{
    
private:

    double con0_par;
    double con1_par;
    double con2_par;
    double con0_perp;
    double con1_perp;
    double con2_perp;
    
    double sigmaT;
    double sigmaPV;
    double sigmaSV;
    double sigma_par;
    double sigma_perp;

public:
    
    Detector();
    
    void Set_con0_par(double User_con0_par);
    double Get_con0_par();
    
    void Set_con1_par(double User_con1_par);
    double Get_con1_par();
    
    void Set_con2_par(double User_con2_par);
    double Get_con2_par();
    
    void Set_con0_perp(double User_con0_perp);
    double Get_con0_perp();
    
    void Set_con1_perp(double User_con1_perp);
    double Get_con1_perp();
    
    void Set_con2_perp(double User_con2_perp);
    double Get_con2_perp();
    
    void Set_sigmaT(double User_sigmaT);
    double Get_sigmaT();
    
    void Set_sigmaPV(double User_sigmaPV);
    double Get_sigmaPV();
    
    void Set_sigmaSV(double User_sigmaSV);
    double Get_sigmaSV();
    
    double Get_Sigma_Perp(const TLorentzVector& sys);
    double Get_Sigma_Par(const TLorentzVector& sys);
    
    TVector3 Smear_MET(const TVector3& inv);
    
    TVector3 Smear_MET_Mag(TVector3& inv);
    
    TVector3 Smear_MET_Dir(TVector3& inv);
    
    TLorentzVector Smear_Electron(const TLorentzVector& Electron_MC_TLV);
    
    double GetMuonResolution(const TLorentzVector& Muon_MC_TLV);
    
    TLorentzVector Smear_Muon(const TLorentzVector& Muon_MC_TLV, double Resolution);
    
    TLorentzVector Smear_Muon(const TLorentzVector& Muon_MC_TLV);
    
    Vertex Smear_PV(Vertex PV);
    
    Vertex Smear_SV(Vertex SV);
    
    TVector3 Smear_PV(const TVector3& PV);
    
    TVector3 Smear_SV(const TVector3& SV);
    
    double Smear_ToF(double ToF);
    
    TVector3 Smear_Beta(Vertex User_PV, Vertex User_SV);
    
    TVector3 Smear_Beta_Mag(Vertex User_PV, Vertex User_SV);
    
    Vertex Smear_Vertex(Vertex User_Vertex, double sigmaV, double sigmaT);

};
#endif

#define Detector_cxx
inline Detector::Detector()
{
    con0_par=14.7467;
    con1_par=1.68788;
    con2_par=-2.95569e-07;
    con0_perp=15.4667;
    con1_perp=1.41597;
    con2_perp=-6.37947e-07;
    sigmaT=0.03;
    sigmaPV=12.0;
    sigmaSV=65.0;
}

inline void Detector::Set_con0_par(double User_con0_par)
{
    con0_par=User_con0_par;
}
inline double Detector::Get_con0_par()
{
    return con0_par;
}

inline void Detector::Set_con1_par(double User_con1_par)
{
    con1_par=User_con1_par;
}
inline double Detector::Get_con1_par()
{
    return con1_par;
}

inline void Detector::Set_con2_par(double User_con2_par)
{
    con2_par=User_con2_par;
}
inline double Detector::Get_con2_par()
{
    return con2_par;
}

inline void Detector::Set_con0_perp(double User_con0_perp)
{
    con0_perp = User_con0_perp;
}
inline double Detector::Get_con0_perp()
{
    return con0_perp;
}

inline void Detector::Set_con1_perp(double User_con1_perp)
{
    con1_perp = User_con1_perp;
}
inline double Detector::Get_con1_perp()
{
    return con1_perp;
}

inline void Detector::Set_con2_perp(double User_con2_perp)
{
    con2_perp = User_con2_perp;
}
inline double Detector::Get_con2_perp()
{
    return con2_perp;
}

inline double Detector::Get_Sigma_Perp(const TLorentzVector& sys)
{
    sigma_perp = (sys.Pt()*sqrt((con0_perp*con0_perp)/(sys.Pt()*sys.Pt())+(con1_perp*con1_perp)/sys.Pt()+con2_perp*con2_perp));
    return sigma_perp;
}

inline double Detector::Get_Sigma_Par(const TLorentzVector& sys)
{
    sigma_par = (sys.Pt()*sqrt((con0_par*con0_par)/(sys.Pt()*sys.Pt())+(con1_par*con1_par)/sys.Pt()+con2_par*con2_par));
    return sigma_par;
}

// smear MET
inline TVector3 Detector::Smear_MET(const TVector3& inv)
{
    TVector3 zhat;
    zhat.SetXYZ(0.0,0.0,1.0);
    
    TVector3 parhat  = inv.Unit();
    TVector3 perphat = inv.Cross(zhat).Unit();
    
    TVector3 MET = inv + gRandom->Gaus(0.,sigma_par)*parhat + gRandom->Gaus(0.,sigma_perp)*perphat;

    return MET;
}

inline TVector3 Detector::Smear_MET_Mag(TVector3& inv)
{
    TVector3 parhat  = inv.Unit();
    TVector3 MET_Mag = inv + gRandom->Gaus(0.,sigma_par)*parhat;
    return MET_Mag;
}

inline TVector3 Detector::Smear_MET_Dir(TVector3& inv)
{
    TVector3 zhat;
    zhat.SetXYZ(0.0,0.0,1.0);
    TVector3 perphat = inv.Cross(zhat).Unit();
    TVector3 MET_Dir = gRandom->Gaus(0.,sigma_perp)*perphat;
    return MET_Dir;
}

// smear Electron Given TLV of Electron
inline TLorentzVector Detector::Smear_Electron(const TLorentzVector& Electron_MC_TLV)
{
    double Smear_Electron_E;
    double Sigma = 0.0;
    
    if(fabs(Electron_MC_TLV.Eta()) <= 1.5)
    {
        Sigma = 0.028;
    }
    else if(fabs(Electron_MC_TLV.Eta()) > 1.5 && fabs(Electron_MC_TLV.Eta()) <= 1.75)
    {
        Sigma = 0.037;
    }
    else if(fabs(Electron_MC_TLV.Eta()) > 1.75 && fabs(Electron_MC_TLV.Eta()) <= 2.15)
    {
        Sigma = 0.038;
    }
    else if(fabs(Electron_MC_TLV.Eta()) > 2.15 && fabs(Electron_MC_TLV.Eta()) <= 3.00)
    {
        Sigma = 0.044;
    }
    else if(fabs(Electron_MC_TLV.Eta()) > 3.00 && fabs(Electron_MC_TLV.Eta()) <= 4.00)
    {
        Sigma = 0.10;
    }
    Smear_Electron_E = Electron_MC_TLV.E()+gRandom->Gaus(0.0,Sigma);
    //double Smear_Electron_PT = sqrt(Smear_Electron_E*Smear_Electron_E - Electron_MC_TLV.M()*Electron_MC_TLV.M())*TMath::Cosh(Electron_MC_TLV.Theta());
    TLorentzVector Electron_Smear;
    //Electron_Smear.SetPtEtaPhiE(Smear_Electron_PT,Electron_MC_TLV.Eta(),Electron_MC_TLV.Phi(),Smear_Electron_E);
    Electron_Smear.SetPtEtaPhiE(Smear_Electron_E/(TMath::CosH(Electron_MC_TLV.Eta())),Electron_MC_TLV.Eta(),Electron_MC_TLV.Phi(),Smear_Electron_E);
    return Electron_Smear;
}

double LogNormal(double mean, double sigma)
{
    double a, b;
    if(mean > 0.0)
    {
        b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
        a = TMath::Log(mean) - 0.5*b*b;
        
        return TMath::Exp(a + b*gRandom->Gaus(0.0, 1.0));
    }
    else
    {
        return 0.0;
    }
}

inline TLorentzVector Detector::Smear_Muon(const TLorentzVector& Muon_MC_TLV, double Resolution)
{
    TLorentzVector Muon_Smear;
    double Smear_Muon_PT;
    double PT_True = Muon_MC_TLV.Pt();
    double Eta_True = fabs(Muon_MC_TLV.Eta());
    //Smear_Muon_PT = PT_True*gRandom->Gaus(1.0,Resolution);
    Smear_Muon_PT = LogNormal(PT_True,Resolution*PT_True);
    Muon_Smear.SetPtEtaPhiE(Smear_Muon_PT,Muon_MC_TLV.Eta(),Muon_MC_TLV.Phi(),Smear_Muon_PT*TMath::CosH(Muon_MC_TLV.Eta()));
    return Muon_Smear;
}

inline TLorentzVector Detector::Smear_Muon(const TLorentzVector& Muon_MC_TLV)
{
    TLorentzVector Muon_Smear;
    double Smear_Muon_PT;
    double Resolution=GetMuonResolution(Muon_MC_TLV);
    double PT_True = Muon_MC_TLV.Pt();
    double Eta_True = fabs(Muon_MC_TLV.Eta());
    //Smear_Muon_PT = PT_True*gRandom->Gaus(1.0,Resolution);
    Smear_Muon_PT = LogNormal(PT_True,Resolution*PT_True);
    Muon_Smear.SetPtEtaPhiE(Smear_Muon_PT,Muon_MC_TLV.Eta(),Muon_MC_TLV.Phi(),Smear_Muon_PT*TMath::CosH(Muon_MC_TLV.Eta()));
    return Muon_Smear;
}

// smear Muon Given TLV of Muon
inline double Detector::GetMuonResolution(const TLorentzVector& Muon_MC_TLV)
{
    double Resolution = 0.0;
    double PT_True = Muon_MC_TLV.Pt();
    double Eta_True = fabs(Muon_MC_TLV.Eta());
    
    if(PT_True > 5.0 && PT_True < 10.0)
    {
        Resolution = 0.00519;
        if(Eta_True > 0.2 && Eta_True < 0.4) {Resolution+=0.00575;}
        else if(Eta_True > 0.4 && Eta_True < 0.6) {Resolution+=0.00629;}
        else if(Eta_True > 0.6 && Eta_True < 0.8) {Resolution+=0.00626;}
        else if(Eta_True > 0.8 && Eta_True < 1.0) {Resolution+=0.00722;}
        else if(Eta_True > 1.0 && Eta_True < 1.2) {Resolution+=0.00968;}
        else if(Eta_True > 1.2 && Eta_True < 1.4) {Resolution+=0.00983;}
        else if(Eta_True > 1.4 && Eta_True < 1.6) {Resolution+=0.0102;}
        else if(Eta_True > 1.6 && Eta_True < 1.8) {Resolution+=0.0126;}
        else if(Eta_True > 1.8 && Eta_True < 2.0) {Resolution+=0.0178;}
        else if(Eta_True > 2.0 && Eta_True < 2.2) {Resolution+=0.0248;}
        else if(Eta_True > 2.2 && Eta_True < 2.4) {Resolution+=0.0347;}
        else if(Eta_True > 2.4 && Eta_True < 2.5) {Resolution+=0.0359;}
        else if(Eta_True > 2.5 && Eta_True < 2.6) {Resolution+=0.0412;}
        else if(Eta_True > 2.6 && Eta_True < 2.7) {Resolution+=0.0422;}
        else if(Eta_True > 2.7 && Eta_True < 2.8) {Resolution+=0.0496;}
    }
    else if(PT_True > 10.0 && PT_True < 20.0)
    {
        Resolution = 0.00533;
        if(Eta_True > 0.2 && Eta_True < 0.4) {Resolution+=0.00592;}
        else if(Eta_True > 0.4 && Eta_True < 0.6) {Resolution+=0.00602;}
        else if(Eta_True > 0.6 && Eta_True < 0.8) {Resolution+=0.00614;}
        else if(Eta_True > 0.8 && Eta_True < 1.0) {Resolution+=0.00736;}
        else if(Eta_True > 1.0 && Eta_True < 1.2) {Resolution+=0.0104;}
        else if(Eta_True > 1.2 && Eta_True < 1.4) {Resolution+=0.00968;}
        else if(Eta_True > 1.4 && Eta_True < 1.6) {Resolution+=0.00991;}
        else if(Eta_True > 1.6 && Eta_True < 1.8) {Resolution+=0.0106;}
        else if(Eta_True > 1.8 && Eta_True < 2.0) {Resolution+=0.0158;}
        else if(Eta_True > 2.0 && Eta_True < 2.2) {Resolution+=0.0207;}
        else if(Eta_True > 2.2 && Eta_True < 2.4) {Resolution+=0.0327;}
        else if(Eta_True > 2.4 && Eta_True < 2.5) {Resolution+=0.0369;}
        else if(Eta_True > 2.5 && Eta_True < 2.6) {Resolution+=0.039;}
        else if(Eta_True > 2.6 && Eta_True < 2.7) {Resolution+=0.0488;}
        else if(Eta_True > 2.7 && Eta_True < 2.8) {Resolution+=0.0543;}
    }
    else if(PT_True > 20.0 && PT_True < 40.0)
    {
        Resolution = 0.00559;
        if(Eta_True > 0.2 && Eta_True < 0.4) {Resolution+=0.00586;}
        else if(Eta_True > 0.4 && Eta_True < 0.6) {Resolution+=0.00614;}
        else if(Eta_True > 0.6 && Eta_True < 0.8) {Resolution+=0.00654;}
        else if(Eta_True > 0.8 && Eta_True < 1.0) {Resolution+=0.00759;}
        else if(Eta_True > 1.0 && Eta_True < 1.2) {Resolution+=0.0103;}
        else if(Eta_True > 1.2 && Eta_True < 1.4) {Resolution+=0.00963;}
        else if(Eta_True > 1.4 && Eta_True < 1.6) {Resolution+=0.0102;}
        else if(Eta_True > 1.6 && Eta_True < 1.8) {Resolution+=0.011;}
        else if(Eta_True > 1.8 && Eta_True < 2.0) {Resolution+=0.0142;}
        else if(Eta_True > 2.0 && Eta_True < 2.2) {Resolution+=0.0197;}
        else if(Eta_True > 2.2 && Eta_True < 2.4) {Resolution+=0.0315;}
        else if(Eta_True > 2.4 && Eta_True < 2.5) {Resolution+=0.0409;}
        else if(Eta_True > 2.5 && Eta_True < 2.6) {Resolution+=0.0471;}
        else if(Eta_True > 2.6 && Eta_True < 2.7) {Resolution+=0.0579;}
        else if(Eta_True > 2.7 && Eta_True < 2.8) {Resolution+=0.0679;}
    }
    else if(PT_True > 40.0 && Eta_True < 2.8)
    {
        Resolution = 0.00656;
        if(Eta_True > 0.2 && Eta_True < 0.4) {Resolution+=0.007;}
        else if(Eta_True > 0.4 && Eta_True < 0.6) {Resolution+=0.00733;}
        else if(Eta_True > 0.6 && Eta_True < 0.8) {Resolution+=0.00752;}
        else if(Eta_True > 0.8 && Eta_True < 1.0) {Resolution+=0.00851;}
        else if(Eta_True > 1.0 && Eta_True < 1.2) {Resolution+=0.0105;}
        else if(Eta_True > 1.2 && Eta_True < 1.4) {Resolution+=0.0105;}
        else if(Eta_True > 1.4 && Eta_True < 1.6) {Resolution+=0.0111;}
        else if(Eta_True > 1.6 && Eta_True < 1.8) {Resolution+=0.0119;}
        else if(Eta_True > 1.8 && Eta_True < 2.0) {Resolution+=0.0165;}
        else if(Eta_True > 2.0 && Eta_True < 2.2) {Resolution+=0.0242;}
        else if(Eta_True > 2.2 && Eta_True < 2.4) {Resolution+=0.0359;}
        else if(Eta_True > 2.4 && Eta_True < 2.5) {Resolution+=0.0474;}
        else if(Eta_True > 2.5 && Eta_True < 2.6) {Resolution+=0.052;}
        else if(Eta_True > 2.6 && Eta_True < 2.7) {Resolution+=0.0639;}
        else if(Eta_True > 2.7 && Eta_True < 2.8) {Resolution+=0.0763;}
    }
    if(Eta_True >= 2.8 && Eta_True < 3.0)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.03253153;}
        else if(PT_True >= 1.0 && PT_True < 10.0) {Resolution = 0.032532 + (PT_True-1.0)*0.000382;}
        else if(PT_True >= 10.0 && PT_True < 100.0) {Resolution = 0.035969 + (PT_True-10.0)*0.000673;}
        else if(PT_True >= 100.0) {Resolution = (0.096568*PT_True)/100.0;}
    }
    else if(Eta_True >= 3.0 && Eta_True < 3.2)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.04551714;}
        else if(PT_True >= 1.0 && PT_True < 10.0) {Resolution = 0.045517 + (PT_True-1.0)*0.002188;}
        else if(PT_True >= 10.0 && PT_True < 100.0) {Resolution = 0.065207 + (PT_True-10.0)*0.003254;}
        else if(PT_True >= 100.0) {Resolution = (0.358089*PT_True)/100.0;}
    }
    else if(Eta_True >= 3.2 && Eta_True < 3.4)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.04724756;}
        else if(PT_True >= 1.0 && PT_True < 10.0) {Resolution = 0.047248 + (PT_True-1.0)*0.003322;}
        else if(PT_True >= 10.0 && PT_True < 100.0) {Resolution = 0.077142 + (PT_True-10.0)*0.004516;}
        else if(PT_True >= 100.0) {Resolution = (0.483608*PT_True)/100.0;}
    }
    else if(Eta_True >= 3.4 && Eta_True < 3.6)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.04493020;}
        else if(PT_True >= 1.0 && PT_True < 10.0) {Resolution = 0.044930 + (PT_True-1.0)*0.002917;}
        else if(PT_True >= 10.0 && PT_True < 100.0) {Resolution = 0.071183 + (PT_True-10.0)*0.004186;}
        else if(PT_True >= 100.0) {Resolution = (0.447915*PT_True)/100.0;}
    }
    else if(Eta_True >= 3.6 && Eta_True < 3.8)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.05464199;}
        else if(PT_True >= 1.0 && PT_True < 10.0) {Resolution = 0.054642 + (PT_True-1.0)*0.003456;}
        else if(PT_True >= 10.0 && PT_True < 100.0) {Resolution = 0.085748 + (PT_True-10.0)*0.004995;}
        else if(PT_True >= 100.0) {Resolution = (0.535262*PT_True)/100.0;}
    }
    else if(Eta_True >= 3.8 && Eta_True < 4.0)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.06551185;}
        else if(PT_True >= 1.0 && PT_True < 10.0) {Resolution = 0.065512 + (PT_True-1.0)*0.006574;}
        else if(PT_True >= 10.0 && PT_True < 100.0) {Resolution = 0.124674 + (PT_True-10.0)*0.009560;}
        else if(PT_True >= 100.0) {Resolution = (0.0985040*PT_True)/100.0;}
    }
    else if(Eta_True >= 0.0 && Eta_True < 0.2)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.00467469;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.004675 + (PT_True-1.0)*0.00056;}
    }
    else if(Eta_True >= 0.2 && Eta_True < 0.4)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.00515885;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.005159 + (PT_True-1.0)*0.00049;}
    }
    else if(Eta_True >= 0.4 && Eta_True < 0.6)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.00509775;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.005098 + (PT_True-1.0)*0.00055;}
    }
    else if(Eta_True >= 0.6 && Eta_True < 0.8)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.00568785;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.005688 + (PT_True-1.0)*0.00061;}
    }
    else if(Eta_True >= 0.8 && Eta_True < 1.0)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.00668287;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.006683 + (PT_True-1.0)*0.00065;}
    }
    else if(Eta_True >= 1.0 && Eta_True < 1.2)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.01047734;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.010477 + (PT_True-1.0)*0.00059;}
    }
    else if(Eta_True >= 1.2 && Eta_True < 1.4)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.01430653;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.014307 + (PT_True-1.0)*0.00038;}
    }
    else if(Eta_True >= 1.4 && Eta_True < 1.6)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.0171924;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.017192 + (PT_True-1.0)*0.00020;}
    }
    else if(Eta_True >= 1.6 && Eta_True < 1.8)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.01783940;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.017839 + (PT_True-1.0)*0.00045;}
    }
    else if(Eta_True >= 1.8 && Eta_True < 2.0)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.01787758;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.017878 + (PT_True-1.0)*0.00154;}
    }
    else if(Eta_True >= 2.0 && Eta_True < 2.2)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.01825549;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.018255 + (PT_True-1.0)*0.000270;}
    }
    else if(Eta_True >= 2.2 && Eta_True < 2.4)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.01803308;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.018033 + (PT_True-1.0)*0.00220;}
    }
    else if(Eta_True >= 2.4 && Eta_True < 2.6)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.02156195;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.0212562 + (PT_True-1.0)*0.000225;}
    }
    else if(Eta_True >= 2.6 && Eta_True < 2.8)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.02691276;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.026913 + (PT_True-1.0)*0.000357;}
    }
    else if(Eta_True >= 2.8 && Eta_True < 3.0)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.03253153;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.032532 + (PT_True-1.0)*0.00382;}
    }
    else if(Eta_True >= 3.0 && Eta_True < 3.2)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.04551714;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.045517 + (PT_True-1.0)*0.002188;}
    }
    else if(Eta_True >= 3.2 && Eta_True < 3.4)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.04724756;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.047248 + (PT_True-1.0)*0.003322;}
    }
    else if(Eta_True >= 3.4 && Eta_True < 3.6)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.04493020;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.044930 + (PT_True-1.0)*0.002917;}
    }
    else if(Eta_True >= 3.6 && Eta_True < 3.8)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.05464199;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.054642 + (PT_True-1.0)*0.003456;}
    }
    else if(Eta_True >= 3.8 && Eta_True < 4.0)
    {
        if(PT_True >= 0.0 && PT_True < 1.0) {Resolution = 0.06551185;}
        else if(PT_True >= 1.0 && PT_True < 5.0) {Resolution = 0.065512 + (PT_True-1.0)*0.006574;}
    }
    return Resolution;
}

//smear Primary Vertex
inline TVector3 Detector::Smear_PV(const TVector3& PV)
{
    TVector3 Smeared_PV;
    Smeared_PV.SetXYZ(gRandom->Gaus(PV.X(),sigmaPV),gRandom->Gaus(PV.Y(),sigmaPV),gRandom->Gaus(PV.Z(),sigmaPV));
    return Smeared_PV;
}

//smear Secondary Vertex
inline TVector3 Detector::Smear_SV(const TVector3& SV)
{
    TVector3 Smeared_SV;
    Smeared_SV.SetXYZ(gRandom->Gaus(SV.X(),sigmaSV),gRandom->Gaus(SV.Y(),sigmaSV),gRandom->Gaus(SV.Z(),sigmaSV));
    return Smeared_SV;
}

//smear Primary Vertex
inline Vertex Detector::Smear_PV(Vertex PV)
{
    Vertex Smeared_PV;
    Smeared_PV.SetXYZTPos(gRandom->Gaus(PV.GetXPos(),sigmaPV),gRandom->Gaus(PV.GetYPos(),sigmaPV),gRandom->Gaus(PV.GetZPos(),sigmaPV),gRandom->Gaus(PV.GetTPos(),sigmaT));
    return Smeared_PV;
}

//smear Secondary Vertex
inline Vertex Detector::Smear_SV(Vertex SV)
{
    Vertex Smeared_SV;
    Smeared_SV.SetXYZTPos(gRandom->Gaus(SV.GetXPos(),sigmaSV),gRandom->Gaus(SV.GetYPos(),sigmaSV),gRandom->Gaus(SV.GetZPos(),sigmaSV),gRandom->Gaus(SV.GetTPos(),sigmaT));
    return Smeared_SV;
}

inline double Detector::Smear_ToF(double ToF)
{
    return (gRandom->Gaus(ToF,sigmaT));
}

inline void Detector::Set_sigmaT(double User_sigmaT)
{
    sigmaT = User_sigmaT;
}

inline double Detector::Get_sigmaT()
{
    return sigmaT;
}

inline void Detector::Set_sigmaPV(double User_sigmaPV)
{
    sigmaPV = User_sigmaPV;
}

inline double Detector::Get_sigmaPV()
{
    return sigmaPV;
}

inline void Detector::Set_sigmaSV(double User_sigmaSV)
{
    sigmaSV = User_sigmaSV;
}

inline double Detector::Get_sigmaSV()
{
    return sigmaSV;
}

inline TVector3 Detector::Smear_Beta(Vertex User_PV, Vertex User_SV)
{
    TVector3 Smeared_Beta = (1./30./(User_SV.GetTPos()-User_PV.GetTPos()))*(User_SV.GetXYZPos()-User_PV.GetXYZPos());
    return Smeared_Beta;
}

inline TVector3 Detector::Smear_Beta_Mag(Vertex User_PV, Vertex User_SV)
{
    TVector3 Beta_True = (1./30./(User_SV.GetTPos()-User_PV.GetTPos()))*(User_SV.GetXYZPos()-User_PV.GetXYZPos());
    
    TVector3 RefVec(1.0,1.0,1.0);
    
    TVector3 parhat  = Beta_True.Unit();
    TVector3 perphat = Beta_True.Cross(RefVec).Unit();
    double sigma_par_beta = Beta_True.Mag()*.1;
    //double sigma_perp = ;
    
    //NOTE: Set divide PerpHat by 10 to turn down the direction smearing
    TVector3 Beta = Beta_True + gRandom->Gaus(0.,sigma_par_beta)*parhat;// + gRandom->Gaus(0.,sigma_perp)*perphat;
    
    return Beta;
}

inline Vertex Detector::Smear_Vertex(Vertex User_Vertex, double User_sigmaV, double User_sigmaT)
{
    Vertex Smeared;
    Smeared.SetXYZPos(gRandom->Gaus(User_Vertex.GetXPos(),User_sigmaV),gRandom->Gaus(User_Vertex.GetYPos(),User_sigmaV),gRandom->Gaus(User_Vertex.GetZPos(),User_sigmaV));
    Smeared.SetTPos(gRandom->Gaus(User_Vertex.GetTPos(),User_sigmaT));
    return Smeared;
}
