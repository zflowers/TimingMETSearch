#ifndef RESOLUTION_H
#define RESOLUTION_H
//#include "Physics.hh"
//#include "Detector.hh"

class Resolution
{
    
private:
    
    Detector m_Detector;

public:
    Resolution();
    Resolution(const Resolution& old_resolution);
    ~Resolution();
    Resolution(Detector u_Detector);
    double GetAngle(TVector3 TLV1, TVector3 TLV2);
    double GetAngle(TLorentzVector TLV1, TLorentzVector TLV2);
    double GetAngle(TVector3 TLV1, TLorentzVector TLV2);
    double GetAngle(TLorentzVector TLV1, TVector3 TLV2);
    double GetAngleError(TVector3 TLV1, TVector3 TLV2, TVector3 TLV1_RECO, TVector3 TLV2_RECO);
    TVector3 GetBeta(Vertex PV, Vertex SV);
    double Energy_Z_Parent(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2);
    double Energy_Z_Parent(Vertex PV, Vertex SV, TLorentzVector Lepton1, TLorentzVector Lepton2);
    double Energy_Z_Parent_Resolution(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2, double sigma_E_lepton1, double sigma_E_lepton2, double sigma_cos_theta1, double sigma_cos_theta2, double sigma_Beta_Mag);
    double Energy_Z_Parent_Resolution(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2, double sigma_E_lepton, double sigma_cos_theta, double sigma_Beta_Mag);
    double Cos_Resolution(TVector3 V1, TVector3 V2, double sigma_V1, double sigma_V2);
    double Mass_Parent(TVector3 Decay, TVector3 Beta);
    double Mass_Parent_Resolution(TVector3 Beta, TVector3 MET, TLorentzVector L1, TLorentzVector L2, double sigma_MET, double sigma_L1, double sigma_L2, double sigma_Beta_Mag);
    double Mass_Invisible_Resolution(TVector3 Beta, TVector3 MET, TLorentzVector L1, TLorentzVector L2, double sigma_MET, double sigma_Beta_Mag);
};
#endif

#define Resolution_cxx
inline Resolution::Resolution()
{
    m_Detector.Set_con0_par(14.7467);
    m_Detector.Set_con1_par(1.68788);
    m_Detector.Set_con2_par(-2.95569e-07);
    m_Detector.Set_con0_perp(15.4667);
    m_Detector.Set_con1_perp(1.41597);
    m_Detector.Set_con2_perp(-6.37947e-07);
    m_Detector.Set_sigmaT(0.03/sqrt(2.));
    m_Detector.Set_sigmaPV(20.0/10000.0);
    m_Detector.Set_sigmaSV(65.0/10000.0);
}

inline Resolution::Resolution(Detector u_Detector)
{
    m_Detector.Set_con0_par(u_Detector.Get_con0_par());
    m_Detector.Set_con1_par(u_Detector.Get_con1_par());
    m_Detector.Set_con2_par(u_Detector.Get_con2_par());
    m_Detector.Set_con0_perp(u_Detector.Get_con0_perp());
    m_Detector.Set_con1_perp(u_Detector.Get_con1_perp());
    m_Detector.Set_con2_perp(u_Detector.Get_con2_perp());
    m_Detector.Set_sigmaT(u_Detector.Get_sigmaT());
    m_Detector.Set_sigmaPV(u_Detector.Get_sigmaPV());
    m_Detector.Set_sigmaSV(u_Detector.Get_sigmaSV());
}

inline Resolution::Resolution(const Resolution& old_resolution)
{
    m_Detector = old_resolution.m_Detector;
}

inline Resolution::~Resolution()
{
}

inline double Resolution::GetAngle(TVector3 TLV1, TVector3 TLV2){
    return TMath::ACos((TLV1.Dot(TLV2))/(TLV1.Mag()*TLV2.Mag()));
}

inline double Resolution::GetAngle(TLorentzVector TLV1, TLorentzVector TLV2){
    return GetAngle(TLV1.Vect(),TLV2.Vect());
}

inline double Resolution::GetAngle(TVector3 TLV1, TLorentzVector TLV2){
    return GetAngle(TLV1,TLV2.Vect());
}

inline double Resolution::GetAngle(TLorentzVector TLV1, TVector3 TLV2){
    return GetAngle(TLV1.Vect(),TLV2);
}

inline double Resolution::GetAngleError(TVector3 TLV1, TVector3 TLV2, TVector3 TLV1_RECO, TVector3 TLV2_RECO){
    double costheta_true = (TLV1.Dot(TLV2))/(TLV1.Mag()*TLV2.Mag());
    double costheta_reco = (TLV1_RECO.Dot(TLV2_RECO))/(TLV1_RECO.Mag()*TLV2_RECO.Mag());
    return costheta_true-costheta_reco;
}

inline TVector3 Resolution::GetBeta(Vertex PV, Vertex SV)
{
    return (1.0/30.0/(SV.GetTPos()-PV.GetTPos())*(SV.GetXYZPos()-PV.GetXYZPos()));
}

inline double Resolution::Energy_Z_Parent(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2)
{
    double gamma = 1.0/sqrt(1.0-Beta.Mag2());
    TLorentzVector Z = Lepton1 + Lepton2;
    return gamma*(Z.E() - Z.Vect().Dot(Beta));
}

inline double Resolution::Energy_Z_Parent(Vertex PV, Vertex SV, TLorentzVector Lepton1, TLorentzVector Lepton2)
{
    TVector3 Beta = GetBeta(PV,SV);
    double gamma = 1.0/sqrt(1.0-Beta.Mag2());
    TLorentzVector Z = Lepton1 + Lepton2;
    return gamma*(Z.E() - Z.Vect().Dot(Beta));
}

inline double Resolution::Energy_Z_Parent_Resolution(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2, double sigma_E_lepton1, double sigma_E_lepton2, double sigma_cos_theta1, double sigma_cos_theta2, double sigma_Beta_Mag)
{
    double gamma = 1.0/sqrt(1.0-Beta.Mag2());
    double gamma3 = gamma*gamma*gamma;
    double a = sigma_E_lepton1*gamma*(1.0-Beta.Mag()*TMath::Cos(GetAngle(Beta,Lepton1)));
    double b = sigma_E_lepton2*gamma*(1.0-Beta.Mag()*TMath::Cos(GetAngle(Beta,Lepton2)));
    double c = sigma_cos_theta1*gamma3*Lepton1.E()*TMath::Cos(GetAngle(Beta,Lepton1));
    double d = sigma_cos_theta2*gamma3*Lepton2.E()*TMath::Cos(GetAngle(Beta,Lepton2));
    double e = sigma_Beta_Mag*gamma3*(Lepton1.E()*Beta.Mag()+Lepton2.E()*Beta.Mag()-Lepton1.E()*TMath::Cos(GetAngle(Beta,Lepton1))-Lepton2.E()*TMath::Cos(GetAngle(Beta,Lepton2)));
    return sqrt(a*a + b*b + e*e);
}

inline double Resolution::Cos_Resolution(TVector3 V1, TVector3 V2, double sigma_V1, double sigma_V2)
{
    double a = sigma_V1*V2.Mag2()/(sqrt((V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())))*(V1.X()*(V1.Y()*V2.Y()+V1.Z()*V2.Z())-V2.X()*(V1.Y()*V1.Y()+V2.Z()*V2.Z()));
    double b = sigma_V1*V2.Mag2()/(sqrt((V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())))*(V1.Y()*(V1.Z()*V2.Z()+V1.X()*V2.X())-V2.Y()*(V1.Z()*V1.Z()+V2.X()*V2.X()));
    double c = sigma_V1*V2.Mag2()/(sqrt((V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())))*(V1.Z()*(V1.X()*V2.X()+V1.Y()*V2.Y())-V2.Z()*(V1.X()*V1.X()+V2.Y()*V2.Y()));
    double d = sigma_V2*V1.Mag2()/(sqrt((V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())))*(V1.X()*(V2.Y()*V2.Y()+V2.Z()*V2.Z())-V2.X()*(V1.Y()*V2.Y()+V1.Z()*V2.Z()));
    double e = sigma_V2*V1.Mag2()/(sqrt((V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())))*(V1.Y()*(V2.Z()*V2.Z()+V2.X()*V2.X())-V2.Y()*(V1.Z()*V2.Z()+V1.X()*V2.X()));
    double f = sigma_V2*V1.Mag2()/(sqrt((V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())*(V1.Mag2()*V2.Mag2())))*(V1.Z()*(V2.X()*V2.X()+V2.Y()*V2.Y())-V2.Z()*(V1.X()*V2.X()+V1.Y()*V2.Y()));
    return sqrt(a*a+b*b+c*c+d*d+e*e+f*f);
}

inline double Resolution::Mass_Parent(TVector3 Decay, TVector3 Beta)
{
    return Beta.Dot(Decay)/(1.0/sqrt(1.0-Beta.Mag2())*Beta.Pt()*Beta.Pt());
}

inline double Resolution::Mass_Parent_Resolution(TVector3 Beta, TVector3 MET, TLorentzVector L1, TLorentzVector L2, double sigma_MET, double sigma_L1, double sigma_L2, double sigma_Beta_Mag){
    
    double gamma = 1.0/sqrt(1.0-Beta.Mag2());
    //Beta.SetZ(0.0);
    double sigmaCos_BetaMET = Cos_Resolution(Beta,MET,sigma_Beta_Mag,sigma_MET);
    double sigmaCos_BetaL1 = Cos_Resolution(Beta,L1.Vect(),sigma_Beta_Mag,sigma_L1);
    double sigmaCos_BetaL2 = Cos_Resolution(Beta,L2.Vect(),sigma_Beta_Mag,sigma_L2);
    
    double a = sigma_MET*(1.0/(gamma*Beta.Pt()))*TMath::Cos(Beta.DeltaPhi(MET));
    double b = sigma_L1*(1.0/(gamma*Beta.Pt()))*TMath::Cos(Beta.DeltaPhi(L1.Vect()));
    double c = sigma_L2*(1.0/(gamma*Beta.Pt()))*TMath::Cos(Beta.DeltaPhi(L2.Vect()));
    double d = sigma_Beta_Mag*(-gamma/(Beta.Pt()*Beta.Pt()*Beta.Mag()))*(Beta.Dot(MET+L1.Vect()+L2.Vect()));
    double e = sigmaCos_BetaMET*MET.Mag()/(gamma*Beta.Mag());
    double f = sigmaCos_BetaL1*L1.Vect().Mag()/(gamma*Beta.Mag());
    double g = sigmaCos_BetaL2*L2.Vect().Mag()/(gamma*Beta.Mag());
    
    return sqrt(a*a + b*b + c*c + d*d);
    //return sqrt(a*a + b*b + c*c + d*d + e*e + f*f + g*g);
}

inline double Resolution::Mass_Invisible_Resolution(TVector3 Beta, TVector3 MET, TLorentzVector L1, TLorentzVector L2, double sigma_MET, double sigma_Beta_Mag)
{
    TLorentzVector L1t = L1;
    TLorentzVector L2t = L2;
    L1t.SetZ(0.0);
    L2t.SetZ(0.0);
    double J = Beta.Unit().Dot(MET+L1t.Vect()+L2t.Vect());
    double gamma = 1.0/sqrt(1.0-Beta.Mag2());
    double Mass_Visible = (L1+L2).M();
    double Mass_P = Mass_Parent(L1t.Vect()+L2t.Vect()+MET,Beta);
    double Mass_Invisible = sqrt(Mass_P*Mass_P-2.0*Mass_P*Energy_Z_Parent(Beta,L1,L2)+Mass_Visible*Mass_Visible);
    
    double a = (sigma_Beta_Mag/Mass_Invisible)*((-J*gamma)/(Beta.Pt()*Beta.Pt())*(Mass_P-Energy_Z_Parent(Beta,L1,L2))+Mass_P*gamma*gamma*gamma*(Beta.Mag()*L1.E()+Beta.Mag()*L2.E()-L1.E()*TMath::Cos(GetAngle(Beta,L1))-L2.E()*TMath::Cos(GetAngle(Beta,L2))));
    double b = (sigma_MET/Mass_Invisible)*(TMath::Cos(Beta.DeltaPhi(MET)))/(gamma*Beta.Pt())*(Mass_P-Energy_Z_Parent(Beta,L1,L2));
    
    return sqrt(a*a + b*b);
}
