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
    double Energy_Z_Parent(Vertex PV, Vertex SV, TLorentzVector Lepton1, TLorentzVector Lepton2);
    double Energy_Z_Parent_Resolution(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2, double sigma_E_lepton1, double sigma_E_lepton2, double sigma_cos_theta1, double sigma_cos_theta2, double sigma_Beta_Mag);
    double Energy_Z_Parent_Resolution(TVector3 Beta, TLorentzVector Lepton1, TLorentzVector Lepton2, double sigma_E_lepton, double sigma_cos_theta, double sigma_Beta_Mag);
    double Cos_Resolution(TVector3 V1, TVector3 V2, double sigma_V1, double sigma_V2);
    double Mass_Parent_Resolution(TVector3 Beta, TVector3 MET, TLorentzVector L1, TLorentzVector L2, double sigma_MET, double sigma_L1, double sigma_L2, double sigma_Beta_MagT);
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
    return sqrt(a*a + b*b + c*c + d*d + e*e);
}

inline double Resolution::Cos_Resolution(TVector3 V1, TVector3 V2, double sigma_V1, double sigma_V2)
{
    
    double a1 = sigma_V1*V2.X();
    double a2 = sigma_V1*V2.Y();
    double a3 = sigma_V2*V1.X();
    double a4 = sigma_V2*V1.Y();
    double a = (1.0/(V1.Mag()*V2.Mag()))*sqrt(a1*a1+a2*a2+a3*a3+a4*a4);
    double b = sigma_V1*V1.Dot(V2)*(1.0/(V1.Mag2()*V2.Mag()));
    double c = sigma_V2*V1.Dot(V2)*(1.0/(V1.Mag()*V2.Mag2()));
    
    
    
    return sqrt(a*a+b*b+c*c);
}

inline double Resolution::Mass_Parent_Resolution(TVector3 Beta, TVector3 MET, TLorentzVector L1, TLorentzVector L2, double sigma_MET, double sigma_L1, double sigma_L2, double sigma_Beta_MagT){
    double gamma = 1.0/sqrt(1.0-Beta.Mag2());
    Beta.SetZ(0.0);
    double sigmaCos_BetaMET = Cos_Resolution(Beta,MET,sigma_Beta_MagT,sigma_MET);
    double sigmaCos_BetaL1 = Cos_Resolution(Beta,L1.Vect(),sigma_Beta_MagT,sigma_L1);
    double sigmaCos_BetaL2 = Cos_Resolution(Beta,L2.Vect(),sigma_Beta_MagT,sigma_L2);
    
    double a = sigma_MET*(1.0/(gamma*Beta.Mag()))*TMath::Cos(GetAngle(Beta,MET));
    double b = sigma_L1*(1.0/(gamma*Beta.Mag()))*TMath::Cos(GetAngle(Beta,L1.Vect()));
    double c = sigma_L2*(1.0/(gamma*Beta.Mag()))*TMath::Cos(GetAngle(Beta,L2.Vect()));
    double d = sigma_Beta_MagT*(-gamma/(Beta.Mag()*Beta.Mag()))*(Beta.Unit().Dot(MET+L1.Vect()+L2.Vect()));
    double e = sigmaCos_BetaMET*MET.Mag()/(gamma*Beta.Mag());
    double f = sigmaCos_BetaL1*L1.Vect().Mag()/(gamma*Beta.Mag());
    double g = sigmaCos_BetaL2*L2.Vect().Mag()/(gamma*Beta.Mag());
    
    return sqrt(a*a + b*b + c*c + d*d);
    //return sqrt(a*a + b*b + c*c + d*d + e*e + f*f + g*g);
}
