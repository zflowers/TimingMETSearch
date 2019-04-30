#ifndef RESOLUTION_H
#define RESOLUTION_H

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
    double Par_Resolution(TVector3 MET, TVector3 L1a, TVector3 L2a, TVector3 L1b, TVector3 L2b, double sigma_MET, double sigma_L1a, double sigma_L2a, double sigma_L1b, double sigma_L2b);
    double Visible_Resolution(TVector3 L1a, TVector3 L2a, TVector3 L1b, TVector3 L2b, double sigma_L1a, double sigma_L2a, double sigma_L1b, double sigma_L2b);
    double Mass_Parents2(TVector3 Par, TVector3 Betaa, TVector3 Betab);
    double Mass_Parents2(TVector3 MET_Mag, TVector3 MET_Dir, TVector3 Vis, TVector3 Betaa, TVector3 Betab);
    double Mass_Parents2_Resolution(TVector3 MET_Mag, TVector3 MET_Dir, TVector3 Vis, TVector3 Betaa, TVector3 Betab, double sigma_Beta_Maga, double sigma_MET_Mag, double sigma_MET_Dir, double sigma_Vis);
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
    
    double a = sigma_MET*(1.0/(gamma*Beta.Pt()))*TMath::Cos(Beta.DeltaPhi(MET));
    double b = sigma_L1*(1.0/(gamma*Beta.Pt()))*TMath::Cos(Beta.DeltaPhi(L1.Vect()));
    double c = sigma_L2*(1.0/(gamma*Beta.Pt()))*TMath::Cos(Beta.DeltaPhi(L2.Vect()));
    double d = sigma_Beta_Mag*(-gamma/(Beta.Pt()*Beta.Pt()*Beta.Mag()))*(Beta.Dot(MET+L1.Vect()+L2.Vect()));
    return sqrt(a*a + b*b + c*c + d*d);
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
    double E_Z_P = Energy_Z_Parent(Beta,L1,L2);
    double Mass_Invisible = sqrt(Mass_P*Mass_P-2.*Mass_P*E_Z_P+Mass_Visible*Mass_Visible);
    
    double dMp_dB = -J*gamma/(Beta.Pt()*Beta.Pt());
    double a = (sigma_Beta_Mag/Mass_Invisible)*(Mass_P*dMp_dB-E_Z_P*dMp_dB-Mass_P*gamma*gamma*gamma*(L1.E()*Beta.Mag()+L2.E()*Beta.Mag()-L1.E()*TMath::Cos(GetAngle(Beta,L1))-L2.E()*TMath::Cos(GetAngle(Beta,L2))));
    double dMp_dM = TMath::Cos(Beta.DeltaPhi(MET))/(gamma*Beta.Pt());
    double b = (sigma_MET/Mass_Invisible)*dMp_dM*(Mass_P-E_Z_P);
    
    return sqrt(a*a + b*b);
}

inline double Resolution::Par_Resolution(TVector3 MET, TVector3 L1a, TVector3 L2a, TVector3 L1b, TVector3 L2b, double sigma_MET, double sigma_L1a, double sigma_L2a, double sigma_L1b, double sigma_L2b)
{
    double Par = (MET+L1a+L2a+L1b+L2b).Mag();
    double a = (sigma_MET/Par)*(MET.Mag()+L1a.Mag()*TMath::Cos(GetAngle(MET,L1a))+L2a.Mag()*TMath::Cos(GetAngle(MET,L2a))+L1b.Mag()*TMath::Cos(GetAngle(MET,L1b))+L2b.Mag()*TMath::Cos(GetAngle(MET,L2b)));
    double b = (sigma_L1a/Par)*(L1a.Mag()+MET.Mag()*TMath::Cos(GetAngle(MET,L1a))+L2a.Mag()*TMath::Cos(GetAngle(L1a,L2a))+L1b.Mag()*TMath::Cos(GetAngle(L1a,L1b))+L2b.Mag()*TMath::Cos(GetAngle(L1a,L2b)));
    double c = (sigma_L2a/Par)*(L2a.Mag()+L1a.Mag()*TMath::Cos(GetAngle(L2a,L1a))+MET.Mag()*TMath::Cos(GetAngle(MET,L2a))+L1b.Mag()*TMath::Cos(GetAngle(L2a,L1b))+L2b.Mag()*TMath::Cos(GetAngle(L2a,L2b)));
    double d = (sigma_L1b/Par)*(L1b.Mag()+L1a.Mag()*TMath::Cos(GetAngle(L1b,L1a))+L2a.Mag()*TMath::Cos(GetAngle(L1b,L2a))+MET.Mag()*TMath::Cos(GetAngle(MET,L1b))+L2b.Mag()*TMath::Cos(GetAngle(L1b,L2b)));
    double e = (sigma_L2b/Par)*(L2b.Mag()+L1a.Mag()*TMath::Cos(GetAngle(L2b,L1a))+L2a.Mag()*TMath::Cos(GetAngle(L2b,L2a))+L1b.Mag()*TMath::Cos(GetAngle(L2b,L1b))+MET.Mag()*TMath::Cos(GetAngle(MET,L2b)));
    
    return sqrt(a*a+b*b+c*c+d*d+e*e);
}

inline double Resolution::Visible_Resolution(TVector3 L1a, TVector3 L2a, TVector3 L1b, TVector3 L2b, double sigma_L1a, double sigma_L2a, double sigma_L1b, double sigma_L2b)
{
    TVector3 Vis = L1a+L2a+L1b+L2b;
    double a = sigma_L1a*(L1a.Mag()+L1a.Dot(L2a)/L1a.Mag()+L1a.Dot(L1b)/L1a.Mag()+L1a.Dot(L2b)/L1a.Mag())/Vis.Mag();
    double b = sigma_L2a*(L2a.Mag()+L2a.Dot(L1a)/L2a.Mag()+L2a.Dot(L1b)/L2a.Mag()+L2a.Dot(L2b)/L2a.Mag())/Vis.Mag();
    double c = sigma_L1b*(L1b.Mag()+L1b.Dot(L2a)/L1b.Mag()+L1b.Dot(L1a)/L1b.Mag()+L1b.Dot(L2b)/L1b.Mag())/Vis.Mag();
    double d = sigma_L2b*(L2b.Mag()+L2b.Dot(L2a)/L2b.Mag()+L2b.Dot(L1b)/L2b.Mag()+L2b.Dot(L1a)/L2b.Mag())/Vis.Mag();
    return sqrt(a*a+b*b+c*c+d*d);
}

inline double Resolution::Mass_Parents2(TVector3 Par, TVector3 Betaa, TVector3 Betab)
{
    TVector3 Zhat(0.0,0.0,1.0);
    double gamma = 1.0/sqrt(1.0-Betaa.Mag2());
    Betaa.SetZ(0.);
    Betab.SetZ(0.);
    return (Betab.Unit().Cross(Par)).Dot(Zhat)*Par.Mag2()/(gamma*Betaa.Pt()*((Betaa.Unit().Dot(Par)*Betab.Unit().Cross(Par))-Betab.Unit().Dot(Par)*Betaa.Unit().Cross(Par)).Dot(Zhat));
}

inline double Resolution::Mass_Parents2(TVector3 MET_Mag, TVector3 MET_Dir, TVector3 Vis, TVector3 Betaa, TVector3 Betab)
{
    TVector3 Zhat(0.0,0.0,1.0);
    double gamma = 1.0/sqrt(1.0-Betaa.Mag2());
    TVector3 BetaaT = Betaa;
    TVector3 BetabT = Betab;
    BetaaT.SetZ(0.);
    BetabT.SetZ(0.);
    //new angles
    //m: MET_Mag
    //d: MET_Dir
    //v: Vis
    
    double cmv = MET_Mag.Unit().Dot(Vis.Unit());
    double cdv = MET_Dir.Unit().Dot(Vis.Unit());
    double sbm = MET_Mag.Unit().Cross(BetabT.Unit()).Dot(Zhat);
    double sbd = MET_Dir.Unit().Cross(BetabT.Unit()).Dot(Zhat);
    double sbv = Vis.Unit().Cross(BetabT.Unit()).Dot(Zhat);
    double cam = MET_Mag.Unit().Dot(BetaaT.Unit());
    double cad = MET_Dir.Unit().Dot(BetaaT.Unit());
    double cav = Vis.Unit().Dot(BetaaT.Unit());
    double cbm = MET_Mag.Unit().Dot(BetabT.Unit());
    double cbd = MET_Dir.Unit().Dot(BetabT.Unit());
    double cbv = Vis.Unit().Dot(BetabT.Unit());
    double sam = MET_Mag.Unit().Cross(BetaaT.Unit()).Dot(Zhat);
    double sad = MET_Dir.Unit().Cross(BetaaT.Unit()).Dot(Zhat);
    double sav = Vis.Unit().Cross(BetaaT.Unit()).Dot(Zhat);

    
    double den1 = MET_Mag.Mag2()*(cam*sbm-cbm*sam);
    double den2 = MET_Dir.Mag2()*(cad*sbd-cbd*sad);
    double den3 = Vis.Mag2()*(cav*sbv-cbv*sav);
    double den4 = MET_Mag.Mag()*MET_Dir.Mag()*(cam*sbd+cad*sbm-cbm*sad-cbd*sam);
    double den5 = MET_Mag.Mag()*Vis.Mag()*(cam*sbv+cav*sbm-cbm*sav-cbv*sam);
    double den6 = MET_Dir.Mag()*Vis.Mag()*(cad*sbv+cav*sbd-cbd*sav-cbv*sad);
    
    double den = gamma*Betaa.Pt()*(den1+den2+den3+den4+den5+den6);

    double num = MET_Mag.Mag2()*MET_Mag.Mag()*sbm+MET_Mag.Mag2()*MET_Dir.Mag()*sbd+2.*MET_Mag.Mag2()*Vis.Mag()*cmv*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*Vis.Mag()*cdv*sbm+2.*MET_Mag.Mag()*Vis.Mag2()*cmv*sbv+MET_Mag.Mag()*Vis.Mag2()*sbm+MET_Dir.Mag2()*MET_Dir.Mag()*sbd+2.*MET_Dir.Mag2()*Vis.Mag()*cdv*sbd+MET_Dir.Mag2()*Vis.Mag()*sbv+2.*MET_Dir.Mag()*Vis.Mag2()*cdv*sbv+MET_Dir.Mag()*Vis.Mag2()*sbd+Vis.Mag2()*Vis.Mag()*sbv+MET_Mag.Mag2()*Vis.Mag()*sbv+MET_Mag.Mag()*MET_Dir.Mag2()*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*Vis.Mag()*cmv*sbd;
    
    return num/den;
}

inline double Resolution::Mass_Parents2_Resolution(TVector3 MET_Mag, TVector3 MET_Dir, TVector3 Vis, TVector3 Betaa, TVector3 Betab, double sigma_Beta_Maga, double sigma_MET_Mag, double sigma_MET_Dir, double sigma_Vis)
{
    TVector3 Zhat(0.0,0.0,1.0);
    double gamma = 1.0/sqrt(1.0-Betaa.Mag2());
    TVector3 BetaaT = Betaa;
    TVector3 BetabT = Betab;
    BetaaT.SetZ(0.);
    BetabT.SetZ(0.);
    //Angles
    //m: MET_Mag
    //d: MET_Dir
    //v: Vis
    
    double cmv = MET_Mag.Unit().Dot(Vis.Unit());
    double cdv = MET_Dir.Unit().Dot(Vis.Unit());
    double sbm = MET_Mag.Unit().Cross(BetabT.Unit()).Dot(Zhat);
    double sbd = MET_Dir.Unit().Cross(BetabT.Unit()).Dot(Zhat);
    double sbv = Vis.Unit().Cross(BetabT.Unit()).Dot(Zhat);
    double cam = MET_Mag.Unit().Dot(BetaaT.Unit());
    double cad = MET_Dir.Unit().Dot(BetaaT.Unit());
    double cav = Vis.Unit().Dot(BetaaT.Unit());
    double cbm = MET_Mag.Unit().Dot(BetabT.Unit());
    double cbd = MET_Dir.Unit().Dot(BetabT.Unit());
    double cbv = Vis.Unit().Dot(BetabT.Unit());
    double sam = MET_Mag.Unit().Cross(BetaaT.Unit()).Dot(Zhat);
    double sad = MET_Dir.Unit().Cross(BetaaT.Unit()).Dot(Zhat);
    double sav = Vis.Unit().Cross(BetaaT.Unit()).Dot(Zhat);
    
    double num = MET_Mag.Mag2()*MET_Mag.Mag()*sbm+MET_Mag.Mag2()*MET_Dir.Mag()*sbd+2.*MET_Mag.Mag2()*Vis.Mag()*cmv*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*Vis.Mag()*cdv*sbm+2.*MET_Mag.Mag()*Vis.Mag2()*cmv*sbv+MET_Mag.Mag()*Vis.Mag2()*sbm+MET_Dir.Mag2()*MET_Dir.Mag()*sbd+2.*MET_Dir.Mag2()*Vis.Mag()*cdv*sbd+MET_Dir.Mag2()*Vis.Mag()*sbv+2.*MET_Dir.Mag()*Vis.Mag2()*cdv*sbv+MET_Dir.Mag()*Vis.Mag2()*sbd+Vis.Mag2()*Vis.Mag()*sbv+MET_Mag.Mag2()*Vis.Mag()*sbv+MET_Mag.Mag()*MET_Dir.Mag2()*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*Vis.Mag()*cmv*sbd;
    double den1 = MET_Mag.Mag2()*(cam*sbm-cbm*sam);
    double den2 = MET_Dir.Mag2()*(cad*sbd-cbd*sad);
    double den3 = Vis.Mag2()*(cav*sbv-cbv*sav);
    double den4 = MET_Mag.Mag()*MET_Dir.Mag()*(cam*sbd+cad*sbm-cbm*sad-cbd*sam);
    double den5 = MET_Mag.Mag()*Vis.Mag()*(cam*sbv+cav*sbm-cbm*sav-cbv*sam);
    double den6 = MET_Dir.Mag()*Vis.Mag()*(cad*sbv+cav*sbd-cbd*sav-cbv*sad);
    double den = (den1+den2+den3+den4+den5+den6);

    double dnum_dMag = 3.*MET_Mag.Mag2()*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*sbd+4.*MET_Mag.Mag()*Vis.Mag()*cmv*sbm+MET_Dir.Mag2()*sbm+2.*MET_Dir.Mag()*Vis.Mag()*cmv*sbd+2.*MET_Mag.Mag()*Vis.Mag()*sbv+2.*MET_Dir.Mag()*Vis.Mag()*cdv*sbm+2.*Vis.Mag2()*cmv*sbv+Vis.Mag2()*sbm;
    double dden_dMag = 2.*MET_Mag.Mag()*(cam*sbm-cbm*sam)+MET_Dir.Mag()*(cam*sbd+cad*sbm-cbm*sad-cbd*sam)+Vis.Mag()*(cam*sbv+cav*sbm-cbm*sav-cbv*sam);
    double dnum_dDir = 3.*MET_Dir.Mag2()*sbd+MET_Mag.Mag2()*sbd+2.*MET_Mag.Mag()*Vis.Mag()*cdv*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*sbm+2.*MET_Mag.Mag()*Vis.Mag()*cmv*sbd+4.*MET_Dir.Mag()*Vis.Mag()*cdv*sbd+2.*MET_Dir.Mag()*Vis.Mag()*sbv+2.*Vis.Mag2()*cdv*sbv+Vis.Mag2()*sbd;
    double dden_dDir = 2.*MET_Dir.Mag()*(cad*sbd-cbd*sad)+MET_Mag.Mag()*(cam*sbd+cad*sbm-cbm*sad-cbd*sam)+Vis.Mag()*(cad*sbv+cav*sbd-cbd*sav-cbv*sad);
    double dnum_dVis = 2.*MET_Mag.Mag2()*cmv*sbm+2.*MET_Mag.Mag()*MET_Dir.Mag()*cdv*sbm+4.*MET_Mag.Mag()*Vis.Mag()*cmv*sbv+2.*MET_Mag.Mag()*Vis.Mag()*sbm+2.*MET_Dir.Mag2()*cdv*sbd+MET_Dir.Mag2()*sbv+4.*MET_Dir.Mag()*Vis.Mag()*cdv*sbv+2.*MET_Dir.Mag()*Vis.Mag()*sbd+3.*Vis.Mag2()*sbv+MET_Mag.Mag2()*sbv+2.*MET_Mag.Mag()*MET_Dir.Mag()*cmv*sbd;
    double dden_dVis = 2.*Vis.Mag()*(cav*sbv-cbv*sav)+MET_Mag.Mag()*(cam*sbv+cav*sbm-cbm*sav-cbv*sam)+MET_Dir.Mag()*(cad*sbv+cav*sbd-cbd*sav-cbv*sad);
    
    double f_MET_Mag = (1./num)*(dnum_dMag)-(1./den)*(dden_dMag);
    double f_MET_Dir = (1./num)*(dnum_dDir)-(1./den)*(dden_dDir);
    double f_Vis = (1./num)*(dnum_dVis)-(1./den)*(dden_dVis);
    
    double MET_MAG_RES = sigma_MET_Mag*Mass_Parents2(MET_Mag,MET_Dir,Vis,Betaa,Betab)*f_MET_Mag;
    double MET_DIR_RES = sigma_MET_Dir*Mass_Parents2(MET_Mag,MET_Dir,Vis,Betaa,Betab)*f_MET_Dir;
    double VIS_RES = sigma_Vis*Mass_Parents2(MET_Mag,MET_Dir,Vis,Betaa,Betab)*f_Vis;
    double BETA_RES = sigma_Beta_Maga*Mass_Parents2(MET_Mag,MET_Dir,Vis,Betaa,Betab)/Betaa.Mag();
    
    return sqrt(BETA_RES*BETA_RES + MET_MAG_RES*MET_MAG_RES + MET_DIR_RES*MET_DIR_RES + VIS_RES*VIS_RES);
}
