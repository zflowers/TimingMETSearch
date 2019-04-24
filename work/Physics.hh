#ifndef PHYSICS_H
#define PHYSICS_H
#include <TRandom.h>
#include <TH3.h>
#include "TLorentzVector.h"
#include "TVector3.h"

class Physics
{
    
private:
    
    TH3* m_histPtvsEta;

public:
    Physics();
    Physics(const Physics& old_physics);
    ~Physics();
    void SetEtaPtMCM(const TH3& hist);
    void GetEtaPtMCM(double& PT, double& eta, double& M);
    Vertex Get_PV(double User_PV_X, double User_PV_Y, double User_PV_Z, double User_PV_T);
    Vertex Get_SV(double User_ctau, TLorentzVector P);
    double Get_Beta(TLorentzVector P);
    TVector3 GetBeta(Vertex PV, Vertex SV);
    TVector3 Get_vBeta(double User_ToF, TLorentzVector P);
    double Get_Gamma(TLorentzVector P);
    double Get_Gamma(double Beta_Mag);
    double Get_ToF(double User_ctau, TLorentzVector& P);
};
#endif

#define Physics_cxx
inline Physics::Physics()
{
    m_histPtvsEta = nullptr;
}

inline Physics::Physics(const Physics& old_physics)
{
    m_histPtvsEta = (TH3*)old_physics.m_histPtvsEta->Clone("copy");
}

inline Physics::~Physics()
{
    if(m_histPtvsEta != nullptr)
        delete m_histPtvsEta;
}

inline void Physics::SetEtaPtMCM(const TH3& hist)
{
    string name = "m_";
    string histname = hist.GetName();
    name+=histname;
    m_histPtvsEta=(TH3*)hist.Clone(name.c_str());
    m_histPtvsEta->SetDirectory(0);
}

inline void Physics::GetEtaPtMCM(double& eta, double& PT, double& M)
{
    m_histPtvsEta->GetRandom3(eta,PT,M);
}

inline Vertex Physics::Get_PV(double User_PV_X, double User_PV_Y, double User_PV_Z, double User_PV_T)
{
    Vertex PV;
    PV.SetXYZTPos(User_PV_X,User_PV_Y,User_PV_Z,User_PV_T);
    return PV;
}

inline Vertex Physics::Get_SV(double User_ToF, TLorentzVector P)
{
    TVector3 vBeta = P.BoostVector();
    vBeta = 30.*User_ToF*vBeta;
    Vertex SV(vBeta.X(),vBeta.Y(),vBeta.Z(),User_ToF);
    return SV;
}

inline double Physics::Get_Beta(TLorentzVector P)
{
    return (P.BoostVector().Mag());
}

//get the velocity from two 4D Vertices
inline TVector3 Physics::GetBeta(Vertex PV, Vertex SV)
{
    return (1.0/30.0/(SV.GetTPos()-PV.GetTPos())*(SV.GetXYZPos()-PV.GetXYZPos()));
}

inline TVector3 Physics::Get_vBeta(double User_ToF, TLorentzVector P)
{
    TVector3 vBeta;
    vBeta = P.BoostVector();
    vBeta = User_ToF*30.*vBeta;
    return vBeta;
}

inline double Physics::Get_Gamma(TLorentzVector P)
{
    double beta = P.BoostVector().Mag();
    return (1./sqrt(1.-beta*beta));
}

inline double Physics::Get_Gamma(double Beta_Mag)
{
    return (1./sqrt(1.-Beta_Mag*Beta_Mag));
}

inline double Physics::Get_ToF(double User_ctau, TLorentzVector& P)
{
    double beta = P.BoostVector().Mag();
    return -1.*1./sqrt(1.-beta*beta)*User_ctau/30.*log((1.-gRandom->Rndm()));
}

