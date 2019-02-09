#ifndef VERTEX_H
#define VERTEX_H
#include <TVector3.h>

class Vertex
{
private:
    TVector3 m_positionXYZ;
    double m_time;
    
public:
    Vertex();
    Vertex(double X, double Y, double Z, double T);
    void SetXPos(double X);
    void SetYPos(double Y);
    void SetZPos(double Z);
    void SetTPos(double T);
    double GetXPos();
    double GetYPos();
    double GetZPos();
    double GetTPos();
    TVector3 GetXYZPos();
    void SetXYZPos(double X, double Y, double Z);
    void SetXYZTPos(double X, double Y, double Z, double T);
};

#endif

#define Vertex_cxx
inline Vertex::Vertex()
{
    m_positionXYZ.SetXYZ(0.0,0.0,0.0);
    m_time = 0.0;
}

inline Vertex::Vertex(double X, double Y, double Z, double T)
{
    m_positionXYZ.SetXYZ(X,Y,Z);
    m_time = T;
}

inline void Vertex::SetXPos(double X)
{
    m_positionXYZ.SetX(X);
}

inline void Vertex::SetYPos(double Y)
{
    m_positionXYZ.SetY(Y);
}

inline void Vertex::SetZPos(double Z)
{
    m_positionXYZ.SetZ(Z);
}

inline void Vertex::SetTPos(double T)
{
    m_time = T;
}

inline double Vertex::GetXPos()
{
    return m_positionXYZ.X();
}

inline double Vertex::GetYPos()
{
    return m_positionXYZ.Y();
}

inline double Vertex::GetZPos()
{
    return m_positionXYZ.Z();
}

inline double Vertex::GetTPos()
{
    return m_time;
}

inline TVector3 Vertex::GetXYZPos()
{
    return m_positionXYZ;
}

inline void Vertex::SetXYZPos(double X, double Y, double Z)
{
    m_positionXYZ.SetXYZ(X,Y,Z);
}

inline void Vertex::SetXYZTPos(double X, double Y, double Z, double T)
{
    m_positionXYZ.SetXYZ(X,Y,Z);
    m_time = T;
}
