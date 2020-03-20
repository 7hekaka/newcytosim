// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matrix34.h"
#include "vector2.h"
#include "random.h"

Vector3 Matrix34::rotationAxis() const
{
    return Vector3(val[6]-val[9], val[8]-val[2], val[1]-val[4]);
}

real Matrix34::rotationAngle() const
{
    return acos(0.5*(1-trace()));
}


void Matrix34::getEulerAngles(real& a, real& b, real& c) const
{
    real cb = sqrt(val[0] * val[0] + val[1] * val[1]);
    
    b = atan2(-val[2], cb);
    
    if ( cb != 0 ) {
        a = atan2(val[1], val[0]);
        c = atan2(val[6], val[10]);
    }
    else {
        a = 0;
        c = atan2(-val[5], val[6]);
    }
}




/// returns a rotation of angle PI around axis Z
Matrix34 Matrix34::rotation180()
{
    return Matrix34(-1, 0, 0, 0, -1, 0, 0, 0, 1);
}


Matrix34 Matrix34::rotationAroundX(const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    return Matrix34(1, 0, 0, 0, c, s, 0, -s, c);
}

Matrix34 Matrix34::rotationAroundY(const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    return Matrix34(c, 0, -s, 0, 1, 0, s, 0, c);
}

Matrix34 Matrix34::rotationAroundZ(const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    return Matrix34(c, s, 0, -s, c, 0, 0, 0, 1);
}


Matrix34 Matrix34::rotationAroundPrincipalAxis(index i, const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    
    i %= 3;
    index j = (i+1)%3;
    index k = (i+2)%3;
    
    Matrix34 res(0, 1);
    res(j,j) = c;
    res(k,j) = s;
    res(j,k) = -s;
    res(k,k) = c;
    return res;
}


Matrix34 Matrix34::rotationFromAngles(const real a[3])
{
    real ca = cos(a[0]), sa = sin(a[0]);
    real cb = cos(a[1]), sb = sin(a[1]);
    real cc = cos(a[2]), sc = sin(a[2]);
    
    Matrix34 res;

    res(0,0) =  ca*cb;
    res(1,0) =  sa*cb;
    res(2,0) = -sb;
    
    res(0,1) =  ca*sb*sc - sa*cc;
    res(1,1) =  sa*sb*sc + ca*cc;
    res(2,1) =  cb*sc;
    
    res(0,2) =  ca*sb*cc + sa*sc;
    res(1,2) =  sa*sb*cc - ca*sc;
    res(2,2) =  cb*cc;

    return res;
}


Matrix34 Matrix34::rotationAroundAxisEuler(const real a[3])
{
    real ca = cos(a[0]), sa = sin(a[0]), ca1 = 1 - ca;
    real cb = cos(a[1]), sb = sin(a[1]);
    real cc = cos(a[2]), sc = sin(a[2]);
    
    real sacc      = sa * cc,         sasc    = sa * sc;
    real saccsb    = sacc * sb,       sacccb  = sacc * cb;
    real ccccca1   = cc * cc * ca1,   ccscca1 = cc * sc * ca1;
    real cbccccca1 = cb * ccccca1;
    
    Matrix34 res;
    res(0,0) = cb * cbccccca1 + ca;
    res(0,1) = sb * cbccccca1 - sasc;
    res(0,2) = cb * ccscca1   + saccsb;
    
    res(1,0) = sb * cbccccca1 + sasc;
    res(1,1) = ca - cb * cbccccca1 + ccccca1;
    res(1,2) = sb * ccscca1   - sacccb;
    
    res(2,0) = cb * ccscca1 - saccsb;
    res(2,1) = sb * ccscca1 + sacccb;
    res(2,2) = 1 - ccccca1;

    return res;
}


Matrix34 Matrix34::randomRotation()
{
    //James Arvo, Fast random rotation matrices. in Graphics Gems 3.
    real u2 = M_PI * RNG.sreal();
    real u3 = RNG.preal();
    Vector3 V( cos(u2)*sqrt(u3), sin(u2)*sqrt(u3), sqrt(1-u3) );
    return householder(V) * rotationAroundZ(M_PI*RNG.sreal());
}


Matrix34 Matrix34::randomRotation(real angle)
{
    return rotationAroundAxis(Vector3::randU(), cos(angle), sin(angle));
}


Matrix34 Matrix34::rotationToVector(const Vector3& vec)
{
    Matrix34 res;
    Vector3 X, Y, Z = normalize(vec);
    Z.orthonormal(X, Y);
    res.setColumns(Z, X, Y);
    return res;
}


Matrix34 Matrix34::randomRotationToVector(const Vector3& vec)
{
    Matrix34 res;
    Vector3 X, Y, Z = normalize(vec);
    Z.orthonormal(X, Y);
#if ( 0 )
    real a = M_PI * RNG.sreal();
    real c = cos(a), s = sin(a);
    res.setColumns(Z, X*c+Y*s, Y*c-X*s);
#else
    const Vector2 V = Vector2::randU();
    res.setColumns(Z, X*V.XX+Y*V.YY, Y*V.XX-X*V.YY);
#endif
    return res;
}

