// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "space_cylinder.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


SpaceCylinder::SpaceCylinder(SpaceProp const* p)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinder is only valid in 3D: use rectangle instead");
    half_ = 0;
    radius_ = 0;
}


void SpaceCylinder::resize(Glossary& opt)
{
    real len = half_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    if ( opt.set(len, "length") )
        len *= 0.5;

    if ( len < 0 )
        throw InvalidParameter("cylinder:length must be >= 0");

    if ( rad < 0 )
        throw InvalidParameter("cylinder:radius must be >= 0");
    
    half_ = len;
    radius_ = rad;
}


void SpaceCylinder::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_,-radius_,-radius_);
    sup.set( half_, radius_, radius_);
}


bool SpaceCylinder::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.normYZSqr();
    return ( abs_real(W.XX) <= half_  &&  RT <= square(radius_) );
#elif ( DIM > 1 )
    return ( abs_real(W.XX) <= half_  &&  abs_real(W.YY) <= radius_ );
#else
    return false;
#endif
}


bool SpaceCylinder::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    const real RT = W.normYZSqr();
    return ( abs_real(W.XX) + rad <= half_  &&  RT <= square(radius_-rad) );
#elif ( DIM > 1 )
    return ( abs_real(W.XX) + rad <= half_  &&  abs_real(W.YY) <= radius_-rad );
#else
    return false;
#endif
}


Vector SpaceCylinder::place() const
{
#if ( DIM >= 3 )
    const Vector2 V = Vector2::randB(radius_);
    return Vector(half_*RNG.sreal(), V.XX, V.YY);
#elif ( DIM > 1 )
    return Vector(half_*RNG.sreal(), radius_*RNG.sreal());
#else
    return Vector(half_*RNG.sreal());
#endif
}

//------------------------------------------------------------------------------
Vector SpaceCylinder::project(Vector const& W) const
{
    Vector P(W);
#if ( DIM >= 3 )
    bool in = true;
    if ( abs_real(W.XX) > half_ )
    {
        P.XX = std::copysign(half_, W.XX);
        in = false;
    }
    
    real n = W.normYZ();
    
    if ( n > radius_ )
    {
        n = radius_ / n;
        P.YY = n * W.YY;
        P.ZZ = n * W.ZZ;
    }
    else if ( in )
    {
        if ( half_ - abs_real(W.XX) < radius_ - n )
        {
            P.XX = std::copysign(half_, W.XX);
        }
        else
        {
            n = radius_ / n;
            P.YY = n * W.YY;
            P.ZZ = n * W.ZZ;
        }
    }
#endif
    return P;
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical part and the caps.
 */
void SpaceCylinder::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca,
                                   real stiff, const real len, const real rad)
{
    bool cap = ( abs_real(pos.XX) > len );
    bool cyl = false;
    real X = std::copysign(len, pos.XX);
    
#if ( DIM > 2 )
    
    real dis = pos.YY*pos.YY + pos.ZZ*pos.ZZ;
    
    if ( rad*rad < dis )
    {
        // outside cylinder in YZ plane
        cyl = true;
    }
    else if ( ! cap )
    {
        // inside cylinder in YZ plane and also inside in X:
        //if ( abs_real( pos.XX - p ) > rad - std::sqrt(dis) )
        if ( dis > square( rad - abs_real(pos.XX-X) ) )
            cyl = true;
        else
            cap = true;
    }
    
#endif

    if ( cap )
        meca.addPlaneClampX(mp, X, stiff);
  
    if ( cyl )
        meca.addCylinderClampX(mp, rad, stiff);
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinder::setConfinement(Vector const& pos, Mecapoint const& mp,
                                   real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);
    real L = max_real(0, half_ - rad);
    
    setConfinement(pos, mp, meca, stiff, L, R);
}

//------------------------------------------------------------------------------

void SpaceCylinder::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LR");
    out.writeUInt16(2);
    out.writeFloat(half_);
    out.writeFloat(radius_);
}


void SpaceCylinder::setLengths(const real len[8])
{
    half_ = len[0];
    radius_ = len[1];
}

void SpaceCylinder::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LR");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_view.h"

void SpaceCylinder::draw3D() const
{
    const float L(half_);
    const float R(radius_);

    gym::stretchAlignZX(-L, L, R);
    gle::tube1();
    gle::discBottom1();
    gle::discTop1();
}

#else

void SpaceCylinder::draw3D() const {}

#endif

