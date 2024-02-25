// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "space_ring.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceRing::SpaceRing(SpaceProp const* p)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("ring is only valid in 3D: use rectangle instead");
    half_ = 0;
    radius_ = 0;
}


void SpaceRing::resize(Glossary& opt)
{
    real len = half_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    if ( opt.set(len, "length") )
        len *= 0.5;

    if ( len < 0 )
        throw InvalidParameter("ring:length must be > 0");
    if ( rad < 0 )
        throw InvalidParameter("ring:radius must be >= 0");

    half_ = len;
    radius_ = rad;
}


void SpaceRing::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_,-radius_,-radius_);
    sup.set( half_, radius_, radius_);
}


Vector SpaceRing::place() const
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
bool SpaceRing::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.normYZSqr();
    return ( abs_real(W.XX) <= half_  &&  RT <= square(radius_) );
#else
    return false;
#endif
}

bool SpaceRing::allInside(Vector const& W, const real rad ) const
{
    assert_true( rad >= 0 );

#if ( DIM > 2 )
    const real RT = W.normYZSqr();
    return ( abs_real(W.XX) + rad <= half_  &&  RT <= square(radius_-rad) );
#else
    return false;
#endif
}

//------------------------------------------------------------------------------
/**
 Project always on the surface of the cylinder
 */
Vector SpaceRing::project(Vector const& W) const
{
    Vector P;
    P.XX = min_real(max_real(W.XX, -half_), half_);
    
#if ( DIM > 2 )
    real n = W.normYZ();
    
    if ( n > 0 )
    {
        n = radius_ / n;
        P.YY = n * W.YY;
        P.ZZ = n * W.ZZ;
    }
    else
    {
        P.YY = radius_;
        P.ZZ = 0;
    }
#endif
    return P;
}

//------------------------------------------------------------------------------

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setConfinement(Vector const& pos, Mecapoint const& mp,
                               Meca& meca, real stiff) const
{
    if ( abs_real(pos.XX) > half_ )
        meca.addPlaneClampX(mp, std::copysign(half_, pos.XX), stiff);
    
    meca.addCylinderClampX(mp, radius_, stiff);
}

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setConfinement(Vector const& pos, Mecapoint const& mp,
                               real rad, Meca& meca, real stiff) const
{
    if ( radius_ > rad )
    {
        real len = max_real(half_-rad, 0);
        if ( abs_real(pos.XX) > len )
            meca.addPlaneClampX(mp, std::copysign(len, pos.XX), stiff);
        
        meca.addCylinderClampX(mp, radius_-rad, stiff);
    }
    else
    {
        meca.addLineClampX(mp, stiff);
        std::cerr << "object is too big to fit in SpaceRing\n";
    }
}

//------------------------------------------------------------------------------

void SpaceRing::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LR");
    out.writeUInt16(2);
    out.writeFloat(half_);
    out.writeFloat(radius_);
}


void SpaceRing::setLengths(const real len[8])
{
    half_ = len[0];
    radius_ = len[1];
}


void SpaceRing::read(Inputter& in, Simul&, ObjectTag)
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

void SpaceRing::draw3D() const
{
    const float L(half_);
    const float R(radius_);

    gym::stretchAlignZX(-L, L, R);
    gle::tube1();
}

#else

void SpaceRing::draw3D() const {}

#endif

