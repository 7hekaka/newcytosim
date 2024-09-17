// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "space_cylinderP.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


SpaceCylinderP::SpaceCylinderP(SpaceProp const* p)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinderP is only valid in 3D: use strip instead");
    half_ = 0;
    radius_ = 0;
}

void SpaceCylinderP::resize(Glossary& opt)
{
    real len = half_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");

    if ( opt.set(len, "length") )
        len *= 0.5;

    if ( rad < 0 )
        throw InvalidParameter("cylinderP:radius must be >= 0");

    if ( len <= 0 )
        throw InvalidParameter("cylinderP:length must be > 0");
    
    half_ = len;
    radius_ = rad;

    update();
}


void SpaceCylinderP::update()
{
#if ENABLE_PERIODIC_BOUNDARIES
    modulo_.reset();
    modulo_.enablePeriodic(0, 2*half_);
#endif
}


void SpaceCylinderP::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_,-radius_,-radius_);
    sup.set( half_, radius_, radius_);
}


bool SpaceCylinderP::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.normYZSqr();
    return ( RT <= square(radius_) );
#elif ( DIM > 1 )
    return ( abs_real(W.YY) <= radius_ );
#else
    return false;
#endif
}


bool SpaceCylinderP::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    const real RT = W.normYZSqr();
    return ( RT <= square(radius_-rad) );
#elif ( DIM > 1 )
    return ( abs_real(W.YY) <= radius_-rad );
#else
    return false;
#endif
}


Vector SpaceCylinderP::place() const
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


Vector SpaceCylinderP::normalToEdge(Vector const& pos) const
{
#if ( DIM >= 3 )
    real n = 1.0 / pos.normYZ();
    return Vector(0, n * pos.YY, n * pos.ZZ);
#elif ( DIM >= 2 )
    return Vector(0, sign_real(pos.YY), 0);
#endif
    return Vector(0, 0, 0);  // intentionally invalid!
}


Vector SpaceCylinderP::placeOnEdge(real) const
{
#if ( DIM >= 3 )
    real C, S;
    RNG.urand2(C, S, radius_);
    return Vector(half_*RNG.sreal(), C, S);
#endif
    return Vector(half_*RNG.sreal(), radius_*RNG.sflip(), 0);
}


//------------------------------------------------------------------------------
Vector SpaceCylinderP::project(Vector const& W) const
{
    Vector P(W);
    
#if ( DIM > 2 )
    real n = W.normYZ();
    if ( n > REAL_EPSILON )
    {
        P.YY = W.YY * ( radius_ / n );
        P.ZZ = W.ZZ * ( radius_ / n );
    }
    else
    {
        real C, S;
        RNG.urand2(C, S, radius_);
        P.YY = C;
        P.ZZ = S;
    }
#endif
    return P;
}


Vector SpaceCylinderP::bounce(Vector const& pos) const
{
    /* This is not correct and may lead to artifacts for small radius */
    if ( !inside(pos) )
        return bounceOnEdges(pos);
    
    Vector P = pos;
#if ENABLE_PERIODIC_BOUNDARIES
    P.XX = modulo_.fold_(pos.XX, 0);
#endif
    return P;
}

//------------------------------------------------------------------------------

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setConfinement(Vector const& pos, Mecapoint const& mp,
                                    Meca& meca, real stiff) const
{
    meca.addCylinderClampX(mp, radius_, stiff);
}

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setConfinement(Vector const& pos, Mecapoint const& mp,
                                    real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);

    meca.addCylinderClampX(mp, R, stiff);
}

//------------------------------------------------------------------------------

void SpaceCylinderP::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LR");
    out.writeUInt16(2);
    out.writeFloat(half_);
    out.writeFloat(radius_);
}


void SpaceCylinderP::setLengths(const real len[8])
{
    half_ = len[0];
    radius_ = len[1];
    update();
}


void SpaceCylinderP::read(Inputter& in, Simul&, ObjectTag)
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

void SpaceCylinderP::draw3D() const
{
    const float L(half_);
    const float R(radius_);

    gym::stretchAlignZX(-L, L, R);
    gle::tube1();
    float WIDTH = 1;
    if ( 1 )
    {
        // mark the edge of the periodicity with dotted lines
        gym::transAlignZX(-L, R, -R);
        gle::dottedCircle(WIDTH);
        gym::transAlignZX(L, R, -R);
        gle::dottedCircle(WIDTH);
    }
}

#else

void SpaceCylinderP::draw3D() const {}

#endif

