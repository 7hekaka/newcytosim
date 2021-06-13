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
    modulo_.reset();
    modulo_.enable(0, 2*half_);
}


void SpaceCylinderP::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_,-radius_,-radius_);
    sup.set( half_, radius_, radius_);
}


real SpaceCylinderP::volume() const
{
    return 2 * M_PI * half_ * square(radius_);
}


bool SpaceCylinderP::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.YY * W.YY + W.ZZ * W.ZZ;
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
    const real RT = W.YY * W.YY + W.ZZ * W.ZZ;
    return ( RT <= square(radius_-rad) );
#elif ( DIM > 1 )
    return ( abs_real(W.YY) <= radius_-rad );
#else
    return false;
#endif
}


Vector SpaceCylinderP::randomPlace() const
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


Vector SpaceCylinderP::randomPlaceOnEdge(real) const
{
#if ( DIM >= 3 )
    const Vector2 YZ = Vector2::randU(radius_);
    return Vector(half_*RNG.sreal(), YZ.XX, YZ.YY);
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
        const Vector2 V = Vector2::randU();
        P.YY = radius_ * V.XX;
        P.ZZ = radius_ * V.YY;
    }
#endif
    return P;
}


void SpaceCylinderP::bounce(Vector& pos) const
{
    if ( !inside(pos) )
        bounceOnEdges(pos);

    pos.XX = fold_real(pos.XX, modulo_.period_[0]);
}

//------------------------------------------------------------------------------

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setConfinement(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    meca.addCylinderClampX(pe, radius_, stiff);
}

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setConfinement(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);

    meca.addCylinderClampX(pe, R, stiff);
}

//------------------------------------------------------------------------------

void SpaceCylinderP::write(Outputter& out) const
{
    writeShape(out, "cylinderP");
    out.writeUInt16(2);
    out.writeFloat(half_);
    out.writeFloat(radius_);
}


void SpaceCylinderP::setLengths(const real len[])
{
    half_ = len[0];
    radius_ = len[1];
    update();
}


void SpaceCylinderP::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "cylinderP");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"

void SpaceCylinderP::draw3D() const
{
    GLfloat L(half_);
    GLfloat R(radius_);

    glPushMatrix();
    gle::stretchAlignZX(-L, L, R);
    gle::tube1();
    glPopMatrix();

    if ( 1 )
    {
        // mark the edge of the periodicity with dotted lines
        glLineStipple(1, 0x000F);
        glEnable(GL_LINE_STIPPLE);
        glPushMatrix();
        gle::transAlignZX(-L, R, -1);
        gle::circle();
        glPopMatrix();
        glPushMatrix();
        gle::transAlignZX(L, R, -1);
        gle::circle();
        glPopMatrix();
        glDisable(GL_LINE_STIPPLE);
    }
}

#else

void SpaceCylinderP::draw3D() const {}

#endif

