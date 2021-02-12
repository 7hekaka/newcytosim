// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
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
    length_ = 0;
    radius_ = 0;
}


void SpaceRing::resize(Glossary& opt)
{
    real len = length_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    if ( opt.set(len, "length") )
        len *= 0.5;

    if ( len < 0 )
        throw InvalidParameter("ring:length must be > 0");
    if ( rad < 0 )
        throw InvalidParameter("ring:radius must be >= 0");

    length_ = len;
    radius_ = rad;
}


void SpaceRing::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_,-radius_,-radius_);
    sup.set( length_, radius_, radius_);
}


real  SpaceRing::volume() const
{
    return 2 * M_PI * length_ * radius_ * radius_;
}


Vector SpaceRing::randomPlace() const
{
#if ( DIM >= 3 )
    const Vector2 V = Vector2::randB(radius_);
    return Vector(length_*RNG.sreal(), V.XX, V.YY);
#elif ( DIM > 1 )
    return Vector(length_*RNG.sreal(), radius_*RNG.sreal());
#else
    return Vector(length_*RNG.sreal());
#endif
}


//------------------------------------------------------------------------------
bool SpaceRing::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.YY * W.YY + W.ZZ * W.ZZ;
    return ( abs_real(W.XX) <= length_  &&  RT <= square(radius_) );
#else
    return false;
#endif
}

bool SpaceRing::allInside(Vector const& W, const real rad ) const
{
    assert_true( rad >= 0 );

#if ( DIM > 2 )
    const real RT = W.YY * W.YY + W.ZZ * W.ZZ;
    return ( abs_real(W.XX) + rad <= length_  &&  RT <= square(radius_-rad) );
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
    Vector P(W);
    if ( W.XX >  length_ )
        P.XX =  length_;
    else if ( W.XX < -length_ )
        P.XX = -length_;
    else
        P.XX = W.XX;
    
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
void SpaceRing::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff, const real len, const real rad)
{
    if ( abs_real(pos.XX) > len )
        meca.addPlaneClampX(pe, std::copysign(len, pos.XX), stiff);
    
    meca.addCylinderClampX(pe, rad, stiff);
}


/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length_, radius_);
}

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length_, radius_);
}

//------------------------------------------------------------------------------

void SpaceRing::write(Outputter& out) const
{
    writeShape(out, "ring");
    out.writeUInt16(2);
    out.writeFloat(length_);
    out.writeFloat(radius_);
}


void SpaceRing::setLengths(const real len[])
{
    length_ = len[0];
    radius_ = len[1];
}


void SpaceRing::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "ring");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

void SpaceRing::draw3D() const
{
    GLfloat L(length_);
    GLfloat R(radius_);

    const GLenum glp = GL_CLIP_PLANE5;
    GLdouble plane[] = { -1, 0, 0, L };
    glEnable(glp);
    glClipPlane(glp, plane);
    glPushMatrix();
    gle::transAlignZ(Vector(-L,0,0), R, Vector(1,0,0));
    gle::halfTube4();
    glPopMatrix();
    glDisable(glp);
}

#else

void SpaceRing::draw3D() const {}

#endif

