// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "space_capsule.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


SpaceCapsule::SpaceCapsule(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("capsule is only defined for DIM = 2 and 3");
}


void SpaceCapsule::resize(Glossary& opt)
{
    real len = half_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    // total length is specified:
    if ( opt.set(len, "length") )
        len = ( len - 2 * rad ) * 0.5;

    if ( len < 0 )
        throw InvalidParameter("capsule:length must be >= 2 * radius");
    if ( rad < 0 )
        throw InvalidParameter("capsule:radius must be >= 0");
    
    half_ = len;
    radius_ = rad;
}


void SpaceCapsule::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_-half_,-radius_,-radius_);
    sup.set( radius_+half_, radius_, radius_);
}


real SpaceCapsule::volume() const
{
#if ( DIM >= 3 )
    return (2.0*M_PI) * half_ * square(radius_) + (4.0*M_PI/3.0) * cube(radius_);
#else
    return 4 * half_ * radius_ + M_PI * square(radius_);
#endif
}


real SpaceCapsule::surface() const
{
#if ( DIM >= 3 )
    return (4.0*M_PI) * half_ * radius_ + (4.0*M_PI) * square(radius_);
#else
    return 4 * half_ + (2.0*M_PI) * radius_;
#endif
}


bool SpaceCapsule::inside(Vector const& W) const
{
    real n = W.normYZSqr() + square(max_real(0, abs_real(W.XX)-half_));
    return ( n <= square(radius_) );
}


bool SpaceCapsule::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
    real n = W.normYZSqr() + square(max_real(0, abs_real(W.XX)-half_));
    
    return ( n <= square(radius_-rad) );
}

//------------------------------------------------------------------------------
Vector SpaceCapsule::project(Vector const& pos) const
{
    // calculate the projection on the axis, within boundaries:
    real X = min_real(half_, max_real(pos.XX, -half_));
    
#if ( DIM >= 3 )
    real n = square(pos.XX-X) + pos.normYZSqr();
    if ( n > 0 )
    {
        n = radius_ / std::sqrt(n);
        return Vector( X + n * ( pos.XX - X ), n * pos.YY, n * pos.ZZ);
    }
    Vector2 YZ = Vector2::randU(radius_);
    return Vector(0, YZ.XX, YZ.YY);
#elif ( DIM >= 2 )
    real n = square(pos.XX-X) + square(pos.YY);
    if ( n > 0 )
    {
        n = radius_ / std::sqrt(n);
        return Vector(X + n * ( pos.XX - X ), n * pos.YY, 0);
    }
#endif
    return Vector(0, radius_*RNG.sflip(), 0);
}


Vector SpaceCapsule::randomPlace() const
{
#if ( DIM >= 3 )
    // volume elements divided by 4 * M_PI * radius_
    const real V0 = half_;
    const real V1 = radius_;  // spherical caps
    const real P = RNG.preal() * ( V0 + V1 );
    if ( P < V0 )
    {
        Vector2 YZ = Vector2::randB(radius_);
        return Vector(RNG.sreal()*half_, YZ.XX, YZ.YY);
    }
#elif ( DIM >= 2 )
    // surface elements divided by radius_
    const real V0 = 4 * half_;
    const real V1 = M_PI * radius_;  // spherical caps
    const real P = RNG.preal() * ( V0 + V1 );
    if ( P < V0 )
        return Vector(RNG.sreal()*half_, RNG.sreal()*radius_, 0);
#endif
    // a position on the caps:
    Vector vec = Vector::randB(radius_);
    vec.XX += std::copysign(half_, vec.XX);
    return vec;
}


Vector SpaceCapsule::normalToEdge(Vector const& pos) const
{
    real X = min_real(half_, max_real(pos.XX, -half_));
#if ( DIM >= 3 )
    real n = square(pos.XX-X) + pos.normYZSqr();
    if ( n > 0 )
    {
        n = 1.0 / std::sqrt(n);
        return Vector(n * ( pos.XX - X ), n * pos.YY, n * pos.ZZ);
    }
    Vector2 YZ = Vector2::randU();
    return Vector(0, YZ.XX, YZ.YY);
#elif ( DIM >= 2 )
    real n = square(pos.XX-X) + square(pos.YY);
    if ( n > 0 )
    {
        n = 1.0 / std::sqrt(n);
        return Vector(n * ( pos.XX - X ), n * pos.YY, 0);
    }
#endif
    return Vector(0, RNG.sflip(), 0);
}


Vector SpaceCapsule::randomPlaceOnEdge(real) const
{
#if ( DIM >= 3 )
    // surface elements divided by 4 * M_PI * radius_
    const real S0 = half_;
    const real S1 = radius_;  // spherical caps
    const real P = RNG.preal() * ( S0 + S1 );
    if ( P < S0 )
    {
        Vector2 YZ = Vector2::randU(radius_);
        return Vector(RNG.sreal()*half_, YZ.XX, YZ.YY);
    }
#else
    // length elements divided by 2
    const real S0 = half_;
    const real S1 = M_PI * radius_;  // spherical caps
    const real P = RNG.preal() * ( S0 + S1 );
    if ( P < S0 )
        return Vector(RNG.sreal()*half_, RNG.sflip()*radius_, 0);
#endif
    // a position on the caps:
    Vector vec = Vector::randU(radius_);
    vec.XX += std::copysign(half_, vec.XX);
    return vec;
}


//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff, const real len, const real rad)
{
    if ( abs_real(pos.XX) > len )
    {
        Vector cen(std::copysign(len, pos.XX),0,0);
        meca.addSphereClamp(pos-cen, pe, cen, rad, stiff);
    }
    else
        meca.addCylinderClampX(pe, rad, stiff);
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, half_, radius_);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    if ( rad < radius_ )
        setInteraction(pos, pe, meca, stiff, half_, radius_-rad);
    else
        setInteraction(pos, pe, meca, stiff, half_, 0);
}

//------------------------------------------------------------------------------

void SpaceCapsule::write(Outputter& out) const
{
    writeShape(out, "capsule");
    out.writeUInt16(2);
    out.writeFloat(half_);
    out.writeFloat(radius_);
}


void SpaceCapsule::setLengths(const real len[])
{
    half_ = len[0];
    radius_ = len[1];
}

void SpaceCapsule::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "capsule");
    setLengths(len);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
#include "gle_flute.h"

void SpaceCapsule::draw2D() const
{
    const float L(half_);
    const real R(radius_);
    constexpr size_t fin = 16 * gle::finesse;
    float* arc = gle::mapFloatBuffer(4*fin+4);
    float* cra = arc + 2*fin + 2;
    gle::compute_arc(fin, arc, R, -M_PI_2, M_PI,  L, 0);
    gle::compute_arc(fin, cra, R,  M_PI_2, M_PI, -L, 0);
    gle::unmapFloatBuffer(2);
    glDrawArrays(GL_LINE_LOOP, 0, 2*fin+2);
}


void SpaceCapsule::draw3D() const
{
    const real L(half_);
    const float R(radius_);
    
    const GLenum glp = GL_CLIP_PLANE5;
    GLdouble plane[] = { 1, 0, 0, 0 };
    glEnable(glp);

    //right side:
    glPushMatrix();
    glClipPlane(glp, plane);
    gle::transAlignZ(Vector(L,0,0), R, Vector(-1,0,0));
    gle::halfTube1();
    gle::hemisphere4();
    gle::arrowStrip(0.5, 2);
    glPopMatrix();

    //left side:
    glPushMatrix();
    plane[0] = -1;
    glClipPlane(glp, plane);
    gle::transAlignZ(Vector(-L,0,0), R, Vector(1,0,0));
    gle::halfTube1();
    gle::hemisphere4();
    gle::arrowStrip(0.5, 2);
    glPopMatrix();
    
    glDisable(glp);
}

#else

void SpaceCapsule::draw2D() const {}
void SpaceCapsule::draw3D() const {}

#endif
