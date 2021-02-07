// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
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
    length_ = 0;
    radius_ = 0;
}


void SpaceCylinder::resize(Glossary& opt)
{
    real len = length_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    if ( opt.set(len, "length") )
        len *= 0.5;

    if ( len < 0 )
        throw InvalidParameter("cylinder:length must be >= 0");

    if ( rad < 0 )
        throw InvalidParameter("cylinder:radius must be >= 0");
    
    length_ = len;
    radius_ = rad;
}


void SpaceCylinder::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_,-radius_,-radius_);
    sup.set( length_, radius_, radius_);
}


real SpaceCylinder::volume() const
{
    return 2 * M_PI * length_ * radius_ * radius_;
}


bool SpaceCylinder::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.YY * W.YY + W.ZZ * W.ZZ;
    return ( abs_real(W.XX) < length_  &&  RT <= radius_ * radius_ );
#elif ( DIM > 1 )
    return ( abs_real(W.XX) < length_  &&  abs_real(W.YY) <= radius_ );
#else
    return false;
#endif
}


bool SpaceCylinder::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    const real RT = W.YY * W.YY + W.ZZ * W.ZZ;
    return ( abs_real(W.XX) + rad < length_  &&  RT <= square(radius_-rad) );
#elif ( DIM > 1 )
    return ( abs_real(W.XX) + rad < length_  &&  abs_real(W.YY) <= radius_-rad );
#else
    return false;
#endif
}


Vector SpaceCylinder::randomPlace() const
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
Vector SpaceCylinder::project(Vector const& W) const
{
    Vector P(W);
#if ( DIM >= 3 )
    bool in = true;
    if ( abs_real(W.XX) > length_ )
    {
        P.XX = std::copysign(length_, W.XX);
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
        if ( length_ - abs_real(W.XX) < radius_ - n )
        {
            P.XX = std::copysign(length_, W.XX);
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
void SpaceCylinder::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca,
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
        meca.addPlaneClampX(pe, X, stiff);
  
    if ( cyl )
        meca.addCylinderClampX(pe, rad, stiff);
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinder::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length_, radius_);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinder::setInteraction(Vector const& pos, Mecapoint const& pe,
                                   real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);
    real L = max_real(0, length_ - rad);
    
    setInteraction(pos, pe, meca, stiff, L, R);
}

//------------------------------------------------------------------------------

void SpaceCylinder::write(Outputter& out) const
{
    writeShape(out, "cylinder");
    out.writeUInt16(2);
    out.writeFloat(length_);
    out.writeFloat(radius_);
}


void SpaceCylinder::setLengths(const real len[])
{
    length_ = len[0];
    radius_ = len[1];
}

void SpaceCylinder::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "cylinder");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

void SpaceCylinder::draw3D() const
{
    const size_t fin = 512;
    GLfloat cir[2*fin+2];
    gle::compute_circle(fin, cir, 1);

    GLfloat L(length_);
    GLfloat R(radius_);
    
    glBegin(GL_TRIANGLE_STRIP);
    for ( size_t i = 0; i <= fin; ++i )
    {
        GLfloat c = cir[2*i], s = cir[1+2*i];
        glNormal3f(0, c, s);
        glVertex3f(+L, R*c, R*s);
        glVertex3f(-L, R*c, R*s);
    }
    glEnd();
    
    // draw the cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( +1, 0, 0 );
    glVertex3f( +L, 0, 0 );
    for ( size_t i = 0; i <= fin; ++i )
        glVertex3f(+L, R*cir[2*i], R*cir[1+2*i]);
    glEnd();
    
    // draw the cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( -1, 0, 0 );
    glVertex3f( -L, 0, 0 );
    for ( size_t i = 0; i <= fin; ++i )
        glVertex3f(-L, -R*cir[2*i], R*cir[1+2*i]);
    glEnd();
}

#else

void SpaceCylinder::draw3D() const {}

#endif

