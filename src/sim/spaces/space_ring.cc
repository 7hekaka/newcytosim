// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_ring.h"
#include "mecapoint.h"
#include "exceptions.h"
#include "meca.h"


SpaceRing::SpaceRing(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1]), radiusSqr(mLengthSqr[1])
{
    if ( DIM < 3 )
        throw InvalidParameter("ring is only valid in 3D: use rectangle instead");
}


void SpaceRing::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length,-radius,-radius);
    sup.set( length, radius, radius);
}


real  SpaceRing::volume() const
{
    return 2 * M_PI * length * radius * radius;
}


Vector SpaceRing::randomPlace() const
{
#if ( DIM >= 3 )
    Vector2 sec = Vector2::randB(radius);
    return Vector(length*RNG.sreal(), sec.XX, sec.YY);
#elif ( DIM > 1 )
    return Vector(length*RNG.sreal(), radius*RNG.sreal());
#else
    return Vector(length*RNG.sreal());
#endif
}


//------------------------------------------------------------------------------
bool  SpaceRing::inside(Vector const& w) const
{
    if ( fabs(w.XX) > length )
        return false;

#if ( DIM > 2 )
    return ( w.YY*w.YY + w.ZZ*w.ZZ <= radiusSqr );
#else
    return false;
#endif
}

bool  SpaceRing::allInside(Vector const& w, const real rad ) const
{
    assert_true( rad >= 0 );
    
    if ( fabs(w.XX) + rad > length )
        return false;

#if ( DIM > 2 )
    return ( w.YY*w.YY + w.ZZ*w.ZZ <= square(radius-rad) );
#else
    return false;
#endif
}

//------------------------------------------------------------------------------
/**
 Project always on the surface of the cylinder
 */
Vector SpaceRing::project(Vector const& w) const
{
    Vector p;
    if ( w.XX >  length )
        p.XX =  length;
    else if ( w.XX < -length )
        p.XX = -length;
    else
        p.XX = w.XX;
    
#if ( DIM > 2 )
    real n = sqrt( w.YY*w.YY+ w.ZZ*w.ZZ );
    
    if ( n > 0 )
    {
        n = radius / n;
        p.YY = n * w.YY;
        p.ZZ = n * w.ZZ;
    }
    else
    {
        p.YY = radius;
        p.ZZ = 0;
    }
#endif
    return p;
}

//------------------------------------------------------------------------------

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff, const real len, const real rad)
{
    const Matrix::index_t inx = DIM * pe.matIndex();

    if ( pos.XX > len )
    {
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    += stiff * len;
    }
    else if ( pos.XX < -len )
    {
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    -= stiff * len;
    }
    
    meca.addCylinderClampX(pe, rad, stiff);
}


/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

/**
 This applies a force directed to the surface of the cylinder
 */
void SpaceRing::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceRing::draw() const
{
#if ( DIM > 2 )

    const size_t fin = 512;
    GLfloat c[fin+1], s[fin+1];
    gle::circle(fin, c, s, GLfloat(radius));

    GLfloat L = GLfloat(length);
    
    glBegin(GL_TRIANGLE_STRIP);
    for ( size_t n = 0; n <= fin; ++n )
    {
        glNormal3f( 0, c[n], s[n]);
        glVertex3f(+L, c[n], s[n]);
        glVertex3f(-L, c[n], s[n]);
    }
    glEnd();
    
#endif
    return true;
}

#else

bool SpaceRing::draw() const
{
    return false;
}

#endif

