// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinder.h"
#include "mecapoint.h"
#include "exceptions.h"
#include "meca.h"


SpaceCylinder::SpaceCylinder(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1])
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinder is only valid in 3D: use rectangle instead");
}


void SpaceCylinder::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length,-radius,-radius);
    sup.set( length, radius, radius);
}


real SpaceCylinder::volume() const
{
    return 2 * M_PI * length * radius * radius;
}


bool SpaceCylinder::inside(Vector const& w) const
{
    if ( fabs(w.XX) > length )
        return false;
    
#if ( DIM > 2 )
    return ( w.YY*w.YY + w.ZZ*w.ZZ <= radius * radius );
#elif ( DIM > 1 )
    return ( fabs(w.YY) <= radius );
#else
    return false;
#endif
}


bool SpaceCylinder::allInside(Vector const& w, const real rad) const
{
    assert_true( rad >= 0 );
    
    if ( fabs(w.XX) + rad > length )
        return false;
    
#if ( DIM > 2 )
    return ( w.YY*w.YY + w.ZZ*w.ZZ <= square(radius-rad) );
#elif ( DIM > 1 )
    return ( fabs(w.YY) <= radius-rad );
#else
    return false;
#endif
}


Vector SpaceCylinder::randomPlace() const
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
Vector SpaceCylinder::project(Vector const& w) const
{
    Vector p = w;
#if ( DIM >= 3 )
    bool inX = 1;
    
    if ( fabs(w.XX) > length )
    {
        p.XX = std::copysign(length, w.XX);
        inX = 0;
    }
    
    real n = w.normYZ();
    
    if ( n > radius )
    {
        n = radius / n;
        p.YY = n * w.YY;
        p.ZZ = n * w.ZZ;
    }
    else
    {
        if ( inX )
        {
            if ( length - fabs(w.XX) < radius - n )
            {
                p.XX = std::copysign(length, w.XX);
            }
            else
            {
                n = radius / n;
                p.YY = n * w.YY;
                p.ZZ = n * w.ZZ;
            }
        }
    }
#endif
    return p;
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical part and the caps.
 */
void SpaceCylinder::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca,
                                   real stiff, const real len, const real rad)
{
    bool cap = ( fabs(pos.XX) > len );
    bool cyl = false;
    real p = std::copysign(len, pos.XX);
    
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
        if ( fabs( pos.XX - p ) > rad - sqrt(dis) )
            cyl = true;
        else
            cap = true;
    }
    
#endif

    if ( cap )
    {
        const Matrix::index_t inx = DIM * pe.matIndex();
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    += stiff * p;
    }
  
    if ( cyl )
        meca.addCylinderClampX(pe, rad, stiff);
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinder::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinder::setInteraction(Vector const& pos, Mecapoint const& pe,
                                   real rad, Meca & meca, real stiff) const
{
    real eRadius = radius - rad;
    if ( eRadius < 0 ) eRadius = 0;
    real eLength = length - rad;
    if ( eLength < 0 ) eLength = 0;
    
    setInteraction(pos, pe, meca, stiff, eLength, eRadius);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCylinder::draw() const
{
#if ( DIM > 2 )

    const size_t fin = 512;

    GLfloat L = (GLfloat)length;
    GLfloat R = (GLfloat)radius;
    
    GLfloat c[fin+1], s[fin+1];
    gle::circle(fin, c, s, 1);
    
    glBegin(GL_TRIANGLE_STRIP);
    for ( size_t sc = 0; sc <= fin; ++sc )
    {
        GLfloat ca = c[sc], sa = s[sc];
        glNormal3f( 0, ca, sa );
        glVertex3f( +L, R*ca, R*sa );
        glVertex3f( -L, R*ca, R*sa );
    }
    glEnd();
    
    // draw the cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( +1, 0, 0 );
    glVertex3f( +L, 0, 0 );
    for ( size_t sc = 0; sc <= fin; ++sc )
        glVertex3f( +L, R*c[sc], R*s[sc] );
    glEnd();
    
    // draw the cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f( -1, 0, 0 );
    glVertex3f( -L, 0, 0 );
    for ( size_t sc = 0; sc <= fin; ++sc )
        glVertex3f( -L,-R*c[sc], R*s[sc] );
    glEnd();
    
#endif
    return true;
}

#else

bool SpaceCylinder::draw() const
{
    return false;
}

#endif

