// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinderP.h"
#include "mecapoint.h"
#include "exceptions.h"
#include "meca.h"


SpaceCylinderP::SpaceCylinderP(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1])
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinderP is only valid in 3D: use strip instead");
}

void SpaceCylinderP::setModulo(Modulo& mod) const
{
    mod.enable(0, length);
}

void SpaceCylinderP::resize()
{
    Space::checkLengths(2, false);
    
    if ( length <= 0 )
        throw InvalidParameter("length of cylinderP must be > 0");
}


void SpaceCylinderP::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length,-radius,-radius);
    sup.set( length, radius, radius);
}



real  SpaceCylinderP::volume() const
{
    return 2 * M_PI * length * radius * radius;
}


bool  SpaceCylinderP::inside(Vector const& w) const
{
#if ( DIM > 2 )
    return ( w.YY*w.YY + w.ZZ*w.ZZ <= radius * radius );
#elif ( DIM > 1 )
    return ( fabs(w.YY) <= radius );
#else
    return false;
#endif
}


bool SpaceCylinderP::allInside(Vector const& w, const real rad ) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    return ( w.YY*w.YY + w.ZZ*w.ZZ <= square(radius-rad) );
#elif ( DIM > 1 )
    return ( fabs(w.YY) <= radius-rad );
#else
    return false;
#endif
}


Vector SpaceCylinderP::randomPlace() const
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
Vector SpaceCylinderP::project(Vector const& w) const
{
    Vector p;
    p.XX = w.XX;
    
#if ( DIM > 2 )
    real n = w.normYZ();
    if ( n > REAL_EPSILON )
    {
        p.YY = w.YY * ( radius / n );
        p.ZZ = w.ZZ * ( radius / n );
    }
    else
    {
        Vector2 yz = Vector2::randU();
        p.YY = radius * yz.XX;
        p.ZZ = radius * yz.YY;
    }
#endif
    return p;
}

//------------------------------------------------------------------------------

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    meca.addCylinderClampX(pe, radius, stiff);
}

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    real eRadius = radius - rad;
    if ( eRadius < 0 ) eRadius = 0;
    
    meca.addCylinderClampX(pe, eRadius, stiff);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCylinderP::draw() const
{
#if ( DIM > 2 )

    const size_t fin = 512;
    GLfloat c[fin+1], s[fin+1];
    gle::circle(fin, c, s, 1);

    GLfloat L = (GLfloat)length;
    GLfloat R = (GLfloat)radius;

    glBegin(GL_TRIANGLE_STRIP);
    for ( size_t n = 0; n <= fin; ++n )
    {
        glNormal3f( 0, c[n], s[n] );
        glVertex3f( +L, R*c[n], R*s[n] );
        glVertex3f( -L, R*c[n], R*s[n] );
    }
    glEnd();
    
    if ( 1 )
    {
        //draw dotted-rings to indicate periodicity
        glLineStipple(1, 0x000F);
        glEnable(GL_LINE_STIPPLE);
        glPushMatrix();
        glTranslatef(L, 0, 0);
        glScalef(R, R, R);
        glRotated(90, 0, 1, 0);
        gle::gleCircle();
        glTranslatef(0, 0, -2*L/R);
        glRotated(180, 0, 1, 0);
        gle::gleCircle();
        glPopMatrix();
        glDisable(GL_LINE_STIPPLE);
    }

#endif
    return true;
}

#else

bool SpaceCylinderP::draw() const
{
    return false;
}

#endif

