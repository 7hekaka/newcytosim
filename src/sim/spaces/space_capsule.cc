// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_capsule.h"
#include "mecapoint.h"
#include "exceptions.h"
#include "meca.h"



SpaceCapsule::SpaceCapsule(const SpaceProp* p)
: Space(p), length(mLength[0]), radius(mLength[1]), radiusSqr(mLengthSqr[1])
{
    if ( DIM == 1 )
        throw InvalidParameter("capsule is only defined for DIM = 2 or 3");
}


void SpaceCapsule::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius-length,-radius,-radius);
    sup.set( radius+length, radius, radius);
}


real SpaceCapsule::volume() const
{
#if ( DIM >= 3 )
    return ( length + (2/3.0) * radius ) * radiusSqr * ( 2 * M_PI );
#else
    return 4 * length * radius + M_PI * radiusSqr;
#endif
}


bool SpaceCapsule::inside(Vector const& w) const
{
#if ( DIM > 2 )
    real nrm = w.YY * w.YY + w.ZZ * w.ZZ;
#elif ( DIM > 1 )
    real nrm = w.YY * w.YY;
#else
    real nrm = 0;
#endif
    real x = fabs(w.XX);
    
    if ( x > length )
        nrm += square( x - length );
    
    return ( nrm <= radiusSqr );
}


bool SpaceCapsule::allInside(Vector const& w, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    real nrm = w.YY * w.YY + w.ZZ * w.ZZ;
#elif ( DIM > 1 )
    real nrm = w.YY * w.YY;
#else
    real nrm = 0;
#endif
    real x = fabs(w.XX);
    
    if ( x > length )
        nrm += square( x - length );
    
    return ( nrm <= square(radius-rad) );
}

//------------------------------------------------------------------------------
Vector SpaceCapsule::project(Vector const& w) const
{
    Vector p;
    real nrm = w.normYZ();
    
    //calculate the projection on the axis, within boundaries:
    if ( fabs(w.XX) > length )
    {
        real L = std::copysign(length, w.XX);
        nrm  += square( w.XX - L );
        //normalize from this point on the axis
        if ( nrm > 0 ) nrm = radius / sqrt( nrm );
        
        p.XX = length + nrm * ( w.XX - L );
    }
    else
    {
        //normalize from this point on the axis
        if ( nrm > 0 ) nrm = radius / sqrt( nrm );
        
        p.XX = w.XX;
    }
    
    if ( nrm > 0 )
    {
#if ( DIM > 1 )
        p.YY = nrm * w.YY;
#endif
#if ( DIM >= 3 )
        p.ZZ = nrm * w.ZZ;
#endif
    }
    else
    {
        //we project on a arbitrary point on the cylinder
#if ( DIM > 1 )
        p.YY = radius;
#endif
#if ( DIM >= 3 )
        p.ZZ = 0;
#endif
    }
    return p;
}



Vector SpaceCapsule::randomPlace() const
{
    size_t nb_trials = 1<<13;
    size_t ouf = 0;
    Vector res;

    do {
        
#if ( DIM == 1 )
        res.set((length+radius)*RNG.sreal());
#elif ( DIM == 2 )
        res.set((length+radius)*RNG.sreal(), radius*RNG.sreal());
#else
        Vector2 sec = Vector2::randB(radius);
        res.set((length+radius)*RNG.sreal(), sec.XX, sec.YY);
#endif
        
        if ( ++ouf > nb_trials )
        {
            std::clog << "placement failed in SpaceCapsule::randomPlace()" << std::endl;
            return Vector(0,0,0);
        }
        
    } while ( ! inside(res) );
    
    return res;
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff, const real len, const real rad)
{
    if ( fabs(pos.XX) > len )
        meca.addSphereClamp(pos, pe, Vector(std::copysign(len, pos.XX),0,0), rad, stiff);
    else
        meca.addCylinderClampX(pe, rad, stiff);
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length, radius);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCapsule::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    if ( rad < radius )
        setInteraction(pos, pe, meca, stiff, length, radius-rad);
    else
        setInteraction(pos, pe, meca, stiff, length, 0);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCapsule::draw() const
{
    //number of sections in the quarter-circle
    constexpr size_t fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    
    GLfloat c[4*fin+1], s[4*fin+1];
    gle::circle(4*fin, c, s, 1);
    
    GLfloat L = (GLfloat)length;
    GLfloat R = (GLfloat)radius;
    
#if ( DIM <= 2 )
    
    //display a loop in X/Y plane
    glBegin(GL_LINE_LOOP);
    
    for ( size_t n = 0;     n <= 2*fin; ++n )
        glVertex2f(R*s[n]+L, R*c[n]);
    
    for ( size_t n = 2*fin; n <= 4*fin; ++n )
        glVertex2f(R*s[n]-L, R*c[n]);
    
    glEnd();
    
#else
    
    //display strips along the side of the volume:
    for ( size_t sc = 0; sc < 4*fin; ++sc )
    {
        //compute the transverse angles:
        GLfloat ctb  = c[sc  ],   stb  = s[sc  ];
        GLfloat cta  = c[sc+1],   sta  = s[sc+1];
        GLfloat ctbR = R*ctb,     stbR = R*stb;
        GLfloat ctaR = R*cta,     staR = R*sta;
        
        //draw one strip of the oval:
        glBegin(GL_TRIANGLE_STRIP);
        for ( size_t ii=0; ii <= fin; ++ii )
        {
            GLfloat ca = c[ii], sa = s[ii];
            glNormal3f( ca, cta*sa, sta*sa );
            glVertex3f( +L+R*ca, ctaR*sa, staR*sa );
            glNormal3f( ca, ctb*sa, stb*sa );
            glVertex3f( +L+R*ca, ctbR*sa, stbR*sa );
        }
        for ( int ii=fin; ii >= 0; --ii)
        {
            GLfloat ca = -c[ii], sa = s[ii];
            glNormal3f( ca, cta*sa, sta*sa );
            glVertex3f( -L+R*ca, ctaR*sa, staR*sa );
            glNormal3f( ca, ctb*sa, stb*sa );
            glVertex3f( -L+R*ca, ctbR*sa, stbR*sa );
        }
        glEnd();
    }
    
    if ( 1 )
    {
        //draw 2 rings on the surface
        glPushMatrix();
        glTranslatef(L, 0, 0);
        glScalef(R, R, R);
        glRotated(90, 0, 1, 0);
        gle::gleArrowedBand(24, 0.25);
        glTranslatef(0, 0, -2*L/R);
        glRotated(180, 0, 1, 0);
        gle::gleArrowedBand(24, 0.25);
        glPopMatrix();
    }

#endif
    return true;
}

#else

bool SpaceCapsule::draw() const
{
    return false;
}

#endif
