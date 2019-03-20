// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_dice.h"
#include "exceptions.h"


SpaceDice::SpaceDice(const SpaceProp* p)
: Space(p), radius(mLength[3]), radiusSqr(mLengthSqr[3])
{
    if ( DIM == 1 )
        throw InvalidParameter("dice is not usable in 1D");
}


void SpaceDice::resize()
{
    Space::checkLengths(4, false);
    
    for ( int d = 0; d < 3; ++d )
        if ( length(d) < radius )
            throw InvalidParameter("all dice's dimensions must be >= radius");
}


/**
 The `dice` is included in the rectangle
 */
void SpaceDice::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length(0),-length(1),-length(2));
    sup.set( length(0), length(1), length(2));
}


/**
 If `radius==0`, the volume should be the volume of a rectangle
 */
real SpaceDice::volume() const
{
#if ( DIM == 1 )
    return 2 * length(0);
#elif ( DIM == 2 )
    return 4 * length(0)*length(1) + (M_PI-4)*radius*radius;
#else
    return 8 * length(0)*length(1)*length(2)
    + 2 * (M_PI-4) * ( length(0) + length(1) + length(2) - 3 * radius ) * radius * radius
    + (4/3.0 * M_PI - 8) * radius * radius * radius;
#endif
}


//------------------------------------------------------------------------------

bool  SpaceDice::inside(Vector const& w) const
{
    real dis = 0;
    for ( int d = 0; d < DIM; ++d )
    {
        real a = fabs(w[d]) - length(d);
        if ( a > 0 )
            return false;
        a += radius;
        if ( a > 0 )
            dis += a * a;
    }
    return ( dis <= radiusSqr );
}



//------------------------------------------------------------------------------

#if ( DIM == 1 )

Vector SpaceDice::project(Vector const& w) const
{
    return Vector(std::copysign(length(0), w.XX), 0, 0);
}

#else

Vector SpaceDice::project(Vector const& w) const
{
    Vector p = w;
    bool in = true;
    
    //calculate projection on the inner cube obtained by subtracting radius
    for ( int d = 0; d < DIM; ++d )
    {
        real test = length(d) - radius;
        if ( fabs(w[d]) > test )
        {
            p[d] = std::copysign(test, w[d]);
            in = false;
        }
    }
    
    if ( in )
    {
        // find the dimensionality corresponding to the closest face
        real d0 = length(0) - fabs(w.XX);
        real d1 = length(1) - fabs(w.YY);
#if ( DIM > 2 )
        real d2 = length(2) - fabs(w.ZZ);
        if ( d2 < d1 )
        {
            if ( d2 < d0 )
                p.ZZ = std::copysign(length(2), w.ZZ);
            else
                p.XX = std::copysign(length(0), w.XX);
        }
        else
#endif
        {
            if ( d1 < d0 )
                p.YY = std::copysign(length(1), w.YY);
            else
                p.XX = std::copysign(length(0), w.XX);
        }
        return p;
    }

    //normalize to radius(), and add to p to get the real projection
    real dis = radius / sqrt((w-p).normSqr());
    for ( int d = 0; d < DIM; ++d )
        p[d] += dis * ( w[d] - p[d] );
    
    return p;
}
#endif

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceDice::draw() const
{
#if ( DIM > 2 )
    
    const real X = length(0) - radius;
    const real Y = length(1) - radius;
    const real Z = length(2) - radius;
 
    const real XR = length(0);
    const real YR = length(1);
    const real ZR = length(2);

    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  XR,  Y, -Z );
    gleVertex(  XR,  Y,  Z );
    gleVertex(  XR, -Y, -Z );
    gleVertex(  XR, -Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -XR, -Y, -Z );
    gleVertex( -XR, -Y,  Z );
    gleVertex( -XR,  Y, -Z );
    gleVertex( -XR,  Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  YR, -Z );
    gleVertex( -X,  YR, -Z );
    gleVertex(  X,  YR,  Z );
    gleVertex( -X,  YR,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X, -YR,  Z );
    gleVertex( -X, -YR,  Z );
    gleVertex(  X, -YR, -Z );
    gleVertex( -X, -YR, -Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y,  ZR );
    gleVertex( -X,  Y,  ZR );
    gleVertex(  X, -Y,  ZR );
    gleVertex( -X, -Y,  ZR );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, -ZR );
    gleVertex(  X, -Y, -ZR );
    gleVertex( -X,  Y, -ZR );
    gleVertex( -X, -Y, -ZR );
    glEnd();
    
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    drawSection( 0, -X, 0.01 );
    drawSection( 0,  X, 0.01 );
    drawSection( 1, -Y, 0.01 );
    drawSection( 1,  Y, 0.01 );
    drawSection( 2, -Z, 0.01 );
    drawSection( 2,  Z, 0.01 );
    glDisable(GL_LINE_STIPPLE);
    glPopAttrib();
    
#else

    drawSection( 2, 0, 0.01 );

#endif
    
    return true;
}

#else

bool SpaceDice::draw() const
{
    return false;
}

#endif


