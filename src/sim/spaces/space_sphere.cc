// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_sphere.h"
#include "exceptions.h"
#include "random.h"
#include "meca.h"


SpaceSphere::SpaceSphere(const SpaceProp* p)
: Space(p), radius(mLength[0]), radiusSqr(mLengthSqr[0])
{
}


void SpaceSphere::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius,-radius,-radius);
    sup.set( radius, radius, radius);
}


real SpaceSphere::volume() const
{
#if ( DIM == 1 )
    return 2 * radius;
#elif ( DIM == 2 )
    return M_PI * radius * radius;
#else
    return 4/3.0 * M_PI * radius * radius * radius;
#endif
}

bool SpaceSphere::inside(Vector const& pos) const
{
    return pos.normSqr() <= radiusSqr;
}

Vector SpaceSphere::project(Vector const& pos) const
{
    real n = pos.normSqr();
    
    if ( n > 0 ) {
        return pos * ( radius / sqrt(n) );
    }
    else {
        //select a random point on the surface
        return radius * Vector::randU();
    }
}

//------------------------------------------------------------------------------

void SpaceSphere::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    meca.addSphereClamp( pos, pe, Vector(0,0,0), radius, stiff );
}


void SpaceSphere::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    if ( radius > rad )
        meca.addSphereClamp( pos, pe, Vector(0,0,0), radius-rad, stiff );
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceSphere\n";
    }
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"

bool SpaceSphere::draw() const
{

#if ( DIM <= 2 )
 
    //number of sections in the quarter-circle
    constexpr size_t fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    
    GLfloat cir[2*fin+2];
    gle::circle(fin, cir, (GLfloat)radius);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);
    glDisableClientState(GL_VERTEX_ARRAY);

#else
    
    GLfloat R = (GLfloat)radius;
    glPushMatrix();
    glScalef(R, R, R);
    gle::gleSphere8B();
    gle::gleThreeBands(128);
    glPopMatrix();
    
#endif
    
    return true;
}

#else

bool SpaceSphere::draw() const
{
    return false;
}


#endif
