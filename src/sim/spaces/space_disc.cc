// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "space_disc.h"
#include "exceptions.h"
#include "random.h"
#include "meca.h"


SpaceDisc::SpaceDisc(const SpaceProp* p)
: Space(p), radius(mLength[0])
{
    if ( DIM != 2 )
        throw InvalidParameter("disc is only usable in 2D");
    rForce = 0;
}


void SpaceDisc::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius,-radius,-radius);
    sup.set( radius, radius, radius);
}


#if ( DIM != 2 )


real SpaceDisc::volume() const
{
    return 0;
}

bool SpaceDisc::inside(Vector const& pos) const
{
    return false;
}

Vector SpaceDisc::project(Vector const&) const
{
    return Vector(0, 0, 0);
}

#else


real SpaceDisc::volume() const
{
    return M_PI * radius * radius;
}

bool SpaceDisc::inside(Vector const& pos) const
{
    return pos.normSqr() <= radius * radius;
}

Vector SpaceDisc::project(Vector const& pos) const
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

#endif

//------------------------------------------------------------------------------

/// add interactions to a Meca
void SpaceDisc::setInteractions(Meca &, FiberSet const&) const
{
    rForce = 0;
}


void SpaceDisc::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    meca.addSphereClamp(pos, pe, Vector(0,0,0), radius, stiff);
    rForce += stiff * ( pos.norm() - radius );
}


void SpaceDisc::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    if ( radius > rad )
    {
        meca.addSphereClamp(pos, pe, Vector(0,0,0), radius-rad, stiff);
        rForce += stiff * ( rad + pos.norm() - radius );
    }
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceDisc\n";
        rForce += 2 * stiff * ( rad - radius );
    }
}


void SpaceDisc::step()
{
    real dr = prop->mobility_dt * rForce;
    //std::clog << "SpaceDisc:  radius " << std::setw(12) << radius << " force " << rForce << " delta_radius " << dr << "\n";
    radius += dr;
}


#ifdef DISPLAY

#include "opengl.h"
#include "gle.h"

bool SpaceDisc::draw() const
{
#if ( DIM <= 2 )

    constexpr size_t fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    GLfloat cir[2*fin+2];
    gle::circle(fin, cir, (GLfloat)radius);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);
    glDisableClientState(GL_VERTEX_ARRAY);
    
#endif
    
    return true;
}

#else

bool SpaceDisc::draw() const
{
    return false;
}


#endif
