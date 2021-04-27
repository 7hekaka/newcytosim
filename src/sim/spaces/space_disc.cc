// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "space_disc.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"


SpaceDisc::SpaceDisc(SpaceDynamicProp const* p)
: Space(p),prop(p)
{
    if ( DIM != 2 )
        throw InvalidParameter("disc is only usable in 2D");
    radius_ = 0;
    force_  = 0;
}


void SpaceDisc::resize(Glossary& opt)
{
    real rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");

    if ( rad < 0 )
        throw InvalidParameter("disc:radius must be >= 0");
    
    radius_ = rad;
}


void SpaceDisc::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_,-radius_,-radius_);
    sup.set( radius_, radius_, radius_);
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
    return M_PI * radius_ * radius_;
}

bool SpaceDisc::inside(Vector const& pos) const
{
    return pos.normSqr() <= radius_ * radius_;
}

Vector SpaceDisc::project(Vector const& pos) const
{
    real n = pos.normSqr();
    
    if ( n > 0 ) {
        return pos * ( radius_ / std::sqrt(n) );
    }
    else {
        //select a random point on the surface
        return radius_ * Vector::randU();
    }
}

#endif

//------------------------------------------------------------------------------

/// add interactions to a Meca
void SpaceDisc::setInteractions(Meca&) const
{
    force_ = 0;
}


void SpaceDisc::setConfinement(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    meca.addSphereClamp(pos, pe, Vector(0,0,0), radius_, stiff);
    force_ += stiff * ( pos.norm() - radius_ );
}


void SpaceDisc::setConfinement(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    if ( radius_ > rad )
    {
        meca.addSphereClamp(pos, pe, Vector(0,0,0), radius_-rad, stiff);
        force_ += stiff * ( rad + pos.norm() - radius_ );
    }
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceDisc\n";
        force_ += 2 * stiff * ( rad - radius_ );
    }
}


void SpaceDisc::step()
{
    real dr = prop->mobility_dt * force_;
    //std::clog << "SpaceDisc:  radius " << std::setw(12) << radius_ << " force " << force_ << " delta_radius " << dr << "\n";
    radius_ += dr;
}

//------------------------------------------------------------------------------

void SpaceDisc::write(Outputter& out) const
{
    writeShape(out, "disc");
    out.writeUInt16(2);
    out.writeFloat(radius_);
    out.writeFloat(force_);
}


void SpaceDisc::setLengths(const real len[])
{
    radius_ = len[0];
    force_  = len[1];
}

void SpaceDisc::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "disc");
    setLengths(len);
}


#ifdef DISPLAY

#include "opengl.h"
#include "gle.h"

void SpaceDisc::draw2D() const
{
    GLfloat R(radius_);
    glPushMatrix();
    glScalef(R, R, R);
    glPopMatrix();
}

void SpaceDisc::draw3D() const
{
    draw2D(); // unfinished
}

#else

void SpaceDisc::draw2D() const {}
void SpaceDisc::draw3D() const {}

#endif
