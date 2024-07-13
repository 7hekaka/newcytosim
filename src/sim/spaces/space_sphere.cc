// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_sphere.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"

SpaceSphere::SpaceSphere(SpaceProp const* p)
: Space(p), radius_(0)
{
}


void SpaceSphere::resize(Glossary& opt)
{
    real rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else if ( opt.set(rad, "surface") )
        rad = std::sqrt(rad*0.25*M_1_PI);
    else opt.set(rad, "radius");
    
    if ( rad < 0 )
        throw InvalidParameter(prop->name()+":radius must be >= 0");
    
    radius_ = rad;
}

void SpaceSphere::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_,-radius_,-radius_);
    sup.set( radius_, radius_, radius_);
}


real SpaceSphere::volume() const
{
#if ( DIM == 1 )
    return 2 * radius_;
#elif ( DIM == 2 )
    return M_PI * square(radius_);
#else
    return ( 4 * M_PI / 3.0 ) * cube(radius_);
#endif
}


real SpaceSphere::surface() const
{
#if ( DIM == 1 )
    return 2;
#elif ( DIM == 2 )
    return ( 2 * M_PI ) * radius_;
#else
    return ( 4 * M_PI ) * square(radius_);
#endif
}


bool SpaceSphere::inside(Vector const& pos) const
{
    return pos.normSqr() <= square(radius_);
}


bool SpaceSphere::allInside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
    real R = max_real(0, radius_-rad);  // remaining radius
    return pos.normSqr() <= R * R;
}


Vector SpaceSphere::project(Vector const& pos) const
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

//------------------------------------------------------------------------------

void SpaceSphere::setConfinement(Vector const& pos, Mecapoint const& mp,
                                 Meca& meca, real stiff) const
{
    meca.addSphereClamp(pos, mp, Vector(0,0,0), radius_, stiff);
}


void SpaceSphere::setConfinement(Vector const& pos, Mecapoint const& mp,
                                 real rad, Meca& meca, real stiff) const
{
    if ( radius_ > rad )
        meca.addSphereClamp(pos, mp, Vector(0,0,0), radius_-rad, stiff);
    else {
        meca.addPointClamp(mp, Vector(0,0,0), stiff);
        std::cerr << "object is too big to fit in SpaceSphere\n";
    }
}


//------------------------------------------------------------------------------

void SpaceSphere::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "R");
    out.writeUInt16(2);
    out.writeFloat(radius_);
    out.writeFloat(0.f);
}


void SpaceSphere::setLengths(const real len[8])
{
    radius_ = len[0];
}


void SpaceSphere::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "R");
    setLengths(len);
}


//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "point_disp.h"
#include "gym_view.h"

void SpaceSphere::draw2D(float width) const
{
    if ( prop->disp && prop->disp->visible & 2 )
    {
        gym::color(prop->disp->color2);
        gle::disc1();
    }

    gym::color(prop->disp->color);
    gle::circle(radius_, width);
}

void SpaceSphere::draw3D() const
{
    gym::scale(radius_);
    gle::sphere();
    if ( prop->disp && prop->disp->visible & 4 )
        gle::threeArrowStrip(0.5, 1);
}

#else

void SpaceSphere::draw2D(float) const {}
void SpaceSphere::draw3D() const {}

#endif
