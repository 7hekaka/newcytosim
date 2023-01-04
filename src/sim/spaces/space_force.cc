// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_force.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


SpaceForce::SpaceForce(SpaceProp const* p)
: Space(p)
{
    force.reset();
    center.reset();
    stiffness = 0;
}

void SpaceForce::resize(Glossary& opt)
{
    opt.set(force, "force");
    opt.set(center, "center");
    opt.set(stiffness, "stiffness");
}


void SpaceForce::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-1, -1, -1);
    sup.set( 1,  1,  1);
}


real SpaceForce::volume() const
{
    throw InvalidParameter("invalid use of space `force'");
    return -1;
}


Vector SpaceForce::project(Vector const&) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setConfinement(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteractions(Meca& meca, Simul const&) const
{
    if ( stiffness > 0 )
        meca.addPointClampToAll(center, stiffness);
    else
        meca.addForceToAll(force);
}


//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"

void SpaceForce::draw3D() const
{
    gle::drawArrow(center, center+force, 0.1 * force.norm());
}

#else

void SpaceForce::draw3D() const {}

#endif


