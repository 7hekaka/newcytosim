// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "space_dynamic_sphere.h"
#include "meca.h"


SpaceDynamicSphere::SpaceDynamicSphere(const SpaceProp* p)
: SpaceSphere(p)
{
    rForce = 0;
}


void SpaceDynamicSphere::setInteractions(Meca &, FiberSet const&) const
{
    rForce = 0;
}


void SpaceDynamicSphere::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    meca.addSphereClamp(pos, pe, Vector(0,0,0), radius, stiff);
    rForce += stiff * ( pos.norm() - radius );
}


void SpaceDynamicSphere::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    if ( radius > rad )
    {
        meca.addSphereClamp(pos, pe, Vector(0,0,0), radius-rad, stiff);
        rForce += stiff * ( rad + pos.norm() - radius );
    }
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceDynamicSphere\n";
        rForce += 2 * stiff * ( rad - radius );
    }
}


void SpaceDynamicSphere::step()
{
    real dr = prop->mobility_dt * rForce;
    std::clog << "SpaceDynamicSphere:  radius " << std::setw(12) << radius << " force " << rForce << " delta_radius " << dr << "\n";
    radius += dr;
}

