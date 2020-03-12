// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "space_dynamic_sphere.h"
#include "meca.h"


SpaceDynamicSphere::SpaceDynamicSphere(DynamicSpaceProp const* p)
: SpaceSphere(p),prop(p)
{
    force_ = 0;
}


void SpaceDynamicSphere::setInteractions(Meca&) const
{
    force_ = 0;
}


void SpaceDynamicSphere::setInteraction(Vector const& pos, Mecapoint const& pe,
                                        Meca& meca, real stiff) const
{
    meca.addSphereClamp(pos, pe, Vector(0,0,0), radius_, stiff);
    force_ += stiff * ( pos.norm() - radius_ );
}


void SpaceDynamicSphere::setInteraction(Vector const& pos, Mecapoint const& pe,
                                        real rad, Meca& meca, real stiff) const
{
    if ( radius_ > rad )
    {
        meca.addSphereClamp(pos, pe, Vector(0,0,0), radius_-rad, stiff);
        force_ += stiff * ( rad + pos.norm() - radius_ );
    }
    else {
        meca.addPointClamp( pe, Vector(0,0,0), stiff );
        std::cerr << "object is too big to fit in SpaceDynamicSphere\n";
        force_ += 2 * stiff * ( rad - radius_ );
    }
}


void SpaceDynamicSphere::step()
{
    real dr = prop->mobility_dt * force_;
    std::clog << "SpaceDynamicSphere:  radius " << std::setw(12) << radius_;
    std::clog << " force " << force_ << " delta_radius " << dr << "\n";
    radius_ += dr;
}

