// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fork.h"
#include "fork_prop.h"
#include "meca.h"


Fork::Fork(ForkProp const* p, Vector const& w)
: Couple(p, w), prop(p)
{
#if ( DIM == 2 )
    sinus = prop->angle_dir.YY;
#endif
}


Fork::~Fork()
{
    prop = nullptr;
}


void Fork::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
    meca.addLink(pt1, pt2, prop->stiffness);
    
#if ( DIM == 2 )
    if ( prop->flip )
    {
        Vector2 dir = prop->angle_dir;
        // flip the angle to match the current configuration of the bond
        sinus = std::copysign(dir.YY, cross(pt1.diff(), pt2.diff()));
        dir.YY = sinus;
        meca.addTorque(pt1, pt2, dir, prop->angular_stiffness);
    }
    else
    {
        meca.addTorque(pt1, pt2, prop->angle_dir, prop->angular_stiffness);
    }
    //meca.addTorquePoliti(pt1, pt2, dir, prop->angular_stiffness);
#elif ( DIM == 3 )
    meca.addTorque(pt1, pt2, prop->angle_dir, prop->angular_stiffness);
#endif
}

