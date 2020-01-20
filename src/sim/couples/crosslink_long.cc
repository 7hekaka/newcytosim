// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "crosslink_long.h"
#include "crosslink_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "meca.h"

extern Modulo const* modulo;


CrosslinkLong::CrosslinkLong(CrosslinkProp const* p, Vector const& w)
: Crosslink(p, w), mArm1(nullTorque), mArm2(nullTorque)
{
}


CrosslinkLong::~CrosslinkLong()
{
}

/*
 Note that, since `mArm` is calculated by setInteraction(),
 the result of sidePos() will be incorrect if 'solve=0'
*/
Vector CrosslinkLong::sidePos() const
{
#if ( DIM > 1 )
    return cHand1->pos() + cross(mArm1, cHand1->dirFiber());
#else
    return cHand1->pos();
#endif
}


/**
 Calculates the force for the interSideLink()
 */
Vector CrosslinkLong::force() const
{
    Vector d = cHand2->pos() - CrosslinkLong::sidePos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


/**
 This uses Meca::addSideSideLink(), which is fully symmetric.
 */
void CrosslinkLong::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    /*
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    real len = 0.5 * prop->length;

#if ( DIM == 2 )
    
    Vector dir = pt2.pos() - pt1.pos();
    mArm1 = std::copysign(len, cross(pt1.diff(), dir));
    mArm2 = std::copysign(len, cross(dir, pt2.diff()));
    meca.addSideSideLink2D(pt1, pt2, mArm1, mArm2, prop->stiffness);

#elif ( DIM >= 3 )

# if FIBER_HAS_FAMILY
    Vector rad1 = fiber1()->radialDir(abscissa1());
    Vector rad2 = fiber2()->radialDir(abscissa2());
    mArm1 = cross(pt1.diff(), rad1).normalized(len);
    mArm2 = cross(pt2.diff(), rad2).normalized(len);
# else
    Vector dir = ptB.pos() - ptA.pos();
    mArm1 = cross(ptA.diff(), dir).normalized(len);
    mArm2 = cross(dir, ptB.diff()).normalized(len);
# endif
    meca.addSideSideLink3D(pt1, pt2, mArm1, mArm2, prop->stiffness);
    
#endif
}

