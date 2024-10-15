// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "crosslink_long.h"
#include "crosslink_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "meca.h"


CrosslinkLong::CrosslinkLong(CrosslinkProp const* p, Vector const& w)
: Crosslink(p, w), mArm1(nullTorque), mArm2(nullTorque)
{
}


CrosslinkLong::~CrosslinkLong()
{
}


void CrosslinkLong::stepAA()
{
    Vector f = CrosslinkLong::force();
    real mag = f.norm();
    
    if ( cHand1->checkKramersDetachment(mag) )
        cHand1->detach();
    else
        cHand1->stepLoaded( f);
    
    if ( cHand2->checkKramersDetachment(mag) )
        cHand2->detach();
    else
        cHand2->stepLoaded(-f);
}

//------------------------------------------------------------------------------

/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector CrosslinkLong::sidePos1() const
{
#if ( DIM > 1 )
    return cHand1->pos() + cross(mArm1, cHand1->dirFiber());
#else
    return cHand1->pos();
#endif
}


/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector CrosslinkLong::sidePos2() const
{
#if ( DIM > 1 )
    return cHand2->pos() + cross(mArm2, cHand2->dirFiber());
#else
    return cHand2->pos();
#endif
}


/**
 Calculates the force from stiffness * ( sidePos2() - sidePos1() )
 */
Vector CrosslinkLong::force() const
{
    Vector d = CrosslinkLong::sidePos2() - CrosslinkLong::sidePos1();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop()->stiffness * d;
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
    real len = 0.5 * prop()->length;

#if ( DIM == 2 )
    
    Vector dir = pt2.pos() - pt1.pos();
    mArm1 = std::copysign(len, cross(pt1.diff(), dir));
    mArm2 = std::copysign(len, cross(dir, pt2.diff()));
    meca.addSideSideLink(pt1, mArm1, pt2, mArm2, prop()->stiffness);

#elif ( DIM >= 3 )
    
    Vector dir = pt2.pos() - pt1.pos();
    if ( modulo )
        modulo->fold(dir);
    
    Vector off1 = cross(pt1.diff(), dir);
    real n1 = off1.norm();
    if ( n1 > REAL_EPSILON )
        mArm1 = off1 * ( len / n1 );
    else
        mArm1 = pt1.dir().randOrthoU(len);
    
    Vector off2 = cross(dir, pt2.diff());
    real n2 = off2.norm();
    if ( n2 > REAL_EPSILON )
        mArm2 = off2 * ( len / n2 );
    else
        mArm2 = pt1.dir().randOrthoU(len);

    //mArm1 = cross(pt1.diff(), dir).normalized(len);
    //mArm2 = cross(dir, pt2.diff()).normalized(len);
    meca.addSideSideLink(pt1, mArm1, pt2, mArm2, prop()->stiffness);
    
#endif
}

