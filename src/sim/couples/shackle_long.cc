// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "shackle_long.h"
#include "shackle_prop.h"
#include "modulo.h"
#include "meca.h"


ShackleLong::ShackleLong(ShackleProp const* p, Vector const& w)
: Shackle(p, w), mArm(nullTorque)
{
}

void ShackleLong::stepAA()
{
    Vector f = ShackleLong::force();
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

Torque ShackleLong::calcArm(Interpolation const& pt, Vector const& pos, real len)
{
    Vector off = pt.pos1() - pos;
    if ( modulo )
        modulo->fold(off);
#if ( DIM >= 3 )
    off = cross(off, pt.diff());
    real n = off.norm();
    if ( n > REAL_EPSILON )
        return off * ( len / n );
    else
        return pt.dir().randOrthoU(len);
#else
    return std::copysign(len, cross(off, pt.diff()));
#endif
}


void ShackleLong::afterAttachment(Hand const* ha)
{
    Couple::afterAttachment(ha);

#if ( DIM > 1 )
    if ( cHand1->attached() && cHand2->attached() )
    {
        Interpolation const& pt1 = cHand1->interpolation();
        Interpolation const& pt2 = cHand2->interpolation();
        mArm = calcArm(pt1, pt2.pos(), prop()->length);
    }
#endif
}


/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector ShackleLong::sidePos1() const
{
#if ( DIM > 1 )
    return cHand1->pos() + cross(mArm, cHand1->dirFiber());
#else
    return cHand1->pos();
#endif
}

/**
 Calculates the force as stiffness * ( cHand2->pos() - sidePos1() )
 */
Vector ShackleLong::force() const
{
    Vector d = cHand2->pos() - ShackleLong::sidePos1();
        
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop()->stiffness * d;
}


/**
 The interaction is slipery on hand1
 */
void ShackleLong::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    real iseg = cHand1->fiber()->segmentationInv();

#if ( DIM == 2 )
    
    mArm = calcArm(pt1, pt2.pos(), prop()->length);
    meca.addSideSlidingLink2D(pt1, mArm*iseg, pt2, cHand1->dirFiber(), prop()->stiffness);
    
#elif ( DIM >= 3 )
    
    mArm = calcArm(pt1, pt2.pos(), prop()->length);
    meca.addSideSlidingLink3D(pt1, mArm*iseg, pt2, cHand1->dirFiber(), prop()->stiffness);
    
#endif
}

