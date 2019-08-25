// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "shackle_long.h"
#include "shackle_prop.h"
#include "modulo.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
ShackleLong::ShackleLong(ShackleProp const* p, Vector const& w)
: Shackle(p, w), mArm(nullTorque)
{
}


//------------------------------------------------------------------------------
#if ( DIM == 2 )

/**
 Returns -len or +len
 */
real ShackleLong::calcArm(const Interpolation & pt, Vector const& pos, real len)
{
    Vector off = pt.pos1() - pos;
    if ( modulo )
        modulo->fold(off);
    return std::copysign(len, cross(off, pt.diff()));
}

#elif ( DIM >= 3 )

/**
 Return a vector of norm `len`, perpendicular to the Fiber referenced by `pt` and aligned with the link.
 */
Vector ShackleLong::calcArm(const Interpolation & pt, Vector const& pos, real len)
{
    Vector off = pt.pos1() - pos;
    if ( modulo )
        modulo->fold(off);
    //return cross(off, pt.diff()).normalized(len);
    off = cross(off, pt.diff());
    real n = off.norm();
    if ( n > REAL_EPSILON )
        return off * ( len / n );
    else
        return pt.diff().randOrthoU(len);
}

#endif

//------------------------------------------------------------------------------
Vector ShackleLong::posSide() const
{
#if ( DIM > 1 )
    return cHand1->pos() + cross(mArm, cHand1->dirFiber());
#else
    return cHand1->pos();
#endif
}

/**
 Calculates the force for the interSideLink()
 */
Vector ShackleLong::force() const
{
    Vector d = cHand2->pos() - ShackleLong::posSide();
        
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}

//------------------------------------------------------------------------------
/**
 The interaction is slipery on hand1
 */
void ShackleLong::setInteractions(Meca & meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();

#if ( DIM == 2 )
    
    mArm = calcArm(pt1, pt2.pos(), prop->length);
    meca.addSideSlidingLink2D(pt1, pt2, mArm, prop->stiffness);
    
#elif ( DIM >= 3 )
    
    mArm = calcArm(pt1, pt2.pos(), prop->length);
    meca.addSideSlidingLink3D(pt1, pt2, mArm, prop->stiffness);
    
#endif
}

