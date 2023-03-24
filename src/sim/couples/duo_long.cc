// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "duo_long.h"
#include "duo_prop.h"
#include "modulo.h"
#include "random.h"
#include "meca.h"


DuoLong::DuoLong(DuoProp const* p, Vector const& w)
: Duo(p, w), mArm(nullTorque)
{
}


DuoLong::~DuoLong()
{
}


void DuoLong::stepAA()
{
    Vector f = DuoLong::force();
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

Torque DuoLong::calcArm(Interpolation const& pt, Vector const& pos, real len)
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
        return pt.diff().randOrthoU(len);
#else
    return std::copysign(len, cross(off, pt.diff()));
#endif
}


/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector DuoLong::sidePos1() const
{
#if ( DIM > 1 )
    return cHand1->pos() + cross(mArm, cHand1->dirFiber());
#else
    return cHand1->pos();
#endif
}


/**
 Calculates the force for the addSideLink()
 */
Vector DuoLong::force() const
{
    Vector d = cHand2->pos() - DuoLong::sidePos1();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop()->stiffness * d;
}


/**
 This uses addSideLink2D() or addSideLink3D().
 
 Another possibility would be addSideSideLink, which is fully symmetric.
 */
void DuoLong::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
    //meca.addSideSideLink(pt1, pt2, prop()->length, prop()->stiffness);
    /*
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    
#if ( DIM == 2 )
    
    mArm = calcArm(pt1, pt2.pos(), prop()->length);
    meca.addSideLink2D(pt1, pt2, mArm, prop()->stiffness);
    
#elif ( DIM >= 3 )
    
    mArm = calcArm(pt1, pt2.pos(), prop()->length);
    meca.addSideLink3D(pt1, pt2, mArm, prop()->stiffness);
    
#endif
    
#if NEW_DUO_HAS_TORQUE
#if ( DIM == 2 )
    if ( prop()->flip )
    {
        Vector2 dir = prop()->rest_dir;
        // flip the angle to match the current configuration of the bond
        sine = std::copysign(dir.YY, cross(pt1.diff(), pt2.diff()));
        dir.YY = sine;
        meca.addTorque(pt1, pt2, dir, prop()->angular_stiffness);
    }
    else
    {
        meca.addTorque(pt1, pt2, prop()->rest_dir, prop()->angular_stiffness);
    }
    //meca.addTorquePoliti(pt1, pt2, dir, prop()->angular_stiffness);
#elif ( DIM == 3 )
    meca.addTorque(pt1, pt2, prop()->rest_dir, prop()->angular_stiffness);
#endif
#endif
}

