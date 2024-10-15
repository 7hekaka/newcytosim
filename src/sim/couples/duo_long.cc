// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "duo_long.h"
#include "duo_prop.h"
#include "modulo.h"
#include "random.h"
#include "object_set.h"
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
    if ( active_ && prop()->vulnerable )
        tryDeactivate();

    if ( active_ )
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
    else
    {
        cHand1->detach();
        cHand2->detach();
        if ( prop()->deactivation_mode )
            recycle();
        return;
    }
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
        return pt.dir().randOrthoU(len);
#else
    return std::copysign(len, cross(off, pt.diff()));
#endif
}



void DuoLong::afterAttachment(Hand const* ha)
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
}

