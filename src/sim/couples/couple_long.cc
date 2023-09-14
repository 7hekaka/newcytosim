// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "couple_long.h"
#include "couple_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "meca.h"


CoupleLong::CoupleLong(CoupleProp const* p, Vector const& w)
: Couple(p, w), mArm(nullTorque)
{
}


CoupleLong::~CoupleLong()
{
}


void CoupleLong::stepAA()
{
    Vector f = CoupleLong::force();
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

Torque CoupleLong::calcArm(Interpolation const& pt, Vector const& pos, real len)
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


void CoupleLong::afterAttachment(Hand const* ha)
{
    Couple::afterAttachment(ha);

#if ( DIM > 1 )
    if ( cHand1->attached() && cHand2->attached() )
    {
        Interpolation const& pt1 = cHand1->interpolation();
        Interpolation const& pt2 = cHand2->interpolation();
        mArm = calcArm(pt1, pt2.pos(), prop->length);
    }
#endif
}


/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector CoupleLong::sidePos1() const
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
Vector CoupleLong::force() const
{
    Vector d = cHand2->pos() - CoupleLong::sidePos1();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


/**
 This uses addSideLink2D() or addSideLink3D().
 
 Another possibility would be addSideSideLink, which is fully symmetric.
 */
void CoupleLong::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
    //meca.addSideSideLink(pt1, pt2, prop()->length, prop()->stiffness);
    
    /*
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    
#if ( DIM == 2 )
    
    mArm = calcArm(pt1, pt2.pos(), prop->length);
    meca.addSideLink2D(pt1, pt2, mArm, prop->stiffness);
    
#elif ( DIM >= 3 )

    mArm = calcArm(pt1, pt2.pos(), prop->length);
    meca.addSideLink3D(pt1, pt2, mArm, prop->stiffness);
    
#endif
}


