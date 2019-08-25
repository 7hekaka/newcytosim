// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "wrist_long.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"


extern Modulo const* modulo;


WristLong::WristLong(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Wrist(sp, mec, pti), mArm(nullTorque)
{
}


WristLong::~WristLong()
{
}

//------------------------------------------------------------------------------

Torque WristLong::calcArm(const Interpolation & pt, Vector const& pos, real len)
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


Vector WristLong::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - WristLong::posSide();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


#if ( 1 )

Vector WristLong::posSide() const
{
#if ( DIM > 1 )
    return sHand->pos() + cross(mArm, sHand->dirFiber());
#else
    return sHand->pos();
#endif
}

/**
 Using a Meca::interSideLink()
 */
void WristLong::setInteractions(Meca & meca) const
{
    Interpolation const& pt = sHand->interpolation();
    
    /* 
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    
#if ( DIM == 2 )
    
    mArm = calcArm(pt, posFoot(), prop->length);
    meca.addSideLink2D(pt, anchor.point(), mArm, prop->stiffness);
    
#elif ( DIM >= 3 )
    
    mArm = calcArm(pt, posFoot(), prop->length);
    meca.addSideLink3D(pt, anchor.point(), mArm, prop->stiffness);
    
#endif
}

#else

Vector WristLong::posSide() const
{
    return sHand->pos() + mArm;
}

/** 
 This uses Meca::interLongLink() 
 */
void WristLong::setInteractions(Meca & meca) const
{
    Interpolation const& pt = sHand->interpolation();
    mArm = ( sBase.pos() - sHand->pos() ).normalized(prop->length);
    meca.addLongLink(pt, sBase, prop->length, prop->stiffness);
}

#endif


