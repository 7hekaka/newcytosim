// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "picket_long.h"
#include "simul_part.h"
#include "meca.h"
#include "modulo.h"



PicketLong::PicketLong(SingleProp const* p, Vector const& w)
: Picket(p, w), mArm(nullTorque)
{
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter(name()+":diffusion cannot be > 0 if activity=fixed");
#endif
}


PicketLong::~PicketLong()
{
    //std::clog<<"~PicketLong("<<this<<")\n";
}

//------------------------------------------------------------------------------

/** This will modify 'pos' if periodic boundary conditions are used */
Torque PicketLong::calcArm(Interpolation const& pt, Vector& pos, real len)
{
    Vector off = pt.pos1() - pos;
    if ( modulo )
        pos += modulo->foldOffset(off);
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


void PicketLong::afterAttachment(Hand const* ha)
{
    Single::afterAttachment(ha);

#if ( DIM > 1 )
    Interpolation const& ipt = sHand->interpolation();
    mArm = calcArm(ipt, sPos, prop->length);
#endif
}

/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector PicketLong::sidePos() const
{
#if ( DIM > 1 )
    return sHand->pos() + cross(mArm, sHand->dirFiber());
#else
    return sHand->pos();
#endif
}


Vector PicketLong::force() const
{
    assert_true( sHand->attached() );
    Vector d = sPos - PicketLong::sidePos();
 
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void PicketLong::stepA()
{
    assert_true( sHand->attached() );
    assert_true( hasLink() );

    Vector f = PicketLong::force();
    if ( sHand->checkKramersDetachment(f.norm()) )
        sHand->detach();
    else
        sHand->stepLoaded(f);
}


void PicketLong::setInteractions(Meca& meca) const
{
#if ( DIM == 1 )
    meca.addPointClamp(sHand->interpolation(), sPos, prop->stiffness);
#else
    Interpolation const& ipt = sHand->interpolation();
    
    /* 
     The 'arm' is recalculated every time, but in 2D at least,
     this may not be necessary, as flipping should occur rarely.
     */
    Vector pos = sPos;
    mArm = calcArm(ipt, pos, prop->length);
    
#if ( DIM == 2 )
    meca.addSidePointClamp2D(ipt, pos, mArm, prop->stiffness);
#elif ( DIM >= 3 )
    meca.addSidePointClamp3D(ipt, pos, mArm, prop->stiffness);
#endif

#endif
}


