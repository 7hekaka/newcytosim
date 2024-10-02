// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#include "wrist_long.h"
#include "simul_part.h"
#include "meca.h"
#include "modulo.h"



WristLong::WristLong(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Wrist(sp, mec, pti)
{
#if WRIST_USES_LONGLINK
    mArm = Vector(0,0,0);
#else
    mArm = nullTorque;
#endif
}


WristLong::~WristLong()
{
}

//------------------------------------------------------------------------------

Vector WristLong::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - WristLong::sidePos();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void WristLong::stepA()
{
    assert_true( sHand->attached() );
    Vector f = WristLong::force();
    if ( sHand->checkKramersDetachment(f.norm()) )
        sHand->detach();
    else
        sHand->stepLoaded(f);
}


#if WRIST_USES_LONGLINK

Vector WristLong::calcArm(Interpolation const& itp, Vector const& pos, real len)
{
    Vector off = pos - itp.pos();
    if ( modulo )
        modulo->fold(off);
    return off.normalized(len);
}


void WristLong::afterAttachment(Hand const* ha)
{
    Single::afterAttachment(ha);

#if ( DIM > 1 )
    Interpolation const& ipt = sHand->interpolation();
    mArm = calcArm(ipt, posFoot(), prop->length);
#endif
}


Vector WristLong::sidePos() const
{
    return sHand->pos() + mArm;
}

/**
 This uses Meca::addLongLink()
 */
void WristLong::setInteractions(Meca& meca) const
{
    Interpolation const& itp = sHand->interpolation();
    mArm = calcArm(itp, posFoot(), prop->length);
    
    // 'mArm' is not used to calculate the interaction here:
    if ( base_.rank() == 1 )
        meca.addLongLink(base_.vertex0(), itp, prop->length, prop->stiffness);
    else
    {
        const index_t off = base_.vertex0().matIndex0();
        real const* coef = base_.coef();
        meca.addLongLink4(itp, off, coef[0], coef[1], coef[2], coef[3], prop->length, prop->stiffness);
    }
}

#else

Torque WristLong::calcArm(Interpolation const& itp, Vector const& pos, real len)
{
    Vector off = itp.pos1() - pos;
    if ( modulo )
        modulo->fold(off);
#if ( DIM >= 3 )
    off = cross(off, itp.diff());
    real n = off.norm();
    if ( n > REAL_EPSILON )
        return off * ( len / n );
    else
        return itp.dir().randOrthoU(len);
#else
    return std::copysign(len, cross(off, itp.diff()));
#endif
}

/*
 Note that, since `mArm` is calculated by setInteractions(),
 the result will be incorrect if 'solve=0'
*/
Vector WristLong::sidePos() const
{
#if ( DIM > 1 )
    return sHand->pos() + cross(mArm, sHand->dirFiber());
#else
    return sHand->pos();
#endif
}

/**
 Using Meca::addSideLink2D() or Meca::addSideLink3D()
 */
void WristLong::setInteractions(Meca& meca) const
{
    Interpolation const& itp = sHand->interpolation();
    
    /* 
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    
#if ( DIM == 2 )
    
    mArm = calcArm(itp, posFoot(), prop->length);
    if ( base_.rank() == 1 )
        meca.addSideLink2D(itp, base_.vertex0(), mArm, prop->stiffness);
    else
        throw InvalidParameter("unfinished WristLong::setInteractions(length>0, Interpolation4)");

#elif ( DIM >= 3 )
    
    mArm = calcArm(itp, posFoot(), prop->length);
    if ( base_.rank() == 1 )
        meca.addSideLink3D(itp, base_.vertex0(), mArm, prop->stiffness);
    else
        throw InvalidParameter("unfinished WristLong::setInteractions(length>0, Interpolation4)");
#endif
}

#endif


