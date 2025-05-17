// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "wrist.h"
#include "meca.h"
#include "modulo.h"
#include "single_set.h"



Wrist::Wrist(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Single(sp)
{
    // 'mec' can be Null when reading from file
    rebase(mec, pti);
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter(name()+":diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::~Wrist()
{
}


Vector Wrist::stretch() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - sHand->pos();
    
    if ( modulo )
        modulo->fold(d);
    
    return d;
}


Vector Wrist::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - sHand->pos();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void Wrist::stepF()
{
    assert_false( sHand->attached() );

    sHand->stepUnattached(simul(), posFoot());
}


void Wrist::stepA()
{
    assert_true( sHand->attached() );
    Vector f = Wrist::force();
    if ( sHand->checkKramersDetachment(f.norm()) )
        sHand->detach();
    else
        sHand->stepLoaded(f);
}


void Wrist::setInteractions(Meca& meca) const
{
    Interpolation i = sHand->interpolation();
    base_.addLink(meca, i, prop->stiffness);
}


void Wrist::write(Outputter& out) const
{
    writeMarker(out, WRIST_TAG);
    sHand->writeHand(out);
    base_.write(out);
}


void Wrist::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    sHand->readHand(in, sim);
    
#if BACKWARD_COMPATIBILITY < 47
    if ( in.formatID() < 47 )
    {
        Mecapoint base;
        base.read(in, sim);
        base_.set(base.mecable(), base.point());
    }
    else
#endif
        base_.read(in, sim);
}

