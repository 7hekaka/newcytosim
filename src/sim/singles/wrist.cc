// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
#include "wrist.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"


extern Modulo const* modulo;


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, const size_t pti)
: Single(sp)
{
    assert_true(mec);
    rebase(mec, pti);
    assert_false(base_.bad());
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter(name()+":diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, size_t a, size_t b, real c)
: Single(sp)
{
    assert_true(mec);
    rebase(mec, a, b, c);
    assert_false(base_.bad());
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter(name()+":diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::Wrist(SingleProp const* sp, Mecable const* mec, size_t ref, Vector pos)
: Single(sp)
{
    assert_true(mec);
    rebase(mec, ref, pos);
    assert_false(base_.bad());
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
    Vector f = force();
    sHand->stepLoaded(f, f.norm());
}


void Wrist::setInteractions(Meca& meca) const
{
    base_.addLink(meca, sHand->interpolation(), prop->stiffness);
}


void Wrist::write(Outputter& out) const
{
    sHand->write(out);
    base_.write(out);
}


void Wrist::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    const bool s = attached();
    sHand->read(in, sim);
    
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 47 )
    {
        Mecapoint base;
        base.read(in, sim);
        base_.set(base.mecable(), base.point());
    }
    else
#endif
        base_.read(in, sim);
    
    /*
     Because the SingleSet has 2 lists where Single are stored depending
     on their bound/unbound state, we need to unlink and relink here, in
     case the state stored on file is different from the current state.
     */
    if ( s != attached() )
    {
        SingleSet * set = static_cast<SingleSet*>(objset());
        if ( set )
        {
            if ( s )
                set->relinkD(this);
            else
                set->relinkA(this);
        }
    }
}

