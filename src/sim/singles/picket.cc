// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "picket.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"



Picket::Picket(SingleProp const* p, Vector const& w)
: Single(p, w)
{
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter(name()+":diffusion cannot be > 0 if activity=fixed");
#endif
}


Picket::~Picket()
{
    //std::clog<<"~Picket("<<this<<")\n";
}


void Picket::beforeDetachment(Hand const*)
{
    assert_true( attached() );

    SingleSet * set = static_cast<SingleSet*>(objset());
    if ( set )
        set->relinkD(this);
}


void Picket::stepF()
{
    assert_false( sHand->attached() );

    sHand->stepUnattached(simul(), sPos);
}


void Picket::stepA()
{
    assert_true( sHand->attached() );
    
    Vector f = Picket::force();
    if ( sHand->checkKramersDetachment(f.norm()) )
        sHand->detach();
    else
        sHand->stepLoaded(f);
}


/**
 This calculates the force corresponding to addPointClamp()
 */
Vector Picket::stretch() const
{
    assert_true( sHand->attached() );
    Vector d = sPos - posHand();
    
    if ( modulo )
        modulo->fold(d);
    
    return d;
}

/**
 This calculates the force corresponding to addPointClamp()
 */
Vector Picket::force() const
{
    assert_true( sHand->attached() );
    Vector d = sPos - posHand();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void Picket::setInteractions(Meca& meca) const
{
    assert_true( prop->length == 0 );
    meca.addPointClamp(sHand->interpolation(), sPos, prop->stiffness);
    //meca.addLineClamp(sHand->interpolation(), sPos, sHand->dirFiber(), prop->stiffness);
}


