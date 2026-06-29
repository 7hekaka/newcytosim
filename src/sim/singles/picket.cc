// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "picket.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"
#include "random.h"
#include "space.h"
#include <cmath>

namespace
{
    static inline Vector radial_xy(Vector const& p)
    {
        const real r = std::sqrt(p.XX * p.XX + p.YY * p.YY);
        if ( r > REAL_EPSILON )
            return Vector(p.XX / r, p.YY / r, 0);
        return Vector(1, 0, 0);
    }

    static inline void slide_surface_anchor(SingleProp const* prop, Simul const& sim, Vector& pos)
    {
        if ( prop->anchor_mode != "slide" || prop->anchor_D <= 0 )
            return;
        if ( prop->confine != CONFINE_ON || !prop->confine_space )
            return;

        const real dt = sim.prop.time_step;
        if ( dt <= 0 )
            return;

        const Vector n = radial_xy(pos);
        const Vector et(-n.YY, n.XX, 0);
        const Vector ez(0, 0, 1);
        const real sigma = std::sqrt(2.0 * prop->anchor_D * dt);

        pos += sigma * ( RNG.gauss() * et + RNG.gauss() * ez );
        pos = prop->confine_space->project(pos);
    }
}



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

    slide_surface_anchor(prop, simul(), sPos);
    sHand->stepUnattached(simul(), sPos);
}


void Picket::stepA()
{
    assert_true( sHand->attached() );
    
    slide_surface_anchor(prop, simul(), sPos);
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

