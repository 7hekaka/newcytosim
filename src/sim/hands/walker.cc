// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "digit.h"
#include "walker.h"
#include "walker_prop.h"
#include "glossary.h"
#include "messages.h"
#include "lattice.h"
#include "simul_part.h"
#include "hand_monitor.h"


Walker::Walker(WalkerProp const* p, HandMonitor* h)
: Digit(p,h)
{
}


void Walker::attach(FiberSite const& s)
{
    Digit::attach(s);
    nextAct = RNG.exponential();
    
#if ( 0 )
    // this allows for step size being a multiple of lattice site
    unsigned n = std::nearbyint( prop()->step_size / lattice()->unit() );
    prop()->stride = std::copysign(n, prop()->unloaded_speed);
#else
    // here digit::step_size must be equal to fiber:step_size
    if ( lattice() && lattice()->unit() != prop()->step_size  )
        throw InvalidParameter("walker:step_size must be equal to fiber:lattice_unit");
    prop()->stride = (int)sign_real(prop()->unloaded_speed);
#endif
}


/**
 Currently, the Walker only makes forward steps, but backward steps exist as well.
 \todo simulate occurence of backward steps
 */
void Walker::stepUnloaded()
{
    assert_true( attached() );
    
    float R = prop()->walking_rate_dt;
    
    nextAct -= max_float(0, R);

    while ( nextAct <= 0 )
    {
        // test detachment due to stepping
        if ( RNG.test(prop()->unbinding_chance) )
            return detach();
        
        lati_t s = site() + prop()->stride;
        int out = outsideMP(s);
        
        if ( out )
        {
            if ( RNG.test_not(prop()->hold_growing_end[out-1]) )
                return detach();
        }
        else if ( !valLattice(s) )
            hopLattice(s);
    
        nextAct += RNG.exponential();
    }
}


/**
 Currently, antagonistic force only reduces the rate of forward stepping.
 However, force is also known to increase the rate of backward steps.
 \todo simulate occurence of backward steps in Walker
 */
void Walker::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // evaluate displacement, given the load parallel to filament:
    float R = prop()->walking_rate_dt + dot(force, dirFiber()) * prop()->var_rate_dt;

    nextAct -= max_float(0, R);
    
    while ( nextAct <= 0 )
    {
        // test detachment due to stepping
        if ( RNG.test(prop()->unbinding_chance) )
            return detach();

        lati_t s = site() + prop()->stride;
        int out = outsideMP(s);
        
        if ( out )
        {
            if ( RNG.test_not(prop()->hold_growing_end[out-1]) )
                return detach();
        }
        else if ( !valLattice(s) )
            hopLattice(s);
        
        nextAct += RNG.exponential();
    }
}

