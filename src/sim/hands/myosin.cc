// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "myosin.h"
#include "myosin_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "lattice.h"
#include "simul_part.h"
#include "hand_monitor.h"


Myosin::Myosin(MyosinProp const* p, HandMonitor* h)
: Digit(p,h)
{
    ABORT_NOW("unfinished class");
}


void Myosin::attach(FiberSite const& s)
{
    Digit::attach(s);
    nextAct = RNG.exponential();
}


/**
 \todo simulate occurence of backward steps
 */
void Myosin::stepUnloaded()
{
    assert_true( attached() );
    
    nextAct -= prop()->walking_rate_dt;
    
    while ( nextAct <= 0 )
    {
        assert_true( attached() );
        lati_t s = site() + 1;
        if ( outsideMP(s) ) //immediately detach at the end of the Fiber:
            return detach();
        if ( !valLattice(s) )
            hopLattice(s);
        nextAct += RNG.exponential();
    }
}


/**
 Currently, antagonistic force only reduced the rate of forward stepping.
 However, force is also known to increase the rate of backward steps.
 \todo simulate occurence of backward steps
 */
void Myosin::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // calculate displacement, dependent on the load along the desired direction of displacement
    float R = prop()->walking_rate_dt + dot(force, dirFiber()) * prop()->var_rate_dt;

    nextAct -= max_float(0, R);

    while ( nextAct <= 0 )
    {
        assert_true( attached() );
        lati_t s = site() + 1;
        if ( outsideMP(s) )  //immediately detach at the end of the Fiber:
            return detach();
        if ( !valLattice(s) )
            hopLattice(s);
        nextAct += RNG.exponential();
    }
}

