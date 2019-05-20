// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "kinesin.h"
#include "kinesin_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "lattice.h"
#include "simul.h"


Kinesin::Kinesin(KinesinProp const* p, HandMonitor* h)
: Digit(p,h), prop(p)
{
    nextStep = RNG.exponential();
    ABORT_NOW("unfinished class");
}


/**
 \todo simulate occurence of backward steps
 */
void Kinesin::stepUnloaded()
{
    assert_true( attached() );
    
    nextStep -= prop->stepping_rate_dt;
    
    while ( nextStep <= 0 )
    {
        assert_true( attached() );
        site_t s = site() + 1;
        if ( edgy(s) )
        {
            nextStep = RNG.exponential();
            //immediately detach at the end of the Fiber:
            detach();
            return;
        }
        if ( vacant(s) )
            hop(s);
        nextStep += RNG.exponential();
    }
    
    testDetachment();
}


/**
 Currently, antagonistic force only reduced the rate of forward stepping.
 However, force is also known to increase the rate of backward steps.
 \todo simulate occurence of backward steps
 */
void Kinesin::stepLoaded(Vector const& force, real force_norm)
{
    assert_true( attached() );
    
    // calculate displacement, dependent on the load along the desired direction of displacement
    real rate_step = prop->stepping_rate_dt + dot(force, dirFiber()) * prop->var_rate_dt;

    nextStep -= rate_step;
    
    while ( nextStep <= 0 )
    {
        assert_true( attached() );
        site_t s = site() + 1;
        if ( edgy(s) )
        {
            nextStep = RNG.exponential();
            //immediately detach at the end of the Fiber:
            detach();
            return;
        }
        if ( vacant(s) )
            hop(s);
        nextStep += RNG.exponential();
    }
    
    assert_true( nextDetach >= 0 );
    testKramersDetachment(force_norm);
}

