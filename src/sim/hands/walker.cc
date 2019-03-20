// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "digit.h"
#include "walker.h"
#include "walker_prop.h"
#include "glossary.h"
#include "lattice.h"
#include "simul.h"


Walker::Walker(WalkerProp const* p, HandMonitor* h)
: Digit(p,h), nextStep(0), prop(p)
{
}


void Walker::attach(FiberSite const& fb)
{
    Digit::attach(fb);
    nextStep = RNG.exponential();
}


int  Walker::stepForward()
{
    if ( prop->unloaded_speed > 0 )
        return Digit::stepP();
    else
        return Digit::stepM();
}


int  Walker::stepBackward()
{
    if ( prop->unloaded_speed > 0 )
        return Digit::stepM();
    else
        return Digit::stepP();
}


/**
 Currently, the Walker only makes forward steps, but backward steps exist as well.
 \todo simulate occurence of backward steps
 */
void Walker::stepUnloaded()
{
    assert_true( attached() );
    
    nextStep   -= prop->stepping_rate_dt;
    
    while ( nextStep <= 0 )
    {
        assert_true( attached() );
        // test detachment due to stepping
        if ( RNG.test(prop->unbinding_chance) )
        {
            detach();
            return;
        }
        
        if ( stepForward() == 2 )
        {
            // we have reached the tip of the fiber
            if ( RNG.test_not(prop->dangling_chance) )
            {
                detach();
                return;
            }
        }
        nextStep += RNG.exponential();
    }
    
    testDetachment();
}


/**
 Currently, antagonistic force only reduces the rate of forward stepping.
 However, force is also known to increase the rate of backward steps.
 \todo simulate occurence of backward steps in Walker
 */
void Walker::stepLoaded(Vector const& force, real force_norm)
{
    assert_true( attached() );
    
    // calculate displacement, dependent on the load along the desired direction of displacement
    real rate_step = prop->stepping_rate_dt + dot(force, dirFiber()) * prop->var_rate_dt;

    nextStep -= rate_step;
    
    while ( nextStep <= 0  &&  attached() )
    {
        // test detachment due to stepping
        if ( RNG.test(prop->unbinding_chance) )
        {
            detach();
            return;
        }
        
        if ( stepForward() == 2 )
        {
            // we have reached the tip of the fiber
            if ( RNG.test_not(prop->dangling_chance) )
            {
                detach();
                return;
            }
        }
        nextStep += RNG.exponential();
    }
    
    assert_true( nextDetach >= 0 );
    if ( prop->unbinding_force_inv > 0 )
        testKramersDetachment(force_norm);
    else
        testDetachment();
}

