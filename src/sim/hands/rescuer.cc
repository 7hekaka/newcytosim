// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "rescuer.h"
#include "rescuer_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"


Rescuer::Rescuer(RescuerProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
}


/**
 Warning:
 This will only work if the time step is small such that only one Hand is
 affected at any time step by the shrinkage of the Fiber.
 Otherwise, the order in which the Hand are considered is random,
 and a distal Hand might be detached, even if the Fiber is rescued prior to this.
 
 This condition should be satisfied for Microtubule systems, since:
 - Shrinking speed ~ -0.1 um/second
 - time_step = 0.010 seconds
 -> shrinkage = 1 nm / time_step
 which would work for the density of 13 binding sites / 8 nm of microtubule.
 */
void Rescuer::handleDisassemblyM()
{
    Fiber * fib = fiber();

    if ( hAbs < fib->abscissaM() && fib->freshAssemblyM() < 0 )
    {
        if ( RNG.test(prop->rescue_prob) )
        {
            //revert the last disassembly step:
            fib->growM(-fib->freshAssemblyM());
            // induce rescue:
            fib->setEndStateM(STATE_GREEN);
        }
        else
            detach();
    }
}

/**
 Warning:
 This will work best if the time step is small such that only one Hand is
 affected at any time step by the shrinkage of the Fiber.
 Otherwise, the order in which the Hand are considered is random,
 and a distal Hand might be detached, even if the Fiber is rescued prior to this.
 
 This condition should be satisfied for Microtubule systems, since:
 - Shrinking speed ~ -0.1 um/second
 - time_step = 0.010 seconds
 -> shrinkage = 1 nm / time_step
 With a density of 13 binding sites / 8 nm of microtubule, this is at most
 2 events / time_step
 */
void Rescuer::handleDisassemblyP()
{
    Fiber * fib = fiber();

    /* increase MT length to cover position of Hand after checking that increment
    would be positive, since another Hand might have done the same job already. */
    if ( hAbs > fib->abscissaP() && fib->freshAssemblyP() < 0 )
    {
        if ( RNG.test(prop->rescue_prob) )
        {
            //revert the last disassembly step:
            fib->growP(-fib->freshAssemblyP());
            // induce rescue:
            fib->setEndStateP(STATE_GREEN);
            //assert_true(hAbs <= fib->abscissaP())
        }
        else
            detach();
    }
}


void Rescuer::stepUnloaded()
{
    assert_true( attached() );
}


void Rescuer::stepLoaded(Vector const& force)
{
    assert_true( attached() );
}

