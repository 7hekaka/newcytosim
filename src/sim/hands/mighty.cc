// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "mighty.h"
#include "mighty_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "hand_monitor.h"


Mighty::Mighty(MightyProp const* p, HandMonitor* h)
: Hand(p, h)
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
void Mighty::handleDisassemblyM()
{
    Fiber * fib = modifiableFiber();

    if ( hAbs < fib->abscissaM() && fib->freshAssemblyM() < 0 )
    {
        if ( RNG.test(prop()->rescue_chance) )
        {
            //revert the last disassembly step:
            fib->undoGrowM();
            // induce rescue:
            //fib->setEndStateM(STATE_GREEN);
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
void Mighty::handleDisassemblyP()
{
    Fiber * fib = modifiableFiber();

    /* increase MT length to cover position of Hand after checking that increment
    would be positive, since another Hand might have done the same job already. */
    if ( hAbs > fib->abscissaP() && fib->freshAssemblyP() < 0 )
    {
        Hand const* h = otherHand();
        bool link = h && h->attached();
        if ( link && RNG.test(prop()->rescue_chance) )
        {
            //revert the last disassembly step:
            fib->undoGrowP();
            // induce rescue:
            fib->setEndStateP(STATE_GREEN);
            //assert_true(hAbs <= fib->abscissaP())
        }
        else
            detach();
    }
}


void Mighty::stepUnloaded()
{
    assert_true( attached() );
    
    real a = hAbs + prop()->set_speed_dt;
    
    if ( a < hFiber->abscissaM() )
    {
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = hFiber->abscissaM();
    }
    
    if ( a > hFiber->abscissaP() )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = hFiber->abscissaP();
    }

    moveTo(a);
}


void Mighty::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // the load is the projection of the force on the local direction of Fiber
    real load = dot(force, dirFiber());
    
    // calculate load-dependent displacement:
    real dab = prop()->set_speed_dt + load * prop()->var_speed_dt;
    
    // possibly limit the range of the speed:
    if ( prop()->limit_speed )
    {
        dab = std::max(dab, prop()->min_dab);
        dab = std::min(dab, prop()->max_dab);
    }
    
    real a = hAbs + dab;
    
    if ( a < hFiber->abscissaM() )
    {
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = hFiber->abscissaM();
    }
    
    if ( a > hFiber->abscissaP() )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = hFiber->abscissaP();
    }
    
    moveTo(a);
}

