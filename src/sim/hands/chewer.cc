// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "chewer.h"
#include "chewer_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"


Chewer::Chewer(ChewerProp const* p, HandMonitor* h)
: Hand(p,h)
{
    engaged = NO_END;
}


void Chewer::attach(FiberSite const& s)
{
    engaged = NO_END;
    Hand::attach(s);
}


void Chewer::stepUnloaded()
{
    assert_true( attached() );
    
#if ( NEW_FIBER_END_CHEW & 1 )
    if ( engaged == PLUS_END )
    {
        modifiableFiber()->chew(engaged, prop()->chewing_speed_dt);
        moveToEnd(engaged);
        return;
    }
#endif
#if ( NEW_FIBER_END_CHEW & 2 )
    if ( engaged == MINUS_END )
    {
        modifiableFiber()->chew(engaged, prop()->chewing_speed_dt);
        moveToEnd(engaged);
        return;
    }
#endif

    real a = hAbs + prop()->line_diffusion_dt * RNG.sreal();
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        a = M;
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        engaged = MINUS_END;
    }
    
    if ( a >= P )
    {
        a = P;
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        engaged = PLUS_END;
    }
    
    moveTo(a);
}


void Chewer::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
#if ( NEW_FIBER_END_CHEW & 1 )
    if ( engaged == PLUS_END )
    {
        modifiableFiber()->chew(engaged, prop()->chewing_speed_dt);
        moveToEnd(engaged);
        return;
    }
#endif
#if ( NEW_FIBER_END_CHEW & 2 )
    if ( engaged == MINUS_END )
    {
        modifiableFiber()->chew(engaged, prop()->chewing_speed_dt);
        moveToEnd(engaged);
        return;
    }
#endif

    // the load is the projection of the force on the local direction of Fiber
    real load = dot(force, dirFiber());
    real diff = prop()->line_diffusion_dt * RNG.sreal();
    real a = hAbs + diff + prop()->movability_dt * load;
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        a = M;
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        engaged = MINUS_END;
    }
    
    if ( a >= P )
    {
        a = P;
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        engaged = PLUS_END;
    }

    moveTo(a);
}

