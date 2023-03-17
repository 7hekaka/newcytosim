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
    
    if ( engaged != NO_END )
    {
#if NEW_FIBER_END_CHEW
        hFiber->chew(engaged, prop()->chewing_speed_dt);
        moveToEnd(engaged);
#else
        throw InvalidParameter("fiber:chew is not enabled");
#endif
        return;
    }

    real a = hAbs + prop()->diffusion_dt * RNG.sreal();
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        a = M;
        engaged = MINUS_END;
    }
    
    if ( a >= P )
    {
        a = P;
        engaged = PLUS_END;
    }
    
    if ( engaged && RNG.test_not(prop()->hold_growing_end) )
        return detach();

    moveTo(a);
}


void Chewer::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    if ( engaged != NO_END )
    {
#if NEW_FIBER_END_CHEW
        hFiber->chew(engaged, prop()->chewing_speed_dt);
        moveToEnd(engaged);
#else
        throw InvalidParameter("fiber:chew is not enabled");
#endif
        return;
    }
    
    // the load is the projection of the force on the local direction of Fiber
    real load = dot(force, dirFiber());
    real a = hAbs + prop()->diffusion_dt * RNG.sreal() + prop()->mobility_dt * load;
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        a = M;
        engaged = MINUS_END;
    }
    
    if ( a >= P )
    {
        a = P;
        engaged = PLUS_END;
    }
    
    if ( engaged && RNG.test_not(prop()->hold_growing_end) )
        return detach();

    moveTo(a);
}

