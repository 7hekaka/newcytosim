// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "slider.h"
#include "slider_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"


Slider::Slider(SliderProp const* p, HandMonitor* h)
: Hand(p,h)
{
}


void Slider::stepUnloaded()
{
    assert_true( attached() );
    
    real a = hAbs + prop()->line_diffusion_dt * RNG.sreal();
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = P;
    }

    moveTo(a);
}


void Slider::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    real load = dot(force, dirFiber());
    real diff = prop()->line_diffusion_dt * RNG.sreal();
    real a = hAbs + diff + prop()->movability_dt * load;
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = P;
    }

    moveTo(a);
}

