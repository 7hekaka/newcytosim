// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "slider.h"
#include "slider_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"


Slider::Slider(SliderProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
}


void Slider::stepUnloaded()
{
    assert_true( attached() );
    
    /// diffusion?
}


void Slider::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    real a = hAbs + dot(force, dirFiber()) * prop->mobility_dt;
    
    if ( a < hFiber->abscissaM() )
    {
        if ( RNG.test_not(prop->hold_growing_end) )
            return detach();
        a = hFiber->abscissaM();
    }
    
    if ( a > hFiber->abscissaP() )
    {
        if ( RNG.test_not(prop->hold_growing_end) )
            return detach();
        a = hFiber->abscissaP();
    }

    moveTo(a);
}

