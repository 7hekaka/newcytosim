// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "mighty.h"
#include "mighty_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"


Mighty::Mighty(MightyProp const* p, HandMonitor* h)
: Hand(p, h), prop(p)
{
}


bool Mighty::attachmentAllowed(FiberSite& sit) const
{
    if ( !Hand::attachmentAllowed(sit) )
        return false;
    
    return true;
}


void Mighty::stepUnloaded()
{
    assert_true( attached() );
    
    real a = hAbs + prop->set_speed_dt;
    
    if ( a < hFiber->abscissaM() )
    {
        if ( RNG.test_not(prop->hold_growing_end) )
        {
            detach();
            return;
        }
        a = hFiber->abscissaM();
    }
    
    if ( a > hFiber->abscissaP() )
    {
        if ( RNG.test_not(prop->hold_growing_end) )
        {
            detach();
            return;
        }
        a = hFiber->abscissaP();
    }
    
    // detachment is also induced by displacement:
    assert_true( nextDetach >= 0 );
    nextDetach -= prop->unbinding_density * abs_real(a-hAbs);
    if ( nextDetach <= 0 )
        detach();
    else
        moveTo(a);
}


void Mighty::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // the load is the projection of the force on the local direction of Fiber
    real load = dot(force, dirFiber());
    
    // calculate load-dependent displacement:
    real dab = prop->set_speed_dt + load * prop->var_speed_dt;
    
    // possibly limit the range of the speed:
    if ( prop->limit_speed )
    {
        dab = std::max(dab, prop->min_dab);
        dab = std::min(dab, prop->max_dab);
    }
    
    real a = hAbs + dab;
    
    if ( a < hFiber->abscissaM() )
    {
        if ( RNG.test_not(prop->hold_growing_end) )
        {
            detach();
            return;
        }
        a = hFiber->abscissaM();
    }
    
    if ( a > hFiber->abscissaP() )
    {
        if ( RNG.test_not(prop->hold_growing_end) )
        {
            detach();
            return;
        }
        a = hFiber->abscissaP();
    }
    
    // detachment is also induced by displacement:
    assert_true( nextDetach >= 0 );
    nextDetach -= prop->unbinding_density * abs_real(a-hAbs);
    
    if ( nextDetach <= 0 )
        detach();
    else
        moveTo(a);
}

