// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "kinesin.h"
#include "kinesin_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "lattice.h"
#include "simul_part.h"
#include "hand_monitor.h"


Kinesin::Kinesin(KinesinProp const* p, HandMonitor* h)
: Digit(p,h)
{
}


void Kinesin::attach(FiberSite const& s)
{
    Digit::attach(s);
    nextAct = RNG.exponential();
    nextBack = RNG.exponential();

    // here digit::step_size must be equal to fiber:step_size
    if ( lattice() && lattice()->unit() != prop()->step_size  )
        throw InvalidParameter("kinesin:step_size must be equal to fiber:lattice_unit");
}


void Kinesin::stepUnloaded()
{
    assert_true( attached() );
    int stride = prop()->stepping_stride;

    nextAct -= prop()->forward_rate_dt * 0.5f;
    nextBack -= prop()->backward_rate_dt / 1.1f;

    while ( std::min(nextAct, nextBack) <= 0 )
    {
        // test detachment due to stepping
        if ( RNG.test(prop()->unbinding_chance) )
            return detach();

        int dir = ( nextAct <= nextBack ) - ( nextAct > nextBack );

        lati_t s = site() + dir * stride;
        int out = outsideMP(s);
        
        if ( out )
        {
            if ( RNG.test_not(prop()->hold_growing_end[out-1]) )
                return detach();
        }
        else if ( !valLattice(s) )
            hopLattice(s);
    
        if ( dir == 1 )
            nextAct += RNG.exponential();
        else
            nextBack += RNG.exponential();
    }
}


void Kinesin::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    int stride = prop()->stepping_stride;
    
    real load = dot(force, dirFiber()) * sign_real(stride);
    
    // antagonistic load is negative
    nextAct -= prop()->forward_rate_dt / ( 1 + std::exp(-load*prop()->force_inv) );
    nextBack -= prop()->backward_rate_dt / ( 0.1 + std::exp(load*prop()->force_inv) );

    while ( std::min(nextAct, nextBack) <= 0 )
    {
        // test detachment due to stepping
        if ( RNG.test(prop()->unbinding_chance) )
            return detach();

        int dir = ( nextAct <= nextBack ) - ( nextAct > nextBack );

        lati_t s = site() + dir * stride;
        int out = outsideMP(s);
        
        if ( out )
        {
            if ( RNG.test_not(prop()->hold_growing_end[out-1]) )
                return detach();
        }
        else if ( !valLattice(s) )
            hopLattice(s);
    
        if ( dir == 1 )
            nextAct += RNG.exponential();
        else
            nextBack += RNG.exponential();
    }
}

