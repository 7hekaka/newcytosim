// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "motor.h"
#include "motor_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "hand_monitor.h"
#include "messages.h"

Motor::Motor(MotorProp const* p, HandMonitor* h)
: Hand(p, h)
{
}


void Motor::stepUnloaded()
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


void Motor::stepLoaded(Vector const& force)
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

    //std::cerr << this << " > " << hAbs << "  " << dab << "\n";
    moveTo(a);
}

