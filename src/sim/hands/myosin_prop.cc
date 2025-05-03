// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "myosin.h"
#include "myosin_prop.h"
#include "exceptions.h"
#include "glossary.h"
#include "messages.h"
#include "simul_part.h"


Hand * MyosinProp::newHand(HandMonitor* m) const
{
    return new Myosin(this, m);
}


void MyosinProp::clear()
{
    DigitProp::clear();

    stall_force     = 0;
    unloaded_speed  = 0;
    walking_rate_dt = 0;
    var_rate_dt     = 0;
}


void MyosinProp::read(Glossary& glos)
{
    DigitProp::read(glos);
    
    glos.set(stall_force,    "stall_force", "force");
    glos.set(unloaded_speed, "unloaded_speed", "speed");
#if BACKWARD_COMPATIBILITY < 100
    if ( glos.set(unloaded_speed, "max_speed") )
        Cytosim::warn("'max_speed' is deprecated: use 'unloaded_speed'\n");
#endif
}


void MyosinProp::complete(Simul const& sim)
{
    DigitProp::complete(sim);
   
    if ( primed(sim) && stall_force <= 0 )
        throw InvalidParameter("myosin:stall_force must be > 0");
    
    if ( unloaded_speed < 0 )
        throw InvalidParameter("myosin:unloaded_speed must be >= 0");

    walking_rate_dt = time_step(sim) * abs_real(unloaded_speed) / step_size;
    var_rate_dt = std::copysign(walking_rate_dt/stall_force, unloaded_speed);
}


void MyosinProp::write_values(std::ostream& os) const
{
    DigitProp::write_values(os);
    write_value(os, "stall_force",    stall_force);
    write_value(os, "unloaded_speed", unloaded_speed);
}

