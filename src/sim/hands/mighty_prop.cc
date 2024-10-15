// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "mighty_prop.h"
#include "mighty.h"
#include "simul_part.h"


Hand * MightyProp::newHand(HandMonitor* m) const
{
    return new Mighty(this, m);
}


void MightyProp::clear()
{
    HandProp::clear();
    
    stall_force       = 0;
    unloaded_speed    = 0;
    limit_speed       = true;
    var_speed_dt = 0;
    set_speed_dt = 0;
    rescue_chance = 0;
}


void MightyProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(stall_force,    "stall_force", "force");
    glos.set(unloaded_speed, "unloaded_speed", "speed");
#if BACKWARD_COMPATIBILITY < 100
    if ( glos.set(unloaded_speed, "max_speed") )
        Cytosim::warn("'max_speed' is deprecated: use 'unloaded_speed'\n");
#endif
    glos.set(limit_speed,    "limit_speed");
    glos.set(rescue_chance, "rescue_chance");
}


void MightyProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( primed(sim) && stall_force <= 0 )
        throw InvalidParameter("mighty:stall_force must be > 0");
    
    const real tau = time_step(sim);

    set_speed_dt = tau * unloaded_speed;
    var_speed_dt = abs_real(set_speed_dt) / stall_force;
    
    // The limits for a displacement in one time step apply if ( limit_speed = true )
    if ( unloaded_speed > 0 )
    {
        min_dab = 0;
        max_dab = 2 * tau * unloaded_speed;
    }
    else
    {
        min_dab = 2 * tau * unloaded_speed;
        max_dab = 0;
    }
    
    if ( rescue_chance < 0 )
        throw InvalidParameter("rescuer:rescue_chance must be >= 0");
}


void MightyProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    HandProp::checkStiffness(stiff, len, mul, kT);
    
    /*
     Compare mobility with stiffness: this can induce instability
     */
    real ef = abs_real(set_speed_dt) * stiff * mul / stall_force;
    if ( unloaded_speed != 0  &&  ef > 0.5 )
    {
        Cytosim::warn("simulating `", name(), "' may fail as:\n",
                      PREF, "time_step * stiffness * unloaded_speed / stall_force = ", ef, '\n',
                      PREF, "-> reduce time_step (really!)\n");
        //throw InvalidParameter(oss.str());
    }
    
    /*
     Compare the energy in a link due to the equipartition theorem
     to the maximum force that the motor can sustain before detaching:
     1/2 kT * DIM  <<  1/2 stiffness x^2 ~ 1/2 force^2 / stiffness;
     */
    if ( std::sqrt( DIM * kT * stiff ) > stall_force )
    {
        Cytosim::warn(name(),":stall_force is too small:\n",
                      PREF, "DIM * kT * stiffness > stall_force\n",
                      PREF, "-> reduce stiffness or increase stall_force\n");
    }
    
    /*
     Compare the force created by traveling during the time 1/unbinding_rate,
     and compare to stall_force. This is limit the efficiency of the motor.
     */
    ef = abs_real( stiff * unloaded_speed / ( unbinding_rate * stall_force ));
    if ( unbinding_rate != 0  &&  unloaded_speed != 0  &&  ef < 1 )
    {
        Cytosim::warn("Mighty `", name(), "' is inefficient since\n",
                      PREF, "stiffness * unloaded_speed / unbinding_rate << stall_force",
                      "; ratio = ", ef, "\n");
    }
    
    /*
     Compare detachment rate at stall-force, with detachment rate at rest
     */
    if ( std::exp( stall_force * unbinding_force_inv ) > 100 )
        Cytosim::warn(name(), ":exp( stall_force / unbinding_force ) is greater than 1!\n");
}


void MightyProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "stall_force",       stall_force);
    write_value(os, "unloaded_speed",    unloaded_speed);
    write_value(os, "limit_speed",       limit_speed);
    write_value(os, "rescue_chance", rescue_chance);
}

