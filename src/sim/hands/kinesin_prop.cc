// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "kinesin.h"
#include "kinesin_prop.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul_part.h"


Hand * KinesinProp::newHand(HandMonitor* m) const
{
    return new Kinesin(this, m);
}


void KinesinProp::clear()
{
    DigitProp::clear();

    /** values here are based on Rob Cross measurements of conventional kinesin:
    
        Mechanics of the kinesin step, Carter & Cross, 2005
        http://www.doi.org/10.1038/nature03528
     */
    stepping_force   = 2;
    forward_rate     = 277;
    backward_rate    = 0.34;
    unbinding_chance = 0.01;
    stepping_stride  = 1;
}


void KinesinProp::read(Glossary& glos)
{
    DigitProp::read(glos);
    
    glos.set(stepping_force,   "stepping_force", "force");
    glos.set(forward_rate,     "forward_rate");
    glos.set(backward_rate,    "backward_rate");
    glos.set(unbinding_chance, "unbinding_chance");
    glos.set(stepping_stride,  "stepping_stride", "stride");
}


void KinesinProp::complete(Simul const& sim)
{
    DigitProp::complete(sim);
   
    if ( primed(sim) && stepping_force <= 0 )
        throw InvalidParameter("kinesin:force must be > 0");
    force_inv = 1.f / stepping_force;
    
    if ( forward_rate < 0 )
        throw InvalidParameter("kinesin:forward_rate must be >= 0");
    forward_rate_dt = forward_rate * time_step(sim);
    
    if ( backward_rate < 0 )
        throw InvalidParameter("kinesin:backward_rate must be >= 0");
    backward_rate_dt = backward_rate * time_step(sim);
    
    if ( primed(sim) )
    {
        real stall = log(forward_rate/backward_rate)*0.5/force_inv;
        printf("Kinesin's stall force is %.4f pN\n", stall);
    }
}


void KinesinProp::write_values(std::ostream& os) const
{
    DigitProp::write_values(os);
    write_value(os, "stepping_force", stepping_force);
    write_value(os, "forward_rate", forward_rate);
    write_value(os, "backward_rate", backward_rate);
    write_value(os, "unbinding_chance", unbinding_chance);
    write_value(os, "stepping_stride", stepping_stride);
}

