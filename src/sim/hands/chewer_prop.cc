// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul_prop.h"
#include "simul_part.h"
#include "chewer_prop.h"
#include "chewer.h"


Hand * ChewerProp::newHand(HandMonitor* m) const
{
    return new Chewer(this, m);
}


void ChewerProp::clear()
{
    HandProp::clear();

    chewing_speed = 0;
    line_diffusion = 0;
}


void ChewerProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(chewing_speed, "chewing_rate");
    glos.set(chewing_speed, "chewing_speed");
    glos.set(line_diffusion, "diffusion");
}


void ChewerProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( chewing_speed < 0 )
        throw InvalidParameter(name()+"chewing_speed must be >= 0");

    chewing_speed_dt = chewing_speed * time_step(sim);
    
    if ( line_diffusion < 0 )
        throw InvalidParameter(name()+"diffusion must be >= 0");

    /*
     This is for unidimensional diffusion along the filaments, and we want:
     var(dx) = 2 D time_step, given that we use dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed in [-1, 1], its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D time_step`
     */
    line_diffusion_dt = std::sqrt(6.0 * line_diffusion * time_step(sim));
    
    // use Einstein's relation to get a mobility:
    real movability = line_diffusion / boltzmann(sim);
    movability_dt = movability * time_step(sim);
    
    if ( primed(sim) )
        Cytosim::log("   Chewer `", name(), "' has mobility = ", movability, " um/s\n");
}


/**
 Estimate numerical stability from mobility and stiffness
 */
void ChewerProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    HandProp::checkStiffness(stiff, len, mul, kT);
    
    real e = movability_dt * stiff * mul;
    if ( e > 2.0 )
    {
        std::ostringstream oss;
        oss << "Diffusible chewer `" << name() << "' is unstable since:\n";
        oss << PREF << "time_step * mobility * stiffness = " << e << '\n';
        oss << PREF << "-> reduce diffusion or time_step\n";
        throw InvalidParameter(oss.str());
        //std::clog << oss.str();
    }
    else
        Cytosim::log("   Slider `", name(), "' stability = ", e, "\n");
}


void ChewerProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "chewing_speed", chewing_speed);
    write_value(os, "diffusion", line_diffusion);
}

