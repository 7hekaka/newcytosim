// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "dim.h"
#include "exceptions.h"
#include "glossary.h"
#include "cutter_prop.h"
#include "cutter.h"
#include "simul_part.h"
#include "messages.h"


Hand * CutterProp::newHand(HandMonitor* m) const
{
    return new Cutter(this, m);
}


void CutterProp::clear()
{
    HandProp::clear();

    line_diffusion = 0;
    cutting_rate = 0;
    cutting_range = INFINITY;
    cut_width = 0;
    
    new_end_state[0] = STATE_WHITE;
    new_end_state[1] = STATE_WHITE;
}


void CutterProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(line_diffusion, "diffusion");
    glos.set(cutting_rate, "cutting_rate");
    glos.set(cutting_range, "cutting_range");
    glos.set(cut_width, "cut_width");

    // possible dynamic states of the ends
    Glossary::dict_type<state_t> keys{{"white",     STATE_WHITE},
                                      {"green",     STATE_GREEN},
                                      {"yellow",    STATE_YELLOW},
                                      {"orange",    STATE_ORANGE},
                                      {"red",       STATE_RED},
                                      {"static",    STATE_WHITE},
                                      {"grow",      STATE_GREEN},
                                      {"growing",   STATE_GREEN},
                                      {"shrink",    STATE_RED},
                                      {"shrinking", STATE_RED}};
    
    glos.set(new_end_state[0], "new_end_state", keys);
    glos.set(new_end_state[1], "new_end_state", 1, keys);
}


void CutterProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
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
    
    if ( line_diffusion > 0 && primed(sim) )
        Cytosim::log("   Cutter `", name(), "' has mobility = ", movability, " um/s\n");

    if ( cutting_rate < 0 )
        throw InvalidParameter(name()+"cutting_rate must be >= 0");

    if ( cutting_range < 0 )
        throw InvalidParameter(name()+"cutting_range must be >= 0");
    
    if ( cut_width < 0 )
        throw InvalidParameter(name()+"cut_width must be >= 0");

    cutting_rate_dt = cutting_rate * time_step(sim);
}


/**
 Estimate numerical stability from mobility and stiffness
 */
void CutterProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    HandProp::checkStiffness(stiff, len, mul, kT);
    
    real e = movability_dt * stiff * mul;
    if ( e > 2.0 )
    {
        std::ostringstream oss;
        oss << "Diffusible cutter `" << name() << "' is unstable since:\n";
        oss << PREF << "time_step * mobility * stiffness = " << e << '\n';
        oss << PREF << "-> reduce diffusion or time_step\n";
        throw InvalidParameter(oss.str());
        //std::clog << oss.str();
    }
    else
        Cytosim::log("   Slider `", name(), "' stability = ", e, "\n");
}


void CutterProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "diffusion", line_diffusion);
    write_value(os, "cutting_rate",  cutting_rate);
    write_value(os, "cutting_range", cutting_range);
    write_value(os, "cut_width", cut_width);
    write_value(os, "new_end_state", new_end_state, 2);
}

