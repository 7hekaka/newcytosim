// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "dim.h"
#include "cymdef.h"
#include "exceptions.h"
#include "glossary.h"
#include "nucleator_prop.h"
#include "nucleator.h"
#include "simul.h"


Hand * NucleatorProp::newHand(HandMonitor* m) const
{
    return new Nucleator(this, m);
}


void NucleatorProp::clear()
{
    HandProp::clear();

    nucleation_rate = 1;
    fiber_type = "";
    fiber_spec = "";
    fiber_class = nullptr;
    track_end  = NO_END;
    hold_end   = MINUS_END;
    addictive  = false;
    addictive_state = STATE_RED;
    branch_angle = 0;
    nucleation_limit = 0;
    nucleate_in_plane = "";
    nucleate_space = nullptr;
    branch_direction = BRANCH_PARALLEL;
}


void NucleatorProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(nucleation_rate, "rate");
    glos.set(fiber_type, "fibers");
    glos.set(fiber_spec, "fibers", 1);
    
    glos.set(nucleation_rate, "nucleate", 0);
    glos.set(fiber_type, "nucleate", 1);
    glos.set(fiber_spec, "nucleate", 2);
    
    glos.set(branch_angle, "branch_angle", "nucleation_angle");
    glos.set(nucleation_limit, "nucleation_limit");
    glos.set(nucleate_in_plane, "nucleate_in_plane");
    
#if BACKWARD_COMPATIBILITY < 100
    glos.set(fiber_spec, "nucleation_spec");
    glos.set(fiber_spec, "spec");
#endif
    
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
    
    glos.set(addictive, "addictive");
    glos.set(addictive_state, "addictive", 1, keys);
    
    if ( glos.set(track_end, "track_end", {{"off", NO_END},
        {"minus_end", MINUS_END}, {"plus_end", PLUS_END}}) )
        hold_end = track_end;
    
    glos.set(hold_end, "hold_end", {{"off", NO_END},
        {"minus_end", MINUS_END}, {"plus_end", PLUS_END}});
    
#if BACKWARD_COMPATIBILITY < 100
    glos.set(branch_direction, "specificity", {
        {"off", BRANCH_SPECIFIED},
        {"mostly_parallel", BRANCH_MOSTLY_PARALLEL} });
#endif
    
    glos.set(branch_direction, "branch_direction", {
        {"parallel", BRANCH_PARALLEL},
        {"mostly_parallel", BRANCH_MOSTLY_PARALLEL},
        {"antiparallel", BRANCH_ANTIPARALLEL},
        {"random", BRANCH_RANDOM},
        {"specified", BRANCH_SPECIFIED}} );
}


void NucleatorProp::complete(Simul const& sim)
{
    HandProp::complete(sim);

    if ( fiber_type.empty() )
        throw InvalidParameter(name()+":nucleate[1] (=fiber_type) must be specified if activity=nucleate");

    fiber_class = sim.findProperty<FiberProp>("fiber", fiber_type);
    
    if ( nucleation_rate < 0 )
        throw InvalidParameter(name()+":nucleate (=rate) must be positive");
    
    if ( nucleation_limit < 0 )
        throw InvalidParameter(name()+":nucleation_limit must be >= 0");

    if ( track_end && track_end != hold_end )
        throw InvalidParameter("if set, "+name()+":track_end should be equal to hold_end");
    
    nucleation_rate_dt = nucleation_rate * time_step(sim) * POOL_UNATTACHED;
    
#if BACKWARD_COMPATIBILITY < 100
    // extended meaning of 'nucleate_in_plane' on 27.11.2024
    if ( nucleate_in_plane == "1" )
    {
        nucleate_in_plane = sim.spaces.nameObject(fiber_class->confine_space);
        Cytosim::log(name(), ":nucleate_in_plane <---- ", nucleate_in_plane, "\n");
    }
#endif

    if ( ! nucleate_in_plane.empty() )
    {
        nucleate_space = sim.findSpace(nucleate_in_plane);
        if ( nucleate_space )
            nucleate_in_plane = sim.spaces.nameObject(nucleate_space);
        else if ( sim.primed() )
            throw InvalidParameter(name()+":nucleate_space `"+nucleate_in_plane+"' was not found");
    }
}



void NucleatorProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "nucleate",  nucleation_rate, fiber_type, "("+fiber_spec+")");
    write_value(os, "branch_angle", branch_angle);
    write_value(os, "branch_direction", branch_direction);
    write_value(os, "nucleation_limit", nucleation_limit);
    write_value(os, "nucleate_in_plane", nucleate_in_plane);

    write_value(os, "hold_end",  hold_end);
    write_value(os, "track_end", track_end);
    write_value(os, "addictive", addictive, addictive_state);
}

