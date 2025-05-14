// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "aster_prop.h"
#include "solid_prop.h"
#include "fiber_prop.h"
#include "glossary.h"
#include "simul.h"


void AsterProp::clear()
{
    stiffness[0] = -1;
    stiffness[1] = -1;
    pole = MINUS_END;
    fiber_rate = 0;
    fiber_type = "";
    fiber_spec = "";
}


void AsterProp::read(Glossary& glos)
{
    glos.set(stiffness, 2, "stiffness");
    glos.set(pole, "focus", {{"plus_end", PLUS_END}, {"minus_end", MINUS_END}});
    glos.set(pole, "pole", {{"plus", PLUS_END}, {"minus", MINUS_END}});

    glos.set(fiber_rate, "nucleate", 0);
    glos.set(fiber_type, "nucleate", 1);
    glos.set(fiber_spec, "nucleate", 2);
    
#if BACKWARD_COMPATIBILITY < 100
    glos.set(fiber_rate, "nucleation_rate");
    glos.set(fiber_type, "fibers");
    glos.set(fiber_spec, "fibers", 1);
#endif
}


void AsterProp::complete(Simul const& sim)
{
    if ( stiffness[0] < 0 )
        throw InvalidParameter("aster:stiffness[0] must be specified and >= 0");
    
    if ( stiffness[1] < 0 )
        throw InvalidParameter("aster:stiffness[1] must be specified and >= 0");

    if ( fiber_rate < 0 )
        throw InvalidParameter("aster:nucleation rate (nucleate[0]) must be >= 0");
    
    if ( fiber_rate > 0 )
    {
        if ( fiber_type.empty() )
            throw InvalidParameter("aster:nucleation fiber type (nucleate[1]) must be specified");
 
        // verify that fiber class exists:
        sim.properties.find_or_die("fiber", fiber_type);
    }
 
    fiber_prob = -std::expm1( -fiber_rate * time_step(sim) );
}


void AsterProp::write_values(std::ostream& os) const
{
    write_value(os, "nucleate", fiber_rate, fiber_type, "("+fiber_spec+")");
    write_value(os, "stiffness", stiffness[0], stiffness[1]);
    write_value(os, "pole", pole);
}

