// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "bundle_prop.h"
#include "glossary.h"
#include "simul.h"


void BundleProp::clear()
{
    stiffness  = -1;
    overlap    = -1;
    pole       = MINUS_END;
    bipolar    = true;
    fiber_rate = 0.0;
    fiber_type = "";
    fiber_spec = "";
}


void BundleProp::read(Glossary& glos)
{
    glos.set(stiffness, "stiffness");
    glos.set(overlap, "overlap");
    glos.set(pole, "focus", {{"plus_end", PLUS_END}, {"minus_end", MINUS_END}});
    glos.set(pole, "pole", {{"plus", PLUS_END}, {"minus", MINUS_END}});
    glos.set(bipolar, "bipolar");
    glos.set(fiber_rate, "nucleate");
    glos.set(fiber_type, "nucleate", 1);
    glos.set(fiber_spec, "nucleate", 2);
}


void BundleProp::complete(Simul const& sim)
{
    if ( fiber_rate < 0 )
        throw InvalidParameter("bundle:nucleation_rate (nucleate) must be >= 0");
    
    if ( fiber_rate > 0 )
    {
        if ( fiber_type.empty() )
            throw InvalidParameter("bundle:fibers (nucleate[1]) must be specified");
        
        // verify that fiber class exists:
        sim.properties.find_or_die("fiber", fiber_type);
    }

    fiber_prob = -std::expm1( -fiber_rate * time_step(sim) );

    if ( overlap < 0 )
        throw InvalidParameter("bundle:overlap must be specified and >= 0");
    
    if ( stiffness < 0 )
        throw InvalidParameter("bundle:stiffness must be specified and >= 0");
}


void BundleProp::write_values(std::ostream& os) const
{
    write_value(os, "stiffness", stiffness);
    write_value(os, "overlap",   overlap);
    write_value(os, "pole",      pole);
    write_value(os, "bipolar",   bipolar);
    write_value(os, "nucleate",  fiber_rate, fiber_type, "("+fiber_spec+")");
}

