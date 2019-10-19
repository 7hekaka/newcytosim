// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#include "tubule_prop.h"
#include "property_list.h"
#include "simul_prop.h"
#include "solid_prop.h"
#include "fiber_prop.h"
#include "glossary.h"


void TubuleProp::clear()
{
    stiffness[0] = 0;
    stiffness[1] = 0;
    fiber_type   = "";
}


void TubuleProp::read(Glossary& glos)
{
    glos.set(stiffness, 2, "stiffness");
    glos.set(fiber_type, "fiber");
}


void TubuleProp::complete(Simul const& sim)
{
    if ( stiffness[0] < 0 )
        throw InvalidParameter("tubule:stiffness[0] must be specified and >= 0");
    
    if ( stiffness[1] < 0 )
        throw InvalidParameter("tubule:stiffness[1] must be specified and >= 0");

    if ( fiber_type.empty() )
        throw InvalidParameter("tubule:fiber must be specified");

    sim.properties.find_or_die("fiber", fiber_type);
}


void TubuleProp::write_values(std::ostream& os) const
{
    write_value(os, "fiber",     fiber_type);
    write_value(os, "stiffness", stiffness[0], stiffness[1]);
}

