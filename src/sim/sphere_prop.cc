// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University

#include "sphere_prop.h"
#include "glossary.h"
#include "messages.h"
#include "sphere.h"
#include "space.h"
#include "cymdef.h"

#include "sphere.h"
#include "simul_prop.h"
#include "space_prop.h"
#include "simul.h"


void SphereProp::clear()
{
    point_mobility = -1;
    viscosity     = -1;
    piston_effect = false;
    steric_key    = 0;
    steric_range  = 0;
    
    confine = CONFINE_OFF;
    confine_stiff = 0;
    confine_spec = "first";
    confine_space = nullptr;
    
    display       = "";
    display_fresh = false;
}


void SphereProp::read(Glossary& glos)
{
    glos.set(point_mobility,  "point_mobility");
    glos.set(piston_effect,   "piston_effect");
    glos.set(viscosity,       "viscosity");
    
    glos.set(steric_key,      "steric");
    glos.set(steric_range,    "steric", 1);
 
    glos.set(confine, "confine", {{"off",        CONFINE_OFF},
                                  {"on",         CONFINE_ON},
                                  {"inside",     CONFINE_INSIDE},
                                  {"none",       CONFINE_OFF},
                                  {"surface",    CONFINE_ON},
                                  {"all_inside", CONFINE_ALL_INSIDE}});
    
    glos.set(confine_stiff, "confine", 1);
    glos.set(confine_spec, "confine", 2);

    glos.set(confine_stiff, "confine_stiff");
    glos.set(confine_spec, "confine_spec");

#if BACKWARD_COMPATIBILITY < 50
    if ( confine_spec == "current" )
        confine_spec = "last";

    glos.set(confine, "confined",{{"none",    CONFINE_OFF},
                                  {"inside",  CONFINE_INSIDE},
                                  {"surface", CONFINE_ON}});
    glos.set(confine_stiff, "confined", 1);
#endif
    
    if ( glos.set(display, "display") )
        display_fresh = true;
}


void SphereProp::complete(Simul const& sim)
{
    if ( viscosity < 0 )
        viscosity = sim.prop.viscosity;
        
    if ( viscosity <= 0 )
        throw InvalidParameter("sphere:viscosity or simul:viscosity should be defined > 0");
    
    confine_space = sim.findSpace(confine_spec);
    if ( confine != CONFINE_OFF )
    {
        if ( confine_space )
        {
            if ( confine_spec.empty() )
                confine_spec = sim.spaces.nameObject(confine_space);
        }
        else
        {
            // this condition may occur when the Property is created before the Space
            if ( primed(sim) )
                throw InvalidParameter(name()+":confine_space `"+confine_spec+"' was not found");
        }
    }

    if ( confine && confine_stiff < 0 )
        throw InvalidParameter(name()+":confine_stiff must be >= 0");
    
    if ( primed(sim) && steric_key && !sim.prop.steric_mode )
        Cytosim::warn(name(), ":steric is set but simul:steric = 0\n");

    if ( point_mobility < 0 )
        throw InvalidParameter("sphere:point_mobility must be specified and >= 0");
}


void SphereProp::write_values(std::ostream& os) const
{
    write_value(os, "viscosity",      viscosity);
    write_value(os, "point_mobility", point_mobility);
    write_value(os, "piston_effect",  piston_effect);
    write_value(os, "steric",         steric_key, steric_range);
    write_value(os, "confine",        confine, confine_stiff, confine_spec);
    write_value(os, "display",        "("+display+")");
}

