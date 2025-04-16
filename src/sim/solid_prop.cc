// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "solid_prop.h"
#include "glossary.h"
#include "messages.h"
#include "simul_prop.h"
#include "simul.h"
#include "cymdef.h"


void SolidProp::clear()
{
    drag         = -1;
    viscosity    = -1;
    steric_key   = 0;
    steric_range = 0;
    
    confine = CONFINE_OFF;
    confine_stiff = 0;
    confine_spec = "first";
    confine_space = nullptr;
    
#if NEW_RADIAL_FLOW
    flow_time[0] = 0;
    flow_time[1] = 0;
    flow_center.reset();
#endif
    display = "";
    display_fresh = false;
}


void SolidProp::read(Glossary& glos)
{
    glos.set(drag,           "drag");
    glos.set(viscosity,      "viscosity");
    
    glos.set(steric_key,     "steric");
    glos.set(steric_range,   "steric", 1);
    
    glos.set(confine,        "confine", {{"off",          CONFINE_OFF},
                                         {"on",           CONFINE_ON},
                                         {"inside",       CONFINE_INSIDE},
                                         {"outside",      CONFINE_OUTSIDE},
                                         {"point",        CONFINE_POINT},
                                         {"point_inside", CONFINE_POINT_INSIDE},
#if BACKWARD_COMPATIBILITY < 50
                                         {"none",         CONFINE_OFF},
                                         {"surface",      CONFINE_ON},
#endif
                                         {"all_inside",   CONFINE_ALL_INSIDE}});
    
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
    
#if NEW_RADIAL_FLOW
    glos.set(flow_time,  2, "flow_time");
    glos.set(flow_center,   "flow_center");
    if ( flow_time[0] > flow_time[1] )
        throw InvalidParameter("flow_time[0] should be lower than flow_time[1]");
#endif

    if ( glos.set(display, "display") )
        display_fresh = true;
}


void SolidProp::complete(Simul const& sim)
{
    if ( viscosity < 0 )
        viscosity = sim.prop.viscosity;
    
    if ( viscosity <= 0 )
        throw InvalidParameter("bead:viscosity or simul:viscosity should be defined > 0");
    
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
    
}


void SolidProp::write_values(std::ostream& os) const
{
    if ( drag > 0 )  write_value(os, "drag", drag);
    write_value(os, "viscosity", viscosity);
    write_value(os, "steric",    steric_key, steric_range);
    write_value(os, "confine",   confine, confine_stiff, confine_spec);
#if NEW_RADIAL_FLOW
    write_value(os, "flow_center", flow_center);
    write_value(os, "flow_time",   flow_time, 2);
#endif
    write_value(os, "display",   "("+display+")");
}

