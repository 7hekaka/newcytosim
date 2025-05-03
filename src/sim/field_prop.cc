// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "field_prop.h"
#include "messages.h"
#include "glossary.h"
#include "simul.h"
#include "dim.h"

//------------------------------------------------------------------------------
void FieldProp::clear()
{
    step                  = 0;
    field_space           = "first";
    field_space_ptr       = nullptr;
    field_periodic        = 0;
    slow_diffusion        = 0;
    full_diffusion        = 0;
    boundary_condition    = 0;
    boundary_value        = 0;
    decay_rate            = 0;
    transport_strength    = 0;
    transport_length      = 0;
    cut_fibers            = 0;
    chew_fibers           = 0;
    positive              = 0;
    save                  = true;
    display_scale         = 1;
    visible               = 1;
}


//------------------------------------------------------------------------------
void FieldProp::read(Glossary& glos)
{
    Glossary::dict_type<int> keys{{"flux", 0}, {"edge", 7},
        {"edgeX", 1}, {"edgeY", 2}, {"edgeZ", 4},
        {"edgeXY", 3}, {"edgeYZ", 6}, {"edgeXZ", 5}};
    
    glos.set(step,               "step");
    glos.set(field_periodic,     "periodic");
    glos.set(field_space,        "space");
    glos.set(slow_diffusion,     "diffusion");
    glos.set(full_diffusion,     "full_diffusion");
    glos.set(boundary_condition, "boundary_condition", keys);
    glos.set(boundary_value,     "boundary_value");
    glos.set(boundary_condition, "boundary", keys);
    glos.set(boundary_value,     "boundary", 1);
    glos.set(decay_rate,         "decay_rate");
    
    if ( glos.has_key("bind_fibers") )
        Cytosim::warn("field::bind_fibers is obsolete, use fiber:mesh_binding_rate\n");

    glos.set(transport_strength, "transport", 0);
    glos.set(transport_length,   "transport", 1);
    glos.set(transport_strength, "transport_strength");
    glos.set(transport_length,   "transport_length");
    
    glos.set(cut_fibers,         "cut_fibers");
    glos.set(chew_fibers,        "chew_fibers");

    glos.set(positive,           "positive");
    glos.set(display_scale,      "display_scale");
    glos.set(display_scale,      "scale_max"); // BACKWARD_COMPATIBILITY
    glos.set(visible,            "visible");
    glos.set(save,               "save");
}


//------------------------------------------------------------------------------
void FieldProp::complete(Simul const& sim)
{
    field_space_ptr = sim.findSpace(field_space);
    
    if ( field_space_ptr )
        field_space = sim.spaces.nameObject(field_space_ptr);
    else if ( primed(sim) )
        throw InvalidParameter("field::space must be created before the field");

    if ( step < REAL_EPSILON )
        throw InvalidParameter("field:step must be defined and > 0");
    
    if ( slow_diffusion < 0 )
        throw InvalidParameter("field:diffusion must be >= 0");
    
    real CFL = 2 * DIM * time_step(sim) * (slow_diffusion+full_diffusion) / ( step * step );
    //std::clog << "The diffusion CFL condition for `" << name() << "' is " << CFL << '\n';
    
    if ( primed(sim)  &&  CFL > 0.5 )
    {
        InvalidParameter e("field:diffusion is too fast");
        e << "The CFL condition ( diffusion * time_step / step^2 ) must be below 1/2,";
        e << "\n  but it is " + std::to_string(CFL) + "\n";
        throw e;
    }

    if ( decay_rate < 0 )
        throw InvalidParameter("field:decay_rate must be >= 0");
}


//------------------------------------------------------------------------------

void FieldProp::write_values(std::ostream& os) const
{
    write_value(os, "step",           step);
    write_value(os, "space",          field_space);
    write_value(os, "periodic",       field_periodic);
    write_value(os, "diffusion",      slow_diffusion);
    write_value(os, "full_diffusion", full_diffusion);
    write_value(os, "boundary",       boundary_condition, boundary_value);
    write_value(os, "decay_rate",     decay_rate);
    write_value(os, "transport",      transport_strength, transport_length);
    write_value(os, "cut_fibers",     cut_fibers);
    write_value(os, "chew_fibers",    chew_fibers);
    write_value(os, "positive",       positive);
    write_value(os, "display_scale",  display_scale);
    write_value(os, "visible",        visible);
    write_value(os, "save",           save);
}

