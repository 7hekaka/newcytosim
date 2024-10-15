// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "fiber_prop.h"
#include "cymdef.h"
#include <cmath>
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "property_list.h"
#include "single_prop.h"
#include "simul_prop.h"
#include "simul.h"
#include "fiber.h"

/**
 This is virtualized to return a derived Fiber from classes derived from FiberProp
 */
Fiber* FiberProp::newFiber() const
{
    return new Fiber(this);
}


/**
 Return the length specified by the user in 'opt'
 */
real FiberProp::newFiberLength(Glossary& opt) const
{
    real len = 1.0;
    
    if ( 0 < opt.num_values("fiber_length") )
    {
        opt.set_from_least_used_value(len, "fiber_length");
    }
    else if ( opt.set(len, "length") )
    {
        real var = 0.0;
        if ( opt.value_is("length", 1, "exponential") )
        {
            real L = len;
            // exponential distribution, truncated (Julio M.Belmonte's student):
            do
                len = L * RNG.exponential();
            while (( len < min_length ) | ( len > max_length ));
        }
        else if ( opt.set(var, "length", 1) )
        {
            real L = len;
            // add variability without changing mean:
            do
                len = L + var * RNG.sreal();
            while (( len < min_length ) | ( len > max_length ));
        }
        else
        {
            if ( len < min_length )
                std::cerr << "Warning: "<<name()<<"'s initial length < min_length\n";
            if ( len > max_length )
                std::cerr << "Warning: "<<name()<<"'s initial length > max_length\n";
        }
    }
    else
    {
        len = std::max(len, min_length);
        len = std::min(len, max_length);
    }
    return len;
}


/**
 @addtogroup FiberGroup
 @{
 <hr>
 
 When creating a new Fiber, you may specify:
 - the initial length,
 - the initial state of the plus and minus ends,
 - if the position refers to the center or to the tip of the fiber
 - the shape, using a set of points
 .
 
 Syntax:
 
     new filament
     {
       length = REAL, LENGTH_MODIFIER
       end_state = PLUS_END_STATE, MINUS_END_STATE
       reference = REFERENCE
     }
 
 The optional LENGTH_MODIFIER can be:
 - `exponential`,
 - REAL
 .
 This introduces variability, without changing the mean length.
 The second form generates a flat distribution of width 2*LENGTH_MODIFIER.
 
 The initial states PLUS_END_STATE and MINUS_END_STATE can be:
 - 0 = `white` or `static`
 - 1 = `green` or `grow`
 - 4 = `red`   or `shrink`
 .
 
 Optional reference specificiation:
 - center [default]
 - plus_end
 - minus_end
 .
 
 To use a persistent random walk as initial shape, set:
 
     new filament
     {
       equilibrate = PERSISTENCE_LENGTH
       position = POSITION
       direction = DIRECTION
     }
 
 In this case the shape will be random and different for each filament, and
 characterized by the given persistence length (units of length). The shape will
 be translated to bring its center of gravity at the specified position, and
 rotated to match the direction specified with the average filament direction.
 As usual, if 'position' and 'direction' are not specified, they are random.

 To specify the shape of a Fiber directly, use:
 
     new filament
     {
         shape = POSITION, POSITION, ...
     }
 
 Examples:
 
     new filament
     {
       length = 1
       plus_end = grow
       minus_end = static
     }

 which is equivalent to:

     new filament
     {
       length = 1
       end_state = green, white
     }

     new filament
     {
       position = 0 0 0
       orientation = 1 0 0
       shape = -4 -3 0, -3 0 0, -1 2 0, 1  3 0
     }
 
 @}
 */
Fiber* FiberProp::newFiber(Glossary& opt) const
{
    Fiber * fib = newFiber();
    const real len = newFiberLength(opt);

#if ( 1 )
    // specify the vertices directly:
    if ( opt.has_key("points") )
    {
        index_t nbp = opt.num_values("points");
        fib->setNbPoints(nbp);
        for ( index_t p = 0; p < nbp; ++p )
        {
            Vector vec(0,0,0);
            if ( ! opt.set(vec, "points", p) )
                throw InvalidParameter("fiber:points must be a list of comma-separated vectors");
            fib->setPoint(p, vec);
        }
        if ( opt.has_key("length") )
            fib->imposeLength(len);
    }
    else
#endif
    if ( opt.has_key("shape") )
    {
        if ( opt.value_is("shape", 0, "curved") )
        {
            real rad = 1, off = 0;
            Vector dir(0, 1, 0);
            opt.set(rad, "shape", 1);
            opt.set(dir, "shape", 2);
            opt.set(off, "shape", 3);
            fib->setCurved(dir, rad, len, off);
        }
        else
        {
            size_t nbp = opt.num_values("shape");
            if ( nbp < 2 )
                throw InvalidParameter("fiber:shape must be a list of comma-separated points");
            
            real* tmp = new_real(DIM*nbp);
            for ( size_t p = 0; p < nbp; ++p )
            {
                Vector vec(0,0,0);
                if ( ! opt.set(vec, "shape", p) )
                    throw InvalidParameter("fiber:shape must be a list of comma-separated points");
                vec.store(tmp+DIM*p);
            }
            fib->setShape(tmp, nbp, 0);
            if ( fib->nbPoints() < 2 )
                throw InvalidParameter("the vectors specified in fiber:shape must not overlap");
            free_real(tmp);
        }
    }
    else
    {
        real pl = 0; // persistence length
        // place fiber horizontally with center at the origin:
        if ( opt.set(pl, "equilibrate") && pl > 0 )
            fib->setEquilibrated(len, pl);
        else
            fib->setStraight(Vector(-0.5*len,0,0), Vector(1,0,0), len);
        
        FiberEnd ref = CENTER;
        if ( opt.set(ref, "reference", {{"plus_end", PLUS_END}, {"minus_end", MINUS_END}, {"center", CENTER}}) )
            fib->placeEnd(ref);
    }
    
    // set abscissa of minus-end
    real a = 0;
    if ( opt.set(a, "origin") )
        fib->setOrigin(a);
    
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
    
    // set state of plus ends:
    state_t p = STATE_WHITE;
#if BACKWARD_COMPATIBILITY < 50
    if ( opt.set(p, "plus_end_state") )
    {
        fib->setEndStateP(p);
        Cytosim::warn("use `plus_end = STATE` instead of `plus_end_state = STATE`\n");
    }
#endif
    if ( opt.set(p, "plus_end", keys) || opt.set(p, "end_state", keys) )
        fib->setEndStateP(p);

    // set state of minus ends:
    state_t m = STATE_WHITE;
#if BACKWARD_COMPATIBILITY < 50
    if ( opt.set(m, "minus_end_state") )
    {
        Cytosim::warn("use `minus_end = STATE` instead of `minus_end_state = STATE`\n");
        fib->setEndStateM(m);
    }
#endif
    if ( opt.set(m, "minus_end", keys) || opt.set(m, "end_state", 1, keys) )
        fib->setEndStateM(m);

#if BACKWARD_COMPATIBILITY < 100
    if ( fib->prop->activity != "none" && m == 0 && p == 0 )
        Cytosim::warn("`", fib->prop->name(), "' may not grow since both ends are in state `white`\n");
#endif
    
    fib->updateFiber();
    
#if FIBER_HAS_DENSITY
    if ( fib->densityField().ready() )
    {
        // enable density initialization
        real val = 0;
        if ( opt.set(val, "mesh_value") )
            fib->densityField().clear(val);
    }
#endif

    return fib;
}


//------------------------------------------------------------------------------
void FiberProp::clear()
{
    rigidity      = -1;
    segmentation  = 1;
    min_length    = 0.008;   // suitable for actin/microtubules
    max_length    = INFINITY;
    total_polymer = INFINITY;
    persistent    = false;

    viscosity    = -1;
    drag_radius  = 0.0125;  // radius of a Microtubule
    drag_length  = 5;
    drag_model   = 0;
    drag_gap     = 0;
    
    binding_key  = ~0U;  //all bits at 1

    lattice      = 0;
    lattice_unit = 0;
    save_lattice = 0;
#if FIBER_HAS_DENSITY
    density                = 0;
    density_unit           = 0;
    density_cut_fiber      = 0;
    density_flux_speed     = 0;
    density_binding_rate   = 0;
    density_unbinding_rate = 0;
    density_aging_rate     = 0;
    density_aging_limit    = 1;
#endif
    confine = CONFINE_OFF;
    confine_stiff[0] = 0;
    confine_stiff[1] = 0;
    confine_spec = "first";
    confine_space = nullptr;
    
#if NEW_FIBER_CONFINE2
    confine2 = CONFINE_OFF;
    confine2_stiff[0] = 0;
    confine2_stiff[1] = 0;
    confine2_spec = "first";
    confine2_space = nullptr;
#endif

    steric_key = 0;
    steric_radius = 0;
    steric_range = 0;
    
    field     = "none";
    field_ptr = nullptr;
    
    glue = 0;
    glue_single = "none";
    glue_prop = nullptr;
    
#if NEW_COLINEAR_FORCE
    colinear_force = 0;
#endif
#if NEW_FIBER_END_CHEW
    max_chewing_speed = 1.0;
#endif
#if NEW_FIBER_LOOP
    loop = 0;
#endif
    activity      = "none";
    display       = "";
    display_fresh = false;
    
    used_polymer = 0;
    free_polymer = 1;
    fiber_count = 0;
    
#if NEW_SQUEEZE_FORCE
    squeeze_mode  = 0;
    squeeze_force = 0;
    squeeze_range = 1;
#endif
#if NEW_FIBER_END_FORCE
    end_force_mode = NO_END;
    end_force.set(0,0,0);
#endif

}

//------------------------------------------------------------------------------
void FiberProp::read(Glossary& glos)
{
    glos.set(rigidity,      "rigidity");
    glos.set(segmentation,  "segmentation");
    glos.set(min_length,    "min_length");
    glos.set(max_length,    "max_length");
    glos.set(total_polymer, "total_polymer");
    glos.set(persistent,    "persistent");
#if BACKWARD_COMPATIBILITY < 50
    bool ds;
    if ( glos.set(ds, "delete_stub") )
    {
        persistent = !ds;
        Cytosim::warn("use `persistent=", !ds, "' instead of `delete_stub="<<ds<<"'\n");
    }
#endif
    
    glos.set(viscosity,    "viscosity");
    glos.set(drag_radius,  "drag_radius");
    glos.set(drag_length,  "drag_length");
    glos.set(drag_model,   "drag_model");
    glos.set(drag_gap,     "drag_model", 1);
#if BACKWARD_COMPATIBILITY < 50
    glos.set(drag_radius,  "hydrodynamic_radius");
    glos.set(drag_length,  "hydrodynamic_radius", 1);
    glos.set(drag_model,   "surface_effect");
    glos.set(drag_gap,     "surface_effect", 1);
#endif

    glos.set(binding_key,  "binding_key");
    
    glos.set(lattice,      "lattice");
    glos.set(lattice_unit, "lattice", 1, "lattice_unit", 1);
    glos.set(save_lattice, "save_lattice");
    
#if FIBER_HAS_DENSITY
    if ( glos.set(density, "density") )
    {
        glos.set(density_unit, "density", 1, "density_unit", 0);
        glos.set(density_cut_fiber, "density_cut_fiber");
        glos.set(density_flux_speed, "density_flux_speed");
        glos.set(density_binding_rate, "density_binding_rate");
        glos.set(density_unbinding_rate, "density_unbinding_rate");
        glos.set(density_aging_rate, "density_aging_rate");
        glos.set(density_aging_rate, "density_aging");
        glos.set(density_aging_limit, "density_aging", 1);
        glos.set(density_aging_limit, "density_aging_limit");
    }
#  ifdef BACKWARD_COMPATIBILITY
    // parameters `lattice_*' were renamed `mesh_*' on 28.11.2019 and `density_*' on 25.05.2025
    else if ( glos.set(density, "mesh") )
    {
        glos.set(density_unit, "mesh", 1, "mesh_unit", 0);
        glos.set(density_cut_fiber, "mesh_cut_fiber");
        glos.set(density_flux_speed, "mesh_flux_speed");
        glos.set(density_binding_rate, "mesh_binding_rate");
        glos.set(density_unbinding_rate, "mesh_unbinding_rate");
        glos.set(density_aging_rate, "mesh_aging_rate");
    }
#  endif
#endif
    
    glos.set(confine, "confine", {{"off",       CONFINE_OFF},
                                  {"on",        CONFINE_ON},
                                  {"inside",    CONFINE_INSIDE},
                                  {"outside",   CONFINE_OUTSIDE},
                                  {"in_out",    CONFINE_IN_OUT},
                                  {"none",      CONFINE_OFF},
                                  {"surface",   CONFINE_ON},
                                  {"plus_end",  CONFINE_PLUS_END},
                                  {"minus_end", CONFINE_MINUS_END},
                                  {"both_ends", CONFINE_BOTH_ENDS},
                                  {"minus_out", CONFINE_MINUS_OUT},
                                  {"plus_out",  CONFINE_PLUS_OUT}});
    
    if ( glos.set(confine_stiff[0], "confine", 1) )
        confine_stiff[1] = confine_stiff[0];
    if ( glos.is_number("confine", 2) )
    {
        glos.set(confine_stiff[1], "confine", 2);
        glos.set(confine_spec, "confine", 3);
    }
    else
        glos.set(confine_spec, "confine", 2);
    glos.set(confine_stiff, 2, "confine_stiff");
    glos.set(confine_spec, "confine_spec");

#if NEW_FIBER_CONFINE2
    glos.set(confine2, "confine2", {{"off",    CONFINE_OFF},
                                   {"on",      CONFINE_ON},
                                   {"inside",  CONFINE_INSIDE},
                                   {"outside", CONFINE_OUTSIDE},
                                   {"in_out",  CONFINE_IN_OUT},}, 1);
    
    if ( glos.set(confine2_stiff[0], "confine2", 1) )
        confine2_stiff[1] = confine2_stiff[0];
    if ( glos.is_number("confine2", 2) )
    {
        glos.set(confine2_stiff[1], "confine2", 2);
        glos.set(confine2_spec, "confine2", 3);
    }
    else
        glos.set(confine2_spec, "confine2", 2);

    glos.set(confine2_stiff, 2, "confine2_stiffness");
    glos.set(confine2_spec, "confine2_space");
#endif

    
#if BACKWARD_COMPATIBILITY < 50
    if ( confine_spec == "current" )
        confine_spec = "last";

    glos.set(confine, "confined", {{"none",     CONFINE_OFF},
                                  {"inside",    CONFINE_INSIDE},
                                  {"outside",   CONFINE_OUTSIDE},
                                  {"surface",   CONFINE_ON},
                                  {"minus_end", CONFINE_MINUS_END},
                                  {"plus_end",  CONFINE_PLUS_END}});
    
    glos.set(confine_stiff[0], "confined", 1);
#endif

#if NEW_SQUEEZE_FORCE
    glos.set(squeeze_mode,  "squeeze", {{"off", 0}, {"all", 1}, {"minus_end", 2}});
    glos.set(squeeze_force, "squeeze", 1);
    glos.set(squeeze_range, "squeeze", 2);
#endif
#if NEW_FIBER_END_FORCE
    glos.set(end_force, "end_force", 1);
    glos.set(end_force_mode, "end_force", {{"off", NO_END}, {"plus_end", PLUS_END},
        {"minus_end", MINUS_END}, {"center", CENTER}, {"both_ends", BOTH_ENDS}, {"torque", ORIGIN}});
#endif

    glos.set(steric_key,    "steric");
    glos.set(steric_radius, "steric", 1);
    glos.set(steric_range,  "steric", 2);
    glos.set(steric_radius, "steric_radius");
    glos.set(steric_range,  "steric_range");
    glos.set(field, "field");
#if FIBER_HAS_GLUE
    glos.set(glue,        "glue");
    glos.set(glue_single, "glue", 1);
#endif
#if NEW_COLINEAR_FORCE
    glos.set(colinear_force, "colinear_force");
#endif
#if NEW_FIBER_END_CHEW
    glos.set(max_chewing_speed, "max_chewing_speed");
#endif
#if NEW_FIBER_LOOP
    glos.set(loop, "loop");
#endif
    glos.set(activity, "activity");
    if ( glos.set(display, "display") )
        display_fresh = true;
}


void FiberProp::complete(Simul const& sim)
{
    if ( viscosity < 0 )
        viscosity = sim.prop.viscosity;
    
    if ( viscosity <= 0 )
        throw InvalidParameter("fiber:viscosity or simul:viscosity should be defined > 0");

    /* confine_spec is also used for `glue' and `shrink_outside',
     and needs to be set here */
    confine_space = sim.findSpace(confine_spec);
    if ( confine_space )
    {
        if ( confine_spec.empty() )
        {
            confine_spec = sim.spaces.nameObject(confine_space);
            //std::cerr << ":confine_spec <-- " << confine_spec << "\n";
        }
    }
    else
    {
        // this condition may occur when the Property is created before the Space
        if ( confine != CONFINE_OFF && primed(sim) )
            throw InvalidParameter(name()+":confine_space `"+confine_spec+"' was not found");
    }
    
    if ( confine && confine_stiff[0] < 0 )
        throw InvalidParameter(name()+":confine_stiff must be >= 0");

#if NEW_FIBER_CONFINE2
    if ( confine2 != CONFINE_OFF )
    {
        confine2_space = sim.findSpace(confine2_spec);
        if ( confine2_space )
        {
            if ( confine2_spec.empty() )
                confine2_spec = sim.spaces.nameObject(confine2_space);
        }
        else
        {
            if ( primed(sim) )
                throw InvalidParameter(name()+":confine2_space `"+confine2_spec+"' was not found");
            // this condition occur when the Property is created before the Space
        }
    }
    else
        confine2_space = nullptr;

    if ( confine2 && confine2_stiff[0] < 0 )
        throw InvalidParameter(name()+":confine2_stiffness must be specified and >= 0");
#endif
    
    if ( primed(sim) && steric_key && !sim.prop.steric_mode )
        Cytosim::warn(name(), ":steric is set but simul:steric = 0\n");

    if ( min_length < 0 )
        throw InvalidParameter("fiber:min_length should be >= 0");

    if ( max_length < 0 )
        throw InvalidParameter("fiber:max_length should be >= 0");

#if BACKWARD_COMPATIBILITY < 100
    if ( total_polymer == 0 )
        total_polymer = INFINITY;
#endif
    if ( total_polymer <= 0 )
        throw InvalidParameter("fiber:total_polymer should be > 0 (you can specify 'inf')");

    if ( glue )
    {
        if ( primed(sim) )
            glue_prop = sim.findProperty<SingleProp>("single", glue_single);
    }
    
    if ( field != "none" )
    {
        Property * p = sim.properties.find_or_die("field", field);
        field_ptr = sim.pickField(p);
    }
    
    if ( lattice && primed(sim) )
    {
#if FIBER_HAS_LATTICE
        if ( lattice_unit <= 0 )
            throw InvalidParameter("fiber:lattice_unit (known as lattice[1]) must be specified and > 0");
#else
        throw InvalidParameter("Cytosim cannot handle fiber:lattice. Please recompile after setting FIBER_HAS_LATTICE to 1");
#endif
    }

#if FIBER_HAS_DENSITY
    if ( density && primed(sim) )
    {
        if ( density_unit <= 0 )
            throw InvalidParameter("fiber:density_unit (known as mesh[1]) must be specified and > 0");
        
        if ( density_flux_speed != 0 || density_binding_rate != 0 || density_unbinding_rate != 0 )
        {
            if ( field.empty() )
                throw InvalidParameter("fiber:mesh features require fiber:field to be specified");

            if ( !field_ptr )
                throw InvalidParameter("fiber:field not found");
        }
        
        if ( density_aging_rate < 0 )
            throw InvalidParameter("fiber:density_aging_rate must be >= 0");
        if ( density_aging_limit < 0 )
            throw InvalidParameter("fiber:density_aging_limit must be >= 0");
        if ( density_aging_rate * time_step(sim) > 1 )
            throw InvalidParameter("fiber:density_aging_rate is too high (unstable)");
    }

    if ( density_aging_rate > 0 && !density && primed(sim) )
        throw InvalidParameter("for `density_aging_rate', the mesh must be defined");
#endif

    if ( rigidity < 0 )
        throw InvalidParameter("fiber:rigidity must be specified and >= 0");
    
    if ( segmentation <= 0 )
        throw InvalidParameter("fiber:segmentation must be > 0");

    if ( primed(sim) )
    {
        // Adjust the segmentation of all Fibers having this FiberProp
        for ( Fiber* fib = sim.fibers.first(); fib; fib=fib->next() )
            if ( fib->property() == this )
                fib->adjustSegmentation(segmentation);
    }
    
    if ( steric_key && steric_radius <= 0 )
        throw InvalidParameter("fiber:steric[1] (radius) must be specified and > 0");
    
    if ( drag_radius <= 0 )
        throw InvalidParameter("fiber:drag_radius must be > 0");
    
    if ( drag_length <= 0 )
        throw InvalidParameter("fiber:drag_length must be > 0");

#if NEW_FIBER_END_CHEW
    if ( max_chewing_speed < 0 )
        throw InvalidParameter("fiber:max_chewing_speed must be >= 0");
    max_chewing_speed_dt = max_chewing_speed * time_step(sim);
#endif
    
#if ( 0 )
    //print information relating to the 'numerical stiffness' of the system
    Fiber fib(this);
    fib.setStraight(Vector(-5,0,0), Vector(1,0,0), 10);
    
    fib.setDragCoefficient();
    real mob_dt = time_step(sim) * fib.pointMobility();
    
    real stiffness = 100;
    real coef1 = mob_dt * stiffness;

    Cytosim::log.print("Numerical hardness (stiffness=%.1f): %7.2f\n", stiffness, coef1);

    real rod   = segmentation;
    real coef2 = mob_dt * rigidity / ( rod * rod * rod );
    
    Cytosim::log.print("Numerical hardness (rigidity=%.1f): %7.2f\n", rigidity, coef2);
#endif
}


//------------------------------------------------------------------------------

void FiberProp::write_values(std::ostream& os) const
{
    write_value(os, "rigidity",            rigidity);
    write_value(os, "segmentation",        segmentation);
    write_value(os, "min_length",          min_length);
    write_value(os, "max_length",          max_length);
    write_value(os, "total_polymer",       total_polymer);
    write_value(os, "persistent",          persistent);
    write_value(os, "viscosity",           viscosity);
    write_value(os, "drag_radius",         drag_radius);
    write_value(os, "drag_length",         drag_length);
    write_value(os, "drag_model",          drag_model, drag_gap);
#if NEW_SQUEEZE_FORCE
    write_value(os, "squeeze",             squeeze_mode, squeeze_force, squeeze_range);
#endif
#if NEW_FIBER_END_FORCE
    write_value(os, "end_force",           end_force_mode, end_force);
#endif
    write_value(os, "binding_key",         binding_key);
    write_value(os, "lattice",             lattice, lattice_unit);
#if FIBER_HAS_DENSITY
    write_value(os, "density",                density, density_unit);
    write_value(os, "density_cut_fiber",      density_cut_fiber);
    write_value(os, "density_flux_speed",     density_flux_speed);
    write_value(os, "density_binding_rate",   density_binding_rate);
    write_value(os, "density_unbinding_rate", density_unbinding_rate);
    write_value(os, "density_aging", density_aging_rate, density_aging_limit);
#endif
    write_value(os, "confine", confine, confine_stiff[0], confine_stiff[1], confine_spec);
#if NEW_FIBER_CONFINE2
    write_value(os, "confine2", confine2, confine2_stiff[0], confine2_stiff[1], confine2_spec);
#endif
    write_value(os, "steric",              steric_key, steric_radius, steric_range);
    write_value(os, "field",               field);
    write_value(os, "glue",                glue, glue_single);
#if NEW_COLINEAR_FORCE
    write_value(os, "colinear_force",      colinear_force);
#endif
#if NEW_FIBER_END_CHEW
    write_value(os, "max_chewing_speed",   max_chewing_speed);
#endif
#if NEW_FIBER_LOOP
    write_value(os, "loop",                loop);
#endif
    write_value(os, "activity",            activity);
    write_value(os, "display",             "("+display+")");
}

