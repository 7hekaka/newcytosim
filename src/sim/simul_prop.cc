// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "simul_prop.h"
#include "assert_macro.h"
#include "space_prop.h"
#include "space_set.h"
#include "simul.h"
#include "exceptions.h"
#include "messages.h"
#include "glossary.h"
#include "random.h"


//------------------------------------------------------------------------------
void SimulProp::clear()
{
    viscosity         = 1;
#if NEW_CYTOPLASMIC_FLOW
    flow.reset();
#endif
#if OLD_SPINDLE_FLUX
    flux_speed        = 0;
#endif
    time              = 0;
    time_step         = 0;
    kT                = 0.0042;
    tolerance         = 0.05;
    acceptable_prob   = 0.5;
    precondition      = 0;
    precondition_span = 2;
    random_seed       = 0;
    
    steric_mode       = 0;
    steric_stiff_push[0] = 0;
    steric_stiff_pull[0] = 0;
    steric_stiff_push[1] = 0;
    steric_stiff_pull[1] = 0;

    steric_max_range  = -1;
    binding_grid_step = -1;
    
    verbose = 0;

    config_file      = "config.cym";
    property_file    = "properties.cmp";
    system_file      = TRAJECTORY;
    clear_trajectory = true;
    
    skip_free_single = 0;
    skip_free_couple = 0;

    display       = "";
    display_fresh = false;
}


void SimulProp::read(Glossary& glos)
{
    // a dimensionality can be specified to stop the program from running
    unsigned d = DIM;
    if ( glos.set(d, "dimension", "dim") && d != DIM )
    {
        std::cerr << "Aborting since the config file specifies a different dimensionality:\n";
        std::cerr << "Cytosim was compiled with DIM="<<DIM<<" but the config specifies `dim="<<d<<"'\n";
        exit(0);
        //throw InvalidParameter("dimensionality missmatch");
    }
    
    glos.set(viscosity,         "viscosity");
#if NEW_CYTOPLASMIC_FLOW
    glos.set(flow,              "flow");
#endif
#if OLD_SPINDLE_FLUX
    glos.set(flux_speed,        "flux_speed");
#endif
    glos.set(time,              "time");
    glos.set(time_step,         "time_step");
    
    real T = 0;
    if ( glos.set(T, "temperature" ) )
        kT = 1.38064852e-5 * T;
    glos.set(kT, "kT", "thermal_energy");

    glos.set(tolerance,         "tolerance");
    glos.set(acceptable_prob,   "acceptable_prob");
    glos.set(precondition,      "precondition");
    glos.set(precondition_span, "precondition", 1);
    
    glos.set(steric_mode,          "steric", {{"off", 0}, {"on", 1}});
    glos.set(steric_stiff_push[0], "steric", 1);
    glos.set(steric_stiff_pull[0], "steric", 2);
    glos.set(steric_stiff_push, 2, "steric_stiff_push");
    glos.set(steric_stiff_pull, 2, "steric_stiff_pull");
    glos.set(steric_max_range,     "steric_max_range");

    glos.set(binding_grid_step, "binding_grid_step");
    
    // these parameters are not written:
    glos.set(verbose,           "verbose");
    
    // names of files and path:
    glos.set(config_file,       "config");
    glos.set(config_file,       ".cytosim"); // fullname extension
    glos.set(config_file,       ".cym");
    
    glos.set(property_file,     "property_file");
    glos.set(property_file,     "property");
    
#if BACKWARD_COMPATIBILITY < 100
    glos.set(system_file, "object_file", "trajectory");
    bool a = false;
    if ( glos.set(a, "append_file") )
        clear_trajectory = !a;
#endif

    glos.set(system_file, "system_file", "system");
    glos.set(clear_trajectory, "clear_trajectory");
    
    glos.set(skip_free_single, "skip_free_single");
    glos.set(skip_free_couple, "skip_free_couple");
    glos.set(random_seed, "random_seed");
    
    if ( glos.set(display, "display") )
        display_fresh = true;
}


/**
 If the Global parameters have changed, we update all derived parameters.
 This makes it possible to change the time-step in the middle of a config file.
 */
void SimulProp::complete(Simul const& sim)
{
    // initialize the random number generator:
    if ( !RNG.seeded() )
    {
        if ( random_seed )
            RNG.seed(random_seed);
        else
            random_seed = RNG.seed();
    }

    if ( time_step < 0 )
        throw InvalidParameter("simul:time_step must be > 0");

    if ( sim.primed() )
    {
        if ( viscosity <= 0 )
            throw InvalidParameter("simul:viscosity must be > 0");

        if ( time_step <= 0 )
            throw InvalidParameter("simul:time_step must be specified and > 0");
        
        if ( kT < 0 )
            throw InvalidParameter("simul:kT must be > 0");
        
        if ( kT == 0 && tolerance > 0.1 )
            throw InvalidParameter("if simul:kT==0, simul:tolerance must be defined and small");
        
        if ( steric_stiff_push[0] < 0 )
            throw InvalidParameter("steric stiffness (push, steric[1]) must be >= 0");

        if ( steric_stiff_pull[0] < 0 )
            throw InvalidParameter("steric stiffness (pull steric[2]) must be >= 0");
    }
    /*
     If the Global parameters have changed, we update all derived parameters.
     To avoid an infinite recurence, the main SimulProp (*this) was
     not included in the PropertyList Simul::properties;
     */
    sim.properties.complete(sim);
}

//------------------------------------------------------------------------------

void SimulProp::write_values(std::ostream& os) const
{
    //write_value(os, "time",            time);
    write_value(os, "time_step",       time_step);
    write_value(os, "kT",              kT);
    write_value(os, "viscosity",       viscosity);
#if NEW_CYTOPLASMIC_FLOW
    write_value(os, "flow", flow);
#endif
#if OLD_SPINDLE_FLUX
    write_value(os, "flux_speed", flux_speed);
#endif
    std::endl(os);
    write_value(os, "tolerance",       tolerance);
    write_value(os, "acceptable_prob", acceptable_prob);
    write_value(os, "precondition",    precondition, precondition_span);
    write_value(os, "random_seed",     random_seed);
    std::endl(os);
    write_value(os, "steric", steric_mode, steric_stiff_push[0], steric_stiff_pull[0]);
    write_value(os, "steric_max_range",  steric_max_range);
    write_value(os, "binding_grid_step", binding_grid_step);
    write_value(os, "verbose", verbose);
    std::endl(os);
    write_value(os, "display", "("+display+")");
}

