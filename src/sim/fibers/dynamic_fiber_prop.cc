// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "cymdef.h"
#include "dynamic_fiber_prop.h"
#include "dynamic_fiber.h"
#include "simul_prop.h"
#include "exceptions.h"
#include "glossary.h"
#include "messages.h"
#include "simul.h"


Fiber* DynamicFiberProp::newFiber() const
{
    return new DynamicFiber(this);
}

/**
 This rounds up the value given by FiberProp::newFiberLength()
 to a multiple of the unit_length
 */
real DynamicFiberProp::newFiberLength(Glossary& opt) const
{
    const real uni = unit_length;
    const real len = FiberProp::newFiberLength(opt);
    
    index_t i = std::round( len / uni );

    while ( i * uni > max_length)
        --i;

    while ( i * uni < min_length)
        ++i;
    
    real dif = len - i * uni;
    //printf("length adjusted by: %.6e\n", dif);

    if ( abs_real(dif) > 0.001 )
        Cytosim::log(name(), ":length rounded up to ", i, " x ", uni, " = ", i*uni, "\n");
    
    return i * uni;
}


void DynamicFiberProp::clear()
{
    FiberProp::clear();
    
    // we use the tubulin heterodimer length by default:
    unit_length = 0.008;
    
    for ( int i = 0; i < 2; ++i )
    {
        growing_speed[i]     = 0;
        growing_off_speed[i] = 0;
        growing_force[i]     = INFINITY;
        hydrolysis_rate[i]   = 0;
        shrinking_speed[i]   = 0;
        rebirth_rate[i]      = 0;
        unhydrolyzed_prob[i] = 0;
    }
#if NEW_STALL_OUTSIDE
    stall_outside = 1;
    stall_label = "";
    stall_space = nullptr;
#endif
}


/**
 Estimate the life time and length of fiber using formula from:
 A theory of microtubule catastrophes and their regulation</b>\n
 Brun L, Rupp B, Ward J, Nedelec F\n
 PNAS 106 (50) 21173-21178; 2009\n
 */
static void splashGHU(std::ostream& os, real g, real h, real unit)
{
    real ctime = ( 7*h*h + 12*g*h + 3*g*g ) / ( 3*h*h * ( 2*h + 3*g ) );
    // if h << g, ctime ~ ( g + 4*h ) / ( 3*h*h );
    real len = g * unit * ctime;
    
    std::streamsize p = os.precision(5);
    os << " hydrolysis " << h << "/s growth " << g << "/s";
    os << "  catastrophe_time " << ctime << "s  rate " << 1/ctime << "/s";
    os << "  length " << len << "um \n";
    os.precision(p);
}


void DynamicFiberProp::splash(std::ostream& os) const
{
    std::ostringstream oss;
    oss << std::setw(16) << name() << ":";
    if ( 0 == growing_off_speed[0] )
    {
        splashGHU(oss, growing_speed[0]/unit_length, hydrolysis_rate[0], unit_length);
    }
    else
    {
        // calculate stall force, from:
        // 0 = growing_speed * std::exp(force/growing_force) + growing_off_speed;
        real f = -growing_force[0] * std::log(-growing_off_speed[0]/growing_speed[0]);
        oss << " stall_force " << f;
    }
    if ( oss.str() != splashed )
    {
        splashed = oss.str();
        os << splashed;
    }
}

/** Calculate the Hydrolysis rate resulting in a lifetime == t, given the growth */
static real back_calculate(real g, real t)
{
    const real g2 = g * g;
    /* Use Newton's method to find root of:
     F = ( 7*h*h + 12*g*h + 3*g*g ) / ( 3*h*h * ( 2*h + 3*g ) );
     d = -2*(7*h^3+24*g*h^2+27*g^2*h+9*g^3) / (3*h^3*(2*h+3*g)^2)
     */
    real h = std::sqrt(0.3333/g); //t = g / ( 3*h*h );
    for ( int i = 0; i < 9; ++i )
    {
        real h2 = h * h, hg = 2*h + 3*g;
        real F = ( 7*h2 + 12*g*h + 3*g2 - t * (3*h2 * hg) ) * h * hg;
        real d = -2 * ( h2 * ( 7*h + 24*g ) + 9 * g2 * ( 3*h + g ));
        h -= F / d;
    }
    return h;
}

void DynamicFiberProp::read(Glossary& glos)
{
    FiberProp::read(glos);
    
    glos.set(unit_length,          "unit_length");
    glos.set(growing_speed,     2, "growing_speed");
    glos.set(growing_off_speed, 2, "growing_off_speed");
    glos.set(growing_force,     2, "growing_force");
    glos.set(hydrolysis_rate,   2, "hydrolysis_rate");
    glos.set(shrinking_speed,   2, "shrinking_speed");
    glos.set(rebirth_rate,      2, "rebirth_rate");
    glos.set(unhydrolyzed_prob, 2, "unhydrolyzed_prob");
#if NEW_STALL_OUTSIDE
    glos.set(stall_outside, "stall_outside");
    glos.set(stall_label, "stall_outside", 1);
#endif
#if BACKWARD_COMPATIBILITY < 44
    if ( glos.set(growing_force[0], "dynamic_force") )
        Cytosim::warn ("fiber:dynamic_force was renamed growing_force\n");
    
    int f = 0;
    if ( glos.set(f, "fate", {{"none", 0}, {"destroy", 1}, {"rescue", 2}}))
    {
        Cytosim::warn("fiber:fate is deprecated: use `persistent` and `rebirth_rate`\n");
        persistent = ( f != 1 );
        rebirth_rate[0] = ( f == 2 ? INFINITY : 0 );
    }
#endif
    /*
     If 'length' is specified, we calculate the corresponding hydrolysis_rate
     */
    real len = 0;
    if ( glos.set(len, "length") && 0 == growing_off_speed[0] )
    {
        real g = growing_speed[0]/unit_length;
        real h = back_calculate(g, len/growing_speed[0]);
        splashGHU(std::clog, g, h, unit_length);
    }
}


void DynamicFiberProp::complete(Simul const& sim)
{
    FiberProp::complete(sim);

    fiber_count = sim.fibers.count(this);
    /// print predicted average length in verbose mode:
    if ( primed(sim) && sim.prop.verbose && fiber_count )
        splash(std::clog);

    for ( int i = 0; i < 2; ++i )
    {
        if ( growing_force[i] <= 0 )
            throw InvalidParameter(name()+":growing_force should be > 0");
        growing_force_inv[i] = 1.0 / growing_force[i];

        if ( growing_speed[i] < 0 )
            throw InvalidParameter(name()+":growing_speed should be >= 0");
        growing_rate_dt[i] = time_step(sim) * abs_real(growing_speed[i]) / unit_length;

        if ( growing_off_speed[i] > 0 )
            throw InvalidParameter(name()+":growing_off_speed should be <= 0");
        growing_off_rate_dt[i] = -time_step(sim) * growing_off_speed[i] / unit_length;

        if ( hydrolysis_rate[i] < 0 )
            throw InvalidParameter(name()+":hydrolysis_rate should be >= 0");
        hydrolysis_rate_2dt[i] = 2 * time_step(sim) * hydrolysis_rate[i];

        if ( shrinking_speed[i] > 0 )
            throw InvalidParameter("fiber:shrinking_speed should be <= 0");
        shrinking_rate_dt[i] = time_step(sim) * abs_real(shrinking_speed[i]) / unit_length;
        
        if ( rebirth_rate[i] < 0 )
            throw InvalidParameter("fiber:rebirth_rate should be >= 0");
        rebirth_prob[i] = -std::expm1( -rebirth_rate[i] * time_step(sim) );

        if ( unhydrolyzed_prob[i] < 0 || unhydrolyzed_prob[i] > 1 )
            throw InvalidParameter("fiber:unhydrolyzed_prob should be in [0, 1]");
    }

    if ( min_length <= 0 )
    {
        Cytosim::log.print("fiber:min_length <-- %.3f\n", unit_length);
        min_length = unit_length;
    }
    
#if NEW_STALL_OUTSIDE
    stall_space = sim.findSpace(stall_label);
    if ( stall_space )
        stall_label = sim.spaces.nameObject(stall_space);
    else if ( primed(sim) && stall_outside != 1 )
        throw InvalidParameter("A space must be specified as stall_outside[1]");
#endif
}


void DynamicFiberProp::write_values(std::ostream& os) const
{
    FiberProp::write_values(os);
    
    write_value(os, "unit_length",          unit_length);
    write_value(os, "growing_speed",        growing_speed, 2);
    write_value(os, "growing_off_speed",    growing_off_speed, 2);
    write_value(os, "growing_force",        growing_force, 2);
    write_value(os, "hydrolysis_rate",      hydrolysis_rate, 2);
    write_value(os, "shrinking_speed",      shrinking_speed, 2);
    write_value(os, "rebirth_rate",         rebirth_rate, 2);
    write_value(os, "unhydrolyzed_prob",    unhydrolyzed_prob, 2);
#if NEW_STALL_OUTSIDE
    write_value(os, "stall_outside", stall_outside, stall_label);
#endif
}

