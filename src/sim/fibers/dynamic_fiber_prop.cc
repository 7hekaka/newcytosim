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
#if OLD_DYNAMIC_ZONE
    zone_radius = INFINITY;
    zone_space  = "";
    for ( int i = 0; i < 2; ++i )
        zone_hydrolysis_rate[i] = 0;
#endif
}


/**
 Estimate the life time and length of fiber using formula from:
 A theory of microtubule catastrophes and their regulation</b>\n
 Brun L, Rupp B, Ward J, Nedelec F\n
 PNAS 106 (50) 21173-21178; 2009\n
 */
static void splash(std::ostream& os, real g, real h, real unit)
{
    real ctime = ( 7*h*h + 12*g*h + 3*g*g ) / ( 3*h*h * ( 2*h + 3*g ) );
    real len = g * unit * ctime;
    
    //const real ctime = g / ( 3*h*h );  // that is only true if g >> h
    std::streamsize p = os.precision();
    os.precision(5);
    os << "  DynamicFiber h " << h << " g " << g << " :";
    os << " catastrophe_time " << ctime << "  rate " << 1/ctime;
    os << " length " << len << "\n";
    os.precision(p);
}

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
#if OLD_DYNAMIC_ZONE
    glos.set(zone_space,              "zone_space");
    glos.set(zone_radius,             "zone_radius");
    glos.set(zone_hydrolysis_rate, 2, "zone_hydrolysis_rate");
#endif
#if BACKWARD_COMPATIBILITY < 44
    if ( glos.set(growing_force[0], "dynamic_force") )
        Cytosim::warn << "fiber:dynamic_force was renamed growing_force\n";
    
    int f = 0;
    if ( glos.set(f, "fate", {{"none", 0}, {"destroy", 1}, {"rescue", 2}}))
    {
        Cytosim::warn << "fiber:fate is deprecated: use `persistent` and `rebirth_rate`\n";
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
        splash(std::clog, g, h, unit_length);
    }
}


void DynamicFiberProp::complete(Simul const& sim)
{
    FiberProp::complete(sim);

    for ( int i = 0; i < 2; ++i )
    {
        if ( growing_force[i] <= 0 )
            throw InvalidParameter("fiber:growing_force should be > 0");
        growing_force_inv[i] = 1.0 / growing_force[i];

        if ( growing_speed[i] < 0 )
            throw InvalidParameter("fiber:growing_speed should be >= 0");
        growing_rate_dt[i] = sim.time_step() * abs_real(growing_speed[i]) / unit_length;

        if ( growing_off_speed[i] > 0 )
            throw InvalidParameter("growing_off_speed should be <= 0");
        growing_off_rate_dt[i] = -sim.time_step() * growing_off_speed[i] / unit_length;

        if ( hydrolysis_rate[i] < 0 )
            throw InvalidParameter("fiber:hydrolysis_rate should be >= 0");
        hydrolysis_rate_2dt[i] = 2 * sim.time_step() * hydrolysis_rate[i];
        
#if OLD_DYNAMIC_ZONE
        zone_space_ptr = sim.findSpace(zone_space);
        if ( zone_radius < 0 )
            throw InvalidParameter("fiber:zone_radius should be >= 0");
        if ( zone_hydrolysis_rate[i] < 0 )
            throw InvalidParameter("fiber:zone_hydrolysis_rate should be >= 0");
        zone_hydrolysis_rate_2dt[i] = 2 * sim.time_step() * zone_hydrolysis_rate[i];

        zone_radius_sqr = square(zone_radius);
#endif

        if ( shrinking_speed[i] > 0 )
            throw InvalidParameter("fiber:shrinking_speed should be <= 0");
        shrinking_rate_dt[i] = sim.time_step() * abs_real(shrinking_speed[i]) / unit_length;
        
        if ( rebirth_rate[i] < 0 )
            throw InvalidParameter("fiber:rebirth_rate should be >= 0");
        rebirth_prob[i] = -std::expm1( -rebirth_rate[i] * sim.time_step() );

        if ( unhydrolyzed_prob[i] < 0 || unhydrolyzed_prob[i] > 1 )
			throw InvalidParameter("fiber:unhydrolyzed_prob should be in [0, 1]");
    }

    if ( min_length <= 0 )
        min_length = unit_length;
     
    /// print predicted average length in verbose mode:
    if ( sim.primed() && sim.prop->verbose )
    {
        if ( 0 == growing_off_speed[0] )
        splash(std::clog, growing_speed[0]/unit_length, hydrolysis_rate[0], unit_length);
        
        // calculate stall force, from:
        // 0 = growing_speed * std::exp(force/growing_force) + growing_off_speed;
        real f = -growing_force[0] * std::log(-growing_off_speed[0]/growing_speed[0]);
        std::clog << name() << ":stall_force = " << f << "\n";
    }
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
#if OLD_DYNAMIC_ZONE
    write_value(os, "zone_space",           zone_space);
    write_value(os, "zone_radius",          zone_radius);
    write_value(os, "zone_hydrolysis_rate", zone_hydrolysis_rate, 2);
#endif
}

