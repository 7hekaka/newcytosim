// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "cymdef.h"
#include "messages.h"
#include "glossary.h"
#include "exceptions.h"
#include "simul_prop.h"
#include "digit_prop.h"
#include "digit.h"
#include "fiber.h"

Hand * DigitProp::newHand(HandMonitor* m) const
{
    return new Digit(this, m);
}


void DigitProp::clear()
{
    HandProp::clear();
}


void DigitProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    if ( glos.set(step_size, "step_size") && site_shift == 0 )
        site_shift = 0.5 * step_size;
    glos.set(site_shift, "site_shift");
    
    unsigned long i = 0;
    if ( glos.set(i, "footprint") )
    {
        footprint = i;
        if ( (unsigned long)footprint != i )
            throw InvalidParameter("out of range `footprint` value");
#if FIBER_HAS_LATTICE
        if ( std::is_same<real, FiberLattice::cell_t>::value )
            Cytosim::warn("`footprint` is only valid with Integer-based fiber:lattice\n");
#else
        Cytosim::warn("with no fiber:lattice support, `footprint` is ignored\n");
#endif
    }

#if BACKWARD_COMPATIBILITY < 100
    bool u = true;
    if ( glos.set(u, "use_lattice") && !u )
        throw InvalidParameter("`use_lattice` is deprecated: set footprint=0");
#endif
}


void DigitProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( step_size <= 0 )
        throw InvalidParameter(name()+":step_size must be defined and > 0");
    
    if ( site_shift < 0 || step_size < site_shift )
        throw InvalidParameter(name()+":site_shift must be in [0, step_size]");
}


void DigitProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "step_size", step_size);
    write_value(os, "footprint", std::to_string(footprint));
    write_value(os, "site_shift", site_shift);
}

