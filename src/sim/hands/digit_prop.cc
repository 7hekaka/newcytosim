// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "exceptions.h"
#include "glossary.h"
#include "common.h"
#include "property_list.h"
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

    step_size = 0;
    footprint = 1;
    site_pos = 0.5;
}


void DigitProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(step_size, "step_size");
    if ( !std::is_same<FiberLattice::cell_t, double>::value )
    {
        unsigned long i = 0;
        if ( glos.set(i, "footprint") )
            footprint = i;
    }
    glos.set(site_pos,  "site_pos");
    
#ifdef BACKWARD_COMPATIBILITY
    bool u = true;
    if ( glos.set(u, "use_lattice") && !u )
        throw InvalidParameter("`use_lattice` is deprecated: set footprint=0");
#endif
}


void DigitProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( step_size <= 0 )
        throw InvalidParameter("Digit:step_size must be defined and > 0");
    
    if ( site_pos < 0 || 1 < site_pos )
        throw InvalidParameter("Digit:site_pos must be in [0, 1]");
}


void DigitProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "step_size", step_size);
    write_value(os, "footprint", footprint);
    write_value(os, "site_pos", site_pos);
}

