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


Hand * DigitProp::newHand(HandMonitor* m) const
{
    return new Digit(this, m);
}


void DigitProp::clear()
{
    HandProp::clear();

    step_size = 1;
}


void DigitProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(step_size, "step_size");
}


void DigitProp::complete(Simul const& sim)
{
    HandProp::complete(sim);

    if ( step_size <= 0 )
        throw InvalidParameter("Digit:step_size must be defined and > 0");

}


void DigitProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "step_size", step_size);
}

