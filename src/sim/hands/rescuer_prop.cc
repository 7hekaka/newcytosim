// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "exceptions.h"
#include "glossary.h"
#include "rescuer_prop.h"
#include "rescuer.h"


Hand * RescuerProp::newHand(HandMonitor* m) const
{
    return new Rescuer(this, m);
}


void RescuerProp::clear()
{
    HandProp::clear();

    rescue_chance = 0;
}


void RescuerProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(rescue_chance, "rescue_prob", "rescue_probability");
    glos.set(rescue_chance, "rescue_chance");
}


void RescuerProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( rescue_chance < 0 )
        throw InvalidParameter("rescuer:rescue_chance must be >= 0");
}


void RescuerProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "rescue_chance", rescue_chance);
}

