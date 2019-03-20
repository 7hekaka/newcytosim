// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event_prop.h"
#include "property_list.h"
#include "glossary.h"
#include "event.h"



Event * EventProp::newEvent(Glossary * opt) const
{
    return new Event(this);
}

//------------------------------------------------------------------------------
void EventProp::clear()
{
    rate = 0;
    code = "";
}


void EventProp::read(Glossary& glos)
{
    glos.set(rate,      "rate");
    glos.set(code,      "code");
}


void EventProp::complete(Simul const& sim)
{
    if ( rate < 0 )
        throw InvalidParameter("event:rate must be > 0");
}



void EventProp::write_values(std::ostream& os) const
{
    write_value(os, "rate",          rate);
    write_value(os, "code",          code);
}

