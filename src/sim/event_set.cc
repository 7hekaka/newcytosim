// Cytosim was created by Francois Nedelec. Copyright 2024

#include "event_set.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "event.h"


void EventSet::steps()
{
    Event * obj = first();
    while ( obj )
    {
        Event * nxt = obj->next();
        obj->step(simul_);
        obj = nxt;
    }
    if ( size() > 1 ) shuffle();
}


Property* EventSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    return nullptr;
}


Object * EventSet::newObject(const ObjectTag tag, PropertyID)
{
    if ( tag == Event::TAG )
        return new Event();

    throw InvalidIO("Warning: unknown Event tag `"+std::to_string(tag)+"'");
    return nullptr;
}


/**
 @defgroup NewEvent How to create an Event
 @ingroup NewObject

 Specify a new Event:
 
     new event NAME
     {
         code = CODE;
         rate = POSITIVE_REAL;
         interval = POSITIVE_REAL;
     }
 
  `rate` (inverse of time) or `interval` (time) must be specified but not both.
 */
ObjectList EventSet::newObjects(Property const*, Glossary& opt)
{
    Event * e = new Event(simul_.time(), opt);
    return ObjectList(e);
}


void EventSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}

void EventSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << "\n" << size() << " events:";
        for ( Event * e=first(); e; e=e->next() )
            os << "\n   : " << e->activity;
    }
}
