// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event_set.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "event.h"


Event * EventSet::first() const
{
    return static_cast<Event*>(pool_.front());
}

Event * EventSet::findID(ObjectID n) const
{
    return static_cast<Event*>(inventory_.get(n));
}

void EventSet::step()
{
    for ( Event * e=first(); e; e=e->next() )
        e->step(simul_);
}


Property* EventSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    return nullptr;
}


Object * EventSet::newObject(const ObjectTag tag, size_t num)
{
    if ( tag == Event::TAG )
        return new Event();
    std::cerr << "Warning: unknown Event tag `"+std::string(1,tag)+"' requested\n";
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
ObjectList EventSet::newObjects(const std::string&, Glossary& opt)
{
    ObjectList res;
    Event * e = new Event(simul_.time(), opt);
    res.push_back(e);
    return res;
}


void EventSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.writeLine("\n#section "+title());
        writeObjects(out, pool_);
    }
}

