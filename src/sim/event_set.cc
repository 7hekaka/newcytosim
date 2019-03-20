// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event_set.h"
#include "event_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"

//------------------------------------------------------------------------------

void EventSet::prepare()
{
    for ( Event * f=first(); f; f=f->next() )
    {
        f->prepare();
    }
}


void EventSet::step()
{
    for ( Event * f=first(); f; f=f->next() )
    {
        f->step();
    }
}

//------------------------------------------------------------------------------
#pragma mark -

Property* EventSet::newProperty(const std::string& kd, const std::string& nm, Glossary&) const
{
    if ( kd == "event" )
        return new EventProp(nm);
    return nullptr;
}


Object * EventSet::newObjectT(const ObjectTag tag, unsigned idx)
{
    if ( tag == Event::TAG )
    {
        EventProp * p = simul.findProperty<EventProp>("event", idx);
        return p->newEvent();
    }
    return nullptr;
}


/**
 @defgroup NewEvent How to create an Event
 @ingroup NewObject

 Specify a new Event:
 
     new event NAME
     {

     }
 */
ObjectList EventSet::newObjects(const std::string& name, Glossary& opt)
{
    EventProp * p = simul.findProperty<EventProp>("event", name);
    Event * obj = p->newEvent();

    ObjectList res;
    res.push_back(obj);
    return res;
}

