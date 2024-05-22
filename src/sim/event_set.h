// Cytosim was created by Francois Nedelec. Copyright 2024

#ifndef EVENT_SET_H
#define EVENT_SET_H

#include "object_set.h"
#include "event.h"

class Simul;

/// a list of Event
/**
 */
class EventSet : public ObjectSet
{
public:
    
    /// creator
    EventSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "event"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// write all Objects to file
    void writeSet(Outputter&) const;
        
    /// print a summary of the content (nb of objects, class)
    void report(std::ostream& os) const;

    //--------------------------
    
    /// first object
    Event * first() const { return static_cast<Event*>(pool_.front()); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Event * identifyObject(ObjectID n) const { return static_cast<Event*>(inventory_.get(n)); }
    
    /// Monte-Carlo simulation step for every Object
    void steps();

};


#endif
