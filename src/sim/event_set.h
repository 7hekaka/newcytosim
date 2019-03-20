// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
    std::string title() const { return "event"; }
    
    /// create a new property for class `kind` with given name
    Property *  newProperty(const std::string& kind, const std::string& name, Glossary&) const;
    
    /// create new objects, of class `kind` and type `name`, give the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(ObjectTag, unsigned);
    
    //--------------------------
    
    /// first object
    Event *     first() const
    {
        return static_cast<Event*>(nodes.front());
    }
    
    ///  return pointer to the Object of given ID, or zero if not found
    Event *     findID(ObjectID n) const
    {
        return static_cast<Event*>(inventory.get(n));
    }
    
    /// get ready to do a step()
    void        prepare();
    
    /// Monte-Carlo simulation step for every Object
    void        step();

};


#endif
