// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SET_H
#define SPACE_SET_H

#include "object_set.h"
#include "space.h"
class Simul;

///a list of Space
class SpaceSet : public ObjectSet
{
    /// the master space
    Space const* master_;

public:

    /// constructor
    SpaceSet(Simul& s) : ObjectSet(s), master_(nullptr) {}

    /// return master
    Space const* master() const { return master_; }

    /// change master
    void setMaster(Space const* s);
    
    //--------------------------
    
    /// identifies the property
    static std::string title() { return "space"; }
    
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
    
    /// add Object
    void link(Object *);
    
    /// remove Object
    void unlink(Object *);

    /// erase all objects
    void erase();
    
    /// Monte-Carlo step for every Space
    void steps();
    
    /// first Space
    Space * first() const { return static_cast<Space*>(pool_.front()); }
    
    /// last Space
    Space * last() const { return static_cast<Space*>(pool_.back()); }
    
    /// first Fiber in inventory
    Space * firstID() const { return static_cast<Space*>(inventory_.first()); }

    /// returns Space immediately following 'obj' in inventory
    Space * nextID(Space const* obj) const { return static_cast<Space*>(inventory_.next(obj)); }

    /// return pointer to the Object of given ID, or zero if not found
    Space * identifyObject(ObjectID n) const { return static_cast<Space*>(inventory_.get(n)); }

    /// return max extension over all Spaces
    real maxExtension() const;
    
};


#endif

