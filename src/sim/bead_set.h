// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef BEAD_SET_H
#define BEAD_SET_H

#include "object_set.h"
#include "bead.h"

class Simul;

/// a list of Bead
class BeadSet : public ObjectSet
{
public:
    
    /// creator
    BeadSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "bead"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// write all Objects to file
    void writeSet(Outputter&) const;
    
    /// like defrost() but also delete Wrist attached to deleted Beads
    void defrostMore();

    /// print a summary of the content (nb of objects, class)
    void report(std::ostream& os) const { writeReport(os, title()); }

    //--------------------------
    
    /// remove from the list
    void remove(Object *);
    
    /// first Object
    Bead * first() const { return static_cast<Bead*>(pool_.front()); }
    
    /// first Bead in inventory
    Bead * firstID() const { return static_cast<Bead*>(inventory_.first()); }
    
    /// returns Bead immediately following 'obj' in inventory
    Bead * nextID(Bead const* obj) const { return static_cast<Bead*>(inventory_.next(obj)); }

    /// find object from its Number
    Bead * identifyObject(ObjectID n) const { return static_cast<Bead*>(inventory_.get(n)); }
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;
    
    /// Monte-Carlo simulation step for every Object
    void steps();
};


#endif

