// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#ifndef SOLID_SET_H
#define SOLID_SET_H

#include "object_set.h"
#include "solid.h"

class Simul;

/// a list of Solid
class SolidSet : public ObjectSet
{
public:
    
    /// creator
    SolidSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "solid"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// write all Objects to file
    void writeSet(Outputter&) const;
    
    /// like defrost() but also delete Wrist attached to deleted Solids
    void defrostMore();

    /// print a summary of the content (nb of objects, class)
    void report(std::ostream& os) const { writeReport(os, title()); }

    //--------------------------
    
    /// remove from the list
    void remove(Object *);
    
    /// first Solid
    Solid * first() const { return static_cast<Solid*>(pool_.front()); }
    
    /// last Solid
    Solid * last() const { return static_cast<Solid*>(pool_.back()); }
    
    /// first Solid in inventory
    Solid * firstID() const { return static_cast<Solid*>(inventory_.first()); }

    /// returns Solid immediately following 'obj' in inventory
    Solid * nextID(Solid const* obj) const { return static_cast<Solid*>(inventory_.next(obj)); }

    /// return pointer to the Object of given ID, or zero if not found
    Solid * identifyObject(ObjectID n) const { return static_cast<Solid*>(inventory_.get(n)); }
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;
    
    /// returns Solid, if one of its Sphere covers the given position (`inx` is set by this function)
    Solid * insideSphere(Vector const&, real range, size_t& inx, SolidProp const*) const;
                             
    /// Monte-Carlo simulation step for every Object
    void steps();
};


#endif

