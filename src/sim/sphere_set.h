// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPHERE_SET_H
#define SPHERE_SET_H

#include "object_set.h"
#include "sphere.h"
class Simul;

///a list of Sphere
class SphereSet : public ObjectSet
{
public:
    
    ///creator
    SphereSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "sphere"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// write all Objects to file
    void writeSet(Outputter&) const;
        
    /// print a summary of the content (nb of objects, class)
    void report(std::ostream& os) const { writeReport(os, title()); }

    //--------------------------
   
    /// remove object
    void remove(Object *);

    /// first Object
    Sphere * first() const { return static_cast<Sphere*>(pool_.front()); }
    
    /// first Sphere in inventory
    Sphere * firstID() const { return static_cast<Sphere*>(inventory_.first()); }
    
    /// returns Sphere immediately following 'obj' in inventory
    Sphere * nextID(Sphere const* obj) const { return static_cast<Sphere*>(inventory_.next(obj)); }

    /// return pointer to the Object of given ID, or zero if not found
    Sphere * identifyObject(ObjectID n) const { return static_cast<Sphere*>(inventory_.get(n));}
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;
    
    /// Monte-Carlo simulation step for every Object
    void steps() {}
 };

#endif
