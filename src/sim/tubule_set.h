// Cytosim was created by Francois Nedelec. Copyright 2024

#ifndef TUBULE_SET_H
#define TUBULE_SET_H

#include "object_set.h"
#include "tubule.h"

class Simul;

/// a list of Tubules
/**
 */
class TubuleSet : public ObjectSet
{
public:
    
    /// creator
    TubuleSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "tubule"; }
    
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

    /// first object
    Tubule * first() const { return static_cast<Tubule*>(pool_.front()); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Tubule * findID(ObjectID n) const { return static_cast<Tubule*>(inventory_.get(n)); }

    /// Monte-Carlo simulation step for every Object
    void steps();

};


#endif
