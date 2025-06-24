// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ORGANIZER_SET_H
#define ORGANIZER_SET_H

#include "object_set.h"
#include "organizer.h"

class Aster;

/// a list of Organizer (Aster, Nucleus, Bundle)
class OrganizerSet : public ObjectSet
{
public:
    
    ///creator
    OrganizerSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "organizer"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// write all Objects to file
    void writeSet(Outputter&) const;
        
    /// print a summary of the content (nb of objects, class)
    void report(std::ostream&) const;

    //--------------------------
    
    /// first Organizer
    Organizer * first() const { return static_cast<Organizer*>(pool_.front()); }
    
    /// find object with given ID
    Organizer * identifyObject(ObjectID n) const { return static_cast<Organizer*>(inventory_.get(n)); }
    
    /// find highest ObjectID among Organizers containing given Mecable
    ObjectID findOrganizerID(Mecable const*) const;

    /// Monte-Carlo simulation step for every Object
    void steps();

    //--------------------------

    /// first Aster with this name
    Aster * pickAster(std::string) const;

};


#endif


