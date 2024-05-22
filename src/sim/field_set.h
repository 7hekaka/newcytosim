// Cytosim was created by Francois Nedelec. Copyright 2024

#ifndef FIELD_SET_H
#define FIELD_SET_H

#include "object_set.h"
#include "field.h"

class Simul;

/// a list of Field
/**
 FieldSet contains all the 'Field' in the simulation

 */
class FieldSet : public ObjectSet
{
public:
    
    /// creator
    FieldSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "field"; }
    
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
    Field * first() const { return static_cast<Field*>(pool_.front()); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Field * identifyObject(ObjectID n) const { return static_cast<Field*>(inventory_.get(n)); }
    
    /// get ready to do steps()
    void prepare();
    
    /// Monte-Carlo simulation step for every Object
    void steps();

};


#endif
