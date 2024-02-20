// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_PROP_H
#define SPACE_PROP_H

#include "real.h"
#include "property.h"
#include "vector.h"

class Glossary;
class PointDisp;
class Space;
class Simul;


/// Property for Space
/**
 @ingroup Properties
 */
class SpaceProp : public Property
{
    friend class Space;

public:

    /**
     @defgroup SpacePar Parameters of Space
     @ingroup Parameters
     @{
    */
    
    /// primitive (e.g. `rectangle`)
    std::string shape;
    
    /// name of the file where dimensions are stored
    std::string dimensions;
    
    /// display string (see @ref PointDispPar)
    std::string display;

    /// @}
    
    /// derived variable: flag to indicate that `display` has a new value
    bool display_fresh;

    /// derived variable: parameters extracted from `display`
    PointDisp * disp;

public:

    /// constructor
    SpaceProp(const std::string& n) : Property(n), disp(nullptr) { clear(); }

    /// destructor
    ~SpaceProp() { }
    
    /// create a new, uninitialized, Space
    virtual Space * newSpace() const;

    /// create a new Space according to specifications
    virtual Space * newSpace(Glossary&) const;
    
    /// identifies the property
    std::string category() const { return "space"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// check and derive some parameters
    void complete();
 
    /// check and derive more parameters
    void complete(Simul const&);
            
    /// return a carbon copy of object
    Property* clone() const { return new SpaceProp(*this); }
    
    /// write all values
    virtual void write_values(std::ostream&) const;
    
};

#endif

