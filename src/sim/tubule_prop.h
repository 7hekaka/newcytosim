// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#ifndef TUBULE_PROP_H
#define TUBULE_PROP_H

#include "property.h"
#include "real.h"

class FiberSet;


/// Property for Tubule
/**
 @ingroup Properties
 
 */
class TubuleProp : public Property
{
    friend class Tubule;
    
public:
    
    /**
     @defgroup TubulePar Parameters of Tubule
     @ingroup Parameters
     @{
     */

    /// stiffness of links between the Solid and the Fiber
    /**
     - stiffness[0] is for the link between protofilaments
     - stiffness[1] is for the angular link.
     .
    */
    real          stiffness[2];
    
    /// name of Fiber that make up the Tubule (know as `filament[1]`)
    std::string   fiber_type;

    /// name of Fiber used to make the backbone
    std::string   bone_type;
    
    /// distance between centerline of protofilament and central axis of tubule
    /* defines the distance between the centerline of the protofilaments
     `radius=0.0135` gives a diameter of 25+2nm, where the 2 represents the position
     of the center of molecules that attach to the tubulin lattice
     */
    real          radius;
    
    /// @}
    
public:
    
    /// constructor
    TubuleProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~TubuleProp() { }
    
    /// identifies the property
    std::string category() const { return "tubule"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// check and derive parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new TubuleProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

