// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef BUNDLE_PROP_H
#define BUNDLE_PROP_H

#include "cymdef.h"
#include "real.h"
#include "property.h"

class FiberSet;


/// Property for Bundle
/**
 @ingroup Properties
*/
class BundleProp : public Property
{
    friend class Bundle;
    
public:
    
    /**
     @defgroup BundlePar Parameters of Bundle
     @ingroup Parameters
     @{
     */
 
    /// stiffness of the links that connect the overlapping fibers
    real stiffness;
    
    /// length of the zone where fibers overlap
    real overlap;
    
    /// designates which end of the fiber is towards the center
    FiberEnd pole;
    
    /// if true, the bundle is made of filaments alternating in polarity
    bool bipolar;

    /// rate for creating fiber in empty slots
    real fiber_rate;
    
    /// name of Fiber in the Bundle
    std::string fiber_type;
    
    /// name of Fiber in the Bundle
    std::string fiber_spec;

    /// @}
    
    /// probability of nucleation
    real fiber_prob;

public:
 
    /// constructor
    BundleProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~BundleProp() { }
    
    /// identifies the property
    std::string category() const { return "bundle"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// check and derive parameters
    void complete(Simul const&);
    
    
    /// return a carbon copy of object
    Property* clone() const { return new BundleProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

