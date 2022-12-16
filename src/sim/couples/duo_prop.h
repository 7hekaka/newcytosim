// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef DUO_PROP_H
#define DUO_PROP_H

#include "couple_prop.h"
class Space;
class SolidProp;


/// Additional Property for Duo and DuoLong
/**
 @ingroup Properties
*/
class DuoProp : public CoupleProp
{
    
    friend class Duo;
    
public:
    
    /**
     @defgroup DuoPar Parameters of Duo
     @ingroup Parameters
     Inherits @ref CouplePar.
     @{
     */
    
    /// name of the Space inside which the Duo is activated
    std::string  activation;

    /// rate of deactivation
    real deactivation_rate;
    
    /// type of deactivation: can lead to object deletion
    int deactivation_type;
    
    /// if true, the deactivation clock runs at all time
    bool vulnerable;
    
    /// @}

    /// deactivation_rate * time_step
    real deactivation_rate_dt;
    
    // Space inside which the Duo is activated
    Space const* activation_space;
    
    // Space inside which the Duo is activated
    SolidProp const* activation_beads;

public:
    
    /// constructor
    DuoProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~DuoProp() { }
    
    /// return a Duo or a DuoLong with this property
    Couple * newCouple(Glossary*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new DuoProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;

};

#endif

