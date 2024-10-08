// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DIGIT_PROP_H
#define DIGIT_PROP_H

#include "hand_prop.h"
#include "fiber.h"

/// Additional Property for Digit
/**
 @ingroup Properties
*/
class DigitProp : public HandProp
{
    friend class Digit;
    
public:
    
    /**
     @defgroup DigitPar Parameters of Digit
     @ingroup Parameters
     Inherits @ref HandPar.
     @{
     */
    
    /// @}
    
public:

    /// constructor
    DigitProp(const std::string& n) : HandProp(n)  { clear(); }
    
    /// destructor
    ~DigitProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new DigitProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

