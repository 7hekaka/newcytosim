// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DIGIT_PROP_H
#define DIGIT_PROP_H

#include "hand_prop.h"


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
    
    /// size of one step
    real   step_size;
    
    ///
    FiberLattice::cell_t footprint;
    
    /// specifies the position occupied within the Lattice site
    /**
     `site_pos` should be in [0, 1]:
     - with `0.0`, the attachment position is at the start of the site
     - with `1.0`, the attachment position is at the end of the site
     - with `0.5`, the attachment is exactly midway
     [default = 0.5]
     */
    real   site_pos;
    
    /// @}
    
public:

    /// constructor
    DigitProp(const std::string& n) : HandProp(n)  { footprint=0; clear(); }
    
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

