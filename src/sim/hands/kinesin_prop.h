// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef KINESIN_PROP_H
#define KINESIN_PROP_H

#include "digit_prop.h"


/// Additional Property for Kinesin
/**
 @ingroup Properties
*/
class KinesinProp : public DigitProp
{
    friend class Kinesin;
    
public:
    
    /**
     @defgroup KinesinPar Parameters of Kinesin
     @ingroup Parameters
     Inherits @ref DigitPar.
     @{
     */
    
    ///
    real    force;

    ///
    real    forward_rate;
    
    /// backward rate
    real    backward_rate;
    
    ///
    real    unbinding_chance;
    
    /// directionality ( -1 / +1 )
    int     stride;
    
    /// @}
    
private:
    
    real    forward_rate_dt;
    real    backward_rate_dt;
    real    force_inv;
    
public:

    /// constructor
    KinesinProp(const std::string& n) : DigitProp(n)  { clear(); }
    
    /// destructor
    ~KinesinProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new KinesinProp(*this); }
    
    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

