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
    
    /// force parameter expressing the sensitivity of stepping
    real stepping_force;

    /// rate of forward step
    real forward_rate;
    
    /// backward rate
    real backward_rate;
    
    /// probability to detach per step
    /**
     This probability in [0,1] is tested for every successful step, and thus
     creates an additionaldetachment opportunity that is proportional to the
     distance travelled by the motor, in contrast to `unbinding_rate` which
     gives a contribution that is proportional to the duration of the interaction.
     
     Hence, this does not contributes to the detachment of a stalled motor.
     */
    real unbinding_chance;
    
    /// sign indicates directionality; magnitude the number of steps
    int stepping_stride;
    
    /// @}
    
private:
    
    float forward_rate_dt;
    float backward_rate_dt;
    float force_inv;
    
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

