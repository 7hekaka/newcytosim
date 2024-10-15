// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef WALKER_PROP_H
#define WALKER_PROP_H

#include "digit_prop.h"

/// Additional Property for Walker
/**
 @ingroup Properties
*/
class WalkerProp : public DigitProp
{
    friend class Walker;
    
public:
    
    /**
     @defgroup WalkerPar Parameters of Walker
     @ingroup Parameters
     Inherits @ref DigitPar.
     @{
     */

    /// force at which stepping rate becomes zero
    real stall_force;
    
    /// speed if `force=0` ( `unloaded_speed = rate * step_size` )
    /**
     A positive value specifies a plus-end directed motor.
     A negative value specifies a minus-end directed motor.
     */
    real unloaded_speed;
    
    /// probability to detach per step
    /**
     This probability in [0,1] is tested for every successful step, creating an
     additional detachment opportunity that is proportional to the distance
     travelled by the motor, in contrast to `unbinding_rate` which gives a
     contribution that is proportional to the duration of the interaction.
     
     In particular, `unbinding_chance` does not lead to the detachment of a stalled motor.
     */
    real unbinding_chance;

    /// @}
    
    
private:
    
    /// derived variable
    float var_rate_dt;
    
    /// derived variable
    float walking_rate_dt;
    
    /// number of fiber's step in one 'step'
    mutable int stride;
    
public:

    /// constructor
    WalkerProp(const std::string& n) : DigitProp(n)  { clear(); }
    
    /// destructor
    ~WalkerProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// perform additional tests for the validity of parameters, given the elasticity of the link
    void checkStiffness(real stiff, real len, real mul, real kT) const;
    
    /// return a carbon copy of object
    Property* clone() const { return new WalkerProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
    /// return 'unload_speed' for the Motor class
    real motorSpeed() const { return unloaded_speed; }

};

#endif

