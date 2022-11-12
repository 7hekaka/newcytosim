// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MIGHTY_H
#define MIGHTY_H

#include "hand.h"
#include "mighty_prop.h"

/// A Hand that can move and do other things to a Fiber
/**
 The Mighty is a Hand, and can thus bind and unbind from Fiber.
 
 Mighty is currently a copy of Motor.
 It can be used to implement custom advanced functionalities.
 
 See Examples and the @ref MightyPar.
 @ingroup HandGroup 
 */
class Mighty : public Hand
{
private:
    
    /// disabled default constructor
    Mighty();

public:
    
    /// Property
    MightyProp const* prop() const { return static_cast<MightyProp const*>(Hand::prop); }

    /// constructor
    Mighty(MightyProp const*, HandMonitor*);


    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

