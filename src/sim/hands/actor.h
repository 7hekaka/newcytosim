// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ACTOR_H
#define ACTOR_H

#include "hand.h"
#include "actor_prop.h"

/// A Hand that can act on a Fiber
/**
 The Actor is a Hand, and can thus bind and unbind from Fiber.
 
 The Actor currently does nothing else.
 It exists as a template, or as a class that can be used when new functionalities are needed

 See Examples and the @ref ActorPar.
 @ingroup HandGroup 
 */
class Actor : public Hand
{
private:
    
    /// disabled default constructor
    Actor();

public:
    
    /// Property
    ActorProp const* prop() const { return static_cast<ActorProp const*>(Hand::prop); }

    /// constructor
    Actor(ActorProp const*, HandMonitor*);
    
    
    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

