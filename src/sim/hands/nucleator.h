// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef NUCLEATOR_H
#define NUCLEATOR_H

#include "hand.h"
#include "nucleator_prop.h"

/// A Hand that can nucleate a Fiber
/**
 The Nucleator is a Hand, and thus can bind and unbind from fibers,
 but in most cases, however you would want `binding_rate=0`.
 
 A free Nucleator can create new fibers with a prescibed rate.
 The rate, the type of fiber and the characteristics of the fiber are set
 as three values in property @ref NucleatorPar "nucleate".
 
 If the nucleator is part of a Couple, the parameter "branch_angle"
 can be set to define the direction of the new fiber with respect to the
 preexisting fiber (the one that is bound by the other Hand of the Couple).
 
 By default the nucleator stays attached at the end of the fiber that it has created.
 This can be changed by setting: @ref NucleatorPar "hold_end"

 See Examples and the @ref NucleatorPar.
 @ingroup HandGroup
 
 */
class Nucleator : public Hand
{
private:
    
    /// disabled default constructor
    Nucleator();
    
    /// create a new Fiber
    ObjectList createFiber(Simul&, Vector pos, FiberProp const*, Glossary&);

public:
    
    /// Property
    NucleatorProp const* prop() const { return static_cast<NucleatorProp const*>(Hand::prop); }

    /// constructor
    Nucleator(NucleatorProp const*, HandMonitor*);
    
    
    /// simulate when is not attached
    void stepUnattached(Simul&, Vector const& pos);

    /// simulate when `this` is attached but not under load
    void stepUnloaded();

    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
    /// detach from Fiber
    void detach();

};

#endif

