// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef KINESIN_H
#define KINESIN_H

#include "digit.h"
#include "kinesin_prop.h"


/// A model of the kinesin motor with discrete stepping
/**
 This is a simple Kinesin stochastic model directly inspired from measurements in:
 
     Mechanics of the kinesin step, Carter & Cross, 2005
     http://www.doi.org/10.1038/nature03528
     
 The Kinesin makes discrete jumps along the fiber.
 
 Stepping is stochastic and dependent on force:
    - antagonistic force decreases the rate of forward steps.
    - antagonistic force increases the rate of backward steps.
 .
 
 See Examples and the @ref KinesinPar.
 @ingroup HandGroup
 FJN 19.02.2020
*/
class Kinesin : public Digit
{
private:
    
    /// disabled default constructor
    Kinesin();
    
    /// Gillespie countdown timer for stepping
    float nextBack;
    
public:
    
    /// Property
    KinesinProp const* prop() const { return static_cast<KinesinProp const*>(Hand::prop); }

    /// constructor
    Kinesin(KinesinProp const*, HandMonitor*);
    
    
    /// attach and update variables
    void attach(FiberSite const&);

    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

