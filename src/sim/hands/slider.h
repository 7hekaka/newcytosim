// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SLIDER_H
#define SLIDER_H

#include "hand.h"
#include "slider_prop.h"

/// A Hand that will move on a Fiber passively with viscous resistance
/**
 The Slider is a Hand, and can thus bind and unbind from Fiber.
 
 A bound Slider will move on the fiber it is pulled by external force.
 The Slider is a passive element that may move which staying attached to the fiber.
 
 The @ref SliderPar "mobility" defines the ratio between speed and force
 
     real load = force * direction_of_fiber;
     real displacement = load * mobility * time_step;
 
 See Examples and the @ref SliderPar.
 @ingroup HandGroup 
 */
class Slider : public Hand
{
private:
    
    /// disabled default constructor
    Slider();
    
public:
    
    /// Property
    SliderProp const* prop() const { return static_cast<SliderProp const*>(Hand::prop); }

    /// constructor
    Slider(SliderProp const*, HandMonitor*);

    
    /// simulate when `this` is attached but not under load
    void stepUnloaded();

    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

