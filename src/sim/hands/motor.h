// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MOTOR_H
#define MOTOR_H

#include "hand.h"
#include "motor_prop.h"

/// a Hand that can move smoothly on a Fiber
/**
 The Motor is a Hand, and thus can bind and unbind from Fiber.
 
 A bound Motor can move along its fiber.
 The direction of the motor is set by the sign of @ref MotorPar "unloaded_speed".
 If the speed is positive, the motor attempts to move towards the plus end.
 
 The speed is linearly proportional to the load of the motor.
 The load is the projection of the force vector on the direction of the fiber.
 
     real load = force * direction_of_fiber;
     real speed = unloaded_speed * ( 1 + load / stall_force );
 
 The actual movement depends on the time step:

     real displacement = speed * timestep;
 
 As defined in Hand, detachment increases exponentially with force.

 See Examples and the @ref MotorPar.
 @ingroup HandGroup
 */
class Motor : public Hand
{
private:
    
    /// disabled default constructor
    Motor();

public:

    /// Property
    MotorProp const* prop() const { return static_cast<MotorProp const*>(Hand::prop); }

    /// constructor
    Motor(MotorProp const*, HandMonitor*);
    
    
    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

