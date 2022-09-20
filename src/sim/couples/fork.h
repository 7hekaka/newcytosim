// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FORK_H
#define FORK_H

#include "couple.h"
#include "fork_prop.h"

/// A specialized kind of Couple
/**
 The fork connect two fibers at an angle.
 It creates a torque between the two fibers, with a resting angle defined by 'ForkProp::angle',
 and a rotational stiffness which is 'ForkProp::angular_stiffness'.
 
 Note that the Fork is unfinished and should not be used.
 
 @ingroup CoupleGroup
 */
class Fork : public Couple
{
#if ( DIM == 2 )
    /// sine of angle, with [up, down] sign in space for 2D
    mutable real sine;
#endif
public:
    
    /// constructor
    Fork(ForkProp const*, Vector const& w = Vector(0,0,0));
    
    /// Property
    ForkProp const* prop() const { return static_cast<ForkProp const*>(Couple::prop); }

    /// destructor
    virtual ~Fork();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;

};


#endif

