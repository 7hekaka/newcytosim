// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef CROSSLINK_LONG_H
#define CROSSLINK_LONG_H

#include "crosslink.h"

/// A Crosslink with a non-zero resting length
/**
 The CrosslinkLong adds a non-zero resting length to Crosslink,
 it uses Meca::addSideSideLink() which is a symmetric interaction with non-zero
 resting length.
 
 CrosslinkLong is automatically selected if ( prop:length > 0 )
 @ingroup CoupleGroup
 */
class CrosslinkLong : public Crosslink
{
    /// the side (top/bottom) of the interaction
    mutable Torque mArm1;

    /// the side (top/bottom) of the interaction
    mutable Torque mArm2;
    
public:
    
    /// constructor
    CrosslinkLong(CrosslinkProp const*, Vector const& w = Vector(0,0,0));

    /// destructor
    virtual ~CrosslinkLong();
    
    /// position on the side of fiber1 used for sideInteractions
    Vector sidePos1() const;
    
    /// position on the side of fiber2 used for sideInteractions
    Vector sidePos2() const;

    /// force between hands, essentially: stiffness * ( cHand2->posHand() - cHand1->posHand() )
    Vector force() const;
    
    /// simulation step for a doubly-attached Couple
    void stepAA();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif

