// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef COUPLE_LONG_H
#define COUPLE_LONG_H

#include "couple.h"

/// A Couple with a non-zero resting length
/**
 The CoupleLong adds a non-zero resting length to Couple,
 it uses Meca::addSideLink2D() or Meca::addSideLink3D()

 CoupleLong is automatically selected if ( prop:length > 0 )
 @ingroup CoupleGroup
 */
class CoupleLong : public Couple
{
    /// the side (top/bottom) of the interaction
    mutable Torque mArm;
    
    /// used to calculate `mArm`
    static Torque calcArm(Interpolation const& pt, Vector const& pos, real len);
    
public:
    
    /// constructor
    CoupleLong(CoupleProp const*, Vector const& w = Vector(0,0,0));

    /// destructor
    virtual ~CoupleLong();
    
    /// recalculates `mArm` when making a bridge
    void afterAttachment(Hand const*);

    /// position on the side of fiber1 to which the link is made
    Vector sidePos1() const;

    /// force between hands, essentially: stiffness * ( cHand2->pos() - sidePos1() )
    Vector force() const;
    
    /// simulation step for a doubly-attached Couple
    void stepAA();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif

