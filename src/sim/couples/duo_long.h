// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DUO_LONG_H
#define DUO_LONG_H

#include "duo.h"

/// A Duo with a non-zero resting length
/**
 The DuoLong is a couple that can be active or inactive:
 - it is activated instantly inside a given space,
 - is is deactivated spontaneously with the given rate.
 .
 See Duo

 The DuoLong differs from Duo in that it uses a non-zero resting length,
 and creates its interaction using Meca::addSideLink2D() or addSideLink3D()
 
 DuoLong is automatically selected if ( prop:length > 0 )
 @ingroup CoupleGroup
 */
class DuoLong : public Duo
{    
    /// the side (top/bottom) of the interaction
    mutable Torque mArm;
    
    /// used to calculate `mArm`
    static Torque calcArm(Interpolation const& pt, Vector const& pos, real len);
    
public:
    
    /// constructor
    DuoLong(DuoProp const*, Vector const& w = Vector(0,0,0));

    /// destructor
    virtual ~DuoLong();
    
    /// recalculates `mArm` when making a bridge
    void afterAttachment(Hand const*);

    /// position on the side of fiber1 used for sideInteractions
    Vector sidePos1() const;
 
    /// force between hands, essentially: stiffness * ( cHand2->posHand() - cHand1->posHand() )
    Vector force() const;
    
    /// simulation step for a doubly-attached Couple
    void stepAA();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif

