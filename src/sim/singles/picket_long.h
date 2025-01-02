// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef PICKET_LONG_H
#define PICKET_LONG_H

#include "picket.h"


/// a Picket with a non-zero resting length.
/**
 This single is fixed at its foot position in absolute space.
 It has a non-zero resting length.

 @ingroup SingleGroup
 */
class PicketLong : public Picket
{
    
    /// the side (top/bottom) of the interaction
    mutable Torque mArm;
    
    /// used to recalculate `mArm`
    static Torque calcArm(Interpolation const& pt, Vector& pos, real len);
    
public:

    /// constructor
    PicketLong(SingleProp const*, Vector const& = Vector(0,0,0));

    /// destructor
    ~PicketLong();
    
    /// recalculates `mArm` when making a bridge
    void afterAttachment(Hand const*);

    /// position on the side of fiber used in setInteractions()
    Vector sidePos() const;
    
    /// force = stiffness * ( posFoot() - posHand() )
    Vector force() const;
    
    /// Monte-Carlo step if Hand is attached
    void stepA();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif
