// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef WRIST_LONG_H
#define WRIST_LONG_H

#include "wrist.h"
#include "mecapoint.h"

/// Switch to use Meca::addLongLink of Meca::addSideLink
#define WRIST_USES_LONGLINK 1


/// a Wrist with a non-zero resting length.
/**
 The anchorage is described by Mecapoint sBase:
 - the Mecable is base_.mecable()
 - the index of the vertex on this Mecable is base_.point()
 .
 It has a non-zero resting length.
 
 @ingroup SingleGroup
 */
class WristLong : public Wrist
{
#if WRIST_USES_LONGLINK
    typedef Vector WristArm;
#else
    typedef Torque WristArm;
#endif
    
    /// the side (top/bottom) of the interaction
    mutable WristArm mArm;

    /// used to calculate `mArm`
    static WristArm calcArm(Interpolation const& pt, Vector const& pos, real len);
    
public:
     
    /// constructor
    WristLong(SingleProp const*, Mecable const*, size_t point);

    /// destructor
    ~WristLong();

    //--------------------------------------------------------------------------
    
    /// position on the side of fiber used for sideInteractions
    Vector sidePos() const;
    
    /// force = stiffness * ( posFoot() - posHand() )
    Vector force() const;
        
    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif
