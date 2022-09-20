// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef CROSSLINK_H
#define CROSSLINK_H

#include "couple.h"
#include "crosslink_prop.h"

/// A specialized kind of Couple
/**
 The Crosslink is a simpler kind of Couple, which does not support `trans_activated`
 
 It has a zero resting length, using Meca::addLink()
 
 The related class CrosslinkLong with a non-zero resting length is selected automatically
 @ingroup CoupleGroup
 */
class Crosslink : public Couple
{
public:

    /// constructor
    Crosslink(CrosslinkProp const*, Vector const& w = Vector(0,0,0));
    
    /// Property
    CrosslinkProp const* prop() const { return static_cast<CrosslinkProp const*>(Couple::prop); }

    /// destructor
    virtual ~Crosslink();
    
    /// simulation step for a free Couple: diffusion
    virtual void stepFF();
    
    /// add interactions to a Meca
    void setInteractions(Meca&) const;

};


#endif

