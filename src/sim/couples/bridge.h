// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef BRIDGE_H
#define BRIDGE_H

#include "couple.h"
#include "bridge_prop.h"

/// A Couple with a different mechanical link
/**
 The Bridge differs from CoupleLong in the nature of the mechanical link that it
 creates betwen two filaments.
 The Bridge uses Meca::addLongLink(), in contrast to other sorts of Long Couples
 
 The "Long link" is a finite resting length Hookean spring, which can freely rotate
 at both of its ends. Hence the angle with respect to the filament is uncontrained,
 unlike the "Side link", in which the spring extends orthogonally to  the direction of the filaments.
 
 Because of this the Bridge does not impose a strict separation between a pair of filaments.
 Longitudinal shear on two filaments connected by 'bridges' will likely affect the distance between them.
 
 \image html meca_links.png
 
 The Bridge should have a non-zero resting length.
 For zero-resting length, use Couple or Crosslink
 @ingroup CoupleGroup
 */
class Bridge : public Couple
{
public:
    
    /// constructor
    Bridge(BridgeProp const*, Vector const& w = Vector(0,0,0));
    
    /// Property
    BridgeProp const* prop() const { return static_cast<BridgeProp const*>(Couple::prop); }

    /// destructor
    virtual ~Bridge();
    
    //--------------------------------------------------------------------------
 
    /// force between hands
    Vector force() const;
    
    /// simulation step for a doubly-attached Couple
    void stepAA();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif

