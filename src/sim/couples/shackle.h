// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SHACKLE_H
#define SHACKLE_H

#include "couple.h"
#include "shackle_prop.h"

/// A specialized kind of Couple
/**
 The Shackle creates a connection that is slipery on Hand1
 it uses Meca::addSlidingLink().
 
 Note: this is highly experimental!
 @ingroup CoupleGroup
 */
class Shackle : public Couple
{
    
public:
    
    /// constructor
    Shackle(ShackleProp const* p, Vector const& w = Vector(0,0,0));
    
    /// Property
    ShackleProp const* prop() const { return static_cast<ShackleProp const*>(Couple::prop); }

    /// destructor
    virtual ~Shackle();
    
    //--------------------------------------------------------------------------
    
    /// simulation step if doubly attached
    void stepAA();
    
    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif

