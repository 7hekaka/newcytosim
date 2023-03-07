// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PICKET_H
#define PICKET_H

#include "single.h"


/// a Single attached at a fixed position.
/**
 This Single is fixed at its foot position in absolute space.
 An elastic link is created when the Hand is attached.

 @ingroup SingleGroup
 */
class Picket : public Single
{
public:
    
    /// sPos should never change
    void beforeDetachment(Hand const*);
    
    /// stiffness of the interaction
    real linkStiffness() const { return prop->stiffness; }
    
    /// true if Single creates an interaction
    bool hasLink() const { return true; }

public:

    /// constructor
    Picket(SingleProp const*, Vector const& = Vector(0,0,0));

    /// destructor
    ~Picket();
    
    ///return the position in space of the object
    Vector position() const { return sPos; }

    /// Picket accepts translation
    int mobile() const { return 1; }
    
    /// translate object's position by the given vector
    void translate(Vector const& w)  { sPos += w; }
    
    /// move object to specified position
    void setPosition(Vector const& x) { sPos = x; }
    
    /// stretch of the link = ( posFoot() - posHand() )
    Vector stretch() const;

    /// tension in the link = stiffness * stretch()
    Vector force() const;

    /// Monte-Carlo step if Hand is detached
    void stepF();
    
    /// Monte-Carlo step if Hand is attached
    void stepA();
    
    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
};


#endif
