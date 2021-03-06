// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef WRIST_H
#define WRIST_H

#include "single.h"
#include "interpolation4.h"

/// a Single anchored to a Mecable.
/**
 The Wrist is anchored to a Solid, on a position that is interpolated from the
 Solid's vertices. See class Interpolation4

 @ingroup SingleGroup
 */
class Wrist : public Single
{
protected:
    
    Interpolation4 base_;
    
public:
    
    /// Construct object anchored at one Mecapoint
    Wrist(SingleProp const*, Mecable const*, size_t point);
    
    /// Construct object anchored between two Mecapoint
    Wrist(SingleProp const*, Mecable const*, size_t, size_t, real);
   
    /// Constructor object interpolated over a triad of Mecapoint
    Wrist(SingleProp const*, Mecable const*, size_t ref, Vector pos);

    /// destructor
    ~Wrist();
    
    //--------------------------------------------------------------------------
    
    /// return the position in space of the object
    Vector  position() const { return base_.position(); }
    
    /// Wrist accepts translation
    int     mobile() const { return 0; }
    
    /// translate object's position by the given vector
    void    translate(Vector const&) { }
    
    /// bring object to centered image using periodic boundary conditions
    void    foldPosition(Modulo const*) { }

    //--------------------------------------------------------------------------
    
    /// Object to which this is anchored
    Mecable const* base() const { return base_.mecable(); }
    
    /// attach at one Mecapoint
    void    rebase(Mecable const* mec, size_t pti) { base_.set(mec, pti); }
    
    /// attach between two Mecapoint
    void    rebase(Mecable const* mec, size_t a, size_t b, real c) { base_.set(mec, a, b, c); }
    
    /// attach over a triad of Mecapoint
    void    rebase(Mecable const* mec, size_t ref, Vector pos) { base_.set(mec, ref, pos); }
    
    
    /// stiffness of the interaction
    real    linkStiffness() const { return prop->stiffness; }

    /// the position of the anchoring point
    Vector  posFoot() const { return base_.position(); }

    /// true if Single creates a link
    bool    hasFoot() const { return true; }
    
    /// stretch of the link = ( posFoot() - posHand() )
    Vector  stretch() const;

    /// force = stiffness * ( posFoot() - posHand() )
    Vector  force() const;
    
    
    /// Monte-Carlo step for a free Single
    void    stepF();
    
    /// Monte-Carlo step for a bound Single
    void    stepA();

    /// add interactions to a Meca
    void    setInteractions(Meca&) const;

    //--------------------------------------------------------------------------
    
    /// the Wrist uses a specific TAG to distinguish itself from the Single
    static const ObjectTag TAG = 'w';
    
    /// return unique character identifying the class
    ObjectTag    tag() const { return TAG; }
    
    /// read from file
    void    read(Inputter&, Simul&, ObjectTag);
    
    /// write to file
    void    write(Outputter&) const;
    
};


#endif
