// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef WRIST_H
#define WRIST_H

#include "single.h"
#include "interpolation4.h"

/// a Single anchored to a Mecable.
/**
 The Wrist is anchored to a Solid, on a position that is interpolated from the
 Solid's vertices. See class Interpolation4
 
 The Wrist can be distinguished from other Solid using the `base()` function:
 `Wrist::base()` returns the Mecable onto which the Wrist is based.
 Otherwise, `Single::base()` returns nullptr

 @ingroup SingleGroup
 */
class Wrist : public Single
{
protected:
    
    Interpolation4 base_;
    
public:
    
    /// Construct object anchored at one Mecapoint
    Wrist(SingleProp const*, Mecable const*, unsigned point);

    /// destructor
    ~Wrist();
    
    //--------------------------------------------------------------------------
    
    /// return the position in space of the object
    Vector position() const { return base_.pos(); }
    
    /// Wrist accepts translation
    int mobile() const { return 0; }
    
    /// translate object's position by the given vector
    void translate(Vector const&) { }
    
    /// bring object to centered image using periodic boundary conditions
    void foldPosition(Modulo const*) { }

    //--------------------------------------------------------------------------
    
    /// Object to which this is anchored
    Mecable const* base() const { return base_.mecable(); }

    /// attach at Mecapoint of specified index
    void rebase(Mecable const* mec, unsigned pti) { base_.set(mec, pti); }
    
    /// attach between two Mecapoints of indices `a` and `b`
    void rebase(Mecable const* mec, unsigned a, unsigned b, real c) { base_.set(mec, a, b, c); }
    
    /// attach over a triad of Mecapoints starting at `ref`
    void rebase(Mecable const* mec, unsigned ref, Vector pos) { base_.set(mec, ref, pos); }
    
    
    /// signature of the Solid underlying the Single
    ObjectSignature baseSignature() const { return base()->signature(); }

    /// true if Single creates a link
    bool hasLink() const { return true; }

    /// stiffness of the interaction
    real linkStiffness() const { return prop->stiffness; }

    /// the position of the anchoring point
    Vector posFoot() const { return base_.pos(); }
    
    /// a radial direction at the anchoring point
    Vector dirFoot() const { return base_.normal(); }

    /// stretch of the link = ( posFoot() - posHand() )
    Vector stretch() const;

    /// force = stiffness * stretch()
    Vector force() const;
    
    
    /// Monte-Carlo step if Hand is detached
    void stepF();
    
    /// Monte-Carlo step if Hand is attached
    void stepA();

    /// add interactions to a Meca
    void setInteractions(Meca&) const;

    //--------------------------------------------------------------------------
    
    /// return unique character identifying the class
    ObjectTag tag() const { return WRIST_TAG; }
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to file
    void write(Outputter&) const;
    
    /// export anchoring data to file
    void writeBase(Outputter& out) const { return base_.write(out); }
    
    /// import anchoring data from file
    void readBase(Inputter& in, Simul& sim) { return base_.read(in, sim); }

    /// check validity of base_
    int invalid() const { return base_.invalid(); }
};


#endif
