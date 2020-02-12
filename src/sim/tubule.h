// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#ifndef TUBULE_H
#define TUBULE_H

#include "assert_macro.h"
#include "tubule_prop.h"
#include "glossary.h"
#include "object.h"
#include "buddy.h"
#include "real.h"

class Fiber;
class Simul;
class Meca;

/**
 13 fibers arranged in a tubular configurations into a Microtubule
 
 FJN, Cambridge, Sept--Oct 2019
 */
class Tubule : public Object, private Buddy
{
private:
    
    /// number of protofilaments
    static constexpr size_t NFIL = 13;
    static constexpr size_t FILM = NFIL+2;
    
    /// initial radius of tubule
    /* defines the distance between the centerline of the protofilaments
     This gives a diameter of 25+2nm, where the 2 represents the position
     of the center of molecules that attach to the tubulin lattice
     */
    static constexpr real tube_radius = 0.0135;

    /// central backbone if present
    Fiber* bone_;
    
    /// constitutive filaments
    Fiber* fil_[FILM];
    
    /// the Property of this object
    TubuleProp const* prop;

public:
    
    /// constructor
    Tubule() : prop(nullptr) { reset(); }
    
    /// constructor
    Tubule(TubuleProp * p);

    /// destructor
    ~Tubule();
    
    /// initialize
    void reset();
    
    /// initialize sister[] and brother[]
    void setFamily(Fiber const*);

    /// create filaments
    ObjectList build(Glossary&, Simul&);
    
    /// handles the disapearance of one of the filament
    void       goodbye(Buddy *);

    
    /// a unique character identifying the class
    static const ObjectTag TAG = 't';

    /// an ASCII character identifying the class of this object
    ObjectTag tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return prop; }

    ///
    void step(Simul&);
    
    ///
    void setInteractions(Meca&);
    ///
    void setInteractionsB(Meca&);
    ///
    void setInteractionsC(Meca&);
    
    
    /// a static_cast<> of Node::next()
    Tubule *  next()  const  { return static_cast<Tubule*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Tubule *  prev()  const  { return static_cast<Tubule*>(nPrev); }
    

    /// read
    void      read(Inputter&, Simul&, ObjectTag);
    
    /// write
    void      write(Outputter&) const;
};

#endif
