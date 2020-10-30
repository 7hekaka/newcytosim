// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#ifndef TUBULE_H
#define TUBULE_H

#include "assert_macro.h"
#include "tubule_prop.h"
#include "object.h"
#include "buddy.h"

class Glossary;
class Fiber;
class Simul;
class Meca;

/*
 13 fibers arranged in a tubular configurations into a Microtubule
 
 FJN, Cambridge, Sept--Oct 2019
 */
class Tubule : public Object, private Buddy
{
private:
    
    /// number of protofilaments
    static constexpr size_t NFIL = 13;
    static constexpr size_t FILM = NFIL+2;

    /// central backbone if present
    Fiber* bone_;
    
    /// constitutive filaments
    Fiber* fil_[FILM];
    
    /// offset in abscissa
    real   offset_[FILM];
    
    /// the Property of this object
    TubuleProp const* prop;

public:
    
    /// constructor
    Tubule() : prop(nullptr) { reset(); }
    
    /// constructor
    Tubule(TubuleProp * p);

    /// destructor
    ~Tubule();
    
    /// simulation
    void step();

    /// initialize
    void reset();
    
    /// create filaments
    ObjectList build(Glossary&, Simul&);

    /// position of centerline at distance 'dis' from the MINUS_END
    Vector     posCenterlineM(real dis);

    
    /// initialize sister[] and brother[]
    void setFamily(Fiber const*);

    /// update length of family members
    void salute(Buddy const*);
    
    /// handles the disapearance of one of the filament
    void goodbye(Buddy const*);
    
    ///
    void setInteractions(Meca&) const;
    ///
    void setInteractionsB(Meca&) const;
    ///
    void setInteractionsC(Meca&) const;
    
    
    /// a static_cast<> of Object::next()
    Tubule *  next() const { return static_cast<Tubule*>(nextO); }
    
    /// a static_cast<> of Object::prev()
    Tubule *  prev() const { return static_cast<Tubule*>(prevO); }
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 't';

    /// an ASCII character identifying the class of this object
    ObjectTag tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return prop; }

    /// read
    void      read(Inputter&, Simul&, ObjectTag);
    
    /// write
    void      write(Outputter&) const;
    
    /// debug printout
    void      report(std::ostream&);
};

#endif
