// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#ifndef TUBULE_H
#define TUBULE_H

#include "assert_macro.h"
#include "tubule_prop.h"
#include "glossary.h"
#include "object.h"
#include "real.h"

class Fiber;
class Simul;
class Meca;


class Tubule : public Object
{
private:
    
    /// number of protofilaments
    constexpr static size_t NFIL = 13;
    
    /// initial radius of tubule
    constexpr static real tube_radius = 0.010;
    
    /// distance between protofilaments
    constexpr static real fil_offset = 0.0045;

    
    /// constitutive filaments
    Fiber* fil_[NFIL+2];
    
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
    
    /// create filaments
    ObjectList build(Glossary&, Simul&);
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 't';

    /// an ASCII character identifying the class of this object
    ObjectTag       tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return prop; }

    ///
    void step(Simul&);
    
    ///
    void setInteractions(Meca&);
    
    
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
