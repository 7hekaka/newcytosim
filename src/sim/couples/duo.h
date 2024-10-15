// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef DUO_H
#define DUO_H

#include "couple.h"
#include "duo_prop.h"

/// A specialized kind of Couple
/**
 The Duo is a couple that can be active or inactive:
 - it is activated instantly inside a given space,
 - is is deactivated spontaneously with the given rate.
 .
 See DuoProp
 
 The DuoLong is automatically selected for non-zero-resting length.
 @ingroup CoupleGroup
 */
class Duo : public Couple
{
    friend class DuoLong;
    
    /// Gillespie countdown timer for deactivation event
    float nextDeact;
    
    /// switch activity 'off'
    void deactivate();
    
    /// check for deactivation
    void tryDeactivate();
    
    /// place in reserve
    void recycle();

protected:
    
    /// active/inactive boolean flag
    int active_;
    
public:
    
    /// constructor
    Duo(DuoProp const*, Vector const& w = Vector(0,0,0));
    
    /// Property
    DuoProp const* prop() const { return static_cast<DuoProp const*>(Couple::prop); }

    /// destructor
    virtual ~Duo();
    
    /// switch activity 'on'
    void activate();

    /// activity flag
    int  active() const { return active_; }
    
    /// simulation step for a free Duo
    void stepFF();
    
    /// simulation step for a Duo attached by Hand1
    void stepAF();
    
    /// simulation step for a Duo attached by Hand2
    void stepFA();
    
    /// simulation step for a linking Duo
    void stepAA();

    /// return unique character identifying the class
    ObjectTag tag() const { return Couple::DUO_TAG; }

    /// write to file
    void write(Outputter&) const;
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
};


#endif

