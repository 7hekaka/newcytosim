// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DIGIT_H
#define DIGIT_H

#include "hand.h"
class DigitProp;


/// A Hand that bind to discrete sites along a Fiber
/**
 The Digit is a Hand that can only bind at discrete positions along a Fiber,
 corresponding to the Lattice associated with this Fiber.
 
 Binding will be limited by the occupancy stored in the Lattice:
 - the digit will bind only at a vacant lattice site,
 - upon binding, the digit will occupy the site entirely
 .
 
 See Examples and the @ref DigitPar.
 @ingroup HandGroup 
 
 As defined in Hand, detachment increases exponentially with force.
 */
class Digit : public Hand
{
    /// disabled default constructor
    Digit();
    
public:
    
    /// Property
    DigitProp const* prop;
    
    /// constructor
    Digit(DigitProp const*, HandMonitor*);
    
    /// destructor
    ~Digit() {}
    
    /// check if attachement is possible according to properties
    bool   attachmentAllowed(FiberSite& site) const;

    /// attach and update variables
    void   attach(FiberSite const& site);
    
    /// detach
    void   detach();

    
    /// relocate without checking intermediate sites
    int    jumpTo(site_t npos);
    
    /// relocate without checking intermediate sites
    int    jumpToEndM();
    
    /// relocate without checking intermediate sites
    int    jumpToEndP();

    
    /// attempt one step towards the PLUS_END
    int    stepP();
    
    /// attempt one step towards the MINUS_END
    int    stepM();

    
    /// attempt one step of size `s` towards the PLUS_END
    int    jumpP(int s);
    
    /// attempt one step of size `s` towards the MINUS_END
    int    jumpM(int s);

    
    /// attempt `n` steps towards the PLUS_END, checking all intermediate sites
    int    crawlP(int n);
    
    /// attempt `n` steps towards the MINUS_END, checking all intermediate sites
    int    crawlM(int n);

    
    /// simulate when `this` is attached but not under load
    void   stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void   stepLoaded(Vector const& force, real force_norm);
 
    
    /// this is called when the attachment point is beyond the PLUS_END
    void   handleDisassemblyM();
    
    /// this is called when the attachment point is below the MINUS_END
    void   handleDisassemblyP();
    
};

#endif

