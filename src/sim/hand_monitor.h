// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef HAND_MONITOR
#define HAND_MONITOR

#include "real.h"
#include "vector.h"
#include "inventoried.h"

class Hand;
class Simul;
class Single;
class Couple;
class FiberSite;


/// base class to monitor and control Hand's actions
/**
 The HandMonitor defines an interface that is implemented in both Single and Couple.
 It has two functions:
 1- It allows to inform Single and Couple if their Hand bind or unbind.
 2- It is a mechanism for a Hand to access data from the Single or Couple
  to which it belongs.
 .
 */
class HandMonitor
{
public:

    /// Returning `false` prevents the attachment (this is called before every attempt)
    virtual bool permitAttachment(FiberSite const&, Hand const*) const { return true; }
    
    /// called after attachement
    virtual void afterAttachment(Hand const*) {}
    
    /// called before detachment
    virtual void beforeDetachment(Hand const*) {}
    
    
    /// return the Hand that is not the argument, in a Couple
    virtual Hand const* otherHand(Hand const*) const { return nullptr; }
    
    /// return this if `this` is a Single, and `nullptr` otherwise
    virtual Single const* toSingle() const { return nullptr; }

    /// return this if `this` is a Couple, and `nullptr` otherwise
    virtual Couple const* toCouple() const { return nullptr; }

    
    /// return the distal position of the link attached to this Hand
    /** returns the position of the other Hand, if the Hand is part of a Couple */
    virtual Vector linkFoot(Hand const*) const { return Vector(0,0,0); }
    
    /// return the direction associated with the anchoring
    virtual Vector linkDir(Hand const*) const { return Vector::randU(); }

    /// resting length of the link involving this Hand
    virtual real linkRestingLength() const { return 0; }
    
    /// stiffness of the link involving this Hand
    virtual real linkStiffness() const { return 0; }
    
    /// signature of the Solid underlying the Single
    virtual ObjectSignature baseSignature() const { return 0; }

};


#endif
