// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef SINGLE_H
#define SINGLE_H

#include "dim.h"
#include "vector.h"
#include "movable.h"
#include "object.h"
#include "hand_monitor.h"
#include "mecapoint.h"
#include "single_prop.h"
#include "hand.h"


class Space;
class Fiber;
class PointDisp;


/// A point-like object containing one Hand.
/**
 A Single contains one pointer to Hand, and consequently
 inherit the 2 possible states: `attached` or `free`.
 
 By default:
 - Free Single are diffusing, and try to bind to nearby Fibers,
 - Attached Singles are moving along the Fiber to which their Hand is attached.
 .
 
 However, two derived classes change this behavior:
 -# a Picket is fixed in position and does not diffuse,
 -# a Wrist is attached to one vertex of a Mecable.
 .
 
 Wrist and Picket may exert force on the Fiber to which their Hand attaches.
 For WristLong and PicketLong, this force can have a non-zero resting length.
 For these class, `hasLink()` returns true. The force may still be zero if the
 link stiffness is zero, which is the default value.

 Wrist and Picket can be distinguished with Single::base():
 - for Single and Picket, `base() == nullptr`,
 - for Wrist, `base()` returns the Mecable on which the Wrist is attached.
 .

 @ingroup SingleGroup
 */

class Single : public Object, public HandMonitor
{
private:
    
    /// specialization of HandMonitor
    Vector linkFoot(Hand const*) const { return posFoot(); }
    /// specialization of HandMonitor
    Vector linkDir(Hand const*) const { return dirFoot(); }
    /// specialization of HandMonitor
    real linkRestingLength() const { return prop->length; }
    /// stiffness of the interaction
    real linkStiffness() const { return 0; }

protected:
    
    /// specialization of HandMonitor
    void afterAttachment(Hand const*);
    /// specialization of HandMonitor
    void beforeDetachment(Hand const*);
    /// specialization of HandMonitor
    Single const* toSingle() const { return this; }

    /// the position of the foot
    Vector sPos;
    
    /// the motor domain
    Hand * sHand;

public:
    
    /// property
    SingleProp const* prop;

    /// constructor at specified position
    Single(SingleProp const*, Vector const& = Vector(0,0,0));

    /// destructor
    virtual ~Single();
    
    //--------------------------------------------------------------------------
    
    /// associated Hand
    Hand * hand() const { return sHand; }
    
    /// sHand->attached()
    bool attached() const { return sHand->attached(); }
    
    /// sHand->attached()
    int state() const { return sHand->attached(); }

    /// Fiber to which this is attached
    Fiber const* fiber() const { return sHand->fiber(); }
    
    /// attachment position of Hand along fiber (call is invalid if Hand is not attached)
    real abscissa() const { return sHand->abscissa(); }
    
    /// position of the Hand (call is invalid if Hand is not attached)
    Vector posHand() const { return sHand->pos(); }
    
    /// direction of Fiber at attachment point (call is invalid if Hand is not attached)
    Vector dirFiber() const { return sHand->dirFiber(); }
    
    /// attach Hand at the given site
    void attach(FiberSite s) { if ( sHand->attachmentAllowed(s) ) sHand->attach(s); }
    
    /// attach Hand at given Fiber end
    void attachEnd(Fiber * f, FiberEnd e) { sHand->attachEnd(f, e); }

    /// move Hand at given end
    void moveToEnd(FiberEnd e) { sHand->moveToEnd(e); }
    
    /// detach
    void detach() { sHand->detach(); }

    //--------------------------------------------------------------------------
    
    ///return the position in space of the object
    Vector position() const;
    
    /// Single can be translated only if it is not attached
    int mobile() const { return !sHand->attached(); }
    
    /// translate object's position by the given vector
    void translate(Vector const& x)   { sPos += x; }
    
    /// move object to specified position
    void setPosition(Vector const& x) { sPos = x; }

    /// bring object to centered image using periodic boundary conditions
    void foldPosition(Modulo const*);
    
    /// move to a random position inside the confining space
    void randomizePosition();

    //--------------------------------------------------------------------------
    
    /// the position of the anchoring point
    virtual Vector posFoot() const { return sPos; }
    
    /// the direction at the anchoring point
    virtual Vector dirFoot() const { return Vector::randU(); }

    /// position on the side of fiber used for sideInteractions
    virtual Vector sidePos() const { return sHand->pos(); }
    

    /// the Mecable to which this is anchored, or zero
    virtual Mecable const* base() const { return nullptr; }

    /// detach from Mecable to which this is anchored
    virtual void unbase() {}
    
    /// true if Single creates an interaction
    virtual bool hasLink() const { return false; }
    
    /// stretch = ( position_anchor - position_hand ), or zero for a diffusible Single
    virtual Vector stretch() const { return Vector(0,0,0); }

    /// force = stiffness * ( position_anchor - position_hand ), or zero for a diffusible Single
    virtual Vector force() const { return Vector(0,0,0); }

    /// Monte-Carlo step if Hand is detached
    virtual void stepF();
    
    /// Monte-Carlo step if Hand is attached
    virtual void stepA();

    /// add interactions to a Meca
    virtual void setInteractions(Meca&) const;
    
    //--------------------------------------------------------------------------

    /// set next element
    void next(Single * x) { next_ = x; }

    /// set previous element
    void prev(Single * x) { prev_ = x; }

    /// a static_cast<> of Object::next()
    Single * next() const { return static_cast<Single*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Single * prev() const { return static_cast<Single*>(prev_); }

    //--------------------------------------------------------------------------

    /// a unique character identifying the class
    static const ObjectTag TAG = 's';
   
    /// a unique character identifying the derived class Wrist
    static const ObjectTag WRIST_TAG = 'w';

    /// return unique character identifying the class
    virtual ObjectTag tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// return PointDisp of associated Hand
    PointDisp const* disp() const { return sHand->property()->disp; }

    /// return Property::confine_space
    Space const* confineSpace() const { return prop->confine_space; }
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to file
    void write(Outputter&) const;
    
    /// export anchoring data to file
    virtual void writeBase(Outputter& out) const { ABORT_NOW("writeBase == 0"); }
    
    /// import anchoring data from file
    virtual void readBase(Inputter& out, Simul& sim) { ABORT_NOW("readBase == 0"); }

    /// check validity
    int invalid() const { return !sPos.valid(); }
};


#endif
