// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef COUPLE_H
#define COUPLE_H

#include "object.h"
#include "hand_monitor.h"
#include "couple_prop.h"
#include "hand.h"

class Meca;


/// A set of two Hand linked by an elastic element
/**
 A Couple contains two pointers to Hand:
 - cHand1
 - cHand2
 .
 There are 4 possible states for a Couple:
 - state FF (0): cHand1 and cHand2 are free,
 - state AF (1): cHand1 is bound, cHand2 is free,
 - state FA (2): cHand1 is free, cHand2 is bound,
 - state AA (3): both hands are attached
 .
 The method state() return the state of the Couple in [0-3].

 Generally the Couple behaves according to its state:
 - FF     : the Couple is diffusing and both Hands are trying to bind fibers,
 - AF, FA : the localization is given by the attachement point on the fiber,
 - AA     : the Couple is acting as a Hookean spring between the two fibers.
 .
 
 The default Couple has:
 - a zero resting length, using Meca::addLink()
 - no specificity
 .

 The related class CoupleLong with a non-zero resting length is selected automatically

 @ingroup CoupleGroup
 */
class Couple : public Object, public HandMonitor
{
public:

    /// associated properties
    CoupleProp const* prop;
    
protected:
    
    /// position of complex when detached
    Vector cPos;
    
    /// first Hand
    Hand * cHand1;

    /// second Hand
    Hand * cHand2;
    
    /// specialization of HandMonitor
    bool permitAttachment(FiberSite const&, Hand const*) const;
    /// specialization of HandMonitor
    void afterAttachment(Hand const*);
    /// specialization of HandMonitor
    void beforeDetachment(Hand const*);
    
    /// specialization of HandMonitor
    Hand const* otherHand(Hand const*) const;
    /// specialization of HandMonitor
    Couple const* toCouple() const { return this; }

    /// true if both Hands are attached
    bool hasLink() const { return cHand1->attached() && cHand2->attached(); }
    /// stiffness of the interaction, if the Couple is bridging
    real linkStiffness() const { return hasLink() * prop->stiffness; }
    /// specialization of HandMonitor
    Vector linkFoot(Hand const*) const;
    /// specialization of HandMonitor
    real linkRestingLength() const { return prop->length; }

    /// update position to account for diffusion in one time step
    void diffuse();
    
public:
    
    /// constructor
    Couple(CoupleProp const*, Vector const& w = Vector(0,0,0));

    /// destructor
    virtual ~Couple();

    /// copy operator
    Couple & operator = (Couple const&);
    
    //--------------------------------------------------------------------------
    
    /// change the property and update the two Hands (experimental)
    void changeProperty(CoupleProp *);
    
    /// add interactions to a Meca
    virtual void setInteractions(Meca&) const;
    
    /// add interactions to a Meca (experimental)
    virtual void setInteractionsAF(Meca&) const;
    
    /// add interactions to a Meca (experimental)
    virtual void setInteractionsFA(Meca&) const;
    
    //--------------------------------------------------------------------------
    
    /// the position of the complex, calculated from cPos, cHand1 and cHand2
    Vector position() const;
   
    /// Couple can be displaced only if it is not attached
    int mobile() const { return !cHand1->attached() && !cHand2->attached(); }
    
    /// translate object's position by the given vector
    void translate(Vector const& x) { cPos += x; }
    
    /// move object to specified position
    void setPosition(Vector const& x) { cPos = x; }

    /// bring object to centered image using periodic boundary conditions
    void foldPosition(Modulo const*);
    
    /// move to a random position inside the confining space
    void randomizePosition();
    
    //--------------------------------------------------------------------------
    
    /// activity flag
    virtual int active() const { return 1; }
    
    /// turn activity on
    virtual void activate() {}

    /// the state of the Couple in { 0 ... 3 } representing { FF, AF, FA, AA }
    int state() const { return cHand1->attached() + 2 * cHand2->attached(); }
    
    /// true if Hand1 or Hand2 are attached
    int attached() const { return cHand1->attached() || cHand2->attached(); }

    /// category of link: Parallel; Anti-parallel; X; T+; V+; T-; V-
    int configuration(real len, real max_cos=0.5) const;

    /// return one of the Hand that is attached, or zero if both are detached
    Hand * attachedHand() const;
     
    /// cosine of the angle between the two Fibers attached by the hands
    real cosAngle() const { return dot(cHand1->dirFiber(), cHand2->dirFiber()); }
    
    /// the position of the complex if it is unattached
    Vector posFree() const { return cPos; }
    
    /// offset between hands, essentially: ( cHand2->posHand() - cHand1->posHand() )
    Vector stretch() const;

    /// position on the side of fiber1 used in setInteractions()
    virtual Vector sidePos1() const { return cHand1->pos(); }
    
    /// position on the side of fiber2 used in sideInteractions()
    virtual Vector sidePos2() const { return cHand2->pos(); }
    
    /// force between hands, essentially: stiffness * ( sidePos2 - sidePos1 )
    virtual Vector force() const;

    //--------------------------------------------------------------------------

    /// simulation step for a free Couple: diffusion
    virtual void stepFF();
    
    /// simulation step for a Couple attached by Hand1
    virtual void stepAF();
    
    /// simulation step for a Couple attached by Hand2
    virtual void stepFA();
    
    /// simulation step for a doubly-attached Couple
    virtual void stepAA();

    /// simulate movement and detachement of Hand1, skipping attachment
    void stepHand1();
    
    /// simulate movement and detachement of Hand2, skipping attachment
    void stepHand2();

    //--------------------------------------------------------------------------

    /// pointer to Hand1
    Hand * hand1() { return cHand1; }
   
    /// pointer to constant Hand1
    Hand const* hand1() const { return cHand1; }
    
    /// true if Hand1 is attached
    bool attached1() const { return cHand1->attached(); }
    
    /// Fiber to which Hand1 is attached, or zero if not attached
    Fiber const* fiber1() const { return cHand1->fiber(); }
    
    /// attachment position of Hand1 along fiber (only valid if Hand1 is attached)
    real abscissa1() const { return cHand1->abscissa(); }
    
    /// position of Hand1 when attached (only valid if Hand1 is attached)
    Vector posHand1() const { return cHand1->pos(); }
    
    /// direction of Fiber at attachment point of Hand1 (only valid if Hand1 is attached)
    Vector dirFiber1() const { return cHand1->dirFiber(); }
 
    /// attach Hand1 at the given FiberSite
    void attach1(FiberSite s) { if ( cHand1->attachmentAllowed(s) ) cHand1->attach(s); }
    
    /// attach Hand1 at the given end
    void attachEnd1(Fiber* f, FiberEnd end) { cHand1->attachEnd(f, end); }
    
    /// move Hand1 to given end
    void moveToEnd1(FiberEnd end) { cHand1->moveToEnd(end); }

    //--------------------------------------------------------------------------

    /// pointer to Hand2
    Hand * hand2() { return cHand2; }
    
    /// pointer to constant Hand2
    Hand const* hand2() const { return cHand2; }

    /// true if Hand2 is attached
    bool attached2() const { return cHand2->attached(); }
    
    /// Fiber to which Hand2 is attached, or zero if not attached
    Fiber const* fiber2() const { return cHand2->fiber(); }
    
    /// attachment position of Hand2 along fiber (only valid if Hand2 is attached)
    real abscissa2() const { return cHand2->abscissa(); }
    
    /// position of Hand2 when attached (only valid if Hand2 is attached)
    Vector posHand2() const { return cHand2->pos(); }
    
    /// direction of Fiber at attachment point of Hand2 (only valid if Hand2 is attached)
    Vector dirFiber2() const { return cHand2->dirFiber(); }
    
    /// attach Hand2 at the given FiberSite
    void attach2(FiberSite s) { if ( cHand2->attachmentAllowed(s) ) cHand2->attach(s); }
    
    /// attach Hand2 at the given end
    void attachEnd2(Fiber *f, FiberEnd end) { cHand2->attachEnd(f, end); }
    
    /// move Hand2 to given end
    void moveToEnd2(FiberEnd end) { cHand2->moveToEnd(end); }

    //--------------------------------------------------------------------------

    /// set next element
    void next(Couple * x) { next_ = x; }

    /// set previous element
    void prev(Couple * x) { prev_ = x; }

    /// a static_cast<> of Object::next()
    Couple * next() const { return static_cast<Couple*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Couple * prev() const { return static_cast<Couple*>(prev_); }
    
    //------------------------------ read/write --------------------------------

    /// a unique character identifying the class
    static const ObjectTag TAG = 'c';
    
    /// a unique character identifying the Duo
    static const ObjectTag DUO_TAG = 'y';
    
    /// return unique character identifying the class
    ObjectTag tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// write to file
    void write(Outputter&) const;
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// return PointDisp of Hand1
    PointDisp const* disp1() const { return cHand1->disp(); }
    
    /// return PointDisp of Hand2
    PointDisp const* disp2() const { return cHand2->disp(); }
    
    /// return PointDisp of Hand1 if it is visible, otherwise Hand2
    PointDisp const* disp12() const;
    
    /// return PointDisp of Hand2 if it is visible, otherwise Hand1
    PointDisp const* disp21() const;

};


#endif

