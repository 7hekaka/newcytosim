// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef HAND_H
#define HAND_H

#include "fiber_site.h"
#include "hand_prop.h"

class HandMonitor;
class FiberGrid;
class FiberProp;
class Simul;



/// Simulates the stochastic binding/unbinding of a FiberSite
/**
 A Hand is always part of a larger construct, for example Single or Couple.
 
 Hand is the parent to many class that implement different fiber-related activities.
 Hand provides binding/unbinding capacity to these derived classes.
 
 Attachment occurs with constant rate @ref HandPar "binding_rate" to any fiber located
 at distance  @ref HandPar "binding_range" or less.
 If attachment occurs, it happens on the closest point of the fiber,
 which is either the projection of the current position on the fiber axis, 
 or one of the fiber end.
 
 You can restrict binding to happen only near the ends by setting `bind_only_end`,
 and the associated cutoff distance `bind_end_range`.

 Detachment increases exponentially with force:

     off_rate = unbinding_rate * exp( force.norm() / unbinding_force )

 See @ref HandPar
 @ingroup HandGroup
 
 @todo Include FiberSite site as member variable, instead of derivation
 */
class Hand : public FiberSite
{
    /// a monitor that does nothing
    static HandMonitor dummyMonitor;

private:
    
    /// disabled default constructor
    Hand();
    
    /// Pointer used to build the list of Hands bound to a Fiber
    Hand * hNext;
    
    /// Pointer used to build the list of Hands bound to a Fiber
    Hand * hPrev;

protected:

    /// the monitor associated with this Hand
    HandMonitor * hMonitor;
    
    /// Gillespie normalized time for detachment (must be set at attachment)
    float nextDetach;
    
    /// Gillespie countdown timer used for other activities
    float nextAct;

public:
    
    /// Property is constant, so we do not need to make it private
    HandProp const* prop;

    /// Property
    HandProp const* property() const { return prop; }
    
    /// constructor
    /**
     Hand are created in Single and Couple using HandProp::newHand().
     HandProp is parent to several classes, exactly mirroring the hierarchy of Hands.
     This ensures that the correct class is created and associated with the
     correct HandProp and HandMonitor:
     
         HandProp * hp = HandProp::newProperty(name, opt);
         Hand * h = hp->newHand(this);
     */
    Hand(HandProp const*, HandMonitor*);

    /// destructor
    virtual ~Hand();


    /// return next Hand in Fiber's list
    Hand * next() const { return hNext; }
    
    /// return previous Hand in Fiber's list
    Hand * prev() const { return hPrev; }

    /// set next Hand in Fiber's list
    void next(Hand * h) { hNext = h; }
    
    /// set previous Hand in Fiber's list
    void prev(Hand * h) { hPrev = h; }

    
    /// a random position, at distance `binding_range' on the side of the fiber
    Vector posSide() const;
    
    /// move attached Hand to a different fiber, at the given abscissa
    void relocate(Fiber* f, real a);
    
    /// move to a different fiber, at the same abscissa
    void relocate(Fiber* f) { relocate(f, abscissa()); }

    /// relocate to the specified tip of the current fiber
    void moveToEnd(FiberEnd);

    /// bind at position `a` on Fiber `f`
    void locate(Fiber* f, real a);

    // Check that binding can occur on Fiber, from BITWISE AND of the binding keys
    bool keyMatch(Fiber const* fib) const { return prop->binding_key & fib->prop->binding_key; }
    
    /// return Monitor
    HandMonitor const* monitor() const { return hMonitor; }
    
    
    /// tell if attachment at given site is permitted
    virtual bool attachmentAllowed(FiberSite&) const;
    
    /// bin at position represented by FiberSite
    virtual void attach(FiberSite const&);
    
    /// detach
    virtual void detach();

    /// simulate when the Hand is not attached
    virtual void stepUnattached(Simul&, Vector const& pos);

    /// simulate when the Hand is attached but not under load
    virtual void stepUnloaded();

    /// simulate when the Hand is attached and under load
    virtual void stepLoaded(Vector const& force);

    /// this is called when disassembly occured plus end
    virtual void handleDisassemblyM();
    
    /// this is called when the attachment point is below the minus end
    virtual void handleDisassemblyP();

        
    /// check abscissa against fiber edge, and calls handle functions if necessary.
    void checkFiberRange(real absM, real absP);

    /// attach at specified distance `ab` from FiberEnd (this calls attach(FiberSite))
    void attach(Fiber * f, real a, FiberEnd ref) { locate(f, f->abscissaFrom(a, ref)); }
    
    /// attach at the given end of Fiber (this calls attach(FiberSite))
    void attachEnd(Fiber * f, FiberEnd end) { locate(f, f->abscissaEnd(end)); }

    /// detach, without updating Monitor
    void detachHand();

    /// attach at abscissa of given Fiber (this calls attach(FiberSite))
    void attachTo(Fiber * f, real a) { attach(FiberSite(f, a)); }
    
    /// attach at specified distance `ab` from FiberEnd (this calls attach(FiberSite))
    void attachTo(Fiber * f, real a, FiberEnd ref) { attach(FiberSite(f, f->abscissaFrom(a, ref))); }
    
    /// attach at the given end of Fiber (this calls attach(FiberSite))
    void attachToEnd(Fiber * f, FiberEnd end) { attach(FiberSite(f, f->abscissaEnd(end))); }
    
    
    /// return position of other Hand, if part of a Couple, or position of Single
    Vector linkFoot() const;
    
    /// return stiffness of associated link
    real linkStiffness() const;
    
    /// read from file
    ObjectID readHand(Inputter&, Simul&);
    
    /// write to file
    void writeHand(Outputter&) const;
    
    /// reset Gillespie's counters
    void resetTimers();
    
    /**
     Test for spontaneous detachment using Gillespie counter (@ref Stochastic)
     @return true if the hand should detach
     */
    bool checkDetachment()
    {
        assert_true( nextDetach >= 0 );
        nextDetach -= prop->unbinding_rate_dt;
        
        return ( nextDetach < 0 );
    }
    
    /**
     Test for force-dependent detachment using Gillespie counter (@ref Stochastic)
     @return true if the hand should detach
     */
    bool checkKramersDetachment(const real force)
    {
        assert_true( nextDetach >= 0 );
        /*
         Attention: the exponential term can easily become numerically "infinite",
         which is problematic if 'unbinding_rate==0' and 'unbinding_force' is finite.
         This issue is handled in HandProp::complete()
         */
        //std::cerr << prop->name() << " " << std::exp(force*prop->unbinding_force_inv) << "\n";
        nextDetach -= prop->unbinding_rate_dt * std::exp(force*prop->unbinding_force_inv);
        
        return ( nextDetach < 0 );
    }
};

/// output operator
std::ostream& operator << (std::ostream&, Hand const&);

#endif

