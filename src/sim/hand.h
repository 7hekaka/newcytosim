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
 The class Hand provides binding/unbinding capacity to Fiber.
 It is the parent to many classes that implement different fiber-related activities, such as motors, cutters, etc..
 A Hand is always part of a larger construct, particularly Single or Couple.

 Attachment occurs with constant rate @ref HandPar "binding_rate" to any fiber located
 at distance  @ref HandPar "binding_range" or less.
 If attachment occurs, it happens on the closest point of the fiber,
 which is either the projection of the current position on the fiber axis, 
 or one of the fiber end.
 
 You can restrict binding to occur only near the ends by setting `bind_only_end`,
 specifying which end, and the associated cutoff distance `bind_end_range`.

 Detachment increases exponentially with force:

     off_rate = unbinding_rate * exp( force.norm() / unbinding_force )

 See @ref HandPar
 @ingroup HandGroup
 But the parameter `unbinding_force` is by default `+inf`, meaning that unbinding is not force-sensitive.
 
 @todo Include FiberSite site as member variable, instead of derivation
 */
class Hand : public FiberSite
{
    /// to expose 
    friend class BlinkSortJob<Hand>;
    
    /// a monitor that does nothing
    static HandMonitor dummyMonitor;

private:
    
    /// disabled default constructor
    Hand();
   
protected:
    
    /// Pointer used to build the list of Hands bound to a Fiber
    Hand * next_;
    
    /// Pointer used to build the list of Hands bound to a Fiber
    Hand * prev_;

protected:
    
    /// Property is constant
    HandProp const* prop;

    /// the monitor associated with this Hand
    HandMonitor * hMonitor;
    
    /// Gillespie normalized time for detachment (must be set at attachment)
    float nextDetach;
    
    /// Gillespie countdown timer used for other activities
    float nextAct;

    /// bind at position `a` on Fiber `f`
    void do_attach(Fiber const* f, real a);

public:

    /// unmodifiable Property
    HandProp const* property() const { return prop; }
    
    /// return unmodifiable PointDisp of HandProp
    PointDisp const* disp() const { return prop->disp; }

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
    Hand * next() const { return next_; }
    
    /// return previous Hand in Fiber's list
    Hand * prev() const { return prev_; }

    /// set next Hand in Fiber's list
    void next(Hand * h) { next_ = h; }
    
    /// set previous Hand in Fiber's list
    void prev(Hand * h) { prev_ = h; }

    
    /// a random position, at distance `binding_range' on the side of the fiber
    Vector unbindingPosition() const;
    
    /// move attached Hand to same or different fiber, at the given abscissa
    void relocate(Fiber const* f, real a);

    /// relocate to the specified tip of the current fiber
    void moveToEnd(FiberEnd);

    // Check that binding can occur on Fiber, from BITWISE AND of the binding keys
    bool keyMatch(Fiber const* fib) const { return prop->binding_key & fib->prop->binding_key; }
    
    /// return Monitor
    HandMonitor const* monitor() const { return hMonitor; }
    
    /// the other hand, if this is part of a Couple, or nullptr
    Hand const* otherHand() const;

    
    /// tell if attachment at given site is permitted
    virtual bool attachmentAllowed(FiberSite&) const;
    
    /// bind at position indicated by FiberSite
    virtual void attach(FiberSite const&);
    
    /// detach
    virtual void detach();

    /// simulate when Hand is not attached
    virtual void stepUnattached(Simul&, Vector const& pos);

    /// simulate when Hand is attached but not under load
    virtual void stepUnloaded();

    /// simulate when Hand is attached and under load
    virtual void stepLoaded(Vector const& force);

    /// this is called when disassembly occured plus end
    virtual void handleDisassemblyM();
    
    /// this is called when the attachment point is below the minus end
    virtual void handleDisassemblyP();


    /// check abscissa against fiber edge, and calls handle functions if necessary.
    void checkFiberRange(real absM, real absP);

    /// attach at specified distance `ab` from FiberEnd (this calls attach(FiberSite))
    void attach(Fiber const* f, real a, FiberEnd ref) { do_attach(f, f->abscissaFrom(a, ref)); }
    
    /// attach at the specified end of given Fiber
    void attachEnd(Fiber const* f, FiberEnd end) { do_attach(f, f->abscissaEnd(end)); }

    /// detach, without updating the Monitor
    void detachHand();

    /// attach at abscissa of given Fiber (this calls attach(FiberSite))
    void attachTo(Fiber const* f, real a) { attach(FiberSite(f, a)); }
    
    /// attach at specified distance `ab` from FiberEnd (this calls attach(FiberSite))
    void attachTo(Fiber const* f, real a, FiberEnd ref) { attach(FiberSite(f, f->abscissaFrom(a, ref))); }
    
    /// attach at the given end of Fiber (this calls attach(FiberSite))
    void attachToEnd(Fiber const* f, FiberEnd end) { attach(FiberSite(f, f->abscissaEnd(end))); }
    

#if FIBER_HAS_LATTICE > 0
    
    /// return current Lattice's site bits corresponding to the hand's footprint
    lati_t valLattice(lati_t s) const { return (hLattice->data(s) & prop->footprint); }

    /// return given lattice site bits corresponding to the hand's footprint
    lati_t valLattice(FiberLattice const* lat, lati_t s) const { return (lat->data(s) & prop->footprint); }

    /// flip footprint bits on current site
    void incLattice() const { assert_false(valLattice(hSite)); hLattice->data(hSite) ^= prop->footprint; }

    /// flip footprint bits on current site
    void decLattice() const { hLattice->data(hSite) ^= prop->footprint; assert_false(valLattice(hSite)); }

    /// set FiberSite at index `s` with an abscissa `off` within the site
    void hopLattice(lati_t s)
    {
        assert_true(attached());
        assert_true(hLattice);
        //std::clog << this << ":hop " << hSite << " ---> " << s << "\n";
        decLattice();
        hSite = s;
        incLattice();
        hAbs = s * hLattice->unit() + prop->site_shift;
    }

#elif FIBER_HAS_LATTICE < 0

    /// return current lattice value at given site
    auto valLattice(lati_t s) const { return hLattice->data(s); }

    /// return given lattice value at given site
    auto valLattice(FiberLattice const* lat, lati_t s) const { return lat->data(s); }

    /// add 1.0 to Lattice's site
    void incLattice() const { hLattice->data(hSite) += prop->footprint; }

    /// subtract 1.0 to Lattice's site
    void decLattice() const { hLattice->data(hSite) -= prop->footprint; }
    
    /// set FiberSite at index `s` with an abscissa `off` within the site
    void hopLattice(lati_t s)
    {
        assert_true(attached());
        assert_true(hLattice);
        //std::clog << this << ":hop " << hSite << " ---> " << s << "\n";
        decLattice();
        hSite = s;
        incLattice();
        hAbs = s * hLattice->unit() + prop->site_shift;
    }

#else

    bool valLattice(lati_t) const { return false; }
    void incLattice() const {}
    void decLattice() const {}
    
    /// set FiberSite at index `s` with an abscissa `off` within the site
    void hopLattice(lati_t s)
    {
        hAbs = s * prop->step_size + prop->site_shift;
    }

#endif

    /// return position of other Hand, if part of a Couple, or position of Single
    Vector linkFoot() const;
    
    /// return stiffness of associated link
    real linkStiffness() const;
    
    /// read from file
    ObjectID readHand(Inputter&, Simul&);
    
    /// write to file
    void writeHand(Outputter& out) const { FiberSite::writeFiberSite(out); }
    
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
         Mathematically:
             unbinding_rate_dt * exp( force / unbinding_force )
         is equivalent to
             exp( force / unbinding_force + log(unbinding_rate_dt) )
         and with precalculations:
             exp( force * unbinding_force_inv + log(unbinding_rate_dt) )
         */
        real x = force * prop->unbinding_force_inv + prop->unbinding_rate_log;
        //std::clog << prop->name() << " " << x << "   " << std::exp(x) << "\n";
        nextDetach -= (float)std::exp(x);
        
        return ( nextDetach < 0 );
    }
};

/// output operator
std::ostream& operator << (std::ostream&, Hand const&);

#endif

