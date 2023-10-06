// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef DIGIT_H
#define DIGIT_H

#include "hand.h"
#include "digit_prop.h"


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
    DigitProp const* prop() const { return static_cast<DigitProp const*>(Hand::prop); }

    /// constructor
    Digit(DigitProp const*, HandMonitor*);
    
    //--------------------------------------------------------------------------

#if FIBER_HAS_LATTICE

    /// true if given Lattice's site is outside Lattice's range
    int outsideMP(lati_t s) const { return hLattice->outsideMP(s); }
    
    /// true if abscissa is below minus end
    bool belowM(lati_t s) const { return hLattice->belowM(s); }
    
    /// true if abscissa is above plus end
    bool aboveP(lati_t s) const { return hLattice->aboveP(s); }

#else
    
    /// ersatz function converting to site when there is no lattice
    lati_t site() const { return (lati_t)std::nearbyint(hAbs/prop()->step_size); }
    
    /// converting site to abscissa when there is no lattice
    real abscissa_(lati_t s) const { return s * prop()->step_size + prop()->site_shift; }

    int outsideMP(lati_t s) const { return hFiber->outsideMP(abscissa_(s)); }
    
    /// true if abscissa is above abscissaM
    bool belowM(lati_t s) const { return hFiber->belowM(abscissa_(s)); }
    
    /// true if abscissa is below abscissaP
    bool aboveP(lati_t s) const { return hFiber->aboveP(abscissa_(s)); }

#endif
#if FIBER_HAS_LATTICE > 0

    /// true if none of the Lattice's site bits matches the footprint
    bool vacant(lati_t s) const { return 0 == (hLattice->data(s) & prop()->footprint); }

    /// flip footprint bits on current site
    void inc() const { assert_true(vacant(hSite)); hLattice->data(hSite) ^= prop()->footprint; }

    /// flip footprint bits on current site
    void dec() const { hLattice->data(hSite) ^= prop()->footprint; assert_true(vacant(hSite)); }
    
#elif FIBER_HAS_LATTICE < 0

    /// true if given Lattice's site is zero
    bool vacant(lati_t s) const { return hLattice->data(s) == 0.0; }

    /// add 1.0 to Lattice's site
    void inc() const { hLattice->data(hSite) += 1.0; }

    /// subtract 1.0 to Lattice's site
    void dec() const { hLattice->data(hSite) -= 1.0; }
    
#else

    bool vacant(lati_t) const { return true; }
    void inc() const {}
    void dec() const {}
    
#endif

    /// check if attachement is possible according to properties
    bool attachmentAllowed(FiberSite&) const;

    /// attach and update variables
    void attach(FiberSite const&);
    
    /// detach
    void detach();

    
    /// transfer to given site
    void hop(lati_t);

    /// transfer to given site if it is vacant
    void jumpTo(lati_t p) { if ( vacant(p) ) hop(p); }
    
    /// relocate without checking intermediate sites
    void jumpToEndM() { jumpTo(lattice()->entry()); }

    /// relocate without checking intermediate sites
    void jumpToEndP() { jumpTo(lattice()->fence()); }

    
    /// attempt one step towards the plus end
    void stepP()      { jumpTo(site()+1); }
    
    /// attempt one step towards the minus end
    void stepM()      { jumpTo(site()-1); }

    
    /// attempt one step of size `s` towards the plus end
    void jumpP(int s) { jumpTo(site()+s); }
    
    /// attempt one step of size `s` towards the minus end
    void jumpM(int s) { jumpTo(site()-s); }

    
    /// attempt `n` steps towards the plus end, checking all intermediate sites
    void crawlP(int n);
    
    /// attempt `n` steps towards the minus end, checking all intermediate sites
    void crawlM(int n);

    
    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
 
    
    /// this is called when the attachment point is beyond the plus end
    void handleDisassemblyM();
    
    /// this is called when the attachment point is below the minus end
    void handleDisassemblyP();
    
    /// Promote a Digit class that is not bound to the lattice
    void promote();
};

/// output operator
std::ostream& operator << (std::ostream&, Digit const&);

#endif

