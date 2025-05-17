// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef FIBER_SITE_H
#define FIBER_SITE_H

#include "assert_macro.h"
#include "interpolation.h"
#include "fiber.h"
#include "cymdef.h"


/// FiberSite indicates a location on a Fiber by its curvilinear abscissa
/**
 The key variable is a pointer to a Fiber, `hFiber`, which is NULL
 in the `unattached` state.
 
 In the `attached` state, the location on the Fiber is recorded using the
 curvilinear abscissa `hAbs`, measured along the fiber, from a reference
 that is fixed on the Fiber, called the Fiber's origin. `hAbs' is a signed
 continuous quantity that increases from Minus to Plus ends. The origin is
 virtual and may reside outside the Fiber ends.
 
 In this way, the value of abscissa is independent from the vertices used to
 represent the Fiber's position, and also unaffected by assembly/disassembly
 at the tips of the Fiber.
 
 If the Fiber has a Lattice, The `FiberSite` also supports binding at discrete
 positions, and in this case uses a pointer `hLattice' and a signed integer
 `hSite` to keep track of the position. The lattice uses the same origin as
 the abscissa scale, such that abscissa always corresponds to `unit * site'.
*/
class FiberSite
{
public:
    
    /// propagate Lattice cell index type
    typedef FiberLattice::lati_t lati_t;

protected:
    
    /// the Fiber of interest, or NULL
    Fiber const* hFiber;

    /// the abscissa from the origin of the Fiber
    real hAbs;
    
    /// interpolation coefficient: position = (1-inter_) * Vertex1 + inter_ * Vertex2
    mutable real inter_;

    /// index of segment where site is located: Vertex1 = hFiber::Vertex(segix_)
    mutable index_t segix_;

#if FIBER_HAS_LATTICE
    /// index in the Fiber's Lattice (a signed integer)
    lati_t hSite;
    
    /// pointer to the Lattice of the Fiber, or NULL if not in use
    FiberLattice * hLattice;
#endif

public:

#if FIBER_HAS_LATTICE
    /// default constructor
    FiberSite() : hFiber(nullptr), hAbs(0), hSite(0), hLattice(nullptr), segix_(0) {}
#else
    FiberSite() : hFiber(nullptr), hAbs(0), segix_(0) {}
#endif

    /// construct at the given distance from the origin (i.e. abscissa)
    FiberSite(Fiber const*, real a);

    /// make destructor non-virtual
    ~FiberSite() {}
    
    /// clear member variables (used in debug mode)
    void clear();

#if FIBER_HAS_LATTICE
    
    /// return associated FiberLattice if enabled
    FiberLattice* lattice() const { return hLattice; }
    
    /// index of Lattice's site
    lati_t site() const { return hSite; }
    
    /// set abscissa
    void setAbscissa(real a) { hAbs = a; }

    /// set FiberSite at index `s` and abscissa `abs`
    void engageLattice(lati_t s, real abs)
    {
        hSite = s;
        hAbs = abs;
        assert_true(hFiber->abscissaM() < hAbs + REAL_EPSILON);
        assert_true(hAbs < hFiber->abscissaP() + REAL_EPSILON);
    }

#else
    
    /// return nullptr since FiberLattice if disabled
    FiberLattice* lattice() const { return nullptr; }

#endif

    //--------------------------------------------------------------------------

    /// return the interpolation
    Interpolation interpolation() const { assert_false(bad()); return Interpolation(hFiber, inter_, segix_); }
    
    /// set Interpolation to argument
    void reinterpolate(Interpolation const& i) const { assert_true(i.mecable()==hFiber); inter_=i.coef1(); segix_=i.point1(); }
    
    /// adjust Interpolation without changing Fiber
    void reinterpolate(real a, unsigned i) const { inter_=a; segix_=i; }

    /// update the Interpolation
    void reinterpolate() const {
#if 0
        reinterpolate(hFiber->interpolate(hAbs));
#else
        real a = max_real(hFiber->segmentationInv()*(hAbs-hFiber->abscissaM()), 0);
        segix_ = std::min((unsigned)a, (unsigned)hFiber->lastSegment());
        inter_ = std::min(a-segix_, real(1));
#endif
    }
    
    /// continuous movement to given abscissa on the current fiber
    void moveTo(real a) { hAbs = a; reinterpolate(); assert_true(nullptr==lattice()); }

    /// relocate to minus end of current fiber
    void relocateM();
    
    /// relocate to plus end of current fiber
    void relocateP();

    //--------------------------------------------------------------------------
    
    /// true if not attached
    bool unattached() const { return !hFiber; }

    /// true if attached
    bool attached() const { return hFiber; }
    
    /// Fiber to which this is attached, or zero if not attached
    Fiber const* fiber() const { return hFiber; }
    
    /// Fiber to which this is attached, or zero if not attached
    Fiber * modifiableFiber() const { return const_cast<Fiber*>(hFiber); }
    
    /// position in space (using current interpolation)
    Vector pos() const { assert_false(bad()); return hFiber->midPoint(segix_, inter_); }
    
#if FIBER_HAS_FAMILY
    /// the position around which attachment is seeked
    Vector outerPos() const;
#else
    /// the position around which attachment is seeked
    Vector outerPos() const { assert_false(bad()); return hFiber->midPoint(segix_, inter_); }
#endif
    
    /// position at abscissa shifted by 'x'
    Vector posHand(real x) const { return hFiber->pos(hAbs+x); }

    /// direction of Fiber obtained by normalization
    Vector dir() const { assert_false(bad()); return hFiber->dirSegment(segix_); }
    
    /// the direction of the Fiber at the point of attachment
    Vector dirFiber() const { assert_false(bad()); return hFiber->dirSegment(segix_); }
    
    /// the abscissa, from the origin of the Fiber
    real abscissa() const { return hAbs; }

    /// abscissa, counted from the minus end
    real abscissaFromM() const { return hAbs - hFiber->abscissaM(); }

    /// abscissa, counted from the center of the fiber
    real abscissaFromC() const { return hAbs - hFiber->abscissaC(); }

    /// inverted abscissa counted from the plus end, positive if ( abscissa < abscissa(plus end) )
    real abscissaFromP() const { return hFiber->abscissaP() - hAbs; }

    /// abscissa, counted from the specified FiberEnd (in reversed direction for the plus end)
    real abscissaFrom(FiberEnd ref) const;
            
    /// nearest end to the current attachment point
    FiberEnd nearestEnd() const;
    
    /// distance to specified fiber tip
    real distanceToEnd(FiberEnd) const;
    
    /// distance to the closest fiber tip
    real distanceToNearestEnd() const;
    
    /// true if abscissa is below minus end
    bool belowM() const { return hFiber->belowM(hAbs); }

    /// true if abscissa is above minus end
    bool aboveM() const { return hFiber->aboveM(hAbs); }

    /// true if abscissa is below plus end
    bool belowP() const { return hFiber->belowP(hAbs); }
    
    /// true if abscissa is above plus end
    bool aboveP() const { return hFiber->aboveP(hAbs); }
    
    /// true if abscissa is not within the fiber's boundaries
    int outsideMP() const { return hFiber->outsideMP(hAbs); }
    
    //--------------------------------------------------------------------------
    
    /// read from file
    ObjectID readFiberSite(Inputter&, Simul&);
    
    /// write to file
    void writeFiberSite(Outputter&) const;
 
    /// Human friendly ouput
    void print(std::ostream&) const;
    
    //---------------------------------------------------------------------
    
    /// check that hAbs is within Chain::abscissaM() and Chain::abscissaP()
    int checkAbscissa() const;
    
    /// check validity of the interpolation (debuging purposes)
    int bad() const;
};

/// output operator for debugging purpose
std::ostream& operator << (std::ostream&, FiberSite const&);

/// a variable-size list of fiber sites
typedef Array<FiberSite, 0> FiberSiteList;

#endif

