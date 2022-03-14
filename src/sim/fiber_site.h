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

    /// allowing update of the interpolation
    friend void Fiber::updateHands();

protected:
    
    /// the interpolation on the Fiber's vertices; `bad()` is used to check validity
    mutable Interpolation hTerp;

protected:
    
    /// the Fiber of interest, or NULL
    Fiber * hFiber;

#if FIBER_HAS_LATTICE
    /// pointer to the Lattice of the Fiber, or NULL if not in use
    FiberLattice * hLattice;
    
    /// index in the Fiber's Lattice (a signed integer)
    lati_t hSite;
#endif
    
    /// the abscissa from the origin of the Fiber
    real hAbs;

public:

#if FIBER_HAS_LATTICE
    /// default constructor
    FiberSite() : hFiber(nullptr), hLattice(nullptr), hSite(0), hAbs(0) {}
#else
    FiberSite() : hFiber(nullptr), hAbs(0) {}
#endif

    /// construct at the given distance from the origin
    FiberSite(Fiber*, real a);

    /// make destructor non-virtual
    ~FiberSite() {}
    
#if FIBER_HAS_LATTICE
    
    /// return Lattice if engaged
    FiberLattice* lattice() const { return hLattice; }
    
    /// index of Lattice's site
    lati_t site() const { return hSite; }
    
    /// set FiberSite at index `s` with an abscissa `off` within the site
    void engageLattice(FiberLattice* l, lati_t s, real off)
    {
        hLattice = l;
        hSite    = s;
        hAbs     = l->unit() * s + off;
        //assert_true(hFiber->abscissaM() < hAbs + REAL_EPSILON);
        //assert_true(hAbs < hFiber->abscissaP() + REAL_EPSILON);
    }

#else
    
    FiberLattice* lattice() const { return nullptr; }

#endif
    //--------------------------------------------------------------------------

    /// return the interpolation
    const Interpolation& interpolation() const { assert_false(bad()); return hTerp; }
    
    /// update the Interpolation
    void reinterpolate() const { hTerp = hFiber->interpolate(hAbs); }
    
    /// move to a different abscissa on the current fiber
    void moveTo(real a) { hAbs = a; reinterpolate(); }

    /// relocate to MINUS_END of current fiber
    void relocateM();
    
    /// relocate to PLUS_END of current fiber
    void relocateP();

    //--------------------------------------------------------------------------
    
    /// true if not attached
    bool unattached() const { return !hFiber; }

    /// true if attached
    bool attached() const { return hFiber; }
    
    /// Fiber to which this is attached, or zero if not attached
    Fiber* fiber() const { return hFiber; }
    
    /// position in space (using current interpolation)
    Vector pos() const { assert_false(bad()); return hTerp.pos(); }
    
#if FIBER_HAS_FAMILY
    /// the position around which attachment is seeked
    Vector outerPos() const;
#else
    /// the position around which attachment is seeked
    Vector outerPos() const { assert_false(bad()); return hTerp.pos(); }
#endif
    
    /// position at abscissa shifted by 'x'
    Vector posHand(real x) const { return hFiber->pos(hAbs+x); }

    /// direction of Fiber obtained by normalization
    Vector dir() const { assert_false(bad()); return hTerp.dir(); }
    
    /// the direction of the Fiber at the point of attachment
    Vector dirFiber() const { assert_false(bad()); return hFiber->dirSegment(hTerp.point1()); }
    
    /// the abscissa, from the origin of the Fiber
    real abscissa() const { return hAbs; }

    /// abscissa, counted from the MINUS_END
    real abscissaFromM() const { return hAbs - hFiber->abscissaM(); }

    /// inverted abscissa counted from the PLUS_END, positive if ( abscissa < abscissa(PLUS_END) )
    real abscissaFromP() const { return hFiber->abscissaP() - hAbs; }

    /// abscissa, counted from the specified FiberEnd (in reversed direction for the PLUS_END)
    real abscissaFrom(FiberEnd ref) const;
            
    /// nearest end to the current attachment point
    FiberEnd nearestEnd() const;
    
    /// distance to the closest fiber tip
    real distanceToEnd(FiberEnd) const;

    /// true if abscissa is below abscissaP
    bool belowP() const { return hFiber->belowP(hAbs); }
    
    /// true if abscissa is above abscissaM
    bool aboveM() const { return hFiber->aboveM(hAbs); }
    
    /// true if abscissa is not within the fiber's boundaries
    bool outsideMP() const { return hFiber->outsideMP(hAbs); }
    
    //--------------------------------------------------------------------------
    
    /// read from file
    void read(Inputter&, Simul&);
    
    /// write to file
    void write(Outputter&) const;
 
    /// Human friendly ouput
    void print(std::ostream&) const;
    
    //---------------------------------------------------------------------
    
    /// check that hAbs is within Fiber::abscissaM() and Fiber::abscissaP()
    int checkAbscissa() const;
    
    /// check validity of the interpolation (debuging purposes)
    int bad() const;
};

/// output operator for debugging purpose
std::ostream& operator << (std::ostream&, FiberSite const&);


#endif

