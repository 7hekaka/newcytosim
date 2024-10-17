// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef FIBER_H
#define FIBER_H

#include <set>
#include <cstdint>
#include "mecafil.h"
#include "fiber_prop.h"
#include "hand_list.h"
#include "lattice.h"
#include "cymdef.h"


class Hand;
class Field;
class Single;
class FiberSet;
class FiberSegment;
class LineDisp;


/// record birth time of fibers, used for display mostly
#define FIBER_HAS_BIRTHTIME 0

/// Switch to add a Lattice of integers/floats to each Fiber {0, 1, -1}
#define FIBER_HAS_LATTICE 0

/// Flag to allow `family` member variable to control Couple's binding {0, 1}
#define FIBER_HAS_FAMILY 0

/// Flag to allow dynamic Single creation/binding at fiber's ends {0, 1}
#define FIBER_HAS_GLUE 0

/// Flag to enable the sorting of targets for attachment of Hands {0, 1}
#define BIND_CLOSEST_FIBER 1


/**
 The type of Lattice associated with each Fiber is defined here:
 */
#if FIBER_HAS_LATTICE > 0
/// Lattice cells contain integers, appropriate for discrete occupancy tracking
typedef Lattice<uint8_t> FiberLattice;
#else
/// Lattice cells contain floating-point values
typedef Lattice<real> FiberLattice;
#endif


/// type of lattice that will be displayed in play:
#if FIBER_HAS_DENSITY
typedef Lattice<real> VisibleLattice;
#else
typedef FiberLattice VisibleLattice;
#endif


/// a Mecafil to which Hands may bind
/**
 The Fiber extends the Mecafil (itself build on Chain), adding in particular
 methods that are necessary to simulate the attachment/detachment of Hand.
 It also adds a Lattice object and a FiberProp to hold parameters.
 
 - `FiberProp * prop` points to the physical properties (ie. parameters) of the Fiber.
 - `FiberDisp * disp` points to display parameters (not used in sim).
 - `fHands` keeps track of all attached Hands.
 .
 
 The Fiber may have a Lattice of integers, used by Digit and derived Hands.
 It can also have a Lattice of reals, for other features.
 
 Fibers are stored in a FiberSet.
 @todo Fiber should be called Filament
 */
class Fiber: public Mecafil
{
private:
    
    /// Disabled copy constructor
    Fiber(Fiber const&);
    
    /// disabled assignment operator
    Fiber& operator = (const Fiber&);

    /// Stores the information needed to sever a Fiber
    class CutFacts
    {
    public:
        real abscissa;    ///< abscissa of the cut, from the reference
        real cutwidth;    ///< amount removed
        state_t stateM;   ///< state of the new minus end
        state_t stateP;   ///< state of the new plus end
        
        /// constructor (abscissa, new_plus_end_state, new_minus_end_state)
        CutFacts(real a, real w, state_t p, state_t m) { abscissa=a; cutwidth=w; stateP=p; stateM=m; }
        
        /// sort from plus end to minus end, i.e. with decreasing abscissa
        bool operator < (CutFacts const&b) const { return abscissa > b.abscissa; }
    };

    /// list of bound Hands
    mutable HandList fHands;
    
#if FIBER_HAS_LATTICE
    /// Associated Lattice used for occupancy of Digit
    mutable FiberLattice fLattice;
#endif
#if FIBER_HAS_DENSITY
    /// Associated Lattice of reals
    mutable Lattice<real> fDensity;
#endif
#if FIBER_HAS_GLUE
    /// a grafted used to immobilize the Fiber
    Single * fGlue;
#endif
#if FIBER_HAS_BIRTHTIME
    /// simulation time when initialized
    real fBirthTime;
#endif

protected:
    
#if NEW_FIBER_END_CHEW
    /// stored chewing at the end
    real fChew[2];
#endif
    
    /// ordered list of future severing positions
    std::set<CutFacts> pendingCuts;

    
    /// cut Fiber at point `pti`, return section `[ pti, plus_end ]`
    virtual Fiber* severJoint(index_t pti);
    
    /// return index of point where there is a kink with ( std::cos(angle) < max_cos )
    index_t hasKink(real max_cos) const;

    
    /// viscous drag coefficient for an ellipsoid moving in an infinite volume of fluid
    static real dragCoefficientEllipsoid(real len, FiberProp const*);
    
    /// viscous drag coefficient for a cylinder moving in an infinite volume of fluid
    static real dragCoefficientCylinder(real len, FiberProp const*);
    
    /// viscous drag coefficient for a cylinder moving close to a surface
    static real dragCoefficientSurface(real len, FiberProp const*);
    
    /// add confinement interactions to a Meca
    void setFiberConfinement(Meca&, Confinement, Space const*, real stiff, real stiff2) const;
    
    /// cut fiber at distance `abs` from the minus end; returns section `[ abs, plus_end ]`
    Fiber * severSegment(real abs1, real abs2);

    /// calculate the edges for a cut of width `w` around `a` (arguments used for input/output)
    virtual void findSeverEdges(real& a, real& w);

    /// cut fiber at abscissa `[abs1, abs2]`, returning section `[ abs2, plus_end ]`
    Fiber * severNow(real abs1, real abs2, const real min);

    /// cut fiber, delete sections shorter than `min`, and set states of the newly created ends
    void severNow(real abs1, real abs2, const real min, state_t P, state_t M);

    /// perform all the cuts registered by severSoon()
    void severNow();

public:
    
#if FIBER_HAS_FAMILY
    /// if set, no connection can be made to another fiber of the same `family`
    /** This option limits the binding of Hands that are part of a Couple
     A Hand may not bind to a fiber, if the other Hand of the Couple is already
     attached to a fiber with the same value of `family`, if ( family > 0 ).
     */
    Fiber const* family_;
    Fiber const* sister_;
    Fiber const* brother_;

    /// position of a point specified by distance from the minus end
    Vector displayPosM(real a) const;

#else
    
    /// position of a point specified by distance from the minus end
    Vector displayPosM(real a) const { return posM(a); }

#endif

    /// the Property of this object
    FiberProp const * prop;
    
    /// the display parameters
    LineDisp * disp;

    //--------------------------------------------------------------------------

    /// constructor
    Fiber(FiberProp const*);
    
    /// mark for deletion
    void adieu();
    
    /// destructor
    virtual ~Fiber();

    /// prepare for Meca
    void prepareMecable();

    /// calculate viscous drag coefficient
    void setDragCoefficient();
    
    /// add force elements relevant for this Fiber to a Meca
    void setInteractions(Meca&) const;
    
    //--------------------------------------------------------------------------

    /// adjust abscissa of Hands by applying mirror image around fiber midpoint
    void flipHandsPolarity();
    
    /// flip Vertices and Hands, exchanging Plus and Minus ends
    void flipPolarity();

    /// remove the portion of size `len` that includes the minus end
    void cutM(real len);
    
    /// remove the portion of size `len` that includes the plus end
    void cutP(real len);
    
    /// Cut all segments intersecting the plane defined by <em> n.pos + a = 0 </em>
    void planarCut(Vector const& n, real a, state_t stateP, state_t stateM, real min_len);

    /// register a cut at abscissa `a` from the ORIGIN, with `m` and `p` the states of the new ends
    void severSoon(real a, real w, state_t p, state_t m) { pendingCuts.insert(CutFacts(a, w, p, m)); }

    /// call Chain::join(), and transfer Hands (caller should delete `fib`).
    virtual void join(Fiber *);
    
    /// implements some growth/shrinkage at the ends as part of a simulation step
    bool updateLength(real, real, bool split = true);
    
    /// update Lattice and Density ranges
    void updateRange(Field*);

    /// should be called if any Fiber tip has elongated or shortened
    void updateFiber();
    
    /// simulation step
    virtual void step();

    //--------------------------------------------------------------------------

    /// the energy due to bending elasticity: 1/2 * rigidity * sum( curvature(s)^2 ds ),
    real bendingEnergy() const { return bendingEnergy0() * prop->rigidity; }
    
    /// return abscissa of the closest point to `w`, and set `dis2` to the square of the distance
    real projectPoint(Vector const& w, real & dis2) const;
    
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of minus end
    virtual state_t endStateM() const { return STATE_WHITE; }

    /// return assembly/disassembly state of plus end
    virtual state_t endStateP() const { return STATE_WHITE; }

    /// return assembly/disassembly state of the FiberEnd
    state_t endState(FiberEnd) const;

    
    /// change state of minus end
    virtual void setEndStateM(state_t) {}

    /// change state of plus end
    virtual void setEndStateP(state_t) {}

    /// change state of FiberEnd
    void setEndState(FiberEnd, state_t);
    
    //--------------------------------------------------------------------------
    
    /// register a new Hands that attached to this Fiber
    void addHand(Hand* h) const { fHands.add(h); }
    
    /// unregister bound Hands (which has detached)
    void removeHand(Hand* h) const { fHands.remove(h); }
    
    /// update all Hands bound to this
    void updateHands() const;
    
    /// sort Hands by order of increasing abscissa
    void sortHands() const { fHands.sort(); }
    
    /// return Hand bound to this fiber (use ->next() to access all other Hands)
    Hand * firstHand() const { return fHands.front(); }
   
    /// number of attached Hands
    size_t nbAttachedHands() const { return fHands.count(); }
    
    /// count attached Hands fitting a given criteria
    long nbAttachedHands(int (*func)(Hand const*)) const { return fHands.count(func); }

    /// number of Hands attached within a range of abscissa
    size_t nbHandsInRange(real abs_min, real abs_max, FiberEnd ref) const;
    
    /// number of Hands attached at a distance less than 'len' from the specified FiberEnd
    size_t nbHandsNearEnd(real len, FiberEnd) const;
    
    /// create Hands attached to this fiber
    void makeAttachedHands(ObjectList&, std::string const&, size_t, Glossary&, std::string const&, Simul&);
    
    //--------------------------------------------------------------------------

#if NEW_FIBER_END_CHEW
    /// register a chewing effect of magnitude 'x'
    void chew(FiberEnd e, const real x) { assert_true(e>0); fChew[e-1] += x; }
#endif
#if FIBER_HAS_BIRTHTIME
    /// returns simulation time at which Fiber was created
    real birthTime() const { return fBirthTime; }

    /// set birth time
    void birthTime(double t) { fBirthTime = t; }

    /// returns current age of the fiber
    double age() const;
#else
    /// set birth time
    void birthTime(double) {}

    /// returns simulation time at which Fiber was created
    real birthTime() const { return 0; }

    /// returns current age of the fiber
    double age() const { return 0; }
#endif

    //--------------------------------------------------------------------------
    
#if FIBER_HAS_LATTICE
    /// modifiable reference to Fiber's Lattice
    FiberLattice * lattice() { return &fLattice; }
    
    /// modifiable reference to Fiber's Lattice
    FiberLattice * lattice() const { return &fLattice; }
    
    /// recalculate all cell values given current list of bound digits
    bool resetLattice(bool);
#else
    /// does nothing
    bool resetLattice(bool) { return false; }
#endif
    
    /// record minium, maximum and sum of lattice values
    void infoLattice(size_t& cnt, size_t& vac, real& sum, real& mn, real& mx) const;

    /// print Lattice data (for debugging purpose)
    void printLattice(std::ostream&) const;

#if FIBER_HAS_DENSITY

    /// modifiable reference to Fiber's density object
    Lattice<real>& densityField() { return fDensity; }

    /// value of the density at given abscissa
    real density(real a) const { if ( fDensity.ready() ) return fDensity.cell(a); return 0; }

#endif
    
    /// initialize lattice sites to represent a constant linear density
    void setDensity(Lattice<real>&, real density) const;

    /// transfer all lattice substance to the Field
    void releaseDensity(Lattice<real>&, Field*) const;

    /// update lattice values as `value <- cst + fac * value`
    void evolveDensity(Lattice<real>&, real cst, real fac) const;

    /// transfer from Field to Lattice at rate `on` and back at rate `off`
    void equilibrateDensity(Lattice<real>&, Field*, real on, real off) const;
    
    /// transfer from Field to Lattice at rate `on`
    void bindDensity(Lattice<real>&, Field*, real rate) const;
    
    /// transfer from Field to Lattice at rate `on`
    void fluxDensity(Lattice<real>&, Field*, real speed) const;
    
    /// sever fiber proportionally to the quantity stored in the Lattice
    void cutFiberByDensity(Lattice<real>&);

    /// find minium, maximum and sum of density values
    void infoDensity(real& len, size_t&, size_t&, real& sm, real& mn, real& mx, bool density) const;

    /// lattice to be displayed
    VisibleLattice const* visibleLattice() const;
    
    //--------------------------------------------------------------------------
    
    /// set Space glue for pure pushing
    void setGlue1(Single* glue, FiberEnd, Space const*);
    
    /// set Space glue for pure pulling
    void setGlue2(Single* glue, FiberEnd, Space const*);
    
    /// set Space glue for pushing and pulling
    void setGlue3(Single* glue, Space const*);
    
    /// set Solid glue type 4
    void setGlueG(Single* glue, FiberEnd);
    
    /// set Solid glue type 5
    void setGlueE(Single* glue, FiberEnd);

    /// a setGlue to rule them all
    void setGlue(Single*& glue, FiberEnd, int mode);
    
    /// create a Single that can be used as glue
    void makeGlue(Single*& glue);
    
    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Fiber * next() const { return static_cast<Fiber*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Fiber * prev() const { return static_cast<Fiber*>(prev_); }

    //--------------------------------------------------------------------------
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'f';
    
    /// identifies angle data format
    static const ObjectTag COMPACT_TAG = 'g';
    
    /// identifies Age and Chiasma info (must be upper case as this is meta-data)
    static const ObjectTag FIBINFO_TAG = 'G';

    /// identifies data for dynamic ends of fibers (must be upper case)
    static const ObjectTag DYNAMIC_TAG = 'F';
    
    /// identifies FiberLattice data (was 'l' before 23/06/2021; must be upper case)
    static const ObjectTag LATTICE_TAG = 'T';
    
    /// identifies Lattice<real> data (was 'L' before 23/06/2021; must be upper case)
    static const ObjectTag FIBMESH_TAG = 'M';

    /// return unique character identifying the class
    ObjectTag tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// return specification of fiber class
    virtual std::string activity() const { return "none"; }
    
    /// convert pointer to Fiber* if the conversion seems valid; returns 0 otherwise
    static Fiber* toFiber(Object * obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Fiber*>(obj);
        return nullptr;
    }
    
    /// convert pointer to Fiber* if the conversion seems valid; returns 0 otherwise
    static Fiber const* toFiber(Object const* obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Fiber const*>(obj);
        return nullptr;
    }

    //--------------------------------------------------------------------------

    /// write to file
    void write(Outputter&) const;
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// check data integrity
    int bad() const;

};

#endif

