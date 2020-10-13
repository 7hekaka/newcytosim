// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_H
#define FIBER_H

#include <set>
#include <stdint.h>
#include "mecafil.h"
#include "fiber_prop.h"
//#include "node_list.h"
#include "hand_list.h"
#include "lattice.h"
#include "sim.h"


class Hand;
class Field;
class Single;
class FiberSet;
class FiberSegment;
class LineDisp;


/// Flag to add a Lattice of integers to each Fiber {0, 1}
#define FIBER_HAS_LATTICE 1

/// Flag to add a Lattice of reals to each Fiber {0, 1}
#define FIBER_HAS_MESH 0

/// Flag to allow `family` member variable to control Couple's binding {0, 1}
#define FIBER_HAS_FAMILY 0

/// Flag to allow dynamic Single creation/binding at fiber's ends {0, 1}
#define FIBER_HAS_GLUE 1

/// Flag to enable the sorting of targets for attachment of Hands {0, 1}
#define BIND_CLOSEST_FIBER 1


/**
 The type of Lattice associated with each Fiber is defined here:
 */
#if FIBER_HAS_LATTICE > 0
// Lattice composed of integers, appropriate for discrete occupancy
typedef Lattice<uint16_t> FiberLattice;
#else
// Lattice composed of floating point values, for continuous values
typedef Lattice<real> FiberLattice;
#endif


/// type of lattice that will be displayed in play:
#if FIBER_HAS_MESH
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
 - `frHands` keeps track of all attached Hands.
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
    Fiber& operator =(const Fiber&);

    /// Stores the information needed to sever a Fiber
    class SeverPos
    {
    public:
        real    abs;      ///< abscissa of the cut, from the reference
        state_t stateM;   ///< state of the new MINUS_END
        state_t stateP;   ///< state of the new PLUS_END
        
        /// constructor (abscissa, new_plus_end_state, new_minus_end_state)
        SeverPos(real a, state_t p, state_t m) { abs=a; stateP=p; stateM=m; }
        
        /// sort from PLUS_END to MINUS_END, i.e. with decreasing abscissa
        bool operator < (SeverPos const&b) const { return abs > b.abs; }
    };

    /// list of bound Hands
    mutable HandList    frHands;
    
#if FIBER_HAS_LATTICE
    /// Associated Lattice used for occupancy of Digit
    FiberLattice        frLattice;
#endif
#if FIBER_HAS_MESH
    /// Associated Lattice of reals
    Lattice<real>       frMesh;
#endif
#if FIBER_HAS_GLUE
    /// a grafted used to immobilize the Fiber
    Single *            frGlue;
#endif

protected:
#if NEW_FIBER_CHEW
    /// stored chewing at the end
    real                frChewM, frChewP;
#endif
    
    /// ordered list of future severing positions
    std::set<SeverPos>  pendingCuts;

    
    /// cut Fiber at point `pti`, return section `[ pti - PLUS_END ]`
    virtual Fiber* severPoint(size_t pti);
    
    /// return index of point where there is a kink with ( std::cos(angle) < max_cos )
    size_t         hasKink(real max_cos) const;

    
    /// viscous drag coefficient for an ellipsoid moving in an infinite volume of fluid
    static real    dragCoefficientEllipsoid(real len, FiberProp const*);
    
    /// viscous drag coefficient for a cylinder moving in an infinite volume of fluid
    static real    dragCoefficientCylinder(real len, FiberProp const*);
    
    /// viscous drag coefficient for a cylinder moving close to a surface
    static real    dragCoefficientSurface(real len, FiberProp const*);
    
public:
    
#if FIBER_HAS_FAMILY
    /// if set, no connection can be made to another fiber of the same `family`
    /** This option limits the binding of Hands that are part of a Couple
     A Hand may not bind to a fiber, if the other Hand of the Couple is already
     attached to a fiber with the same value of `family`, if ( family > 0 ).
     */
    Fiber const*  family_;
    Fiber const*  sister_;
    Fiber const*  brother_;

    /// radial direction at the specified distance from the MINUS_END
    Vector  radialDirM(real a) const { assert_true(this!=family_); return posM(a) - family_->posM(a); }

    /// radial direction at the specified abscissa
    Vector  radialDir(real a) const { return radialDirM(a-abscissaM()); }
    
    /// position of a point specified by distance from the MINUS_END
    Vector  displayPosM(real a) const;

#else
    
    /// position of a point specified by distance from the MINUS_END
    Vector  displayPosM(real a) const { return posM(a); }

#endif

    /// the Property of this object
    FiberProp const*    prop;
    
    /// the display parameters
    LineDisp mutable*   disp;

    //--------------------------------------------------------------------------

    /// constructor
    Fiber(FiberProp const*);
    
    /// destructor
    virtual ~Fiber();

    //--------------------------------------------------------------------------
    
    /// prepare for Meca
    void           prepareMecable();

    /// calculate viscous drag coefficient
    void           setDragCoefficient();
    
    /// add interactions to a Meca
    void           setInteractions(Meca&) const;
    

    /// invert polarity and adjust abscissa of Hands to keep them at the same place
    void           flipHandsPolarity();
    
    /// remove the portion of size `len` that includes the MINUS_END
    void           cutM(real len);
    
    /// remove the portion of size `len` that includes the PLUS_END
    void           cutP(real len);
    
    /// Cut all segments intersecting the plane defined by <em> n.pos + a = 0 </em>
    void           planarCut(Vector const& n, real a, state_t stateP, state_t stateM);
    
    /// cut fiber at distance `abs` from the MINUS_END; returns section `[ abs - PLUS_END ]`
    Fiber *        severP(real abs);

    /// cut fiber at abscissa `abs`; returns section `[ abs - PLUS_END ]`
    Fiber *        severNow(real abs) { return severP(abs-abscissaM()); }

    /// register a cut at abscissa `a` from the ORIGIN, with `m` and `p` the states of the new ends
    void           sever(real a, state_t p, state_t m) { pendingCuts.insert(SeverPos(a, p, m)); }
    
    /// perform all the cuts registered by sever()
    void           severNow();

#if NEW_FIBER_CHEW
    /// register a chewing quantity
    void           chew(const real x, FiberEnd end) { if ( end == PLUS_END ) frChewP += x; else frChewM += x; }
#endif

    /// call Chain::join(), and transfer Hands (caller should delete `fib`).
    virtual void   join(Fiber * fib);
    
    /// simulation step
    virtual void   step();
    
    /// called if a Fiber tip has elongated or shortened
    void           updateFiber();
    
    //--------------------------------------------------------------------------

    /// the energy due to bending rigidity: 1/2 * rigidity * sum( curvature(s)^2 ds ),
    real           bendingEnergy() const { return bendingEnergy0() * prop->rigidity; }
    
    /// return the abscissa of the closest position to `w` on this Fiber, and set `dis` to the square of the distance
    real           projectPoint(Vector const& w, real & dis) const;
    
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of MINUS_END
    virtual state_t endStateM() const { return STATE_WHITE; }

    /// return assembly/disassembly state of PLUS_END
    virtual state_t endStateP() const { return STATE_WHITE; }

    /// return assembly/disassembly state of the FiberEnd
    state_t         endState(FiberEnd end) const;

    
    /// change state of MINUS_END
    virtual void   setEndStateM(state_t) {}

    /// change state of PLUS_END
    virtual void   setEndStateP(state_t) {}

    /// change state of FiberEnd `end` to `s`
    void           setEndState(FiberEnd end, state_t s);
    
    
    /// the length of freshly assembled polymer at the MINUS_END during the last time step
    virtual real   freshAssemblyM() const { return 0; }

    /// the length of freshly assembled polymer at the PLUS_END during the last time step
    virtual real   freshAssemblyP() const { return 0; }

    /// the length of freshly assembled polymer during the last time step
    real           freshAssembly(FiberEnd end) const;
    
    
    /// true if the tip `end` has grown in the last time step ( freshAssembly(which) > 0 )
    bool           isGrowing(FiberEnd end) const { return freshAssembly(end) > 0; }
    
    /// true if the tip `end` has shrunk in the last time step ( freshAssembly(which) < 0 )
    bool           isShrinking(FiberEnd end) const { return freshAssembly(end) < 0; }
    
    //--------------------------------------------------------------------------
    
    /// register a new Hands that attached to this Fiber
    void           addHand(Hand* h) const { frHands.add(h); }
    
    /// unregister bound Hands (which has detached)
    void           removeHand(Hand* h) const  { frHands.remove(h); }
    
    /// update all Hands bound to this
    void           updateHands() const { frHands.update(); }

    /// detach all Hands
    void           detachHands() const { frHands.detachAll(); }
    
    /// sort Hands by order of increasing abscissa
    void           sortHands() const { frHands.sort(); }
    
    /// return Hand bound to this fiber (use ->next() to access all other Hands)
    Hand *         firstHand() const { return frHands.front(); }
   
    /// number of attached Hands
    size_t         nbHands() const { return frHands.count(); }
    
    /// a function to count Hands using a custom criteria
    int            nbHands(int (*func)(Hand const*)) const { return frHands.count(func); }

    /// number of Hands attached within a range of abscissa
    size_t         nbHandsInRange(real abs_min, real abs_max, FiberEnd ref) const;
    
    /// number of Hands attached at a distance less than 'len' from the specified FiberEnd
    size_t         nbHandsNearEnd(real len, FiberEnd end) const;
    
    //--------------------------------------------------------------------------
#if FIBER_HAS_LATTICE
    /// modifiable reference to Fiber's Lattice
    FiberLattice&  lattice() { return frLattice; }
    
    /// const reference to Fiber's Lattice
    FiberLattice const& lattice() const { return frLattice; }
        
    /// recalculate occupancy lattice from bound Digits
    void           resetLattice();
#else
    /// does nothing
    void           resetLattice() {}
#endif
    
    /// record minium, maximum and sum of lattice values
    void           infoLattice(real& len, size_t&, real& sm, real& mn, real& mx) const;

    /// print Lattice data (for debugging purpose)
    void           printLattice(std::ostream&) const;

    
#if FIBER_HAS_MESH

    /// modifiable reference to Fiber's Lattice
    Lattice<real> const&  mesh() const { return frMesh; }

    /// value of the frMesh at given abscissa
    real           meshValue(real a) const { if ( frMesh.ready() ) return frMesh.cell(a); return 0; }

#endif
    
    /// initialize lattice sites to represent a constant linear density
    void           setMeshValues(Lattice<real>&, real density) const;

    /// transfer all lattice substance to the Field
    void           releaseMeshValues(Lattice<real>&, Field*) const;

    /// update lattice values as `value <- cst + fac * value`
    void           evolveMeshValues(Lattice<real>&, real cst, real fac) const;

    /// transfer from Field to Lattice at rate `on` and back at rate `off`
    void           equilibrateMesh(Lattice<real>&, Field*, real on, real off) const;
    
    /// transfer from Field to Lattice at rate `on`
    void           bindMesh(Lattice<real>&, Field*, real rate) const;
    
    /// transfer from Field to Lattice at rate `on`
    void           fluxMesh(Lattice<real>&, Field*, real speed) const;
    
    /// sever fiber proportionally to the quantity stored in the Lattice
    void           cutFiberMesh(Lattice<real>&);

    /// find minium, maximum and sum of mesh values
    void           infoMesh(real& len, size_t&, real& sm, real& mn, real& mx, bool density) const;

    /// lattice to be displayed
    VisibleLattice const* visibleLattice() const;
    
    //--------------------------------------------------------------------------
    
    /// set Space glue for pure pushing
    void           setGlue1(Single* glue, FiberEnd, Space const*);
    
    /// set Space glue for pure pulling
    void           setGlue2(Single* glue, FiberEnd, Space const*);
    
    /// set Space glue for pushing and pulling
    void           setGlue3(Single* glue, Space const*);
    
    /// set Solid glue
    void           setGlueG(Single* glue, FiberEnd);

    /// a setGlue to rule them all
    void           setGlue(Single*& glue, FiberEnd);
    
    /// create a Single that can be used as glue
    void           makeGlue(Single*& glue);
    
    //--------------------------------------------------------------------------

    /// a static_cast<> of Node::next()
    Fiber *  next() const { return static_cast<Fiber*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Fiber *  prev() const { return static_cast<Fiber*>(nPrev); }

    //--------------------------------------------------------------------------
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'f';
    
    /// identifies data for dynamic ends of fibers
    static const ObjectTag TAG_DYNAMIC = 'F';
    
    /// identifies FiberLattice data
    static const ObjectTag TAG_LATTICE = 'l';
    
    /// identifies Lattice<real> data
    static const ObjectTag TAG_FIBMESH = 'L';

    /// return unique character identifying the class
    ObjectTag       tag() const { return TAG; }
    
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
    void        write(Outputter&) const;
    
    /// read from file
    void        read(Inputter&, Simul&, ObjectTag);

};

#endif

