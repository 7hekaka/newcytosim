// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef CHAIN_H
#define CHAIN_H

#include "cymdef.h"
#include "vector.h"
#include "mecable.h"

class Interpolation;
class Mecapoint;
class Glossary;

/// include a normal used for fancy display of fibers as helices
#define FIBER_HAS_NORMAL 0


/// Mecable with linear geometry
/**
 This class describes a thin flexible filament that is longitudinally incompressible.
 The curvilinear length of the filament can be changed by growP(), growM(), cutP() and cutM().
 
 \par Number of points:
 
 The best number of points to describe a Chain is automatically calculated:
 It is the integer `number_of_points` that minimizes:
 
    abs_real( length() / number_of_points - FiberProp::segmentation )
 
 where FiberProp::segmentation is a parameter of the Fiber class.
 All the segments in a fiber all have the same length

    Chain::segmentation() = length() / ( number_of_points - 1 )

 Note that Chain::segmentation() is not always equal to FiberProp::segmentation.
 If the fibers have various length, their segmentation() will be different,
 even though they all share the same value of FiberProp::segmentation.

 See related functions: length(), nbPoints() and segmentation().
 
 \par Longitudinal incompressibility:
 
 Successive vertices are kept at a constant distance via constrained dynamics:

    norm( posPoint(N+1)-posPoint(N) ) == Chain::segmentation()
 
 \par Origin:
 
 An abscissa is a curvilinear distance taken along the Fiber,
 and the Chain provides an origin to make this independent of the vertices. 
 Thus even if the fiber lengthen from its ends, a position described by an abscissa will
 stay associated with the same local lattice site.
 
 Functions are provided in Chain to convert abscissa measured from different references,
 and to obtain positions of the fiber for a given abcissa.

 \par Derived classes:
 
 The class FiberSite keeps track of its position using an abscissa from the origin,
 and all Hand objects are built from this class.
 The class Fiber keeps track of the FiberSite that are attached to itself.
 
*/
class Chain : public Mecable
{
public:
    
    /// the ideal number of points for ratio = length / segmentation
    static index_t bestNumberOfPoints(real ratio);

    /// calculate length of given string of points
    static real contourLength(const real pts[], size_t n_pts);
    
private:
        
    /// actual section length: distance between consecutive points, and inverse
    real fnCut, iCut;
    
    /// target segmentation (copy of 'FiberProp::segmentation')
    real fnSegmentation;
    
    /// abscissa of the plus-end (equal to length initially)
    real fnAbscissaP;

    /// abscissa of the minus-end (equal to zero initially)
    real fnAbscissaM;

#if FIBER_HAS_NORMAL
    /// vector orthogonal to backbone at the origin, used for display only
    mutable Vector3 fnNormal;
#endif
    
protected:
    
    /// length increment at minus end during last time-step
    real cDeltaM;
    
    /// length increment at plus end during last time-step
    real cDeltaP;

    /// flag to update
    bool needUpdate;

    /// this signals that update is needed, to be called after a change in length
    void postUpdate() { needUpdate = true; }
    
    /// called if a Fiber tip has elongated or shortened
    virtual void updateFiber() {}

    /// restore the distance between two points
    static void reshape_two(const real*, real*, real cut);

    /// oldest method to restore the distance between successive vertices
    static void reshape_global(index_t, const real*, real*, real cut);

    /// apply the forces movements needed to the distance between two points
    static void reshape_apply(index_t, const real*, const real*, real*);

    /// iterative method to restore the distance between successive vertices
    static int reshape_calculate(index_t, real, index_t max_iter, real const*, real const*, real const*, real*, size_t);

    /// iterative method to restore the distance between successive vertices
    static int reshape_local(index_t, const real*, real*, real cut, index_t max_iter, real* tmp, size_t);

    /// change segmentation
    void setSegmentation(real c) { fnCut = std::max(c, REAL_EPSILON); iCut = real(1) / fnCut; }
    
public:
    
    /// Constructor
    Chain();
    
    /// Destructor
    ~Chain() {}
    
    /// Number of segments = nbPoints() - 1
    index_t nbSegments()  const { return nPoints - 1; }
    
    /// Index of the last segment = nbPoints() - 2
    index_t lastSegment() const { return nPoints - 2; }
    
    /// return P where segment [ P, P+1 [ contains point at distance `a` from the minus end
    /** returns 0 if `a < 0` and last point index if `a > lastSegment()` */
    index_t indexSegmentM(const real a) const { return std::min(index_t(std::max(a,(real)0)/fnCut), lastSegment()); }

    //---------------------

    /// set position of minus end with given direction (length and Nb of points are not modified)
    /** dir does not need to be normalized */
    void setStraight(Vector const& pos, Vector const& dir);
    
    /// set position of minus end, direction and length of Fiber
    void setStraight(Vector const& pos, Vector const& dir, real len);
    
    /// set on the surface of a sphere of radius 'rad', in direction 'dir'
    void setCurved(Vector dir, real rad, real len, real off);

    /// translate Fiber to place 'ref' at the position where the CENTER is located
    void placeEnd(FiberEnd ref);
    
    /// set shape with `np` points from the given array of size DIM*n_pts
    void setShape(const real pts[], size_t n_pts, index_t np);

    /// set shape as a random walk with given parameters
    void setEquilibrated(real length, real persistence_length);

    /// change the current segmentation to force `length()==len` (normally not needed)
    void imposeLength(real len) { setSegmentation(len/real(nbSegments())); fnAbscissaP = fnAbscissaM + len; }
    
    /// Number of distance constraints applied to the movements of vertices
    index_t nbConstraints() const { return nPoints - 1; }
    
    /// change Lagrange multipliers (do not use: this is done by computeTensions)
    virtual void setTensions(const real*) {}

    //---------------------
    
    /// return Mecapoint of given end
    Mecapoint exactEnd(FiberEnd) const;

    /// interpolation of minus end
    Interpolation interpolateEndM() const;

    /// interpolation of plus end
    Interpolation interpolateEndP() const;

    /// interpolation of given end
    Interpolation interpolateEnd(FiberEnd) const;

    /// interpolation of the mid-point between the two ends
    Interpolation interpolateCenter() const;
    
    /// interpolation of the site specified from the minus end
    Interpolation interpolateM(real ab) const;
    
    /// interpolation of a site specified by its distance from a FiberEnd
    Interpolation interpolateFrom(real ab, FiberEnd ref) const;
    
    /// interpolation of the site specified by its distance from the ORIGIN
    Interpolation interpolateAbs(real ab) const;

    //---------------------
    
    /// length of the Fiber, estimated from the difference of abscissa at the ends
    /** This is the quantity that defines the length that the filament should have! */
    real length()                const { return fnAbscissaP - fnAbscissaM; }

    /// length of the Fiber, estimated from the segmentation and number of segments
    /** This should be equal to length(), if computers did exact calculus */
    real length1()               const { return nPoints * fnCut - fnCut; }
    
    /// the sum of the distance between consecutive vertices (used for debugging)
    real contourLength()         const { return contourLength(pPos, nPoints); }
    
    /// true if `( abscissaM() <= a ) AND ( a <= abscissaP() )`
    bool betweenMP(const real a) const { return abscissaM() <= a + REAL_EPSILON && a <= abscissaP() + REAL_EPSILON; }
    
    /// return 2 if site `i` is partly or entirely below the minus end, 1 if above the plus end
    int outsideMP(const real a) const { return 2 * ( a < abscissaM() ) | ( abscissaP() < a ); }
    
    /// true if `a` is below the abscissa of minus end, and thus outside
    bool belowM(const real a) const { return a < abscissaM(); }
    
    /// true if `a` is above the abscissa of minus end
    bool aboveM(const real a) const { return abscissaM() <= a; }

    /// true if `a` is below the abscissa of plus end
    bool belowP(const real a) const { return a <= abscissaP(); }
    
    /// true if `a` is above the abscissa of plus end, and thus outside
    bool aboveP(const real a) const { return abscissaP() < a; }

    /// calculate the domain in which ab is located (near a FiberEnd, or central)
    FiberEnd whichEndDomain(real a, real lambda) const;

    //---------------------
    
    /// displace the ORIGIN of abscissa at distance `a` from the minus end
    void setOrigin(real a) { fnAbscissaM = -a; fnAbscissaP = fnCut * real(nbSegments()) - a; }

    /// signed distance from ORIGIN to minus end (abscissa of minus end)
    real abscissaM() const { return fnAbscissaM; }
    
    /// abscissa of center, midway between minus and plus ends
    real abscissaC() const { return 0.5 * (fnAbscissaM + fnAbscissaP); }

    /// signed distance from ORIGIN to plus end (abscissa of plus end)
    real abscissaP() const { return fnAbscissaP; }

    /// signed distance from ORIGIN to vertex specified with index (or intermediate position)
    real abscissaPoint(const real n) const { return fnAbscissaM + fnCut * n; }

    /// distance from minus end to a vertex specified as index (or intermediate position)
    real distancePointM(const real n) const { return fnCut * n; }

    /// distance from plus end to a vertex specified as index (or intermediate position)
    real distancePointP(const real n) const { return fnCut * ( lastPoint() - n ); }

    /// distance from CENTER to a vertex specified as index (or intermediate position)
    real distancePointC(const real n) const { return fnCut * ( n - 0.5 * lastPoint() ); }

    /// signed distance from the ORIGIN to the specified FiberEnd
    real abscissaEnd(FiberEnd) const;
    
    /// converts distance from the specified FiberEnd, to abscissa from the ORIGIN
    real abscissaFrom(real dis, FiberEnd ref) const;
    
    /// return abscissa specified with `dis, ref, mod`
    real someAbscissa(real dis, FiberEnd ref, int mod, real alpha) const;

    //---------------------

#if ( DIM == 1 )
    /// position at distance `ab` from the minus end
    Vector posM(real ab) const { return Vector(pPos[0]+std::copysign(ab, pPos[1]-pPos[0])); }
#else
    /// position at distance `ab` from the minus end
    Vector posM(real ab) const;
#endif
    /// position of a point specified by abscissa from the ORIGIN
    Vector pos(real ab) const { return posM(ab-fnAbscissaM); }

    /// position of a point specified by abscissa `ab` from reference `ref`
    Vector posFrom(real ab, FiberEnd ref) const { return pos(abscissaFrom(ab, ref)); }

    /// position of the point taken mid-way along the curve
    Vector posMiddle() const { if ( nPoints&1 ) return posPoint(nPoints/2); return midPoint(nPoints/2-1, 0.5); }
    
    /// position of a FiberEnd
    Vector posEnd(FiberEnd) const;
    
    /// position of minus end
    Vector posEndM() const { return Vector(pPos); }

    /// position of plus end
    Vector posEndP() const { return Vector(pPos+DIM*nPoints-DIM); }
    
    /// external force acting on minus end
    Vector netForceEndM() const { return netForce(0); }
    
    /// external force acting on plus end
    Vector netForceEndP() const { return netForce(nPoints-1); }

    //---------------------

#if ( 1 )
    /// unit tangent vector to the fiber within segment [p, p+1]
    /** Using iCut, expected to be the inverse of the distance between vertices */
    Vector dirSegment(index_t p) const { return diffPoints(p) * iCut; }
#else
    /// normalized tangent vector to the fiber within segment [p, p+1]
    Vector dirSegment(index_t p) const { return normalize(diffPoints(p)); }
#endif
#if ( DIM == 1 )
    /// direction at distance `ab` from the minus end
    Vector dirM(real ab) const { return Vector(sign_real(pPos[1]-pPos[0])); }
#else
    /// direction at distance `ab` from the minus end
    Vector dirM(real ab) const;
#endif
    /// normalized tangent vector to the fiber at given abscissa from the origin
    Vector dir(real ab) const { return dirM(ab-fnAbscissaM); }

    /// normalized tangent vector to the fiber at given abscissa from given reference
    Vector dir(real ab, FiberEnd ref) const { return dirM(abscissaFrom(ab, ref)); }
    
    /// normalized tangent vector to the fiber at given end
    Vector dirEnd(FiberEnd) const;
    
    /// normalized tangent vector at minus end, orientated towards the plus end
    Vector dirEndM() const { return dirSegment(0); }
    
    /// normalized tangent vector to the fiber at plus end
    Vector dirEndP() const { return dirSegment(lastSegment()); }

    /// force on the minus end projected on the direction of elongation
    real projectedForceEndM() const;

    /// force on the plus end projected on the direction of elongation
    real projectedForceEndP() const;

    /// dot-product (force at the end of the Fiber).(direction of Fiber growth)
    real projectedForceEnd(FiberEnd) const;
    
    /// direction averaged over the entire length
    Vector direction() const { return normalize(posEndP()-posEndM()); }
    
    /// return updated `normal` that is orthogonal to `d` (used for fake 3D display)
    Vector3 adjustedNormal(Vector3 const& d) const;

    //--------------------- Segmentation / discrete representation
    
    /// set desired segmentation
    void targetSegmentation(real c) { assert_true(c>0); fnSegmentation = c; }
    
    /// return desired segmentation (this is not the length of the segments)
    real targetSegmentation() const { return fnSegmentation; }
    
    /// the current segment length (distance between successive vertices)
    real segmentation() const { return fnCut; }
    
    /// should return 1.0 / segmentation()
    real segmentationInv() const { return iCut; }
    
    /// reinterpolate vertices and adjust fiber to have `ns` segments
    void resegment(index_t ns);
    
    /// automatically select the number of points if needed, and resegment the fiber
    void adjustSegmentation();
    
    /// change the target segmentation, and adjust number of points if needed
    void adjustSegmentation(real);

    /// change all vertices to given array of coordinates
    void getPoints(real const*);
    
    /// restore the distance between successive vertices
    void reshape() { getPoints(pPos); }
    
    //--------------------- Info
    
    /// calculate average and variance of the segment length
    static void computeMeanVar(index_t cnt, real const* ptr, real, real&, real&);

    /// calculate the minimum and maximum segment length, for `cnt` segments
    static void computeMinMax(index_t cnt, real const* ptr, real&, real&);
    
    /// calculate the minimum and maximum segment length
    void segmentMinMax(real& n, real& x) const { computeMinMax(nbSegments(), pPos, n, x); }

    /// curvature calculated at joint `p`, where `0 < p < nbPoints()-1`
    real curvature(index_t p) const;
    
    /// normalized energy associated with bending
    real bendingEnergy0() const;

    /// the cosine of the maximum segment angle: indicate the errors due to curvature
    real minCosine() const;
    
    /// number of joints at which ( cosine(angle) < threshold )
    index_t nbKinks(real threshold = 0) const;
    
    /// calculate intersection between segment `s` and the plane defined by <em> n.pos + a = 0 </em>
    real planarIntersect(index_t s, Vector const& n, const real a) const;

    //--------------------- Growing/Shrinking

    /// increase/decrease length of Fiber by `delta`, at the minus end
    void growM(real delta);
    /// increase/decrease length of Fiber by `delta`, at the plus end
    void growP(real delta);

    /// add a segment of length segmentation() at the minus end
    void addSegmentM();
    /// add a segment of length segmentation() at the plus end
    void addSegmentP();

    /// remove a portion of length `delta` including the minus end
    void cutM(real delta);
    /// remove a portion of length `delta` including the plus end
    void cutP(real delta);

    /// Discard vertices in [ 0, P-1 ] and keep [ P, end ]
    virtual void truncateM(index_t P);
    /// Keep vertices [ 0, P ] and discard the others
    virtual void truncateP(index_t P);

    /// length of polymer made in last timestep, at the minus end (negative for shrinking)
    real freshAssemblyM() const { return cDeltaM; }
    /// length of polymer made in last timestep, at the plus end (negative for shrinking)
    real freshAssemblyP() const { return cDeltaP; }
    
    /// undo last length change at minus end
    void undoGrowM() { growM(-cDeltaM); cDeltaM = 0; }
    /// undo last length change at plus end
    void undoGrowP() { growP(-cDeltaP); cDeltaP = 0; }

    //---------------------

    /// the length of freshly assembled polymer during the last time step
    real freshAssembly(FiberEnd) const;
    
    /// true if the tip `end` has grown ( freshAssembly(which) > 0 )
    bool isGrowing(FiberEnd e) const { return freshAssembly(e) > 0; }
    
    /// true if the tip `end` has shrunk ( freshAssembly(which) < 0 )
    bool isShrinking(FiberEnd e) const { return freshAssembly(e) < 0; }

    /// grow at specified end (plus end or minus end)
    void grow(FiberEnd, real delta);
    
    /// shorten or lengthen Fiber without changing the position of `end`
    void adjustLength(real len, FiberEnd);
    
    /// merge two fibers by attaching given Chain at the plus end of `this`
    void join(Chain const*);

    /// invert polarity (swap PLUS end MINUS ends in place)
    void flipChainPolarity();

    //---------------------

    /// print info such as length and segmentation
    void briefdoc(std::ostream&, real, real, real, real) const;

    /// print info such as length and segmentation
    void document(std::ostream&, real, real, real, real) const;
    
    /// print info such as length and segmentation
    void document(std::ostream&, real const* ptr) const;
    
    /// return string with info such as length and segmentation
    std::string document(real const* ptr) const;

    /// check length and segmentation and write if suspicious
    int checkLength(real const* ptr, std::ostream&, real len) const;

    /// write to Outputter
    void write(Outputter&) const;
    
    /// read from Inputter
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to Outputter
    void writeAngles(Outputter&) const;

    /// read from Inputter, using compact binary mode
    void readAngles(Inputter&, Simul&, ObjectTag);

};


#endif
