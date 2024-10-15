// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef MECABLE_H
#define MECABLE_H

#include "dim.h"
#include "object.h"
#include "buddy.h"
#include "matfull.h"
#include "cymdef.h"

class Meca;
class Space;


/**
 Add correction terms to the projection in constrainted dynamics
 The effect is to stabilize fibers under traction, at some modest CPU cost.
 Set ADD_PROJECTION_DIFF = 7 to enable validation codes in makeProjectionDiff()
*/
#define ADD_PROJECTION_DIFF 1


/// Can be simulated using a Meca.
/**
 A Mecable is an Object made of points that can can be simulated in a Meca.
 
 Mecable defines an interface that is implemented in Bead, Fiber, Sphere and Solid.
 It derives from Movable, and thus can be translated or rotated.
 Mecable is also a Buddy, and can thus be part of an Organizer.
 
 This is one of the fundamental class of Cytosim, working closely with Meca.
 
 Acquired PointSet on 20/01/2018.
 */
class Mecable : public Object, public Buddy
{
public:
    
    /// to save memory, SIZE_T could be defined here to use only 2 bytes.
    /* The limit imposed on the size of the Mecable is then 65535 vertices */
    typedef unsigned short SIZE_T;

protected:

    /// array of size DIM*pAllocated contains DIM*nPoints coordinates
    /**
     The coordinates are arranged as follows:
     X1,       X2,       etc. for DIM==1
     X1 Y1,    X2 Y2,    etc. for DIM==2
     X1 Y1 Z1, X2 Y2 Z2, etc. for DIM==3
     */
    real * pPos;

private:
    
    /// Array containing force-coordinates which is allocated in Meca
    real const* pForce;

    /// Matrix block used for preconditionning in Meca
    real * pBlock;
    
    /// pointer to allocated memory for pivot indices used in matrix factorization
    int * pPivot;
    
    /// Index that Object coordinates occupy in the matrices and vectors of Meca
    index_t pIndex;
    
    /// Allocated size of pBlock[], capable of size^2
    index_t pBlockAlc;

protected:

    /// Number of points in the Mecable
    SIZE_T nPoints;

private:

    /// Currently allocated size: pPos[] can hold DIM*pAllocated scalars
    SIZE_T pAllocated;
    
    /// type of preconditionner
    SIZE_T pBlockType;
    
public:

    /// The constructor resets the pointers to memory
    Mecable();
    
    /// Destructor should release memory
    virtual ~Mecable() { release(); }

    /// Copy constructor
    Mecable(Mecable const&);
    
    /// Copy assignment
    Mecable& operator = (Mecable const&);
    
    //--------------------------------------------------------------------------
    
    /// Set the number of points of the object
    void setNbPoints(const index_t n);
    
    /// Returns number of points
    index_t nbPoints() const { return nPoints; }
    
    /// Index of the last point = nbPoints() - 1
    index_t lastPoint() const { return nPoints - 1; }
    
    /// size currently allocated
    index_t allocated() const { return pAllocated; }
    
    /// Number of distance constraints applied to the movements of vertices
    virtual index_t nbConstraints() const { return 0; }

    //--------------------------------------------------------------------------
    
    /// Position of vertex number 'p' (indices starting at zero)
    Vector posPoint(index_t p) const { assert_true(pPos && p<nPoints); return Vector(pPos+DIM*p); }
    
    /// Position of point 'p' of the object
    /** this is identical to posPoint(), it exists for historical reasons*/
    Vector posP(index_t p)     const { return Vector(pPos+DIM*p); }
    
    /// modifiable address of coordinate array
    real * addrPoints() { return pPos; }

    /// Address of coordinate array
    const real * addrPoints() const { return pPos; }

    /// Address of point `p`
    const real * addrPoint(index_t p) const { assert_true(p<nPoints); return pPos + DIM*p; }

    /// Set position of point `i` to `x`
    void setPoint(index_t i, Vector const& x) { assert_true(i<nPoints); x.store(pPos+DIM*i); }
    
    /// Shift point at index `i` by `x`
    void movePoint(index_t i, Vector const& x) { assert_true(i<nPoints); x.add_to(pPos+DIM*i); }
    
    /// replace current coordinates by values from the given array
    virtual void getPoints(real const*);
    
    /// Set to `n_pts` points copied from `pts[]`
    void setPoints(const real pts[], index_t n_pts);
    
    /// copy current vertex coordinates to given array, assuming it is suitably allocated
    void putPoints(real*) const;

    /// Copy current vertex coordinates to `ptr[]`, already allocated to hold `sup` scalars
    void putPoints(float ptr[], index_t sup) const;

    /// Add a point and expand the object, returning the array index that was used
    index_t addPoint(Vector const& w);
    
    /// Remove `nbp` points starting from index `inx`
    void removePoints(index_t inx, index_t nbp);
    
    /// Remove all points
    void clearPoints() { nPoints = 0; }
    
    /// Shift `nbp` points starting from index `inx`
    void shiftPoints(index_t inx, index_t nbp);
    
    /// Remove all points with indices [ 0, p-1 ], keep [ p, nbPoints() ]
    virtual void truncateM(index_t p);
    
    /// Keep points [ 0, p ], remove other points
    virtual void truncateP(index_t p);
    
    /// Set all coordinates to zero (nicer for debug/testing)
    void resetPoints();
    
    /// Add random noise uniformly to all coordinate (used for testing purposes)
    void addNoise(real amount);
    
    /// calculate first and second momentum of point coordinates
    void calculateMomentum(Vector& avg, Vector& dev);
    
    //--------------------------------------------------------------------------
    
    /// Difference of two consecutive points: (P+1) - (P)
    Vector diffPoints(const index_t P) const
    {
        assert_true( P+1 < nPoints );
        Vector vec;
        vec.load_diff(pPos+DIM*P);
        return vec;
    }
    
    /// Difference of two points = Q - P = vector PQ
    Vector diffPoints(const index_t P, const index_t Q) const
    {
        assert_true( P < nPoints );
        assert_true( Q < nPoints );
        Vector vec;
        vec.load_diff(pPos+DIM*Q, pPos+DIM*P);
        return vec;
    }
    
    /// intermediate position between P and Q=P+1 : P + A * ( Q - P )
    /** returns P if A=0; returns P+1 if A=1, and returns middle of segment if A = 0.5*/
    Vector midPoint(const index_t P, const real A) const
    {
        assert_true( P < nPoints );
        assert_true( A == 0 || P+1 < nPoints );
        //assert_true( 0 <= A && A <= 1 );
        return Vector::interpolated(pPos+DIM*P, A, pPos+DIM*P+DIM);
    }
    
    /// middle position between P and P+1 = 0.5 * ( P + Q )
    Vector midPoint(const index_t P) const
    {
        assert_true( P+1 < nPoints );
        return 0.5 * ( Vector(pPos+DIM*P) + Vector(pPos+DIM*P+DIM) );
    }

    /// Calculate intermediate position = P + a * ( Q - P )
    Vector interpolatePoints(const index_t P, const index_t Q, const real A) const
    {
        assert_true( P < nPoints );
        assert_true( Q < nPoints );
        //assert_true( 0 <= A && A <= 1 );
        return Vector::interpolated(pPos+DIM*P, A, pPos+DIM*Q);
    }
    
    /// interpolate 'rank' points starting from 'ref' with coefficients
    Vector interpolatePoints(index_t ref, const real coef[], index_t rank) const;

    //--------------------------------------------------------------------------

    /// Allocate memory to store given number of vertices and extra vectors
    real * allocateMemory(index_t, index_t);

    /// Allocate memory to store given number of vertices
    virtual void allocateMecable(index_t n) { allocateMemory(n, 0); }
    
    /// free allocated memory
    void release();
    
    /// prepare the Mecable to solve the mechanics in Meca
    /**
     This should prepare necessary variables to solve the system:
     - set bending elasticity coefficients, for addRigidity() to work properly
     - set drag mobility, for projectForces() to work,
     - set matrix/variables necessary for constrained dynamics
     .
     */
    virtual void prepareMecable() = 0;

    /// Calculate the mobility coefficient
    virtual void setDragCoefficient() = 0;
    
    /// The total drag coefficient of the object ( force = drag * speed )
    virtual real dragCoefficient() const = 0;

    /// The mobility of a model vertex ( speed = mobility * point_force )
    virtual real pointMobility() const = 0;
    
    /// Add Brownian noise terms to a force vector (alpha = kT / timestep)
    virtual real addBrownianForces(real const* fce, real alpha, real* rhs) const { return INFINITY; }
    
    /// add the interactions (for example due to confinements)
    virtual void setInteractions(Meca&) const {}
    
    //--------------------------------------------------------------------------
    
    /// Store the index where coordinates are located in Meca
    void setIndex(index_t i) { pIndex = i; }
    
    /// Index of the first vertex in the isotropic matrix (Meca::mISO)
    /**
     The coordinates are stored consecutively in Meca's vectors:
     [X1, Y1, Z1] are stored startig at DIM*matIndex()
     [X2, Y2, Z2] is store at DIM*(1+matIndex()) ...
     The index can be used directly to address the isotropic matrix (Meca::mISO)
     and should be multiplied by DIM to address the general matrix (Meca::mFUL)
     */
    index_t matIndex() const { return pIndex; }
    
    /// set size of preconditionner block, allocating memory for 'alc' scalars
    void blockSize(index_t, index_t alc, index_t pivot);
    
    /// Returns allocated size of preconditionner block
    index_t blockLimit() const { return pBlockAlc; }

    /// True if preconditionner block is 'in use'
    SIZE_T blockType()  const { return pBlockType; }

    /// Returns address of memory allocated for preconditionning
    real * pblock()     const { return pBlock; }
    
    /// Returns address of memory available to store pivoting indices
    int * pivot() const { return pPivot; }
 
    /// Type of block: 0=identity; 1=full; 2=band; 3=custom
    void blockType(SIZE_T t) { pBlockType = t; }

    //--------------------------------------------------------------------------
    
    /// returns the force on point `p` calculated at the previous Meca's solve
    Vector netForce(const index_t p) const;
    
    /// replace current forces by the ones provided as argument
    virtual void getForces(const real* ptr) { pForce = ptr; }
    
    //--------------------------------------------------------------------------

    /// return coefficient for filament's bending elasticity
    virtual real jointRigidity() const { return 0; }

    /// Add bending elasticity terms Y <- Y + Rigidity * X
    /**
        Rigidity can be any force acting internally to the objects
     in particular, the bending elasticity of Fibers.
     This function is used in Meca to calculate force from positions.
     */
    virtual void addRigidity(const real* X, real* Y) const {}

    /// Calculate speeds for given forces: Y <- forces(X)
    /**
     The function calculates the 'legal' forces with constraints applied.
     It may or may not scale by the object's mobility coefficient, but one may
     derive the speeds in conjunction with `leftoverMobility()`:
     
         speed = leftoverMobility() * projectForces(forces)
     
     Note that:
     - The input `X` and output `Y` must be vectors of size `DIM * nPoints`
     - `X` and `Y` may point to the same address
     
     The default implementation ( Y <- 0 ) makes the object immobile
     */
    virtual void projectForces(const real* X, real* Y) const { zero_real(DIM*nPoints, Y); }

    /// Return drag coefficient that was not applied by projectForces()
    virtual real leftoverMobility() const { return 1.0; }

    //--------------------------------------------------------------------------
    
#if ADD_PROJECTION_DIFF
    /// true if projectionDiff is used
    bool useProjectionDiff;

    /// set terms derived from the Projection operator, from the given forces
    virtual void makeProjectionDiff(const real* force) {}
    
    /// add terms from projection correction terms: Y <- Y + diffP * X;
    virtual void addProjectionDiff(const real* X, real* Y) const {}
    
    /// add terms from projection correction matrix: mat <- diffP
    virtual void addProjectionDiff(real* mat) const {}
    
    /// true if addProjectionDiff() does something
    bool hasProjectionDiff() const { return useProjectionDiff; }
#endif
    
    //--------------------------------------------------------------------------
    //           Position-related functions derived from Movable
    //--------------------------------------------------------------------------
    
    /// Position of center of gravity
    virtual Vector position() const;
    
    /// Mecable accepts translation and rotation
    virtual int mobile() const { return 3; }
    
    /// Translate object (moves all the points by the same vector)
    virtual void translate(Vector const&);
    
    /// Rotate object by given rotation
    virtual void rotate(Rotation const&);
    
    /// bring object to centered image using periodic boundary conditions
    virtual void foldPosition(Modulo const*);
    
    /// true if all points are inside Space
    bool allPointsInside(Space const*) const;

    //--------------------------------------------------------------------------
    
    /// Write to file
    void write(Outputter&) const;
    
    /// Read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// Human friendly ouput
    void print(std::ostream&, real const*) const;
    
    /// return (validated) index encoded in `str`
    index_t point_index(std::string const& str) const;

    /// check for NaNs in the position vector
    int invalid() const;
};


/// output operator:
std::ostream& operator << (std::ostream& os, Mecable const&);

#endif
