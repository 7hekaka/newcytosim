// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef MECA_H
#define MECA_H

#include "dim.h"
#include "array.h"
#include "vector.h"
//#include "sparmat.h"
#include "sparmatsym1.h"
#include "sparmatsymblk.h"
#include "sparmatblk.h"
#include "allocator.h"


class Modulo;
class Mecable;
class Mecapoint;
class Interpolation;
class SimulProp;
class Simul;

/**
 Set to 1 to distribute the Matrix Vector-multiplication in multithreaded code.
 This will select a type for matrix mC specifically built for that purpose.
 This option can only be beneficial if NUM_THREADS > 1
 Do not enable this option for sequential code.
 */
#define PARALLELIZE_MATRIX 0


// known Matrix block types:
class Matrix11;
class Matrix22;
class Matrix33;
class Matrix34;

/// MatrixBlock is an alias to a matrix class of size DIM * DIM
/**
 MatrixBlock is used to update the matrix mC in 'meca_inter.cc',
 and should match the class used for the blocks of mC.
 */
#if ( DIM == 1 )
typedef Matrix11 MatrixBlock;
#elif ( DIM == 2 )
typedef Matrix22 MatrixBlock;
#elif PARALLELIZE_MATRIX
typedef Matrix34 MatrixBlock;
#else
typedef Matrix33 MatrixBlock;
#endif


/// set TRUE to use matrix mB and mC (the traditional way)
/** This option should be 0 if PARALLELIZE_MATRIX == 1 */
#define USE_ISO_MATRIX 0

/**
 Option to allow the user to see Links made by Meca in 'play'.
 This option affects display speed since it requires two calls to setInteractions()
 This option is normally OFF. Only supported with Xcode compilation.
 */
#define DRAW_MECA_LINKS 0


/// A class to calculate the motion of objects in Cytosim
/**
Meca solves the motion of objects defined by points (i.e. Mecable),
using an equation that includes terms for each interaction between Objects,
and also forces that are internal to an object, for instance bending elasticity
for Fibers, and external forces such as confinements.
The equation is formulated using linear-algebra:
 
    d vPTS/dt = mobility * mP * ( Force + mdiffP * vPTS )
 
 with
 
    Force = vBAS + ( mB + mC + mR ) * vPTS
 
 The equation is solved for a small increment of time `time_step`, in the presence
 of Brownian motion, and at low Reynolds number, ie. a regime in which inertial
 forces that are proportional to mass are negligible.
 
 The equation contains `DIM * nbPoints()` degrees of freedom, where `nbPoints()`
 is the total number of points in the system. It contains vectors and matrices.
 The different  terms of the equation are:
 
 - Vector vPTS containing all the Mecable coordinates (x, y, z):
   Fiber, Sphere, Solid and other Mecable. 
 
 - Vector vBAS is of same size as vPTS, and includes the constant part obtained by
   linearization of the forces. It includes for instance the positions of Single,
   calibrated random forces simulating Brownian motion, and also offsets for periodic
   boundary conditions.
 
 - Matrix mB is the isotropic part obtained after linearization of the forces.
   It operates similarly and independently on the different dimension X, Y and Z.
   mB is square of size nbPoints(), symmetric and sparse.
 
 - Matrix mC is the non-isotropic part obtained after linearization of the forces.
   mC is square of size DIM*nbPoints(), symmetric and sparse.
 .
 
 Typically, mB and mC will inherit the stiffness coefficients of the interactions, 
 while vBAS will get forces (stiffness * position). They are set by the member functions
 addLink(), addLongLink(), addSideLink(), addSlidingLink(), etc.

 - mR add the bending elasticity for Mecafil, or other internal forces.
   mR is symmetric of size DIM*nbPoints(), diagonal by blocks, each block corresponding to a Fiber.
 
 - mP applies the projection due to constrained dynamics.
   For Mecafil, this maintains the distance between neighboring points (longitudinal incompressibility). 
   mP is symmetric of size DIM*nbPoints(), diagonal by blocks, each block corresponding to a Fiber.
   mP is not actually calculated as a matrix:
   its application on each block is done by Mecable::projectForces()
 
 - mdiffP is a term coming from the derivative of the projection P.
   It can provide better numerical stability in some situations where the filament are stretched.
   You can however define ADD_PROJECTION_DIFF = 0 in meca.cc to remove mdiffP.
 .
 
 
 Note: All Links have no effect if the given Mecapoint or Interpolation have a 
 point in common, because the matrix elements would not be calcuated correctly 
 in that case. Generally, such interactions are anyway not desirable. It would 
 correspond for example to a link between two point of the same segment, without 
 effect since the segment is straight, or between two successive segments on the
 same Fiber, which at best would fold it in a non-physical way.

 */

class Meca
{
public:
    
    /// verbose level
    int             doNotify;

    /// enables graphical display of all interactions
    bool            drawLinks;
    
    /// used for recording CPU cycles
    mutable unsigned long long cycles_;

private:
    
    /// list of Mecable containing points to simulate
    Array<Mecable*> mecables;
    
    /// local copy of time step
    real            time_step;

    /// flag to indicate that result is available
    int             ready_;
    
    /// total number of points in the system
    size_t          nPoints_;
    
    /// size of the currently allocated memory
    size_t          allocated_;
    
    /// residual threshold used for iterative solver
    real            tolerance_;

    //--------------------------------------------------------------------------
    // Vectors of size DIM * nbPoints()
    
    real*  vPTS;         ///< coordinates of Mecable points
    real*  vSOL;         ///< coordinates after the dynamics has been solved
    real*  vBAS;         ///< part of the force that is independent of positions
    real*  vRND;         ///< vector of Gaussian random numbers
    real*  vRHS;         ///< right hand side of the dynamic system
    real*  vFOR;         ///< the calculated forces, with Brownian components
    real*  vTMP;         ///< intermediate of calculus
    
    //--------------------------------------------------------------------------

    /// working memory allocator for BCGS and GMRES used in solve()
    LinearSolvers::Allocator allocator;
    
    /// secondary memory allocator for GMRES
    LinearSolvers::Allocator temporary;
    
    /// Matrices used for GMRES
    //LinearSolvers::Matrix mH, mV;

private:
#if USE_ISO_MATRIX    
    /// true if the matrix mC is non-zero
    bool   useMatrixC;

    /// isotropic symmetric part of the dynamic
    /** 
     This is a symmetric square matrix of size `nbPoints()`
     It contains terms which have identical coefficients on the X, Y, Z subspaces, such as addLink()
    */
    SparMatSym1  mB;
#endif
    
    /// non-isotropic symmetric part of the dynamic
    /** 
     This is a symmetric square matrix of size `DIM*nbPoints()`
     It contains terms which are different in the X, Y, Z subspaces,
     arising from addSideLink() addSideSlidingLink(), etc.
    */
#if PARALLELIZE_MATRIX
    SparMatBlk     mC;
#else
    SparMatSymBlk  mC;
#endif
    
public:

    /// return address of vector where positions are stored
    real const* addrPTS() const { return vPTS; }

    /// position interpolated from two points in vPTS[]
    Vector  position1(const size_t inx) const;

    /// position interpolated from two points in vPTS[]
    Vector  position2(const size_t inx[2], const real coef[2]) const;

    /// position interpolated from three points in vPTS[]
    Vector  position3(const size_t inx[3], const real coef[3]) const;

    /// position interpolated from four points in vPTS[]
    Vector  position4(const size_t inx[4], const real coef[4]) const;
    
    /// position interpolated from five points in vPTS[]
    Vector  position5(const size_t inx[5], const real coef[5]) const;
    
    /// position interpolated from six points in vPTS[]
    Vector  position6(const size_t inx[6], const real coef[6]) const;

private:
    
    /// add block 'T' to mC at position (i, j)
    void add_block(size_t i, size_t j, MatrixBlock const& T);
 
    /// add block 'alpha*T' to mC at position (i, j)
    void add_block(size_t i, size_t j, real alpha, MatrixBlock const& T);
    
    /// subtract block 'T' to mC at position (i, j)
    void sub_block(size_t i, size_t j, MatrixBlock const& T);

    /// subtract block 'alpha*T' to mC at position (i, j)
    void sub_block(size_t i, size_t j, real alpha, MatrixBlock const& T);

    /// add block 'T' to mC at position (i, i)
    void add_block_diag(size_t i, MatrixBlock const& T);
    
    /// subtract block 'T' to mC at position (i, i)
    void sub_block_diag(size_t i, MatrixBlock const& T);
    
    /// add block 'alpha*T' to mC at position (i, i)
    void add_block_diag(size_t i, real alpha, MatrixBlock const& T);

    /// add value to mB at position (i, j)
    void add_iso(size_t i, size_t j, real val);

    /// add value to mB at position (i, j)
    void sub_iso(size_t i, size_t j, real val);

    /// add value to vBAS at index `i`
    void add_base(size_t i, Vector const&);

    /// add value to vBAS at index `i`
    void add_base(size_t i, Vector const&, real);

    /// sub value to vBAS at index `i`
    void sub_base(size_t i, Vector const&);

private:
    
    /// allocate memory
    void allocate(size_t);
    
    /// release memory
    void release();
    
    /// prepare matrices for 'solve'
    void prepareMatrices();
    
    /// calculate forces for one Mecable
    void multiply1(const Mecable*, const real* X, real* Y) const;
    
    /// implements multiply() followed by precondition() for one Mecable
    void multiply_precondition1(const Mecable*, const real*, real*) const;

    /// implements multiply() followed by precondition() for one Mecable
    void multiply_precondition1(const Mecable*, const real*, real*, real*) const;

    /// calculate the linear part of forces:  Y <- B + ( mB + mC ) * X
    void calculateForces(const real* X, const real* B, real* Y) const;

    /// add forces due to bending elasticity
    void addAllRigidity(const real* X, real* Y) const;

    /// extract the matrix on-diagonal block corresponding to a Mecable
    void getBlock(const Mecable*, real* mat) const;
    
    /// extract the 5-bands symmetric on-diagonal block corresponding to a Mecable
    void getBandedBlock(const Mecable*, real* mat) const;

    /// extract the istropic projection of the on-diagonal block corresponding to a Mecable
    void getIsoBlock(const Mecable*, real* mat) const;

    /// DEBUG: extract the matrix on-diagonal block corresponding to a Mecable using 'multiply()'
    void extractBlock(const Mecable*, real* mat) const;
    
    /// DEBUG: compare `blk` with block extracted using extractBlockSlow()
    void verifyBlock(const Mecable*, const real*);
    
    /// DEBUG: test if `blk` is inverse of block extracted using extractBlockSlow()
    void checkBlock(const Mecable*, const real*);
    
    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondAlt(Mecable*, real*, real*, size_t);

    /// compute all blocks of the preconditionner
    void computePrecondAlt();
    
    /// compute the preconditionner block corresponding to given Mecable
    void renewPreconditionner(Mecable*, int, real*, int*, real*, size_t);
    
    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondBand(Mecable*);
    
    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondIsoS(Mecable*);

    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondIsoP(Mecable*);

    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondFull(Mecable*);
    
    /// compute all blocks of the preconditionner (method=1)
    void computePreconditionner(int, int);
    
    /// compute all blocks of the preconditionner
    void renewPreconditionner(int);

public:
    
    /// constructor
    Meca();
    
    /// destructor
    ~Meca() { release(); }
    
    /// Add a Mecable to the list of objects to be simulated
    void     addMecable(Mecable* p) { mecables.push_back(p); }
    
    /// Number of Mecable
    size_t   nbMecables() const { return mecables.size(); }
    
    /// Number of points in the Mecable that has the most number of points
    size_t   largestMecable() const;

    /// true if system does not contain any object
    bool     empty() const { return nPoints_ == 0; }
    
    /// number of points in the system
    size_t   nbVertices() const { return nPoints_; }
    
    /// Implementation of LinearOperator::size()
    size_t   dimension() const { return DIM * nPoints_; }
    
    /// calculate Y <- M*X, where M is the matrix associated with the system
    void multiply(const real* X, real* Y) const;

    /// apply preconditionner: Y <- P*X (note that X maybe equal to Y)
    void precondition(const real* X, real* Y) const;
    
    /// calculates Y <- P*M*X, for left-sided preconditinning
    void multiply_precondition(const real* X, real* Y) const;
    
    /// equivalent to: multiply(X, T); precondition(T, Y);
    void multiply_precondition(const real* X, real* T, real* Y) const;

    //---------------------- EXPLICIT FORCE ELEMENTS ---------------------------

    /// Add a constant force on Mecapoint
    void addForce(Mecapoint const&, Vector const& force);
    
    /// Add a constant force on Interpolated point
    void addForce(Interpolation const&, Vector const& force);
    
    /// Add a constant force to every points
    void addForceToAll(Vector const& force);
    
    /// Add a torque to the segment indicated by Interpolation
    void addTorque(Interpolation const&, Torque const& torque);
    
    /// Add a torque to constrain the segment to be oriented in direction `dir`
    void addTorqueClamp(Interpolation const&, Vector const& dir, real weight);
    
    /// Add an explicit torque to constrain two segments to be parallel
    void addTorqueExplicit(Interpolation const&, Interpolation const&, real weight);

    /// Add an explicit torque to constrain two segments to an angle defined by (sinus, cosinus)
    void addTorqueExplicit(Interpolation const&, Interpolation const&, real cosinus, real sinus, real weight);
    
    //------------------------- IMPLICIT ELEMENTS ------------------------------

    /// Add a torque to constrain two segments to an angle defined by (sinus, cosinus)
    static MatrixBlock torqueMatrix(real weight, Torque const& axis, real cosinus, real sinus);

    /// this has been replaced by interTorque()
    void addTorquePoliti(Interpolation const&, Interpolation const&, real cosinus, real sinus, real weight);
    
    /// Add a torque to constrain two segments to an angle defined by (sinus, cosinus)
    void addTorque(Interpolation const&, Interpolation const&, MatrixBlock const&, real weight);

    /// Add a torque to constrain two segments to an angle defined by (sinus, cosinus)
    void addTorque(Interpolation const&, Interpolation const&, real cosinus, real sinus, real weight);

    /// Add a torque on 3 points with equilibrium angle defined by (sinus, cosinus)
    void addTorque(Mecapoint const&, Mecapoint const&, Mecapoint const&, MatrixBlock const&, real weight);

    /// Add a torque on 3 points with equilibrium angle defined by (sinus, cosinus)
    void addTorquePlane(Mecapoint const&, Mecapoint const&, Mecapoint const&, Torque const&, real cosinus, real sinus, real weight);

    /// Add a torque on 3 points with equilibrium angle defined by (sinus, cosinus), add LongLink on two points
    void addTorqueLong(Mecapoint const&, Mecapoint const&, Mecapoint const&, MatrixBlock const&, real weight, real len, real weightL);

    /// Link of stiffness `weight` from fixed position
    void addPointClamp(Mecapoint const&, Vector, real weight);
    
    /// Link of stiffness `weight` from fixed position
    void addPointClamp(Interpolation const&, Vector, real weight);
    
    /// Link of stiffness `weight` from fixed position, in the XY plane
    void addPointClampXY(Mecapoint const&, Vector, real weight);

    /// A Hookean force linking all vertices to `cen`
    void addPointClampToAll(Vector const& cen, real weight);

    /// Link of stiffness `weight` and sphere of radius `rad` and center `cen`
    void addSphereClamp(Vector const& off, Mecapoint const&, Vector const& cen, real rad, real weight);
    
    /// Link of stiffness `weight` and sphere of radius `rad` and center `cen`
    void addSphereClamp(Vector const& off, Interpolation const&, Vector const& cen, real rad, real weight);

    /// Link of stiffness `weight` and sphere of radius `rad` and center `cen`
    void addSphereClamp(Mecapoint const&, Vector cen, real rad, real weight);
    
    /// Link of stiffness `weight` and sphere of radius `rad` and center `cen`
    void addSphereClamp(Interpolation const&, Vector  cen, real rad, real weight);
    
    /// Link of stiffness `weight` with cylinder of axis X and radius `len`
    void addCylinderClampX(Mecapoint const&, real rad, real weight);
    
    /// Link of stiffness `weight` with cylinder of axis T and radius `len`
    void addCylinderClampY(Mecapoint const&, real rad, real weight);

    /// Link of stiffness `weight` with cylinder of axis Z and radius `len`
    void addCylinderClampZ(Mecapoint const&, real rad, real weight);
    
    /// Link of stiffness `weight` with cylinder of axis X and radius `len`
    void addCylinderClamp(Mecapoint const&, Vector const&, Vector const&, real rad, real weight);

#if ( DIM == 2 )
    /// Link of stiffness `weight` and resting length `len`, on the side of first segment
    void addSidePointClamp2D(Interpolation const&, Vector, real arm, real weight);
#endif
    /// Link of stiffness `weight` and resting length `len`, on the side of first segment
    void addSidePointClamp3D(Interpolation const&, Vector, Torque const& arm, real weight);

    /// Link of stiffness `weight` with fixed position `pos`, on the side of the segment
    void addSidePointClamp(Interpolation const&, Vector const& pos, real len, real weight);
    
    /// Link of stiffness `weight` with a line defined by `pos` and its tangent `dir`
    void addLineClamp(Mecapoint const&, Vector const& pos, Vector const& dir, real weight);
    
    /// Link of stiffness `weight` with a line defined by `pos` and its tangent `dir`
    void addLineClamp(Interpolation const&, Vector const& pos, Vector const& dir, real weight);
    
    /// Link of stiffness `weight` to coordinate corresponding to `inx`
    void addPlaneClamp(size_t inx, real off, real weight);
    
    /// Link of stiffness `weight` with a plane parallel to YZ offset by `off`
    void addPlaneClampX(Mecapoint const&, real off, real weight);
    
    /// Link of stiffness `weight` with a plane parallel to XZ offset by `off`
    void addPlaneClampY(Mecapoint const&, real off, real weight);
    
    /// Link of stiffness `weight` with a plane parallel to XY offset by `off`
    void addPlaneClampZ(Mecapoint const&, real off, real weight);

    
    /// Link of stiffness `weight` with a plane defined by `pos` and its normal `dir`
    void addPlaneClamp(Mecapoint const&, Vector const& pos, Vector const& dir, real weight);

    /// Link of stiffness `weight` with a plane defined by `pos` and its normal `dir`
    void addPlaneClamp(Interpolation const&, Vector const& pos, Vector const& dir, real weight);

    //------------ ZERO-RESTING LENGTH ELEMENTS LINKING POINTS -----------------
    
    /// Link of stiffness `weight` between two vertices
    void addLink(Mecapoint const&, Mecapoint const&, real weight);
    
    /// Link of stiffness `weight` (use the other one)
    void addLink(Interpolation const&, Mecapoint const&, real weight);
    
    /// Link of stiffness `weight` between a vertex and a interpolated point
    void addLink(Mecapoint const&, Interpolation const&, real weight);
    
    /// Link of stiffness `weight` between two interpolated points
    void addLink(Interpolation const&, Interpolation const&, real weight);
    
    
    /// Link of stiffness `weight` between vertex and interpolated point
    void addLink2(Mecapoint const&, const size_t[], const real[], real weight);
    
    /// Link of stiffness `weight` between vertex and interpolated point
    void addLink3(Mecapoint const&, const size_t[], const real[], real weight);

    /// Link of stiffness `weight` between vertex and interpolated point
    void addLink4(Mecapoint const&, const size_t[], const real[], real weight);
    
    
    /// Link of stiffness `weight` between Interpolation and vertex
    void addLink1(Interpolation const&, size_t, real weight);

    /// Link of stiffness `weight` between Interpolation and interpolated point
    void addLink2(Interpolation const&, const size_t[], const real[], real weight);
    
    /// Link of stiffness `weight` between Interpolation and interpolated point
    void addLink3(Interpolation const&, const size_t[], const real[], real weight);

    /// Link of stiffness `weight` between Interpolation and interpolated point
    void addLink4(Interpolation const&, const size_t[], const real[], real weight);

    //----------------------- ELEMENTS LINKING POINTS --------------------------

    /// Link of stiffness `weight` and resting length `len`
    void addLongLink(Mecapoint const&, Mecapoint const&, real len, real weight);
    
    /// Link of stiffness `weight` and resting length `len`
    void addLongLink(Mecapoint const&, Interpolation const&, real len, real weight);
    
    /// Link of stiffness `weight` and resting length `len`
    void addLongLink(Interpolation const&, Interpolation const&, real len, real weight);

#if ( DIM == 2 )
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment
    void addSideLink2D(Interpolation const&, Mecapoint const&, real arm, real weight);
#endif
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment
    void addSideLink3D(Interpolation const&, Mecapoint const&, Torque const& arm, real weight);

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment
    void addSideLink(Interpolation const&, Mecapoint const&, real arm, real weight);

    
#if ( DIM == 2 )
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment
    void addSideLink2D(Interpolation const&, Interpolation const&, real arm, real weight);
#endif
    /// Link of stiffness `weight`, at distance `arm` on the side of segment supporting first argument
    void addSideLink3D(Interpolation const&, Interpolation const&, Torque const& arm, real weight);
    
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment
    void addSideLink(Interpolation const&, Interpolation const&, real len, real weight);
    
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment
    void testSideLink(Interpolation const&, Mecapoint const&, Torque const& arm, real weight);

#if ( DIM == 2 )
    /// Link of stiffness `weight` and resting length `arm1+arm2`, on the sides of both fibers
    void addSideSideLink2D(Interpolation const&, Interpolation const&, real arm1, real arm2, real weight);
    void addSideSideLink2Dalt(Interpolation const&, Interpolation const&, real arm1, real arm2, real weight);
#endif
    /// Link of stiffness `weight` and resting length `arm1+arm2`, on the sides of both fibers
    void addSideSideLink3D(Interpolation const&, Interpolation const&, Torque const& arm1, Torque const& arm2, real weight);

    /// Link of stiffness `weight` and resting length `arm`, on the sides of both fibers
    void addSideSideLink(Interpolation const&, Interpolation const&, real arm, real weight);

    /// Link of stiffness `weight` and perpendicular to first segment
    void addSlidingLink(Interpolation const&, Mecapoint const&, real weight);
    
    /// Link of stiffness `weight` and perpendicular to first segment
    void addSlidingLink(Interpolation const&, Interpolation const&, real weight);

    
#if ( DIM == 2 )
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink2D(Interpolation const&, Mecapoint const&, real arm, real weight);
    
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLinkS(Interpolation const&, Mecapoint const&, real arm, real weight);
#elif ( DIM >= 3 )
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLinkS(Interpolation const&, Mecapoint const&, Vector const& arm, real weight);
#endif
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink3D(Interpolation const&, Mecapoint const&, Torque const& arm, real weight);

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink(Interpolation const&, Mecapoint const&, real len, real weight);
    
    
#if ( DIM == 2 )
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink2D(Interpolation const&, Interpolation const&, real arm, real weight);
    
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLinkS(Interpolation const&, Interpolation const&, real arm, real weight);
#elif ( DIM >= 3 )
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLinkS(Interpolation const&, Interpolation const&, Torque const& arm, real weight);
#endif
    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink3D(Interpolation const&, Interpolation const&, Torque const&, real weight);

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink(Interpolation const&, Interpolation const&, real len, real weight);
    
    
    /// Create a 3-way link with given weights on each branch
    void addTriLink(Interpolation const& pt1, real w1, Interpolation const& pt2, real w2, Interpolation const& pt3, real w3);

    /// Linearized Coulomb repulsive force (experimental)
    void addCoulomb(Mecapoint const&, Mecapoint const&, real weight);
    
    //-------------------------- COMPUTING METHODS -----------------------------

    /// Allocate the memory necessary to solve(). This must be called after the last add()
    void prepare(Simul const*);
    
    /// Calculate motion of all Mecables in the system; returns number of step of the iterative solver
    size_t solve(SimulProp const*, unsigned precondition);
    
    /// transfer newly calculated point coordinates back to Mecables
    void apply();

    /// calculate Forces on Mecables and Lagrange multipliers for Fiber, without thermal motion
    void computeForces();
    
    //----------------------- EXPORT/DEBUG FUNCTIONS ---------------------------
    
    /// set Mecable:flag() according to connectivity defined by matrix elements
    void flagClusters() const;
    
    /// Count number of non-zero entries in the entire system
    size_t nbNonZeros(real threshold) const;

    /// Extract the complete dynamic matrix in column-major format in a C-array
    void getMatrix(size_t, real * matrix) const;
    
    /// Save complete matrix in Matrix Market format
    void saveMatrix(FILE *, real threshold) const;
    
    /// Save right-hand-side vector
    void saveRHS(FILE *) const;
    
    /// Output vectors and matrices, in a format that can be imported in MATLAB
    void saveSystem(const char dirname[]) const;

    /// Save complete matrix in binary format
    void dumpMatrix(FILE *) const;
    
    /// Save elasticity matrix in binary format
    void dumpElasticity(FILE *) const;
    
    /// Save mobility/projection matrix in binary format
    void dumpMobility(FILE *) const;
    
    /// Save preconditionner in binary format
    void dumpPreconditionner(FILE *) const;
    
    /// Save drag coefficients associated with each degree of freedom in binary format
    void dumpDrag(FILE *) const;
    
    /// Save the object ID associated with each degree of freedom
    void dumpObjectID(FILE *) const;
    
    /// Output vectors and matrices, in a format that can be imported in MATLAB
    void dump() const;
    
    /// Output vectors and matrices, in a format that can be imported in MATLAB
    void dump(const char dirname[]) const;

    /// Output vectors and matrices in various files (for debugging)
    void dumpSparse();

};

#endif

