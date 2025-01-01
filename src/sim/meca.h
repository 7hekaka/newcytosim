// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef MECA_H
#define MECA_H

#include "dim.h"
#include "array.h"
#include "vector.h"
//#include "sparmat.h"
#include "sparmatsym1.h"
#include "sparmatsymblk.h"
#include "sparmatsymblkdiag.h"
#include "sparmatblk.h"
#include "allocator.h"
#include "point_grid.h"
#include "locus_grid.h"


class Modulo;
class Mecable;
class Mecapoint;
class Interpolation;
class Interpolation4;
class FiberSegment;
class SimulProp;
class Simul;


/// set 1 to use matrix mISO and mFUL (the traditional way)
/**
 USE_ISO_MATRIX should not affect the results in any way, but the speed of execution.
 This makes the simulation a lot faster on isotropic systems (eg. self.cym), that
 are using isotropic force elements such as the Hookean link (Meca::addLink().
 It is useless for a purely non-isotropic system and causes a bit of overhead.
 */
#define USE_ISO_MATRIX 1


typedef SparMatSymBlkDiag BlockMatrixType;


/// MatrixBlock is an alias to a matrix class of size DIM * DIM
typedef BlockMatrixType::Block MatrixBlock;

/**
 Option to allow 'play' to display Meca links graphically.
 This option affects sim and play speed with some code added to all Meca::addLink()
 This option is normally OFF.
 */
#define DRAW_MECA_LINKS 0

// this should be defined to use uniform_flow
#define NEW_CYTOPLASMIC_FLOW 0


/// A class to calculate the motion of objects in Cytosim
/**
Meca solves the motion of objects defined by points (i.e. Mecable),
using an equation that includes terms for each interaction between Objects,
and also forces that are internal to an object, for instance bending elasticity
for Fibers, and external forces such as confinements.
The equation is formulated using linear-algebra:
 
    d vPTS/dt = mobility * mP * ( Force + mDiffP * vPTS )
 
 with
 
    Force = vBAS + ( mISO + mFUL + mR ) * vPTS
 
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
 
 - Matrix mISO is the isotropic part obtained after linearization of the forces.
   It operates similarly and independently on the different dimension X, Y and Z.
   mISO is square of size nbPoints(), symmetric and sparse.
 
 - Matrix mFUL is the non-isotropic part obtained after linearization of the forces.
   mFUL is square of size DIM*nbPoints(), symmetric and sparse.
 .
 
 Typically, mISO and mFUL will inherit the stiffness coefficients of the interactions, 
 while vBAS will get forces (stiffness * position). They are set by the member functions
 addLink(), addLongLink(), addSideLink(), addSlidingLink(), etc.

 - mR add the bending elasticity for Mecafil, or other internal forces.
   mR is symmetric of size DIM*nbPoints(), diagonal by blocks, each block corresponding to a Fiber.
 
 - mP applies the projection due to constrained dynamics.
   For Mecafil, this maintains the distance between neighboring points (longitudinal incompressibility). 
   mP is symmetric of size DIM*nbPoints(), diagonal by blocks, each block corresponding to a Fiber.
   mP is not actually calculated as a matrix:
   its application on each block is done by Mecable::projectForces()
 
 - mDiffP is a term coming from the derivative of the projection P.
   It can provide better numerical stability in some situations where the filament are stretched.
   You can however define ADD_PROJECTION_DIFF=0 to remove mDiffP.
 .
 
 
 Note: addLinks() calls have no effect if the given Mecapoint or Interpolation have a
 point in common, because the matrix elements would not be calcuated correctly 
 in that case. Generally, such interactions are anyway not desirable. It would 
 correspond for example to a link between two point of the same segment, without 
 effect since the segment is straight, or between two successive segments on the
 same Fiber, which at best would fold it in a non-physical way.

 */

class Meca
{
    friend class Simul;
    friend class Display;
    
private:
    
    /// time step for Brownian Mechanics = copy of simul:time_step
    real tau_;
    
    /// magnitude of Brownian motion = 2 * simul:kT / simul:time_step
    real alpha_;
    
    /// accepted residual threshold when solving linear system
    real tolerance_;
    
    /// total number of points in the system
    index_t nPoints_;
    
    /// size of the currently allocated memory
    index_t allocated_;

    /// number of preconditionner blocks that could not be factorized
    index_t bump_;
    
    /// flag to indicate that result is available
    int ready_;
    
    /// flag to include steric interactions
    int steric_;
    
    /// preconditionning mode
    int precond_;
    
    /// verbosity level
    int verbose_;
    
    /// list of Mecable containing points to simulate
    Array<Mecable*> mecables;

#if NEW_CYTOPLASMIC_FLOW
    Vector uniform_flow_dt_;
#endif
    //--------------------------------------------------------------------------
    // Vectors of size DIM * nbPoints()
    
    real * vPTS;         ///< coordinates of Mecable points
    real * vSOL;         ///< coordinates after the dynamics has been solved
    real * vBAS;         ///< part of the force that is independent of positions
    real * vRHS;         ///< right hand side of the dynamic system
    real * vFOR;         ///< the calculated forces, with Brownian components
    real * vTMP;         ///< intermediate of calculus
    
    //--------------------------------------------------------------------------

    /// working memory allocator for BCGS
    LinearSolvers::Allocator allocator_;
    
    /// grid used for steric interaction between Fiber/Solid/Bead/Sphere
    PointGrid pointGrid;
    
    /// alternative grid used for steric interaction, instead of pointGrid
    LocusGrid locusGrid;

    /// Matrices used for GMRES
    //LinearSolvers::Matrix mH, mV;

    /// address of force base (for pycytosim)
    real * base() { return vBAS; }
    
    /// address of points vector (for pycytosim)
    real * points() { return vPTS; }
    
    /// address of force vector (for pycytosim)
    real * force() { return vFOR; }

public:
    
    /// record of time (CPU cycles)
    mutable unsigned long long cycles_;

    /// verbose level
    int doNotify;

    /// enables graphical display of all interactions
    int drawLinks;

private:
#if USE_ISO_MATRIX    
    /// true if the matrix mFUL is non-zero
    bool useFullMatrix;

    /// isotropic symmetric part of the dynamic
    /** 
     This is a symmetric square matrix of size `nbPoints()`, acting
     identically and separately on the X, Y, Z subspaces, such as addLink()
    */
    SparMatSym1 mISO;
#endif
    
    /// non-isotropic symmetric part of the dynamic
    /** 
     This is a symmetric square matrix of size `DIM * nbPoints()`
     It contains terms which are different in the X, Y, Z subspaces,
     arising from addSideLink() addSideSlidingLink(), etc.
    */
    BlockMatrixType mFUL;
    
public:

    /// return address of vector where positions are stored
    real const* addrPTS() const { return vPTS; }

private:
    
    /// add block 'T' to mFUL at position (i, j)
    void add_block(index_t i, index_t j, MatrixBlock const& T);
 
    /// add block 'alpha*T' to mFUL at position (i, j)
    void add_block(index_t i, index_t j, real alpha, MatrixBlock const& T);
    
    /// add block '-T' to mFUL at position (i, j)
    void sub_block(index_t i, index_t j, MatrixBlock const& T);

    /// add block '-alpha*T' to mFUL at position (i, j)
    void sub_block(index_t i, index_t j, real alpha, MatrixBlock const& T);

    /// add block 'T' to mFUL at position (i, i)
    void add_block_diag(index_t i, MatrixBlock const& T);
    
    /// add block '-T' to mFUL at position (i, i)
    void sub_block_diag(index_t i, MatrixBlock const& T);
    
    /// add block 'alpha*T' to mFUL at position (i, i)
    void add_block_diag(index_t i, real alpha, MatrixBlock const& T);
    
    /// add block 'alpha * (T + dia * Id)' to mFUL at position (i, i)
    void add_block_diag(index_t i, real alpha, MatrixBlock const& T, real dia);

    /// add isotropic stiffness at position (i, j)
    void add_iso(index_t i, index_t j, real val);
    
    /// add isotropic stiffness on diagonal at position (i, i)
    void add_iso_diag(index_t i, real val);

    /// subtract isotropic stiffness at position (i, j)
    void sub_iso(index_t i, index_t j, real val);
    
    /// subtract isotropic stiffness on diagonal at position (i, i)
    void sub_iso_diag(index_t i, real val);

    /// add vector to vBAS at index `i`
    void add_base(index_t, Vector const&) const;

    /// add scaled vector to vBAS at index `i`
    void add_base(index_t, Vector const&, real) const;

    /// sub vector to vBAS at index `i`
    void sub_base(index_t, Vector const&) const;

    /// sub vector to vBAS at index `i`
    void sub_base(index_t, Vector const&, real) const;

private:
    
    /// allocate memory
    void allocate(index_t);
    
    /// release memory
    void release();
    
    /// prepare matrices for 'solve'
    void prepareMatrices();
    
    /// calculate forces for one Mecable
    void multiply1(const Mecable*, const real* X, real* Y) const;

    /// calculate the linear part of forces:  Y <- B + ( mISO + mFUL ) * X
    void calculateForces(const real* X, const real* B, real* Y) const;

    /// redraw new values of the noise, and update system's RHS
    void renewBrownianForces();
    
    /// add forces due to bending elasticity
    void addAllRigidity(const real* X, real* Y) const;
    
    /// extract the matrix on-diagonal block corresponding to a Mecable
    void getBandedBlock(const Mecable*, real* mat, index_t ldd, index_t rank) const;

    /// extract the matrix on-diagonal block corresponding to a Mecable
    void getHalfBlock(const Mecable*, real* mat) const;
 
    /// extract the matrix on-diagonal block corresponding to a Mecable
    void getFullBlock(const Mecable*, real* mat) const;

    /// extract the 5-bands symmetric on-diagonal block corresponding to a Mecable
    void getIsoBandedBlock(const Mecable*, real* mat, index_t kd, index_t ldd) const;

    /// extract the istropic projection of the on-diagonal block corresponding to a Mecable
    void getIsoBlock(const Mecable*, real* mat) const;

    /// DEBUG: extract the matrix on-diagonal block corresponding to a Mecable using 'multiply()'
    void extractBlock(const Mecable*, real* mat) const;
    
    /// DEBUG: compare `blk` with block returned by extractBlock()
    void verifyBlock(const Mecable*, const real* blk);
    
    /// DEBUG: test if `blk` is inverse of block returned by extractBlock()
    void checkBlock(const Mecable*, const real* blk);
    
    /// compute the preconditionner block corresponding to given Mecable
    void renewPreconditionner(Mecable*, int, real*, int*, real*, index_t);
    
    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondIsoB(Mecable*, real*);
    
    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondIsoS(Mecable*);

    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondIsoP(Mecable*, real*);
    
    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondBand(Mecable*, real*);

    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondHalf(Mecable*, real*);

    /// compute the preconditionner block corresponding to given Mecable
    void computePrecondFull(Mecable*, real*);
    
    /// compute all blocks of the preconditionner (method=1)
    void computePreconditionner();
    
    /// compute all blocks of the preconditionner
    void renewPreconditionner(int);

public:
    
    /// constructor
    Meca();
    
    /// destructor
    ~Meca() { release(); }
    
    /// Add a Mecable to the list of objects to be simulated
    void addMecable(Mecable* p) { mecables.push_back(p); }
    
    /// true if system does not contain any object
    bool empty() const { return nPoints_ == 0; }

    /// Number of points in the Mecable that has the most number of points
    index_t largestMecable() const;
    
    /// Number of points in the Mecable that has the least number of points
    index_t smallestMecable() const;

    /// Number of distance constraints applied to the movements of vertices
    index_t nbConstraints() const;
    
    /// Number of Mecable
    index_t nbMecables() const { return (index_t)mecables.size(); }

    /// number of points in the system
    index_t nbVertices() const { return nPoints_; }
    
    /// Implementation of LinearOperator::size()
    index_t dimension() const { return DIM * nPoints_; }
    
    /// total allocated memory size for preconditionner
    size_t preconditionnerSize() const;
    
    /// returns the norm of `RHS - MATRIX * SOL`
    real residualNorm() const;

    /// calculate Y <- M*X, where M is the matrix associated with the system
    void multiply(const real* X, real* Y) const;

    /// apply preconditionner: Y <- P*X (this works even if X == Y)
    void precondition(const real* X, real* Y) const;

    //---------------------- EXPLICIT FORCE ELEMENTS ---------------------------

    /// Add a constant force on Mecapoint
    void addForce(Mecapoint const&, Vector const& force);
    
    /// Add a constant force on Mecapoint
    void addForce(Mecable const*, index_t inx, Vector const& force);
    
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

    /// Add an explicit torque to constrain two segments to an angle defined by ang = (cosine, sine)
    void addTorqueExplicit(Interpolation const&, Interpolation const&, Vector2 const& ang, real weight);
    
    //------------------------- IMPLICIT ELEMENTS ------------------------------

    /// Add a torque to constrain two segments to an angle defined by ang = (cosine, sine)
    static MatrixBlock torqueMatrix(real weight, Torque const& axis, Vector2 const& ang);

    /// old code that has been replaced by interTorque()
    void addTorquePoliti(Interpolation const&, Interpolation const&, Vector2 const& ang, real weight);
    
    /// Add a torque to constrain two segments to be aligned
    void addTorque(Interpolation const&, Interpolation const&, real weight);

    /// Add a torque to constrain two segments to an angle defined by MatrixBlock
    void addTorque(Interpolation const&, Interpolation const&, MatrixBlock const&, real weight);

    /// Add a torque to constrain two segments to an angle defined by ang = (cosine, sine)
    void addTorque(Interpolation const&, Interpolation const&, Vector2 const& ang, real weight);
    
    /// Add a 'bending elasticity' torque on 3 points 
    void addTorque3(Mecapoint const&, Mecapoint const&, Mecapoint const&, real scale, real weight);

    /// Add a torque on 3 points with equilibrium angle defined by MatrixBlock
    void addTorque3(Mecapoint const&, Mecapoint const&, Mecapoint const&, MatrixBlock const&, real weight);

    /// Add a torque on 3 points with equilibrium angle defined by ang = (cosine, sine)
    void addTorque3Plane(Mecapoint const&, Mecapoint const&, Mecapoint const&, Torque const&, Vector2 const& ang, real weight);

    /// Add a torque on 3 points with equilibrium angle defined by ang = (cosine, sine), add LongLink on two points
    void addTorque3Long(Mecapoint const&, Mecapoint const&, Mecapoint const&, MatrixBlock const&, real weight, real len, real weightL);
    
    /// Add a torque on 4 points to align AB with CD
    void addTorque4(index_t indexA, index_t indexB, index_t indexC, index_t indexD, real weight);
    
    /// Add a torque on 4 points to align AB with CD
    void addTorque4(Mecapoint const&, Mecapoint const&, real weight);

    /// Link of stiffness `weight` from fixed position
    void addPointClamp(Mecapoint const&, Vector, real weight);
    
    /// Link of stiffness `weight` from fixed position
    void addPointClamp(Interpolation const&, Vector, real weight);

    /// Link of stiffness `weight` from fixed position, in the X plane only
    void addPointClampX(Mecapoint const&, real x_pos, real weight);

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
    
    /// Link of stiffness `weight` with cylinder of axis X and radius `rad`
    void addCylinderClampX(Mecapoint const&, real rad, real weight);
    
    /// Link of stiffness `weight` with cylinder of axis T and radius `rad`
    void addCylinderClampY(Mecapoint const&, real rad, real weight);

    /// Link of stiffness `weight` with cylinder of axis Z and radius `rad`
    void addCylinderClampZ(Mecapoint const&, real rad, real weight);
    
    /// Link of stiffness `weight` with cylinder of axis X and radius `rad`
    void addCylinderClamp(Mecapoint const&, Vector const&, Vector const&, real rad, real weight);

#if ( DIM == 2 )
    /// Link of stiffness `weight` and resting length `arm`, on the side of first segment
    void addSidePointClamp2D(Interpolation const&, Vector, real arm, real weight);
#endif
    /// Link of stiffness `weight` and resting length `arm`, on the side of first segment
    void addSidePointClamp3D(Interpolation const&, Vector, Torque const& arm, real weight);

    /// Link of stiffness `weight` with fixed position `pos`, on the side of the segment
    void addSidePointClamp(Interpolation const&, Vector pos, real len, real weight);
    
    /// Link of stiffness `weight` with the X-axis
    void addLineClampX(Mecapoint const&, real weight);

    /// Link of stiffness `weight` with a line defined by `pos` and its tangent `dir`
    void addLineClamp(Mecapoint const&, Vector const& pos, Vector const& dir, real weight);
    
    /// Link of stiffness `weight` with a line defined by `pos` and its tangent `dir`
    void addLineClamp(Interpolation const&, Vector const& pos, Vector const& dir, real weight);

    
    /// Link of stiffness `weight` orthogonal to one of the principal plane and offset by `off`
    void addPlaneClampXYZ(Mecapoint const& P, int xyz, real off, real weight);

    /// Link of stiffness `weight` with a plane parallel to YZ and offset by `off`
    void addPlaneClampX(Mecapoint const& P, real off, real weight) { addPlaneClampXYZ(P, 0, off, weight); }
    
    /// Link of stiffness `weight` with a plane parallel to XZ and offset by `off`
    void addPlaneClampY(Mecapoint const& P, real off, real weight) { addPlaneClampXYZ(P, 1, off, weight); }
    
    /// Link of stiffness `weight` with a plane parallel to XY and offset by `off`
    void addPlaneClampZ(Mecapoint const& P, real off, real weight) { addPlaneClampXYZ(P, 2, off, weight); }

    
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
    void addLink2(Mecapoint const&, const index_t, real, real, real weight);
    
    /// Link of stiffness `weight` between vertex and interpolated point
    void addLink3(Mecapoint const&, const index_t, real, real, real, real weight);

    /// Link of stiffness `weight` between vertex and interpolated point
    void addLink4(Mecapoint const&, const index_t, real, real, real, real, real weight);
    
    
    /// Link of stiffness `weight` between Interpolation and vertex
    void addLink1(Interpolation const&, index_t, real weight);

    /// Link of stiffness `weight` between Interpolation and interpolated point
    void addLink2(Interpolation const&, const index_t, real, real, real weight);
    
    /// Link of stiffness `weight` between Interpolation and interpolated point
    void addLink3(Interpolation const&, const index_t, real, real, real, real weight);

    /// Link of stiffness `weight` between Interpolation and interpolated point
    void addLink4(Interpolation const&, const index_t, real, real, real, real, real weight);

    //----------------------- ELEMENTS LINKING POINTS --------------------------

    /// Link of stiffness `weight` and resting length `len`, tweaked version
    void addLongLink1(Mecapoint const&, Mecapoint const&, Vector const&, real ab2, real len, real weight);

    /// Link of stiffness `weight` and resting length `len`, tweaked version
    void addLongLink2(Mecapoint const&, Mecapoint const&, Vector const&, real ab2, real len, real weight);

    /// Link of stiffness `weight` and resting length `len`
    void addLongLink(Mecapoint const&, Mecapoint const&, real len, real weight);
    
    /// Link of stiffness `weight` and resting length `len`
    void addLongLink(Mecapoint const&, Interpolation const&, real len, real weight);
    
    /// Link of stiffness `weight` and resting length `len`
    void addLongLink(Interpolation const&, Interpolation const&, real len, real weight);
    
    /// Link of stiffness `weight` and resting length `len`
    void addLongLink4(Interpolation const&, const index_t pts, real, real, real, real, real len, real weight);
    
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

    /// Link of stiffness `weight` and resting length `arm1+arm2`, on the sides of both fibers
    void addSideSideLink2D(Interpolation const&, real arm1, Interpolation const&, real arm2, real weight);

    /// Link of stiffness `weight` and resting length `arm1+arm2`, on the sides of both fibers
    void addSideSideLink(Interpolation const&, Torque const& arm1, Interpolation const&, Torque const& arm2, real weight);

    /// Link of stiffness `weight` and resting length `arm`, on the sides of both fibers
    void addSideSideLink(Interpolation const&, Interpolation const&, real arm, real weight);

    /// Link of stiffness `weight` and perpendicular to first segment
    void addSlidingLink(Interpolation const&, Mecapoint const&, real weight);
    
    /// Link of stiffness `weight` and perpendicular to first segment
    void addSlidingLink(Interpolation const&, Interpolation const&, real weight);

    
#if ( DIM == 2 )
    /// Link on the side of first argument, using rotation `leg`, with the force along `dir` removed
    void addSideSlidingLink2D(Interpolation const&, real leg, Mecapoint const&, Vector const& dir, real weight);
#endif

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLinkS(Interpolation const&, Mecapoint const&, Torque const& arm, real weight);

    /// Link on the side of first argument, using rotation `leg`, with the force along `dir` removed
    void addSideSlidingLink3D(Interpolation const&, Torque const& leg, Mecapoint const&, Vector const& dir, real weight);

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink(FiberSegment const&, real, Mecapoint const&, real len, real weight);
    
    
#if ( DIM == 2 )
    /// Link on the side of first argument, using rotation `leg`, with the force along `dir` removed
    void addSideSlidingLink2D(Interpolation const&, real leg, Interpolation const&, Vector const& dir, real weight);
#endif

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLinkS(Interpolation const&, Interpolation const&, Torque const& arm, real weight);
    
    /// Link on the side of first argument, using rotation `leg`, with the force along `dir` removed
    void addSideSlidingLink3D(Interpolation const&, Torque const& leg, Interpolation const&, Vector const& dir, real weight);

    /// Link of stiffness `weight`, at distance `arm` on the side of first segment and perpendicular to this segment
    void addSideSlidingLink(FiberSegment const&, real, Interpolation const&, real len, real weight);
    
    
    /// Create a 3-way link with given weights on each branch
    void addTriLink(Interpolation const& pt1, real w1, Interpolation const& pt2, real w2, Interpolation const& pt3, real w3);

    //-------------------------- COMPUTING METHODS -----------------------------

    /// select which engine will be used, and ready it
    void selectStericEngine(Simul const&, SimulProp const&);

    /// add steric interactions between spheres, solids and fibers to Meca
    static void addStericInteractions(PointGrid&, Simul const&);

    /// add steric interactions between spheres, solids and fibers to Meca
    static void addStericInteractions(LocusGrid&, Simul const&);

    //-------------------------- COMPUTING METHODS -----------------------------
    
    /// import points and attributes indices
    void readyMecables();
    
    /// import useful parameters
    void importParameters(SimulProp const&);

    /// Allocate the memory necessary to `solve`. This must be called after the last addMecable
    void getReady(Simul const&);
    
    /// Calculate motion of all Mecables in the system; returns number of step of the iterative solver
    unsigned solve();
    
    /// transfer newly calculated point coordinates back to Mecables
    void apply();

    /// calculate Forces on Mecables and Lagrange multipliers for Fiber, without thermal motion
    void calculateForces();

    //----------------------- EXPORT/DEBUG FUNCTIONS ---------------------------

    /// member function pointer
    using MultiplyFuncPtr = void (Meca::*)(const real*, real*) const;

    /// linear operator, corresponding to the elasticity part of the dynamic matrix
    void multiplyElasticity(const real* X, real* Y) const;

    /// Count number of non-zero entries in the full system matrix
    size_t countTerms(real threshold) const;

    /// Extract matrix defined by template function in a preallocated C-array
    template < MultiplyFuncPtr >
    void getMatrix(real * mat, index_t lda) const;

    /// Save matrix defined by template function in binary format
    template < MultiplyFuncPtr >
    void dumpMatrix(FILE *, bool nat=true) const;
    
    /// Save mobility/projection matrix in binary format
    void dumpProjection(FILE *, bool nat=true) const;
    
    /// Save preconditionner in binary format
    void dumpPreconditionner(FILE *, bool nat=true) const;
    
    /// Save drag coefficients associated with each degree of freedom in binary format
    void dumpMobility(FILE *, bool nat=true) const;
    
    /// Save the object ID associated with each degree of freedom
    void dumpObjectID(FILE *) const;
    
    /// Output vectors and matrices, in a format that can be imported in MATLAB
    void dumpSystem(bool nat=true) const;
    
    
    /// Save the object ID associated with each degree of freedom
    void saveObjectID(FILE *) const;

    /// Save drag coefficients associated with each degree of freedom in binary format
    void saveMobility(FILE *) const;

    /// Save matrix defined by template function in Matrix Market format
    template < MultiplyFuncPtr >
    void saveMatrix(FILE *, size_t dim, real threshold) const;
    
    /// Output vectors and matrices, in a format that can be imported in MATLAB
    void saveSystem() const;

    /// Output vectors and matrices in various files (for debugging)
    void exportSystem() const;

    /// set Mecable:flag() according to connectivity defined by matrix elements
    void flagClusters() const;
    
    /// export bitmap images to reveal the matrices' sparsity patterns
    void saveMatrixBitmaps(const char[], unsigned inc = 1) const;
    
    /// save image of connectivity between Mecables
    void saveConnectivityBitmap() const;

};

#endif

