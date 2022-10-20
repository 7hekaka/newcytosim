// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University

/**
 * -----------------------------------------------------------------------------
 *                     -- Meca is the heart of Cytosim --
 * -----------------------------------------------------------------------------
 *             It solves the equations of motion for the Mecables,
 *      using implicit integration and iterative methods with sparse matrix
 * -----------------------------------------------------------------------------
 * @todo See if Lagrangian dynamics could work better than constrainted dynamics
 * -----------------------------------------------------------------------------
 */

#include <fstream>
#include "dim.h"

#include "assert_macro.h"
#include "blas.h"
#include "lapack.h"
#include "cytoblas.h"

#include "meca.h"
#include "mecable.h"
#include "messages.h"
#include "simul_prop.h"
#include "exceptions.h"
#include "vecprint.h"
#include "filepath.h"
#include "bicgstab.h"
#include "gmres.h"
#include "timer.h"
#include "simul.h"

/**
Set SEPARATE_RIGIDITY_TERMS to chose how Rigidity term are calculated:
   0. Rigidity terms are added to 'mISO' or 'mFUL'
   1. Rigidity values are calculated on the fly using `Mecable::addRigidity()`
.
With a sequential simulation, the second option is usually faster.
(in any case, with DIM==1, this should be 0)
 */
#define SEPARATE_RIGIDITY_TERMS ( DIM > 1 )

/// number of threads running in parallel
#define NUM_THREADS 1


#if NUM_THREADS > 1
/*
 Parallelization uses Intel's OpenMP.
 This requires a specific flag for the compiler, so adjust the makefile.inc
 CXXFLG := -std=gnu++14 -fopenmp
 */
#include <omp.h>
#endif

/// this define will enable explicit integration (should be off)
#define EXPLICIT_INTEGRATION 0

// shortcut
#if ( DIM == 1 )
#   define VECMULADDISO vecMulAdd
#elif ( DIM == 2 )
#   define VECMULADDISO vecMulAddIso2D
#elif ( DIM == 3 )
#   define VECMULADDISO vecMulAddIso3D
#endif

//------------------------------------------------------------------------------

#include "meca_inter.cc"
#include "meca_steric.cc"
#include "meca_rigidity.cc"
#include "meca_math.cc"
#include "meca_precond.cc"
#include "meca_util.cc"

//------------------------------------------------------------------------------
#pragma mark - Allocate

Meca::Meca()
: mecables(32, 32), pointGrid(*this), locusGrid(*this)
{
    tau_ = 0;
    alpha_ = 0;
    tolerance_ = 0;
    nPoints_ = 0;
    allocated_ = 0;
    bump_ = 0;
    ready_ = -1;
    steric_ = 0;
    precond_ = 0;
    verbose_ = 0;
#if NEW_CYTOPLASMIC_FLOW
    uniform_flow_dt_.reset();
#endif

    vPTS = nullptr;
    vSOL = nullptr;
    vBAS = nullptr;
    vRND = nullptr;
    vRHS = nullptr;
    vFOR = nullptr;
    vTMP = nullptr;
#if USE_ISO_MATRIX
    useFullMatrix = false;
#endif
    doNotify = 0;
    drawLinks = 0;
}


void allocate_vector(size_t s, real *& ptr, bool reset)
{
    free_real(ptr);
    ptr = new_real(s);
    if ( reset )
        zero_real(s, ptr);
}

void free_vector(real *& ptr)
{
    free_real(ptr);
    ptr = nullptr;
}

void Meca::allocate(size_t alc)
{
    if ( alc > allocated_ )
    {
        // make a multiple of chunk to keep pointers aligned:
        allocated_ = chunk_real(alc);
        
        alc = DIM * allocated_;
        allocate_vector(alc, vPTS, 1);
        allocate_vector(alc, vSOL, 1);
        allocate_vector(alc, vBAS, 0);
        allocate_vector(alc, vRND, 1);
        allocate_vector(alc, vRHS, 1);
        vFOR = vRHS; //allocate_vector(alc, vFOR, 1);
        vTMP = vSOL; //allocate_vector(alc, vTMP, 0);
        //std::clog << "Meca::allocate(" << allocated_ << ")\n";
    }
}


void Meca::release()
{
    //std::clog << "Meca::release()\n";
    free_vector(vPTS);
    free_vector(vSOL);
    free_vector(vBAS);
    free_vector(vRND);
    free_vector(vRHS);
    //free_vector(vFOR);
    //free_vector(vTMP);
}


size_t Meca::largestMecable() const
{
    size_t res = 0;
    for ( Mecable * mec : mecables )
        res = std::max(res, mec->nbPoints());
    return res;
}

size_t Meca::nbConstraints() const
{
    size_t res = 0;
    for ( Mecable * mec : mecables )
        res += mec->nbConstraints();
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Multiply

/**
 calculate the forces into `F`, given the Mecable coordinates `X`:
 
     F <- B + mISO * X + mFUL * X

 If `B == 0`, this term is ommited. With `B = vBAS` and `X = vPTS`, this
 function calculates the forces in the system in `F`:
 
     F <- vBAS + mISO * vPTS + mFUL * vPTS

 */
void Meca::calculateForces(const real* X, real const* B, real* F) const
{
    assert_true( empty() || ( X != F && X != B && F != B ));
    
#if USE_ISO_MATRIX
    if ( useFullMatrix )
        mFUL.vecMul(X, F);    // F <- mFUL * X
    else
        zero_real(dimension(), F);
    // F <- F + mISO * X
    mISO.VECMULADDISO(X, F);
#else
    // F <- mFUL * X
    mFUL.vecMul(X, F);
#endif
    
    // F <- F + B
    blas::xadd(dimension(), B, F);
}


void Meca::addAllRigidity(const real* X, real* Y) const
{
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
        mec->addRigidity(X+inx, Y+inx);
    }
}


#if PARALLELIZE_MATRIX

/**
calculate the matrix vector product corresponding to 'mec'

    Y <- X + alpha * speed( Y + P' * X );

*/
void Meca::multiply1(Mecable const* mec, const real* X, real* Y) const
{
    const size_t inx = DIM * mec->matIndex();

    // multiply the lines corresponding to this Mecable:
    mFUL.vecMul(X, Y, mec->matIndex(), mec->matIndex()+mec->nbPoints());

#if SEPARATE_RIGIDITY_TERMS
    mec->addRigidity(X+inx, Y+inx);
#endif

#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(X+inx, Y+inx);
#endif

    mec->projectForces(Y+inx, Y+inx);
    // Y <- X + beta * Y
    const real beta = -tau_ * mec->leftoverMobility();
    blas::xpay(DIM*mec->nbPoints(), X+inx, beta, Y+inx);
}


/**
 calculate the matrix product needed for the conjugate gradient algorithm
 
     Y <- X - time_step * speed( mISO + mFUL + P' ) * X;
 
 */
void Meca::multiply(const real* X, real* Y) const
{
#if NUM_THREADS > 1
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable * mec : mecables )
        multiply1(mec, X, Y);
#else
    mFUL.vecMul(X, Y);

    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
#  if SEPARATE_RIGIDITY_TERMS
        mec->addRigidity(X+inx, Y+inx);
#  endif
#  if ADD_PROJECTION_DIFF
        if ( mec->hasProjectionDiff() )
            mec->addProjectionDiff(X+inx, Y+inx);
#  endif
        mec->projectForces(Y+inx, Y+inx);
        // Y <- X + beta * Y
        const real beta = -tau_ * mec->leftoverMobility();
        blas::xpay(DIM*mec->nbPoints(), X+inx, beta, Y+inx);
    }
#endif
}

#else  // PARALLELIZE_MATRIX

/// Y <- X - time_step * speed( mISO + mFUL + P' ) * X;
void Meca::multiply(const real* X, real* Y) const
{
#if USE_ISO_MATRIX
    // Y <- mFUL * X
    if ( useFullMatrix )
        mFUL.vecMul(X, Y);
    else
        zero_real(dimension(), Y);
    // Y <- Y + mISO * X
    mISO.VECMULADDISO(X, Y);
#else
    mFUL.vecMul(X, Y);
#endif
    
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
#if SEPARATE_RIGIDITY_TERMS
        mec->addRigidity(X+inx, Y+inx);
#endif
#if ADD_PROJECTION_DIFF
        if ( mec->hasProjectionDiff() )
            mec->addProjectionDiff(X+inx, Y+inx);
#endif
        mec->projectForces(Y+inx, Y+inx);
        // Y <- X + beta * Y
        const real beta = -tau_ * mec->leftoverMobility();
        blas::xpay(DIM*mec->nbPoints(), X+inx, beta, Y+inx);
    }
}

#endif  // PARALLELIZE_MATRIX


//------------------------------------------------------------------------------
#pragma mark - Solve


/// qsort function comparing number of points of Mecables
static int compareMecables(const void * A, const void * B)
{
    size_t a = (*static_cast<Mecable *const*>(A))->nbPoints();
    size_t b = (*static_cast<Mecable *const*>(B))->nbPoints();
    return ( a < b ) - ( a > b );
}


/**
 This is currently limited to forces generated by Couple and Single
 attached to Mecables registered in 'meca'
 */
void Meca::setSomeInteractions()
{
    for ( Mecable const* mec : mecables )
    {
        mec->setInteractions(*this);
        
        Fiber const* fib = Fiber::toFiber(mec);
        for ( Hand * h = fib->firstHand(); h; h = h->next() )
        {
            HandMonitor const* m = h->monitor();
            Hand const* oh = m->otherHand(h);
            if ( oh > h  &&  oh->attached() )
                static_cast<Couple const*>(m)->setInteractions(*this);
            else if ( !oh )
                static_cast<Single const*>(m)->setInteractions(*this);
        }
    }
    
    if ( steric_ == 2 )
        addSomeStericInteractions();
}

/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::readyMecables()
{
    ready_ = 0;
    /*
     Attributes a position in the vector/matrix to each Mecable
     */
    size_t cnt = 0;
    for ( Mecable * mec : mecables )
    {
        mec->setIndex(cnt);
        cnt += mec->nbPoints();
    }
    nPoints_ = cnt;
    // allocate extra to allow some SIMD instruction burr
    allocate(cnt+1);
    
    // allocate sparse matrices:
#if USE_ISO_MATRIX
    mISO.resize(cnt);
    mISO.reset();
#endif
    mFUL.resize(DIM*cnt);
    mFUL.reset();
    
    // reset base:
    zero_real(DIM*cnt, vBAS);
    
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable * mec : mecables )
    {
        mec->putPoints(vPTS+DIM*mec->matIndex());
        mec->prepareMecable();
#if ( DIM > 1 ) && !SEPARATE_RIGIDITY_TERMS
        if ( mec->hasRigidity() )
        {
#   if USE_ISO_MATRIX
            addBendingRigidityMatrix(mISO, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   else
            addBendingRigidityBlockMatrix<DIM>(mFUL, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   endif
        }
#endif
    }
    //fprintf(stderr, "Meca::prepare() isnan %lu\n", has_nan(dimension(), vPTS));
}


void Meca::importParameters(SimulProp const& prop)
{
    tau_ = prop.time_step;
    alpha_ = prop.kT / tau_;
    tolerance_ = prop.tolerance;
    precond_ = prop.precondition;
    verbose_ = prop.verbose;
#if NEW_CYTOPLASMIC_FLOW
    uniform_flow_dt_ = prop.uniform_flow * prop.time_step;
#endif
}


void Meca::getReady(Simul const& sim)
{
    mecables.clear();
    for ( Fiber  * f= sim.fibers.first(); f; f=f->next() )
        addMecable(f);
    for ( Solid  * s= sim.solids.first(); s; s=s->next() )
        addMecable(s);
    for ( Sphere * o=sim.spheres.first(); o; o=o->next() )
        addMecable(o);
    for ( Bead   * b=  sim.beads.first(); b; b=b->next() )
        addMecable(b);
    
#if 0
    /*
     Sorting Mecables can improve multithreaded performance by distributing
     the work more equally between threads. Note that his operation is not free
     and for large systems random partitionning may not be so bad. Moreover for
     homogeneous systems (if all filaments have the same length) this is useless.
    */
    mecables.sort(compareMecables);
    
    /*
    for ( Mecable const* mec : mecables )
        std::clog << mec->reference() << " sorted " << mec->nbPoints() << "\n";
     */
#endif
    importParameters(sim.prop);
    selectStericEngine(sim);
    readyMecables();
}

/**
 Prepare matrices mISO and mFUL for multiplication
 This should be called after setInteractions()
 */
inline void Meca::prepareMatrices()
{
#if USE_ISO_MATRIX
    mISO.prepareForMultiply(DIM);
    useFullMatrix = mFUL.prepareForMultiply(1);
#else
    mFUL.prepareForMultiply(1);
#endif
}


/**
 Calculates forces due to external links, without adding Thermal motion,
 and also excluding bending elasticity of Fibers.
 
 Mecable::getForces will also sets the Lagrange multipliers for the Fiber.
 
 The function will not change the positions of any Mecable.
 */
void Meca::calculateForces()
{
    prepareMatrices();
    
    // vFOR <- external forces
    calculateForces(vPTS, vBAS, vFOR);

    for ( Mecable * mec : mecables )
    {
        mec->getForces(vFOR+DIM*mec->matIndex());
    }
}


/**
This updates the right-hand-side:
 
    rhs <- tau * Projection * ( rhs + alpha * rnd )
 
 Also prepare Projection diff is requested

 Vector 'rnd' is input, a set of independent Gaussian random numbers
 Vector 'rhs' is both input and output.
*/
real brownian1(Mecable* mec, real const* rnd, const real alpha, real tau, real* rhs)
{        
    real n = mec->addBrownianForces(rnd, alpha, rhs);

#if ADD_PROJECTION_DIFF == 2
    /* This uses the force to calculate the Lagrange multipliers */
    mec->makeProjectionDiff(rhs);
#endif

    // Calculate the right-hand-side of the system:
    mec->projectForces(rhs, rhs);
    
#if ADD_PROJECTION_DIFF
    /* assumes that the Lagrange multipliers were set correctly in the
     previous call to projectForces(); */
    mec->makeProjectionDiff(nullptr);
#endif

    // rhs <- tau * rhs, resulting in time_step * P * fff:
    blas::xscal(DIM*mec->nbPoints(), tau*mec->leftoverMobility(), rhs, 1);

    /*
     At this stage, `fff` contains the external forces in each vertex but also
     internal force such as bending elasticity terms, and the Lagrange multipliers
     do not represent the true tension in the filaments.
     Hence we do not call 'computeTensions(fff)' here
     */
    
    return n;
}


//------------------------------------------------------------------------------
#pragma mark - Solve & Apply

/**
 Meca::solve() solves the equation of motion with all Mecables:
 
     drag * ( Xnew - Xold ) / time_step = P * Force + Noise
 
 Where X is large a vector containing all the coordinates.
 P is the projection associated with constrains in the dynamics: P*P = P
 The projection P and scaling by `mobility = 1/drag` are implemented together via
 
     Mecable::projectForces()
     Mecable::leftoverMobility()
 
 We note here `mobP` the combination: mobP * X = ( 1/drag ) * P * X.
 To calculate Xnew, explicit integration would be:
 
     Xnew = Xold + time_step * mobP * ( Force + Noise )
 
 For a semi-implicit integration, we use a linearization of the force:
 
     Force(X) = M * X + B
 
 where M is a matrix and B a vector. The linearization is performed by the
 functions that update the matrix M, such as Meca::addLink() in `meca_inter.cc`.
 The force is usually linearized around the positions of equilibrium of that force,
 but it is then used around Xold, so we write:
 
     Force(X) = M * ( X - Xold ) + F
 
 where F = M * Xold + B = Force(Xold), leading to:
 
     ( I - time_step * mobP * M ) ( Xnew - Xold ) = time_step * mobP * ( F + Noise )
 
 with:
 
     Noise = std::sqrt(2*kT*time_step*mobility) * Gaussian(0,1)
 
 With implicit integration a large time step can be used.
 The matrix ( I - time_step * mobP * M ) remains definite positive.
 Moreover, both mobP and M are sparse, such that the matrix-vector product
 is calculated as follows in Meca::multiply():
 
 ( I - time_step * mobP * M ) * X = X - time_step * ( mobP * ( M * X ) )
 
 Further, M is not formed and instead we keep separate components:
 
     M = mISO + mFUL + Rigidity
 
 Where mISO is isotropic: it applies similarly in the X, Y and Z subspaces, while
 mFUL can accept crossterms between different subspaces. Using mISO is optional.
 In this way when calculating M * X, components can be processed in parallel.
 
 Normally, Meca::solve is called after:

     'mISO', 'mFUL' and 'B=vBAS' are set in Meca::setAllInteractions()
     'vPTS = Xold' is set from Mecables' points in Meca::prepare()
 
 The outline of the calculation is:
 
     'vRND' <- calibrated Gaussian random terms ~N(0,1)
     'vFOR' <- F = M * Xold + B
     'vRHS' <- set right-hand-side: time_step * mobP * F + vRND
     Solve the linear system ( I - time_step * mob * P * M ) * vSOL = vRHS
     'vSOL' <- solution to the linear system of equations
     'vPTS' <- calculate new positions: 'Xnew = vPTS + vSOL'
     'vFOR' <- calculate force with new positions: 'M * Xnew + B'
 
 The function Meca::apply() sends 'VPTS' and 'vFOR' back to the Mecable.
 
 
 Note: We currently solve ( 1 - time_step * P * M ) * X = Y
 Since both M and P are symmetric, following Woodbury's identity we have:
         X = Y + time_step * P * inverse( 1 - time_step * M * P ) * M * Y
 This adds 2 MAT.vec, but swaps M and P for the iterative solver.
 */
size_t Meca::solve()
{
    assert_true(ready_==0);

    prepareMatrices();
    
    // calculate external forces in vRHS:
    calculateForces(vPTS, vBAS, vRHS);
    
#if SEPARATE_RIGIDITY_TERMS
    addAllRigidity(vPTS, vRHS);
#endif
    
    /* 
     Fill `vRND` with Gaussian random numbers 
     This operation can be done in parallel, in a separate thread
     */
    RNG.gauss_set(vRND, dimension());
    
    /*
     Add Brownian motions to 'vRHS', and calculate vRHS by multiplying by mobilities.
     As Brownian terms are added, we record the magnitude of the typical smallest
     scalar contribution in `noiseLevel`. The dynamics will later be solved with 
     a residual that is proportional to this level:
     SimulProp::tolerance * noiseLevel
     As long as SimulProp::tolerance is smaller than 1, this should allow for a
     level of numerical error is small with respect to the Brownian noise in
     the system, and the results should be physically appropriate.
     */
    
    real noiseLevel = INFINITY;

    /*
     Add Brownian contributions and calculate Minimum value of it
      vRHS <- vRHS + mobility_coefficient * vRND
      vRHS <- tau * P * vRHS:
     */
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        real local = INFINITY;
        #pragma omp for
        for ( Mecable * mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            real n = brownian1(mec, vRND+inx, alpha_, tau_, vRHS+inx);
            local = std::min(local, n);
            //printf("thread %i min: %f\n", omp_get_thread_num(), local);
        }
        #pragma omp critical
        noiseLevel = std::min(noiseLevel, local);
    }

    // scale minimum noise level to serve as a measure of required precision
    noiseLevel *= tau_;
    
    //printf("noiseLeveld = %8.2e   variance(vRHS) / estimate = %8.4f\n",
    //       noiseLevel, blas::nrm2(dimension(), vRHS) / (noiseLevel * std::sqrt(dimension())) );

#if NEW_CYTOPLASMIC_FLOW
    /**
     Includes a constant fluid flow displacing all the objects along
     */
    if ( uniform_flow_dt_.norm() > REAL_EPSILON )
    {
        LOG_ONCE("NEW_CYTOPLASMIC_FLOW code enabled\n");
        for ( size_t p = 0; p < nbVertices(); ++p )
            uniform_flow_dt_.add_to(vRHS+DIM*p);
    }
#endif
    
#if EXPLICIT_INTEGRATION
    /*
     This implements the forward Euler integration, for testing purposes.
     The result is very inefficient, since we have built the stiffness matrix,
     which is not necessary for this explicit scheme.
     */
    blas::xadd(dimension(), vRHS, vPTS);
    ready_ = 1;
    return 1;
#endif

    // compute preconditionner:
    auto start = timer();
    computePreconditionner();
    auto factor = timer() - start;
    cycles_ = 0;

    /*
     Choose the initial guess for the solution of the system (Xnew - Xold):
     we could use the solution at the previous step, or a vector of zeros.
     Using the previous solution could be advantageous if the speed were
     somehow continuous. However, the system is without inertia. In addition,
     objects are considered in a random order to build the linear system, such
     that the blocks from two consecutive iterations do not match.
     From this, using zero for the initial guess seems safer:
     */
    zero_real(dimension(), vSOL);

    /*
     We now solve the system MAT * vSOL = vRHS  by an iterative method:
     the convergence tolerance is scaled to the contribution of Brownian motions
     contained in vRHS. Since we collected in 'noiseLevel' the minimul level
     of the Brownian contribution, this should work well if tolerance << 1
     */
    
    // tolerance is normally relative to the level of noise
    if ( noiseLevel > 0 )
        tolerance_ *= noiseLevel;
    else
    {
        // tolerance will be understood as an absolute quantity
        if ( alpha_ > 0 )
            Cytosim::log << "Warning: all Brownian terms are zero?\n";
    }
    
    /*
     With exact arithmetic, biConjugate Gradient should converge at most
     in a number of iterations equal to the size of the linear system,
     with each BCGGS iteration involving 2 matrix-vector multiplications.
     This limit is however too large, and we set an arbitrary limit in practice.
     */
    size_t max_iter = std::min(1111UL, 2*dimension());
    LinearSolvers::Monitor monitor(max_iter, tolerance_);

    //fprintf(stderr, "System size %6lu  limit %6lu  tolerance %f precondition %i\n", dimension(), max_iter, tolerance_, precond);

    /*
     GMRES may converge faster than BCGGS, but has overheads and uses more memory
     hence for very large systems, BCGGS is often advantageous.
     */

    //------- call the iterative solver:
    if ( precond_ )
    {
        // change initial condition to be `P * RHS`:
        precondition(vRHS, vSOL);
        LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
        //fprintf(stderr, "    BCGS     count %4u  residual %.3e\n", monitor.count(), monitor.residual());
        //LinearSolvers::GMRES(*this, vRHS, vSOL, 32, monitor, allocator_, mH, mV, temporary_);
    }
    else
    {
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
        //LinearSolvers::GMRES(*this, vRHS, vSOL, 64, monitor, allocator_, mH, mV, temporary_);
    }
    
    //fprintf(stderr, "    BCGS%u    count %4i  residual %.3e\n", precond, monitor.count(), monitor.residual());

#if ( 0 )
    // enable this to compare with GMRES using different restart parameters
    for ( int RS : {8, 16, 32} )
    {
        monitor.reset();
        zero_real(dimension(), vSOL);
        LinearSolvers::GMRES(*this, vRHS, vSOL, RS, monitor, allocator_, mH, mV, temporary_);
        fprintf(stderr, "    GMRES-%i  count %4i  residual %.3e\n", RS, monitor.count(), monitor.residual());
    }
#endif
#if ( 0 )
    // enable this to compare BCGS and GMRES
    fprintf(stderr, "    BCGS     count %4i  residual %.3e\n", monitor.count(), monitor.residual());
    monitor.reset();
    zero_real(dimension(), vSOL);
    LinearSolvers::GMRES(*this, vRHS, vSOL, 64, monitor, allocator_, mH, mV, temporary_);
    fprintf(stderr, "    GMRES-64 count %4i  residual %.3e\n", monitor.count(), monitor.residual());
#endif
#if ( 0 )
    // enable this to compare with another implementation of biconjugate gradient stabilized
    monitor.reset();
    zero_real(dimension(), vSOL);
    LinearSolvers::bicgstab(*this, vRHS, vSOL, monitor, allocator_);
    fprintf(stderr, "    bcgs     count %4i  residual %.3e\n", monitor.count(), monitor.residual());
#endif
    
    if ( !monitor.converged() )
    {
        Cytosim::out("Failed with size %lu precond %i flag %u count %4u residual %.3e (%.3f)",
            dimension(), precond_, monitor.flag(), monitor.count(), monitor.residual(), monitor.residual()/tolerance_);
        
        // in case the solver did not converge, we try other methods:
        monitor.reset();
#if !SAFER_CONVERGENCE
        // try with a different seed
        precondition(vRHS, vSOL);
#endif
        if ( precond_ )
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
        else
            LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
        
        Cytosim::out(" --> restarted: count %4i residual %.3e\n", monitor.count(), monitor.residual());

        // relax the convergence criteria a bit
        if ( monitor.residual() > 1.4142 * tolerance_ )
        {
            // try with our strongest preconditioner
            precond_ = 6;
            computePreconditionner();
            monitor.reset();
#if !SAFER_CONVERGENCE
            zero_real(dimension(), vSOL);
#endif
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
            Cytosim::out(" --> restarted precond 6: count %4i residual %.3e\n", monitor.count(), monitor.residual());
        }

#if SAFER_CONVERGENCE
        // relax the convergence criteria a bit more
        if ( monitor.residual() > 1.4142 * tolerance_ )
        {
            // try with different seed and strongest preconditioner
            monitor.reset();
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
            Cytosim::out(" -> final: count %4i residual %.3e\n", monitor.count(), monitor.residual());
        }
#endif

        // if the solver did not converge, its result cannot be used!
        if ( monitor.residual() > 1.4142 * tolerance_ )
            throw Exception("no convergence, residual ", monitor.residual(),
                            " achieved ", monitor.residual()/tolerance_);
    }

    //printf("\n   /sol "); VecPrint::print(std::cerr, dimension(), vSOL, 3);
    //printf("\n   >pts "); VecPrint::print(std::cerr, dimension(), vPTS, 3);

    auto solve = cycles_;
    cycles_ = timer() - start;
    
    //add the solution (the displacement) to update the Mecable's vertices
    blas::xadd(dimension(), vSOL, vPTS);

    ready_ = 1;

    // report on the matrix type and size, sparsity, and the number of iterations
    if (( 0 < doNotify ) || ( verbose_ & 1 ))
    {
        --doNotify;
        std::stringstream oss;
        oss << "\tsize " << DIM << "*" << nbVertices() << " kern " << largestMecable();
        //oss << " constraints " << nbConstraints();
#if USE_ISO_MATRIX
        oss << " " << mISO.what();
        if ( useFullMatrix )
#endif
        oss << " " << mFUL.what();
        oss << " precond " << precond_ << " (" << preconditionnerSize() << ")";
        oss << " count " << std::setw(4) << monitor.count();
        oss << " residual " << std::setw(11) << std::left << monitor.residual();
        size_t dim = dimension();
        if ( verbose_ & 8 )
        {
            // calculate true residual: tmp = rhs - A * x
            real * tmp = allocator_.bind(0);
            multiply(vSOL, tmp);
            blas::xsub(dim, vRHS, tmp);
            oss << ": " << std::setw(11) << std::left << blas::nrm8(dim, tmp);
            oss << " dx " << std::setw(11) << std::left << blas::nrm8(dim, vSOL);
        }
        if ( verbose_ & 4 )
        {
            unsigned cnt = std::max(1U, monitor.count());
            oss << "  cycles " << precond_ << "T " << std::setw(8) << cycles_;
            oss << " F " << std::setw(8) << factor << std::setw(6) << factor/cnt;
            oss << " S " << std::setw(8) << solve << std::setw(6) << solve/cnt;
            oss << " R " << std::setw(6) << ( cycles_ - factor - solve ) / cnt;
        }
        Cytosim::out << oss.str() << std::endl;
    }
    
    return monitor.count();
}


/**
 This transfers coordinates calculated in Meca::solve() back to the Mecables
 It also calculates the corresponding Forces and transfer them back.
 */
void Meca::apply()
{
    if ( ready_ )
    {
        /*
         Re-calculate forces with the new coordinates, excluding bending elasticity.
         In this way the forces returned to the fibers do not sum-up to zero, and
         are appropriate for example to calculate the effect of force on assembly.
         */
        calculateForces(vPTS, vBAS, vFOR);
        
        #pragma omp parallel for num_threads(NUM_THREADS)
        for ( Mecable * mec : mecables )
        {
            const size_t off = DIM * mec->matIndex();
            mec->addBrownianForces(vRND+off, alpha_, vFOR+off);
            //fprintf(stderr, "\n  "); VecPrint::print(stderr, DIM*mec->nbPoints(), vFOR+off, 2, DIM);
            if ( 1 )
            {
                // check validity of results:
                const size_t dim = DIM * mec->nbPoints();
                size_t a = has_nan(dim, vPTS+off);
                size_t b = has_nan(dim, vFOR+off);
                //fprintf(stderr, "Meca::solve isnan %i %i\n", a, b);
                if ( a | b )
                {
                    size_t c = has_nan(dim, vRND+off);
                    fprintf(stderr, "invalid results for Mecable %s isnan %lu %lu %lu\n", mec->reference().c_str(), a, b, c);
                    continue;
                }
            }
            // transfer new coordinates to Mecable:
            mec->getForces(vFOR+off);
            mec->getPoints(vPTS+off);
        }
    }
    else
    {
        // if !ready_, the result is not usable
        //printf("superfluous call to Meca::apply()\n");
    }
}

