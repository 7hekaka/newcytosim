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
: mecables(32), pointGrid(*this), locusGrid(*this)
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

void Meca::allocate(index_t alc)
{
    if ( alc > allocated_ )
    {
        // make a multiple of chunk to keep pointers aligned:
        allocated_ = chunk_real(alc);
        
        alc = DIM * allocated_;
        allocate_vector(alc, vPTS, 1);
        allocate_vector(alc, vSOL, 1);
        allocate_vector(alc, vBAS, 0);
        vFOR = vSOL; //allocate_vector(alc, vFOR, 1);
        allocate_vector(alc, vRHS, 1);
        //allocate_vector(alc, vTMP, 0);
        //std::clog << "Meca::allocate(" << allocated_ << ")\n";
    }
}


void Meca::release()
{
    //std::clog << "Meca::release()\n";
    free_vector(vPTS);
    free_vector(vSOL);
    free_vector(vBAS);
    //free_vector(vFOR);
    free_vector(vRHS);
    //free_vector(vTMP);
}


index_t Meca::largestMecable() const
{
    index_t res = 0;
    for ( Mecable * mec : mecables )
        res = std::max(res, mec->nbPoints());
    return res;
}

index_t Meca::smallestMecable() const
{
    index_t res = ~0U;
    for ( Mecable * mec : mecables )
        res = std::min(res, mec->nbPoints());
    return res;
}

index_t Meca::nbConstraints() const
{
    index_t res = 0;
    for ( Mecable * mec : mecables )
        res += mec->nbConstraints();
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Multiply

/**
 calculate the forces into `F`, given the Mecable coordinates `X`:
 
     F <- B + mISO * X + mFUL * X

 If `B == 0`, this term is omitted. With `B = vBAS` and `X = vPTS`, this
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
    #pragma omp parallel for
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
        mec->addRigidity(X+inx, Y+inx);
    }
}


/// Y <- X - time_step * speed( mISO + mFUL + diffP ) * X;
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
        const size_t cnt = DIM * mec->nbPoints();
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
        blas::xpay(cnt, X+inx, beta, Y+inx);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/// qsort function comparing number of vertices in Mecables
[[maybe_unused]]
static int compareMecables(const void * A, const void * B)
{
    auto a = (*static_cast<Mecable *const*>(A))->nbPoints();
    auto b = (*static_cast<Mecable *const*>(B))->nbPoints();
    return ( a < b ) - ( a > b );
}


/**
 Allocate and reset matrices and vectors necessary for `Meca::solve`,
 copy coordinates of Mecables into vPTS[]
 */
void Meca::readyMecables()
{
    ready_ = 0;
    /*
     Attributes a position in the vector/matrix to each Mecable
     */
    index_t cnt = 0;
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
    
    #pragma omp parallel for
    for ( Mecable * mec : mecables )
    {
        mec->putPoints(vPTS+DIM*mec->matIndex());
        mec->prepareMecable();
#if ( DIM > 1 ) && !SEPARATE_RIGIDITY_TERMS
        real jR = mec->jointRigidity();
        if ( jR != 0 )
        {
#   if USE_ISO_MATRIX
            addBendingRigidityMatrix(mISO, mec->matIndex(), mec->nbPoints(), jR);
#   else
            addBendingRigidityBlockMatrix<DIM>(mFUL, mec->matIndex(), mec->nbPoints(), jR);
#   endif
        }
#endif
    }
    //fprintf(stderr, "Meca::prepare() %lu NaNs\n", has_nan(dimension(), vPTS));
}


void Meca::importParameters(SimulProp const& prop)
{
    tau_ = prop.time_step;
    alpha_ = 2 * prop.kT / tau_;
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
    mecables.quick_sort(compareMecables);
    
    /*
    for ( Mecable const* mec : mecables )
        std::clog << mec->reference() << " sorted " << mec->nbPoints() << "\n";
     */
#endif
    importParameters(sim.prop);
    selectStericEngine(sim, sim.prop);
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
 
 Also prepares ProjectionDiff if the feature is enabled

 Vector 'fce' is input
 Vector 'rhs' is both input and output:
 - on entry, it should contain a set of independent Gaussian random numbers
 - on exit, it will return the result
*/
real brownian1(Mecable* mec, real const* fce, const real alpha, real tau, real* rhs)
{
    real n = mec->addBrownianForces(fce, alpha, rhs);

#if ADD_PROJECTION_DIFF == 7
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



/*
 Change the vector of noise, and recalculate the RHS of the system.
 This is a fail-safe option in case of failed convergence... outside the normal path
 
 recalculating vFOR, because it is equal to vSOL and therefore erased in `solve`.
 */
void Meca::renewBrownianForces()
{
    if ( vSOL == vFOR )
    {
        // calculate elastic forces in vFOR:
        calculateForces(vPTS, vBAS, vFOR);
        
#if SEPARATE_RIGIDITY_TERMS
        addAllRigidity(vPTS, vFOR);
#endif
    }
    
    //fprintf(stderr, "\n/"); VecPrint::print(stderr, dimension(), vRHS, 2, DIM);
    RNG.gauss_set(vRHS, dimension());

    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
        brownian1(mec, vFOR+inx, alpha_, tau_, vRHS+inx);
    }
    //fprintf(stderr, "\nL"); VecPrint::print(stderr, dimension(), vRHS, 2, DIM);
    //fprintf(stderr, "\n");
}

//------------------------------------------------------------------------------
#pragma mark - Solve & Apply

/// calculate the norm of `RHS - MATRIX * SOL`
real Meca::residualNorm() const
{
    size_t dim = dimension();
    real * tmp = allocator_.bind(0);
    multiply(vSOL, tmp);
    blas::xsub(dim, vRHS, tmp);
    real res = blas::nrm8(dim, tmp);
#if 0
    blas::xscal(dim, 1.0/res, tmp, 1);
    fprintf(stderr, "\n: "); VecPrint::print(stderr, dim, tmp, 2);
    fprintf(stderr, " R %f\n", res);
#endif
    return res;
}


/**
 `Meca::solve` solves the equation of motion with all Mecables:
 
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
 
     'vFOR' <- F = M * Xold + B (elastic force components)
     'vRHS' <- calibrated Gaussian random terms ~N(0,1)
     'vRHS' <- time_step * mobP * F + scale * vRHS (elastic forces+noise)
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
unsigned Meca::solve()
{
    const index_t dim = dimension();
    assert_true(ready_==0);

    prepareMatrices();
    
    // calculate elastic forces in vFOR:
    calculateForces(vPTS, vBAS, vFOR);
    
#if SEPARATE_RIGIDITY_TERMS
    addAllRigidity(vPTS, vFOR);
#endif
    
    /*
     Fill `vRHS` with Gaussian random numbers
     This operation could be done in parallel, in a separate thread
     */
    RNG.gauss_set(vRHS, dim);
    
    /*
     Add Brownian motions to 'vRHS', and calculate vRHS by multiplying by mobilities.
     As Brownian terms are added, we record their magnitude in `noiseLevel`,
     which will represent the magnitude of displacement due to Brownian motion.
     The dynamics will later be solved with a residual that is proportional:
     SimulProp::tolerance * noiseLevel
     As long as SimulProp::tolerance is smaller than 1, this should allow for a
     level of numerical error that is small with respect to the Brownian noise in
     the system, and the results should be physically appropriate.
     */
    
    real noiseLevel = INFINITY;

    /*
     Add Brownian contributions and estimate the magnitude of the noise
     on entry, vRHS contains Gaussian random numbers.
              vRHS <- vFOR + brownian_scale * vRHS
     the `brownian_scale` is calculated from the Mecable's mobility
              vRHS <- tau * P * vRHS:
     */
    #pragma omp parallel
    {
        real local = INFINITY;
        #pragma omp for
        for ( Mecable * mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            real n = brownian1(mec, vFOR+inx, alpha_, tau_, vRHS+inx);
            local = std::min(local, n);
            //printf("thread %i min: %f\n", omp_get_thread_num(), local);
        }
        #pragma omp critical
        noiseLevel = std::min(noiseLevel, local);
    }

    // scale minimum noise level to serve as a measure of required precision
    noiseLevel *= tau_;
    
#if NEW_CYTOPLASMIC_FLOW
    /**
     Includes a constant fluid flow displacing all the objects along
     */
    if ( uniform_flow_dt_.norm() > REAL_EPSILON )
    {
        LOG_ONCE("NEW_CYTOPLASMIC_FLOW code enabled\n");
        for ( index_t p = 0; p < nbVertices(); ++p )
            uniform_flow_dt_.add_to(vRHS+DIM*p);
    }
#endif
    
#if EXPLICIT_INTEGRATION
    /*
     This implements the forward Euler integration, for testing purposes.
     The method is quite inefficient, since we have constructed the stiffness matrix,
     which is not necessary for this explicit scheme. Force would be sufficient.
     */
    copy_real(dim, vRHS, vSOL);
    ready_ = 1;
    return 1;
#endif

    // compute preconditionner:
    auto start = machine_time();
    computePreconditionner();
    auto factorize = machine_time() - start;
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
    zero_real(dim, vSOL);

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
        if ( alpha_ > 0 )
        {
            // this is likely an error, since kT > 0:
            Cytosim::log("Warning: all Brownian terms are zero?\n");
        }
        // if kT==0, the tolerance will be understood as an absolute quantity
    }
    
    /*
     With exact arithmetic, biConjugate Gradient should converge within a number
     of iterations equal to the size of the linear system, with each BCGGS
     iteration involving 2 matrix-vector multiplications.
     Here we set the limit to the theoretical maximum:
     */
    index_t max_iter = 2 * dim;
    LinearSolvers::Monitor monitor(max_iter, tolerance_);

    //fprintf(stderr, "\nSystem size %6lu  limit %6lu  tolerance %f\n", dim, max_iter, tolerance_);

    //------- call the iterative solver:
    if ( precond_ )
    {
        LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
        //fprintf(stderr, "    BCGS     count %4u  residual %.3e\n", monitor.count(), monitor.residual());
    }
    else
    {
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
    }
    
#if ( 0 )
    // enable this to compare with another implementation of biconjugate gradient stabilized
    monitor.reset();
    zero_real(dim, vSOL);
    LinearSolvers::bicgstab(*this, vRHS, vSOL, monitor, allocator_);
    fprintf(stderr, "    bcgs     count %4i  residual %.3e\n", monitor.count(), monitor.residual());
#endif
    
    real resid = residualNorm();
    //fprintf(stderr, "bCGS%u %4i resid %.3e (%.3e)\n", precond_, monitor.count(), resid, monitor.residual());
    
    if ( resid > tolerance_ )
    {
        doNotify = 1;
        Cytosim::out.print("Failed size %lu precond %i flag %u count %4u residual %.3e (%.3e)",
                           dim, precond_, monitor.flag(), monitor.count(), resid, monitor.residual());
        
        // in case the solver did not converge, continue with same method/vectors:
        monitor.reset();
        if ( precond_ )
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
        else
            LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
        
        resid = residualNorm();
        Cytosim::out.print(" --> prolong: count %4i residual %.3e\n", monitor.count(), resid);
        
        if ( resid > tolerance_ )
        {
            // recalculate solution:
            monitor.reset();
            renewBrownianForces();
            zero_real(dim, vSOL);
            if ( precond_ )
                LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
            else
                LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
            resid = residualNorm();
            Cytosim::out.print(" --> new noise count %4i residual %.3e", monitor.count(), resid);
        }
        
        if ( resid > tolerance_ )
        {
            // recalculate solution:
            monitor.reset();
            zero_real(dim, vSOL);
            if ( precond_ )
            {
                LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
                resid = residualNorm();
                Cytosim::out.print(" --> extended: count %4i residual %.3e\n", monitor.count(), resid);
            }
            else
            {
                // try with our strongest preconditioner
                precond_ = 6;
                computePreconditionner();
                LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
                resid = residualNorm();
                Cytosim::out.print(" --> restarted precond 6: count %4i residual %.3e\n", monitor.count(), resid);
            }
        }

        // stop if the solver did not converge:
        if ( resid > 2 * tolerance_ )
            throw Exception("no convergence, residual ", resid, " achieved fraction ", resid/tolerance_);
    }

    //printf("\n   /sol "); VecPrint::print(std::cerr, dim, vSOL, 3);
    //printf("\n   >pts "); VecPrint::print(std::cerr, dim, vPTS, 3);

    auto apply = cycles_;
    cycles_ = machine_time() - start;

#if 0
    // print displacement of Mecable that has moved the most
    real dis = 0;
    static real mean = 0;
    Mecable * mac = nullptr;
    for ( Mecable * mec : mecables )
    {
        const index_t off = DIM * mec->matIndex();
        const index_t nbc = DIM * mec->nbPoints();
        real d = blas::nrm2(nbc, vSOL+off);
        if ( d > dis )
        {
            dis = d;
            mac = mec;
        }
    }
    if ( mac )
    {
        std::stringstream oss;
        oss << "max disp. (" << mac->reference() << ") = " << std::setprecision(6) << dis;
        oss << "  residual " << residualNorm();
        mean = 0.875 * mean + 0.125 * dis;
        if ( dis > 2 * mean )
            oss << "  ****";
        std::cout << oss.str() << std::endl;
    }
#endif
    
    ready_ = 1;
    
    // report on the matrix type and size, sparsity, and the number of iterations
    if (( 0 < doNotify ) || ( verbose_ & 1 ))
    {
        --doNotify;
        std::stringstream oss;
        if ( bump_ > 0 )
            oss << "\t\tfailed to compute " << bump_ << " / " << mecables.size() << " preconditionner blocks\n";
        oss << "\tsize " << DIM << "*" << nbVertices() << " kern " << largestMecable();
        //oss << " constraints " << nbConstraints();
#if USE_ISO_MATRIX
        oss << " " << mISO.what();
        if ( useFullMatrix )
#endif
        oss << " " << mFUL.what();
        //oss << " noise " << noiseLevel << "  " << blas::nrm8(dim, vRND);
        if ( precond_ )
            oss << " precond " << precond_ << " (" << preconditionnerSize() << ")" << CHOUCROUTE;
        oss << " count " << std::setw(4) << monitor.count();
        oss << " residual " << std::setw(11) << std::left << monitor.residual();
        if ( verbose_ & 4 )
        {
            unsigned cnt = std::max(1U, monitor.count());
            oss << "  cycles " << precond_ << "T " << std::setw(8) << cycles_;
            oss << " F " << std::setw(8) << factorize << std::setw(6) << factorize/cnt;
            oss << " S " << std::setw(8) << apply << std::setw(6) << apply/cnt;
            oss << " M " << std::setw(6) << ( cycles_ - factorize - apply ) / cnt;
        }
        Cytosim::out << oss.str() << std::endl;
        //std::clog << oss.str() << std::endl;
    }
    
    return monitor.count();
}


/**
 This transfers coordinates calculated in `Meca::solve` back to the Mecables
 It also calculates the corresponding Forces and transfer them back.
 */
void Meca::apply()
{
    if ( ready_ )
    {
#if 0
        // print vertices of Mecables
        for ( Mecable * mec : mecables )
        {
            const size_t off = DIM * mec->matIndex();
            const size_t len = DIM * std::min(4U, mec->nbPoints());
            VecPrint::print(mec->reference().c_str(), len, vSOL+off, 5, DIM);
        }
        fprintf(stderr, "\n");
#endif
        // add calculated displacement to obtain vertices positions:
        blas::xadd(dimension(), vSOL, vPTS);

        /*
         Re-calculate forces with the new coordinates, excluding bending elasticity.
         In this way the forces returned to the fibers do not sum-up to zero, and
         are appropriate for example to calculate the effect of force on assembly.
         */
        calculateForces(vPTS, vBAS, vFOR);
        /*
         It is unclear if we should add Brownian values to the force returned to Mecables:
         The Brownian components will scale like sqrt( kT * drag_coefficient / time_step ),
         So it will be large for small time_step and large objects. They are not added.
         */

        #pragma omp parallel for
        for ( Mecable * mec : mecables )
        {
            const index_t off = DIM * mec->matIndex();
            const index_t len = DIM * mec->nbPoints();
#ifndef __FAST_MATH__
            // check validity of results:
            size_t a = has_nan(len, vPTS+off);
            size_t b = has_nan(len, vFOR+off);
            if ( a | b )
            {
                fprintf(stderr, "invalid results for %s : %lu %lu NaNs / %u\n", mec->reference().c_str(), a, b, len);
                continue;
            }
#endif
            // transfer new coordinates to Mecable:
            mec->getForces(vFOR+off);
            mec->getPoints(vPTS+off);
        }
    }
    else
    {
        // if `ready_ == false`, the results are not usable
        //printf("superfluous call to Meca::apply()\n");
    }
}

