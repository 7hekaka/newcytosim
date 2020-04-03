// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 * ------------------------------------------------------------------------------
 *                   -- Meca is the heart of Cytosim --
 * ------------------------------------------------------------------------------
 *             It solves the equations of motion for the Mecables,
 *      using implicit integration and iterative methods with sparse matrix
 * ------------------------------------------------------------------------------
 * @todo See if Lagrangian dynamics could work better than constrainted dynamics
 * @todo Implement the PARDISO sparse matrix format
 * @todo Check if IDR(s) would perform better than BCGS or GMRES
 * ------------------------------------------------------------------------------
 */

#include <fstream>

#include "meca.h"
#include "mecable.h"
#include "messages.h"
#include "simul_prop.h"
#include "blas.h"
#include "lapack.h"
#include "exceptions.h"
#include "vecprint.h"
#include "filepath.h"
#include "tictoc.h"
#include "bicgstab.h"
#include "gmres.h"

#include "meca_inter.cc"
#include "meca_rigidity.cc"
#include "meca_math.cc"

/**
 Add correction term to the constrainted dynamics
 The effect is to stabilize fibers under traction, at some modest CPU cost.
*/
#define ADD_PROJECTION_DIFF 1


/**
 Assumes that Matrix Vector-multiplication can be distributed
 */
#define PARALLELIZE_MATRIX 0


/**
Set SEPARATE_RIGIDITY_TERMS to chose how Rigidity term are calculated:
   0. Rigidity terms are added to Matrix 'mB' or 'mC'
   1. Rigidity values are calculated on the fly using `Mecable::addRigidity()`
.
With a sequential simulation, the second option is usually faster.
(in any case, with DIM==1, this should be 0)
 */
#define SEPARATE_RIGIDITY_TERMS ( DIM > 1 )

/// this define will enable explicit integration (should be off)
#define EXPLICIT_INTEGRATION 0


/// use this to generate code for validation
#define DEBUG_MECA 0


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


//------------------------------------------------------------------------------
#pragma mark - Allocate

Meca::Meca()
: mecables(32)
{
    ready_ = 0;
    nbPts = 0;
    allocated_ = 0;
    vPTS = nullptr;
    vSOL = nullptr;
    vBAS = nullptr;
    vRND = nullptr;
    vRHS = nullptr;
    vFOR = nullptr;
    vTMP = nullptr;
#if USE_ISO_MATRIX
    useMatrixC = false;
#endif
    drawLinks = false;
    time_step = 0;
}


void allocate_vector(size_t s, real *& ptr, bool reset)
{
    free_real(ptr);
    ptr = new_real(s);
    if ( reset )
        zero_real(s, ptr);
}


void Meca::allocate(size_t alc)
{
    //allocate the vectors
    if ( alc > allocated_ )
    {
        // make a multiple of chunk to align pointers:
        allocated_ = chunk_real(alc);
        
        // pad with 4 doubles to allow SIMD instruction burr
        alc = DIM * allocated_ + 4;
        allocate_vector(alc, vPTS, 1);
        allocate_vector(alc, vSOL, 1);
        allocate_vector(alc, vBAS, 0);
        allocate_vector(alc, vRND, 1);
        allocate_vector(alc, vRHS, 1);
        allocate_vector(alc, vFOR, 1);
        allocate_vector(alc, vTMP, 0);
    }
}


void Meca::release()
{
    //std::clog << "Meca::release()\n";
    free_real(vPTS);
    free_real(vSOL);
    free_real(vBAS);
    free_real(vRND);
    free_real(vRHS);
    free_real(vFOR);
    free_real(vTMP);
    vPTS = nullptr;
    vSOL = nullptr;
    vBAS = nullptr;
    vRND = nullptr;
    vRHS = nullptr;
    vFOR = nullptr;
    vTMP = nullptr;
}


size_t Meca::largestMecable() const
{
    size_t res = 0;
    for ( Mecable * mec : mecables )
        res = std::max(res, mec->nbPoints());
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Multiply


// shortcut

#if ( DIM == 1 )
#   define VECMULADDISO vecMulAdd
#elif ( DIM == 2 )
#   define VECMULADDISO vecMulAddIso2D
#elif ( DIM == 3 )
#   define VECMULADDISO vecMulAddIso3D
#endif

/**
 calculate the forces into `F`, given the Mecable coordinates `X`:
 
     F <- B + mB * X + mC * X

 If `B == 0`, this term is ommited. With `B = vBAS` and `X = vPTS`, this
 function calculates the forces in the system in `F`:
 
     F <- vBAS + mB * vPTS + mC * vPTS

 */
void Meca::calculateForces(const real* X, real const* B, real* F) const
{
    assert_true( empty() || ( X != F && X != B && F != B ));

    // F <- B  of F <- 0
    if ( B )
        copy_real(dimension(), B, F);      //blas::xcopy(dimension(), B, 1, F, 1);
    else
        zero_real(dimension(), F);
    
    // F <- F + mB * X
    
#if USE_ISO_MATRIX
    mB.VECMULADDISO(X, F);

    if ( useMatrixC )
#endif
    {
        // F <- F + mC * X
        mC.vecMulAdd(X, F);
    }
}


void Meca::addAllRigidity(const real* X, real* Y) const
{
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            const size_t inx = DIM * (*mci)->matIndex();
            (*mci)->addRigidity(X+inx, Y+inx);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
        mec->addRigidity(X+inx, Y+inx);
    }
#endif
}


#if PARALLELIZE_MATRIX

/**
calculate the matrix vector product corresponding to 'mec'

    Y <- X + alpha * speed( Y + P' * X );

*/
void Meca::multiply1(Mecable const* mec, const real* X, real* Y) const
{
    const size_t inx = DIM * mec->matIndex();
    const size_t bks = DIM * mec->nbPoints();

    // multiply the columns corresponding to this Mecable:
    mC.vecMul(X, Y, inx, inx+bks);

#if SEPARATE_RIGIDITY_TERMS
    mec->addRigidity(X+inx, Y+inx);
#endif

#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(X+inx, Y+inx);
#endif

    mec->projectForces(Y+inx, Y+inx);

    // Y <- X + alpha * Y
    blas::xpay(bks, X+inx, -time_step*mec->leftoverMobility(), Y+inx);
}


/**
 calculate the matrix product needed for the conjugate gradient algorithm
 
     Y <- X - time_step * speed( mB + mC + P' ) * X;
 
 */
void Meca::multiply(const real* X, real* Y) const
{
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            multiply1(*mci, X, Y);
            mci += NUM_THREADS;
        }
    }
#elif ( 0 )
    for ( Mecable * mec : mecables )
        multiply1(mec, X, Y);
#else
    mC.vecMul(X, Y);

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
        // Y <- X + alpha * Y
        blas::xpay(DIM*mec->nbPoints(), X+inx, -time_step*mec->leftoverMobility(), Y+inx);
    }
#endif
}

/**
 This is equivalent to
     multiply(X, T);
     precondition(T, Y);
 for the block corresponding to 'mec'
 */
void Meca::multiply_precondition1(Mecable const* mec, real const* X, real* T, real* Y) const
{
    const size_t inx = DIM * mec->matIndex();
    const size_t bks = DIM * mec->nbPoints();

    mC.vecMul(X, T, inx, inx+bks);
    
#if SEPARATE_RIGIDITY_TERMS
    mec->addRigidity(X+inx, T+inx);
#endif
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(X+inx, T+inx);
#endif
    
    if ( Y == X )
    {
        // T <- P * T
        mec->projectForces(T+inx, T+inx);
        // X <- X + alpha * T = X + alpha * P * FORCE
        blas::xaxpy(bks, -time_step*mec->leftoverMobility(), T+inx, 1, Y+inx, 1);
    }
    else
    {
        // Y <- P * T
        mec->projectForces(T+inx, Y+inx);
        // Y <- X + alpha * Y = X + alpha * P * FORCE
        blas::xpay(bks, X+inx, -time_step*mec->leftoverMobility(), Y+inx);
    }

    if ( mec->useBlock() )
    {
        // Y <- PRECONDITIONNER_BLOCK * Y
        int info = 0;
        lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y+inx, bks, &info);
        assert_true(info==0);
    }
}


/**
 This is equivalent to
 
     multiply(X, T);       // T <- M*X
     precondition(T, Y);   // Y <- P*T
 
 This is used for left-sided preconditionning
 */
void Meca::multiply_precondition(real const* X, real* T, real* Y) const
{
#if NUM_THREADS > 1
#pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            multiply_precondition1(*mci, X, T, Y);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : mecables )
        multiply_precondition1(mec, X, T, Y);
#endif
}


void Meca::multiply_precondition1(Mecable const* mec, real const* X, real* Y) const
{
    assert_true(X!=Y);
    
    const size_t inx = DIM * mec->matIndex();
    const size_t bks = DIM * mec->nbPoints();

    mC.vecMul(X, Y, inx, inx+bks);
    
#if SEPARATE_RIGIDITY_TERMS
    mec->addRigidity(X+inx, Y+inx);
#endif
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(X+inx, Y+inx);
#endif
    
    // Y <- P * Y
    mec->projectForces(Y+inx, Y+inx);
    // Y <- X + alpha * Y = X + alpha * P * FORCE
    blas::xpay(bks, X+inx, -time_step*mec->leftoverMobility(), Y+inx);
    
    if ( mec->useBlock() )
    {
        // Y <- PRECONDITIONNER_BLOCK * Y
        int info = 0;
        lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y+inx, bks, &info);
        assert_true(info==0);
    }
}

/**
 This calculates `Y <- P * M * X`, and is equivalent to
 
     multiply(X, TEMP);
     precondition(TEMP, Y);

 as needed for left-sided preconditionning.
 But it does so using only one pass, when multi-threaded
 */
void Meca::multiply_precondition(real const* X, real* Y) const
{
#if NUM_THREADS > 1
#pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            multiply_precondition1(*mci, X, Y);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : mecables )
        multiply_precondition1(mec, X, Y);
#endif
}

#else  // PARALLELIZE_MATRIX

/// Y <- X - time_step * speed( mB + mC + P' ) * X;
void Meca::multiply(const real* X, real* Y) const
{
#if USE_ISO_MATRIX
    // Y <- ( mB + mC ) * X
    calculateForces(X, nullptr, Y);
#else
    mC.vecMul(X, Y);
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
        // Y <- X + alpha * Y
        blas::xpay(DIM*mec->nbPoints(), X+inx, -time_step*mec->leftoverMobility(), Y+inx);
    }
}

#endif  // PARALLELIZE_MATRIX

//------------------------------------------------------------------------------
#pragma mark - Precondition


/**
 Get the total diagonal block corresponding to an Object, which is:
 
     I - time_step * P ( mB + mC + P' )
 
 The result is constructed by using functions from mB and mC
 This block is square but not symmetric!
 */
void Meca::getBlock(real* res, const Mecable * mec) const
{
    const size_t np = mec->nbPoints();
    const size_t bs = DIM * np;
    
    zero_real(bs*bs, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    if ( mec->hasRigidity() )
    {
        addRigidityLower(res, bs, mec->nbPoints(), mec->fiberRigidity());
        //std::clog<<"Rigidity block " << mec->reference() << "\n";
        //VecPrint::print(std::clog, bs, bs, res, bs, 0);
    }
#endif
#if USE_ISO_MATRIX
    mB.addTriangularBlock(res, bs, mec->matIndex(), np, DIM);
#endif
    expand_lower_matrix(bs, res);
    
#if USE_ISO_MATRIX
    if ( useMatrixC )
#endif
        mC.addDiagonalBlock(res, bs, DIM*mec->matIndex(), bs);
    
#if ( 0 )
    std::clog<<"mB+mC block:\n";
    VecPrint::print(std::clog, bs, bs, res, bs);
#endif
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
    {
        // Include the corrections P' in preconditioner, vector by vector.
        real* tmp = vTMP + DIM * mec->matIndex();
        zero_real(bs, tmp);
        for ( size_t ii = 0; ii < bs; ++ii )
        {
            tmp[ii] = 1.0;
            mec->addProjectionDiff(tmp, res+bs*ii);
            tmp[ii] = 0.0;
        }
        //std::clog<<"dynamic with P'\n";
        //VecPrint::print(std::clog, bs, bs, res, bs);
    }
#endif
    
    //compute the projection, by applying it to each column vector:
    /*
     This could be vectorized by having projectForces()
     accept multiple vectors as arguments, using SIMD instructions
     */
    for ( size_t i = 0; i < bs; ++i )
        mec->projectForces(res+bs*i, res+bs*i);
    
    // scale
    real beta = -time_step * mec->leftoverMobility();
    //blas::xscal(bs*bs, beta, res, 1);
    for ( size_t n = 0; n < bs*bs; ++n )
        res[n] = beta * res[n];
    // add Identity matrix:
    for ( size_t i = 0; i < bs; ++i )
        res[i+bs*i] += 1.0;
}


/**
 This version builds the diagonal block indirectly using Meca:multiply().
 This is a slow method that calls 'multiply()' n-times, where
 'n' is the size of the block.
 
 This should be used for validation only.
*/
void Meca::extractBlock(real* res, const Mecable* mec) const
{
    const size_t dim = dimension();
    const size_t bks = DIM * mec->nbPoints();
    const size_t off = DIM * mec->matIndex();
    
    assert_true( off+bks <= dim );
    real * vec = new_real(dim);
    real * tmp = new_real(dim);
    
    zero_real(dim, vec);
    //zero_real(bs*bs, res);
    
    // proceed column by column:
    for ( size_t jj = 0; jj < bks; ++jj )
    {
        vec[jj+off] = 1.0;
        multiply(vec, tmp);
        vec[jj+off] = 0.0;
        blas::xcopy(bks, tmp+off, 1, res+jj*bks, 1);
    }
    
    free_real(vec);
    free_real(tmp);
}


// DEBUG: compare `blk` with block extracted using extractBlockSlow()
void Meca::verifyBlock(const Mecable * mec, const real* blk)
{
    const size_t bks = DIM * mec->nbPoints();
    real * wrk = new_real(bks*bks);
    
    extractBlock(wrk, mec);
    
    blas::xaxpy(bks*bks, -1.0, blk, 1, wrk, 1);
    real err = blas::nrm2(bks*bks, wrk);
 
    std::clog << "verifyBlock ";
    std::clog << std::setw(10) << mec->reference() << " " << std::setw(6) << bks;
    std::clog << "  | B - K | = " << std::setprecision(3) << err << std::endl;
    
    if ( err > bks * bks * REAL_EPSILON )
    {
        VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        
        size_t s = std::min(16UL, bks);
        extractBlock(wrk, mec);
        std::clog << " blockS\n";
        VecPrint::print(std::clog, s, s, wrk, bks, 3);
        
        std::clog << " block \n";
        VecPrint::print(std::clog, s, s, blk, bks, 3);
    }
    
    free_real(wrk);
}


/**
 Multiply here `blk` with the dynamic block extracted by extractBlockSlow()
 and check that we recover the identity matrix
 */
void Meca::checkBlock(const Mecable * mec, const real* blk)
{
    const size_t bks = DIM * mec->nbPoints();
    
    std::clog << "  checkBlock " << mec->useBlock() << " ";
    std::clog << std::setw(10) << mec->reference() << " " << std::setw(6) << bks;
 
    if ( !mec->useBlock() )
    {
        std::clog << std::endl;
        return;
    }
    
    real * wrk = new_real(bks*bks);
    real * mat = new_real(bks*bks);
    real * vec = new_real(bks);
    
    extractBlock(wrk, mec);   // wrk <- PRECONDITIONNER_BLOCK
   
    copy_real(bks*bks, wrk, mat);  // mat <- wrk
    if ( mec->useBlock() )
    {
        int info = 0;
        // mat <- PRECONDITIONNER_BLOCK * mat
        lapack::xgetrs('N', bks, bks, blk, bks, mec->pivot(), mat, bks, &info);
        assert_true(info==0);
    }
    
    for ( size_t k=0; k < bks*bks; k += 1+bks )
        mat[k] -= 1.0;
    
    real err = blas::nrm2(bks*bks, mat) / bks;
    std::clog << " | 1 - PM | = " << std::setprecision(3) << err;
    
    if ( 1 )
    {
        // chose initial vector for power iteration
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        //this is not valid for triangular preconditionner
        real eig = -1;
        if ( mec->useBlock() )
            eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, mat);
        std::clog << "  eigen(1-PM) = " << eig;
    }
    
    std::clog << std::endl;

    if ( err > 1 )
    {
        // print preconditionner block for visual inspection:
        unsigned s = std::min(16UL, bks);
        std::clog << " matrix: \n";
        VecPrint::print(std::clog, s, s, wrk, bks);
        std::clog << "\nprecond: \n";
        VecPrint::print(std::clog, s, s, blk, bks);
        std::clog << "\nprecond * matrix:\n";
        VecPrint::print(std::clog, s, s, mat, bks);
        std::clog << "\n";
    }
    free_real(vec);
    free_real(mat);
    free_real(wrk);
}


/**
 Arrays 'tmp' and 'wrk' should be of size (nb_points*DIM)^2 or more
 'method' in {0, 2}
 
 Method 2 attempts to keep the existing block, and evaluates it fitness
 by iterations to identify the largest eigenvalue, which may or may not perform well
 */
void Meca::computePreconditionnerAlt(Mecable* mec, real* tmp, real* wrk, size_t wrksize)
{
    const size_t bks = DIM * mec->nbPoints();
    assert_true( bks*bks <= wrksize );

    bool may_keep = false;
    
    if ( mec->blockSize() == bks )
        may_keep = true;
    else
        mec->allocateBlock();
    
    real* blk = mec->block();
    real* vec = vTMP + DIM * mec->matIndex();

    if ( may_keep )
    {
        // extract diagonal matrix block corresponding to this Mecable:
        getBlock(wrk, mec);
        
        // chose initial vector for power iteration
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        
        // power iterations scale like (bs^2)
        real eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, tmp);
 
        //std::clog << "mec " << std::setw(4) << mec->identity();
        if ( eig < 1.0 )
        {
            //std::clog << "    keep     eigen(1-PM) = " << eig << std::endl;
            mec->useBlock(1);
            return;
        }
        else
        {
            //std::clog << "    RENEW    eigen(1-PM) = " << eig << std::endl;
            blas::xcopy(bks*bks, wrk, 1, blk, 1);
        }
    }
    else
    {
        // extract diagonal matrix block corresponding to this Mecable:
        getBlock(blk, mec);
    }
    
    //verifyBlock(mec, blk);

    int info = 0;
    
    /**
     We invert here the full block using LAPACK general matrix functions
     the workload scales as SIZE ^ 3, as a function of the size of the block
     */
    //TEST truncate_matrix(bs, blk, 3*DIM-1);
    
    // calculate LU factorization:
    lapack::xgetf2(bks, bks, blk, bks, mec->pivot(), &info);
    
    if ( info )
    {
        //failed to factorize matrix !!!
        std::clog << "Meca::computePreconditionner failed (lapack::xgetf2, info = " << info << ")\n";
        return;
    }
    
    mec->useBlock(1);
    //checkBlock(mec, blk);
}


// Compute sequentially all the blocks of the preconditionner
void Meca::computePreconditionnerAlt()
{
    const size_t vecsize = DIM * largestMecable();
    const size_t wrksize = vecsize * vecsize;
 
    // allocate memory:
    real * tmp = new_real(wrksize);
    real * wrk = new_real(wrksize);
    
    for ( Mecable * mec : mecables )
        computePreconditionnerAlt(mec, tmp, wrk, wrksize);
    
    free_real(wrk);
    free_real(tmp);
}


/**
Compute block of the preconditionner corresponding to 'mec'
 */
void Meca::computePreconditionner(Mecable* mec)
{
    const size_t bks = DIM * mec->nbPoints();

    mec->allocateBlock();
    real* blk = mec->block();
 
    // extract diagonal matrix block corresponding to this Mecable:
    getBlock(blk, mec);

    //verifyBlock(mec, blk);

    // calculate LU factorization:
    int info = 0;
    lapack::xgetf2(bks, bks, blk, bks, mec->pivot(), &info);
    
    if ( info == 0 )
    {
        mec->useBlock(1);
        //checkBlock(mec, blk);
        //std::clog << "Meca::computePreconditionner(" << mec->reference() << ")\n";
    }
    else
    {
        assert_true(mec->useBlock() == 0);
        std::clog << "Meca::computePreconditionner failed (lapack::xgetf2, info " << info << ")\n";
    }
}


/// Compute all the blocks of the preconditionner
/**
 This is method 1, that should perform well and can be multithreaded
 */
void Meca::computePreconditionner()
{
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            computePreconditionner(*mci);
            mci += NUM_THREADS;
        }
        //printf("thread %i complete %i\n", omp_get_thread_num(), TicToc::microseconds());
    }
#else
    for ( Mecable * mec : mecables )
        computePreconditionner(mec);
#endif
}


void Meca::precondition(const real* X, real* Y) const
{    
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            Mecable const* mec = *mci;
            const size_t inx = DIM * mec->matIndex();
            const size_t bks = DIM * mec->nbPoints();
            if ( Y != X )
                blas::xcopy(bks, X+inx, 1, Y+inx, 1);  // Y <- X
            if ( mec->useBlock() )
            {
                int info = 0;
                // Y <- PRECONDITIONNER_BLOCK * Y
                lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y+inx, bks, &info);
                assert_true(info==0);
            }
            mci += NUM_THREADS;
        }
    }
#else
    if ( Y != X )
        blas::xcopy(dimension(), X, 1, Y, 1);
    for ( Mecable const* mec : mecables )
    {
        if ( mec->useBlock() )
        {
            int info = 0;
            const size_t inx = DIM * mec->matIndex();
            const size_t bks = DIM * mec->nbPoints();
            // Y <- PRECONDITIONNER_BLOCK * Y
            lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y+inx, bks, &info);
        }
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/// function to sort Mecables
int smaller_mecable(const void * ap, const void * bp)
{
    Mecable const** a = (Mecable const**)(ap);
    Mecable const** b = (Mecable const**)(bp);

    if ( (*a)->nbPoints() > (*b)->nbPoints() ) return -1;
    if ( (*a)->nbPoints() < (*b)->nbPoints() ) return  1;
    return 0;
}


/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::prepare()
{
#if NUM_THREADS > 1
    /*
     Sorting Mecables can improve multithreaded performance by distributing
     the work more equally between threads. Note that his operation is not free
     and for large systems random partitionning may not be so bad. Moreover for
     homogeneous systems (if all filaments have the same length) this is useless.
    */
    mecables.sort(smaller_mecable);
    
    /*
    for ( Mecable const* mec : mecables )
        std::clog << mec->reference() << " " << mec->nbPoints() << "\n";
     */
#endif
    
    /*
     Attributes a position in the vector/matrix to each Mecable
     */
    size_t cnt = 0;
    for ( Mecable * mec : mecables )
    {
        mec->matIndex(cnt);
        cnt += mec->nbPoints();
    }
    nbPts = cnt;
    allocate(cnt);
    
#if USE_ISO_MATRIX
    //allocate the sparse matrices:
    mB.resize(cnt);
    mB.reset();
#endif
    
    mC.resize(DIM*cnt);
    mC.reset();
    
    // reset base:
    zero_real(DIM*cnt, vBAS);
    
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            Mecable * mec = *mci;
            mec->putPoints(vPTS+DIM*mec->matIndex());
            mec->prepareMecable();
#if ( DIM > 1 ) && !SEPARATE_RIGIDITY_TERMS
            if ( mec->hasRigidity() )
            {
#   if USE_ISO_MATRIX
                addRigidityMatrix(mB, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   else
                addRigidityBlockMatrix(mC, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   endif
            }
#endif
            mec->useBlock(0);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : mecables )
    {
        mec->putPoints(vPTS+DIM*mec->matIndex());
        mec->prepareMecable();
#if ( DIM > 1 ) && !SEPARATE_RIGIDITY_TERMS
        if ( mec->hasRigidity() )
        {
#   if USE_ISO_MATRIX
            addRigidityMatrix(mB, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   else
            addRigidityBlockMatrix(mC, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   endif
        }
#endif
        mec->useBlock(0);
    }
#endif
}


/**
 Prepare matrices mB and mC for multiplication
 This should be called after setInteractions()
 */
void Meca::prepareMatrices()
{
#if USE_ISO_MATRIX
    mB.prepareForMultiply(DIM);
    useMatrixC = mC.prepareForMultiply(1);
#else
    mC.prepareForMultiply(1);
#endif
}


/**
 Calculates forces due to external links, without adding Thermal motion,
 and also excluding bending elasticity of Fibers.
 
 This also sets the Lagrange multipliers for the Fiber.
 
 The function will not change the position of the Mecables.
 */
void Meca::computeForces()
{
    prepareMatrices();
    
    // calculate forces in vFOR, but without Brownian noise:
    calculateForces(vPTS, vBAS, vFOR);
    copy_real(dimension(), vFOR, vTMP);

    // add rigidity, and calculate Lagrange Multipliers:
    for ( Mecable * mec : mecables )
    {
        real * fff = vFOR + DIM * mec->matIndex();
        real * ttt = vTMP + DIM * mec->matIndex();
        
        // add bending rigidity:
#if SEPARATE_RIGIDITY_TERMS
        real * xxx = vPTS + DIM * mec->matIndex();
        mec->addRigidity(xxx, ttt);
#endif
        // calculate Lagrange multipliers for Fibers (irrelevant for other object)
        mec->projectForces(ttt, ttt);
        mec->storeTensions(ttt);
        // register force (all objects)
        mec->getForces(fff);
    }
}


/**
This preforms:
 
    vFOR <- vFOR + Noise
    vRHS <- time_step * P * vFOR:
 
 Also prepare Projection diff is requested

 'rhs' and 'fff' are output. Input 'rnd' is a set of independent random numbers
*/
real brownian1(Mecable* mec, real const* rnd, real alpha, real* fff, real beta, real* rhs)
{
    real n = mec->addBrownianForces(rnd, alpha, fff);
    
    // Calculate the right-hand-side of the system in vRHS:
    mec->projectForces(fff, rhs);
    
    // rhs <- beta * rhs, resulting in time_step * P * vFOR:
    blas::xscal(DIM*mec->nbPoints(), beta*mec->leftoverMobility(), rhs, 1);

    /*
     In this case, `fff` contains the true force in each vertex of the system
     and the Lagrange multipliers will represent the tension in the filaments
     */
    mec->storeTensions(fff);
    
#if ADD_PROJECTION_DIFF
    mec->makeProjectionDiff(fff);
#endif
    
    return n;
}


/**
 This solves the equation:
 
     ( Xnew - Xold ) / time_step = P * Force + Noise
 
 Explicit integration would lead to:
 
     Xnew = Xold + time_step * P * Force + Noise
 
 For implicit integration, we use a linearization of the force:
 
     Force(X) = M * X + B
 
 where M is a matrix and B is a vector, leading to:
 
     ( I - time_step * P * M ) ( Xnew - Xold ) = time_step * P * F + Noise
 
 where:
 
     F = M * Xold + B
 
     Noise = sqrt(2*kT*time_step*mobility) * Gaussian(0,1)
 
 Implicit integration ensures numerical stability. In the code,
 the matrix M is decomposed as:
 
     M = mB + mC + the rigidity terms of Mecables
 
 and the vectors are:
 
     'vPTS' is Xold
     'vBAS' is B
     'vRND' is Noise, a vector of calibrated Brownian terms
     'vFOR' is the force ( M * Xold + B ) and then ( M * Xnew + B )
     'vRHS' is the right-hand-side of the system (time_step * P * F + vRND)
     'vSOL' is `Xnew - Xold`, the solution to the system
 
 */
void Meca::solve(SimulProp const* prop, const unsigned precond)
{
    // get global time step
    time_step = prop->time_step;

    prepareMatrices();
    
    // calculate forces before constraints in vFOR:
    calculateForces(vPTS, vBAS, vFOR);
    
#if SEPARATE_RIGIDITY_TERMS
    addAllRigidity(vPTS, vFOR);
#endif
 
    /* 
     Fill `vRND` with Gaussian random numbers 
     This operation can be done in parallel, in a separate thread
     */
    RNG.gauss_set(vRND, dimension());
    
    /*
     Add Brownian motions to 'vFOR', and calculate vRHS by multiplying by mobilities.
     As Brownian terms are added, we record the magnitude of the typical smallest
     scalar contribution in `noiseLevel`. The dynamics will later be solved with 
     a residual that is proportional to this level:
     SimulProp::tolerance * noiseLevel
     As long as SimulProp::tolerance is smaller than 1, this should allow for a
     level of numerical error is small with respect to the Brownian noise in
     the system, and the results should be physically appropriate.
     */
    
    real noiseLevel = INFINITY;
    const real alpha = prop->kT/time_step;
    
    /*
     Add Brownian contributions and calculate Minimum value of it
      vFOR <- vFOR + Noise
      vRHS <- P * vFOR:
     */
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        real local = INFINITY;
        Mecable ** mci = mecables.begin() + omp_get_thread_num();
        while ( mci < mecables.end() )
        {
            const size_t inx = DIM * (*mci)->matIndex();
            real n = brownian1(*mci, vRND+inx, alpha, vFOR+inx, time_step, vRHS+inx);
            local = std::min(local, n);
            mci += NUM_THREADS;
        }
        //printf("thread %i min: %f\n", omp_get_thread_num(), local);
    #pragma omp critical
        noiseLevel = std::min(noiseLevel, local);
    }
#else
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
        real n = brownian1(mec, vRND+inx, alpha, vFOR+inx, time_step, vRHS+inx);
        noiseLevel = std::min(noiseLevel, n);
    }
#endif

    // scale minimum noise level to serve as a measure of required precision
    noiseLevel *= time_step;
    
    //printf("noiseLeveld = %8.2e   variance(vRHS) / estimate = %8.4f\n",
    //       noiseLevel, blas::nrm2(dimension(), vRHS) / (noiseLevel * sqrt(dimension())) );

#if NEW_CYTOPLASMIC_FLOW
    /**
     Includes a constant fluid flow displacing all the objects along
     */
    if ( prop->flow.norm() > REAL_EPSILON )
    {
        LOG_ONCE("NEW_CYTOPLASMIC_FLOW code enabled\n");
        Vector flow_dt = prop->flow * time_step;
        for ( int p = 0; p < dimension(); ++p )
            flow_dt.add_to(vRHS+DIM*p);
    }
#endif
    
#if EXPLICIT_INTEGRATION
    /*
     This implements the forward Euler integration, for testing purposes.
     The result is very inefficient, since we have built the entire stiffness matrix,
     which is not necessary for this simple explicit scheme.
     */
    blas::xaxpy(dimension(), 1.0, vRHS, 1, vPTS, 1);
    
    for ( Mecable * mec : mecables )
    {
        mec->getPoints(vPTS+DIM*mec->matIndex());
        mec->getForces(vFOR+DIM*mec->matIndex());
    }
    return;
#endif
    
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
     the tolerance is in scaled to the contribution of Brownian
     motions contained in vRHS, assuming that the error is equally spread 
     along all degrees of freedom, this should work for tolerance << 1
     here a printf() can be used to check that the estimate is correct:
    */ 
    
    // NOTE: the tolerance to solve the system should be such that the solution
    // found does not depend on the initial guess.
    
    real abstol = noiseLevel * prop->tolerance;
    
    if ( prop->kT == 0 )
        abstol = prop->tolerance;
    
    /*
     With exact arithmetic, biConjugate Gradient should converge at most
     in a number of iterations equal to the size of the linear system, with
     each iteration involving 2 matrix-vector multiplications.
     We set here the max limit to the number of matrix-vector multiplication:
     */
    LinearSolvers::Monitor monitor(2*dimension(), abstol);

    //------- call the iterative solver:

    if ( precond == 1 )
        computePreconditionner();
    else if ( precond == 2 )
        computePreconditionnerAlt();

    //fprintf(stderr, "Solve precond %i size %6i absolute tolerance %f\n", precond, dimension(), abstol);

    /*
     GMRES may converge faster than BCGGS, but has overheads and uses more memory
     hence for large systems, biCGGS is often advantageous.
     */

    if ( precond )
    {
        LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
        //LinearSolvers::GMRES(*this, vRHS, vSOL, 32, monitor, allocator, mH, mV, temporary);
    }
    else
    {
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
        //LinearSolvers::GMRES(*this, vRHS, vSOL, 64, monitor, allocator, mH, mV, temporary);
    }
    
#if ( 0 )
    //size_t mem = dimension() * (64+2) * sizeof(real);
    //std::clog << "GMRES mem " << (mem >> 20) << "MB\n";
    fprintf(stderr, "System size %6i precondition %i", dimension(), precond);
    fprintf(stderr, "    Solver count %4lu  residual %.3e\n", monitor.count(), monitor.residual());
#endif
#if ( 0 )
    // enable this to compare with GMRES using different restart parameters
    for ( int RS : {8, 16, 32} )
    {
        monitor.reset();
        zero_real(dimension(), vSOL);
        LinearSolvers::GMRES(*this, vRHS, vSOL, RS, monitor, allocator, mH, mV, temporary);
        fprintf(stderr, "    GMRES-%i  count %4lu  residual %.3e\n", RS, monitor.count(), monitor.residual());
    }
#endif
#if ( 0 )
    // enable this to compare BCGS and GMRES
    fprintf(stderr, "    BCGS     count %4lu  residual %.3e\n", monitor.count(), monitor.residual());
    monitor.reset();
    zero_real(dimension(), vSOL);
    LinearSolvers::GMRES(*this, vRHS, vSOL, 64, monitor, allocator, mH, mV, temporary);
    fprintf(stderr, "    GMRES-64 count %4lu  residual %.3e\n", monitor.count(), monitor.residual());
#endif
#if ( 0 )
    // enable this to compare with another implementation of biconjugate gradient stabilized
    monitor.reset();
    zero_real(dimension(), vSOL);
    LinearSolvers::bicgstab(*this, vRHS, vSOL, monitor, allocator);
    fprintf(stderr, "    bcgs      count %4lu  residual %.3e\n", monitor.count(), monitor.residual());
#endif
    
    //------- in case the solver did not converge, we try other methods:
    
    if ( !monitor.converged() )
    {
        Cytosim::out("Solver failed: precond %i flag %i count %4lu residual %.3e\n",
            precond, monitor.flag(), monitor.count(), monitor.residual());
        
        // try with different initial seed: vRHS
        monitor.reset();
        copy_real(dimension(), vRHS, vSOL);
        LinearSolvers::GMRES(*this, vRHS, vSOL, 255, monitor, allocator, mH, mV, temporary);
        Cytosim::out("     seed: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
        
        if ( !monitor.converged() )
        {
            monitor.reset();
            zero_real(dimension(), vSOL);
            
            if ( precond )
            {
                // try with the other method:
                LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                Cytosim::out("    BCGSP: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
            }
            else
            {
                // try with a preconditioner
                computePreconditionner();
                LinearSolvers::GMRES(*this, vRHS, vSOL, 127, monitor, allocator, mH, mV, temporary);
                Cytosim::out("    GMRES: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
                if ( !monitor.converged() )
                {
                    // try with other method:
                    monitor.reset();
                    zero_real(dimension(), vSOL);
                    LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                }
            }
            
            if ( !monitor.converged() )
            {
                // try with different GMRES parameters:
                monitor.reset();
                zero_real(dimension(), vSOL);
                LinearSolvers::GMRES(*this, vRHS, vSOL, 255, monitor, allocator, mH, mV, temporary);
                Cytosim::out("    GMRES(256): count %4lu residual %.3e\n", monitor.count(), monitor.residual());
                
                if ( !monitor.converged() )
                {
                    // no method could converge... this is really bad!
                    Exception e("convergence failure");
                    e << monitor.count() << " iterations, residual " << monitor.residual() << '\n';
                    throw e;
                }
            }
        }
    }
    
#ifndef NDEBUG
    
    //check validity of the data:
    for ( size_t i = 0; i < dimension(); ++i )
    {
        if ( std::isnan(vSOL[i]) )
        {
            fprintf(stderr, "Meca::solve produced NaN\n");
            abort();
        }
    }
    
#endif
    
    //add the solution of the system (=dPTS) to the points coordinates
    blas::add(dimension(), vSOL, vPTS);
    
    /*
     Re-calculate forces with the new coordinates, excluding bending elasticity
     and Brownian terms on the vertices.
     In this way the result returned to the fibers does not sum-up to zero,
     and is appropriate for example to calculate the effect of force on assembly.
     */
    calculateForces(vPTS, vBAS, vFOR);
    ready_ = 1;
    
    // report on the matrix type and size, sparsity, and the number of iterations
    if ( prop->verbose )
    {
        std::stringstream oss;
        oss << "Meca " << DIM << "*" << nbPts;
        oss << " kern " << largestMecable();
#if USE_ISO_MATRIX
        oss << " " << mB.what();
        if ( useMatrixC )
#endif
        oss << " " << mC.what();
        oss << " precond " << precond;
        oss << " count " << monitor.count();
        //oss << " flag " << monitor.flag();
        oss << " residual = " << monitor.residual() << "\n";
        Cytosim::out << oss.str();
        if ( prop->verbose > 1 )
            std::clog << oss.str();
    }
}


// transfer newly calculated point coordinates back to Mecables
void Meca::apply()
{
    if ( ready_ )
    {
        ready_ = 0;
        
#if NUM_THREADS > 1
#pragma omp parallel num_threads(NUM_THREADS)
        {
            Mecable ** mci = mecables.begin() + omp_get_thread_num();
            while ( mci < mecables.end() )
            {
                Mecable * mec = *mci;
                mec->getForces(vFOR+DIM*mec->matIndex());
                mec->getPoints(vPTS+DIM*mec->matIndex());
                mci += NUM_THREADS;
            }
        }
#else
        for ( Mecable * mec : mecables )
        {
            mec->getForces(vFOR+DIM*mec->matIndex());
            mec->getPoints(vPTS+DIM*mec->matIndex());
        }
#endif
    }
    else
    {
        //write(2, "extra Meca::apply() calls\n", 26);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Debug/Output Functions


/**
 Count number of non-zero entries in the entire system
 */
size_t Meca::nbNonZeros(real threshold) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * dst = new_real(dim);
    
    zero_real(dim, src);
    
    size_t cnt = 0;
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1.0;
        multiply(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            cnt += ( std::abs(dst[i]) > threshold );
        src[j] = 0.0;
    }
    
    free_real(dst);
    free_real(src);
    return cnt;
}

/**
 Extract the full matrix associated with matVect, in `mat[]`
 The matrix should be preallocated of size `dim`, which should be equal
 to the system's dimension().
 */
void Meca::getMatrix(size_t dim, real * mat) const
{
    if ( dim != dimension() )
        throw InvalidIO("wrong matrix dimension");
    
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    zero_real(dim, res);
    
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1.0;
        multiply(src, res);
        blas::xcopy(dim, res, 1, mat+j*dim, 1);
        src[j] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save a sparse matrix in Matrix Market format
 https://math.nist.gov/MatrixMarket/formats.html
 This is a Sparse text format
 */
void Meca::saveMatrix(FILE * file, real threshold) const
{
    fprintf(file, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(file, "%% This is a matrix produced by Cytosim\n");
    fprintf(file, "%% author: FJ Nedelec\n");
    fprintf(file, "%% kind: biological cell simulation (cytoskeleton)\n");

    const size_t dim = dimension();
    real * src = new_real(dim);
    real * dst = new_real(dim);
    zero_real(dim, src);
    
    const size_t cnt = nbNonZeros(threshold);
    fprintf(file, "%lu %lu %lu\n", dim, dim, cnt);
    
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1.0;
        multiply(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            if ( std::abs(dst[i]) > threshold )
                fprintf(file, "%lu %lu %f\n", i, j, dst[i]);
        src[j] = 0.0;
    }
    
    free_real(dst);
    free_real(src);
}


void Meca::saveRHS(FILE * file) const
{
    fprintf(file, "%% This is a vector produced by Cytosim\n");
    fprintf(file, "%% author: FJ Nedelec\n");
    fprintf(file, "%% kind: biological cell simulation (cytoskeleton)\n");
    
    const size_t dim = dimension();
    fprintf(file, "%lu\n", dim);
    
    for ( size_t i = 0; i < dim; ++i )
        fprintf(file, "%f\n", vRHS[i]);
}


/**
 Save the full matrix associated with multiply(), in binary format
 */
void Meca::dumpMatrix(FILE * file) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        multiply(src, res);
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the elasticity matrix
 */
void Meca::dumpElasticity(FILE * file) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        
        calculateForces(src, nullptr, res);
        
#if SEPARATE_RIGIDITY_TERMS
        addAllRigidity(src, res);
#endif

#if ADD_PROJECTION_DIFF
        for ( Mecable const* mec : mecables )
        {
            if ( mec->hasProjectionDiff() )
            {
                const size_t inx = DIM * mec->matIndex();
                mec->addProjectionDiff(src+inx, res+inx);
            }
        }
#endif
        
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the projection matrix multiplied by the mobility
 */
void Meca::dumpMobility(FILE * file) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t i = 0; i < dim; ++i )
    {
        src[i] = 1.0;
        zero_real(dim, res); // this should not be necessary
        
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            // this includes the mobility, but not the time_step:
            mec->projectForces(src+inx, res+inx);
            blas::xscal(DIM*mec->nbPoints(), mec->leftoverMobility(), res+inx, 1);
        }
        // write column to file directly:
        fwrite(res, sizeof(real), dim, file);
        src[i] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save matrix associated with the preconditionner, in binary format
 This relies on `Meca::precondition()`, which may apply a dummy preconditionner
 */
void Meca::dumpPreconditionner(FILE * file) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        precondition(src, res);
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


void Meca::dumpObjectID(FILE * file) const
{
    real * vec = new_real(largestMecable());
    
    for ( size_t ii = 0; ii < mecables.size(); ++ii )
    {
        const size_t nbp = mecables[ii]->nbPoints();
        for ( size_t p=0; p < nbp; ++p )
            vec[p] = ii;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), nbp, file);
    }
    
    free_real(vec);
}


void Meca::dumpDrag(FILE * file) const
{
    real * vec = new_real(largestMecable());
    
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = mec->nbPoints();
        const real drag = mec->dragCoefficient() / nbp;
        for ( size_t p=0; p < nbp; ++p )
            vec[p] = drag;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), nbp, file);
    }
    
    free_real(vec);
}


/**
 This dump the total matrix and some vectors in binary files.
 
 This MATLAB code should read the output:
 
     ord = load('ord.txt');
     time_step = load('stp.txt');
     obj = fread(fopen('obj.bin'), ord, 'double');
     drg = fread(fopen('drg.bin'), ord, 'double');
     sys = fread(fopen('sys.bin'), [ord, ord], 'double');
     ela = fread(fopen('ela.bin'), [ord, ord], 'double');
     mob = fread(fopen('mob.bin'), [ord, ord], 'double');
     con = fread(fopen('con.bin'), [ord, ord], 'double');
     pts = fread(fopen('pts.bin'), ord, 'double');
     rhs = fread(fopen('rhs.bin'), ord, 'double');
     sol = fread(fopen('sol.bin'), ord, 'double');
 
 To display the matrices:

     imshow(abs(sys))
     imshow(abs(ela))
 
 You can then compare the results with matlab's own iterative method,
 and compare the result using a scatter plot:
 
     x = bicgstab(sys, rhs, 0.001, ord);
     plot(x, sol, '.');
 
 */
void Meca::dump() const
{
    FILE * f = fopen("ord.txt", "w");
    fprintf(f, "%lu\n", dimension());
    fclose(f);
    
    f = fopen("stp.txt", "w");
    fprintf(f, "%f\n", time_step);
    fclose(f);
 
    f = fopen("drg.bin", "wb");
    dumpDrag(f);
    fclose(f);
    
    f = fopen("obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    f = fopen("rhs.bin", "wb");
    fwrite(vRHS, sizeof(real), dimension(), f);
    fclose(f);
    
    f = fopen("sol.bin", "wb");
    fwrite(vSOL, sizeof(real), dimension(), f);
    fclose(f);
    
    f = fopen("pts.bin", "wb");
    fwrite(vPTS, sizeof(real), dimension(), f);
    fclose(f);
    
    f = fopen("sys.bin", "wb");
    dumpMatrix(f);
    fclose(f);
    
    f = fopen("ela.bin", "wb");
    dumpElasticity(f);
    fclose(f);
    
    f = fopen("mob.bin", "wb");
    dumpMobility(f);
    fclose(f);
    
    f = fopen("con.bin", "wb");
    dumpPreconditionner(f);
    fclose(f);
}


/**
 output vectors and matrices in a text-based sparse format
 */
void Meca::dumpSparse()
{
#if SEPARATE_RIGIDITY_TERMS
    std::clog << "incorrect dump since SEPARATE_RIGIDITY_TERMS is defined\n";
#endif
    std::clog << "dumping matrices in binary format\n";
        
    FILE * f = fopen("d_drg.bin", "wb");
    dumpDrag(f);
    fclose(f);
    
    f = fopen("d_obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    std::ofstream os("d_sol.txt");
    VecPrint::dump(os, dimension(), vPTS);
    os.close();
    
    os.open("d_rhs.txt");
    VecPrint::dump(os, dimension(), vRHS);
    os.close();
    
#if USE_ISO_MATRIX
    os.open("d_matB.txt");
    mB.printSparse(os, 0);
    os.close();
#endif
    
    os.open("d_matC.txt");
    mC.printSparse(os, 0);
    os.close();
        
    size_t alc = 0;
    for ( Mecable const* mec : mecables )
        alc = std::max(alc, mec->nbPoints());

    real * tmp1 = new_real(DIM*alc);
    real * tmp2 = new_real(DIM*DIM*alc*alc);
    
    os.open("diagonal.txt");
    
    for ( Mecable * mec : mecables )
    {
        const size_t bks = DIM * mec->nbPoints();
        extractBlock(tmp2, mec);
        VecPrint::sparse_off(os, bks, bks, tmp2, bks, DIM*mec->matIndex());
    }
    os.close();
    
    free_real(tmp1);
    free_real(tmp2);
}

