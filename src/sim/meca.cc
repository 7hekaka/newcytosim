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

#include "assert_macro.h"
#include "blas.h"
#include "lapack.h"
#include "cytoblas.h"
#include "xtbsv.h"
#include "xtrsm.h"

#include "meca.h"
#include "mecable.h"
#include "messages.h"
#include "simul_prop.h"
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


/// use this to generate code for validation
#define DEBUG_MECA 0


/// if you like this famous dish from Alsace, set this to a prime number
#define CHOUCROUTE 7

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
: mecables(32, 32)
{
    ready_ = -1;
    nPoints_ = 0;
    allocated_ = 0;
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
    tau_ = 0;
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
        
        // pad with 4 doubles to allow some SIMD instruction burr
        alc = DIM * allocated_ + 4;
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
 
     F <- B + mISO * X + mFUL * X

 If `B == 0`, this term is ommited. With `B = vBAS` and `X = vPTS`, this
 function calculates the forces in the system in `F`:
 
     F <- vBAS + mISO * vPTS + mFUL * vPTS

 */
void Meca::calculateForces(const real* X, real const* B, real* F) const
{
    assert_true( empty() || ( X != F && X != B && F != B ));
    
    // F <- B
    copy_real(dimension(), B, F);      //blas::xcopy(dimension(), B, 1, F, 1);
    
#if USE_ISO_MATRIX
    // F <- F + mISO * X
    mISO.VECMULADDISO(X, F);

    if ( useFullMatrix )
#endif
    {
        // F <- F + mFUL * X
        mFUL.vecMulAdd(X, F);
    }
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


// the leading dimension of the banded matrix used for preconditionning
constexpr size_t BAND_LDD = 3;

/// apply preconditionner block in band storage
inline void applyPrecondBand(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
#if CHOUCROUTE
    alsatian_xpbtrsL<DIM>(nbp, mec->block(), BAND_LDD, Y);
#else
    /*
     we cannot call lapack::DPBTRS('L', bks, KD, 1, mec->block(), KD+1, Y, bks, &info)
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     But calling DTBSV gets the required work done.
     */
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', nbp, 2, mec->block(), BAND_LDD, Y+d, DIM);
        blas::xtbsv('L', 'T', 'N', nbp, 2, mec->block(), BAND_LDD, Y+d, DIM);
    }
#endif
}


/// apply isotropic preconditionner block in full storage
inline void applyPrecondIsoS(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
    /*
     we cannot call lapack::DPOTRS('L', nbp, mec->block(), nbp, Y, DIM, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     */
#if 0
    real * tmp = new_real(DIM*nbp);
    copy_real(DIM*nbp, Y, tmp);
    iso_xpotrsL<DIM>(nbp, mec->block(), nbp, tmp);
    size_t S = std::min(16, nbp);
    std::clog << "\n "; VecPrint::print(std::clog, S, tmp, 3, 100.0);
    free_real(tmp);
#endif
#if CHOUCROUTE
    alsatian_xpotrsL<DIM>(nbp, mec->block(), nbp, Y);
#else
    iso_xpotrsL<DIM>(nbp, mec->block(), nbp, Y);
#endif
    //std::clog << "\nL"; VecPrint::print(std::clog, S, Y, 3, 100.0);
}


/// apply isotropic preconditionner block in full storage
inline void applyPrecondIsoP(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
    /*
     we cannot call lapack::DPOTRS('L', nbp, mec->block(), nbp, Y, DIM, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     */
    //iso_xgetrsN<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
    alsatian_xgetrsN<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
}


/// apply preconditionner block in full storage
inline void applyPrecondFull(Mecable const* mec, real* Y)
{
    int bks = mec->blockSize();
    int info = 0;
    lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y, bks, &info);
    //lapack_xgetrsN(bks, mec->block(), bks, mec->pivot(), Y);
    assert_true(info==0);
}


inline void applyPreconditionner(Mecable const* mec, real* Y)
{
    switch ( mec->blockType() )
    {
        case 0: break;
        case 1: applyPrecondBand(mec, Y); break;
        case 2: applyPrecondIsoS(mec, Y); break;
        case 3: applyPrecondIsoP(mec, Y); break;
        case 4: applyPrecondFull(mec, Y); break;
        default: ABORT_NOW("unknown Mecable::blockType()"); break;
    }
}


inline void applyPreconditionner(Mecable const* mec, real const* X, real* Y)
{
    assert_true( Y != X );
    //if ( Y != X ) blas::xcopy(bks, X, 1, Y, 1);
    switch ( mec->blockType() )
    {
        case 0: break;
        case 1: applyPrecondBand(mec, Y); break;
        case 2: applyPrecondIsoS(mec, Y); break;
        case 3: applyPrecondIsoP(mec, Y); break;
        case 4: applyPrecondFull(mec, Y); break;
#if EXPERIMENTAL_PRECONDITIONNERS
        case 7: mec->blockMultiply(X, Y); break;
#endif
        default: ABORT_NOW("unknown Mecable::blockType()"); break;
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
    const size_t bks = DIM * mec->nbPoints();

    // multiply the lines corresponding to this Mecable:
    mFUL.vecMul(X, Y, inx, inx+bks);

#if SEPARATE_RIGIDITY_TERMS
    mec->addRigidity(X+inx, Y+inx);
#endif

#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(X+inx, Y+inx);
#endif

    mec->projectForces(Y+inx, Y+inx);

    // Y <- X + alpha * Y
    blas::xpay(bks, X+inx, -tau_*mec->leftoverMobility(), Y+inx);
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
    #if SEPARATE_RIGIDITY_TERMS
        mec->addRigidity(X+inx, Y+inx);
    #endif
    #if ADD_PROJECTION_DIFF
        if ( mec->hasProjectionDiff() )
            mec->addProjectionDiff(X+inx, Y+inx);
    #endif
        mec->projectForces(Y+inx, Y+inx);
        // Y <- X + alpha * Y
        blas::xpay(DIM*mec->nbPoints(), X+inx, -tau_*mec->leftoverMobility(), Y+inx);
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

    mFUL.vecMul(X, T, inx, inx+bks);
    
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
        blas::xaxpy(bks, -tau_*mec->leftoverMobility(), T+inx, 1, Y+inx, 1);
    }
    else
    {
        // Y <- P * T
        mec->projectForces(T+inx, Y+inx);
        // Y <- X + alpha * Y = X + alpha * P * FORCE
        blas::xpay(bks, X+inx, -tau_*mec->leftoverMobility(), Y+inx);
    }

    // Y <- PRECONDITIONNER_BLOCK * Y
    applyPreconditionner(mec, Y+inx);
}


/**
 This is equivalent to
 
     multiply(X, T);       // T <- M*X
     precondition(T, Y);   // Y <- P*T
 
 This is used for left-sided preconditionning
 */
void Meca::multiply_precondition(real const* X, real* T, real* Y) const
{
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable * mec : mecables )
        multiply_precondition1(mec, X, T, Y);
}


void Meca::multiply_precondition1(Mecable const* mec, real const* X, real* Y) const
{
    assert_true( X != Y );
    
    const size_t inx = DIM * mec->matIndex();
    const size_t bks = DIM * mec->nbPoints();

    mFUL.vecMul(X, Y, inx, inx+bks);
    
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
    blas::xpay(bks, X+inx, -tau_*mec->leftoverMobility(), Y+inx);
    
    applyPreconditionner(mec, Y+inx);
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
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable * mec : mecables )
        multiply_precondition1(mec, X, Y);
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
        // Y <- X + alpha * Y
        const real beta = -tau_ * mec->leftoverMobility();
        blas::xpay(DIM*mec->nbPoints(), X+inx, beta, Y+inx);
    }
}

#endif  // PARALLELIZE_MATRIX

//------------------------------------------------------------------------------
#pragma mark - Precondition


void printBlock(Mecable* mec, size_t sup)
{
    const size_t bks = DIM * mec->nbPoints();
    size_t S = std::min(bks, sup);
    real * blk = mec->block();

    std::clog << "Diagonal block " << bks << "\n";
    VecPrint::print(std::clog, S, S, mec->block(), bks, 3);
    
    std::clog << "S ";
    char str[32];
    for ( size_t i = 0; i < S; ++i )
    {
        real sum = 0;
        for ( size_t j = 0; j < bks; ++j )
            sum += blk[j+bks*i];
        snprintf(str, sizeof(str), " %8.3f", sum);
        std::clog << str;
    }
    std::clog << "\n";
}


/**
 Build a banded block that can be used as a preconditionner:
 
     I - time_step * mobility * Rigidity
 
 This block is square, symmetric, definite positive and totally predictable
 */
void Meca::getBandedBlock(const Mecable * mec, real* res) const
{
    const size_t nbp = mec->nbPoints();

    const real beta = -tau_ * mec->leftoverMobility();
    
    if ( mec->hasRigidity() )
    {
        if ( BAND_LDD != 3 )
            zero_real(BAND_LDD*nbp, res);
        setRigidityBanded(res, BAND_LDD, nbp, beta*mec->fiberRigidity());
    }
    else
        zero_real(BAND_LDD*nbp, res);

    // add ones to the diagonal:
    for ( size_t i = 0; i < nbp; ++i )
        res[i*BAND_LDD] += 1.0;
    
    //std::clog << "\nrigidity band " << nbp << "\n";
    //VecPrint::print(std::clog, 3, std::min(nbp, 16UL), res, BAND_LDD, 1);

#if USE_ISO_MATRIX
    mISO.addTriangularBlockBanded(beta, res, BAND_LDD, mec->matIndex(), nbp, 2);
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalTraceBanded(beta/DIM, res, BAND_LDD, DIM*mec->matIndex(), DIM*nbp, 2);
}


/**
 Get the total diagonal block corresponding to an Object, which is:
 
     I - time_step * P ( mISO + mFUL + P' )
 
 The result is constructed by using functions from mISO and mFUL
 This block is square but not symmetric!
 */
void Meca::getIsoBlock(const Mecable * mec, real* res) const
{
    const size_t nbp = mec->nbPoints();
    
    zero_real(nbp*nbp, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    if ( mec->hasRigidity() )
        //addRigidityLower<1>(res, nbp, mec->nbPoints(), mec->fiberRigidity());
        addRigidity<1>(res, nbp, mec->nbPoints(), mec->fiberRigidity());
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, nbp, mec->matIndex(), nbp);
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalTrace(1.0/DIM, res, nbp, DIM*mec->matIndex(), DIM*nbp);
#if ( 0 )
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
    {
        // Include the corrections P' in preconditioner, vector by vector.
        real* tmp = vTMP + DIM * mec->matIndex();
        zero_real(bks, tmp);
        for ( size_t i = 0; i < bks; ++i )
        {
            tmp[i] = 1.0;
            mec->addProjectionDiff(tmp, res+bks*i);
            tmp[i] = 0.0;
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
    
    for ( size_t i = 0; i < bks; ++i )
        mec->projectForces(res+bks*i, res+bks*i);
#endif

    // scale
    const real beta = -tau_ * mec->leftoverMobility();
    //blas::xscal(bs*bs, beta, res, 1);
    for ( size_t n = 0; n < nbp*nbp; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    size_t bs1 = nbp + 1;
    for ( size_t i = 0; i < nbp*nbp; i += bs1 )
        res[i] += 1.0;
}


/**
 Get the total diagonal block corresponding to an Object, which is:
 
     I - time_step * P ( mISO + mFUL + P' )
 
 The result is constructed by using functions from mISO and mFUL
 This block is square but not symmetric!
 */
void Meca::getFullBlock(const Mecable * mec, real* res) const
{
    const size_t nbp = mec->nbPoints();
    const size_t bks = DIM * nbp;
    
    zero_real(bks*bks, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    if ( mec->hasRigidity() )
    {
        addRigidityLower<DIM>(res, bks, mec->nbPoints(), mec->fiberRigidity());
        //std::clog<<"Rigidity block " << mec->reference() << "\n";
        //VecPrint::print(std::clog, bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
#endif
    expand_lower_matrix<DIM>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, DIM*mec->matIndex(), bks);
    
#if ( 0 )
    std::clog<<"mISO+mFUL block:\n";
    VecPrint::print(std::clog, bks, bks, res, bks);
#endif
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
    {
        // Include the corrections P' in preconditioner, vector by vector.
        real* tmp = vTMP + DIM * mec->matIndex();
        zero_real(bks, tmp);
        for ( size_t i = 0; i < bks; ++i )
        {
            tmp[i] = 1.0;
            mec->addProjectionDiff(tmp, res+bks*i);
            tmp[i] = 0.0;
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
    
    for ( size_t i = 0; i < bks; ++i )
        mec->projectForces(res+bks*i, res+bks*i);
    
    // scale
    const real beta = -tau_ * mec->leftoverMobility();
    //blas::xscal(bs*bs, beta, res, 1);
    for ( size_t n = 0; n < bks*bks; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    size_t bs1 = bks + 1;
    for ( size_t i = 0; i < bks*bks; i += bs1 )
        res[i] += 1.0;
}


/**
 This version builds the diagonal block indirectly using Meca:multiply().
 This is a slow method that calls 'multiply()' n-times, where
 'n' is the size of the block.
 
 This should be used for validation only.
*/
void Meca::extractBlock(const Mecable* mec, real* res) const
{
    const size_t dim = dimension();
    const size_t bks = DIM * mec->nbPoints();
    const size_t off = DIM * mec->matIndex();
    
    assert_true( off+bks <= dim );
    real * vec = new_real(dim);
    real * tmp = new_real(dim);
    
    zero_real(dim, vec);
    //zero_real(bks*bks, res);
    
    // proceed column by column:
    for ( size_t jj = 0; jj < bks; ++jj )
    {
        vec[jj+off] = 1.0;
        multiply(vec, tmp);
        vec[jj+off] = 0.0;
        copy_real(bks, tmp+off, res+jj*bks);
    }
    
    free_real(vec);
    free_real(tmp);
}


/**
 DEBUG: compare `blk` with block extracted using extractBlock()
 */
void Meca::verifyBlock(const Mecable * mec, const real* blk)
{
    const size_t bks = DIM * mec->nbPoints();
    real * wrk = new_real(bks*bks);
    
    extractBlock(mec, wrk);
    
    blas::xaxpy(bks*bks, -1.0, blk, 1, wrk, 1);
    real err = blas::nrm2(bks*bks, wrk);
 
    std::clog << "verifyBlock ";
    std::clog << std::setw(8) << mec->reference() << "  size " << std::setw(4) << bks;
    std::clog << " | B - K | = " << std::setprecision(3) << err << '\n';
    
    if ( err > bks * bks * REAL_EPSILON )
    {
        //VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        extractBlock(mec, wrk);
        
        const size_t S = std::min(9UL, bks);
        std::clog << " extracted:\n"; VecPrint::print(std::clog, S, S, wrk, bks, 3);
        std::clog << " computed:\n";  VecPrint::print(std::clog, S, S, blk, bks, 3);
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
    
    if ( 0 == mec->blockType() )
        return;

    std::clog << "  checkBlock " << mec->blockType() << " ";
    std::clog << std::setw(10) << mec->reference() << " " << std::setw(6) << bks;
    
    real * wrk = new_real(bks*bks);
    real * mat = new_real(bks*bks);
    real * vec = new_real(bks);
    
    extractBlock(mec, wrk);   // wrk <- MEC_BLOCK
   
    copy_real(bks*bks, wrk, mat);  // mat <- wrk
    for ( size_t i = 0; i < bks; ++i )
    {
        applyPreconditionner(mec, wrk+bks*i, mat+bks*i);
        mat[i+bks*i] -= 1.0;
    }
    real err = blas::nrm2(bks*bks, mat) / bks;
    std::clog << " | 1 - PM | = " << std::setprecision(3) << err;
    
    if ( 1 )
    {
        // chose initial vector for power iteration
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        real eig = -1;
        if ( 1 == mec->blockType() )
            eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, mat);
        std::clog << "  eigen(1-PM) = " << eig;
    }
    
    std::clog << '\n';

    if ( err > 1 )
    {
        // print preconditionner block for visual inspection:
        const size_t S = std::min(16UL, bks);
        std::clog << " matrix: \n";
        VecPrint::print(std::clog, S, S, wrk, bks);
        std::clog << "\nprecond: \n";
        VecPrint::print(std::clog, S, S, blk, bks);
        std::clog << "\nprecond * matrix:\n";
        VecPrint::print(std::clog, S, S, mat, bks);
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
void Meca::computePrecondAlt(Mecable* mec, real* tmp, real* wrk, size_t wrksize)
{
    const size_t bks = DIM * mec->nbPoints();
    assert_true( bks*bks <= wrksize );

    bool may_keep = false;
    
    if ( mec->blockSize() == bks )
        may_keep = true;
    else
        mec->blockSize(bks, bks*bks, bks);
    
    real* blk = mec->block();
    real* vec = vTMP + DIM * mec->matIndex();

    if ( may_keep )
    {
        // extract diagonal matrix block corresponding to this Mecable:
        getFullBlock(mec, wrk);
        
        // chose initial vector for power iteration
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        
        // power iterations scale like (bs^2)
        real eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, tmp);
 
        //std::clog << "mec " << std::setw(4) << mec->identity();
        if ( eig < 1.0 )
        {
            //std::clog << "    keep     eigen(1-PM) = " << eig << '\n';
            mec->blockType(4);
            return;
        }
        else
        {
            //std::clog << "    RENEW    eigen(1-PM) = " << eig << '\n';
            blas::xcopy(bks*bks, wrk, 1, blk, 1);
        }
    }
    else
    {
        // extract diagonal matrix block corresponding to this Mecable:
        getFullBlock(mec, blk);
    }
    
    //verifyBlock(mec, blk);

    int info = 0;
    
    /**
     We invert here the full block using LAPACK general matrix functions
     the workload scales as SIZE ^ 3, as a function of the size of the block
     */
    //TEST truncate_matrix(bks, blk, 3*DIM-1);
    
    // calculate LU factorization:
    lapack::xgetf2(bks, bks, blk, bks, mec->pivot(), &info);
    
    if ( info )
    {
        //failed to factorize matrix !!!
        std::clog << "Meca::computePrecond failed (lapack::xgetf2, info = " << info << ")\n";
        return;
    }
    
    mec->blockType(4);
    //checkBlock(mec, blk);
}


/**
 Compute sequentially all the blocks of the preconditionner
 */
void Meca::computePrecondAlt()
{
    const size_t vecsize = DIM * largestMecable();
    const size_t wrksize = vecsize * vecsize;
 
    // allocate memory:
    real * tmp = new_real(wrksize);
    real * wrk = new_real(wrksize);
    
    for ( Mecable * mec : mecables )
        computePrecondAlt(mec, tmp, wrk, wrksize);
    
    free_real(wrk);
    free_real(tmp);
}


/**
 Compute a preconditionner block corresponding to 'mec'
 The dimension is reduced by DIM and banded with diagonal + 2 off-diagonals
 This block is symmetric definite positive, and is factorized by Cholesky's method
 */
void Meca::computePrecondBand(Mecable* mec)
{
    assert_true(BAND_LDD>2);
    const size_t nbp = mec->nbPoints();
    mec->blockSize(DIM*nbp, BAND_LDD*nbp, 0);
    getBandedBlock(mec, mec->block());
    
    //std::clog << "banded preconditionner " << nbp << "\n";
    //VecPrint::print(std::clog, 3, std::min(nbp, 16UL), mec->block(), BAND_LDD, 1);
    
    /**
     Factorize banded matrix with Andre-Louis Cholesky's method
     born 15.10.1875 in Montguyon, France
     died 31.08.1918 in Bagneux, following wounds received in battlefield.
     */
    int info = 0;
#if CHOUCROUTE
    alsatian_xpbtf2L<2>(nbp, mec->block(), BAND_LDD, &info);
#else
    lapack::xpbtf2('L', nbp, 2, mec->block(), BAND_LDD, &info);
#endif
    if ( 0 == info )
    {
        mec->blockType(1);
        //std::clog<<"factorized banded preconditionner: " << nbp << "\n";
        //VecPrint::print(std::clog, 3, std::min(nbp, 16UL), mec->block(), BAND_LDD, 1);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute Band Preconditionner block\n";
        ++bump_;
    }
}


/**
 Compute a preconditionner block corresponding to 'mec':
 Block of dimension reduced by DIM and without projection
 This block is symmetric definite positive, and is factorized by Cholesky's method
 */
void Meca::computePrecondIsoS(Mecable* mec)
{
    const size_t nbp = mec->nbPoints();
    //const size_t S = std::min(nbp, 16UL);
#if 0
    const size_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, bks);
    getFullBlock(mec, mec->block());
    project_matrix<DIM>(nbp, mec->block(), bks, mec->block(), nbp);
    std::clog<<"projected: " << bks << "\n";
    VecPrint::print(std::clog, S, S, mec->block(), nbp, 2);
#endif

    mec->blockSize(DIM*nbp, nbp*nbp, 0);
    
    //getBandedBlock(mec, mec->block());
    //std::clog << "banded preconditionner " << nbp << "\n";
    //VecPrint::print(std::clog, 3, S, mec->block(), BAND_LDD, 1);

    getIsoBlock(mec, mec->block());
    
    //std::clog<<"iso symmetric preconditionner: " << nbp << "\n";
    //VecPrint::print(std::clog, S, S, mec->block(), nbp, 2);

    // calculate Cholesky factorization:
    int info = 0;
#if CHOUCROUTE
    alsatian_xpotf2L(nbp, mec->block(), nbp, &info);
#else
    lapack::xpotf2('L', nbp, mec->block(), nbp, &info);
#endif
    
    if ( 0 == info )
    {
        mec->blockType(2);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute IsoS Preconditionner block\n";
        ++bump_;
    }
}

/**
Compute a preconditionner block corresponding to 'mec':
 Block of dimension reduced by DIM, with projection
 The block is not symmetric and is factorized by LU decomposition
 */
void Meca::computePrecondIsoP(Mecable* mec)
{
    const size_t nbp = mec->nbPoints();
    int info = 0;

    const size_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, nbp);
    getFullBlock(mec, mec->block());
    project_matrix<DIM>(nbp, mec->block(), bks, mec->block(), nbp);

    //const size_t S = std::min(nbp, 16UL);
    //std::clog<<"isop preconditionner: " << nbp << "\n";
    //VecPrint::print(std::clog, S, S, mec->block(), nbp, 2);

    // calculate LU factorization:
    lapack::xgetf2(nbp, nbp, mec->block(), nbp, mec->pivot(), &info);

    if ( 0 == info )
    {
        mec->blockType(3);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute IsoP Preconditionner block\n";
        ++bump_;
    }
}

/**
Compute preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondFull(Mecable* mec)
{
    const size_t bks = DIM * mec->nbPoints();
    mec->blockSize(bks, bks*bks, bks);
    
    getFullBlock(mec, mec->block());
    //verifyBlock(mec, mec->block());
    
    // calculate LU factorization:
    int info = 0;
    lapack::xgetf2(bks, bks, mec->block(), bks, mec->pivot(), &info);
    
    if ( 0 == info )
    {
        mec->blockType(4);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute full Preconditionner block\n";
        ++bump_;
    }
}


/**
 Initialize Mecable::blockMatrix() as the preconditionner block,
 using the factorization given by 'blk' and 'piv'
 */
void convertPreconditionner(Mecable* mec, real* blk, int* piv, real* wrk)
{
#if EXPERIMENTAL_PRECONDITIONNERS
    const size_t bks = DIM * mec->nbPoints();
    // create Identity matrix:
    init_matrix(bks, wrk, 1, 0);
    int info = 0;
    // calculate matrix inverse
    lapack::xgetrs('N', bks, bks, blk, bks, piv, wrk, bks, &info);
    if ( info == 0 )
    {
        MatrixFull& mat = mec->blockMatrix();
        mat.resize(bks);
        mat.importMatrix(bks, wrk, bks);
        mec->blockSize(bks, 0, 0);
        mec->blockType(7);
        //std::clog << "R";
    }
#endif
}


/**
renew preconditionner block corresponding to 'mec'
 */
void Meca::renewPreconditionner(Mecable* mec, int span, real* blk, int* piv, real* wrk, size_t wrksize)
{
#if EXPERIMENTAL_PRECONDITIONNERS
    const size_t bks = DIM * mec->nbPoints();
    const int age = mec->blockAge();

    if (( mec->blockType() != 7 ) | ( mec->blockSize() != bks ) | ( age >= span ))
    {
        // recalculate block!
        getFullBlock(mec, blk);
        int info = 0;
        // factorize block
        lapack::xgetf2(bks, bks, blk, bks, piv, &info);
        if ( info == 0 )
            convertPreconditionner(mec, blk, piv, wrk);
        else
            mec->blockType(0);
    }
    else
        mec->blockAge(age+1);
#endif
}



// Compute sequentially all the blocks of the preconditionner
void Meca::renewPreconditionner(int span)
{
    const size_t vecsize = DIM * largestMecable();
    const size_t wrksize = vecsize * vecsize;
 
    // allocate memory:
    real * tmp = new_real(wrksize);
    real * wrk = new_real(wrksize);
    int* pivot = new int[vecsize];

    for ( Mecable * mec : mecables )
        renewPreconditionner(mec, span, tmp, pivot, wrk, wrksize);
    
    delete[] pivot;
    free_real(wrk);
    free_real(tmp);
}


void Meca::computePreconditionner(int precond, int span)
{
    bump_ = 0;
    switch( precond )
    {
        case 0:
            for ( Mecable * mec : mecables )
                mec->blockType(0);
            break;
        case 1:
            for ( Mecable * mec : mecables )
                computePrecondBand(mec);
            break;
        case 2:
            for ( Mecable * mec : mecables )
                computePrecondIsoS(mec);
            break;
        case 3:
            for ( Mecable * mec : mecables )
                computePrecondIsoP(mec);
            break;
        case 4:
            for ( Mecable * mec : mecables )
                computePrecondFull(mec);
            break;
        case 5:
            renewPreconditionner(span);
            break;
        default:
            throw InvalidParameter("unknown `precondition' value");
            break;
    }
    if ( bump_ > 0 )
        Cytosim::log << "failed to compute " << bump_ << " / " << mecables.size() << " preconditionner blocks\n";
}


void Meca::precondition(const real* X, real* Y) const
{
    auto rdt = __rdtsc();
    if ( Y != X )
    {
        blas::xcopy(dimension(), X, 1, Y, 1);
        #pragma omp parallel for num_threads(NUM_THREADS)
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            applyPreconditionner(mec, X+inx, Y+inx);
        }
    }
    else
    {
        #pragma omp parallel for num_threads(NUM_THREADS)
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            applyPreconditionner(mec, Y+inx);
        }
    }
    cycles_ += __rdtsc() - rdt;
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/// function to sort Mecables
int compareMecables(const void * A, const void * B)
{
    size_t a = (*static_cast<Mecable *const*>(A))->nbPoints();
    size_t b = (*static_cast<Mecable *const*>(B))->nbPoints();
    return ( a < b ) - ( a > b );
}


/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::prepare(Simul const* sim)
{
    ready_ = 0;
    mecables.clear();

    for ( Fiber  * f= sim->fibers.first(); f ; f=f->next() )
        addMecable(f);
    for ( Solid  * s= sim->solids.first(); s ; s=s->next() )
        addMecable(s);
    for ( Sphere * o=sim->spheres.first(); o ; o=o->next() )
        addMecable(o);
    for ( Bead   * b=  sim->beads.first(); b ; b=b->next() )
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
    allocate(cnt);
    
#if USE_ISO_MATRIX
    //allocate the sparse matrices:
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
            addRigidityMatrix(mISO, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   else
            addRigidityBlockMatrix<DIM>(mFUL, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   endif
        }
#endif
    }
    //fprintf(stderr, "Meca::prepare() isnan %i\n", isnan(dimension(), vPTS));
}


/**
 Prepare matrices mISO and mFUL for multiplication
 This should be called after setInteractions()
 */
void Meca::prepareMatrices()
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
void Meca::computeForces()
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
This preforms:
 
    fff <- fff + Brownian
    rhs <- time_step * P * fff:
 
 Also prepare Projection diff is requested

 'rhs' and 'fff' are output. Input 'rnd' is a set of independent random numbers
*/
real brownian1(Mecable* mec, real const* rnd, const real alpha, real* fff, real tau, real* rhs)
{
    real n = mec->addBrownianForces(rnd, alpha, fff);
    
    // Calculate the right-hand-side of the system in vRHS:
    mec->projectForces(fff, rhs);
    
    // rhs <- tau * rhs, resulting in time_step * P * fff:
    blas::xscal(DIM*mec->nbPoints(), tau*mec->leftoverMobility(), rhs, 1);

    /*
     At this stage, `fff` contains the external forces in each vertex but also
     internal force such as bending elasticity terms, and the Lagrange multipliers
     do not represent the true tension in the filaments.
     Hence we do not call 'computeTensions(fff)' here
     */

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
 
     ( I - time_step * P * M ) ( Xnew - Xold ) = time_step * P * Force + Noise
 
 where:
 
     Force = M * Xold + B
 
     Noise = std::sqrt(2*kT*time_step*mobility) * Gaussian(0,1)
 
 Implicit integration ensures numerical stability. In the code,
 the sparse matrix M is decomposed as:
 
     M = mISO + mFUL + Rigidity_of_Mecables
 
 Where mISO is isotropic: it applies similarly in the X, Y and Z subspaces, and
 has no crossterms between X and Y or X and Z or Y ans Z subspaces.
 All crossterms go into mFUL.
 The terms should already be set:

     'mISO', 'mFUL' and 'B=vBAS' are set in Meca::setAllInteractions()
     'vPTS = Xold' is set from Mecables' points in Meca::prepare()
     
 The outline of the calculation is:
 
     'vRND' <- calibrated Gaussian random terms ~N(0,1)
     'vFOR' <- force 'M * Xold + B'
     'vRHS' <- time_step * P * F + vRND (right-hand-side of the system)
     'vSOL' <- solution to the linear system of equations is calculated
     'vSOL' <- 'Xnew = Xold + vSOL'
     'vFOR' <- force 'M * Xnew + B'
 
 */
size_t Meca::solve(SimulProp const* prop, const unsigned precond)
{
    assert_true(ready_==0);
    tau_ = prop->time_step;

    prepareMatrices();
    
    // calculate external forces in vFOR:
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
    const real alpha = prop->kT / tau_;
    
    /*
     Add Brownian contributions and calculate Minimum value of it
      vFOR <- vFOR + mobility_coefficient * vRND
      vRHS <- P * vFOR:
     */
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        real local = INFINITY;
        #pragma omp for
        for ( Mecable * mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            real n = brownian1(mec, vRND+inx, alpha, vFOR+inx, tau_, vRHS+inx);
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
    if ( prop->flow.norm() > REAL_EPSILON )
    {
        LOG_ONCE("NEW_CYTOPLASMIC_FLOW code enabled\n");
        Vector flow_dt = prop->flow * tau_;
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
    return 1;
#endif

    // compute preconditionner:
    auto start = __rdtsc();
    computePreconditionner(precond, prop->precondition_span);
    auto factor = __rdtsc() - start;
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
    
    if ( noiseLevel > 0 )
        tolerance_ = noiseLevel * prop->tolerance;
    else
    {
        if ( alpha > 0 )
            Cytosim::log << "Warning: all Brownian terms are zero\n";
        // when temperature == 0, use tolerance as an absolute quantity:
        tolerance_ = prop->tolerance;
    }
    
    /*
     With exact arithmetic, biConjugate Gradient should converge at most
     in a number of iterations equal to the size of the linear system,
     with each BCGGS iteration involving 2 matrix-vector multiplications.
     We set here this theoretical limit to the number multiplications:
     */
    size_t max_iter = std::min(2048UL, 2*dimension());
    LinearSolvers::Monitor monitor(max_iter, tolerance_);

    //fprintf(stderr, "System size %6lu  limit %6lu  tolerance %f precondition %i\n", dimension(), max_iter, tolerance_, precond);

    /*
     GMRES may converge faster than BCGGS, but has overheads and uses more memory
     hence for very large systems, BCGGS is often advantageous.
     */

    //------- call the iterative solver:
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
    
    //fprintf(stderr, "    BCGS%u    count %4lu  residual %.3e\n", precond, monitor.count(), monitor.residual());

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
    fprintf(stderr, "    bcgs     count %4lu  residual %.3e\n", monitor.count(), monitor.residual());
#endif
    
    //------- in case the solver did not converge, we try other methods:
    
    if ( !monitor.converged() )
    {
        Cytosim::out("  no convergence: size %lu precond %i flag %i count %4lu residual %.3e",
            dimension(), precond, monitor.flag(), monitor.count(), monitor.residual());

        monitor.reset();
#if !SAFER_CONVERGENCE_FAILURE
        // try with a different seed
        copy_real(dimension(), vRHS, vSOL);
#endif
        if ( precond )
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
        else
            LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
        
        Cytosim::out(" -> restarted: count %4lu residual %.3e\n", monitor.count(), monitor.residual());

        if ( !monitor.converged() )
        {
            monitor.reset();
            zero_real(dimension(), vSOL);
            if ( precond )
            {
                // try without preconditioner
                LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
                //LinearSolvers::GMRES(*this, vRHS, vSOL, 255, monitor, allocator, mH, mV, temporary);
                Cytosim::out("    BCGS  : count %4lu residual %.3e\n", monitor.count(), monitor.residual());
            }
            else
            {
                // try with strongest preconditioner
                computePreconditionner(4, 0);
                LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                Cytosim::out("    BCGSP4: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
            }
        }

        if ( !monitor.converged() )
        {
            monitor.reset();
            // try with different seed and strongest preconditioner
            copy_real(dimension(), vRND, vSOL);
            computePreconditionner(4, 0);
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
            Cytosim::out(" reseeded: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
        }
        
        if ( !monitor.converged() )
        {
            // if the solver did not converge, its result cannot be used!
            if ( monitor.residual() > 4*tolerance_ )
                throw Exception("no convergence after ",monitor.count()," iterations, residual ",monitor.residual());
            return 0;
        }
    }
    
    const auto total = ( __rdtsc() - start ) >> 10;

    //add the solution (the displacement) to update the Mecable's vertices
    blas::add(dimension(), vSOL, vPTS);

    /*
     Re-calculate forces with the new coordinates, excluding bending elasticity.
     In this way the forces returned to the fibers do not sum-up to zero, and
     are appropriate for example to calculate the effect of force on assembly.
     */
    calculateForces(vPTS, vBAS, vFOR);
    
    // add Brownian terms:
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
        mec->addBrownianForces(vRND+inx, alpha, vFOR+inx);
        //fprintf(stderr, "\n  "); VecPrint::print(stderr, DIM*mec->nbPoints(), vFOR+inx, 2, DIM);
    }

    ready_ = 1;

    // report on the matrix type and size, sparsity, and the number of iterations
    if (( 0 < doNotify ) || ( prop->verbose & 1 ))
    {
        --doNotify;
        std::stringstream oss;
        oss << "\tsize " << DIM << "*" << nbVertices() << " kern " << largestMecable();
#if USE_ISO_MATRIX
        oss << " " << mISO.what();
        if ( useFullMatrix )
#endif
        oss << " " << mFUL.what();
        oss << " precond " << precond;
        oss << " count " << std::setw(3) << monitor.count();
        oss << " residual " << std::setw(11) << std::left << monitor.residual();
        if ( prop->verbose & 4 )
        {
            int cnt = std::max(1UL, monitor.count());
            factor = factor >> 10;
            auto solve = cycles_ >> 10;
            oss << "  cycles " << precond << "T " << std::setw(8) << total;
            oss << " F " << std::setw(8) << factor << std::setw(6) << factor/cnt;
            oss << " S " << std::setw(8) << solve << std::setw(6) << solve/cnt;
            oss << " R " << std::setw(6) << ( total - factor - solve ) / cnt;
        }
        Cytosim::out << oss.str() << "\n";
        if ( prop->verbose & 2 )
            std::clog << oss.str() << "\n";
    }
    
    cycles_ = total;
    return monitor.count();
}


// transfer newly calculated point coordinates back to Mecables
void Meca::apply()
{
    if ( ready_ )
    {
        if ( 1 )
        {
            //check validity of the data:
            bool a = isnan(dimension(), vPTS);
            bool b = isnan(dimension(), vFOR);
            //fprintf(stderr, "Meca::solve isnan %i %i\n", a, b);
            if ( a | b )
            {
                fprintf(stderr, "Meca::solve failed (not-a-number %i %i):\n", a, b);
                for ( Mecable * mec : mecables )
                {
                    b = isnan(DIM*mec->nbPoints(), vPTS+DIM*mec->matIndex());
                    fprintf(stderr, "Mecable %s isnan %i\n", mec->reference().c_str(), b);
                }
                abort();
            }
        }
        #pragma omp parallel for num_threads(NUM_THREADS)
        for ( Mecable * mec : mecables )
        {
            size_t off = DIM * mec->matIndex();
            mec->getPoints(vPTS+off);
            mec->getForces(vFOR+off);
        }
    }
    else
    {
        // if !ready_, the result is not usable
        //printf("superfluous call to Meca::apply()\n");
    }
}


//------------------------------------------------------------------------------
#pragma mark - Analysis Functions

/** Assuming that Mecable::flag() have been set already */
void Meca::flagClusters() const
{
    const size_t MAX = dimension() / DIM;
    Mecable ** which = new Mecable*[MAX];
    
    //for ( size_t i = 0; i < MAX; ++i )
    //    which[i] = nullptr;

    for ( Mecable * mec : mecables )
    {
        //mec->matchFlagIdentity();
        const size_t inx = mec->matIndex();
        const size_t end = mec->nbPoints() + inx;
        assert_true( end <= MAX );
        for ( size_t i = inx; i < end; ++i )
            which[i] = mec;
    }
    
#if USE_MATRIX_BLOCK
    /// equalize flags for any non-zero matrix element between Mecables:
    for ( size_t jj = 0; jj < MAX; ++jj )
    {
        Mecable const* A = which[jj];
        auto const& col = mFUL.column(DIM*jj);
        for ( size_t n = 0; n < col.size_; ++n )
        {
            size_t ii = col.inx_[n]/DIM;
            assert_true( ii < MAX );
            Mecable const* B = which[ii];
            if ( A->flag() != B->flag() )
            {
                ObjectFlag a = std::min(A->flag(), B->flag());
                ObjectFlag b = std::max(A->flag(), B->flag());
                // replace b -> a everywhere:
                for ( Mecable * mec : mecables )
                {
                    if ( mec->flag() == b )
                        mec->flag(a);
                }
            }
        }
    }
#endif
#if USE_ISO_MATRIX
    /// equalize flags for any non-zero matrix element between Mecables:
    for ( size_t jj = 0; jj < MAX; ++jj )
    {
        Mecable const* A = which[jj];
        auto const* col = mISO.column(jj);
        for ( size_t n = 0; n < mISO.column_size(jj); ++n )
        {
            size_t ii = col[n].inx;
            assert_true( ii < MAX );
            Mecable const* B = which[ii];
            if ( A->flag() != B->flag() )
            {
                ObjectFlag a = std::min(A->flag(), B->flag());
                ObjectFlag b = std::max(A->flag(), B->flag());
                // replace b -> a everywhere:
                for ( Mecable * mec : mecables )
                {
                    if ( mec->flag() == b )
                        mec->flag(a);
                }
            }
        }
    }
#endif
    delete[] which;
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
            cnt += ( abs_real(dst[i]) > threshold );
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
            if ( abs_real(dst[i]) > threshold )
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


void Meca::saveSystem(const char dirname[]) const
{
    std::string cwd = FilePath::get_cwd();
    FilePath::change_dir(dirname, true);
    FILE * f = fopen("matrix.mtx", "w");
    if ( f && ~ferror(f) )
    {
        saveMatrix(f, 0);
        fclose(f);
    }
    f = fopen("rhs.mtx", "w");
    if ( f && ~ferror(f) )
    {
        saveRHS(f);
        fclose(f);
    }
    fprintf(stderr, "Cytosim saved its matrix in `%s'\n", dirname);
    FilePath::change_dir(cwd);
}


//------------------------------------------------------------------------------
#pragma mark - MATLAB export Functions

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
        
        mFUL.vecMul(src, res);
        mISO.VECMULADDISO(src, res);
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
    
    int ii = 0;
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = mec->nbPoints();
        for ( size_t p = 0; p < nbp; ++p )
            vec[p] = ii;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), nbp, file);
        ++ii;
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
    fprintf(f, "%lu %i\n", dimension(), DIM);
    fclose(f);
    
    f = fopen("stp.txt", "w");
    fprintf(f, "%f %f\n", tau_, tolerance_);
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


void Meca::dump(const char dirname[]) const
{
    std::string cwd = FilePath::get_cwd();
    FilePath::change_dir(dirname, true);
    dump();
    FilePath::change_dir(cwd);
    fprintf(stderr, "Cytosim exported a %iD system of size %lu in `%s'\n", DIM, dimension(), dirname);
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
    mISO.printSparse(os, 0);
    os.close();
#endif
    
    os.open("d_matC.txt");
    mFUL.printSparse(os, 0);
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
        extractBlock(mec, tmp2);
        VecPrint::sparse_off(os, bks, bks, tmp2, bks, DIM*mec->matIndex());
    }
    os.close();
    
    free_real(tmp1);
    free_real(tmp2);
}

