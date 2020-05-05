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
#include "assert_macro.h"
#include "blas.h"
#include "lapack.h"
#include "cytoblas.h"
#include "xtbsv.h"
#include "xtrsm.h"
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
#define ADD_PROJECTION_DIFF 0


/**
 Set to 1 if the Matrix Vector-multiplication can be distributed.
 This will work only if the matrix mC is specifically build for that purpose
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
: mecables(32)
{
    ready_ = 0;
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
    useMatrixC = false;
#endif
    doNotify = 0;
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
        allocate_vector(alc, vFOR, 1);
        allocate_vector(alc, vTMP, 0);
        
        //std::clog << "Meca::allocate(" << allocated_ << ") " << vFOR << '\n';
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
    free_vector(vFOR);
    free_vector(vTMP);
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
#  if ( DIM == 3 )
    alsatian_xtbsvLN_3D(nbp, mec->block(), BAND_LDD, Y);
    alsatian_xtbsvLT_3D(nbp, mec->block(), BAND_LDD, Y);
#  elif ( DIM == 2 )
    alsatian_xtbsvLN_2D(nbp, mec->block(), BAND_LDD, Y);
    alsatian_xtbsvLT_2D(nbp, mec->block(), BAND_LDD, Y);
#  elif ( DIM == 1 )
    alsatian_xtbsvLN(nbp, mec->block(), BAND_LDD, Y);
    alsatian_xtbsvLT(nbp, mec->block(), BAND_LDD, Y);
#  endif
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
#  if ( DIM == 3 )
    alsatian_xtrsmLLNN_3D(nbp, mec->block(), nbp, Y);
    alsatian_xtrsmLLTN_3D(nbp, mec->block(), nbp, Y);
#  elif ( DIM == 2 )
    alsatian_xtrsmLLNN_2D(nbp, mec->block(), nbp, Y);
    alsatian_xtrsmLLTN_2D(nbp, mec->block(), nbp, Y);
#  elif ( DIM == 1 )
    alsatian_xtrsmLLNN_1D(nbp, mec->block(), nbp, Y);
    alsatian_xtrsmLLTN_1D(nbp, mec->block(), nbp, Y);
#  endif
#else
    //iso_xpotrsL<DIM>(nbp, mec->block(), nbp, Y);
    iso_xtrsmLLN<DIM>('N', nbp, mec->block(), nbp, Y);
    iso_xtrsmLLT<DIM>('N', nbp, mec->block(), nbp, Y);
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
    iso_xgetrsL<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
}


/// apply preconditionner block in full storage
inline void applyPrecondFull(Mecable const* mec, real* Y)
{
    int bks = mec->blockSize();
    int info = 0;
    lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y, bks, &info);
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
        default: throw InvalidParameter("unknown precondition block type"); break;
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
        case 7: mec->blockMultiply(X, Y); break;
        default: throw InvalidParameter("unknown precondition block type"); break;
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
    #pragma omp parallel for num_threads(NUM_THREADS)
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
void getBand(const Mecable * mec, real* res, real time_step)
{
    const size_t nbp = mec->nbPoints();

    if ( mec->hasRigidity() )
    {
        real beta = -time_step * mec->leftoverMobility() * mec->fiberRigidity();
        setRigidityBanded(res, BAND_LDD, nbp, beta);
    }
    else
        zero_real(BAND_LDD*nbp, res);

    // add ones to the diagonal:
    for ( size_t i = 0; i < nbp; ++i )
        res[i*BAND_LDD] += 1.0;

    //std::clog<<"banded preconditionner " << nbp << "\n";
    //VecPrint::print(std::clog, 3, std::min(nbp, 16UL), res, BAND_LDD, 1);
}


/**
 Get the total diagonal block corresponding to an Object, which is:
 
     I - time_step * P ( mB + mC + P' )
 
 The result is constructed by using functions from mB and mC
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
    mB.addTriangularBlock(res, nbp, mec->matIndex(), nbp, 1);
    if ( useMatrixC )
#endif
        mC.addDiagonalTrace(1.0/DIM, res, nbp, DIM*mec->matIndex(), DIM*nbp);
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
    real beta = -time_step * mec->leftoverMobility();
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
 
     I - time_step * P ( mB + mC + P' )
 
 The result is constructed by using functions from mB and mC
 This block is square but not symmetric!
 */
void Meca::getBlock(const Mecable * mec, real* res) const
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
    mB.addTriangularBlock(res, bks, mec->matIndex(), nbp, DIM);
#endif
    expand_lower_matrix<DIM>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useMatrixC )
#endif
        mC.addDiagonalBlock(res, bks, DIM*mec->matIndex(), bks);
    
#if ( 0 )
    std::clog<<"mB+mC block:\n";
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
    real beta = -time_step * mec->leftoverMobility();
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
    std::clog << " | B - K | = " << std::setprecision(3) << err << std::endl;
    
    if ( err > bks * bks * REAL_EPSILON )
    {
        VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        
        size_t s = std::min(16UL, bks);
        extractBlock(mec, wrk);
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
        getBlock(mec, wrk);
        
        // chose initial vector for power iteration
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        
        // power iterations scale like (bs^2)
        real eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, tmp);
 
        //std::clog << "mec " << std::setw(4) << mec->identity();
        if ( eig < 1.0 )
        {
            //std::clog << "    keep     eigen(1-PM) = " << eig << std::endl;
            mec->blockType(4);
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
        getBlock(mec, blk);
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
Compute banded preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondBand(Mecable* mec)
{
    int info = 0;
    const size_t nbp = mec->nbPoints();
    mec->blockSize(DIM*nbp, BAND_LDD*nbp, 0);
    getBand(mec, mec->block(), time_step);
#if CHOUCROUTE
    alsatian_xpbtf2L(nbp, 2, mec->block(), BAND_LDD, &info);
#else
    // cholesky decomposition involves calculating SQRT(diagonal terms)
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
        std::clog << "failed to compute Preconditionner block\n";
    }
}


/**
Compute preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondIsoS(Mecable* mec)
{
    const size_t nbp = mec->nbPoints();
#if 0
    const size_t S = std::min(nbp, 16UL);
    const size_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, bks);
    getBlock(mec, mec->block());
    project_matrix<DIM>(nbp, mec->block(), bks, mec->block(), nbp);
    std::clog<<"projected: " << bks << "\n";
    VecPrint::print(std::clog, S, S, mec->block(), bks, 2);
#endif
    int info = 0;

    mec->blockSize(DIM*nbp, nbp*nbp, 0);
    
    getIsoBlock(mec, mec->block());
    
    //std::clog<<"iso symmetric preconditionner: " << nbp << "\n";
    //VecPrint::print(std::clog, S, S, mec->block(), nbp, 2);

// calculate LU factorization:
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
        std::clog << "failed to compute Preconditionner block\n";
    }
}

/**
Compute preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondIsoP(Mecable* mec)
{
    const size_t nbp = mec->nbPoints();
    int info = 0;

    const size_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, nbp);
    getBlock(mec, mec->block());
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
        std::clog << "failed to compute Preconditionner block\n";
    }
}

/**
Compute preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondFull(Mecable* mec)
{
    const size_t bks = DIM * mec->nbPoints();
    mec->blockSize(bks, bks*bks, bks);
    
    getBlock(mec, mec->block());
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
        std::clog << "failed to compute Preconditionner block\n";
    }
}


/**
 Initialize Mecable::blockMatrix() as the preconditionner block,
 using the factorization given by 'blk' and 'piv'
 */
void convertPreconditionner(Mecable* mec, real* blk, int* piv, real* wrk)
{
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
}


/**
renew preconditionner block corresponding to 'mec'
 */
void Meca::renewPreconditionner(Mecable* mec, int span, real* blk, int* piv, real* wrk, size_t wrksize)
{
    const size_t bks = DIM * mec->nbPoints();
    const int age = mec->blockAge();

    if (( mec->blockType() != 7 ) | ( mec->blockSize() != bks ) | ( age >= span ))
    {
        // recalculate block!
        getBlock(mec, blk);
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
    accum_ += __rdtsc() - rdt;
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/// function to sort Mecables
int compareMecables(const void * ap, const void * bp)
{
    size_t a = (*((Mecable const**)ap))->nbPoints();
    size_t b = (*((Mecable const**)bp))->nbPoints();
    return ( a < b ) - ( a > b );
}


/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::prepare(Simul* sim)
{
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
        mec->matIndex(cnt);
        cnt += mec->nbPoints();
    }
    nPoints_ = cnt;
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
    
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable * mec : mecables )
    {
        mec->putPoints(vPTS+DIM*mec->matIndex());
        mec->prepareMecable();
#if ( DIM > 1 ) && !SEPARATE_RIGIDITY_TERMS
        if ( mec->hasRigidity() )
        {
#   if USE_ISO_MATRIX
            addRigidityMatrix<DIM>(mB, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   else
            addRigidityBlockMatrix<DIM>(mC, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   endif
        }
#endif
    }
    //fprintf(stderr, "Meca::prepare() isnan %i\n", isnan(dimension(), vPTS));
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
size_t Meca::solve(SimulProp const* prop, const unsigned precond)
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
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        real local = INFINITY;
        #pragma omp for
        for ( Mecable * mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            real n = brownian1(mec, vRND+inx, alpha, vFOR+inx, time_step, vRHS+inx);
            local = std::min(local, n);
            //printf("thread %i min: %f\n", omp_get_thread_num(), local);
        }
        #pragma omp critical
        noiseLevel = std::min(noiseLevel, local);
    }

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
    return 1;
#endif

    // compute preconditionner:
    auto start_ = __rdtsc();
    computePreconditionner(precond, prop->precondition_span);
    auto precond_ = __rdtsc() - start_;
    start_ = __rdtsc();
    accum_ = 0;

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
     We set here this theoretical limit to the number multiplications:
     */

    LinearSolvers::Monitor monitor(2*dimension(), abstol);

    //fprintf(stderr, "Solve precond %i size %6i absolute tolerance %f\n", precond, dimension(), abstol);
    
    /*
     GMRES may converge faster than BCGGS, but has overheads and uses more memory
     hence for large systems, biCGGS is often advantageous.
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
        Cytosim::out("  no convergence: size %lu precond %i flag %i count %4lu residual %.3e\n",
            dimension(), precond, monitor.flag(), monitor.count(), monitor.residual());

        monitor.reset();
        // try with different seed
        copy_real(dimension(), vRHS, vSOL);
        if ( precond )
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
        else
            LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
        Cytosim::out("  new seed: count %4lu residual %.3e\n", monitor.count(), monitor.residual());

        if ( !monitor.converged() )
        {
            monitor.reset();
            zero_real(dimension(), vSOL);
            if ( precond )
            {
                // try with different method
                LinearSolvers::GMRES(*this, vRHS, vSOL, 255, monitor, allocator, mH, mV, temporary);
                Cytosim::out("    GMRES: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
            }
            else
            {
                // try with a preconditioner
                computePreconditionner(1, 0);
                LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                Cytosim::out("    BCGSP: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
            }
        }

        if ( !monitor.converged() )
        {
            monitor.reset();
            // try with different seed
            copy_real(dimension(), vRND, vSOL);
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
            Cytosim::out(" another seed: count %4lu residual %.3e\n", monitor.count(), monitor.residual());
        }
        
        if ( !monitor.converged() )
        {
            // no method could converge... this is really bad!
            if ( monitor.residual() > 4*abstol )
                throw Exception("no convergence after ",monitor.count()," iterations, residual ",monitor.residual());
            return 0;
        }
    }
    
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
    if ( 0 < doNotify || prop->verbose )
    {
        --doNotify;
        std::stringstream oss;
        oss << "\tsize " << DIM << "*" << nb_points() << " kern " << largestMecable();
#if USE_ISO_MATRIX
        oss << " " << mB.what();
        if ( useMatrixC )
#endif
        oss << " " << mC.what();
        oss << " precond " << precond;
        oss << " count " << std::setw(3) << monitor.count();
        oss << " residual " << std::setw(11) << std::left << monitor.residual();
        if ( prop->verbose & 4 )
        {
            auto total = ( __rdtsc() - start_ ) >> 10;
            int cnt = std::max(1UL, monitor.count());
            precond_ = precond_ >> 10;
            accum_ = accum_ >> 10;
            oss << "  cycles T " << std::setw(8) << total;
            oss << " F " << std::setw(8) << precond_ << std::setw(6) << precond_/cnt;
            oss << " S " << std::setw(8) << accum_ << std::setw(6) << accum_/cnt;
            oss << " R " << std::setw(6) << ( total - accum_ ) / cnt;
        }
        Cytosim::out << oss.str() << "\n";
        if ( prop->verbose & 2 )
            std::clog << oss.str() << "\n";
    }
    
    return monitor.count();
}


// transfer newly calculated point coordinates back to Mecables
void Meca::apply()
{
    if ( ready_ )
    {
        ready_ = 0;
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
        //#pragma omp parallel for num_threads(NUM_THREADS)
        for ( Mecable * mec : mecables )
        {
            mec->getForces(vFOR+DIM*mec->matIndex());
            mec->getPoints(vPTS+DIM*mec->matIndex());
        }
    }
    else
    {
        printf("superfluous Meca::apply() call\n");
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
        extractBlock(mec, tmp2);
        VecPrint::sparse_off(os, bks, bks, tmp2, bks, DIM*mec->matIndex());
    }
    os.close();
    
    free_real(tmp1);
    free_real(tmp2);
}

