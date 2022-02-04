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
#include "xtbsv.h"
#include "xtrsm.h"

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

#include "meca_inter.cc"
#include "meca_steric.cc"
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
    // Y <- X + beta * Y
    const real beta = -tau_ * mec->leftoverMobility();
    blas::xpay(bks, X+inx, beta, Y+inx);
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
#pragma mark - Apply Preconditionner


// the leading dimension of the banded matrix used for preconditionning
constexpr size_t ISOB_LDD = 3;
/*
 We use a band preconditionner with 2*DIM off-diagonals to include
 the 'diagonal' terms from the blocks that are offset by 2 from the diagonal
 */
constexpr size_t BAND_NUD = 2*DIM;
constexpr size_t BAND_LDD = BAND_NUD+DIM;

/// apply banded symmetric isotropic preconditionner block
inline void applyPrecondIsoB(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();

#if CHOUCROUTE
    alsatian_xpbtrsL<DIM>(nbp, mec->block(), ISOB_LDD, Y);
#else
    /*
     we cannot call lapack::DPBTRS('L', bks, KD, 1, mec->block(), KD+1, Y, bks, &info)
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     But calling DTBSV gets the required work done.
     */
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', nbp, 2, mec->block(), ISOB_LDD, Y+d, DIM);
        blas::xtbsv('L', 'T', 'N', nbp, 2, mec->block(), ISOB_LDD, Y+d, DIM);
    }
#endif
}


/// apply symmetric isotropic preconditionner block
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
    std::clog << "\n "; VecPrint::print(S, tmp, 3, 100.0);
    free_real(tmp);
#endif
#if CHOUCROUTE
    alsatian_xpotrsL<DIM>(nbp, mec->block(), nbp, Y);
#else
    iso_xpotrsL<DIM>(nbp, mec->block(), nbp, Y);
#endif
    //std::clog << "\nL"; VecPrint::print(S, Y, 3, 100.0);
}


/// apply non-symmetric but isotropic preconditionner block
inline void applyPrecondIsoP(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
    //iso_xgetrsN<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
    alsatian_xgetrsN<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
}


/// apply banded symmetric preconditionner block
inline void applyPrecondBand(Mecable const* mec, real* Y)
{
    const int bks = mec->blockSize();
    assert_true( (int)BAND_NUD < bks );
#if CHOUCROUTE
    alsatian_xpbtrsLK<BAND_NUD>(bks, mec->block(), BAND_LDD, Y);
#elif 1
    blas_xtbsvLN<'N'>(bks, BAND_NUD, mec->block(), BAND_LDD, Y, 1);
    blas_xtbsvLT<'N'>(bks, BAND_NUD, mec->block(), BAND_LDD, Y, 1);
#elif 1
    blas::xtbsv('L', 'N', 'N', bks, BAND_NUD, mec->block(), BAND_LDD, Y, 1);
    blas::xtbsv('L', 'T', 'N', bks, BAND_NUD, mec->block(), BAND_LDD, Y, 1);
#else
    int info = 0;
    lapack::xpbtrs('L', bks, BAND_NUD, 1, mec->block(), BAND_LDD, Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply symmetric preconditionner block in full storage
inline void applyPrecondHalf(Mecable const* mec, real* Y)
{
    const int bks = mec->blockSize();
#if CHOUCROUTE
    // assuming that diagonal terms of the preconditionner block have been inverted:
    alsatian_xpotrsL(bks, mec->block(), bks, Y);
#elif 1
    iso_xpotrsL<1>(bks, mec->block(), bks, Y);
#else
    int info = 0;
    lapack::xpotrs('L', bks, 1, mec->block(), bks, Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply non-symmetric preconditionner block in full storage
inline void applyPrecondFull(Mecable const* mec, real* Y)
{
    const int bks = mec->blockSize();
#if CHOUCROUTE && defined(__SSE3__)
    // assuming that diagonal terms of the preconditionner block have been inverted:
    alsatian_xgetrsN_SSE(bks, mec->block(), bks, mec->pivot(), Y);
#elif CHOUCROUTE
    alsatian_xgetrsN(bks, mec->block(), bks, mec->pivot(), Y);
#elif 1
    // translated LAPACK's reference code:
    lapack_xgetrsN(bks, mec->block(), bks, mec->pivot(), Y);
#else
    // using LAPACK's library
    int info = 0;
    lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply preconditionner block corresponding to Mecable
inline void applyPreconditionner(Mecable const* mec, real* Y)
{
    switch ( mec->blockType() )
    {
        case 0: break;
        case 1: applyPrecondIsoB(mec, Y); break;
        case 2: applyPrecondIsoS(mec, Y); break;
        case 3: applyPrecondIsoP(mec, Y); break;
        case 4: applyPrecondBand(mec, Y); break;
        case 5: applyPrecondHalf(mec, Y); break;
        case 6: applyPrecondFull(mec, Y); break;
        default: ABORT_NOW("unknown Mecable::blockType()"); break;
    }
}


/// apply preconditionner to entire system: Y <- Preconditionner * X
/** This can be done in parallel if the preconditionner is block-diagonal */
void Meca::precondition(const real* X, real* Y) const
{
    auto rdt = timer();
    if ( Y != X )
        copy_real(dimension(), X, Y);
    
#pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable const* mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
#if EXPERIMENTAL_PRECONDITIONNERS
        if ( mec->blockType() == 7 )
            mec->blockMultiply(X+inx, Y+inx);
        else
#endif
            applyPreconditionner(mec, Y+inx);
    }
    cycles_ += timer() - rdt;
}


/// total allocated memory size for preconditionner
size_t Meca::preconditionnerSize() const
{
    size_t res = 0;
    for ( Mecable const* mec : mecables )
        res += mec->blockAllocated();
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner

void printBlock(Mecable* mec, size_t sup)
{
    const size_t bks = DIM * mec->nbPoints();
    size_t S = std::min(bks, sup);
    real * blk = mec->block();

    VecPrint::full("Diagonal block"+std::to_string(bks), S, S, mec->block(), bks, 3);
    
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
 
 This block is square, symmetric, definite positive and predictable
 */
void Meca::getIsoBBlock(const Mecable * mec, real* res, size_t ldd) const
{
    const size_t nbp = mec->nbPoints();

    const real beta = -tau_ * mec->pointMobility();
    
    if ( mec->hasRigidity() )
    {
        if ( ldd != 3 )
            zero_real(ldd*nbp, res);
        setBendingRigidity<1>(res, ldd-1, nbp, beta*mec->fiberRigidity());
    }
    else
        zero_real(ldd*nbp, res);
    
    //std::clog << "\nrigidity band " << nbp << "\n";
    //VecPrint::full(3, std::min(nbp, 16UL), res, ldd, 1);

    /*
     The matrix `res` is stored in 'packed symmetric banded storage':
     usually, mat(i, j) is stored in mat[i+ldd*j]
     but with banded storage, mat(i, j) is stored in mat[i-j+ldd*j] for i > j
     So we use below ldd-1
     */

#if USE_ISO_MATRIX
    mISO.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, 2);
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalTrace(beta/DIM, res, ldd-1, DIM*mec->matIndex(), DIM*nbp, DIM*2, false);

    // add Identity matrix to band storage:
    for ( size_t i = 0; i < nbp; ++i )
        res[ldd*i] += 1;
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
        //addBendingRigidityLower<1>(res, nbp, mec->nbPoints(), mec->fiberRigidity());
        addBendingRigidity<1>(res, nbp, mec->nbPoints(), mec->fiberRigidity());
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, nbp, mec->matIndex(), nbp);
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalTrace(1.0/DIM, res, nbp, DIM*mec->matIndex(), DIM*nbp, DIM*nbp, true);

#if ( 0 )
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
    {
        // Include the corrections P' in preconditioner, vector by vector.
        real* tmp = vTMP + DIM * mec->matIndex();
        zero_real(bks, tmp);
        for ( size_t i = 0; i < bks; ++i )
        {
            tmp[i] = 1;
            mec->addProjectionDiff(tmp, res+bks*i);
            tmp[i] = 0;
        }
        //VecPrint::full("dynamic with P'", bs, bs, res, bs);
    }
#endif
    
    //include the projection, by applying it to each column vector:
    for ( size_t i = 0; i < bks; ++i )
        mec->projectForces(res+bks*i, res+bks*i);

    const real beta = -tau_ * mec->leftoverMobility();
#else
    // the projection is not called, so we scale by mobility
    const real beta = -tau_ * mec->pointMobility();
#endif

    //blas::xscal(bs*bs, beta, res, 1);
    for ( size_t n = 0; n < nbp*nbp; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    for ( size_t i = 0; i < nbp*nbp; i += nbp+1 )
        res[i] += 1;
}


/**
 Get a diagonal block corresponding to an Object, which is:
 
     I - time_step * mob * ( mISO + mFUL )
 
 The result is constructed by using functions from mISO and mFUL, and then
 multiplied by the vertex mobility, but projection is not applied.
 This block is banded and symmetric, and can be factorized by Cholesky's method!
 */
void Meca::getBandedBlock(const Mecable * mec, real* res, size_t ldd, size_t rank) const
{
    const size_t nbp = mec->nbPoints();
    const size_t bks = DIM * nbp;

    zero_real(ldd*bks, res);

    // multiply by mobility coefficient, skipping projection
    const real beta = -tau_ * mec->pointMobility();
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    if ( mec->hasRigidity() )
    {
        setBendingRigidity<DIM>(res, ldd-1, nbp, beta*mec->fiberRigidity());
        //std::clog<<"Rigidity block " << mec->reference() << "\n";
        //VecPrint::full(bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, rank/DIM);
#endif
    // VecPrint::full("\niso", ldd, std::min(bks, 24ul), res, ldd);
    copy_lower_subspace<DIM>(bks, res, ldd-1, rank);
    // VecPrint::full("\ncopy_subspace", ldd, std::min(bks, 24ul), res, ldd);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addLowerBand(beta, res, ldd-1, DIM*mec->matIndex(), bks, rank);
    
    // add Identity matrix to band storage:
    for ( size_t i = 0; i < bks; ++i )
        res[ldd*i] += 1;
}

/**
 Get a diagonal block corresponding to an Object, which is:
 
     I - time_step * mob * ( mISO + mFUL )
 
 The result is constructed by using functions from mISO and mFUL, and then
 multiplied by the vertex mobility, to approximate the dynamics.
 This block is square and symmetric, and can be factorized by Cholesky's method!
 */
void Meca::getHalfBlock(const Mecable * mec, real* res) const
{
    const size_t nbp = mec->nbPoints();
    const size_t bks = DIM * nbp;
    
    zero_real(bks*bks, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    if ( mec->hasRigidity() )
    {
        addBendingRigidityLower<DIM>(res, bks, mec->nbPoints(), mec->fiberRigidity());
        //VecPrint::full("Rigidity block", bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
#endif
    copy_lower_subspace<DIM, true>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, DIM*mec->matIndex(), bks);
    
    // multiply by mobility coefficient, skipping projection
    const real beta = -tau_ * mec->pointMobility();
    
    //blas::xscal(bs*bs, beta, res, 1);
    for ( size_t n = 0; n < bks*bks; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    for ( size_t i = 0; i < bks*bks; i += bks+1 )
        res[i] += 1;
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
        addBendingRigidityLower<DIM>(res, bks, mec->nbPoints(), mec->fiberRigidity());
        //VecPrint::full("Rigidity block", bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
#endif
    copy_lower_subspace<DIM, true>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, DIM*mec->matIndex(), bks);
    
    //VecPrint::full("mISO+mFUL block", bks, bks, res, bks);
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
    {
        // Include the corrections P' in preconditioner, vector by vector.
        real* tmp = vTMP + DIM * mec->matIndex();
        zero_real(bks, tmp);
        for ( size_t i = 0; i < bks; ++i )
        {
            tmp[i] = 1;
            mec->addProjectionDiff(tmp, res+bks*i);
            tmp[i] = 0;
        }
        //VecPrint::full("dynamic with P'", bs, bs, res, bs);
    }
#endif
    
    // include the projection, by applying it to each column vector:
    /* This could be vectorized */
    for ( size_t i = 0; i < bks; ++i )
        mec->projectForces(res+bks*i, res+bks*i);
    
    // scale
    const real beta = -tau_ * mec->leftoverMobility();
    //blas::xscal(bs*bs, beta, res, 1);
    for ( size_t n = 0; n < bks*bks; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    for ( size_t i = 0; i < bks*bks; i += bks+1 )
        res[i] += 1;
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
        vec[jj+off] = 1;
        multiply(vec, tmp);
        vec[jj+off] = 0;
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
    std::clog << " | B - K | = " << err << '\n';

    if ( err > bks * bks * REAL_EPSILON )
    {
        //VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        extractBlock(mec, wrk);
        
        const size_t S = std::min(9UL, bks);
        VecPrint::full("extracted", S, S, wrk, bks, 3);
        VecPrint::full("computed", S, S, blk, bks, 3);
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
        applyPreconditionner(mec, mat+bks*i);
        mat[i+bks*i] -= 1;
    }
    real err = blas::nrm2(bks*bks, mat) / bks;
    std::clog << " | 1 - PM | = " << err;
    
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
        VecPrint::full("matrix", S, S, wrk, bks);
        VecPrint::full("precond", S, S, blk, bks);
        VecPrint::full("precond * matrix", S, S, mat, bks);
        std::cout << "\n";
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
            mec->blockType(6);
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
    
    mec->blockType(6);
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
 This block is usually symmetric definite positive, and is factorized by Cholesky's method
 */
void Meca::computePrecondIsoB(Mecable* mec)
{
    assert_true(ISOB_LDD>2);
    const size_t nbp = mec->nbPoints();
    mec->blockSize(DIM*nbp, ISOB_LDD*nbp, 0);
    
    //std::clog << "banded preconditionner " << nbp << "\n";
    //VecPrint::full(3, std::min(nbp, 16UL), mec->block(), ISOB_LDD, 1);
    
    /**
     Factorize banded matrix with Andre-Louis Cholesky's method
     born 15.10.1875 in Montguyon, France
     died 31.08.1918 in Bagneux, following wounds received in battle.
     */
    int bt, info = 0;
    if ( 1 ) //ISOB_LDD <= nbp )
    {
        getIsoBBlock(mec, mec->block(), ISOB_LDD);
        // calculate Banded Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpbtf2L<2>(nbp, mec->block(), ISOB_LDD, &info);
#else
        lapack::xpbtf2('L', nbp, 2, mec->block(), ISOB_LDD, &info);
#endif
        bt = 1;
    }
    else
    {
        getIsoBlock(mec, mec->block());
        // calculate Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpotf2L(nbp, mec->block(), nbp, &info);
#else
        lapack::xpotf2('L', nbp, mec->block(), nbp, &info);
#endif
        bt = 2;
    }

    if ( 0 == info )
    {
        mec->blockType(bt);
        //std::clog<<"factorized banded preconditionner: " << nbp << "\n";
        //VecPrint::full(3, std::min(nbp, 16UL), mec->block(), ISOB_LDD, 1);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute Band Preconditionner block of size " << nbp << "\n";
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
    VecPrint::full(S, S, mec->block(), nbp, 2);
#endif

    mec->blockSize(DIM*nbp, nbp*nbp, 0);
    
    //getIsoBBlock(mec, mec->block(), ISOB_LDD);
    //std::clog << "banded preconditionner " << nbp << "\n";
    //VecPrint::full(3, S, mec->block(), ISOB_LDD, 1);

    getIsoBlock(mec, mec->block());
    
    //std::clog<<"iso symmetric preconditionner: " << nbp << "\n";
    //VecPrint::full(S, S, mec->block(), nbp, 2);

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
        //std::clog << "failed to compute IsoS Preconditionner block of size " << nbp << "\n";
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
    //VecPrint::full(S, S, mec->block(), nbp, 2);

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
        //std::clog << "failed to compute IsoP Preconditionner block of size " << nbp << "\n";
        ++bump_;
    }
}

/**
 Compute banded symmetric preconditionner block corresponding to 'mec',
 factorized by Cholesky's method.
 */
void Meca::computePrecondBand(Mecable* mec)
{
    assert_true(BAND_NUD < BAND_LDD);
    const size_t bks = DIM * mec->nbPoints();
   
    int bt, info = 0;

    if ( BAND_LDD < bks )
    {
        mec->blockSize(bks, BAND_LDD*bks, 0);
        getBandedBlock(mec, mec->block(), BAND_LDD, BAND_NUD);
        
#if ( 0 )
        size_t S = std::min(bks, 24UL);
        VecPrint::full("band block "+std::to_string(bks), BAND_LDD, S, mec->block(), BAND_LDD);
        mec->blockSize(bks, bks*bks, 0);
        getHalfBlock(mec, mec->block());
        VecPrint::full("half block", S, S, mec->block(), bks);
#endif
        
        // calculate Cholesky factorization for band storage:
#if CHOUCROUTE
        alsatian_xpbtf2L<BAND_NUD>(bks, mec->block(), BAND_LDD, &info);
#else
        lapack::xpbtf2('L', bks, BAND_NUD, mec->block(), BAND_LDD, &info);
#endif
        bt = 4;
    }
    else
    {
        mec->blockSize(bks, bks*bks, 0);
        getHalfBlock(mec, mec->block());
        // calculate Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpotf2L(bks, mec->block(), bks, &info);
#else
        lapack::xpotf2('L', bks, mec->block(), bks, &info);
#endif
        bt = 5;
    }
    
    if ( 0 == info )
    {
        mec->blockType(bt);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute half Preconditionner block of size " << bks << "\n";
        ++bump_;
    }
}

/**
 Compute preconditionner block corresponding to 'mec'
 This block is symmetric, and factorized by Cholesky's method.
 */
void Meca::computePrecondHalf(Mecable* mec)
{
    const size_t bks = DIM * mec->nbPoints();
    mec->blockSize(bks, bks*bks, 0);
    getHalfBlock(mec, mec->block());
    
#if ( 0 )
    size_t S = std::min(bks, 16UL);
    VecPrint::full("half block "+std::to_string(bks), S, S, mec->block(), bks);
#endif

    int info = 0;
    // calculate Cholesky factorization:
#if CHOUCROUTE
    alsatian_xpotf2L(bks, mec->block(), bks, &info);
#else
    lapack::xpotf2('L', bks, mec->block(), bks, &info);
#endif

    if ( 0 == info )
    {
        mec->blockType(5);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute half Preconditionner bloc of size " << bks << "\n";
        ++bump_;
    }
}

/**
Compute preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondFull(Mecable* mec)
{
    const size_t bks = DIM * mec->nbPoints();
    
#if CHOUCROUTE && REAL_IS_DOUBLE
    mec->blockSize(bks, 4+bks*bks/2, bks);
    // use temporary memory to build matrix block:
    double* blk = new_real(4+bks*bks);
#else
    mec->blockSize(bks, bks*bks, bks);
    real * blk = mec->block();
#endif
    
    getFullBlock(mec, blk);
    //verifyBlock(mec, blk);
    
    // calculate LU factorization:
    int info = 0;
#if CHOUCROUTE
    alsatian_xgetf2(bks, blk, bks, mec->pivot(), &info);
#else
    lapack::xgetf2(bks, bks, blk, bks, mec->pivot(), &info);
#endif
    
    if ( 0 == info )
    {
        mec->blockType(6);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute full Preconditionner block of size " << bks << "\n";
        ++bump_;
    }
    
#if CHOUCROUTE && REAL_IS_DOUBLE
    convert_to_floats(bks*bks, blk, (float*)mec->block());
    free_real(blk);
#endif
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


/*
 Here the preconditionner blocks are calculated,
 according to Simul::precondition
 */
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
                computePrecondIsoB(mec);
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
                computePrecondBand(mec);
            break;
        case 5:
            for ( Mecable * mec : mecables )
                computePrecondHalf(mec);
            break;
        case 6:
            for ( Mecable * mec : mecables )
                computePrecondFull(mec);
            break;
        case 7:
            renewPreconditionner(span);
            break;
        default:
            throw InvalidParameter("unknown `precondition' value");
            break;
    }
    if ( bump_ > 0 )
        Cytosim::log << "failed to compute " << bump_ << " / " << mecables.size() << " preconditionner blocks\n";
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/// qsort function comparing number of points of Mecables
static int compareMecables(const void * A, const void * B)
{
    size_t a = (*static_cast<Mecable *const*>(A))->nbPoints();
    size_t b = (*static_cast<Mecable *const*>(B))->nbPoints();
    return ( a < b ) - ( a > b );
}


void Meca::pickMecables(Simul const& sim)
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
}

/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::getReady()
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
            addBendingRigidityMatrix(mISO, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   else
            addBendingRigidityBlockMatrix<DIM>(mFUL, mec->matIndex(), mec->nbPoints(), mec->fiberRigidity());
#   endif
        }
#endif
    }
    //fprintf(stderr, "Meca::prepare() isnan %i\n", has_nan(dimension(), vPTS));
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
size_t Meca::solve(SimulProp const* prop, const unsigned precond)
{
    assert_true(ready_==0);
    tau_ = prop->time_step;

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
    alpha_ = prop->kT / tau_;
    
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
        mec->getForces(vRHS+DIM*mec->matIndex());
    }
    return 1;
#endif

    // compute preconditionner:
    auto start = timer();
    computePreconditionner(precond, prop->precondition_span);
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
    
    if ( noiseLevel > 0 )
        tolerance_ = noiseLevel * prop->tolerance;
    else
    {
        if ( alpha_ > 0 )
            Cytosim::log << "Warning: all Brownian terms are zero\n";
        // when temperature == 0, use tolerance as an absolute quantity:
        tolerance_ = prop->tolerance;
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
    if ( precond )
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
            dimension(), precond, monitor.flag(), monitor.count(), monitor.residual(), monitor.residual()/tolerance_);
        
        // in case the solver did not converge, we try other methods:
        monitor.reset();
#if !SAFER_CONVERGENCE
        // try with a different seed
        precondition(vRHS, vSOL);
#endif
        if ( precond )
            LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator_);
        else
            LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
        
        Cytosim::out(" --> restarted: count %4i residual %.3e\n", monitor.count(), monitor.residual());

        // relax the convergence criteria a bit
        if ( monitor.residual() > 1.4142 * tolerance_ )
        {
            // try with our strongest preconditioner
            computePreconditionner(6, 0);
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
    if (( 0 < doNotify ) || ( prop->verbose & 1 ))
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
        oss << " precond " << precond << " (" << preconditionnerSize() << ")";
        oss << " count " << std::setw(4) << monitor.count();
        oss << " residual " << std::setw(11) << std::left << monitor.residual();
        size_t dim = dimension();
        if ( prop->verbose & 8 )
        {
            // calculate true residual: tmp = rhs - A * x
            real * tmp = allocator_.bind(0);
            multiply(vSOL, tmp);
            blas::xsub(dim, vRHS, tmp);
            oss << ": " << std::setw(11) << std::left << blas::nrm8(dim, tmp);
            oss << " dx " << std::setw(11) << std::left << blas::nrm8(dim, vSOL);
        }
        if ( prop->verbose & 4 )
        {
            unsigned cnt = std::max(1U, monitor.count());
            oss << "  cycles " << precond << "T " << std::setw(8) << cycles_;
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
        
        // add Brownian terms:
        for ( Mecable * mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            mec->addBrownianForces(vRND+inx, alpha_, vFOR+inx);
            //fprintf(stderr, "\n  "); VecPrint::print(stderr, DIM*mec->nbPoints(), vFOR+inx, 2, DIM);
        }
        
        if ( 1 )
        {
            //check validity of the data:
            bool a = has_nan(dimension(), vPTS);
            bool b = has_nan(dimension(), vFOR);
            //fprintf(stderr, "Meca::solve isnan %i %i\n", a, b);
            if ( a | b )
            {
                fprintf(stderr, "Meca::solve failed (not-a-number %i %i):\n", a, b);
                for ( Mecable * mec : mecables )
                {
                    b = has_nan(DIM*mec->nbPoints(), vPTS+DIM*mec->matIndex());
                    fprintf(stderr, "Mecable %s isnan %i\n", mec->reference().c_str(), b);
                }
                abort();
            }
        }

        #pragma omp parallel for num_threads(NUM_THREADS)
        for ( Mecable * mec : mecables )
        {
            size_t off = DIM * mec->matIndex();
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


//------------------------------------------------------------------------------
#pragma mark - Connectivity Analysis


/// equalize flags for any existing matrix element between Mecables
template < typename MatrixClass >
void computeClusters(Array<Mecable*> const mecables, const size_t MAX, Mecable** table,
                     MatrixClass const& MAT, size_t ORD)
{
    for ( size_t j = 0; j < MAX; ++j )
    {
        size_t jj = ORD * j;
        Mecable const* A = table[j];
        for ( size_t n = 0; n < MAT.column_size(jj); ++n )
        {
            // we do not check the value here, but just having a block
            size_t i = MAT.column_index(jj, n) / ORD;
            assert_true( i < MAX );
            Mecable const* B = table[i];
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
}


/** Assuming that Mecable::flag() have been set already */
void Meca::flagClusters() const
{
    const size_t MAX = nbVertices();
    Mecable ** table = new Mecable*[MAX]{nullptr};
    
    for ( Mecable * mec : mecables )
    {
        //mec->matchFlagIdentity();
        const size_t inx = mec->matIndex();
        const size_t end = mec->nbPoints() + inx;
        assert_true( end <= MAX );
        for ( size_t i = inx; i < end; ++i )
            table[i] = mec;
    }
    
#if USE_MATRIX_BLOCK
    computeClusters(mecables, MAX, table, mFUL, DIM);
#endif
#if USE_ISO_MATRIX
    computeClusters(mecables, MAX, table, mISO, 1);
#endif
    delete[] table;
}


//------------------------------------------------------------------------------
#pragma mark - Matrix Extraction

/**
 Count number of non-zero entries in the full system matrix
 */
size_t Meca::countTerms(const real threshold) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * dst = new_real(dim);
    zero_real(dim, src);
    
    size_t cnt = 0;
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        multiply(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            cnt += ( abs_real(dst[i]) >= threshold );
        src[j] = 0;
    }
    
    free_real(dst);
    free_real(src);
    return cnt;
}

/**
 Extract the full matrix associated with Meca::multiply().
 The array `mat[]` should be preallocated to hold `dim*lda` real scalars,
 with `dim >= Meca::dimension()`, and `lda >= dim` the leading dimension of
 the array.
 */
void Meca::getMatrix(real * mat, size_t lda) const
{
    size_t dim = dimension();
    if ( lda < dim )
        throw InvalidIO("invalid matrix dimensions");
    real * src = new_real(dim);
    zero_real(dim, src);
    
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        multiply(src, mat+j*lda);
        src[j] = 0;
    }
    
    free_real(src);
}

//------------------------------------------------------------------------------
#pragma mark - Text Export

static void saveVector(FILE * fp, size_t dim, real const* VEC)
{
    fprintf(fp, "%% This is a vector produced by Cytosim\n");
    fprintf(fp, "%% author: Francois J. Nedelec\n");
    fprintf(fp, "%% kind: biological cell simulation (cytoskeleton)\n");
    
    fprintf(fp, "%lu\n", dim);
    for ( size_t i = 0; i < dim; ++i )
        fprintf(fp, "%f\n", VEC[i]);
}


void Meca::saveObjectID(FILE * fp) const
{
    int i = 1;
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = DIM * mec->nbPoints();
        for ( size_t p = 0; p < nbp; ++p )
            fprintf(fp, "%if\n", i);
        ++i;
    }
}

void Meca::saveMobility(FILE * fp) const
{
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = mec->nbPoints();
        const real val = mec->pointMobility();
        for ( size_t p = 0; p < DIM * nbp; ++p )
            fprintf(fp, "%f\n", val);
    }
}

/**
 Save a sparse matrix in Matrix Market format
 https://math.nist.gov/MatrixMarket/formats.html
 This is a Sparse text format
 */
void Meca::saveMatrix(FILE * fp, real threshold) const
{
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%% This is a matrix produced by Cytosim\n");
    fprintf(fp, "%% author: Francois J. Nedelec\n");
    fprintf(fp, "%% kind: biological cell simulation (cytoskeleton)\n");

    const size_t dim = dimension();
    real * src = new_real(dim);
    real * dst = new_real(dim);
    zero_real(dim, src);
    
    fprintf(fp, "%lu %lu ", dim, dim);

    fpos_t pos;
    fgetpos(fp, &pos);
    size_t cnt = 0;
    fprintf(fp, "%10lu\n", cnt);

    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        multiply(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            if ( abs_real(dst[i]) > threshold )
            {
                fprintf(fp, "%3lu %3lu %f\n", i, j, dst[i]);
                ++cnt;
            }
        src[j] = 0;
    }
    
    fsetpos(fp, &pos);
    fprintf(fp, "%10lu\n", cnt);

    free_real(dst);
    free_real(src);
}


/**
 Save Matrix and Right-hand-side Vector
 */
void Meca::saveSystem() const
{
    FILE * f = FilePath::open_file("matrix.mtx", "w");
    saveMatrix(f, 0);
    fclose(f);
    
    f = FilePath::open_file("vector.mtx", "w");
    saveVector(f, dimension(), vRHS);
    fclose(f);
}


/**
 save vectors and matrices in a text-based sparse formats
 */
void Meca::exportSystem() const
{
#if SEPARATE_RIGIDITY_TERMS
    std::clog << "incorrect dump since SEPARATE_RIGIDITY_TERMS is defined\n";
#endif
    FILE * f = FilePath::open_file("ord.txt", "w");
    fprintf(f, "%lu %i %lu\n", dimension(), DIM, sizeof(real));
    fclose(f);
    
    f = FilePath::open_file("stp.txt", "w");
    fprintf(f, "%f %f\n", tau_, tolerance_);
    fclose(f);
    
    f = FilePath::open_file("mob.txt", "w");
    saveMobility(f);
    fclose(f);
    
    f = FilePath::open_file("obj.txt", "w");
    saveObjectID(f);
    fclose(f);
    
    std::ofstream os("sol.txt");
    VecPrint::dump(os, dimension(), vPTS);
    os.close();
    
    os.open("rhs.txt");
    VecPrint::dump(os, dimension(), vRHS);
    os.close();
    
#if USE_ISO_MATRIX
    os.open("iso.txt");
    mISO.printSparse(os, 0);
    os.close();
#endif
    
    os.open("full.txt");
    mFUL.printSparse(os, 0);
    os.close();
        
    size_t alc = 0;
    for ( Mecable const* mec : mecables )
        alc = std::max(alc, mec->nbPoints());

    real * tmp1 = new_real(DIM*alc);
    real * tmp2 = new_real(DIM*DIM*alc*alc);
    
    os.open("diag.txt");
    
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


//------------------------------------------------------------------------------
#pragma mark - Binary Export

static void dumpVector(FILE * fp, size_t dim, real* vec, bool nat)
{
    static float * low = nullptr;
    static size_t alc = 0;
    
    if ( !fp )
    {
        delete[] low;
        low = nullptr;
        return;
    }
    if ( !nat && std::is_same<real, double>::value )
    {
        if ( dim > alc )
        {
            delete[] low;
            low = new float[dim];
            alc = dim;
        }
        copy_real(dim, vec, low);
        fwrite(low, sizeof(float), dim, fp);
    }
    else
        fwrite(vec, sizeof(real), dim, fp);
}


void Meca::dumpObjectID(FILE * fp) const
{
    uint32_t * vec = new uint32_t[largestMecable()];
    
    uint32_t i = 1;
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = mec->nbPoints();
        for ( size_t p = 0; p < nbp; ++p )
            vec[p] = i;
        for ( int d = 0; d < DIM; ++d )
            fwrite(vec, sizeof(uint32_t), nbp, fp);
        ++i;
    }
    
    delete[](vec);
}


void Meca::dumpMobility(FILE * fp, bool nat) const
{
    real * vec = new_real(largestMecable());
    
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = mec->nbPoints();
        const real val = mec->pointMobility();
        for ( size_t p=0; p < nbp; ++p )
            vec[p] = val;
        for ( int d = 0; d < DIM; ++ d )
            dumpVector(fp, nbp, vec, nat);
    }
    
    free_real(vec);
}


/**
 Save the full matrix associated with multiply(), in binary format
 */
void Meca::dumpMatrix(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1;
        multiply(src, res);
        dumpVector(fp, dim, res, nat);
        src[ii] = 0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the elasticity matrix, in binary format
 */
void Meca::dumpElasticity(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1;
        
        mFUL.vecMul(src, res);
#if USE_ISO_MATRIX
        mISO.VECMULADDISO(src, res);
#endif
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
        
        dumpVector(fp, dim, res, nat);
        src[ii] = 0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the projection matrix multiplied by the mobility, in binary format
 */
void Meca::dumpProjection(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * vec = new_real(dim);
        
    for ( size_t i = 0; i < dim; ++i )
    {
        zero_real(dim, vec);
        vec[i] = 1;
        
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            // this includes the mobility, but not the time_step:
            mec->projectForces(vec+inx, vec+inx);
            blas::xscal(DIM*mec->nbPoints(), mec->leftoverMobility(), vec+inx, 1);
        }
        // write column to fp directly:
        dumpVector(fp, dim, vec, nat);
    }
    
    free_real(vec);
}


/**
 Save matrix associated with the preconditionner, in binary format
 This relies on `Meca::precondition()`, which may apply a dummy preconditionner
 */
void Meca::dumpPreconditionner(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * vec = new_real(dim);
    
    for ( size_t i = 0; i < dim; ++i )
    {
        zero_real(dim, vec);
        vec[i] = 1;
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            applyPreconditionner(mec, vec+inx);
        }
        dumpVector(fp, dim, vec, nat);
    }
    
    free_real(vec);
}


/**
 This dump the total matrix and some vectors in binary files.
 
 This MATLAB code should read the output:
 
     ord = load('ord.txt');
     time_step = load('stp.txt');
     precision = 'double' % or float?
     obj = fread(fopen('obj.bin'), ord, 'uint32');
     drg = fread(fopen('drg.bin'), ord, precision);
     sys = fread(fopen('sys.bin'), [ord, ord], precision);
     ela = fread(fopen('ela.bin'), [ord, ord], precision);
     mob = fread(fopen('mob.bin'), [ord, ord], precision);
     con = fread(fopen('con.bin'), [ord, ord], precision);
     pts = fread(fopen('pts.bin'), ord, precision);
     rhs = fread(fopen('rhs.bin'), ord, precision);
     sol = fread(fopen('sol.bin'), ord, precision);
 
 To display the matrices:

     imshow(abs(sys))
     imshow(abs(ela))
 
 You can then compare the results with matlab's own iterative method,
 and compare the result using a scatter plot:
 
     x = bicgstab(sys, rhs, 0.001, ord);
     plot(x, sol, '.');
 
 */
void Meca::dumpSystem(bool nat) const
{
    FILE * f = FilePath::open_file("ord.txt", "w");
    fprintf(f, "%lu %i %lu\n", dimension(), DIM, sizeof(real));
    fclose(f);
    
    f = FilePath::open_file("stp.txt", "w");
    fprintf(f, "%.12f %.12f\n", tau_, tolerance_);
    fclose(f);
    
    f = FilePath::open_file("mob.bin", "wb");
    dumpMobility(f, nat);
    fclose(f);
    
    f = FilePath::open_file("obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    f = FilePath::open_file("rhs.bin", "wb");
    dumpVector(f, dimension(), vRHS, nat);
    fclose(f);
    
    f = FilePath::open_file("sol.bin", "wb");
    dumpVector(f, dimension(), vSOL, nat);
    fclose(f);
    
    f = FilePath::open_file("pts.bin", "wb");
    dumpVector(f, dimension(), vPTS, nat);
    fclose(f);
    
    f = FilePath::open_file("sys.bin", "wb");
    dumpMatrix(f, nat);
    fclose(f);
    
    f = FilePath::open_file("ela.bin", "wb");
    dumpElasticity(f, nat);
    fclose(f);
    
    f = FilePath::open_file("prj.bin", "wb");
    dumpProjection(f, nat);
    fclose(f);
    
    f = FilePath::open_file("con.bin", "wb");
    dumpPreconditionner(f, nat);
    fclose(f);
    
    dumpVector(nullptr, 0, nullptr, nat);
}


//------------------------------------------------------------------------------
#pragma mark - Matrix Bitmap Export

#include "../base/save_bitmap.cc"
static const uint16_t BitsPerPixel = 1;
static const uint8_t white = 1;
static const uint8_t gray = 1;


template < typename MatrixClass >
void setMatrixBitmap(uint8_t* bytes, size_t bpr, size_t nbv, MatrixClass const& MAT, size_t ORD)
{
    memset(bytes, 0, bpr*nbv);
    
    for ( size_t j = 0; j < nbv; ++j )
    {
        size_t jj = j * ORD;
        for ( size_t n = 0; n < MAT.column_size(jj); ++n )
        {
            size_t i = MAT.column_index(jj, n) / ORD;
            // swap i,j and flip i to display matrix properly on screen
            set_bitmap(bytes, bpr, BitsPerPixel*j, nbv-1-i, white);
        }
    }
}

void Meca::saveMatrixBitmaps() const
{
    static size_t cnt = 0;
    const size_t nbv = nbVertices();
    const size_t ims = bitmap_size(nbv, nbv, BitsPerPixel);
    const size_t bpr = bytes_per_row(nbv, BitsPerPixel);
    char str[32] = { 0 };
    
    uint8_t* bytes = new uint8_t[ims];
#if USE_ISO_MATRIX
    snprintf(str, sizeof(str), "iso%08lu.bmp", cnt);
    FILE * f = fopen(str, "w");
    if ( f ) {
        if ( !ferror(f) ) {
            setMatrixBitmap(bytes, bpr, nbv, mISO, 1);
            save_bitmap(f, bytes, nbv, nbv, BitsPerPixel);
        }
        fclose(f);
    }
#endif
    snprintf(str, sizeof(str), "ful%08lu.bmp", cnt);
    FILE * g = fopen(str, "w");
    if ( g ) {
        if ( !ferror(g) ) {
            setMatrixBitmap(bytes, bpr, nbv, mFUL, DIM);
            // mark pixels on the diagonals to indicate where each block starts:
            for ( Mecable * mec : mecables )
            {
                size_t i = mec->matIndex();
                set_bitmap(bytes, bpr, BitsPerPixel*i, nbv-1-i, gray);
            }
            save_bitmap(g, bytes, nbv, nbv, BitsPerPixel);
        }
        fclose(g);
    }
    delete[] bytes;
    ++cnt;
}

