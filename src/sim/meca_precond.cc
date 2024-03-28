// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "xtbsv.h"
#include "xtrsm.h"

#ifdef __APPLE__
#  include <mach/mach_time.h>
unsigned long machine_time()
{
    // the units of this counter is not specified, but seems close to nanosec
    return mach_absolute_time() >> 20;
}
#else
#  include <sys/time.h>
unsigned long machine_time()
{
    timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return 1000 * (unsigned long)tv.tv_sec + tv.tv_nsec / 1000;
}
#endif

/// if you like Alsatian specialities, set this to a prime number
#define CHOUCROUTE 7
#define SAUERKRAUT 7

/// leading dimension of the banded matrix used for iso symmetric blocks
constexpr size_t ISOB_KD = 2;
constexpr size_t ISOB_LDD = 3;

/*
 number of off-diagonals and leading dimension for the non-isotropic banded matrix
 We use a band preconditionner with 2*DIM off-diagonals to include the near-
 diagonal terms from the blocks that are offset by 2 from the matrix diagonal
 */
constexpr size_t BAND_NUD = 6;

// should allocate to also hold the true diagonal: BAND_LDD > BAND_NUD
constexpr size_t BAND_LDD = BAND_NUD+2;


//------------------------------------------------------------------------------
#pragma mark - Apply Blocks

/// apply banded symmetric isotropic preconditionner block
static inline void applyPrecondIsoB(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();

#if CHOUCROUTE
    assert_true( ISOB_KD == 2 );
    alsatian_iso_xpbtrsL<DIM>(nbp, mec->pblock(), ISOB_LDD, Y);
#else
    /*
     we cannot call lapack::DPBTRS('L', bks, KD, 1, mec->pblock(), KD+1, Y, bks, &info)
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     But calling DTBSV gets the required work done.
     */
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', nbp, ISOB_KD, mec->pblock(), ISOB_LDD, Y+d, DIM);
        blas::xtbsv('L', 'T', 'N', nbp, ISOB_KD, mec->pblock(), ISOB_LDD, Y+d, DIM);
    }
#endif
}


/// apply symmetric isotropic preconditionner block
static inline void applyPrecondIsoS(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
    /*
     we cannot call lapack::DPOTRS('L', nbp, mec->pblock(), nbp, Y, DIM, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     */
#if 0
    real * tmp = new_real(DIM*nbp);
    copy_real(DIM*nbp, Y, tmp);
    iso_xpotrsL<DIM>(nbp, mec->pblock(), nbp, tmp);
    std::clog << "\n "; VecPrint::print(nbp, tmp, 3, 100.0);
    free_real(tmp);
#endif
#if CHOUCROUTE
    alsatian_xpotrsL<DIM>(nbp, mec->pblock(), nbp, Y);
#else
    iso_xpotrsL<DIM>(nbp, mec->pblock(), nbp, Y);
#endif
    //std::clog << "\nL"; VecPrint::print(nbp, Y, 3, 100.0);
}


/// apply non-symmetric but isotropic preconditionner block
static inline void applyPrecondIsoP(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
#if SAUERKRAUT
    alsatian_xgetrsN<DIM>(nbp, mec->pblock(), nbp, mec->pivot(), Y);
#else
    iso_xgetrsN<DIM>(nbp, mec->pblock(), nbp, mec->pivot(), Y);
#endif
}


/// apply banded symmetric preconditionner block
static inline void applyPrecondBand(Mecable const* mec, real* Y)
{
    const int bks = DIM * mec->nbPoints();
    assert_true( (int)BAND_NUD < bks );
#if SAUERKRAUT
    alsatian_xpbtrsLK<BAND_NUD>(bks, mec->pblock(), BAND_LDD, Y);
#elif 1
    blas_xtbsvLN<'N'>(bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
    blas_xtbsvLT<'N'>(bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
#elif 1
    blas::xtbsv('L', 'N', 'N', bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
    blas::xtbsv('L', 'T', 'N', bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
#else
    int info = 0;
    lapack::xpbtrs('L', bks, BAND_NUD, 1, mec->pblock(), BAND_LDD, Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply full size symmetric preconditionner block
static inline void applyPrecondHalf(Mecable const* mec, real* Y)
{
    const int bks = DIM * mec->nbPoints();
#if SAUERKRAUT
    // assuming that diagonal terms of the preconditionner block have been inverted:
    alsatian_xpotrsL(bks, mec->pblock(), bks, Y);
#elif 1
    iso_xpotrsL<1>(bks, mec->pblock(), bks, Y);
#else
    int info = 0;
    lapack::xpotrs('L', bks, 1, mec->pblock(), bks, Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply full size non-symmetric preconditionner block
static inline void applyPrecondFull(Mecable const* mec, real* Y)
{
    const int bks = DIM * mec->nbPoints();
#if CHOUCROUTE && USE_SIMD
    // assuming that diagonal terms of the preconditionner block have been inverted:
    alsatian_xgetrsN_SSE(bks, mec->pblock(), bks, mec->pivot(), Y);
#elif CHOUCROUTE
    alsatian_xgetrsN(bks, mec->pblock(), bks, mec->pivot(), Y);
#elif 1
    // translated LAPACK's reference code:
    lapack_xgetrsN(bks, mec->pblock(), bks, mec->pivot(), Y);
#else
    // using LAPACK's library
    int info = 0;
    lapack::xgetrs('N', bks, 1, mec->pblock(), bks, mec->pivot(), Y, bks, &info);
    assert_true(info==0);
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Apply Preconditionner


/// apply preconditionner block corresponding to Mecable
static inline void applyPreconditionner(Mecable const* mec, real* Y)
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
    unsigned long rdt = machine_time();
    if ( Y != X )
        copy_real(dimension(), X, Y);
    
    #pragma omp parallel for num_threads(NUM_THREADS)
    for ( Mecable const* mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
#if RECYCLED_PRECONDITIONNER
        if ( mec->blockType() == 7 )
            mec->blockMultiply(X+inx, Y+inx);
        else
#endif
            applyPreconditionner(mec, Y+inx);
    }
    // CPU time is recorded, but only for info, so this code can be removed
    cycles_ += machine_time() - rdt;
}


/// total allocated memory size for preconditionner
size_t Meca::preconditionnerSize() const
{
    size_t res = 0;
    for ( Mecable const* mec : mecables )
        res += mec->blockLimit();
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Debug Preconditionner


static void printPreconditionnerBlock(Mecable* mec, size_t sup)
{
    const size_t bks = DIM * mec->nbPoints();
    size_t S = std::min(bks, sup);
    real * blk = mec->pblock();

    VecPrint::full("Diagonal block"+std::to_string(bks), S, S, mec->pblock(), bks, 3);
    
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
    
    real mag = blas::nrm2(bks*bks, blk);
    blas::xaxpy(bks*bks, -1.0, blk, 1, wrk, 1);
    real err = blas::nrm2(bks*bks, wrk);
 
    std::clog << "verifyBlock ";
    std::clog << std::setw(8) << mec->reference() << "  " << std::setw(4) << bks;
    std::clog << " | B - K | = " << err << ' ';

    if ( err > mag * 1e-6 )
    {
        //VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        extractBlock(mec, wrk);
        
        std::cout << std::endl;
        VecPrint::full("extracted", bks, bks, wrk, bks, 3);
        VecPrint::full("computed", bks, bks, blk, bks, 3);
    }
    
    free_real(wrk);
}


/**
 Multiply here `blk` with the dynamic block extracted by extractBlock()
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
    
    if ( 1 == mec->blockType() )
    {
        // use power iterations to estimate largest eigenvalue
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        real eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, mat);
        std::clog << "  eigen(1-PM) = " << eig;
    }
    
    std::clog << '\n';

    if ( err > 1 )
    {
        std::cout << "\n";
        // print preconditionner block for visual inspection:
        VecPrint::full("matrix", bks, bks, wrk, bks);
        VecPrint::full("precond", bks, bks, blk, bks);
        VecPrint::full("precond * matrix", bks, bks, mat, bks);
    }
    free_real(vec);
    free_real(mat);
    free_real(wrk);
}


//------------------------------------------------------------------------------
#pragma mark - Extract Diagonal Block

/**
 Build a banded block that can be used as a preconditionner:
 
     I - time_step * mobility * Rigidity
 
 This block is square, symmetric, definite positive and well-behaved
 */
void Meca::getIsoBandedBlock(const Mecable * mec, real* res, size_t kd, size_t ldd) const
{
    const size_t nbp = mec->nbPoints();

    const real beta = -tau_ * mec->pointMobility();
    real jR = mec->jointRigidity();
    if ( jR )
    {
        if ( ldd != 3 )
            zero_real(ldd*nbp, res);
        setBendingRigidity<1>(res, ldd-1, nbp, beta*jR);
    }
    else
        zero_real(ldd*nbp, res);
    
    //VecPrint::full("rigidity band", 3, nbp, res, ldd, 1);

    /*
     The matrix `res` is stored in 'packed symmetric banded storage':
     usually, mat(i, j) is stored in mat[i+ldd*j]
     but with banded storage, mat(i, j) is stored in mat[i-j+ldd*j] for i > j
     So we use below ldd-1
     */

#if USE_ISO_MATRIX
    mISO.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, 1, kd);
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalTrace(beta/DIM, res, kd, mec->matIndex(), nbp, DIM, kd, false);

    // add Identity matrix to band storage:
    for ( size_t i = 0; i < nbp; ++i )
        res[ldd*i] += 1;
}


/**
 Extract block of reduced dimension that does not include projection:
 
     I - time_step * point_mobility ( mISO + mFUL )
 
 The result is constructed by using efficient methods from mISO and mFUL
 The block is symmetric and can be factorized by Cholesky, which may fail
 */
void Meca::getIsoBlock(const Mecable * mec, real* res) const
{
    const size_t nbp = mec->nbPoints();
    
    zero_real(nbp*nbp, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    real jR = mec->jointRigidity();
    if ( jR )
        //addBendingRigidityLower<1>(res, nbp, mec->nbPoints(), mec->fiberRigidity());
        addBendingRigidity<1>(res, nbp, mec->nbPoints(), jR);
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, nbp, mec->matIndex(), nbp, 1, 1);
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalTrace(1.0/DIM, res, nbp, mec->matIndex(), nbp, DIM, nbp, true);

    // the projection is not called, so we scale by mobility
    const real beta = -tau_ * mec->pointMobility();

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
    real jR = mec->jointRigidity();
    if ( jR )
    {
        setBendingRigidity<DIM>(res, ldd-1, nbp, beta*jR);
        //VecPrint::full("Rigidity block ", bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, 1, rank/DIM);
#endif
    // VecPrint::full("\niso", ldd, bks, res, ldd);
    copy_lower_subspace<DIM>(bks, res, ldd-1, rank);
    // VecPrint::full("\ncopy_subspace", ldd, bks, res, ldd);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, DIM, rank);
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
    real jR = mec->jointRigidity();
    if ( jR )
    {
        addBendingRigidityLower<DIM>(res, bks, mec->nbPoints(), jR);
        //VecPrint::full("Rigidity block", bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, 1, DIM);
#endif
    copy_lower_subspace<DIM, true>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
    
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
 
     I - time_step * P ( mISO + mFUL + diffP )
 
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
    real jR = mec->jointRigidity();
    if ( jR )
    {
        addBendingRigidityLower<DIM>(res, bks, mec->nbPoints(), jR);
        //VecPrint::full("Rigidity block", bks, bks, res, bks, 0);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, 1, DIM);
#endif
    copy_lower_subspace<DIM, true>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
    
    //VecPrint::full("mISO+mFUL block", bks, bks, res, bks);
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(res);
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

//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner Block

/**
 Compute a preconditionner block corresponding to 'mec'
 The dimension is reduced by DIM and banded with diagonal + 2 off-diagonals
 This block is usually symmetric definite positive, and is factorized by Cholesky's method
 */
void Meca::computePrecondIsoB(Mecable* mec)
{
    const size_t nbp = mec->nbPoints();
    mec->blockSize(DIM*nbp, ISOB_LDD*nbp, 0);
    
    /**
     Factorize banded matrix with Andre-Louis Cholesky's method
     born 15.10.1875 in Montguyon, France
     died 31.08.1918 in Bagneux, following wounds received in battle.
     */
    int bt, info = 0;
    if ( ISOB_LDD <= nbp )
    {
#if 0
        // visual comparison with getIsoBlock()
        mec->blockSize(DIM*nbp, std::max(ISOB_LDD, nbp)*nbp, 0);
        getIsoBlock(mec, mec->pblock());
        VecPrint::full("iso symmetric", nbp, nbp, mec->pblock(), nbp, 2);
#endif
        getIsoBandedBlock(mec, mec->pblock(), ISOB_KD, ISOB_LDD);
        //VecPrint::full("banded", ISOB_LDD, nbp, mec->pblock(), ISOB_LDD, 2);

        // calculate Banded Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpbtf2L(nbp, ISOB_KD, mec->pblock(), ISOB_LDD, &info);
#else
        lapack::xpbtf2('L', nbp, ISOB_KD, mec->pblock(), ISOB_LDD, &info);
#endif
        bt = 1;
    }
    else
    {
        getIsoBlock(mec, mec->pblock());
        // calculate Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpotf2L(nbp, mec->pblock(), nbp, &info);
#else
        lapack::xpotf2('L', nbp, mec->pblock(), nbp, &info);
#endif
        //VecPrint::full("isoB", nbp, nbp, mec->pblock(), nbp, 2);
        bt = 2;
    }

    if ( 0 == info )
    {
        mec->blockType(bt);
        //VecPrint::full("factorized", 3, nbp, mec->pblock(), ISOB_LDD, 1);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute Band Preconditionner block of size " << nbp << "\n";
        ++bump_;
        //computePrecondIsoS(mec);
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
#if 0
    const size_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, bks);
    getFullBlock(mec, mec->pblock());
    project_matrix<DIM>(nbp, mec->pblock(), bks, mec->pblock(), nbp);
    VecPrint::full("projected: ", nbp, nbp, mec->pblock(), nbp, 2);
#endif

    mec->blockSize(DIM*nbp, nbp*nbp, 0);
    
    //getIsoBandedBlock(mec, mec->pblock(), ISOB_LDD);
    //VecPrint::full("banded preconditionner ", 3, nbp, mec->pblock(), ISOB_LDD, 1);

    getIsoBlock(mec, mec->pblock());
    
    //VecPrint::full("iso symmetric", nbp, nbp, mec->pblock(), nbp, 2);

    // calculate Cholesky factorization:
    int info = 0;
#if CHOUCROUTE
    alsatian_xpotf2L(nbp, mec->pblock(), nbp, &info);
#else
    lapack::xpotf2('L', nbp, mec->pblock(), nbp, &info);
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
 Block of dimension reduced by DIM, including the projection
 The block is not symmetric and is factorized by LU decomposition
 */
void Meca::computePrecondIsoP(Mecable* mec)
{
    const size_t nbp = mec->nbPoints();
    int info = 0;

    const size_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, bks);
    double * blk = mec->pblock();

    // we claim too much memory here, but this preconditionner is not good anyhow:
    getFullBlock(mec, blk);
    project_matrix<DIM>(nbp, blk, bks, blk, nbp);

    //VecPrint::full("isop preconditionner", nbp, nbp, mec->pblock(), nbp, 2);

    // calculate LU factorization:
#if SAUERKRAUT
    alsatian_xgetf2(nbp, blk, nbp, mec->pivot(), &info);
#else
    lapack::xgetf2(nbp, nbp, blk, nbp, mec->pivot(), &info);
#endif

    if ( 0 == info )
    {
        mec->blockType(3);
        //checkBlock(mec, blk);
#if SAUERKRAUT && REAL_IS_DOUBLE
    convert_to_floats(bks*bks, blk, (float*)blk);
#endif
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
    mec->blockSize(bks, std::min(BAND_LDD, bks)*bks, 0);
    double * blk = mec->pblock();

    int bt, info = 0;

    if ( BAND_LDD < bks )
    {
        getBandedBlock(mec, blk, BAND_LDD, BAND_NUD);
        
#if ( 0 )
        VecPrint::full("band block "+std::to_string(bks), BAND_LDD, bks, blk, BAND_LDD);
        mec->blockSize(bks, bks*bks, 0);
        getHalfBlock(mec, blk);
        VecPrint::full("half block", bks, bks, blk, bks);
#endif
        
        // calculate Cholesky factorization for band storage:
#if SAUERKRAUT
        alsatian_xpbtf2L(bks, BAND_NUD, blk, BAND_LDD, &info);
#else
        lapack::xpbtf2('L', bks, BAND_NUD, blk, BAND_LDD, &info);
#endif
        bt = 4;
    }
    else
    {
        getHalfBlock(mec, blk);
        // calculate Cholesky factorization:
#if SAUERKRAUT
        alsatian_xpotf2L(bks, blk, bks, &info);
#else
        lapack::xpotf2('L', bks, blk, bks, &info);
#endif
        bt = 5;
    }
    
    if ( 0 == info )
    {
        mec->blockType(bt);
        //checkBlock(mec, blk);
#if SAUERKRAUT && REAL_IS_DOUBLE
        convert_to_floats(bks*bks, blk, (float*)blk);
#endif
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute band Preconditionner block of size " << bks << "\n";
        ++bump_;
    }
}

/**
 Compute preconditionner block corresponding to 'mec'
 This block is symmetric, and factorized by Cholesky's method.
 */
void Meca::computePrecondHalf(Mecable* mec, real* tmp)
{
    const size_t bks = DIM * mec->nbPoints();
#if SAUERKRAUT && REAL_IS_DOUBLE
    mec->blockSize(bks, 4+bks*bks/2, bks);
    // use temporary memory to build matrix block:
    real * blk = tmp;
#else
    mec->blockSize(bks, bks*bks, 0);
    real * blk = mec->pblock();
#endif
    getHalfBlock(mec, blk);
    
    //VecPrint::full("half block "+std::to_string(bks), bks, bks, mec->pblock(), bks);

    int info = 0;
    // calculate Cholesky factorization:
#if SAUERKRAUT
    alsatian_xpotf2L(bks, blk, bks, &info);
#else
    lapack::xpotf2('L', bks, blk, bks, &info);
#endif

    if ( 0 == info )
    {
        mec->blockType(5);
        //checkBlock(mec, blk);
        //VecPrint::full("half", bks, bks, blk, bks);
#if SAUERKRAUT && REAL_IS_DOUBLE
        convert_to_floats(bks*bks, blk, (float*)mec->pblock());
#endif
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
void Meca::computePrecondFull(Mecable* mec, real* tmp)
{
    const size_t bks = DIM * mec->nbPoints();
    
#if CHOUCROUTE && REAL_IS_DOUBLE
    mec->blockSize(bks, 4+bks*bks/2, bks);
    real * blk = tmp;
#else
    mec->blockSize(bks, bks*bks, bks);
    real * blk = mec->pblock();
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
        //if ( bks < 4 ) VecPrint::full("full", bks, bks, blk, bks);
#if CHOUCROUTE && REAL_IS_DOUBLE
        convert_to_floats(bks*bks, blk, (float*)mec->pblock());
#endif
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute full Preconditionner block of size " << bks << "\n";
        ++bump_;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner

/*
 Here the preconditionner blocks are calculated,
 according to Simul::precondition
 */
void Meca::computePreconditionner()
{
    bump_ = 0;
    size_t sup = ( 1 + DIM * largestMecable() ) & ~1;
    real * tmp = new_real(4+sup*sup);

    switch( precond_ )
    {
        case 0:
            for ( Mecable * mec : mecables )
                mec->blockType(0);
            break;
        case 1:
            for ( Mecable * mec : mecables )
            {
                if ( mec->tag() == Solid::TAG )
                    computePrecondFull(mec, tmp);
                else
                    computePrecondIsoB(mec);
            }
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
                computePrecondHalf(mec, tmp);
            break;
        case 6:
            for ( Mecable * mec : mecables )
                computePrecondFull(mec, tmp);
            break;
        case 7:
            for ( Mecable * mec : mecables )
            {
                if ( mec->tag() == Bead::TAG )
                    computePrecondIsoS(mec);
                else
                    computePrecondFull(mec, tmp);
            }
            break;
#if RECYCLED_PRECONDITIONNER
        case 9:
            renewPreconditionner(span);
            break;
#endif
        default:
            throw InvalidParameter("unknown `precondition' value");
            break;
    }
    free_real(tmp);
}



//------------------------------------------------------------------------------
#pragma mark - Experimental Methods

#if RECYCLED_PRECONDITIONNER

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
    
    real* blk = mec->pblock();
    real* vec = new_real(bks);

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
            free_real(vec);
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
    free_real(vec);
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
 Initialize Mecable::blockMatrix() as the preconditionner block,
 using the factorization given by 'blk' and 'piv'
 */
static void convertPreconditionner(Mecable* mec, real* blk, int* piv, real* wrk)
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

    if (( mec->blockType() != 7 ) | ( age >= span ))
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
#endif
