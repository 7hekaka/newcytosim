// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "xtbsv.h"
#include "xtrsm.h"

/// if you like Alsatian specialities, set this to a prime number
#define CHOUCROUTE 7


/// leading dimension of the banded matrix used for iso symmetric blocks
constexpr size_t ISOB_LDD = 3;

/*
 number of off-diagonals and leading dimension for the non-isotropic banded matrix
 We use a band preconditionner with 2*DIM off-diagonals to include the
 'diagonal' terms from the blocks that are offset by 2 from the matrix diagonal
 */
constexpr size_t BAND_NUD = 2*DIM;
constexpr size_t BAND_LDD = BAND_NUD+DIM;


//------------------------------------------------------------------------------
#pragma mark - Apply Blocks

/// apply banded symmetric isotropic preconditionner block
static inline void applyPrecondIsoB(Mecable const* mec, real* Y)
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
static inline void applyPrecondIsoS(Mecable const* mec, real* Y)
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
static inline void applyPrecondIsoP(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
    //iso_xgetrsN<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
    alsatian_xgetrsN<DIM>(nbp, mec->block(), nbp, mec->pivot(), Y);
}


/// apply banded symmetric preconditionner block
static inline void applyPrecondBand(Mecable const* mec, real* Y)
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
static inline void applyPrecondHalf(Mecable const* mec, real* Y)
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
static inline void applyPrecondFull(Mecable const* mec, real* Y)
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
#pragma mark - Debug Preconditionner


static void printPreconditionnerBlock(Mecable* mec, size_t sup)
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
    std::clog << std::setw(8) << mec->reference() << "  " << std::setw(4) << bks;
    std::clog << " | B - K | = " << err << ' ';

    if ( err > bks * bks * REAL_EPSILON )
    {
        //VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        extractBlock(mec, wrk);
        
        const size_t S = std::min(9UL, bks);
        std::cout << std::endl;
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
    
    if ( 1 == mec->blockType() )
    {
        // chose initial vector for power iteration
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        real eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, mat);
        std::clog << "  eigen(1-PM) = " << eig;
    }
    
    std::clog << '\n';

    if ( err > 1 )
    {
        std::cout << "\n";
        // print preconditionner block for visual inspection:
        const size_t S = std::min(16UL, bks);
        VecPrint::full("matrix", S, S, wrk, bks);
        VecPrint::full("precond", S, S, blk, bks);
        VecPrint::full("precond * matrix", S, S, mat, bks);
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

//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner Block

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


//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner

/*
 Here the preconditionner blocks are calculated,
 according to Simul::precondition
 */
void Meca::computePreconditionner()
{
    bump_ = 0;
    switch( precond_ )
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
#if EXPERIMENTAL_PRECONDITIONNERS
        case 7:
            renewPreconditionner(span);
            break;
#endif
        default:
            throw InvalidParameter("unknown `precondition' value");
            break;
    }
    if ( bump_ > 0 )
        Cytosim::log << "failed to compute " << bump_ << " / " << mecables.size() << " preconditionner blocks\n";
}



//------------------------------------------------------------------------------
#pragma mark - Experimental Methods

#if EXPERIMENTAL_PRECONDITIONNERS

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
