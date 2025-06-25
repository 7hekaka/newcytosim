// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN 30.06.2020

#include <sys/time.h>
#include <iostream>
#include <limits>

#define DIM 3

#include "real.h"
#include "timer.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "simd.h"
#include "assert_macro.h"
#include "cytoblas.h"
#include "xtrsm.h"
#include "xtbsv.h"

#define UNFINISHED_PRECONDITIONNER 0

/// convert doubles to floats
void convert_to_floats(size_t cnt, double const* src, float* dst)
{
    #pragma omp simd
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = (float)src[i];
}

template < size_t ORD >
static void clarify_matrix(size_t lin, size_t col, real* mat, size_t ldd)
{
    for ( size_t jj = 0; jj < col; ++jj )
    {
        real * ptr = mat + jj * ldd;
        for ( size_t ii = 0; ii < lin; ++ii )
        {
            if ( ii % ORD )
                ptr[ii] = 0;
        }
    }
}

template < size_t ORD >
static void modify_matrix(size_t lin, size_t col, real* mat, size_t ldd)
{
    real * ptr = mat;
    for ( size_t jj = 0; jj < col; ++jj )
    {
        for ( size_t ii = 0; ii < lin; ++ii )
            *ptr++ = mat[ORD*ii+jj*ldd];
    }
}


/// print 12 scalars from `vec[]` of dimension `num`
inline void print_vector(int num, real const* vec)
{
    if ( num > 12 )
    {
        VecPrint::print(6, vec, 3);
        fprintf(stderr, "...");
        VecPrint::print(6, vec+num-6, 3);
        fprintf(stderr, " |");
        VecPrint::print(4, vec+num, 1);
    }
    else
    {
        VecPrint::print(num, vec, 3);
        fprintf(stderr, " |");
        VecPrint::print(4, vec+num, 1);
    }
    real sum = vec[0];
    for ( int i = 1; i < num; ++i )
        sum += vec[i];
    fprintf(stderr, "  %+22.16f ", sum);
}

void nan_spill(real * dst)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( size_t i = 0; i < 4; ++i )
        dst[i] = n;
}

template < void (*FUNC)(int, real const*, int, real*) >
void check(int N, int ORD, real const* S, real const* AB, int LDA, real* B, char const str[], size_t rep, size_t sub=128)
{
    printf("\n");
    // verification:
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    FUNC(N, AB, LDA, B);
    print_vector(ORD*N, B);
    // performance:
    tick();
    for ( size_t n = 0; n < rep; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < sub; ++u )
            FUNC(N, AB, LDA, B);
    }
    printf(" %-14s cpu %7.0f", str, tock());
}

template < void (*FUNC)(int, real const*, int, real*) >
void multi(int N, int ORD, real const* S, real const* AB, int LDA, real* B, char const str[], size_t rep, size_t sub)
{
    size_t BLK = N * N + 4;
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    tick();
    for ( size_t n = 0; n < rep; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < sub; ++u )
            FUNC(N, AB+u*BLK, LDA, B);
    }
    printf("\n uncached %-14s cpu %7.0f", str, tock());
}

//------------------------------------------------------------------------------
#pragma mark - Isotropic Solvers for interleaved vectors

const size_t KD = 2;

void iso0(int N, real const* AB, int LDA, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', N, KD, AB, LDA, B+d, DIM);
        blas::xtbsv('L', 'T', 'N', N, KD, AB, LDA, B+d, DIM);
    }
}

void iso1(int N, real const* AB, int LDA, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas_xtbsvLN<'C'>(N, KD, AB, LDA, B+d, DIM);
        blas_xtbsvLT<'C'>(N, KD, AB, LDA, B+d, DIM);
    }
}

void iso2(int N, real const* AB, int LDA, real* B)
{
    alsatian_iso_xtbsvLNN<DIM>(N, KD, AB, LDA, B);
    alsatian_iso_xtbsvLTN<DIM>(N, KD, AB, LDA, B);
}

void iso3(int N, real const* AB, int LDA, real* B)
{
    alsatian_iso_xpbtrsL<DIM>(N, AB, LDA, B);
}

void iso5(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 3 ) && defined(__AVX__) && REAL_IS_DOUBLE
    alsatian_iso3_xtbsvLNN2K_AVX(N, AB, LDA, B);
    alsatian_iso3_xtbsvLTN2K_AVX(N, AB, LDA, B);
#elif ( DIM == 3 ) && USE_SIMD
    alsatian_iso3_xtbsvLNN2K_SIMD(N, AB, LDA, B);
    alsatian_iso3_xtbsvLTN2K_SIMD(N, AB, LDA, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_iso2_xtbsvLNN2K_SIMD(N, AB, LDA, B);
    alsatian_iso2_xtbsvLTN2K_SIMD(N, AB, LDA, B);
#elif ( DIM == 1 )
    alsatian_xtbsvLNN2K(N, AB, LDA, B);
    alsatian_xtbsvLTN2K(N, AB, LDA, B);
#else
    alsatian_iso_xtbsvLNN<DIM>(N, 2, AB, LDA, B);
    alsatian_iso_xtbsvLTN<DIM>(N, 2, AB, LDA, B);
#endif
}

void isoLNN(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 1 )
    alsatian_xtbsvLNN2K(N, AB, LDA, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_iso2_xtbsvLNN2K_SIMD(N, AB, LDA, B);
#elif ( DIM == 2 )
    alsatian_iso2_xtbsvLNN2K(N, AB, LDA, B);
#elif ( DIM == 3 ) && defined(__AVX__)
    alsatian_iso3_xtbsvLNN2K_AVX(N, AB, LDA, B);
#elif ( DIM == 3 ) && USE_SIMD
    alsatian_iso3_xtbsvLNN2K_SIMD(N, AB, LDA, B);
#else
    alsatian_iso_xtbsvLNN<DIM>(N, 2, AB, LDA, B);
#endif
}

void isoLTN(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 1 )
    alsatian_xtbsvLTN2K(N, AB, LDA, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_iso2_xtbsvLTN2K_SIMD(N, AB, LDA, B);
#elif ( DIM == 2 )
    alsatian_iso2_xtbsvLTN2K(N, AB, LDA, B);
#elif ( DIM == 3 ) && defined(__AVX__)
    alsatian_iso3_xtbsvLTN2K_AVX(N, AB, LDA, B);
#elif ( DIM == 3 ) && USE_SIMD
    alsatian_iso3_xtbsvLTN2K_SIMD(N, AB, LDA, B);
#else
    alsatian_iso_xtbsvLTN<DIM>(N, 2, AB, LDA, B);
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testISO(int N, size_t rep)
{
    std::cout << DIM << "D rank 2 Cholesky factorization & iso solve " << N << " points --- real " << sizeof(real);
    std::cout << " --- " << __VERSION__;
    const size_t LDAB = 3;

    real * AB = new_real(N*LDAB+4);
    real * S = new_real(N*DIM);
    real * B = new_real(N*DIM+4);

    zero_real(N*LDAB, AB);
    nan_spill(AB+N*LDAB);
    for ( int i = 0; i < N*DIM; ++i )
        S[i] = RNG.sreal();
    for ( int i = 0; i < N; ++i )
    {
        AB[  LDAB*i] = 1+0.5 * RNG.sreal();
        AB[1+LDAB*i] = 0.125 * RNG.sreal();
        AB[2+LDAB*i] = 0.125 * RNG.sreal();
    }
    int info;
    if ( 1 )
    {
        real * ABc = new_real(N*LDAB+4);
        copy_real(N*LDAB+4, AB, ABc);
        lapack::xpbtf2('L', N, 2, ABc, LDAB, &info);
        check<iso0>(N, DIM, S, ABc, LDAB, B, "blas:tbsv", rep);
        free_real(ABc);
    }
    // factorize, Alsatian's way:
    alsatian_xpbtf2L(N, 2, AB, LDAB, &info);
    
    //check<iso0>(NPTS, DIM, S, AB, B, "fail BLAS", rep);
    check<iso1>(N, DIM, S, AB, LDAB, B, "blas_pbtrsL", rep);
    check<iso2>(N, DIM, S, AB, LDAB, B, "alsa_pbtrsL<D>", rep);
    check<iso3>(N, DIM, S, AB, LDAB, B, "alsa_pbtrs", rep);
    check<iso5>(N, DIM, S, AB, LDAB, B, "alsa_pbtrs_SSE", rep);

#if 0 && REAL_IS_DOUBLE
    check<isoLNN>(N, DIM, S, AB, LDAB, B, "tbsvLNN3", rep);
    check<isoLTN>(N, DIM, S, AB, LDAB, B, "tbsvLTN3", rep);
#endif

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - Cholesky factorization

void pot0(int N, real const* AB, int LDA, real* B)
{
    iso_xpotrsL_lapack<DIM>(N, AB, LDA, B);
}

void pot1(int N, real const* AB, int LDA, real* B)
{
    alsatian_iso_xpotrsLref<DIM>(N, AB, LDA, B);
}

void pot2(int N, real const* AB, int LDA, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, AB, LDA, B);
    iso_xtrsmLLT<DIM,'I'>(N, AB, LDA, B);
}

void pot3(int N, real const* AB, int LDA, real* B)
{
    alsatian_iso_xpotrsL<DIM>(N, AB, LDA, B);
}

void pot4(int N, real const* AB, int LDA, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, (float*)AB, LDA, B);
    iso_xtrsmLLT<DIM,'I'>(N, (float*)AB, LDA, B);
}

void testPOTRS(int N, size_t rep)
{
    std::cout << "\n Cholesky factorization of full matrix & iso solve " << N << " points --- real " << sizeof(real);
    std::cout << " --- " << __VERSION__;
    
    const int LDA = ( N + 3 ) & ~3;
    real * AB = new_real(N*LDA+4);
    real * B = new_real(N*DIM+4);
    real * S = new_real(N*DIM);

    for ( int i = 0; i < N*DIM; ++i )
        S[i] = RNG.sreal();
    zero_real(N*LDA, AB);
    for ( int i = 0; i < N; ++i )
    {
        real r = 0.03125 * RNG.sreal();
        AB[LDA*i+i] = 1.5;
        if ( i < N-1 ) AB[LDA*i+i+1] = 0.0625 + r;
        if ( i < N-2 ) AB[LDA*i+i+2] = 0.0625 - r;
    }
    nan_spill(AB+N*LDA);
    //VecPrint::full("AB", N, N, AB, LDA);
    int info = 0;
    if ( 1 )
    {
        real * ABc = new_real(N*LDA+4);
        copy_real(N*LDA+4, AB, ABc);
        lapack::xpotf2('L', N, ABc, LDA, &info);
        check<pot0>(N, DIM, S, ABc, LDA, B, "blas:potrs", rep);
        free_real(ABc);
    }
    
    alsatian_xpotf2L(N, AB, LDA, &info);
    if ( info == 0 )
    {
        check<pot1>(N, DIM, S, AB, LDA, B, "alsa_potrsLref", rep);
        check<pot2>(N, DIM, S, AB, LDA, B, "iso_trsmLL<D>", rep);
        check<pot3>(N, DIM, S, AB, LDA, B, "alsa_potrs", rep);
#if REAL_IS_DOUBLE
        convert_to_floats(N*LDA, AB, (float*)AB);
        check<pot4>(N, DIM, S, AB, LDA, B, "iso_trsm_float", rep);
#endif
    }
    else
    {
        std::cout << "\n ERROR: failed factorization (" << info << ")   ";
        //VecPrint::full("AB", N, N, AB, LDA);
    }
    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - Banded matrices

const int RANK = 6;

// this assumes that the diagonal terms are not inverted
void uni0(int N, real const* AB, int LDA, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, RANK, AB, LDA, B, 1);
    blas::xtbsv('L', 'T', 'N', N, RANK, AB, LDA, B, 1);
}

void uni1(int N, real const* AB, int LDA, real* B)
{
    blas_xtbsvLN<'C'>(N, RANK, AB, LDA, B);
    blas_xtbsvLT<'C'>(N, RANK, AB, LDA, B);
}

void uni2(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN(N, RANK, AB, LDA, B);
    alsatian_xtbsvLTN(N, RANK, AB, LDA, B);
}

void uni3(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNNK<RANK>(N, AB, LDA, B);
    alsatian_xtbsvLTNK<RANK>(N, AB, LDA, B);
}

void uni6(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN6K(N, AB, LDA, B);
    alsatian_xtbsvLTN6K(N, AB, LDA, B);
}

void uni7(int N, real const* AB, int LDA, real* B)
{
    alsatian_xpbtrsLK<RANK>(N, AB, LDA, B);
}

// this gives wrong results
void uniLNB(int N, real const* AB, int LDA, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, RANK, AB, LDA, B, 1);
    //blas::xtbsv('L', 'T', 'N', N, RANK, AB, LDA, B);
}

void uniLN0(int N, real const* AB, int LDA, real* B)
{
    blas_xtbsvLN<'C'>(N, RANK, AB, LDA, B);
    //blas_xtbsvLT<'C'>(N, RANK, AB, LDA, B);
}

void uniLN1(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN(N, RANK, AB, LDA, B);
    //alsatian_xtbsvLTN(N, RANK, AB, LDA, B);
}

void uniLN2(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNNK<RANK>(N, AB, LDA, B);
    //alsatian_xtbsvLTNK<RANK>(N, AB, LDA, B);
}

void uniLN3(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN6K(N, AB, LDA, B);
}

#if REAL_IS_DOUBLE && USE_SIMD
void uniLN4(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN6K_SSE(N, AB, LDA, B);
}
#endif

void uniLN5(int N, real const* AB, int LDA, real* B)
{
    alsatian_striped_xtbsvLNN6K_3D(N, AB, 3, B);
}

// this gives wrong results
void uniLTB(int N, real const* AB, int LDA, real* B)
{
    //blas::xtbsv('L', 'N', 'N', N, RANK, AB, LDA, B);
    blas::xtbsv('L', 'T', 'N', N, RANK, AB, LDA, B, 1);
}

void uniLT0(int N, real const* AB, int LDA, real* B)
{
    //blas_xtbsvLN<'C'>(N, RANK, AB, LDA, B);
    blas_xtbsvLT<'C'>(N, RANK, AB, LDA, B);
}

void uniLT1(int N, real const* AB, int LDA, real* B)
{
    //alsatian_xtbsvLNN(N, RANK, AB, LDA, B);
    alsatian_xtbsvLTN(N, RANK, AB, LDA, B);
}

void uniLT2(int N, real const* AB, int LDA, real* B)
{
    //alsatian_xtbsvLNNK<RANK>(N, AB, LDA, B);
    alsatian_xtbsvLTNK<RANK>(N, AB, LDA, B);
}

void uniLT3(int N, real const* AB, int LDA, real* B)
{
    //alsatian_xtbsvLNN6K(N, AB, LDA, B);
    alsatian_xtbsvLTN6K(N, AB, LDA, B);
}

#if REAL_IS_DOUBLE && USE_SIMD
void uniLT4(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLTN6K_SSE(N, AB, LDA, B);
}
#endif

void uniLT5(int N, real const* AB, int LDA, real* B)
{
    alsatian_striped_xtbsvLTN6K_3D(N, AB, 3, B);
}

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testTBSV(int N, size_t rep)
{
    std::cout << "\nTBSV Banded matrix Cholesky factorization ---- rank ";
    std::cout << RANK << " ---- " << N << " points --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__;

    /// rank of diagonal matrices:
    const int BLDD = RANK+2;

    real * AB = new_real(N*BLDD+4);
    real * MC = new_real(N*BLDD+4);
    real * S = new_real(N);
    real * B = new_real(N+4);

    for ( int i = 0; i < N; ++i )
        S[i] = RNG.sreal();
    
    zero_real(N*BLDD, AB);
    nan_spill(AB+N*BLDD);
    for ( int i = 0; i < N; ++i )
    {
        real * col = AB + BLDD * i;
        real R = 5, r = 0.0125;
        for ( size_t j = 1; j <= RANK; ++j )
            col[j] = r * RNG.sreal();
        col[0] = R - r * RANK; //diagonal term
        col[DIM] -= R / 2;
        col[DIM*2] += R / 4;
    }
    //VecPrint::full("mat", RANK+1, N, AB, BLDD, 1);
    copy_real(N*BLDD, AB, MC);
    
    // factorize, Alsatian's way:
    int info = 0;
    alsatian_xpbtf2L(N, RANK, AB, BLDD, &info);
    //VecPrint::full("factorized", RANK+1, N, AB, BLDD, 1);
    //clarify_matrix<DIM>(RANK+1, N, AB, BLDD);
    //VecPrint::full("clarified", RANK+1, N, AB, BLDD, 1);

    if ( 0 )
    {
        lapack::xpbtf2('L', N, RANK, MC, BLDD, &info);
        check<uni0>(N, 1, S, MC, BLDD, B, "blas:tbsv", rep);
        // our implementations:
        check<uni1>(N, 1, S, AB, BLDD, B, "blas_tbsv", rep);
        check<uni2>(N, 1, S, AB, BLDD, B, "tbsvLxN", rep);
        check<uni3>(N, 1, S, AB, BLDD, B, "tbsvLxNK<KD>", rep);
        check<uni6>(N, 1, S, AB, BLDD, B, "LLN6K+LTN6K", rep);
        check<uni7>(N, 1, S, AB, BLDD, B, "alsa_xpbtrsLK", rep);
    }
    if ( 1 )
    {
        std::cout << "\nTBSVLN --- triangular band forward solve --- rank " << RANK;
        //check<uniLNB>(N, 1, S, AB, BLDD, B, "blas::xtbsv", rep);
        check<uniLN0>(N, 1, S, AB, BLDD, B, "blas_xtbsvLN", rep);
        check<uniLN1>(N, 1, S, AB, BLDD, B, "xtbsvLNN", rep);
        check<uniLN2>(N, 1, S, AB, BLDD, B, "LNNK<KD>", rep);
    }
    if ( DIM == 3 )
    {
        check<uniLN3>(N, 1, S, AB, BLDD, B, "LNN6K", rep);
#if REAL_IS_DOUBLE && USE_SIMD
        check<uniLN4>(N, 1, S, AB, BLDD, B, "LNN6K_SSE", rep);
#endif
#if UNFINISHED_PRECONDITIONNER
        copy_real(N*BLDD, AB, MC);
        modify_matrix<3>(3, N, MC, BLDD);
        //VecPrint::full("modified", 3, N, MC, 3, 1);
        check<uniLN5>(N, 1, S, MC, BLDD, B, "LNN6K/3", rep);
#endif
    }
    if ( 1 )
    {
        std::cout << "\nTBSVLT --- triangular band backward solve --- rank " << RANK;
        //check<uniLTB>(N, 1, S, AB, BLDD, B, "blas::tbsv", rep);
        check<uniLT0>(N, 1, S, AB, BLDD, B, "blas_xtbsvLT", rep);
        check<uniLT1>(N, 1, S, AB, BLDD, B, "xtbsvLTN", rep);
        check<uniLT2>(N, 1, S, AB, BLDD, B, "LTNK<KD>", rep);
    }
    if ( DIM == 3 )
    {
        check<uniLT3>(N, 1, S, AB, BLDD, B, "LTN6K", rep);
#if REAL_IS_DOUBLE && USE_SIMD
        check<uniLT4>(N, 1, S, AB, BLDD, B, "LTN6K_SSE", rep);
#endif
#if UNFINISHED_PRECONDITIONNER
        copy_real(N*BLDD, AB, MC);
        modify_matrix<3>(3, N, MC, BLDD);
        check<uniLT5>(N, 1, S, MC, BLDD, B, "LTN6K/3", rep);
#endif
    }
    free_real(B);
    free_real(S);
    free_real(AB);
    free_real(MC);
}

//------------------------------------------------------------------------------
#pragma mark - GETRS

int* pivot = nullptr;

void getrs0(int N, real const* B, int LDB, real* Y)
{
    int info = 0;
    lapack::xgetrs('N', N, 1, B, LDB, pivot, Y, N, &info);
    assert_true(info==0);
}

void getrs1(int N, real const* B, int LDB, real* Y)
{
    lapack_xgetrsN(N, B, LDB, pivot, Y);
}

void getrs2(int N, real const* B, int LDB, real* Y)
{
    alsatian_iso_xgetrsN<1>(N, B, LDB, pivot, Y);
}

void getrs3(int N, real const* B, int LDB, real* Y)
{
    alsatian_xgetrsN(N, (float*)B, LDB, pivot, Y);
}

void getrs4(int N, real const* B, int LDB, real* Y)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot);
    // Solve L*X = Y, overwriting Y with X.
    alsatian_xtrsmLLN1U(N, (float*)B, LDB, Y);
    // Solve U*X = Y, overwriting Y with X.
    alsatian_xtrsmLUN1C(N, (float*)B, LDB, Y);
}

void getrs5(int N, real const* B, int LDB, real* Y)
{
#if USE_SIMD
    alsatian_xgetrsN_SSE(N, (float*)B, LDB, pivot, Y);
#else
    zero_real(N, Y);
#endif
}

void getrs6(int N, real const* B, int LDB, real* Y)
{
#if ( DIM == 3 )
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_3D(N, (float*)B, LDB, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C_3D(N, (float*)B, LDB, Y);
#endif
}

void getrs7(int N, real const* B, int LDB, real* Y)
{
#if USE_SIMD
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_3D_SSE(N, (float*)B, LDB, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C_3D_SSE(N, (float*)B, LDB, Y);
#endif
}

/* initialize a general matrix that would be inversible. */
void setMatrix(size_t N, real* A, size_t LDA)
{
    zero_real(N*N, A);
    nan_spill(A+N*N);
    for ( size_t i = 0; i < N; ++i )
    {
        for ( size_t j = 0; j < N; ++j )
            A[j+LDA*i] = RNG.sreal();
        A[i+LDA*i] = std::sqrt(N); // dominant diagonal term
    }
}

void testGETRS(int N, size_t rep)
{
    std::cout << "\n" << DIM << "D xGETRS " << N << "x" << N << " --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__;

    int info = 0;
    int LDA = ( N + 3 ) & ~3;
    const size_t MULTI = 128;
    size_t BLK = N * LDA + 4;
    
    real * A = new_real(BLK*MULTI);
    real * Y = new_real(N+4);
    real * S = new_real(N+4);
    pivot = new int[N+4];
    
    for ( int i = 0; i < N; ++i )
        S[i] = RNG.sreal();
    
    if ( 1 )
    {
        setMatrix(N, A, LDA);
        real * Ac = A + BLK;
        copy_real(BLK, A, Ac);
        lapack::xgetf2(N, N, Ac, LDA, pivot, &info);
        printf("\n %lu swaps", count_swaps(1, N, pivot));
        if ( info == 0 )
        {
            check<getrs0>(N, 1, S, Ac, LDA, Y, "lapack::xgetrs", rep);
            check<getrs1>(N, 1, S, Ac, LDA, Y, "lapack_xgetrs", rep);
        }
        
        alsatian_xgetf2(N, A, LDA, pivot, &info);
#if REAL_IS_DOUBLE
        convert_to_floats(N*LDA, A, (float*)A);
#endif
        if ( info == 0 )
        {
            check<getrs2>(N, 1, S, A, LDA, Y, "alsa_getrsN<>", rep);
            check<getrs3>(N, 1, S, A, LDA, Y, "alsa_getrs", rep);
            check<getrs4>(N, 1, S, A, LDA, Y, "alsa_getrs_", rep);
            check<getrs5>(N, 1, S, A, LDA, Y, "alsa_getrsSSE", rep);
            if ( DIM == 3 )
            {
                check<getrs6>(N, 1, S, A, LDA, Y, "alsa_getrs_3D", rep);
                check<getrs7>(N, 1, S, A, LDA, Y, "alsa_getrs3SSE", rep);
            }
        }
    }
    if ( 1 )
    {
        std::cout << "\n" << DIM << "D xGETRS " << N << "x" << N << " --- " << MULTI << " matrices";
        for ( size_t u = 0; u < MULTI; ++u )
        {
            real* mat = A + u * BLK;
            setMatrix(N, mat, LDA);
            lapack::xgetf2(N, N, mat, N, pivot, &info);
        }
        multi<getrs1>(N, 1, S, A, LDA, Y, "lapack::xgetrs", rep, MULTI);
        multi<getrs2>(N, 1, S, A, LDA, Y, "lapack_xgetrs", rep, MULTI);
    }
    if ( 1 )
    {
        for ( size_t u = 0; u < MULTI; ++u )
        {
            real* mat = A + u * BLK;
            setMatrix(N, mat, LDA);
            alsatian_xgetf2(N, mat, LDA, pivot, &info);
#if REAL_IS_DOUBLE
            convert_to_floats(N*LDA, mat, (float*)mat);
#endif
        }
        multi<getrs3>(N, 1, S, A, LDA, Y, "alsa_getrs", rep, MULTI);
        multi<getrs4>(N, 1, S, A, LDA, Y, "alsa_getrs_", rep, MULTI);
        multi<getrs5>(N, 1, S, A, LDA, Y, "alsa_getrsSSE", rep, MULTI);
        if ( DIM == 3 )
        {
            multi<getrs6>(N, 1, S, A, LDA, Y, "alsa_getrs_3D", rep, MULTI);
            multi<getrs7>(N, 1, S, A, LDA, Y, "alsa_getrs3SSE", rep, MULTI);
        }
    }
    
    free_real(Y);
    free_real(S);
    free_real(A);
    delete[] pivot;
}


int main(int argc, char* argv[])
{
    int CNT = 37;
    if ( argc > 1 )
        CNT = std::max(1, atoi(argv[1]));
    size_t REP = 1024;
    RNG.seed();
    testISO(CNT, REP);
    testPOTRS(CNT, REP);
    testTBSV(DIM*CNT, REP);
    testGETRS(DIM*CNT, REP);
    printf("\n");
}
