// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN 30.06.2020

#include <sys/time.h>
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

#define DEVELOP_XTBSV 1
#include "xtbsv.h"

namespace U
{
#undef XTBSV_H
#undef DEVELOP_XTBSV
#define DEVELOP_XTBSV 0
#include "xtbsv.h"
}


/// print 12 scalars from `vec[]` of dimension `num`
inline void print_vector(int num, real const* vec)
{
    if ( num > 12 )
    {
        VecPrint::print(6, vec, 3);
        fprintf(stdout, "...");
        VecPrint::print(6, vec+num-6, 3);
        fprintf(stdout, " |");
        VecPrint::print(2, vec+num, 1);
    }
    else
    {
        VecPrint::print(num, vec, 3);
        fprintf(stdout, " |");
        VecPrint::print(2, vec+num, 1);
    }
    real sum = vec[0];
    for ( int i = 1; i < num; ++i )
        sum += vec[i];
    printf("  sum %+18.16f ", sum);
}

void nan_spill(real * dst)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( size_t i = 0; i < 4; ++i )
        dst[i] = n;
}

template < void (*FUNC)(int, real const*, int, real*) >
void check(int N, int ORD, real const* S, real const* AB, int LDA, real* B, char const str[], size_t cnt, size_t sub=128)
{
    printf("\n");
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    FUNC(N, AB, LDA, B);
    print_vector(ORD*N, B);
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < sub; ++u )
            FUNC(N, AB, LDA, B);
    }
    printf(" %-14s cpu %7.0f", str, tock());
}

template < void (*FUNC)(int, real const*, int, real*) >
void multi(int N, int ORD, real const* S, real const* AB, int LDA, real* B, char const str[], size_t cnt, size_t sub)
{
    size_t BLK = N * N + 4;
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < sub; ++u )
            FUNC(N, AB+u*BLK, LDA, B);
    }
    printf("\n uncached %-14s cpu %7.0f", str, tock());
}

//------------------------------------------------------------------------------
#pragma mark -

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
        blas_xtbsvLN<'I'>(N, KD, AB, LDA, B+d, DIM);
        blas_xtbsvLT<'I'>(N, KD, AB, LDA, B+d, DIM);
    }
}

void iso2(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN<DIM>(N, KD, AB, LDA, B);
    alsatian_xtbsvLTN<DIM>(N, KD, AB, LDA, B);
}

void iso3(int N, real const* AB, int LDA, real* B)
{
    U::alsatian_xpbtrsL<DIM>(N, AB, LDA, B);
}

void iso4(int N, real const* AB, int LDA, real* B)
{
    alsatian_xpbtrsL<DIM>(N, AB, LDA, B);
}

void iso5(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 3 ) && defined(__AVX__) && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_AVX(N, AB, LDA, B);
    alsatian_xtbsvLTN2K3_AVX(N, AB, LDA, B);
#elif ( DIM == 3 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_SIMD(N, AB, LDA, B);
    alsatian_xtbsvLTN2K3_SIMD(N, AB, LDA, B);
#elif ( DIM == 3 ) && defined(__SSE3__) && !REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_SSE(N, AB, LDA, B);
    alsatian_xtbsvLTN2K3_SSE(N, AB, LDA, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K2_SIMD(N, AB, LDA, B);
    alsatian_xtbsvLTN2K2_SIMD(N, AB, LDA, B);
#elif ( DIM == 1 )
    alsatian_xtbsvLNN2K1(N, AB, LDA, B);
    alsatian_xtbsvLTN2K1(N, AB, LDA, B);
#else
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDA, B);
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDA, B);
#endif
}

void isoLNN(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 1 )
    alsatian_xtbsvLNN2K1(N, AB, LDA, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K2_SIMD(N, AB, LDA, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLNN2K2(N, AB, LDA, B);
#elif ( DIM == 3 ) && defined(__AVX__)
    alsatian_xtbsvLNN2K3_AVX(N, AB, LDA, B);
#elif ( DIM == 3 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_SIMD(N, AB, LDA, B);
#else
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDA, B);
#endif
}

void isoLTN(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 1 )
    alsatian_xtbsvLTN2K1(N, AB, LDA, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLTN2K2_SIMD(N, AB, LDA, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLTN2K2(N, AB, LDA, B);
#elif ( DIM == 3 ) && defined(__AVX__)
    alsatian_xtbsvLTN2K3_AVX(N, AB, LDA, B);
#elif ( DIM == 3 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLTN2K3_SIMD(N, AB, LDA, B);
#else
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDA, B);
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testISO(int N, size_t cnt)
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
        S[i] = 1.0 + RNG.sreal();
    for ( int i = 0; i < N; ++i )
    {
        AB[  LDAB*i] = 2 + 0.25 * RNG.sreal();
        AB[1+LDAB*i] = -0.125 + 0.25 * RNG.preal();
        AB[2+LDAB*i] = -0.125 + 0.25 * RNG.preal();
    }
    int info;
    alsatian_xpbtf2L<2>(N, AB, LDAB, &info);
    
    //check<iso0>(NPTS, DIM, S, AB, B, "fail BLAS", cnt);
    check<iso1>(N, DIM, S, AB, LDAB, B, "blas_pbtrsL", cnt);
    check<iso2>(N, DIM, S, AB, LDAB, B, "alsa_pbtrsL<D>", cnt);
    check<iso3>(N, DIM, S, AB, LDAB, B, "alsa_pbtrs_U", cnt);
    check<iso4>(N, DIM, S, AB, LDAB, B, "alsa_pbtrs", cnt);
    check<iso5>(N, DIM, S, AB, LDAB, B, "alsa_pbtrs_SSE", cnt);

#if 0 && REAL_IS_DOUBLE
    check<isoLNN>(N, DIM, S, AB, LDAB, B, "tbsvLNN3", cnt);
    check<isoLTN>(N, DIM, S, AB, LDAB, B, "tbsvLTN3", cnt);
#endif

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - Cholesky factorization

void pot1(int N, real const* AB, int LDA, real* B)
{
    alsatian_xpotrsLref<DIM>(N, AB, LDA, B);
}

void pot2(int N, real const* AB, int LDA, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, AB, LDA, B);
    iso_xtrsmLLT<DIM,'I'>(N, AB, LDA, B);
}

#if defined(__AVX__) && REAL_IS_DOUBLE
void pot3(int N, real const* AB, int LDA, real* B)
{
#if ( DIM == 3 )
    alsatian_xtrsmLLN3<'I'>(N, AB, LDA, B);
    alsatian_xtrsmLLT3<'I'>(N, AB, LDA, B);
#elif ( DIM == 2 )
    alsatian_xtrsmLLN2<'I'>(N, AB, LDA, B);
    alsatian_xtrsmLLT2<'I'>(N, AB, LDA, B);
#elif ( DIM == 1 )
    alsatian_xtrsmLLN1<'I'>(N, AB, LDA, B);
    alsatian_xtrsmLLT1<'I'>(N, AB, LDA, B);
#endif
}
#endif

void testPOTRS(int N, size_t cnt)
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
    alsatian_xpotf2L(N, AB, LDA, &info);
    if ( info == 0 )
    {
        check<pot1>(N, DIM, S, AB, LDA, B, "alsa_potrsLref", cnt);
        check<pot2>(N, DIM, S, AB, LDA, B, "iso_trsmLLN<D>", cnt);
#if defined(__AVX__) && REAL_IS_DOUBLE
        check<pot3>(N, DIM, S, AB, LDA, B, "alsa_trsmLLND", cnt);
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
#pragma mark -
const int RANK = 2*DIM;

// this gives wrong results as the diagonal terms are not inverted
void uni0(int N, real const* AB, int LDA, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, RANK, AB, LDA, B, 1);
    blas::xtbsv('L', 'T', 'N', N, RANK, AB, LDA, B, 1);
}

void uni1(int N, real const* AB, int LDA, real* B)
{
    blas_xtbsvLN<'I'>(N, RANK, AB, LDA, B);
    blas_xtbsvLT<'I'>(N, RANK, AB, LDA, B);
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

void uni4(int N, real const* AB, int LDA, real* B)
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
    blas_xtbsvLN<'I'>(N, RANK, AB, LDA, B);
    //blas_xtbsvLT<'I'>(N, RANK, AB, LDA, B);
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
    U::alsatian_xtbsvLNN6K_SSE(N, AB, LDA, B);
}

void uniLN5(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLNN6K_SSE(N, AB, LDA, B);
}
#endif

// this gives wrong results
void uniLTB(int N, real const* AB, int LDA, real* B)
{
    //blas::xtbsv('L', 'N', 'N', N, RANK, AB, LDA, B);
    blas::xtbsv('L', 'T', 'N', N, RANK, AB, LDA, B, 1);
}

void uniLT0(int N, real const* AB, int LDA, real* B)
{
    //blas_xtbsvLN<'I'>(N, RANK, AB, LDA, B);
    blas_xtbsvLT<'I'>(N, RANK, AB, LDA, B);
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
    U::alsatian_xtbsvLTN6K_SSE(N, AB, LDA, B);
}

void uniLT5(int N, real const* AB, int LDA, real* B)
{
    alsatian_xtbsvLTN6K_SSE(N, AB, LDA, B);
}
#endif

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testTBSV(int N, size_t cnt)
{
    std::cout << "\n rank " << RANK << " banded matrix Cholesky factorization " << N << " points --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__;

    /// rank of diagonal matrices:
    const int BLDD = RANK+1;

    real * AB = new_real(N*std::max(N, BLDD)+4);
    real * S = new_real(N);
    real * B = new_real(N+4);

    for ( int i = 0; i < N; ++i )
        S[i] = RNG.sreal();
    zero_real(N*BLDD, AB);
    nan_spill(AB+N*BLDD);
    for ( int i = 0; i < N; ++i )
    {
        real s = 0, r = 0.0625 * RNG.sreal();
        for ( size_t j = 1; j < BLDD; ++j )
        {
            real x = -0.0625;
            AB[j+BLDD*i] = x;
            s += x;
        }
        AB[  BLDD*i] = 1.0 - 2*s; //diagonal term
        AB[1+BLDD*i] += r;
        AB[2+BLDD*i] -= r;
    }
    int info;
    alsatian_xpbtf2L<RANK>(N, AB, BLDD, &info);
    
#if 1
    check<uni0>(N, 1, S, AB, BLDD, B, "fail blas:tbsv", cnt);
    check<uni1>(N, 1, S, AB, BLDD, B, "blas_tbsv", cnt);
    check<uni2>(N, 1, S, AB, BLDD, B, "tbsvLxN", cnt);
    check<uni3>(N, 1, S, AB, BLDD, B, "tbsvLxNK<KD>", cnt);
    check<uni4>(N, 1, S, AB, BLDD, B, "alsa_xpbtrsLK", cnt);
#endif
    if ( 1 )
    {
        std::cout << "\n rank" << RANK << " TBSVLN --- triangular band matrix forward solve";
        //check<uniLNB>(N, S, AB, B, "blas::xtbsv", cnt);
        check<uniLN0>(N, 1, S, AB, BLDD, B, "blas_xtbsvLN", cnt);
        check<uniLN1>(N, 1, S, AB, BLDD, B, "xtbsvLNN", cnt);
        check<uniLN2>(N, 1, S, AB, BLDD, B, "LNNK<KD>", cnt);
    }
    if ( DIM == 3 )
    {
        check<uniLN3>(N, 1, S, AB, BLDD, B, "LLN6K", cnt);
#if REAL_IS_DOUBLE && USE_SIMD
        check<uniLN4>(N, 1, S, AB, BLDD, B, "U:LNN6K_SSE", cnt);
        check<uniLN5>(N, 1, S, AB, BLDD, B, "LNN6K_SSE", cnt);
#endif
    }
    if ( 1 )
    {
        std::cout << "\n rank" << RANK << " TBSVLT --- triangular band matrix backward solve";
        //check<uniLTB>(N, S, AB, B, "blas::tbsv", cnt);
        check<uniLT0>(N, 1, S, AB, BLDD, B, "blas_xtbsvLT", cnt);
        check<uniLT1>(N, 1, S, AB, BLDD, B, "xtbsvLTN", cnt);
        check<uniLT2>(N, 1, S, AB, BLDD, B, "LTNK<KD>", cnt);
    }
    if ( DIM == 3 )
    {
        check<uniLT3>(N, 1, S, AB, BLDD, B, "LTN6K", cnt);
#if REAL_IS_DOUBLE && USE_SIMD
        check<uniLT4>(N, 1, S, AB, BLDD, B, "U:LTN6K_SSE", cnt);
        check<uniLT5>(N, 1, S, AB, BLDD, B, "LTN6K_SSE", cnt);
#endif
    }
    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - GETRS

int* pivot = nullptr;

void getrs1(int N, real const* B, int LDB, real* Y)
{
    int info = 0;
    lapack::xgetrs('N', N, 1, B, LDB, pivot, Y, N, &info);
    assert_true(info==0);
}

void getrs2(int N, real const* B, int LDB, real* Y)
{
    lapack_xgetrsN(N, B, LDB, pivot, Y);
}

void getrs3(int N, real const* B, int LDB, real* Y)
{
    alsatian_xgetrsN(N, B, LDB, pivot, Y);
}

void getrs4(int N, real const* B, int LDB, real* Y)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U(N, (float*)B, LDB, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C(N, (float*)B, LDB, Y);
}

void getrs5(int N, real const* B, int LDB, real* Y)
{
#if USE_SIMD
    alsatian_xgetrsN_SSE(N, B, LDB, pivot, Y);
#else
    zero_real(N, Y);
#endif
}

void getrs6(int N, real const* B, int LDB, real* Y)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_3D(N, (float*)B, LDB, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C_3D(N, (float*)B, LDB, Y);
}

void getrs7(int N, real const* B, int LDB, real* Y)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_3D_SSE(N, (float*)B, LDB, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C_3D_SSE(N, (float*)B, LDB, Y);
}

/// convert doubles to floats
void convert_to_floats(size_t cnt, double const* src, float* dst)
{
    #pragma omp simd
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = (float)src[i];
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
        A[i+LDA*i] = 1.0 + std::sqrt(N); // dominant diagonal term
    }
}

void testGETRS(int N, size_t cnt)
{
    std::cout << "\n" << DIM << "D xGETRS " << N << "x" << N << " --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__;

    int info = 0;
    int LDA = ( N + 3 ) & ~3;
    const size_t MULTI = 128;
    size_t BLK = N * LDA + 4;
    
    real * A = new_real(BLK*MULTI);
    real * M = A + BLK;
    real * Y = new_real(N+4);
    real * S = new_real(N+4);
    pivot = new int[N+4];
    
    for ( int i = 0; i < N; ++i )
        S[i] = RNG.sreal();
    
    if ( 1 )
    {
        setMatrix(N, A, LDA);
        copy_real(BLK, A, M);
        lapack::xgetf2(N, N, M, LDA, pivot, &info);
        if ( info == 0 )
        {
            check<getrs1>(N, 1, S, M, LDA, Y, "lapack::xgetrs", cnt);
            check<getrs2>(N, 1, S, M, LDA, Y, "lapack_xgetrs", cnt);
        }
        alsatian_xgetf2(N, A, LDA, pivot, &info);
#if REAL_IS_DOUBLE
        convert_to_floats(N*LDA, A, (float*)A);
#endif
        if ( info == 0 )
        {
            check<getrs3>(N, 1, S, A, LDA, Y, "alsa_getrs", cnt);
            check<getrs4>(N, 1, S, A, LDA, Y, "alsa_getrs_", cnt);
            check<getrs5>(N, 1, S, A, LDA, Y, "alsa_getrsSSE", cnt);
            if ( DIM == 3 )
            {
                check<getrs6>(N, 1, S, A, LDA, Y, "alsa_getrs_3D", cnt);
                check<getrs7>(N, 1, S, A, LDA, Y, "alsa_getrs3SSE", cnt);
            }
        }
    }
    if ( 1 )
    {
        std::cout << "\n" << DIM << "D xGETRS--- " << MULTI << " matrices";
        for ( size_t u = 0; u < MULTI; ++u )
        {
            real* mat = A + u * BLK;
            setMatrix(N, mat, LDA);
            lapack::xgetf2(N, N, mat, N, pivot, &info);
        }
        multi<getrs1>(N, 1, S, A, LDA, Y, "lapack::xgetrs", cnt, MULTI);
        multi<getrs2>(N, 1, S, A, LDA, Y, "lapack_xgetrs", cnt, MULTI);
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
        multi<getrs3>(N, 1, S, A, LDA, Y, "alsa_getrs", cnt, MULTI);
        multi<getrs4>(N, 1, S, A, LDA, Y, "alsa_getrs", cnt, MULTI);
        multi<getrs5>(N, 1, S, A, LDA, Y, "alsa_getrsSSE", cnt, MULTI);
        if ( DIM == 3 )
        {
            multi<getrs6>(N, 1, S, A, LDA, Y, "alsa_getrs_3D", cnt, MULTI);
            multi<getrs7>(N, 1, S, A, LDA, Y, "alsa_getrs3SSE", cnt, MULTI);
        }
    }
    
    free_real(Y);
    free_real(S);
    free_real(A);
    delete[] pivot;
}


int main(int argc, char* argv[])
{
    int CNT = 13;
    if ( argc > 1 )
        CNT = std::max(1, atoi(argv[1]));
    size_t REP = 1<<10;
    RNG.seed();
    //testISO(CNT, REP);
    //testPOTRS(CNT, REP);
    //testTBSV(CNT, REP);
    testGETRS(DIM*CNT, REP);
    printf("\n");
}
