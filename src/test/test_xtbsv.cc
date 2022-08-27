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


/// rank of diagonal matrices:
const int RANK = 2*DIM;
const int BLDD = 2*DIM+1;


/// print only 16 scalars from given vector
inline void print(unsigned num, real const* vec)
{
    if ( num > 12 )
    {
        VecPrint::print(6, vec, 3);
        fprintf(stdout, "...");
        VecPrint::print(6, vec+num-8, 3);
        fprintf(stdout, " |");
        VecPrint::print(2, vec+num, 1);
    }
    else
    {
        VecPrint::print(num, vec, 3);
        fprintf(stdout, " |");
        VecPrint::print(2, vec+num, 1);
    }
    printf("  nrm %9.7f ", blas::nrm2(num, vec));
}

void nan_spill(real * dst)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( size_t i = 0; i < 4; ++i )
        dst[i] = n;
}

template < void (*FUNC)(int, real const*, real*) >
void check(int N, int ORD, real const* S, real const* AB, real* B, char const str[], size_t cnt, size_t sub=128)
{
    printf("\n");
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    FUNC(N, AB, B);
    print(ORD*N, B);
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < sub; ++u )
            FUNC(N, AB, B);
    }
    printf(" %-14s cpu %7.0f", str, tock());
}

template < void (*FUNC)(int, real const*, real*) >
void multi(int N, int ORD, real const* S, real const* AB, real* B, char const str[], size_t cnt, size_t sub)
{
    size_t BLK = N * N + 4;
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < sub; ++u )
            FUNC(N, AB+u*BLK, B);
    }
    printf("\n uncached %-14s cpu %7.0f", str, tock());
}

//------------------------------------------------------------------------------
#pragma mark -

const size_t KD = 2;
const size_t LDAB = 3;

void iso0(int N, real const* AB, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', N, KD, AB, LDAB, B+d, DIM);
        blas::xtbsv('L', 'T', 'N', N, KD, AB, LDAB, B+d, DIM);
    }
}

void iso1(int N, real const* AB, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas_xtbsvLN<'I'>(N, KD, AB, LDAB, B+d, DIM);
        blas_xtbsvLT<'I'>(N, KD, AB, LDAB, B+d, DIM);
    }
}

void iso2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN<DIM>(N, KD, AB, LDAB, B);
    alsatian_xtbsvLTN<DIM>(N, KD, AB, LDAB, B);
}

void iso3(int N, real const* AB, real* B)
{
    U::alsatian_xpbtrsL<DIM>(N, AB, LDAB, B);
}

void iso4(int N, real const* AB, real* B)
{
    alsatian_xpbtrsL<DIM>(N, AB, LDAB, B);
}

void iso5(int N, real const* AB, real* B)
{
#if ( DIM == 3 ) && defined(__AVX__) && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_AVX(N, AB, LDAB, B);
    alsatian_xtbsvLTN2K3_AVX(N, AB, LDAB, B);
#elif ( DIM == 3 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_SIMD(N, AB, LDAB, B);
    alsatian_xtbsvLTN2K3_SIMD(N, AB, LDAB, B);
#elif ( DIM == 3 ) && defined(__SSE3__) && !REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_SSE(N, AB, LDAB, B);
    alsatian_xtbsvLTN2K3_SSE(N, AB, LDAB, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K2_SIMD(N, AB, LDAB, B);
    alsatian_xtbsvLTN2K2_SIMD(N, AB, LDAB, B);
#elif ( DIM == 1 )
    alsatian_xtbsvLNN2K1(N, AB, LDAB, B);
    alsatian_xtbsvLTN2K1(N, AB, LDAB, B);
#else
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDAB, B);
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDAB, B);
#endif
}

void isoLNN(int N, real const* AB, real* B)
{
#if ( DIM == 1 )
    alsatian_xtbsvLNN2K1(N, AB, LDAB, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K2_SIMD(N, AB, LDAB, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLNN2K2(N, AB, LDAB, B);
#elif ( DIM == 3 ) && defined(__AVX__)
    alsatian_xtbsvLNN2K3_AVX(N, AB, LDAB, B);
#elif ( DIM == 3 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLNN2K3_SIMD(N, AB, LDAB, B);
#else
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDAB, B);
#endif
}

void isoLTN(int N, real const* AB, real* B)
{
#if ( DIM == 1 )
    alsatian_xtbsvLTN2K1(N, AB, LDAB, B);
#elif ( DIM == 2 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLTN2K2_SIMD(N, AB, LDAB, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLTN2K2(N, AB, LDAB, B);
#elif ( DIM == 3 ) && defined(__AVX__)
    alsatian_xtbsvLTN2K3_AVX(N, AB, LDAB, B);
#elif ( DIM == 3 ) && USE_SIMD && REAL_IS_DOUBLE
    alsatian_xtbsvLTN2K3_SIMD(N, AB, LDAB, B);
#else
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDAB, B);
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
    check<iso1>(N, DIM, S, AB, B, "blas_pbtrsL", cnt);
    check<iso2>(N, DIM, S, AB, B, "alsa_pbtrsL<D>", cnt);
    check<iso3>(N, DIM, S, AB, B, "alsa_pbtrs_U", cnt);
    check<iso4>(N, DIM, S, AB, B, "alsa_pbtrs", cnt);
    check<iso5>(N, DIM, S, AB, B, "alsa_pbtrs_SSE", cnt);

#if 0 && REAL_IS_DOUBLE
    check<isoLNN>(N, DIM, S, AB, B, "tbsvLNN3", cnt);
    check<isoLTN>(N, DIM, S, AB, B, "tbsvLTN3", cnt);
#endif

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - Cholesky factorization

void pot1(int N, real const* AB, real* B)
{
    alsatian_xpotrsLref<DIM>(N, AB, N, B);
}

void pot2(int N, real const* AB, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, AB, N, B);
    iso_xtrsmLLT<DIM,'I'>(N, AB, N, B);
}

#if defined(__AVX__) && REAL_IS_DOUBLE
void pot3(int N, real const* AB, real* B)
{
#if ( DIM == 3 )
    alsatian_xtrsmLLN3<'I'>(N, AB, N, B);
    alsatian_xtrsmLLT3<'I'>(N, AB, N, B);
#elif ( DIM == 2 )
    alsatian_xtrsmLLN2<'I'>(N, AB, N, B);
    alsatian_xtrsmLLT2<'I'>(N, AB, N, B);
#elif ( DIM == 1 )
    alsatian_xtrsmLLN1<'I'>(N, AB, N, B);
    alsatian_xtrsmLLT1<'I'>(N, AB, N, B);
#endif
}
#endif

void testPOTRS(int N, size_t cnt)
{
    std::cout << "\n Cholesky factorization of full matrix & iso solve " << N << " points --- real " << sizeof(real);
    std::cout << " --- " << __VERSION__;

    real * AB = new_real(N*N+4);
    real * S = new_real(N*DIM);
    real * B = new_real(N*DIM+4);

    for ( size_t i = 0; i < N*DIM; ++i )
        S[i] = RNG.sreal();
    zero_real(N*N, AB);
    nan_spill(AB+N*N);
    for ( size_t i = 0; i < N; ++i )
    {
        real r = 0.0625 * RNG.sreal();
        AB[i+N*i] = 1.5;
        AB[i+1+N*i] = -0.125 + r;
        AB[i+2+N*i] = -0.125 - r;
    }
    int info;
    alsatian_xpotf2L(N, AB, N, &info);

    check<pot1>(N, DIM, S, AB, B, "alsa_potrsLref", cnt);
    check<pot2>(N, DIM, S, AB, B, "iso_trsmLLN<D>", cnt);
#if defined(__AVX__) && REAL_IS_DOUBLE
    check<pot3>(N, DIM, S, AB, B, "alsa_trsmLLND", cnt);
#endif
    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark -

// this gives wrong results as the diagonal terms are not inverted
void uni0(int N, real const* AB, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, RANK, AB, BLDD, B, 1);
    blas::xtbsv('L', 'T', 'N', N, RANK, AB, BLDD, B, 1);
}

void uni1(int N, real const* AB, real* B)
{
    blas_xtbsvLN<'I'>(N, RANK, AB, BLDD, B);
    blas_xtbsvLT<'I'>(N, RANK, AB, BLDD, B);
}

void uni2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN(N, RANK, AB, BLDD, B);
    alsatian_xtbsvLTN(N, RANK, AB, BLDD, B);
}

void uni3(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNNK<RANK>(N, AB, BLDD, B);
    alsatian_xtbsvLTNK<RANK>(N, AB, BLDD, B);
}

void uni4(int N, real const* AB, real* B)
{
    alsatian_xpbtrsLK<RANK>(N, AB, BLDD, B);
}

// this gives wrong results
void uniLNB(int N, real const* AB, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, RANK, AB, BLDD, B, 1);
    //blas::xtbsv('L', 'T', 'N', N, RANK, AB, BLDD, B);
}

void uniLN0(int N, real const* AB, real* B)
{
    blas_xtbsvLN<'I'>(N, RANK, AB, BLDD, B);
    //blas_xtbsvLT<'I'>(N, RANK, AB, BLDD, B);
}

void uniLN1(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN(N, RANK, AB, BLDD, B);
    //alsatian_xtbsvLTN(N, RANK, AB, BLDD, B);
}

void uniLN2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNNK<RANK>(N, AB, BLDD, B);
    //alsatian_xtbsvLTNK<RANK>(N, AB, BLDD, B);
}

void uniLN3(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN6K(N, AB, BLDD, B);
}

#if REAL_IS_DOUBLE && USE_SIMD
void uniLN4(int N, real const* AB, real* B)
{
    U::alsatian_xtbsvLNN6K_SSE(N, AB, BLDD, B);
}

void uniLN5(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN6K_SSE(N, AB, BLDD, B);
}
#endif

// this gives wrong results
void uniLTB(int N, real const* AB, real* B)
{
    //blas::xtbsv('L', 'N', 'N', N, RANK, AB, BLDD, B);
    blas::xtbsv('L', 'T', 'N', N, RANK, AB, BLDD, B, 1);
}

void uniLT0(int N, real const* AB, real* B)
{
    //blas_xtbsvLN<'I'>(N, RANK, AB, BLDD, B);
    blas_xtbsvLT<'I'>(N, RANK, AB, BLDD, B);
}

void uniLT1(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN(N, RANK, AB, BLDD, B);
    alsatian_xtbsvLTN(N, RANK, AB, BLDD, B);
}

void uniLT2(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNNK<RANK>(N, AB, BLDD, B);
    alsatian_xtbsvLTNK<RANK>(N, AB, BLDD, B);
}

void uniLT3(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN6K(N, AB, BLDD, B);
    alsatian_xtbsvLTN6K(N, AB, BLDD, B);
}

#if REAL_IS_DOUBLE && USE_SIMD
void uniLT4(int N, real const* AB, real* B)
{
    U::alsatian_xtbsvLTN6K_SSE(N, AB, BLDD, B);
}

void uniLT5(int N, real const* AB, real* B)
{
    alsatian_xtbsvLTN6K_SSE(N, AB, BLDD, B);
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
    check<uni0>(N, 1, S, AB, B, "fail blas:tbsv", cnt);
    check<uni1>(N, 1, S, AB, B, "blas_tbsv", cnt);
    check<uni2>(N, 1, S, AB, B, "tbsvLxN", cnt);
    check<uni3>(N, 1, S, AB, B, "tbsvLxNK<KD>", cnt);
    check<uni4>(N, 1, S, AB, B, "alsa_xpbtrsLK", cnt);
#endif
    if ( 1 )
    {
        std::cout << "\n rank" << RANK << " TBSVLN --- triangular band matrix forward solve";
        //check<uniLNB>(N, S, AB, B, "blas::xtbsv", cnt);
        check<uniLN0>(N, 1, S, AB, B, "blas_xtbsvLN", cnt);
        check<uniLN1>(N, 1, S, AB, B, "xtbsvLNN", cnt);
        check<uniLN2>(N, 1, S, AB, B, "LNNK<KD>", cnt);
        check<uniLN3>(N, 1, S, AB, B, "LLN6", cnt);
#if REAL_IS_DOUBLE && USE_SIMD
        check<uniLN4>(N, 1, S, AB, B, "LNN6SSE_U", cnt);
        check<uniLN5>(N, 1, S, AB, B, "LNN6SSE", cnt);
#endif
    }
    if ( 1 )
    {
        std::cout << "\n rank" << RANK << " TBSVLT --- triangular band matrix backward solve";
        //check<uniLTB>(N, S, AB, B, "blas::tbsv", cnt);
        check<uniLT0>(N, 1, S, AB, B, "blas_xtbsvLT", cnt);
        check<uniLT1>(N, 1, S, AB, B, "xtbsvLTN", cnt);
        check<uniLT2>(N, 1, S, AB, B, "LTNK<KD>", cnt);
        check<uniLT3>(N, 1, S, AB, B, "LTN6", cnt);
#if REAL_IS_DOUBLE && USE_SIMD
        check<uniLT4>(N, 1, S, AB, B, "LTN6SSE_U", cnt);
        check<uniLT5>(N, 1, S, AB, B, "LTN6SSE", cnt);
#endif
    }
    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - GETRS

int* pivot = nullptr;

void getrs1(int N, real const* B, real* Y)
{
    int info = 0;
    lapack::xgetrs('N', N, 1, B, N, pivot, Y, N, &info);
    assert_true(info==0);
}

void getrs2(int N, real const* B, real* Y)
{
    lapack_xgetrsN(N, B, N, pivot, Y);
}

void getrs3(int N, real const* B, real* Y)
{
    //alsatian_xgetrsN(N, B, N, pivot, Y);
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U(N, (float*)B, N, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1I(N, (float*)B, N, Y);
}

void getrs4(int N, real const* B, real* Y)
{
#if USE_SIMD
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_SSE(N, (float*)B, N, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1I_SSE(N, (float*)B, N, Y);
#else
    zero_real(N, Y);
#endif
}

void getrs5(int N, real const* B, real* Y)
{
#if USE_SIMD
    alsatian_xgetrsN_SSE(N, B, N, pivot, Y);
#else
    zero_real(N, Y);
#endif
}


/// convert doubles to floats
void convert_to_floats(size_t cnt, double const* src, float* dst)
{
    #pragma ivdep
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = (float)src[i];
}

void setMatrix(size_t N, real* A)
{
    zero_real(N*N, A);
    nan_spill(A+N*N);
    for ( size_t i = 0; i < N; ++i )
    {
        for ( size_t j = 0; j < N; ++j )
            A[j+N*i] = 0.01 * RNG.sreal();
        A[i+N*i] = 1.0; //dominant diagonal term
    }
}

void testGETRS(int N, size_t cnt)
{
    std::cout << "\n" << DIM << "D xGETRS " << N << "x" << N << " --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__;

    int info = 0;
    const size_t MULTI = 128;
    size_t BLK = N * N + 4;
    
    real * A = new_real(BLK*MULTI);
    real * M = A + BLK;
    real * Y = new_real(N+4);
    real * S = new_real(N);
    pivot = new int[N+4];
    
    for ( int i = 0; i < N; ++i )
        S[i] = RNG.sreal();
    
    setMatrix(N, A);
    copy_real(BLK, A, M);
    lapack::xgetf2(N, N, M, N, pivot, &info);
    if ( info == 0 )
    {
        check<getrs1>(N, 1, S, M, Y, "lapack::xgetrs", cnt);
        check<getrs2>(N, 1, S, M, Y, "lapack_xgetrsN", cnt);
    }
    alsatian_xgetf2(N, A, N, pivot, &info);
#if REAL_IS_DOUBLE
    convert_to_floats(N*N, A, (float*)A);
#endif
    if ( info == 0 )
    {
        check<getrs3>(N, 1, S, A, Y, "alsa_getrsN", cnt);
        check<getrs4>(N, 1, S, A, Y, "xtrsmLLN1U_SSE", cnt);
        check<getrs5>(N, 1, S, A, Y, "alsa_getrsNSSE", cnt);
    }
    
    std::cout << "\n" << DIM << "D xGETRS--- " << MULTI << " matrices";
    for ( size_t u = 0; u < MULTI; ++u )
    {
        real* mat = A + u * BLK;
        setMatrix(N, mat);
        lapack::xgetf2(N, N, mat, N, pivot, &info);
    }
    multi<getrs1>(N, 1, S, A, Y, "lapack::xgetrs", cnt, MULTI);
    multi<getrs2>(N, 1, S, A, Y, "lapack_xgetrsN", cnt, MULTI);

    for ( size_t u = 0; u < MULTI; ++u )
    {
        real* mat = A + u * BLK;
        setMatrix(N, mat);
        alsatian_xgetf2(N, mat, N, pivot, &info);
#if REAL_IS_DOUBLE
        convert_to_floats(N*N, mat, (float*)mat);
#endif
    }
    multi<getrs3>(N, 1, S, A, Y, "alsa_getrsN", cnt, MULTI);
    multi<getrs5>(N, 1, S, A, Y, "alsa_getrsNSSE", cnt, MULTI);

    free_real(Y);
    free_real(S);
    free_real(A);
    delete[] pivot;
}


int main(int argc, char* argv[])
{
    int CNT = 113;
    if ( argc > 1 )
        CNT = std::max(1, atoi(argv[1]));
    size_t REP = 1<<10;
    RNG.seed();
    testISO(CNT, REP);
    testPOTRS(CNT, REP);
    testTBSV(CNT, REP);
    testGETRS(CNT, REP);
    printf("\n");
}
