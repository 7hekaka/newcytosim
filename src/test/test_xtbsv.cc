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


/// number of segments:
constexpr size_t NPTS = 113;
constexpr size_t NVAL = DIM * NPTS;

constexpr size_t RANK = 2*DIM;
constexpr size_t BLDD = 2*DIM+1;


/// print only 16 scalars from given vector
inline void print(size_t num, real const* vec)
{
    if ( num > 16 )
    {
        VecPrint::print(8, vec, 3);
        fprintf(stdout, "...");
        VecPrint::print(8, vec+num-8, 3);
        fprintf(stdout, " |");
        VecPrint::print(2, vec+num, 1);
    }
    else
    {
        VecPrint::print(num, vec, 3);
        fprintf(stdout, " |");
        VecPrint::print(2, vec+num, 1);
    }
    printf("  nrm %7.3f ", blas::nrm8(num, vec));
}

void nan_spill(real * dst)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( size_t i = 0; i < 4; ++i )
        dst[i] = n;
}

template < void (*FUNC)(int, real const*, real*) >
void check(int N, int ORD, real const* S, real const* AB, real* B, char const str[], size_t cnt)
{
    copy_real(ORD*N, S, B);
    nan_spill(B+ORD*N);
    FUNC(N, AB, B);
    print(ORD*N, B);
    tic();
    const size_t SUB = 128;
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(ORD*N, S, B);
        for ( size_t u = 0; u < SUB; ++u )
            FUNC(N, AB, B);
    }
    printf(" %-14s cpu %7.0f\n", str, toc(cnt*SUB));
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

void isoLN(int N, real const* AB, real* B)
{
#ifdef __AVX__
#if ( DIM == 1 )
    alsatian_xtbsvLNN1(N, AB, LDAB, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLNN2(N, AB, LDAB, B);
#elif ( DIM == 3 )
    alsatian_xtbsvLNN3(N, AB, LDAB, B);
#endif
#else
    zero_real(N, B);
#endif
}

void isoLT(int N, real const* AB, real* B)
{
#ifdef __AVX__
#if ( DIM == 1 )
    alsatian_xtbsvLTN1(N, AB, LDAB, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLTN2(N, AB, LDAB, B);
#elif ( DIM == 3 )
    alsatian_xtbsvLTN3(N, AB, LDAB, B);
#endif
#else
    zero_real(N, B);
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testISO(size_t cnt)
{
    std::cout << DIM << "D isoTBSV " << NPTS << " points --- real " << sizeof(real);
    std::cout << " --- " << __VERSION__ << "\n";

    real * AB = new_real(NPTS*LDAB+4);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL+4);

    zero_real(NPTS*LDAB, AB);
    nan_spill(AB+NPTS*LDAB);
    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = 1.0 + RNG.sreal();
    for ( size_t i = 0; i < NPTS; ++i )
    {
        real r = 0.0625 * RNG.sreal();
        AB[  LDAB*i] = 1.5;
        AB[1+LDAB*i] = -0.125 + r;
        AB[2+LDAB*i] = -0.125 - r;
    }
    int info;
    alsatian_xpbtf2L<2>(NPTS, AB, LDAB, &info);
    
    //check<iso0>(NPTS, DIM, S, AB, B, "buggy BLAS", cnt);
    check<iso1>(NPTS, DIM, S, AB, B, "blas_pbtrsL", cnt);
    check<iso2>(NPTS, DIM, S, AB, B, "alsa_pbtrsL<D>", cnt);
    check<iso3>(NPTS, DIM, S, AB, B, "alsa_pbtrs_U", cnt);
    check<iso4>(NPTS, DIM, S, AB, B, "alsa_pbtrs", cnt);

#if 0 && REAL_IS_DOUBLE
    check<isoLN>(NPTS, DIM, S, AB, B, "tbsvLNN3", cnt);
    check<isoLT>(NPTS, DIM, S, AB, B, "tbsvLTN3", cnt);
#endif

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark -

void pot1(int N, real const* AB, real* B)
{
    alsatian_xpotrsLref<DIM>(N, AB, N, B);
}

void pot2(int N, real const* AB, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, AB, N, B);
    iso_xtrsmLLT<DIM,'I'>(N, AB, N, B);
}

void pot3(int N, real const* AB, real* B)
{
#if defined(__AVX__) && REAL_IS_DOUBLE
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
#else
    zero_real(DIM*N, B);
#endif
}

void testPOTRS(size_t cnt)
{
    std::cout << DIM << "D xTBSVLN " << NPTS << " points --- real " << sizeof(real);
    std::cout << " --- " << __VERSION__ << "\n";

    real * AB = new_real(NPTS*NPTS+4);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL+4);

    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    zero_real(NPTS*NPTS, AB);
    nan_spill(AB+NPTS*NPTS);
    for ( size_t i = 0; i < NPTS; ++i )
    {
        real r = 0.0625 * RNG.sreal();
        AB[i+NPTS*i] = 1.5;
        AB[i+1+NPTS*i] = -0.125 + r;
        AB[i+2+NPTS*i] = -0.125 - r;
    }
    int info;
    alsatian_xpotf2L(NPTS, AB, NPTS, &info);

    check<pot1>(NPTS, DIM, S, AB, B, "alsa_potrsLref", cnt);
    check<pot2>(NPTS, DIM, S, AB, B, "iso_trsmLLN<D>", cnt);
    check<pot3>(NPTS, DIM, S, AB, B, "alsa_trsmLLND", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark -

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
    alsatian_xtbsvLNN6(N, AB, BLDD, B);
}

void uniLN4(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE && defined(__SSE3__)
    U::alsatian_xtbsvLNN6SSE(N, AB, BLDD, B);
#else
    zero_real(N, B);
#endif
}

void uniLN5(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE && defined(__SSE3__)
    alsatian_xtbsvLNN6SSE(N, AB, BLDD, B);
#else
    zero_real(N, B);
#endif
}


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
    //alsatian_xtbsvLNN6(N, AB, BLDD, B);
    alsatian_xtbsvLTN6(N, AB, BLDD, B);
}

void uniLT4(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE && defined(__SSE3__)
    U::alsatian_xtbsvLTN6SSE(N, AB, BLDD, B);
#else
    zero_real(N, B);
#endif
}

void uniLT5(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE && defined(__SSE3__)
    alsatian_xtbsvLTN6SSE(N, AB, BLDD, B);
#else
    zero_real(N, B);
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testTBSV(size_t cnt)
{
    std::cout << DIM << "D xTBSV " << NPTS << " points --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__ << "\n";

    real * AB = new_real(NVAL*std::max(NVAL, BLDD)+4);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL+4);

    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    zero_real(NVAL*BLDD, AB);
    nan_spill(AB+NPTS*BLDD);
    for ( size_t i = 0; i < NVAL; ++i )
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
    alsatian_xpbtf2L<RANK>(NVAL, AB, BLDD, &info);
    
#if 1
    check<uni0>(NVAL, 1, S, AB, B, "blas::tbsv", cnt);
    check<uni1>(NVAL, 1, S, AB, B, "blas_tbsv", cnt);
    check<uni2>(NVAL, 1, S, AB, B, "tbsvLxN", cnt);
    check<uni3>(NVAL, 1, S, AB, B, "tbsvLxNK<KD>", cnt);
    check<uni4>(NVAL, 1, S, AB, B, "alsa_xpbtrsLK", cnt);
#endif
#if 1
    std::cout << "xTBSVLN ---\n";
    
    //check<uniLNB>(NVAL, S, AB, B, "blas::xtbsv", cnt);
    check<uniLN0>(NVAL, 1, S, AB, B, "blas_xtbsvLN", cnt);
    check<uniLN1>(NVAL, 1, S, AB, B, "xtbsvLNN", cnt);
    check<uniLN2>(NVAL, 1, S, AB, B, "LNNK<KD>", cnt);
    check<uniLN3>(NVAL, 1, S, AB, B, "LLN6", cnt);
    check<uniLN4>(NVAL, 1, S, AB, B, "LNN6SSE_U", cnt);
    check<uniLN5>(NVAL, 1, S, AB, B, "LNN6SSE", cnt);
#endif
    
    std::cout << "xTBSVLT ---\n";

    //check<uniLTB>(NVAL, S, AB, B, "blas::tbsv", cnt);
    check<uniLT0>(NVAL, 1, S, AB, B, "blas_xtbsvLT", cnt);
    check<uniLT1>(NVAL, 1, S, AB, B, "xtbsvLTN", cnt);
    check<uniLT2>(NVAL, 1, S, AB, B, "LTNK<KD>", cnt);
    check<uniLT3>(NVAL, 1, S, AB, B, "LTN6", cnt);
    check<uniLT4>(NVAL, 1, S, AB, B, "LTN6SSE_U", cnt);
    check<uniLT5>(NVAL, 1, S, AB, B, "LTN6SSE", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------
#pragma mark - GETRS

int* pivot = nullptr;

void getrs1(int N, real const* B, real* Y)
{
#if defined(__SSE3__)
    alsatian_xgetrsN_SSE(N, B, N, pivot, Y);
#endif
}

void getrs2(int N, real const* B, real* Y)
{
    //alsatian_xgetrsN(N, B, N, pivot, Y);
    // Apply row interchanges to the right hand side.
    xlaswp1(Y, 1, N, pivot); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U(N, (float*)B, N, Y);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1I(N, (float*)B, N, Y);
}

void getrs3(int N, real const* B, real* Y)
{
    lapack_xgetrsN(N, B, N, pivot, Y);
}

void getrs4(int N, real const* B, real* Y)
{
    int info = 0;
    lapack::xgetrs('N', N, 1, B, N, pivot, Y, N, &info);
    assert_true(info==0);
}

/// convert doubles to floats
void convert_to_floats(size_t cnt, double const* src, float* dst)
{
    #pragma ivdep
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = (float)src[i];
}

void testGETRS(size_t cnt)
{
    std::cout << DIM << "D xGETRS " << NPTS << " points --- real " << sizeof(real);
    std::cout << " --- "  << __VERSION__ << "\n";

    real * A = new_real(NVAL*NVAL+4);
    real * M = new_real(NVAL*NVAL+4);
    real * Y = new_real(NVAL+4);
    real * S = new_real(NVAL);
    pivot = new int[NVAL+4];
    
    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    
    zero_real(NVAL*NVAL, A);
    nan_spill(A+NVAL*NVAL);
    for ( size_t i = 0; i < NVAL; ++i )
    {
        for ( size_t j = 0; j < NVAL; ++j )
            A[j+NVAL*i] = 0.0125 * RNG.sreal();
        A[i+NVAL*i] = 2.0; //diagonal term
    }
    copy_real(NVAL*NVAL, A, M);

    int info = 0;
    lapack::xgetf2(NVAL, NVAL, M, NVAL, pivot, &info);
    if ( info == 0 )
    {
        check<getrs4>(NVAL, 1, S, M, Y, "lapack::xgetrs", cnt);
        check<getrs3>(NVAL, 1, S, M, Y, "lapack_xgetrsN", cnt);
    }
    alsatian_xgetf2(NVAL, A, NVAL, pivot, &info);
#if REAL_IS_DOUBLE
    convert_to_floats(NVAL*NVAL, A, (float*)A);
#endif
    if ( info == 0 )
    {
        check<getrs2>(NVAL, 1, S, A, Y, "alsa_getrsN", cnt);
        check<getrs1>(NVAL, 1, S, A, Y, "alsa_getrsNSSE", cnt);
    }
    free_real(Y);
    free_real(S);
    free_real(A);
    free_real(M);
    delete[] pivot;
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testISO(1<<9);
    testPOTRS(1<<7);
    testTBSV(1<<9);
    testGETRS(1<<9);
}
