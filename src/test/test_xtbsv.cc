// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN 30.06.2020

#include <sys/time.h>

#include "real.h"
#include "tictoc.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "simd.h"
#include "assert_macro.h"
#include "cytoblas.h"
#include "xtbsv.h"
#include "xtrsm.h"

using namespace TicToc;

#define DIM 3

/// number of segments:
constexpr size_t NSEG = 113;
constexpr size_t NVAL = DIM * ( NSEG + 1 );

constexpr size_t BAND_NUD = 2*DIM;
constexpr size_t BAND_LDD = 2*DIM+1;


/// print only 16 scalars from given vector
inline void print(size_t n, real const* vec)
{
    if ( n > 16 )
    {
        VecPrint::print(std::cout, 8, vec, 3);
        fprintf(stdout, " :");
        VecPrint::print(std::cout, 8, vec+n-8, 3);
        fprintf(stdout, " |");
        VecPrint::print(std::cout, 2, vec+n, 1);
    }
    else
    {
        VecPrint::print(std::cout, n, vec, 3);
        fprintf(stdout, " |");
        VecPrint::print(std::cout, 2, vec+n, 1);
    }
}

void nan_spill(size_t N, real * dst)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( size_t i = 0; i < 4; ++i )
        dst[N+i] = n;
}

template < void (*FUNC)(int, real const*, real*) >
void check(int N, real const* S, real const* AB, real* B, char const str[], size_t cnt)
{
    const size_t SUB = 128;
    copy_real(N, S, B);
    nan_spill(NVAL, B);
    FUNC(N, AB, B);
    print(N, B);
    printf("  nrm %7.3f ", blas::nrm8(N, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NVAL, S, B);
        for ( size_t u = 0; u < SUB; ++u )
            FUNC(N, AB, B);
    }
    printf(" %-14s %7.3f\n", str, toc(cnt*SUB));
}

//------------------------------------------------------------------------------

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
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDAB, B);
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDAB, B);
}

void iso3(int N, real const* AB, real* B)
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
    std::cout << "isoTBSV " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << DIM << "D --- " << __VERSION__ << "\n";

    real * AB = new_real(NSEG*LDAB+4);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL+4);

    zero_real(NSEG*LDAB, AB);
    nan_spill(NSEG*LDAB, AB);
    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = 1.0 + RNG.sreal();
    for ( size_t i = 0; i < NSEG; ++i )
    {
        real r = 0.0625 * RNG.sreal();
        AB[  LDAB*i] = 1.5;
        AB[1+LDAB*i] = -0.125 + r;
        AB[2+LDAB*i] = -0.125 - r;
    }
    int info;
    alsatian_xpbtf2L<2>(NSEG, AB, LDAB, &info);
    
    //check<iso0>(NSEG, S, AB, B, "buggy BLAS", cnt);
    check<iso1>(NSEG, S, AB, B, "blas_pbtrsL", cnt);
    check<iso2>(NSEG, S, AB, B, "alsa_pbtrsL<D>", cnt);
    check<iso3>(NSEG, S, AB, B, "alsa_pbtrsL", cnt);

#if 0 && REAL_IS_DOUBLE
    check<isoLN>(NSEG, S, AB, B, "tbsvLNN3", cnt);
    check<isoLT>(NSEG, S, AB, B, "tbsvLTN3", cnt);
#endif

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------

void iso4(int N, real const* AB, real* B)
{
    alsatian_xpotrsLref<DIM>(N, AB, N, B);
}

void iso5(int N, real const* AB, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, AB, N, B);
    iso_xtrsmLLT<DIM,'I'>(N, AB, N, B);
}

void iso6(int N, real const* AB, real* B)
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
    zero_real(N, B);
#endif
}

void testPOTRS(size_t cnt)
{
    std::cout << "xTBSVLN " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << DIM << "D --- " << __VERSION__ << "\n";

    real * AB = new_real(NSEG*NSEG+4);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL+4);

    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    zero_real(NSEG*NSEG, AB);
    nan_spill(NSEG*NSEG, AB);
    for ( size_t i = 0; i < NSEG; ++i )
    {
        real r = 0.0625 * RNG.sreal();
        AB[i+NSEG*i] = 1.5;
        AB[i+1+NSEG*i] = -0.125 + r;
        AB[i+2+NSEG*i] = -0.125 - r;
    }
    int info;
    alsatian_xpotf2L(NSEG, AB, NSEG, &info);

    check<iso4>(NSEG, S, AB, B, "alsa_potrsLref", cnt);
    check<iso5>(NSEG, S, AB, B, "iso_trsmLLN<D>", cnt);
    check<iso6>(NSEG, S, AB, B, "alsa_trsmLLND", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------

void uni0(int N, real const* AB, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, BAND_NUD, AB, BAND_LDD, B, 1);
    blas::xtbsv('L', 'T', 'N', N, BAND_NUD, AB, BAND_LDD, B, 1);
}

void uni1(int N, real const* AB, real* B)
{
    blas_xtbsvLN<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
    blas_xtbsvLT<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
}

void uni2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN(N, BAND_NUD, AB, BAND_LDD, B);
    alsatian_xtbsvLTN(N, BAND_NUD, AB, BAND_LDD, B);
}

void uni3(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNNK<BAND_NUD>(N, AB, BAND_LDD, B);
    alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
}

void uni4(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN(N, BAND_NUD, AB, BAND_LDD, B);
    alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
}

// this gives wrong results
void uniLNB(int N, real const* AB, real* B)
{
    blas::xtbsv('L', 'N', 'N', N, BAND_NUD, AB, BAND_LDD, B, 1);
    //blas::xtbsv('L', 'T', 'N', N, BAND_NUD, AB, BAND_LDD, B);
}

void uniLN0(int N, real const* AB, real* B)
{
    blas_xtbsvLN<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
    //blas_xtbsvLT<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
}

void uniLN1(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN(N, BAND_NUD, AB, BAND_LDD, B);
    //alsatian_xtbsvLTN(N, BAND_NUD, AB, BAND_LDD, B);
}

void uniLN2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNNK<BAND_NUD>(N, AB, BAND_LDD, B);
    //alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
}

void uniLN3(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN6(N, AB, BAND_LDD, B);
    zero_real(N, B);
}

void uniLN4(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE && defined(__SSE__)
    alsatian_xtbsvLNN6SSEone(N, AB, BAND_LDD, B);
#else
    zero_real(N, B);
#endif
}

void uniLN5(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE && defined(__SSE__)
    alsatian_xtbsvLNN6SSE(N, AB, BAND_LDD, B);
#else
    zero_real(N, B);
#endif
}


// this gives wrong results
void uniLTB(int N, real const* AB, real* B)
{
    //blas::xtbsv('L', 'N', 'N', N, BAND_NUD, AB, BAND_LDD, B);
    blas::xtbsv('L', 'T', 'N', N, BAND_NUD, AB, BAND_LDD, B, 1);
}

void uniLT0(int N, real const* AB, real* B)
{
    //blas_xtbsvLN<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
    blas_xtbsvLT<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
}

void uniLT1(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN(N, BAND_NUD, AB, BAND_LDD, B);
    alsatian_xtbsvLTN(N, BAND_NUD, AB, BAND_LDD, B);
}

void uniLT2(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNNK<BAND_NUD>(N, AB, BAND_LDD, B);
    alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
}

void uniLT3(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN6(N, AB, BAND_LDD, B);
    alsatian_xtbsvLTN6(N, AB, BAND_LDD, B);
}

void uniLT4(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE
    alsatian_xtbsvLTN6SSEone(N, AB, BAND_LDD, B);
#else
    zero_real(N, B);
#endif
}

void uniLT5(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE
    alsatian_xtbsvLTN6SSE(N, AB, BAND_LDD, B);
#else
    zero_real(N, B);
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void test(size_t cnt)
{
    std::cout << "xTBSV " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << DIM << "D --- " << __VERSION__ << "\n";

    real * AB = new_real(NVAL*std::max(NVAL, BAND_LDD)+4);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL+4);

    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    zero_real(NVAL*BAND_LDD, AB);
    nan_spill(NSEG*BAND_LDD, AB);
    for ( size_t i = 0; i < NVAL; ++i )
    {
        real s = 0, r = 0.0625 * RNG.sreal();
        for ( size_t j = 1; j < BAND_LDD; ++j )
        {
            real x = -0.0625;
            AB[j+BAND_LDD*i] = x;
            s += x;
        }
        AB[  BAND_LDD*i] = 1.0 - 2*s; //diagonal term
        AB[1+BAND_LDD*i] += r;
        AB[2+BAND_LDD*i] -= r;
    }
    int info;
    alsatian_xpbtf2L<BAND_NUD>(NVAL, AB, BAND_LDD, &info);
    
    //check<uni0>(NVAL, S, AB, B, "blas::", cnt);
    check<uni1>(NVAL, S, AB, B, "blas_tbsv", cnt);
    check<uni2>(NVAL, S, AB, B, "tbsvLxN", cnt);
    check<uni3>(NVAL, S, AB, B, "tbsvLxNK<KD>", cnt);
    check<uni4>(NVAL, S, AB, B, "mixed", cnt);

    std::cout << "xTBSVLN ---\n";
    
    //check<uniLNB>(NVAL, S, AB, B, "blas::xtbsv", cnt);
    check<uniLN0>(NVAL, S, AB, B, "blas_xtbsvLN", cnt);
    check<uniLN1>(NVAL, S, AB, B, "tbsvLNN", cnt);
    check<uniLN2>(NVAL, S, AB, B, "LNNK<KD>", cnt);
    check<uniLN4>(NVAL, S, AB, B, "LNN6SSEone", cnt);
    check<uniLN5>(NVAL, S, AB, B, "LNN6SSE", cnt);

    std::cout << "xTBSVLT ---\n";

    //check<uniLTB>(NVAL, S, AB, B, "blas::tbsv", cnt);
    check<uniLT0>(NVAL, S, AB, B, "blas_tbsvLT", cnt);
    check<uniLT1>(NVAL, S, AB, B, "tbsvLTN", cnt);
    check<uniLT2>(NVAL, S, AB, B, "LTNK<KD>", cnt);
    check<uniLT3>(NVAL, S, AB, B, "LTN6", cnt);
    check<uniLT4>(NVAL, S, AB, B, "LTN6SSEone", cnt);
    check<uniLT5>(NVAL, S, AB, B, "LTN6SSE", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testISO(1<<9);
    testPOTRS(1<<7);
    test(1<<9);
}
