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
const size_t NSEG = 121;
const size_t NVAL = DIM * ( NSEG + 1 );


template < void (*FUNC)(int, real const*, real*) >
void check(int N, real const* AB, real* B, char const str[], size_t cnt)
{
    FUNC(N, AB, B);
    VecPrint::print(std::clog, std::min(16,N), B, 3);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        FUNC(N, AB, B);
    printf(" %-14s %7.3f\n", str, toc(cnt));
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

    real * AB = new_real(NSEG*LDAB);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL);

    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    zero_real(NSEG*LDAB, AB);
    for ( size_t i = 0; i < NSEG; ++i )
    {
        AB[  LDAB*i] = 5.0;
        AB[1+LDAB*i] = 0.5;
        AB[2+LDAB*i] = RNG.sreal();
    }
    int info;
    alsatian_xpbtf2L<2>(NSEG, AB, LDAB, &info);
    
#if ( 0 )
    copy_real(NVAL, S, B);
    check<iso0>(AB, B, "buggy BLAS", cnt);
#endif
    
    copy_real(NVAL, S, B);
    check<iso1>(NSEG, AB, B, "blas_pbtrsL", cnt);
    
    copy_real(NVAL, S, B);
    check<iso2>(NSEG, AB, B, "alsa_pbtrsL<D>", cnt);

    copy_real(NVAL, S, B);
    check<iso3>(NSEG, AB, B, "alsa_pbtrsL", cnt);

#if 0 && REAL_IS_DOUBLE
    copy_real(NVAL, S, B);
    check<isoLN>(NSEG, AB, B, "tbsvLNN3", cnt);
    
    copy_real(NVAL, S, B);
    check<isoLT>(NSEG, AB, B, "tbsvLTN3", cnt);
#endif

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------

void iso4(int N, real const* AB, real* B)
{
#ifdef __AVX__
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

void iso5(int N, real const* AB, real* B)
{
    iso_xtrsmLLN<DIM,'I'>(N, AB, N, B);
    iso_xtrsmLLT<DIM,'I'>(N, AB, N, B);
}

void iso6(int N, real const* AB, real* B)
{
    alsatian_xpotrsLref<DIM>(N, AB, N, B);
}

void testPOTRS(size_t cnt)
{
    std::cout << "xTBSVLN " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << DIM << "D --- " << __VERSION__ << "\n";

    real * AB = new_real(NSEG*NSEG);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL);

    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    zero_real(NSEG*NSEG, AB);
    for ( size_t i = 0; i < NSEG; ++i )
    {
        AB[i+NSEG*i] = 5.0;
        AB[i+1+NSEG*i] = 0.5;
        AB[i+2+NSEG*i] = RNG.sreal();
    }
    int info;
    alsatian_xpotf2L(NSEG, AB, NSEG, &info);

    copy_real(NVAL, S, B);
    check<iso4>(NSEG, AB, B, "alsa_trsmLLN3<", cnt);
    
    copy_real(NVAL, S, B);
    check<iso5>(NSEG, AB, B, "iso_trsmLLN<D>", cnt);
    
    copy_real(NVAL, S, B);
    check<iso6>(NSEG, AB, B, "alsa_potrsLref", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}

//------------------------------------------------------------------------------

constexpr size_t BAND_LDD = 2*DIM+1;
constexpr size_t BAND_NUD = 2*DIM;

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
    blas_xtbsvLN<'I'>(N, BAND_NUD, AB, BAND_LDD, B);
    alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
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
#if REAL_IS_DOUBLE
    alsatian_xtbsvLNN6SSE(N, AB, BAND_LDD, B);
#else
    zero_real(N, B);
#endif
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

    real * AB = new_real(NVAL*std::max(NVAL, BAND_LDD));
    real * S = new_real(NVAL);
    real * B = new_real(NVAL);

    zero_real(NVAL*BAND_LDD, AB);
    for ( size_t i = 0; i < NVAL; ++i )
    {
        S[i] = RNG.sreal();
        AB[BAND_LDD*i] = 8.0; //diagonal term
        for ( size_t j = 1; j < BAND_LDD; ++j )
            AB[j+BAND_LDD*i] = RNG.sreal() / j;
    }
    int info;
    alsatian_xpbtf2L<2>(NVAL, AB, BAND_LDD, &info);
    
#if 0
    copy_real(NVAL, S, B);
    check<uni0>(NVAL, AB, B, "blas::", cnt);
#endif
    copy_real(NVAL, S, B);
    check<uni1>(NVAL, AB, B, "blas_tbsv", cnt);

    copy_real(NVAL, S, B);
    check<uni2>(NVAL, AB, B, "tbsvLNN", cnt);
    
    copy_real(NVAL, S, B);
    check<uni3>(NVAL, AB, B, "tbsvLNNK<KD>", cnt);
    
    copy_real(NVAL, S, B);
    check<uni4>(NVAL, AB, B, "mixed", cnt);

    std::cout << "xTBSVLN ---\n";
    
    copy_real(NVAL, S, B);
    check<uniLN0>(NVAL, AB, B, "blas_xtbsvLN", cnt);

    copy_real(NVAL, S, B);
    check<uniLN1>(NVAL, AB, B, "tbsvLNN", cnt);

    copy_real(NVAL, S, B);
    check<uniLN2>(NVAL, AB, B, "LNNK<KD>", cnt);

    copy_real(NVAL, S, B);
    check<uniLN4>(NVAL, AB, B, "LNN6SSE", cnt);

    std::cout << "xTBSVLT ---\n";

    copy_real(NVAL, S, B);
    check<uniLT0>(NVAL, AB, B, "blas_tbsvLT", cnt);
    
    copy_real(NVAL, S, B);
    check<uniLT1>(NVAL, AB, B, "tbsvLTN", cnt);

    copy_real(NVAL, S, B);
    check<uniLT2>(NVAL, AB, B, "LTNK<KD>", cnt);
    
    copy_real(NVAL, S, B);
    check<uniLT3>(NVAL, AB, B, "LTN6", cnt);
    
    copy_real(NVAL, S, B);
    check<uniLT4>(NVAL, AB, B, "LTN6SSE", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testISO(1<<16);
    testPOTRS(1<<14);
    test(1<<16);
}
