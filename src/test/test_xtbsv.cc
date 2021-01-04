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
#include "xtbsv.h"
#include "cytoblas.h"

using namespace TicToc;

#define DIM 3

/// number of segments:
const size_t NSEG = 124;
const size_t NVAL = NSEG*DIM;


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

void iso4(int N, real const* AB, real* B)
{
#ifdef __AVX__
#if ( DIM == 1 )
    alsatian_xtbsvLNN1(N, AB, LDAB, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLNN2(N, AB, LDAB, B);
#elif ( DIM == 3 )
    alsatian_xtbsvLNN3(N, AB, LDAB, B);
#endif
#endif
}

void iso5(int N, real const* AB, real* B)
{
#ifdef __AVX__
#if ( DIM == 1 )
    alsatian_xtbsvLTN1(N, AB, LDAB, B);
#elif ( DIM == 2 )
    alsatian_xtbsvLTN2(N, AB, LDAB, B);
#elif ( DIM == 3 )
    alsatian_xtbsvLTN3(N, AB, LDAB, B);
#endif
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testISO(size_t cnt)
{
    std::cout << "isoXTBSV " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << DIM << "D --- " << __VERSION__ << "\n";

    real * AB = new_real(NSEG*LDAB);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL);

    for ( size_t i = 0; i < NSEG; ++i )
    {
        AB[  LDAB*i] = 5.0;
        AB[1+LDAB*i] = 0.5;
        AB[2+LDAB*i] = RNG.sreal();
    }
    for ( size_t i = 0; i < NVAL; ++i )
        S[i] = RNG.sreal();
    int info;
    alsatian_xpbtf2L<2>(NSEG, AB, LDAB, &info);
    
#if ( 0 )
    copy_real(NVAL, S, B);
    check<iso0>(AB, B, "buggy BLAS", cnt);
#endif
    
    copy_real(NVAL, S, B);
    check<iso1>(NSEG, AB, B, "xtbsv<'I'>", cnt);
    
    copy_real(NVAL, S, B);
    check<iso2>(NSEG, AB, B, "xtbsv<DIM>", cnt);

    copy_real(NVAL, S, B);
    check<iso3>(NSEG, AB, B, "xtbsvLNN3", cnt);
    
#if 0
    copy_real(NVAL, S, B);
    check<iso4>(NSEG, AB, B, "xtbsvLNN3", cnt);
    
    copy_real(NVAL, S, B);
    check<iso5>(NSEG, AB, B, "xtbsvLTN3", cnt);
#endif
    
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
    blas_xtbsvLN<'I'>(N, BAND_NUD, AB, BAND_LDD, B, 1);
    blas_xtbsvLT<'I'>(N, BAND_NUD, AB, BAND_LDD, B, 1);

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
    alsatian_xtbsvLNNK<BAND_NUD>(N, AB, BAND_LDD, B);
    //alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
}

void uni5(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN6(N, AB, BAND_LDD, B);
    //alsatian_xtbsvLTN6(N, AB, BAND_LDD, B);
}

void uni6(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNNK<BAND_NUD>(N, AB, BAND_LDD, B);
    alsatian_xtbsvLTNK<BAND_NUD>(N, AB, BAND_LDD, B);
}

void uni7(int N, real const* AB, real* B)
{
    //alsatian_xtbsvLNN6(N, AB, BAND_LDD, B);
    alsatian_xtbsvLTN6(N, AB, BAND_LDD, B);
}

void uni8(int N, real const* AB, real* B)
{
#if REAL_IS_DOUBLE
    alsatian_xtbsvLTN6SSE(N, AB, BAND_LDD, B);
#endif
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void test(size_t cnt)
{
    std::cout << "XTBSV " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << DIM << "D --- " << __VERSION__ << "\n";

    real * AB = new_real(NVAL*BAND_LDD);
    real * S = new_real(NVAL);
    real * B = new_real(NVAL);
    
    zero_real(NVAL*BAND_LDD, AB);
    for ( size_t i = 0; i < NVAL; ++i )
    {
        S[i] = RNG.sreal();
        AB[  BAND_LDD*i] = 8.0;
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
    check<uni1>(NVAL, AB, B, "blas_xtbsv", cnt);

    copy_real(NVAL, S, B);
    check<uni2>(NVAL, AB, B, "xtbsvLNN", cnt);
    
    copy_real(NVAL, S, B);
    check<uni3>(NVAL, AB, B, "xtbsvLNNK<KD>", cnt);
    
    std::cout << " ---\n";
    
    copy_real(NVAL, S, B);
    check<uni4>(NVAL, AB, B, "LNNK<KD>", cnt);
    
    copy_real(NVAL, S, B);
    check<uni5>(NVAL, AB, B, "LNN6", cnt);

    copy_real(NVAL, S, B);
    check<uni6>(NVAL, AB, B, "LTNK<KD>", cnt);
    
    copy_real(NVAL, S, B);
    check<uni7>(NVAL, AB, B, "LTN6", cnt);
    
    copy_real(NVAL, S, B);
    check<uni8>(NVAL, AB, B, "LTN6_SSE", cnt);

    free_real(B);
    free_real(S);
    free_real(AB);
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testISO(1<<16);
    test(1<<16);
}
