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

const size_t KD = 2;
const size_t LDAB = 3;

/// number of segments:
const size_t NSEG = 1240;
const size_t DISP = 16UL;


 void alsatian0(int N, real const* AB, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', N, KD, AB, LDAB, B+d, DIM);
        blas::xtbsv('L', 'T', 'N', N, KD, AB, LDAB, B+d, DIM);
    }
}

 void alsatian1(int N, real const* AB, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas_xtbsvLN<'I'>(N, KD, AB, LDAB, B+d, DIM);
        blas_xtbsvLT<'I'>(N, KD, AB, LDAB, B+d, DIM);
    }
}

 void alsatian2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDAB, B);
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDAB, B);
}

 void alsatian3(int N, real const* AB, real* B)
{
    alsatian_xpbtrsL<DIM>(N, AB, LDAB, B);
}

void alsatian4(int N, real const* AB, real* B)
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

void alsatian5(int N, real const* AB, real* B)
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
void testTBSV(size_t cnt)
{
    std::cout << "testTBSV " << NSEG << " segments --- real " << sizeof(real);
    std::cout << " --- " << __VERSION__ << "\n";

    real * AB  = new_real(NSEG*LDAB);
    real * Bs = new_real(NSEG*DIM);
    real * B  = new_real(NSEG*DIM);

    for ( size_t i = 0; i < NSEG; ++i )
    {
        AB[  LDAB*i] = 5.0;
        AB[1+LDAB*i] = 0.5;
        AB[2+LDAB*i] = RNG.sreal();
    }
    for ( size_t i = 0; i < NSEG*DIM; ++i )
        Bs[i] = RNG.sreal();
    int info;
    alsatian_xpbtf2L<2>(NSEG, AB, LDAB, &info);
    
#if ( 0 )
    copy_real(NSEG*DIM, Bs, B);
    alsatian0(NSEG, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NSEG), B, 2);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian0(NSEG, AB, B);
    printf("    wrong BLAS  %7.3f\n", toc(cnt));
#endif
    
    copy_real(NSEG*DIM, Bs, B);
    alsatian1(NSEG, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NSEG), B, 2);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian1(NSEG, AB, B);
    printf("    xtbsv<'I'>  %7.3f\n", toc(cnt));
    
    copy_real(NSEG*DIM, Bs, B);
    alsatian2(NSEG, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NSEG), B, 2);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian2(NSEG, AB, B);
    printf("    xtbsv<DIM>  %7.3f\n", toc(cnt));

    copy_real(NSEG*DIM, Bs, B);
    alsatian3(NSEG, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NSEG), B, 2);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian3(NSEG, AB, B);
    printf("    xtbsvLNN3   %7.3f\n", toc(cnt));
    
#if 0
    copy_real(NSEG*DIM, Bs, B);
    alsatian4(NSEG, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NSEG), B, 2);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian4(NSEG, AB, B);
    printf("    xtbsvLNN3   %7.3f\n", toc(cnt));
    
    copy_real(NSEG*DIM, Bs, B);
    alsatian5(NSEG, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NSEG), B, 2);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian5(NSEG, AB, B);
    printf("    xtbsvLTN3   %7.3f\n", toc(cnt));
#endif
    
    free_real(B);
    free_real(Bs);
    free_real(AB);
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testTBSV(1<<16);
}
