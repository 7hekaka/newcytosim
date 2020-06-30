// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN 30.06.2020

#include <sys/time.h>

#include "real.h"
#include "tictoc.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "assert_macro.h"
#include "xtbsv.h"
#include "cytoblas.h"

#define DIM 3

const size_t KD = 2;
const size_t LDAB = 3;

/// number of segments:
const size_t NBS = 1240;
const size_t DISP = 16UL;


inline void alsatian0(int N, real const* AB, real* B)
{
    for ( int d = 0; d < DIM; ++d )
    {
        blas_xtbsvLN<'I'>(N, KD, AB, LDAB, B+d, DIM);
        blas_xtbsvLT<'I'>(N, KD, AB, LDAB, B+d, DIM);
    }
}

inline void alsatian1(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN<DIM>(N, 2, AB, LDAB, B);
    alsatian_xtbsvLTN<DIM>(N, 2, AB, LDAB, B);
}

inline void alsatian2(int N, real const* AB, real* B)
{
    alsatian_xtbsvLNN3(N, AB, LDAB, B);
    alsatian_xtbsvLTN3(N, AB, LDAB, B);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testTBSV(size_t cnt)
{
    
    std::cout << "testDPTT " << __VERSION__ << "\n";

    real * AB  = new_real(NBS*LDAB);
    real * Bs = new_real(NBS);
    real * B  = new_real(NBS);

    for ( size_t i = 0; i < NBS; ++i )
    {
        AB[  LDAB*i] =  5.0;
        AB[1+LDAB*i] = -0.5;
        AB[2+LDAB*i] =  RNG.sreal();
        Bs[i] = RNG.sreal();
    }
    int info;
    alsatian_xpbtf2L<2>(NBS, AB, LDAB, &info);

    copy_real(NBS, Bs, B);
    alsatian0(NBS, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 2);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian0(NBS, AB, B);
    TicToc::toc("    xtbsv<'I'>");
    
    copy_real(NBS, Bs, B);
    alsatian1(NBS, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 2);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian1(NBS, AB, B);
    TicToc::toc("    xtbsv<DIM>");

    copy_real(NBS, Bs, B);
    alsatian2(NBS, AB, B);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 2);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian2(NBS, AB, B);
    TicToc::toc("    xtbsvLNN3 ");

    free_real(B);
    free_real(Bs);
    free_real(AB);
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testTBSV(1<<16);
    return EXIT_SUCCESS;
}
