// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include <sys/time.h>

#define DIM 2

#include "real.h"
#include "tictoc.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "assert_macro.h"
#include "xpttrf.h"
#include "cytoblas.h"

/// number of segments:
const size_t NBS = 1240;
const size_t DISP = 16UL;

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testDPTT(size_t cnt)
{
    real * D = new_real(NBS);
    real * U = new_real(NBS);
    real * B = new_real(NBS);
    real * Ds = new_real(NBS);
    real * Us = new_real(NBS);
    real * Bs = new_real(NBS);

    for ( size_t i = 0; i < NBS; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -RNG.preal();
        Bs[i] = RNG.sreal();
    }

    int info;
    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    lapack_xpttrf(NBS, D, U, &info);
    lapack_xptts2(NBS, 1, D, U, B, 1);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        lapack_xptts2(NBS, 1, D, U, B, 1);
    TicToc::toc("   clapack");
    
    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    lapack::xpttrf(NBS, D, U, &info);
    lapack::xptts2(NBS, 1, D, U, B, 1);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        lapack::xptts2(NBS, 1, D, U, B, 1);
    TicToc::toc("    lapack");

    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    italian_xpttrf(NBS, D, U, &info);
    italian_xptts2(NBS, 1, D, U, B, 1);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        italian_xptts2(NBS, 1, D, U, B, 1);
    TicToc::toc("   italian");
    
    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    alsatian_xpttrf(NBS, D, U, &info);
    alsatian_xptts2(NBS, 1, D, U, B, 1);
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian_xptts2(NBS, 1, D, U, B, 1);
    TicToc::toc("  alsatian");

    free_real(D);
    free_real(U);
    free_real(B);
    free_real(Ds);
    free_real(Us);
    free_real(Bs);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testThomas(size_t cnt)
{
    real * D = new_real(NBS);
    real * U = new_real(NBS);
    real * B = new_real(NBS);
    real * Ds = new_real(NBS);
    real * Us = new_real(NBS);
    real * Bs = new_real(NBS);

    for ( size_t i = 0; i < NBS; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -RNG.preal();
        Bs[i] = RNG.sreal();
    }
    
    int info = 0;
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        lapack::xpttrf(NBS, D, U, &info);
        lapack::xptts2(NBS, 1, D, U, B, 1);
    }
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::toc("    lapack");

    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        //italian_xpttrf(NBS, D, U, &info);
        //italian_xptts2(NBS, 1, D, U, B, 1);
        italian_thomas(NBS, U, D, U, B);
    }
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::toc("   italian");
    
    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        alsatian_thomas(NBS, D, U, B);
    }
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::toc("  alsatian");

    TicToc::tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        tridiagonal_solve(NBS, U, D, U, B);
    }
    VecPrint::print(std::clog, std::min(DISP,NBS), B, 3);
    TicToc::toc("  tridiag.");

    free_real(D);
    free_real(U);
    free_real(B);
    free_real(Ds);
    free_real(Us);
    free_real(Bs);
}

int main(int argc, char* argv[])
{
    RNG.seed();
    std::cout << "testPTTRS  --- real " << sizeof(real) << " --- " << __VERSION__ << "\n";

    testThomas(1<<14);
    //testDPTT(1<<17);
    
    return EXIT_SUCCESS;
}
