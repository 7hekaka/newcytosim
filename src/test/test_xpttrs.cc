// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include <sys/time.h>

#include "real.h"
#include "timer.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "assert_macro.h"
#include "xpttrf.h"
#include "cytoblas.h"

//------------------------------------------------------------------------------

template < void (*FUNC)(int, real*, real*, int*) >
void test(int N, real const* D, real const* U, real* d, real* u, char const str[], size_t cnt)
{
    int info;
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(N, D, d);
        copy_real(N, U, u);
        FUNC(N, d, u, &info);
    }
    printf(" %12s cpu %5.0f\n", str, tock());
}


/// benchmark division
void divide(int size, real* D, real* E, int*)
{
    for ( int i = 0; i < size-1; ++i )
        E[i] = E[i] / D[i];
}

void italian_factor(int size, real* D, real* E, int* info)
{
    italian_xpttrf(size, D, E, info);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testFactor(int NSEG, size_t cnt)
{
    real * d = new_real(NSEG);
    real * u = new_real(NSEG);
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);

    for ( int i = 0; i < NSEG; ++i )
    {
        D[i] = 2.0;
        U[i] = -0.5 * RNG.preal();
    }
    
    test<divide>(NSEG, D, U, d, u, "divide", cnt);
    test<lapack::xpttrf>(NSEG, D, U, d, u, "lapack", cnt);
    test<lapack_xpttrf>(NSEG, D, U, d, u, "c-lapack", cnt);
    test<italian_factor>(NSEG, D, U, d, u, "italian", cnt);
    test<alsatian_xpttrf>(NSEG, D, U, d, u, "alsatian", cnt);

    free_real(d);
    free_real(u);
    free_real(D);
    free_real(U);
}

//------------------------------------------------------------------------------

template < void (*FACTOR)(int, real*, real*, int*), void (*SOLVE)(int, real const*, real const*, real*) >
void check(int N, real const* D, real const* U, real const* B, real* d, real* u, real* b, real* S, char const str[], size_t cnt)
{
    int info;
    copy_real(N, D, d);
    copy_real(N, U, u);
    copy_real(N, B, b);
    FACTOR(N, d, u, &info);
    SOLVE(N, d, u, b);
    VecPrint::edges(N, b);
    real err = blas::difference(N, S, b);
    printf(" err %f %12s", err, str);
    if ( err < 0.001 )
    {
        tick();
        for ( size_t n = 0; n < cnt; ++n )
            SOLVE(N, d, u, b);
        printf(" cpu %5.0f", tock());
    }
    printf("\n");
}

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testSolve(int NSEG, size_t cnt)
{
    real * d = new_real(NSEG);
    real * u = new_real(NSEG);
    real * b = new_real(NSEG);
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * B = new_real(NSEG);
    real * S = new_real(NSEG);

    for ( int i = 0; i < NSEG; ++i )
    {
        D[i] = 2.0;
        U[i] = -RNG.preal();
        B[i] = RNG.sreal();
    }

    // calculate reference result:
    int info;
    copy_real(NSEG, D, d);
    copy_real(NSEG, U, u);
    copy_real(NSEG, B, S);
    lapack::xpttrf(NSEG, d, u, &info);
    lapack::xptts2(NSEG, d, u, S);

    check<lapack::xpttrf, lapack::xptts2>(NSEG, D, U, B, d, u, b, S, "lapack", cnt);
    check<lapack_xpttrf, lapack_xptts2>(NSEG, D, U, B, d, u, b, S, "c-lapack", cnt);
    check<italian_factor, italian_xptts2>(NSEG, D, U, B, d, u, b, S, "italian", cnt);
    check<alsatian_xpttrf, alsatian_xptts2>(NSEG, D, U, B, d, u, b, S, "alsatian", cnt);

    free_real(d);
    free_real(u);
    free_real(b);
    free_real(D);
    free_real(U);
    free_real(B);
    free_real(S);
}

//------------------------------------------------------------------------------

template < void (*FUNC)(int, real*, real*, real*) >
void verify(int N, real const* D, real const* U, real const* B, real* d, real* u, real* b, char const str[], size_t cnt)
{
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(N, D, d);
        copy_real(N, U, u);
        copy_real(N, B, b);
        FUNC(N, d, u, b);
    }
    VecPrint::edges(N, b, 3);
    printf(" %12s cpu %5.0f\n", str, tock());
}

void solveL(int N, real* D, real* U, real* B)
{
    int info;
    lapack::xpttrf(N, D, U, &info);
    lapack::xptts2(N, D, U, B);
}

void solveC(int N, real* D, real* U, real* B)
{
    int info;
    lapack_xpttrf(N, D, U, &info);
    lapack_xptts2(N, D, U, B);
}

void solveI(int N, real* D, real* U, real* B)
{
    italian_solve(N, U, D, U, B);
}

void solveA(int N, real* D, real* U, real* B)
{
    int info;
    alsatian_xpttrf(N, D, U, &info);
    alsatian_xptts2(N, D, U, B);
}

void solveF(int N, real* D, real* U, real* B)
{
    alsatian_solve(N, D, U, B);
}

void solveW(int N, real* D, real* U, real* B)
{
    wikipedia_solve(N, U, D, U, B);
}

void solveT(int N, real* D, real* U, real* B)
{
    tridiagonal_solve(N, U, D, U, B);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testFused(int NSEG, size_t cnt)
{
    real * d = new_real(NSEG);
    real * u = new_real(NSEG);
    real * b = new_real(NSEG);
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * B = new_real(NSEG);

    for ( int i = 0; i < NSEG; ++i )
    {
        D[i] = 4.0 - RNG.preal();
        U[i] = -RNG.preal();
        B[i] = RNG.sreal();
    }
    
    verify<solveL>(NSEG, D, U, B, d, u, b, "lapack", cnt);
    verify<solveC>(NSEG, D, U, B, d, u, b, "c-lapack", cnt);
    verify<solveI>(NSEG, D, U, B, d, u, b, "italian", cnt);
    verify<solveA>(NSEG, D, U, B, d, u, b, "alsatian", cnt);
    verify<solveF>(NSEG, D, U, B, d, u, b, "alsafused", cnt);
    verify<solveW>(NSEG, D, U, B, d, u, b, "wikipedia", cnt);
    verify<solveT>(NSEG, D, U, B, d, u, b, "tridiagonal", cnt);

    free_real(d);
    free_real(u);
    free_real(b);
    free_real(D);
    free_real(U);
    free_real(B);
}

int main(int argc, char* argv[])
{
    int nbs = 117;
    if ( argc > 1 )
        nbs = std::max(1, atoi(argv[1]));
    
    RNG.seed();
    std::cout << "Tridiagonal positive symmetric matrix --- real " << sizeof(real) << " bytes --- " << __VERSION__ << "\n";
    std::cout << "Factorize\n";
    testFactor(nbs, 1<<20);
    std::cout << "Solve\n";
    testSolve(nbs, 1<<20);
    std::cout << "Factorize & Solve\n";
    testFused(nbs, 1<<18);
}
