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
#include "../math/xptsl.cc"

//------------------------------------------------------------------------------

/// reference code
inline void solveL(int N, real* D, real* U, real* B)
{
    int info;
    lapack::xpttrf(N, D, U, &info);
    lapack::xptts2(N, D, U, B);
}

/// benchmark division
inline void divide(int size, real* D, real* E, int*)
{
    for ( int i = 0; i < size-1; ++i )
        E[i] = E[i] / D[i];
}

inline void italian_factor(int size, real* D, real* E, int* info)
{
    italian_xpttrf(size, D, E, info);
}


template < void (*FUNC)(int, real*, real*, int*) >
void factor(int N, real const* D, real const* U, real* d, real* u, char const str[], size_t cnt)
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

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testFactor(int seg, size_t cnt)
{
    real * d = new_real(seg);
    real * u = new_real(seg);
    real * D = new_real(seg);
    real * U = new_real(seg);

    for ( int i = 0; i < seg; ++i )
    {
        D[i] = 2.0;
        U[i] = -0.5 * RNG.preal();
    }
    
    factor<divide>(seg, D, U, d, u, "divide", cnt);
    factor<lapack::xpttrf>(seg, D, U, d, u, "lapack", cnt);
    factor<lapack_xpttrf>(seg, D, U, d, u, "c-lapack", cnt);
    factor<italian_factor>(seg, D, U, d, u, "italian", cnt);
    factor<alsatian_xpttrf>(seg, D, U, d, u, "alsatian", cnt);
    factor<alsadual_xpttrf>(seg, D, U, d, u, "alsadual", cnt);

    free_real(d);
    free_real(u);
    free_real(D);
    free_real(U);
}

//------------------------------------------------------------------------------
#pragma mark -

template < void (*FACTOR)(int, real*, real*, int*), void (*SOLVE)(int, real const*, real const*, real*) >
void solve(int N, real const* D, real const* U, real const* B, real* d, real* u, real* b, real* S, char const str[], size_t cnt)
{
    int info;
    copy_real(N, D, d);
    copy_real(N, U, u);
    copy_real(N, B, b);
    FACTOR(N, d, u, &info);
    SOLVE(N, d, u, b);
    VecPrint::edges(N, b);
    real err = blas::difference(N, S, b);
    printf("  Err %f %12s", err, str);
    if ( err < 0.001 )
    {
        tick();
        for ( size_t n = 0; n < cnt; ++n )
            SOLVE(N, d, u, b);
        printf(" cpu %5.0f", tock());
    }
    printf("\n");
#if 0
    for ( int i = 0; i < N; ++i )
    {
        zero_real(N, b);
        b[i] = 1;
        SOLVE(N, d, u, b);
        VecPrint::print("t", N, b, 3);
    }
#endif
}

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testSolve(int seg, size_t cnt)
{
    real * d = new_real(seg);
    real * u = new_real(seg);
    real * b = new_real(seg);
    real * D = new_real(seg);
    real * U = new_real(seg);
    real * B = new_real(seg);
    real * S = new_real(seg);

    for ( int i = 0; i < seg; ++i )
    {
        D[i] = 2.0;
        U[i] = -1;
        B[i] = RNG.sreal();
    }

    // calculate reference result:
    copy_real(seg, D, d);
    copy_real(seg, U, u);
    copy_real(seg, B, S);
    solveL(seg, d, u, S);
    solve<lapack::xpttrf, lapack::xptts2>(seg, D, U, B, d, u, b, S, "lapack", cnt);
    solve<lapack_xpttrf, lapack_xptts2>(seg, D, U, B, d, u, b, S, "c-lapack", cnt);
    solve<italian_factor, italian_xptts2>(seg, D, U, B, d, u, b, S, "italian", cnt);
    solve<alsatian_xpttrf, alsatian_xptts2>(seg, D, U, B, d, u, b, S, "alsatian", cnt);
    solve<alsadual_xpttrf, alsadual_xptts2>(seg, D, U, B, d, u, b, S, "alsadual", cnt);

    free_real(d);
    free_real(u);
    free_real(b);
    free_real(D);
    free_real(U);
    free_real(B);
    free_real(S);
}

//------------------------------------------------------------------------------
#pragma mark -

inline void solveC(int N, real* D, real* U, real* B)
{
    int info;
    lapack_xpttrf(N, D, U, &info);
    lapack_xptts2(N, D, U, B);
}

inline void solveI(int N, real* D, real* U, real* B)
{
    italian_solve(N, U, D, U, B);
}

inline void solveA(int N, real* D, real* U, real* B)
{
    int info;
    alsatian_xpttrf(N, D, U, &info);
    alsatian_xptts2(N, D, U, B);
}

inline void solveF(int N, real* D, real* U, real* B)
{
    alsatian_solve(N, D, U, B);
}

inline void solveW(int N, real* D, real* U, real* B)
{
    wikipedia_solve(N, U, D, U, B);
}

inline void solveN(int N, real* D, real* U, real* B)
{
    nr_tridag(N, U, D, U, B);
}

inline void solveJ(int N, real* D, real* U, real* B)
{
    linpack_xptsl(N, D, U, B);
}

inline void solveX(int N, real* D, real* U, real* B)
{
    alsatian_xptsl(N, D, U, B);
}


template < void (*FUNC)(int, real*, real*, real*) >
void fused(int N, real const* D, real const* U, real const* B, real* d, real* u, real* b, real* S, char const str[], size_t cnt)
{
    tick();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(N, D, d);
        copy_real(N, U, u);
        copy_real(N, B, b);
        FUNC(N, d, u, b);
    }
    double cpu = tock();
    VecPrint::edges(N, b, 2);
    real err = blas::difference(N, S, b);
    printf("  Err %f %12s cpu %5.0f\n", err, str, cpu);
}

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testFused(int seg, size_t cnt)
{
    real * d = new_real(seg);
    real * u = new_real(seg);
    real * b = new_real(seg);
    real * D = new_real(seg);
    real * U = new_real(seg);
    real * B = new_real(seg);
    real * S = new_real(seg);

    for ( int i = 0; i < seg; ++i )
    {
        D[i] = 4.0 - RNG.preal();
        U[i] = -RNG.preal();
        B[i] = RNG.sreal();
    }
    // calculate reference result:
    copy_real(seg, D, d);
    copy_real(seg, U, u);
    copy_real(seg, B, S);
    solveL(seg, d, u, S);

#if 1
    fused<solveL>(seg, D, U, B, d, u, b, S, "lapack", cnt);
    fused<solveI>(seg, D, U, B, d, u, b, S, "italian", cnt);
    fused<solveA>(seg, D, U, B, d, u, b, S, "alsatian", cnt);
    fused<solveW>(seg, D, U, B, d, u, b, S, "wikipedia", cnt);
#endif
    fused<solveC>(seg, D, U, B, d, u, b, S, "c-lapack", cnt);
    fused<solveN>(seg, D, U, B, d, u, b, S, "n_recipee", cnt);
    fused<solveJ>(seg, D, U, B, d, u, b, S, "linpack", cnt);
    fused<solveF>(seg, D, U, B, d, u, b, S, "alsafused", cnt);
    fused<solveX>(seg, D, U, B, d, u, b, S, "alsadual", cnt);

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
    const size_t rep = 1<<20;
    if ( argc > 1 )
        nbs = std::max(1, atoi(argv[1]));
    
    RNG.seed();
    std::cout << "Tridiagonal positive symmetric matrix --- real " << sizeof(real) << " bytes --- " << __VERSION__ << "\n";

    std::cout << "Factorize\n";
    testFactor(nbs, rep);
    std::cout << "Solve\n";
    testSolve(nbs, rep);
    testSolve(nbs+1, rep);
    std::cout << nbs << " Factorize & Solve\n";
    testFused(nbs, rep);
    std::cout << nbs+1 << " Factorize & Solve\n";
    testFused(nbs+1, rep);
}
