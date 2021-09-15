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
void test(int N, real const* Ds, real const* Us, real* D, real* U, char const str[], size_t cnt)
{
    int info;
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(N, Ds, D);
        copy_real(N, Us, U);
        FUNC(N, D, U, &info);
    }
    printf(" %12s %5.2f\n", str, toc(N*cnt));
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
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * Ds = new_real(NSEG);
    real * Us = new_real(NSEG);

    for ( int i = 0; i < NSEG; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -0.5 * RNG.preal();
    }
    
    test<divide>(NSEG, Ds, Us, D, U, "divide", cnt);
    test<lapack_xpttrf>(NSEG, Ds, Us, D, U, "clapack", cnt);
    test<lapack::xpttrf>(NSEG, Ds, Us, D, U, "lapack", cnt);
    test<italian_factor>(NSEG, Ds, Us, D, U, "italian", cnt);
    test<alsatian_xpttrf>(NSEG, Ds, Us, D, U, "alsatian", cnt);

    free_real(D);
    free_real(U);
    free_real(Ds);
    free_real(Us);
}

//------------------------------------------------------------------------------

template < void (*FACTOR)(int, real*, real*, int*), void (*SOLVE)(int, real const*, real const*, real*) >
void check(int N, real const* Ds, real const* Us, real const* Bs, real* D, real* U, real* B, real* S, char const str[], size_t cnt)
{
    int info;
    copy_real(N, Ds, D);
    copy_real(N, Us, U);
    copy_real(N, Bs, B);
    FACTOR(N, D, U, &info);
    SOLVE(N, D, U, B);
    VecPrint::edges(N, B);
    real err = blas::difference(N, S, B);
    printf(" err %f %12s", err, str);
    if ( err < 0.001 )
    {
        tic();
        for ( size_t n = 0; n < cnt; ++n )
            SOLVE(N, D, U, B);
        printf(" %5.2f", toc(N*cnt));
    }
    printf("\n");
}

/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testSolve(int NSEG, size_t cnt)
{
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * B = new_real(NSEG);
    real * Ds = new_real(NSEG);
    real * Us = new_real(NSEG);
    real * Bs = new_real(NSEG);
    real * S = new_real(NSEG);

    for ( int i = 0; i < NSEG; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -RNG.preal();
        Bs[i] = RNG.sreal();
    }

    // calculate reference result:
    int info;
    copy_real(NSEG, Ds, D);
    copy_real(NSEG, Us, U);
    copy_real(NSEG, Bs, S);
    lapack::xpttrf(NSEG, D, U, &info);
    lapack::xptts2(NSEG, D, U, S);

    check<lapack_xpttrf, lapack_xptts2>(NSEG, Ds, Us, Bs, D, U, B, S, "clapack", cnt);
    check<lapack::xpttrf, lapack::xptts2>(NSEG, Ds, Us, Bs, D, U, B, S, "lapack", cnt);
    check<italian_factor, italian_xptts2>(NSEG, Ds, Us, Bs, D, U, B, S, "italian", cnt);
    check<alsatian_xpttrf, alsatian_xptts2>(NSEG, Ds, Us, Bs, D, U, B, S, "alsatian", cnt);

    free_real(D);
    free_real(U);
    free_real(B);
    free_real(Ds);
    free_real(Us);
    free_real(Bs);
    free_real(S);
}

//------------------------------------------------------------------------------

template < void (*FUNC)(int, real*, real*, real*) >
void verify(int N, real const* Ds, real const* Us, real const* Bs, real* D, real* U, real* B, char const str[], size_t cnt)
{
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(N, Ds, D);
        copy_real(N, Us, U);
        copy_real(N, Bs, B);
        alsatian_thomas(N, D, U, B);
    }
    VecPrint::edges(N, B);
    printf(" %12s %5.2f\n", str, toc(N*cnt));
}

void solve1(int N, real* D, real* U, real* B)
{
    int info;
    lapack::xpttrf(N, D, U, &info);
    lapack::xptts2(N, D, U, B);
}

void solve2(int N, real* D, real* U, real* B)
{
    italian_thomas(N, U, D, U, B);
}

void solve3(int N, real* D, real* U, real* B)
{
    int info;
    alsatian_xpttrf(N, D, U, &info);
    alsatian_xptts2(N, D, U, B);
}

void solve4(int N, real* D, real* U, real* B)
{
    alsatian_thomas(N, D, U, B);
}

void solve5(int N, real* D, real* U, real* B)
{
    tridiagonal_solve(N, U, D, U, B);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testThomas(int NSEG, size_t cnt)
{
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * B = new_real(NSEG);
    real * Ds = new_real(NSEG);
    real * Us = new_real(NSEG);
    real * Bs = new_real(NSEG);

    for ( int i = 0; i < NSEG; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -RNG.preal();
        Bs[i] = RNG.sreal();
    }
    
    verify<solve1>(NSEG, Ds, Us, Bs, D, U, B, "lapack", cnt);
    verify<solve2>(NSEG, Ds, Us, Bs, D, U, B, "italian", cnt);
    verify<solve3>(NSEG, Ds, Us, Bs, D, U, B, "alsatian2", cnt);
    verify<solve4>(NSEG, Ds, Us, Bs, D, U, B, "alsatian", cnt);
    verify<solve5>(NSEG, Ds, Us, Bs, D, U, B, "tridiag", cnt);

    free_real(D);
    free_real(U);
    free_real(B);
    free_real(Ds);
    free_real(Us);
    free_real(Bs);
}

int main(int argc, char* argv[])
{
    int nbs = 117;
    if ( argc > 1)
        nbs = std::max(1, atoi(argv[1]));
    
    RNG.seed();
    std::cout << "Tridiagonal positive symmetric matrix --- real " << sizeof(real) << " bytes --- " << __VERSION__ << "\n";
    std::cout << "testPTTRF\n";
    testFactor(nbs, 1<<10);
    std::cout << "testPTTS\n";
    testSolve(nbs, 1<<17);
    std::cout << "testThomas\n";
    testThomas(nbs, 1<<15);
}
