// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include <sys/time.h>


#include "real.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "assert_macro.h"
#include "xpttrf.h"
#include "cytoblas.h"


/// keeping time using Intel's cycle counters
unsigned long long rdt = 0;
/// start timer
inline void tic() { rdt = __rdtsc(); }
/// stop timer and print time
inline void toc(const char* str, double num) { printf("  %10s %5.2f\n", str, double(__rdtsc()-rdt)/num); }


/// print only 16 scalars from given vector
inline void print(size_t n, real const* vec)
{
    if ( n > 16 )
    {
        VecPrint::print(std::cout, 8, vec, 3);
        VecPrint::print(std::cout, 8, vec+n-8, 3);
    }
    else
    {
        VecPrint::print(std::cout, n, vec, 3);
    }
}


/// benchmark division
void divide(int size, real const* D, real* E)
{
    for ( int i = 0; i < size-1; ++i )
        E[i] = E[i] / D[i];
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testDPTTF(size_t NSEG, size_t cnt)
{
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * Ds = new_real(NSEG);
    real * Us = new_real(NSEG);

    for ( size_t i = 0; i < NSEG; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -0.5 * RNG.preal();
    }
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        //copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        divide(NSEG, D, U);
    }
    toc("divide", NSEG*cnt);
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        //copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        xpttrf(NSEG, D, U);
    }
    toc("pttrf", NSEG*cnt);

    int info;
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        lapack_xpttrf(NSEG, D, U, &info);
    }
    toc("clapack", NSEG*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        lapack::xpttrf(NSEG, D, U, &info);
    }
    toc("lapack", NSEG*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        italian_xpttrf(NSEG, D, U, &info);
    }
    toc("italian", NSEG*cnt);
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        alsatian_xpttrf(NSEG, D, U, &info);
    }
    toc("alsatian", NSEG*cnt);

    free_real(D);
    free_real(U);
    free_real(Ds);
    free_real(Us);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testDPTTS(size_t NSEG, size_t cnt)
{
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * B = new_real(NSEG);
    real * Ds = new_real(NSEG);
    real * Us = new_real(NSEG);
    real * Bs = new_real(NSEG);
    real * S = new_real(NSEG);

    for ( size_t i = 0; i < NSEG; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -RNG.preal();
        Bs[i] = RNG.sreal();
    }

    int info;
    copy_real(NSEG, Ds, D);
    copy_real(NSEG, Us, U);
    copy_real(NSEG, Bs, B);
    lapack_xpttrf(NSEG, D, U, &info);
    lapack_xptts2(NSEG, 1, D, U, B, 1);
    copy_real(NSEG, B, S);
    print(NSEG, B);
    printf(" err %f", blas::max_diff(NSEG, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        lapack_xptts2(NSEG, 1, D, U, B, 1);
    toc("clapack", NSEG*cnt);

    copy_real(NSEG, Ds, D);
    copy_real(NSEG, Us, U);
    copy_real(NSEG, Bs, B);
    lapack::xpttrf(NSEG, D, U, &info);
    lapack::xptts2(NSEG, 1, D, U, B, 1);
    print(NSEG, B);
    printf(" err %f", blas::max_diff(NSEG, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        lapack::xptts2(NSEG, 1, D, U, B, 1);
    toc("lapack", NSEG*cnt);

    copy_real(NSEG, Ds, D);
    copy_real(NSEG, Us, U);
    copy_real(NSEG, Bs, B);
    italian_xpttrf(NSEG, D, U, &info);
    italian_xptts2(NSEG, 1, D, U, B, 1);
    print(NSEG, B);
    printf(" err %f", blas::max_diff(NSEG, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        italian_xptts2(NSEG, 1, D, U, B, 1);
    toc("italian", NSEG*cnt);
    
    copy_real(NSEG, Ds, D);
    copy_real(NSEG, Us, U);
    copy_real(NSEG, Bs, B);
    alsatian_xpttrf(NSEG, D, U, &info);
    alsatian_xptts2(NSEG, 1, D, U, B, 1);
    print(NSEG, B);
    printf(" err %f", blas::max_diff(NSEG, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian_xptts2(NSEG, 1, D, U, B, 1);
    toc("alsatian", NSEG*cnt);

    free_real(D);
    free_real(U);
    free_real(B);
    free_real(Ds);
    free_real(Us);
    free_real(Bs);
    free_real(S);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testThomas(size_t NSEG, size_t cnt)
{
    real * D = new_real(NSEG);
    real * U = new_real(NSEG);
    real * B = new_real(NSEG);
    real * Ds = new_real(NSEG);
    real * Us = new_real(NSEG);
    real * Bs = new_real(NSEG);

    for ( size_t i = 0; i < NSEG; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -RNG.preal();
        Bs[i] = RNG.sreal();
    }
    
    int info = 0;
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        copy_real(NSEG, Bs, B);
        lapack::xpttrf(NSEG, D, U, &info);
        lapack::xptts2(NSEG, 1, D, U, B, 1);
    }
    print(NSEG, B);
    toc("lapack", NSEG*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        copy_real(NSEG, Bs, B);
        //italian_xpttrf(NSEG, D, U, &info);
        //italian_xptts2(NSEG, 1, D, U, B, 1);
        italian_thomas(NSEG, U, D, U, B);
    }
    print(NSEG, B);
    toc("italian", NSEG*cnt);
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        copy_real(NSEG, Bs, B);
        alsatian_xpttrf(NSEG, D, U, &info);
        alsatian_xptts2(NSEG, 1, D, U, B, 1);
    }
    print(NSEG, B);
    toc("alsatian2", NSEG*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Us, U);
        copy_real(NSEG, Bs, B);
        alsatian_thomas(NSEG, D, U, B);
    }
    print(NSEG, B);
    toc("alsatian", NSEG*cnt);

    tic();
    copy_real(NSEG, Us, U);
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NSEG, Ds, D);
        copy_real(NSEG, Bs, B);
        tridiagonal_solve(NSEG, U, D, U, B);
    }
    print(NSEG, B);
    toc("tridiag", NSEG*cnt);

    free_real(D);
    free_real(U);
    free_real(B);
    free_real(Ds);
    free_real(Us);
    free_real(Bs);
}

int main(int argc, char* argv[])
{
    size_t nbs = 117;
    if ( argc > 1)
        nbs = std::max(1, atoi(argv[1]));
    
    RNG.seed();
    std::cout << "testPTTF   --- real " << sizeof(real) << " --- " << __VERSION__ << "\n";
    testDPTTF(nbs, 1<<10);
    std::cout << "testPTTS\n";
    testDPTTS(nbs, 1<<17);
    std::cout << "testThomas\n";
    testThomas(nbs, 1<<15);
    
    return EXIT_SUCCESS;
}
