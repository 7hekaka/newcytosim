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


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testDPTTF(size_t NBS, size_t cnt)
{
    real * D = new_real(NBS);
    real * U = new_real(NBS);
    real * Ds = new_real(NBS);
    real * Us = new_real(NBS);

    for ( size_t i = 0; i < NBS; ++i )
    {
        Ds[i] = 2.0;
        Us[i] = -0.5 * RNG.preal();
    }

    int info;
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        lapack_xpttrf(NBS, D, U, &info);
    }
    toc("clapack", NBS*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        lapack::xpttrf(NBS, D, U, &info);
    }
    toc("lapack", NBS*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        italian_xpttrf(NBS, D, U, &info);
    }
    toc("italian", NBS*cnt);
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        alsatian_xpttrf(NBS, D, U, &info);
    }
    toc("alsatian", NBS*cnt);

    free_real(D);
    free_real(U);
    free_real(Ds);
    free_real(Us);
}


/**
 Test Lapack and custom implementation of routines used to factorize
 a symmetric tri-diagonal matrix and solve the associated system.
 */
void testDPTTS(size_t NBS, size_t cnt)
{
    real * D = new_real(NBS);
    real * U = new_real(NBS);
    real * B = new_real(NBS);
    real * Ds = new_real(NBS);
    real * Us = new_real(NBS);
    real * Bs = new_real(NBS);
    real * S = new_real(NBS);

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
    copy_real(NBS, B, S);
    print(NBS, B);
    printf(" err %f", blas::max_diff(NBS, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        lapack_xptts2(NBS, 1, D, U, B, 1);
    toc("clapack", NBS*cnt);

    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    lapack::xpttrf(NBS, D, U, &info);
    lapack::xptts2(NBS, 1, D, U, B, 1);
    print(NBS, B);
    printf(" err %f", blas::max_diff(NBS, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        lapack::xptts2(NBS, 1, D, U, B, 1);
    toc("lapack", NBS*cnt);

    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    italian_xpttrf(NBS, D, U, &info);
    italian_xptts2(NBS, 1, D, U, B, 1);
    print(NBS, B);
    printf(" err %f", blas::max_diff(NBS, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        italian_xptts2(NBS, 1, D, U, B, 1);
    toc("italian", NBS*cnt);
    
    copy_real(NBS, Ds, D);
    copy_real(NBS, Us, U);
    copy_real(NBS, Bs, B);
    alsatian_xpttrf(NBS, D, U, &info);
    alsatian_xptts2(NBS, 1, D, U, B, 1);
    print(NBS, B);
    printf(" err %f", blas::max_diff(NBS, S, B));
    tic();
    for ( size_t n = 0; n < cnt; ++n )
        alsatian_xptts2(NBS, 1, D, U, B, 1);
    toc("alsatian", NBS*cnt);

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
void testThomas(size_t NBS, size_t cnt)
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
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        lapack::xpttrf(NBS, D, U, &info);
        lapack::xptts2(NBS, 1, D, U, B, 1);
    }
    print(NBS, B);
    toc("lapack", NBS*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        //italian_xpttrf(NBS, D, U, &info);
        //italian_xptts2(NBS, 1, D, U, B, 1);
        italian_thomas(NBS, U, D, U, B);
    }
    print(NBS, B);
    toc("italian", NBS*cnt);
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        alsatian_xpttrf(NBS, D, U, &info);
        alsatian_xptts2(NBS, 1, D, U, B, 1);
    }
    print(NBS, B);
    toc("alsatian2", NBS*cnt);

    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Us, U);
        copy_real(NBS, Bs, B);
        alsatian_thomas(NBS, D, U, B);
    }
    print(NBS, B);
    toc("alsatian", NBS*cnt);

    tic();
    copy_real(NBS, Us, U);
    for ( size_t n = 0; n < cnt; ++n )
    {
        copy_real(NBS, Ds, D);
        copy_real(NBS, Bs, B);
        tridiagonal_solve(NBS, U, D, U, B);
    }
    print(NBS, B);
    toc("tridiag", NBS*cnt);

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
