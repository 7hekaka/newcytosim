// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#define DIM 3

#include "real.h"
#include "timer.h"
#include "vector.h"
#include "random.h"
#include "blas.h"
#include "cytoblas.h"
#include "vecprint.h"
#include "simd.h"

#define FOR 4

const size_t MSEG = 64;
const size_t ALOC = ( MSEG + 1 ) * DIM;

real * vP = nullptr;
real * vX = nullptr, * vY = nullptr, * vZ = nullptr;

//------------------------------------------------------------------------------

void setFilament(size_t nbs, real* ptr, real seg, real persistence_length)
{
    nbs = std::min(nbs, MSEG);
    real sigma = std::sqrt(2.0*seg/persistence_length);
    
    Vector pos(0,0,0);
    Vector dir(1,0,0);
    
    pos.store(ptr);
    for ( size_t p = 1 ; p <= nbs; ++p )
    {
        pos += seg * dir;
        pos.store(ptr+DIM*p);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        dir = std::cos(a) * dir + dir.randOrthoU(std::sin(a));
    }
}

void new_reals(real*& p, real*& x, real*& y, real*& z, real mag)
{
    p = new_real(ALOC);
    x = new_real(ALOC);
    y = new_real(ALOC);
    z = new_real(ALOC);
    
    if ( mag > 0 )
    {
        for ( size_t i=0; i<ALOC; ++i )
        {
            x[i] = mag * RNG.sreal();
            y[i] = mag * RNG.sreal();
            z[i] = mag * RNG.sreal();
        }
    }
    else {
        zero_real(ALOC, vX);
        zero_real(ALOC, vY);
        zero_real(ALOC, vZ);
    }
}

void free_reals(real* p, real* x, real* y, real* z)
{
    free_real(p);
    free_real(x);
    free_real(y);
    free_real(z);
}


//------------------------------------------------------------------------------
#pragma mark - RIGIDITY

/*
 This is the reference implementation
 */
void add_rigidity0(const int nbp, const real* X, const real R1, real* Y)
{
    const real R2 = 2.0 * R1;
    const real two = 2.0;
    #pragma omp simd
    const int end = DIM * ( nbp - 2 );
    for ( int i = 0; i < end; ++i )
    {
        real f = ( X[i+DIM*2] + X[i] ) - two * X[i+DIM];
        Y[i      ] -= f * R1;
        Y[i+DIM  ] += f * R2;
        Y[i+DIM*2] -= f * R1;
    }
}

/// In this version for 2D, the loop has dependencies preventing unrolling
void add_rigidity2D(const int nbp, const real* X, const real R1, real* Y)
{
    real fx = 0;
    real fy = 0;
    real y0 = Y[0];
    real y1 = Y[1];
    const real two = -2.0;
    real const*const end = X + DIM * ( nbp - 2 );
    while ( X < end )
    {
        real gx = X[2] * two + ( X[4] + X[0] );
        real gy = X[3] * two + ( X[5] + X[1] );
        X += DIM;
        real rx = fx - gx;
        real ry = fy - gy;
        fx = gx;
        fy = gy;
        Y[0] = y0 + R1 * rx;
        Y[1] = y1 + R1 * ry;
        y0 = Y[2] - R1 * rx;
        y1 = Y[3] - R1 * ry;
        Y += DIM;
    }
    Y[0] = R1 * fx + y0;
    Y[1] = R1 * fy + y1;
    Y[2] -= R1 * fx;
    Y[3] -= R1 * fy;
}

/// In this version for 3D, the loop has dependencies preventing unrolling
void add_rigidity3D(const int nbp, const real* X, const real R1, real* Y)
{
    real fx = 0;
    real fy = 0;
    real fz = 0;
    real y0 = Y[0];
    real y1 = Y[1];
    real y2 = Y[2];
    const real two = -2.0;
    real const*const end = X + DIM * ( nbp - 2 );
    while ( X < end )
    {
        real gx = X[3] * two + ( X[0] + X[6] );
        real gy = X[4] * two + ( X[1] + X[7] );
        real gz = X[5] * two + ( X[2] + X[8] );
        X += DIM;
        real rx = fx - gx;
        real ry = fy - gy;
        real rz = fz - gz;
        fx = gx;
        fy = gy;
        fz = gz;
        Y[0] = y0 + R1 * rx;
        Y[1] = y1 + R1 * ry;
        Y[2] = y2 + R1 * rz;
        y0 = Y[3] - R1 * rx;
        y1 = Y[4] - R1 * ry;
        y2 = Y[5] - R1 * rz;
        Y += DIM;
    }
    Y[0] = R1 * fx + y0;
    Y[1] = R1 * fy + y1;
    Y[2] = R1 * fz + y2;
    Y[3] -= R1 * fx;
    Y[4] -= R1 * fy;
    Y[5] -= R1 * fz;
}

#if REAL_IS_DOUBLE && USE_SIMD

void add_rigidity2D_SSE(const int nbp, const real* X, const real rigid, real* Y)
{
    vec2 R = set2(rigid);
    real *const end = Y + DIM * ( nbp - 2 );

    vec2 nn = load2(X+2);
    vec2 oo = mul2(R, sub2(nn, load2(X)));
    vec2 yy = load2(Y);
    vec2 zz = load2(Y+2);
    
    while ( Y < end )
    {
        vec2 mm = load2(X+4);
        X += 2;
        vec2 dd = mul2(R, sub2(mm, nn));
        vec2 ff = sub2(dd, oo);
        oo = dd;
        nn = mm;
        store2(Y, sub2(yy, ff));
        yy = add2(zz, add2(ff, ff));
        zz = sub2(load2(Y+4), ff);
        Y += 2;
    }
    store2(Y, yy);
    store2(Y+2, zz);
}

/// older implementation
void add_rigidity2D_SSO(const int nbp, const real* X, const real rigid, real* Y)
{
    vec2 R = set2(rigid);
    real *const end = Y + DIM * ( nbp - 2 );

    vec2 xx  = load2(X+DIM);
    vec2 d   = sub2(xx, load2(X));
    vec2 df  = setzero2();
    vec2 of  = setzero2();
    vec2 yy  = load2(Y);
    
    X += 2*DIM;
    while ( Y < end )
    {
        vec2 nn = load2(X);
        X += DIM;
        vec2 e = sub2(nn, xx);
        xx = nn;
        
        vec2 f = mul2(R, sub2(e, d));
        d  = e;
        df = sub2(f, of);
        of = f;
        store2(Y, sub2(yy, df));
        yy = add2(load2(Y+DIM), df);
        Y += DIM;
    }
    store2(Y, add2(load2(Y), add2(df, of)));
    store2(Y+DIM, sub2(load2(Y+DIM), of));
}

#endif
#if REAL_IS_DOUBLE && defined(__AVX__)

void add_rigidity2D_AVX(const size_t nbp, const real* X, const real rigid, real* Y)
{
    vec4 R = set4(rigid);
    vec4 two = set4(2.0);
    
    real *const end = Y + DIM * ( nbp - 2 ) - 8;
    
    vec4 xxx = load4(X);
    vec4 eee = setzero4();
#if 1
    // unrolled 2x2
    while ( Y < end )
    {
        vec4 nnn = load4(X+4);
        vec4 iii = catshift2(xxx, nnn);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = load4(X+8);
        X += 8;
        vec4 ppp = catshift2(eee, ddd);
        vec4 jjj = catshift2(nnn, xxx);
        store4(Y, fnmadd4(R, fnmadd4(two, ppp, add4(eee, ddd)), load4(Y)));
        eee = sub4(sub4(xxx, jjj), sub4(jjj, nnn));
        ppp = catshift2(ddd, eee);
        store4(Y+4, fnmadd4(R, fnmadd4(two, ppp, add4(ddd, eee)), load4(Y+4)));
        Y += 8;
    }
#endif
#if 1
    if ( Y < end+4 )
    {
        vec4 nnn = load4(X+4);
        vec4 iii = catshift2(xxx, nnn);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = nnn;
        X += 4;
        vec4 ppp = catshift2(eee, ddd);
        store4(Y, fnmadd4(R, fnmadd4(two, ppp, add4(eee, ddd)), load4(Y)));
        eee = ddd;
        Y += 4;
    }
#endif
    vec2 nn = gethi(xxx);
    vec2 oo = sub2(nn, getlo(xxx));
    vec2 ee = gethi(eee);
    vec2 yy = fnmadd2(getlo(two), ee, getlo(eee));
    while ( Y < end+8 )
    {
        vec2 mm = load2(X+4);
        X += 2;
        vec2 ff = sub2(mm, nn);
        vec2 dd = sub2(ff, oo);
        nn = mm;
        oo = ff;
        store2(Y, fmadd2(getlo(R), add2(yy, dd), load2(Y)));
        yy = fnmadd2(getlo(two), dd, ee);
        ee = dd;
        Y += 2;
    }
    store2(Y  , fnmadd2(getlo(R), yy, load2(Y  )));
    store2(Y+2, fnmadd2(getlo(R), ee, load2(Y+2)));
}

#endif

void add_rigidityF(const int nbp, const real* X, const real R1, real* Y)
{
    const real six = 6.0;
    const real R4 = R1 * 4.0;
    const real R2 = R1 * 2.0;

    const int end = DIM * ( nbp - 2 );
    #pragma omp simd
    for ( int i = DIM*2; i < end; ++i )
        Y[i] = Y[i] + R4 * (X[i-DIM]+X[i+DIM]) - R1 * (six*X[i]+(X[i-DIM*2]+X[i+DIM*2]));
    
    // special cases near the edges:
    real      * Z = Y + DIM * ( nbp - 1 );
    real const* E = X + DIM * ( nbp - 1 );
    if ( nbp > 3 )
    {
        #pragma omp simd
        for ( int d = 0; d < DIM; ++d )
        {
            Y[d+DIM] -= R1 * (X[d+DIM]+X[d+DIM*3]) + R4 * (X[d+DIM]-X[d+DIM*2]) - R2 * X[d];
            Z[d-DIM] -= R1 * (E[d-DIM]+E[d-DIM*3]) + R4 * (E[d-DIM]-E[d-DIM*2]) - R2 * E[d];
        }
    }
    else
    {
        #pragma omp simd
        for ( int d = 0; d < DIM; ++d )
            Y[d+DIM] += R2 * (X[d+DIM*2]+X[d]) - R4 * X[d+DIM];
    }
    #pragma omp simd
    for ( int d = 0; d < DIM; ++d )
    {
        Y[d] -= R1 * (X[d+DIM*2]+X[d]) - R2 * X[d+DIM];
        Z[d] -= R1 * (E[d-DIM*2]+E[d]) - R2 * E[d-DIM];
    }
}

//------------------------------------------------------------------------------
#pragma mark - TEST Rigidity

template < void (*FUNC)(const int, const real*, real, real*) >
void testRigidity(int nbp, size_t cnt, char const* str)
{
    const real alpha = 64.0;
    
    zero_real(ALOC, vX);
    zero_real(ALOC, vY);
    zero_real(ALOC, vZ);
    
    FUNC(nbp, vP, alpha, vX);
    VecPrint::edges(DIM*nbp, vX);
    fprintf(stderr, " |");
    VecPrint::print(DIM, vX+DIM*nbp);
    add_rigidity0(nbp, vP, alpha, vY);
    real err = blas::difference(DIM*nbp, vX, vY);
    fprintf(stderr, " |");

    tick();
    for ( size_t i = 0; i < cnt; ++i )
    {
        FUNC(nbp, vY, alpha, vZ);
        FUNC(nbp, vZ, alpha, vX);
        FUNC(nbp, vX, alpha, vY);
    }
    if ( abs_real(err) > 64*REAL_EPSILON )
        printf(" XXXX %e ", err);
    else
        printf("  --> %e ", err);
    printf(" %4s cpu %5.2f\n", str, tock());
}


void test(int num, size_t cnt)
{
    testRigidity<add_rigidity0>(num, cnt, "0  ");
#if ( DIM == 2 )
    testRigidity<add_rigidity2D>(num, cnt, "2D ");
#else
    testRigidity<add_rigidity3D>(num, cnt, "3D ");
#endif
    testRigidity<add_rigidityF>(num, cnt, "F  ");
#if USE_SIMD & ( DIM == 2 ) & REAL_IS_DOUBLE
    testRigidity<add_rigidity2D_SSO>(num, cnt, "SSO");
    testRigidity<add_rigidity2D_SSE>(num, cnt, "SSE");
#endif
#if defined(__AVX__) & ( DIM == 2 ) & REAL_IS_DOUBLE
    testRigidity<add_rigidity2D_AVX>(num, cnt, "AVX");
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Main

int main(int argc, char* argv[])
{
    RNG.seed();
    new_reals(vP, vX, vY, vZ, 1.0);
    for ( int i : { 0, 1, 2, 3, 5, 7, 10, 12, 17, 31, 47 } )
    {
        setFilament(i, vP, 0.1, 17.0);
        std::cout << "addRigidity" << DIM << "D,  ";
        std::cout << i+1 << " segments,   " << __VERSION__ << "\n";
        test(i, 1<<21);
    }
    free_reals(vP, vX, vY, vZ);
}
