// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#define DIM 2

#include "real.h"
#include "vector.h"
#include "random.h"
#include "blas.h"
#include "cytoblas.h"
#include "vecprint.h"
#include "simd.h"

#define FOR 4

/// number of segments:
const size_t NSEG = 127;
const size_t NVAL = FOR * ( NSEG + 1 );
const size_t ALOC = NVAL + 8;

const size_t DISP = 16UL;

real * vP = nullptr;
real * vX = nullptr, * vY = nullptr, * vZ = nullptr;

//------------------------------------------------------------------------------

void setFilament(size_t nbs, real* ptr, real seg, real persistence_length)
{
    nbs = std::min(nbs, NSEG);
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
    p = new_real(NVAL);
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
    
    print_alignment(p, "p");
    print_alignment(x, "x");
    print_alignment(y, "y");
    print_alignment(y, "z");
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
 This is the simple implementation
 */
void add_rigidity0(const size_t nbt, const real* X, const real rigid, real* Y)
{
    #pragma vector unaligned
    for ( size_t jj = 0; jj < nbt; ++jj )
    {
        real f = rigid * (( X[jj+DIM*2] - X[jj+DIM] ) - ( X[jj+DIM] - X[jj] ));
        Y[jj      ] -=   f;
        Y[jj+DIM  ] += 2*f;
        Y[jj+DIM*2] -=   f;
    }
}

/*
 This an implementation for 2D
 */
void add_rigidity2(const size_t nbt, const real* vec, const real rigid, real* Y)
{
    real fx = 0;
    real fy = 0;
    #pragma vector unaligned
    for ( size_t jj = 0; jj < nbt; jj += 2 )
    {
        real const* X = vec + jj;
        real gx = rigid * ( X[4] - X[2] - X[2] + X[0] );
        real gy = rigid * ( X[5] - X[3] - X[3] + X[1] );
        Y[jj  ] += fx-gx;
        Y[jj+1] += fy-gy;
        Y[jj+2] -= fx-gx;
        Y[jj+3] -= fy-gy;
        fx = gx;
        fy = gy;
    }
    Y[nbt  ] += fx;
    Y[nbt+1] += fy;
    Y[nbt+2] -= fx;
    Y[nbt+3] -= fy;
}

/*
 In this version for 2D or 3D, the loop is unrolled, pointers are used,
 and ( a0 -2*a1 + a2 ) is replaced by (a2-a1)-(a1-a0).
 */
void add_rigidity3(const size_t nbt, const real* X, const real rigid, real* Y)
{
    const real * xn = X + DIM;
    
    real x0 = xn[0];
    real x1 = xn[1];
#if ( DIM == 3 )
    real x2 = xn[2];
#endif
    
    real d0 = x0 - X[0];
    real d1 = x1 - X[1];
#if ( DIM == 3 )
    real d2 = x2 - X[2];
#endif
    
    real df0 = 0, of0 = 0, odf0 = 0;
    real df1 = 0, of1 = 0, odf1 = 0;
#if ( DIM == 3 )
    real df2 = 0, of2 = 0, odf2 = 0;
#endif
    
    xn += DIM;
    
    real * yp = Y;
    real *const end = Y + nbt;
    while ( yp < end )
    {
        real e0 = *xn - x0;
        x0 = *xn;
        ++xn;
        real f0 = rigid * ( e0 - d0 );
        d0      = e0;
        df0     = f0 - of0;
        of0     = f0;
        *yp    += odf0 - df0;
        odf0    = df0;
        ++yp;
        
        real e1 = *xn - x1;
        x1 = *xn;
        ++xn;
        real f1 = rigid * ( e1 - d1 );
        d1      = e1;
        df1     = f1 - of1;
        of1     = f1;
        *yp    += odf1 - df1;
        odf1    = df1;
        ++yp;
        
#if ( DIM == 3 )
        real e2 = *xn - x2;
        x2 = *xn;
        ++xn;
        real f2 = rigid * ( e2 - d2 );
        d2      = e2;
        df2     = f2 - of2;
        of2     = f2;
        *yp    += odf2 - df2;
        odf2    = df2;
        ++yp;
#endif
    }
    
    yp[0]   += df0 + of0;
    yp[1]   += df1 + of1;
#if ( DIM == 3 )
    yp[2]   += df2 + of2;
#endif
    
    yp += DIM;
    
    yp[0] -= of0;
    yp[1] -= of1;
#if ( DIM == 3 )
    yp[2] -= of2;
#endif
}

#if REAL_IS_DOUBLE && defined(__SSE3__)

void add_rigidity2D_SSE(const size_t nbt, const real* X, const real rigid, real* Y)
{
    vec2 R = set2(rigid);
    real *const end = Y + nbt;

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
void add_rigidity2D_SSO(const size_t nbt, const real* X, const real rigid, real* Y)
{
    vec2 R = set2(rigid);
    real *const end = Y + nbt;

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

void add_rigidity2D_AVX(const size_t nbt, const real* X, const real rigid, real* Y)
{
    vec4 R = set4(rigid);
    vec4 two = set4(2.0);
    
    real *const end = Y + nbt - 8;
    
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
        store4(Y, fmadd4(R, fmsub4(two, ppp, add4(eee, ddd)), load4(Y)));
        eee = sub4(sub4(xxx, jjj), sub4(jjj, nnn));
        ppp = catshift2(ddd, eee);
        store4(Y+4, fmadd4(R, fmsub4(two, ppp, add4(ddd, eee)), load4(Y+4)));
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
        store4(Y, fmadd4(R, fmsub4(two, ppp, add4(eee, ddd)), load4(Y)));
        eee = ddd;
        Y += 4;
    }
#endif
    vec2 nn = gethi(xxx);
    vec2 oo = sub2(nn, getlo(xxx));
    vec2 ee = gethi(eee);
    vec2 yy = fmsub2(getlo(two), ee, getlo(eee));
    while ( Y < end+8 )
    {
        vec2 mm = load2(X+4);
        X += 2;
        vec2 ff = sub2(mm, nn);
        vec2 dd = sub2(ff, oo);
        nn = mm;
        oo = ff;
        store2(Y, fmadd2(getlo(R), sub2(yy, dd), load2(Y)));
        yy = fmsub2(getlo(two), dd, ee);
        ee = dd;
        Y += 2;
    }
    store2(Y  ,  fmadd2(getlo(R), yy, load2(Y  )));
    store2(Y+2, fnmadd2(getlo(R), ee, load2(Y+2)));
}

#endif

void add_rigidityF(const size_t nbt, const real* X, const real R1, real* Y)
{
    const real R6 = R1 * 6;
    const real R4 = R1 * 4;
    const real R2 = R1 * 2;
/*
    if ( nbt == DIM )
    {
        for ( size_t d = 0; d < DIM; ++d )
        {
            real x = 2 * X[d+DIM] - ( X[d+DIM*2] + X[d] );
            Y[d      ] += R1 * x;
            Y[d+DIM  ] -= R2 * x;
            Y[d+DIM*2] += R1 * x;
        }
        return;
    }
*/
    const size_t end = nbt;
    #pragma ivdep
    for ( size_t i = DIM*2; i < end; ++i )
        Y[i] += R4 * (X[i-DIM]+X[i+DIM]) - R1 * (X[i-DIM*2]+X[i+DIM*2]) - R6 * X[i];
    
    // special cases near the edges:
    real      * Z = Y + nbt;
    real const* E = X + nbt + DIM;
    #pragma ivdep
    for ( int d = 0; d < DIM; ++d )
    {
        Y[d+DIM] -= R1 * (X[d+DIM]+X[d+DIM*3]) + R4 * (X[d+DIM]-X[d+DIM*2]) - R2 * X[d];
        Z[d    ] -= R1 * (E[d-DIM]+E[d-DIM*3]) + R4 * (E[d-DIM]-E[d-DIM*2]) - R2 * E[d];
        Y[d    ] -= R1 * (X[d+DIM*2]+X[d]) - R2 * X[d+DIM];
        Z[d+DIM] -= R1 * (E[d-DIM*2]+E[d]) - R2 * E[d-DIM];
    }
}

/// only valid if ( nbt > DIM )
void add_rigidityG(const size_t nbt, const real* X, const real R1, real* Y)
{
    const real R6 = R1 * 6;
    const real R4 = R1 * 4;
    const real R2 = R1 * 2;

    const size_t end = nbt;
    #pragma ivdep
    for ( size_t i = DIM*2; i < end; ++i )
        Y[i] += R4 * ((X-DIM)[i]+(X+DIM)[i]) - R1 * ((X-DIM*2)[i]+(X+DIM*2)[i]) - R6 * X[i];
    
    // special cases near the edges:
    real      * Z = Y + nbt;
    real const* E = X + nbt + DIM;

    #pragma ivdep
    for ( size_t d = 0; d < DIM; ++d )
    {
        Y[d+DIM] -= R1 * ((X+DIM)[d]+(X+DIM*3)[d]) + R4 * ((X+DIM)[d]-(X+DIM*2)[d]) - R2 * X[d];
        Z[d    ] -= R1 * ((E-DIM)[d]+(E-DIM*3)[d]) + R4 * ((E-DIM)[d]-(E-DIM*2)[d]) - R2 * E[d];
        Y[d    ] -= R1 * ((X+DIM*2)[d]+X[d]) - R2 * (X+DIM)[d];
        Z[d+DIM] -= R1 * ((E-DIM*2)[d]+E[d]) - R2 * (E-DIM)[d];
    }
}

/// only valid if ( nbt > DIM )
void add_rigidity4(const size_t nbt, const real* X, const real R1, real* Y)
{
    const real R6 = R1 * 6;
    const real R4 = R1 * 4;
    const real R2 = R1 * 2;
    
    const size_t end = FOR * nbt / DIM;
    #pragma ivdep
    for ( size_t i = FOR*2; i < end; ++i )
        Y[i] += R4 * ((X-FOR)[i]+(X+FOR)[i]) - R1 * ((X-FOR*2)[i]+(X+FOR*2)[i]) - R6 * X[i];
    
    // special cases near the edges:
    real      * Z = Y + nbt;
    real const* E = X + nbt + FOR;
    
    #pragma ivdep
    for ( size_t d = 0; d < FOR; ++d )
    {
        Y[d+FOR] -= R1 * ((X+FOR)[d]+(X+FOR*3)[d]) + R4 * ((X+FOR)[d]-(X+FOR*2)[d]) - R2 * X[d];
        Z[d    ] -= R1 * ((E-FOR)[d]+(E-FOR*3)[d]) + R4 * ((E-FOR)[d]-(E-FOR*2)[d]) - R2 * E[d];
        Y[d    ] -= R1 * ((X+FOR*2)[d]+X[d]) - R2 * (X+FOR)[d];
        Z[d+FOR] -= R1 * ((E-FOR*2)[d]+E[d]) - R2 * (E-FOR)[d];
    }
}

//------------------------------------------------------------------------------
#pragma mark - TEST Rigidity

/// keeping time using Intel's cycle counters
unsigned long long rdt = 0;
/// start timer
inline void tic() { rdt = __rdtsc(); }
/// return time since last `tic()`
inline double toc(double num) { return double(__rdtsc()-rdt)/num; }


template < void (*FUNC)(const size_t, const real*, real, real*) >
void testRigidity(size_t cnt, char const* str)
{
    const size_t nbt = DIM * ( NSEG - 1 );
    const real alpha = 64.0;
    
    zero_real(ALOC, vX);
    zero_real(ALOC, vY);
    zero_real(ALOC, vZ);
    
    FUNC(nbt, vP, alpha, vX);
    VecPrint::print(std::min(DISP,nbt), vX);
    std::cout << " |";
    VecPrint::print(DIM, vX+NVAL);
    add_rigidity0(nbt, vP, alpha, vY);
    real err = blas::difference(nbt+2*DIM, vX, vY);

    tic();
    for ( size_t i=0; i<cnt; ++i )
    {
        FUNC(nbt, vY, alpha, vZ);
        FUNC(nbt, vZ, alpha, vX);
        FUNC(nbt, vX, alpha, vY);
    }
    if ( abs_real(err) > 64*REAL_EPSILON )
        printf(" XXXX %e ", err);
    else
        printf("  --> %e ", err);
    printf(" %4s %5.2f cycles\n", str, toc(3*cnt*nbt));
}


void test(size_t cnt)
{
    testRigidity<add_rigidity0>(cnt, "0  ");
#if ( DIM == 2 )
    testRigidity<add_rigidity2>(cnt, "2  ");
#endif
    testRigidity<add_rigidity3>(cnt, "3  ");
    testRigidity<add_rigidityF>(cnt, "F  ");
    testRigidity<add_rigidityG>(cnt, "G  ");
    testRigidity<add_rigidityF>(cnt, "F  ");
    testRigidity<add_rigidityG>(cnt, "G  ");
    testRigidity<add_rigidity4>(cnt, "4  ");
#if defined(__SSE3__) & ( DIM == 2 ) & REAL_IS_DOUBLE
    testRigidity<add_rigidity2D_SSO>(cnt, "SSO");
    testRigidity<add_rigidity2D_SSE>(cnt, "SSE");
#endif
#if defined(__AVX__) & ( DIM == 2 ) & REAL_IS_DOUBLE
    testRigidity<add_rigidity2D_AVX>(cnt, "AVX");
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Main

int main(int argc, char* argv[])
{
    RNG.seed();
    new_reals(vP, vX, vY, vZ, 1.0);
    setFilament(NSEG, vP, 0.1, 17.0);
    std::cout << "addRigidity " << DIM << "D,  " << NSEG;
    std::cout << " segments,   " << __VERSION__ << "\n";
    test(1<<20);
    free_reals(vP, vX, vY, vZ);
}
