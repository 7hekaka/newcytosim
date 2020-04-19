// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "mecafil.h"
#include "blas.h"
#include "lapack.h"
#include "matsym.h"
#include "random.h"
#include "vecprint.h"

//------------------------------------------------------------------------------
Mecafil::Mecafil()
{
    buildProjection();
    iPointMobility = 0;
    iRigidity = 0;
    iDir = nullptr;
    iLag = nullptr;
    iLLG = nullptr;
    iVTP = nullptr;
#if NEW_ANISOTROPIC_FIBER_DRAG
    iAni = nullptr;
#endif
    useProjectionDiff = false;
}


Mecafil::~Mecafil()
{
    destroyProjection();
    free_real(iDir);
    iDir = nullptr;
    iLag = nullptr;
    iLLG = nullptr;
    iVTP = nullptr;
}


//------------------------------------------------------------------------------
Mecafil::Mecafil(Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Mecafil");
}


Mecafil& Mecafil::operator=(Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Mecafil");
}


//------------------------------------------------------------------------------
size_t Mecafil::allocateMecable(const size_t nbp)
{
    size_t ms = Mecable::allocateMecable(nbp);
    /*
     if Mecable::allocateMecable() allocated memory, it will return the 
     size of the new array, and we allocate the same size for other arrays.
     */
    if ( ms )
    {
        //std::clog << "Mecafil::allocatePoints " << ms << std::endl;
        allocateProjection(ms);
        
        // allocate memory:
        free_real(iDir);
        
#if NEW_ANISOTROPIC_FIBER_DRAG
        iDir = new_real(ms*(4*DIM+1));
        iLag  = iDir + ms*DIM;
        iLLG  = iLag + ms;
        iAni  = iLLG + ms*DIM;
        iVTP  = iAni + ms*DIM;
#else
        iDir = new_real(ms*(4*DIM+1));
        iLag  = iDir + ms*DIM;
        iLLG  = iLag + ms;
        iVTP  = iLLG + ms*DIM;
#endif
        
        // reset Lagrange multipliers
        zero_real(ms, iLag);
    }
    return ms;
}

void Mecafil::release()
{
    free_real(iDir);
    iDir = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The argument should be: sc = kT / dt;
 */
real Mecafil::addBrownianForces(real const* rnd, real sc, real* rhs) const
{
    real b = sqrt( 2 * sc / iPointMobility );

    for ( size_t jj = 0; jj < DIM*nPoints; ++jj )
        rhs[jj] += b * rnd[jj];
    
    return b * iPointMobility;
}


//------------------------------------------------------------------------------

/**
 Calculate the normalized difference between successive vertices of the fiber:

     for ( int n = 0; n < DIM*lastPoint(); ++n )
         iDir[n] = ( pPos[n+DIM] - pPos[n] ) / segmentation();

 */

void Mecafil::storeDirections()
{
#if ( 1 )
    //checkSegmentation(0.01);
    /*
     assume here that successive points are correctly separated, which is usally
     the case, such that any error would be small
     */
    const real sc  = 1.0 / segmentation();
    const size_t end = DIM * lastPoint();
    #pragma ivdep
    for ( size_t p = 0; p < end; ++p )
        iDir[p] = sc * ( pPos[p+DIM] - pPos[p] );
#else
    for ( size_t p = 0; p < lastPoint(); ++p )
        normalize(diffPoints(p)).store(iDir+DIM*p);
#endif
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    /*
     Calculate the average filament direction at each vertex of the fiber in iAni[].
     Note that iAni[] is calculated here from iDir[]
     */
    
    // for the extremities, the direction of the nearby segment is used.
    for ( size_t d = 0; d < DIM; ++d )
    {
        iAni[d]     = iDir[d];
        iAni[d+end] = iDir[d+end-DIM];
    }
    
    // for intermediate points, the directions of the flanking segments are averaged
    for ( size_t p = DIM ; p < end; ++p )
        iAni[p] = 0.5 * ( iDir[p-DIM] + iDir[p] );

    //VecPrint::print(std::clog, last+DIM, iAni);
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Project

#if ( DIM > 1 )

#  if PROJECT_WITH_MATRIX
#     include "mecafil_projectmat.cc"
#     warning "projection matrices are built explicitly"
#  else
#     include "mecafil_project.cc"
#  endif

#else

void Mecafil::buildProjection()   {}  //DIM == 1
void Mecafil::makeProjection()    {}  //DIM == 1
void Mecafil::destroyProjection() {}  //DIM == 1
void Mecafil::allocateProjection(size_t) {}  //DIM == 1

void Mecafil::projectForces(const real* X, real* Y) const
{
    real sum = X[0];
    for ( size_t ii = 1; ii < nPoints; ++ii )
        sum += X[ii];
    
    sum = sum / (real)nPoints;
#if NEW_ANISOTROPIC_FIBER_DRAG
    sum *= 2;
#endif
    for ( size_t ii = 0; ii < nPoints; ++ii )
        Y[ii] = sum;
}

void Mecafil::computeTensions(const real*) {} //DIM == 1
void Mecafil::storeTensions(const real*) {} //DIM == 1
void Mecafil::makeProjectionDiff(const real*) {} //DIM == 1
void Mecafil::addProjectionDiff(const real*, real*) const {} //DIM == 1

#endif


void Mecafil::printTensions(std::ostream& os) const
{
    os << "\n" << reference() << " ";
    VecPrint::print(os, nbSegments(), iLag, 2);
    os << " end " << std::fixed << std::setprecision(2) << netForceEndM() << "   " << netForceEndP();
}

//-----------------------------------------------------------------------
#pragma mark -

/*
 This is the reference implementation
 */
void add_rigidity0(const size_t nbt, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    for ( size_t jj = 0; jj < nbt; ++jj )
    {
        real f = rigid * (( X[jj+DIM*2] - X[jj+DIM] ) - ( X[jj+DIM] - X[jj] ));
        Y[jj      ] -=   f;
        Y[jj+DIM  ] += 2*f;
        Y[jj+DIM*2] -=   f;
    }
}

/*
 In this version the loop is unrolled
 */
void add_rigidity2(const size_t nbt, const real* X, const real rigid, real* Y)
{    
    assert_true( X != Y );
    real dx = X[DIM  ] - X[0];
    real dy = X[DIM+1] - X[1];
#if ( DIM > 2 )
    real dz = X[DIM+2] - X[2];
#endif
    
    real * yv = Y;
    const real*const end = X + nbt + DIM;
    
    const real* xv = X+DIM;
    while ( xv < end )
    {
        real d0 = xv[DIM] - xv[0];
        real f0 = rigid * ( d0 - dx );
        dx = d0;
        yv[0    ] -=    f0;
        yv[DIM  ] += f0+f0;
        yv[DIM*2] -=    f0;
        ++yv;
        ++xv;
        
        real d1 = xv[DIM] - xv[0];
        real f1 = rigid * ( d1 - dy );
        dy = d1;
        yv[0    ] -=    f1;
        yv[DIM  ] += f1+f1;
        yv[DIM*2] -=    f1;
        ++yv;
        ++xv;
        
#if ( DIM > 2 )
        real d2 = xv[DIM] - xv[0];
        real f2 = rigid * ( d2 - dz );
        dz = d2;
        Y[0    ] -=    f2;
        Y[DIM  ] += f2*f2;
        Y[DIM*2] -=    f2;
        ++yv;
        ++xv;
#endif
    }
}


/*
 In this version the loop is unrolled, pointers are used
 and further optimization are made by replacing
 ( a0 -2*a1 + a2 ) by (a2-a1)-(a1-a0).
 */
void add_rigidity3(const size_t nbt, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    const real * xn = X + DIM;
    
    real x0 = xn[0];
    real x1 = xn[1];
#if ( DIM >= 3 )
    real x2 = xn[2];
#endif
    
    real d0 = x0 - X[0];
    real d1 = x1 - X[1];
#if ( DIM >= 3 )
    real d2 = x2 - X[2];
#endif
    
    real df0 = 0, of0 = 0, odf0 = 0;
    real df1 = 0, of1 = 0, odf1 = 0;
#if ( DIM >= 3 )
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
        
#if ( DIM >= 3 )
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
#if ( DIM >= 3 )
    yp[2]   += df2 + of2;
#endif
    
    yp += DIM;
    
    yp[0] -= of0;
    yp[1] -= of1;
#if ( DIM >= 3 )
    yp[2] -= of2;
#endif
}


#if ( DIM == 2 ) && REAL_IS_DOUBLE

#include "simd.h"

#ifdef __SSE3__
/**
 2D implemention using SSE 128bit vector instructions with double precision
 */
void add_rigidity_SSE(const size_t nbt, const real* X, const real rigid, real* Y)
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
#endif


#ifdef __AVX__
/**
 2D implemention using AVX 256bit vector instructions with double precision
 FJN 15.09.2018 -- 17.09.2018
 
 Note that the vectors X and Y are not aligned to memory!
 */
void add_rigidity_AVX(const size_t nbt, const real* X, const real rigid, real* Y)
{
    vec4 R = set4(rigid);
    vec4 two = set4(2.0);
    
    real *const end = Y + nbt - 8;
    
    vec4 xxx = loadu4(X);
    vec4 eee = setzero4();

    // process data 8 by 8:
    while ( Y < end )
    {
        vec4 nnn = loadu4(X+4);
        vec4 iii = permute2f128(xxx, nnn, 0x21);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = loadu4(X+8);
        X += 8;
        vec4 ppp = permute2f128(eee, ddd, 0x21);
        vec4 jjj = permute2f128(nnn, xxx, 0x21);
#ifdef __FMA__
        storeu4(Y, fmadd4(R, fmsub4(two, ppp, add4(eee, ddd)), loadu4(Y)));
#else
        storeu4(Y, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y)));
#endif
        eee = sub4(sub4(xxx, jjj), sub4(jjj, nnn));
        ppp = permute2f128(ddd, eee, 0x21);
#ifdef __FMA__
        storeu4(Y+4, fmadd4(R, fmsub4(two, ppp, add4(ddd, eee)), loadu4(Y+4)));
#else
        storeu4(Y+4, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y+4)));
#endif
        Y += 8;
    }

    // process data 4 by 4:
    if ( Y < end+4 )
    {
        vec4 nnn = loadu4(X+4);
        X += 4;
        vec4 iii = permute2f128(xxx, nnn, 0x21);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = nnn;
        vec4 ppp = permute2f128(eee, ddd, 0x21);
#ifdef __FMA__
        storeu4(Y, fmadd4(R, fmsub4(two, ppp, add4(eee, ddd)), loadu4(Y)));
#else
        storeu4(Y, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y)));
#endif
        eee = ddd;
        Y += 4;
    }

    // process data 2 by 2 using SSE instructions:
    vec2 nn = gethi(xxx);
    vec2 oo = sub2(nn, getlo(xxx));
    vec2 ee = gethi(eee);
    vec2 yy = fmsub2(getlo(two), ee, getlo(eee));
    while ( Y < end+8 )
    {
        vec2 mm = loadu2(X+4);
        X += 2;
        vec2 ff = sub2(mm, nn);
        vec2 dd = sub2(ff, oo);
        nn = mm;
        oo = ff;
        storeup(Y, fmadd2(getlo(R), sub2(yy, dd), loadu2(Y)));
#ifdef __FMA__
        yy = fmsub2(getlo(two), dd, ee);
#else
        yy = sub2(add2(dd, dd), ee);
#endif
        ee = dd;
        Y += 2;
    }
    storeup(Y  ,  fmadd2(getlo(R), yy, loadu2(Y  )));
    storeup(Y+2, fnmadd2(getlo(R), ee, loadu2(Y+2)));
}
#endif
#endif

/*
 This is an optimized implementation
 */
void add_rigidityF(const size_t nbt, const real* X, const real R1, real* Y)
{
    assert_true(nbt > DIM);
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;
    const real R6 = R1 * 6;
    
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

/**
 Add rigidity terms between three points {A, B, C}
 Done with Serge DMITRIEFF, 2015
 */
void add_rigidity(size_t A, size_t B, size_t C, const real* X, const real R1, real* Y)
{
#if ( DIM > 1 )
    for ( size_t d = 0; d < DIM; ++ d )
    {
        real x = 2*X[B*DIM+d] - ( X[A*DIM+d] + X[C*DIM+d] );
        Y[A*DIM+d] += x * R1;
        Y[B*DIM+d] -= x * (R1+R1);
        Y[C*DIM+d] += x * R1;
    }
#endif
}

//------------------------------------------------------------------------------

#define CHECK_RIGIDITY 0

/**
 calculates the second-derivative of point's coordinates,
 scale by the rigidity term, and add to vector Y
*/
void Mecafil::addRigidity(const real* X, real* Y) const
{
#if CHECK_RIGIDITY
    // compare to default implementation:
    real * tmp = new_real(DIM*nPoints);
    copy_real(DIM*nPoints, Y, tmp);
    add_rigidity0(nbt, X, iRigidity, tmp);
#endif
    if ( nPoints > 3 )
    {
        size_t nbt = DIM * ( nPoints - 2 );  // number of triplet values

#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
        add_rigidityF(nbt, X, iRigidity, Y);
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__SSE3__)
        add_rigidity_SSE(nbt, X, iRigidity, Y);
#elif ( DIM > 1 )
        add_rigidityF(nbt, X, iRigidity, Y);
#endif
    
#if NEW_FIBER_LOOP
        if ( iRigidityLoop )
        {
            /*
             With Serge DMITRIEFF:
             Link first and last point in the same way as all other points,
             making the fiber mechanically homogeneous and all points equivalent
             */
            const size_t L = lastPoint();
            add_rigidity(L,   0, 1, X, iRigidity, Y);
            add_rigidity(L-1, L, 0, X, iRigidity, Y);
        }
#endif
    }
    else if ( nPoints > 2 )
    {
        add_rigidity(0, 1, 2, X, iRigidity, Y);
    }
    
#if CHECK_RIGIDITY
    static size_t cnt = 0;
    real err = blas::max_diff(DIM*nPoints, tmp, Y);
    if ( err > 1.0e-6 || ++cnt > 1<<14 )
    {
        cnt = 0;
        printf("addRigidity(%u) error %e\n", nPoints, err);
    }
    free_real(tmp);
#endif
}

