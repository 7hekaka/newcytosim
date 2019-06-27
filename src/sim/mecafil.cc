// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "sim.h"
#include "mecafil.h"
#include "cblas.h"
#include "clapack.h"
#include "matsym.h"
#include "random.h"
#include "fiber_prop.h" // needed for NEW_FIBER_LOOP
//#include "vecprint.h"

//------------------------------------------------------------------------------
Mecafil::Mecafil()
{
    buildProjection();
    rfDragPoint = 0;
    rfRigidity  = 0;
    rfDiff = nullptr;
    rfLag  = nullptr;
    rfLLG  = nullptr;
    rfVTP  = nullptr;
#if NEW_ANISOTROPIC_FIBER_DRAG
    rfDir  = 0;
#endif
    useProjectionDiff = false;
}


Mecafil::~Mecafil()
{
    destroyProjection();
    free_real(rfDiff);
    rfDiff = nullptr;
    rfLag  = nullptr;
    rfLLG  = nullptr;
    rfVTP  = nullptr;
}


//------------------------------------------------------------------------------
Mecafil::Mecafil(Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Fiber");
}


Mecafil& Mecafil::operator=(Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Fiber");
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
        free_real(rfDiff);
        
#if NEW_ANISOTROPIC_FIBER_DRAG
        rfDiff = new_real(ms*(4*DIM+1));
        rfLag  = rfDiff + ms*DIM;
        rfLLG  = rfLag + ms;
        rfDir  = rfLLG + ms*DIM;
        rfVTP  = rfDir + ms*DIM;
#else
        rfDiff = new_real(ms*(2*DIM+1));
        rfLag  = rfDiff + ms*DIM;
        rfLLG  = rfLag + ms;
#endif
        
        // reset Lagrange multipliers
        zero_real(ms, rfLag);
    }
    return ms;
}

void Mecafil::release()
{
    free_real(rfDiff);
    rfDiff = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The argument should be: sc = kT / dt;
 */
real Mecafil::addBrownianForces(real const* rnd, real sc, real* rhs) const
{
    real b = sqrt( 2 * sc * rfDragPoint );

    for ( unsigned jj = 0; jj < DIM*nPoints; ++jj )
        rhs[jj] += b * rnd[jj];
    
    return b / rfDragPoint;
}


//------------------------------------------------------------------------------

/**
 Calculate the normalized difference between successive vertices of the fiber:

     for ( int n = 0; n < DIM*lastPoint(); ++n )
         rfDiff[n] = ( pPos[n+DIM] - pPos[n] ) / segmentation();

 */

void Mecafil::storeDirections()
{
#if ( 1 )
    //checkSegmentation(0.01);
    /*
     assume here that successive points are correctly separated, which is usally
     not the case, but the error is usually small
     */
    const real sc  = 1.0 / segmentation();
    const unsigned end = DIM * lastPoint();
    for ( unsigned p = 0; p < end; ++p )
        rfDiff[p] = sc * ( pPos[p+DIM] - pPos[p] );
#else
    for ( unsigned p = 0; p < lastPoint(); ++p )
        normalize(diffPoints(p)).store(rfDiff+DIM*p);
#endif
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    /*
     Calculate the average filament direction at each vertex of the fiber in rfDir[].
     Note that rfDir[] is calculated here from rfDiff[]
     */
    
    // for the extremities, the direction of the nearby segment is used.
    for ( unsigned d = 0; d < DIM; ++d )
    {
        rfDir[d]     = rfDiff[d];
        rfDir[d+end] = rfDiff[d+end-DIM];
    }
    
    // for intermediate points, the directions of the two flanking segments are averaged
    for ( unsigned p = DIM ; p < end; ++p )
        rfDir[p] = 0.5 * ( rfDiff[p-DIM] + rfDiff[p] );

    //VecPrint::print(std::clog, last+DIM, rfDir);
#endif
}


#if NEW_ANISOTROPIC_FIBER_DRAG

/**
 This will perform:

     Y = X + (TT') X

 Where T is the local direction of the fiber given by `dir`.
 This is used to multiply the tangential component of X by a factor 2,
 without changing the orthogonal components
 */
void scaleTangentially(unsigned nbp, const real* X, const real* dir, real * Y)
{
    for ( unsigned p = 0; p < nbp; ++p )
    {
        real const* xxx = X   + DIM * p;
        real const* ddd = dir + DIM * p;
        real      * yyy = Y   + DIM * p;
#if ( DIM == 2 )
        real s = xxx[0] * ddd[0] + xxx[1] * ddd[1];
        yyy[0] = xxx[0] + s * ddd[0];
        yyy[1] = xxx[1] + s * ddd[1];
#elif ( DIM >= 3 )
        real s = xxx[0] * ddd[0] + xxx[1] * ddd[1] + xxx[2] * ddd[2];
        yyy[0] = xxx[0] + s * ddd[0];
        yyy[1] = xxx[1] + s * ddd[1];
        yyy[2] = xxx[2] + s * ddd[2];
#endif
    }
}
#endif

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

void Mecafil::setSpeedsFromForces(const real* X, real alpha, real* Y) const
{
    real sum = X[0];
    for ( unsigned int ii = 1; ii < nPoints; ++ii )
        sum += X[ii];
    
    sum = alpha * sum / ( rfDragPoint * nPoints );
#if NEW_ANISOTROPIC_FIBER_DRAG
    sum *= 2;
#endif
    for ( unsigned int ii = 0; ii < nPoints; ++ii )
        Y[ii] = sum;
}

void Mecafil::computeTensions(const real*) {} //DIM == 1
void Mecafil::storeTensions(const real*) {} //DIM == 1
void Mecafil::makeProjectionDiff(const real*) {} //DIM == 1
void Mecafil::addProjectionDiff(const real*, real*) const {} //DIM == 1

#endif


//-----------------------------------------------------------------------
#pragma mark -

#if ( 0 )

/**
 Set rigidity terms with modulus 'R1' in lower part of `mat`,
 for a filament with 'cnt' points (and thus with 'cnt/DIM-2' triplets).
 */
template<typename MATRIX>
void add_rigidity_matrix(const int cnt, MATRIX& mat, const int inx, const real R1, const int dim)
{
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;

    const int end = dim * ( inx + cnt - 1 );
    
    for ( int i = dim*(inx+1); i < end; ++i )
    {
        mat(i-dim, i-dim) -= R1;
        mat(i    , i-dim) += R2;
        mat(i+dim, i-dim) -= R1;
        mat(i    , i    ) -= R4;
        mat(i+dim, i    ) += R2;
        mat(i+dim, i+dim) -= R1;
    }
}

#else

/**
 Set rigidity terms with modulus 'R1' in upper diagonal of `mat`,
 for a filament with 'cnt' points (and thus with 'cnt-2' triplets).
 */
template<typename MATRIX>
void add_rigidity_matrix(const int cnt, MATRIX& mat, const int inx, const real R1, const int dim)
{
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;
    const real R5 = R1 * 5;
    const real R6 = R1 * 6;

    const int D = dim, T = dim*2, V = dim*3;
    const int s = dim * inx;
    const int e = s + dim * ( cnt - 2 );

    mat(s  , s  ) -= R1;
    mat(s+D, s  ) += R2;
    mat(s+T, s  ) -= R1;
    
    mat(e+D, e+D) -= R1;
    mat(e+D, e  ) += R2;

    if ( 3 < cnt )
    {
        mat(s+D, s+D) -= R5;
        mat(s+T, s+D) += R4;
        mat(s+V, s+D) -= R1;
        mat(e  , e  ) -= R5;
    }
    else
    {
        mat(s+D, s+D) -= R4;
    }
    
    for ( int n = s+T; n < e ; n += D )
    {
        mat(n,   n) -= R6;
        mat(n+D, n) += R4;
        mat(n+T, n) -= R1;
    }
}

#endif

void Mecafil::addRigidityMatrix(MatrixSparseSymmetric1& mat, const int inx, int dim) const
{
    if ( nPoints > 2 )
    {
        add_rigidity_matrix(nPoints, mat, inx, rfRigidity, dim);
#if ( 1 )
        int N = nPoints;
        MatrixSymmetric m(N);
        add_rigidity_matrix(N, m, 0, 1.0, 1);
        VecPrint::print(std::clog, N, N, m.data(), N, 0);
#endif
    }
}


/**
 Set elements of matrix `mat` corresponding to the elastic terms of the Fiber.
 The array `mat` must be square of dimension `dim * this->nPoints`
 Only terms above the diagonal and corresponding to the first subspace are set
 */
void add_rigidity_upper(unsigned cnt, real* mat, unsigned ldd, const real R1)
{
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;
    const real R5 = R1 * 5;
    const real R6 = R1 * 6;

    constexpr unsigned D = DIM, T = DIM*2, V = DIM*3;
    const unsigned e = DIM * ( cnt - 2 );
    const unsigned f = DIM * ( cnt - 1 );
    
    mat[0      ] -= R1;
    mat[  ldd*D] += R2;
    mat[  ldd*T] -= R1;
    
    mat[e+ldd*f] += R2;
    mat[f+ldd*f] -= R1;
    
    if ( 3 < cnt )
    {
        mat[D+ldd*D] -= R5;
        mat[D+ldd*T] += R4;
        mat[D+ldd*V] -= R1;
        mat[e+ldd*e] -= R5;
    }
    else
    {
        mat[D+ldd*D] -= R4;
    }
    
    for ( unsigned n = T; n < e; n += D )
    {
        mat[n+ldd* n   ] -= R6;
        mat[n+ldd*(n+D)] += R4;
        mat[n+ldd*(n+T)] -= R1;
    }
}


void Mecafil::addRigidityUpper(real * mat, unsigned ldd) const
{
    if ( nPoints > 2 )
    {
        add_rigidity_upper(nPoints, mat, ldd, rfRigidity);
#if ( 0 )
        int N = DIM*nPoints;
        MatrixSymmetric m(N);
        add_rigidity_upper(nPoints, m.data(), ldd, 1.0);
        VecPrint::print(std::clog, N, N, m.data(), N, 0);
#endif
    }
}

//------------------------------------------------------------------------------

/*
 This is the reference implementation
 */
void add_rigidity0(const unsigned nbt, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    for ( unsigned jj = 0; jj < nbt; ++jj )
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
void add_rigidity2(const unsigned nbt, const real* X, const real rigid, real* Y)
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
void add_rigidity3(const unsigned nbt, const real* X, const real rigid, real* Y)
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
    
    real df0, of0 = 0, odf0 = 0;
    real df1, of1 = 0, odf1 = 0;
#if ( DIM >= 3 )
    real df2, of2 = 0, odf2 = 0;
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
void add_rigidity_SSE(const unsigned nbt, const real* X, const real rigid, real* Y)
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
void add_rigidity_AVX(const unsigned nbt, const real* X, const real rigid, real* Y)
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
        storeu2(Y, fmadd2(getlo(R), sub2(yy, dd), loadu2(Y)));
#ifdef __FMA__
        yy = fmsub2(getlo(two), dd, ee);
#else
        yy = sub2(add2(dd, dd), ee);
#endif
        ee = dd;
        Y += 2;
    }
    storeu2(Y  ,  fmadd2(getlo(R), yy, loadu2(Y  )));
    storeu2(Y+2, fnmadd2(getlo(R), ee, loadu2(Y+2)));
}
#endif
#endif

/*
 This is an optimized implementation
 */
void add_rigidityF(const unsigned nbt, const real* X, const real R1, real* Y)
{
    assert_true(nbt > DIM);
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;
    const real R6 = R1 * 6;
    
    const int end = nbt;
    #pragma ivdep
    for ( int i = DIM*2; i < end; ++i )
        Y[i] += R4 * (X[i-DIM]+X[i+DIM]) - R1 * (X[i-DIM*2]+X[i+DIM*2]) - R6 * X[i];
    
    // special cases near the edges:
    real      * Z = Y + nbt + DIM;
    real const* E = X + nbt + DIM;
    #pragma ivdep
    for ( int d = 0; d < DIM; ++d )
    {
        Y[d    ] -= R1 * (X[d+DIM*2]+X[d]) - R2 * X[d+DIM];
        Y[d+DIM] -= R1 * (X[d+DIM]+X[d+DIM*3]) + R4 * (X[d+DIM]-X[d+DIM*2]) - R2 * X[d];
        Z[d-DIM] -= R1 * (E[d-DIM]+E[d-DIM*3]) + R4 * (E[d-DIM]-E[d-DIM*2]) - R2 * E[d];
        Z[d    ] -= R1 * (E[d-DIM*2]+E[d]) - R2 * E[d-DIM];
    }
}

/**
 Add rigidity terms between three points {A, B, C}
 Done with Serge DMITRIEFF, 2015
 */
void add_rigidity(unsigned A, unsigned B, unsigned C, const real* X, const real R1, real* Y)
{
    for ( unsigned d = 0; d < DIM; ++ d )
    {
        real x = 2*X[B*DIM+d] - ( X[A*DIM+d] + X[C*DIM+d] );
        Y[A*DIM+d] += x * R1;
        Y[B*DIM+d] -= x * (R1+R1);
        Y[C*DIM+d] += x * R1;
    }
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
    add_rigidity0(nbt, X, rfRigidity, tmp);
#endif
    if ( nPoints > 3 )
    {
        unsigned nbt = DIM * ( nPoints - 2 );  // number of triplet values

#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
        add_rigidityF(nbt, X, rfRigidity, Y);
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__SSE3__)
        add_rigidity_SSE(nbt, X, rfRigidity, Y);
#else
        add_rigidityF(nbt, X, rfRigidity, Y);
#endif
    
#if NEW_FIBER_LOOP
        if ( rfRigidityLoop )
        {
            /*
             With Serge DMITRIEFF:
             Link first and last point in the same way as all other points,
             making the fiber mechanically homogeneous and all points equivalent
             */
            const unsigned L = lastPoint();
            add_rigidity(L,   0, 1, X, rfRigidity, Y);
            add_rigidity(L-1, L, 0, X, rfRigidity, Y);
        }
#endif
    }
    else if ( nPoints > 2 )
    {
        add_rigidity(0, 1, 2, X, rfRigidity, Y);
    }
    
#if CHECK_RIGIDITY
    static int cnt = 0;
    real err = blas::max_diff(DIM*nPoints, tmp, Y);
    if ( err > 1.0e-6 || ++cnt > 1<<14 )
    {
        cnt = 0;
        printf("addRigidity(%u) error %e\n", nPoints, err);
    }
    free_real(tmp);
#endif
}

