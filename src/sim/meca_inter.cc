// Cytosim was created by Francois Nedelec. Copyright 2019 Cambridge University

#include "dim.h"
#include "cymdef.h"
#include "space.h"
#include "modulo.h"
#include "mecapoint.h"
#include "interpolation.h"
#include "matrix11.h"
#include "matrix22.h"
#include "matrix33.h"
#include "matrix44.h"

//#include "vecprint.h"

/// set TRUE to update matrix mFUL using block directives
/** This is significantly faster on machine with the AVX instruction set */
#define USE_MATRIX_BLOCK 1

/// use Mecapoint::pos() instead of interpolating vPTS[]
#define USE_GLOBAL_POSITION 0


#if DRAW_MECA_LINKS
#  include "meca_inter_draw.cc"
#else
#  define DRAW_LINK(...) ((void) 0)
#endif

//------------------------------------------------------------------------------
#pragma mark - Accessory functions

/// true if 'a' is equal to 'b' or 'c'
static inline bool any_equal(const index_t a, const index_t b, const index_t c)
{
    return ( a == c ) | ( a == b );
}

/// true if 'a' is equal to 'b', 'c' or 'd'
static inline bool any_equal(const index_t a, const index_t b,
                             const index_t c, const index_t d)
{
    return ( a == b ) | ( a == c ) | ( a == d );
}

/// true if 'a' is equal to 'b', 'c', 'd' or 'e'
static inline bool any_equal(const index_t a, const index_t b, const index_t c,
                             const index_t d, const index_t e)
{
    return ( a == b ) | ( a == c ) | ( a == d ) | ( a == e );
}

//------------------------------------------------------------------------------
#pragma mark - Linear Interpolation of Vectors
//------------------------------------------------------------------------------

static Vector interpolate1(const real vec[], const index_t inx)
{
    return Vector(vec+DIM*inx);
}

static Vector interpolate2(const real vec[], const index_t inx, real a, real b)
{
    Vector P0(vec+DIM*inx);
    Vector P1(vec+DIM*(inx+1));
    return a * P0 + b * P1;
}

static Vector interpolate3(const real vec[], const index_t inx, const real a, real b, real c)
{
    Vector P0(vec+DIM*inx);
    Vector P1(vec+DIM*(inx+1));
    Vector P2(vec+DIM*(inx+2));
    return ( a * P0 + b * P1 ) + c * P2;
}

static Vector interpolate4(const real vec[], const index_t inx, const real a, real b, real c, real d)
{
    Vector P0(vec+DIM*inx);
    Vector P1(vec+DIM*(inx+1));
    Vector P2(vec+DIM*(inx+2));
    Vector P3(vec+DIM*(inx+3));
    return ( a * P0 + b * P1 ) + ( c * P2 + d * P3 );
}

//------------------------------------------------------------------------------

template < typename T, typename U >
static inline Vector position_diff(T const& A, U const& B)
{
    return B.pos() - A.pos();
}

template < typename T, typename U >
static inline Vector modulo_offset(T const& A, U const& B)
{
    return modulo->offset(position_diff(A, B));
}

//------------------------------------------------------------------------------
#pragma mark - Functions to set matrix elements
//------------------------------------------------------------------------------

#if 0
#  define CHECK_INDEX(I) ((void) 0)
#  define CHECK_INDICES(I,J,C) ((void) 0)
#else
#  define CHECK_INDEX(I) { assert_true(I<nPoints_); }
#  define CHECK_INDICES(I,J,C) { if (J>I) printf(" wrong-sided %s %u %u\n",C,I,J); }
#endif

#define PRINT_BLOCK(I,J,B) ((void) 0)
//#define PRINT_BLOCK(I,J,B) { std::clog<<std::setw(3)<<I<<" "<<std::setw(3)<<J<<" "<<std::setprecision(2)<<std::setw(10)<<B<<'\n'; }


/// add `alpha * T` to the matrix.
inline void Meca::add_block(index_t i, index_t j, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"add");
#if USE_MATRIX_BLOCK
    mFUL.block(i, j).add_full(T);
    PRINT_BLOCK(i,j,T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.element(i,j) += T.value();
#else
    assert_true( i > j );
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = 0; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*j+x) += T(y,x);
#endif
}

/// add `alpha * T` to the matrix.
inline void Meca::add_block(index_t i, index_t j, const real alpha, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"add_alpha");
#if USE_MATRIX_BLOCK
    mFUL.block(i, j).add_full(alpha, T);
    PRINT_BLOCK(i,j,alpha*T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.element(i,j) += alpha * T.value();
#else
    assert_true( i > j );
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = 0; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*j+x) += alpha * T(y,x);
#endif
}

/// subtract `T` to the matrix.
inline void Meca::sub_block(index_t i, index_t j, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"sub");
#if USE_MATRIX_BLOCK
    mFUL.block(i, j).sub_full(T);
    PRINT_BLOCK(i,j,T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.element(i,j) -= T.value();
#else
    assert_true( i > j );
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = 0; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*j+x) -= T(y,x);
#endif
}

/// subtract `alpha * T` to the matrix.
inline void Meca::sub_block(index_t i, index_t j, const real alpha, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"sub_alpha");
#if USE_MATRIX_BLOCK
    mFUL.block(i, j).sub_full(alpha, T);
    PRINT_BLOCK(i,j,alpha*T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.element(i,j) -= alpha * T.value();
#else
    assert_true( i > j );
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = 0; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*j+x) -= alpha * T(y,x);
#endif
}

/// add `T` to the matrix diagonal. `T` should be symmetric
inline void Meca::add_block_diag(index_t i, MatrixBlock const& T)
{
    CHECK_INDEX(i);
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mFUL.diag_block(i).add_half(T);
    PRINT_BLOCK(i,i,T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.diagonal(i) += T.value();
#else
    // add lower part of block
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = x; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*i+x) += T(y,x);
#endif
}

/// add `alpha * T` to the matrix diagonal. `T` should be symmetric
inline void Meca::add_block_diag(index_t i, const real alpha, MatrixBlock const& T)
{
    CHECK_INDEX(i);
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mFUL.diag_block(i).add_half(alpha, T);
    PRINT_BLOCK(i,i,alpha*T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.diagonal(i) += alpha * T.value();
#else
    // add lower part of block
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = x; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*i+x) += alpha * T(y,x);
#endif
}

/// add `alpha * T` to the matrix diagonal. `T` should be symmetric
inline void Meca::add_block_diag(index_t i, const real alpha, MatrixBlock const& T, const real dia)
{
    CHECK_INDEX(i);
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mFUL.diag_block(i).add_half(alpha, T, dia);
    PRINT_BLOCK(i,i,alpha*T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.diagonal(i) += alpha * ( T.value() + dia );
#else
    // add lower part of block
    for ( size_t x = 0; x < DIM; ++x )
    {
        mFUL(DIM*i+y, DIM*i+x) += alpha * ( T(y,x) + dia );
        for ( size_t y = x+1; y < DIM; ++y )
            mFUL(DIM*i+y, DIM*i+x) += alpha * T(y,x);
    }
#endif
}

/// add `-T` to the matrix diagonal. `T` should be symmetric
inline void Meca::sub_block_diag(index_t i, MatrixBlock const& T)
{
    CHECK_INDEX(i);
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mFUL.diag_block(i).sub_half(T);
    PRINT_BLOCK(i,i,T);
#elif ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.diagonal(i) -= T.value();
#else
    // add lower part of block
    for ( index_t x = 0; x < DIM; ++x )
    for ( index_t y = x; y < DIM; ++y )
        mFUL(DIM*i+y, DIM*i+x) -= T(y,x);
#endif
}


/// add `val` to the XYZ-isometric matrix
inline void Meca::add_iso(index_t i, index_t j, const real val)
{
    CHECK_INDICES(i,j,"add_iso");
#if USE_ISO_MATRIX
    mISO.element(i,j) += val;
#elif USE_MATRIX_BLOCK
    mFUL.block(i, j).add_diag(val);
#else
    for ( index_t x = 0; x < DIM; ++x )
        mFUL(DIM*i+x, DIM*j+x) += val;
#endif
}

/// add `-val` to the XYZ-isometric matrix
inline void Meca::sub_iso(index_t i, index_t j, const real val)
{
    CHECK_INDICES(i,j,"sub_iso");
#if USE_ISO_MATRIX
    mISO.element(i,j) -= val;
#elif USE_MATRIX_BLOCK
    mFUL.block(i, j).add_diag(-val);
#else
    for ( index_t x = 0; x < DIM; ++x )
        mFUL(DIM*i+x, DIM*j+x) -= val;
#endif
}

/// add `val` to the XYZ-isometric matrix
inline void Meca::add_iso_diag(index_t i, const real val)
{
    CHECK_INDEX(i);
#if USE_ISO_MATRIX
    mISO.diagonal(i) += val;
#elif USE_MATRIX_BLOCK
    mFUL.diag_block(i).add_diag(val);
#else
    for ( index_t x = 0; x < DIM; ++x )
        mFUL(DIM*i+x, DIM*i+x) += val;
#endif
}

/// add `-val` to the XYZ-isometric matrix
inline void Meca::sub_iso_diag(index_t i, const real val)
{
    CHECK_INDEX(i);
#if USE_ISO_MATRIX
    mISO.diagonal(i) -= val;
#elif USE_MATRIX_BLOCK
    mFUL.diag_block(i).add_diag(-val);
#else
    for ( index_t x = 0; x < DIM; ++x )
        mFUL(DIM*i+x, DIM*i+x) -= val;
#endif
}


/// add `vec` to the base
inline void Meca::add_base(index_t i, Vector const& vec) const
{
    CHECK_INDEX(i);
    vec.add_to(vBAS+DIM*i);
}

/// add `alpha * vec` to the base
inline void Meca::add_base(index_t i, Vector const& vec, const real alpha) const
{
    CHECK_INDEX(i);
    vec.add_to(alpha, vBAS+DIM*i);
}

/// add `-vec` to the base
inline void Meca::sub_base(index_t i, Vector const& vec) const
{
    CHECK_INDEX(i);
    vec.sub_to(vBAS+DIM*i);
}

/// add `-alpha * vec` to the base
inline void Meca::sub_base(index_t i, Vector const& vec, const real alpha) const
{
    CHECK_INDEX(i);
    vec.sub_to(alpha, vBAS+DIM*i);
}

//------------------------------------------------------------------------------
#pragma mark - Explicit (constant) Forces
//------------------------------------------------------------------------------

/**
 Add constant force to a vertex
 */
void Meca::addForce(Mecapoint const& pte, Vector const& force)
{
    add_base(pte.matIndex0(), force);
}

/**
 Add constant force to a vertex
 */
void Meca::addForce(Mecable const* mec, const index_t inx, Vector const& force)
{
    add_base(mec->matIndex() + inx, force);
}


/**
Add constant force to an interpolated position
 */
void Meca::addForce(Interpolation const& pti, Vector const& force)
{
    add_base(pti.matIndex1(), force, pti.coef0());
    add_base(pti.matIndex2(), force, pti.coef1());
}


void Meca::addForceToAll(Vector const& force)
{
    for ( index_t i = 0; i < nbVertices(); ++i )
        force.add_to(vBAS+DIM*i);
}

//------------------------------------------------------------------------------
#pragma mark - Explicit (constant) Torque
//------------------------------------------------------------------------------

/**
 Add constant torque in `pti`:
 
     force = cross(torque, position)
 
 This is explicit and all contributions go in the force vector vBAS[]
*/
void Meca::addTorque(Interpolation const& pti, const Torque & torque)
{
    Vector d = pti.diff();
    Vector f = cross(torque/d.normSqr(), d);
    
    sub_base(pti.matIndex1(), f);
    add_base(pti.matIndex2(), f);
}


/**
 Add an explicit torque to constrain a segment in a given direction `dir`,
 with a given weight:
 
     torque = weigth * cross(normalize(segment), dir)
     force = cross(torque, position)
 
 This code assumes norm(dir) == 1
 This is explicit and all contributions go in the force vector vBAS[]
 \todo update addTorqueClamp to implicit form
*/
void Meca::addTorqueClamp(Interpolation const& pti,
                          Vector const& dir,
                          const real weight)
{
    assert_true( weight >= 0 );
    
    Vector d = pti.diff();
    real n = d.normSqr();

    Torque Tq = cross(d, dir);

#if ( DIM >= 3 )
    const real Tn = Tq.norm();
#else
    const real Tn = abs_real(Tq);
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = std::atan2(Tn, dot(d, dir));
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     knowing that Tn = Tq.norm() = sin(angle) * sqrt(n)
     
     To have a Torque proportional to sin(angle), use:
     real nn = weight / ( n * sqrt(n) );
     */
    real nn = weight * angle / ( n * Tn );
    
    Vector f = cross(Tq * nn, d);
    
    sub_base(pti.matIndex1(), f);
    add_base(pti.matIndex2(), f);
}


/**
 Add an explicit torque to bring two segments parallel to each other,
 with a given weight:
 
     torque = weigth * cross(dirA, dirB)
     forceA =  cross(torque, dirA)
     forceB = -cross(torque, dirB)
 
 This is explicit and all contributions go in the force vector vBAS[]
 */
void Meca::addTorqueExplicit(Interpolation const& ptA,
                             Interpolation const& ptB,
                             const real weight)
{
    assert_true( weight >= 0 );
    
    Vector da = ptA.diff();
    Vector db = ptB.diff();

    real na = da.normSqr();
    real nb = db.normSqr();

    Torque Tq = cross(da, db);
    
#if ( 0 )
    
#if ( DIM >= 3 )
    const real Tn = Tq.norm();
#else
    const real Tn = abs_real(Tq);
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = std::atan2(Tn, dot(da, db));
    
    /*
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     knowing that Tn = Tq.norm() = sin(angle) * sqrt( na * nb )
     */
    na = weight * angle / ( na * Tn );
    nb = weight * angle / ( nb * Tn );
#else
    // To have a Torque proportional to sin(angle)
    real nn = std::sqrt( na * nb ); //or nn = Tn / angle
    na = weight / ( na * nn );
    nb = weight / ( nb * nn );
    
#endif
    
    Vector fa = cross(Tq * na, da);
    Vector fb = cross(db, Tq * nb);

    sub_base(ptA.matIndex1(), fa);
    add_base(ptA.matIndex2(), fa);
    sub_base(ptB.matIndex1(), fb);
    add_base(ptB.matIndex2(), fb);
}


/**
 Add an explicit torque to induce two segments to make an angle
 defined by (cosine, sine) relative to each other:
 
     torque = weigth * cross( dirA , dirB.rotated(angle) )
     forceA =  cross(torque, dirA)
     forceB = -cross(torque, dirB)
 
 The direction of `ptB` is rotated around `axis` defined as cross( dirA, dirB ).
 The calculation is explicit and all contributions go in the force vector vBAS[]
 It is assumed that `ang.norm() = 1`
 Note that if ( sine == 0 ), you can use addTorque(ptA, ptB, weight)
 */
void Meca::addTorqueExplicit(Interpolation const& ptA,
                             Interpolation const& ptB,
                             Vector2 const& ang, const real weight)
{
    assert_true( weight >= 0 );
    assert_small( ang.normSqr() - 1.0 );
    
    Vector da = ptA.diff();
    Vector db = ptB.diff();

    real na = da.normSqr();
    real nb = db.normSqr();
    
#if ( DIM >= 3 )
    
    /*
     in 3D the axis of torque is perpendicular to both `da` and `db`,
     and the angle is only defined between 0 and PI,
     */
    Vector axis = cross(db, da).normalized(sign_real(ang.YY));
    
    // rotate vector `db` around `arm` by angle specified as (cosine, sine):
    Vector rot = ang.XX * db + ang.YY * cross(axis, db);
    
#elif ( DIM == 2 )

    // this correspond to the Z-direction, up or down:
    real dir = sign_real(cross(da, db));

    // rotate vector `db` by angle defined by (cosine, sine) around Z
    Vector rot( db.XX*ang.XX + db.YY*ang.YY*dir, db.YY*ang.XX - db.XX*ang.YY*dir );
    
#else
    
    // this is meaningless but makes compilation possible
    Vector rot(0.0);
    
    throw InvalidParameter("Meca::addTorque is meaningless in 1D");

#endif

    // calculate torque by vector-product:
    Torque Tq = cross(da, rot);
    
#if ( DIM >= 3 )
    const real Tn = Tq.norm();
#else
    const real Tn = abs_real(Tq);
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = std::atan2(Tn, dot(da, rot));
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     but knowing that Tn = Tq.norm() = sin(angle) * sqrt( na * nb )
     
     To have a Torque proportional to sin(angle), use:
     real nn = sqrt( na * nb ); or nn = Tn / angle;
     na = weight / ( na * nn );
     nb = weight / ( nb * nn );
     */
    na = weight * angle / ( na * Tn );
    nb = weight * angle / ( nb * Tn );
 
    // forces are divided appropriately to reach the desired torque:
    Vector fa = cross(Tq * na, da);
    Vector fb = cross(db, Tq * nb);
    
    // explicit contributions in vBAS
    sub_base(ptA.matIndex1(), fa);
    add_base(ptA.matIndex2(), fa);
    sub_base(ptB.matIndex1(), fb);
    add_base(ptB.matIndex2(), fb);
}


//------------------------------------------------------------------------------
#pragma mark - Implicit Torque
//------------------------------------------------------------------------------

#if ( DIM == 2 )
/**
 Add torque between segments AB and CD containing `pt1` and `pt2`.
 Implicit version with linearized force 2D
 Angle is between AB and CD. Force is along normal N_A and N_C pointing to the other filament
 L_AB and L_CD is the length of the segments AB and CD
 force_A = torque_weight * ( Delta angle ) * N_A/L_AB =-force_B
 force_C = torque_weight * ( Delta angle ) * N_C/L_CD =-force_D
 Delta_angle is the difference between actual angle and resting angle between AB and CD
 
 Antonio Politi, 2013
 
 This code is outdated, and one should use addTorque() instead
 */
void Meca::addTorquePoliti(Interpolation const& pt1,
                           Interpolation const& pt2,
                           Vector const& ang,
                           const real weight)
{
    assert_true( weight >= 0 );
    assert_small( ang.normSqr() - 1.0 );

    if ( pt1.overlapping(pt2) )
        return;
    
    // full indices:
    const size_t idx[] = { DIM*pt1.matIndex1(), DIM*pt1.matIndex1()+1,
                           DIM*pt1.matIndex2(), DIM*pt1.matIndex2()+1,
                           DIM*pt2.matIndex1(), DIM*pt2.matIndex1()+1,
                           DIM*pt2.matIndex2(), DIM*pt2.matIndex2()+1 };
    
    //Vectors and points of torque
    Vector ab = pt1.diff();
    Vector cd = pt2.diff();
    Vector a = pt1.pos1();
    Vector b = pt1.pos2();
    Vector c = pt2.pos1();
    Vector d = pt2.pos2();
    const real coord[]={a.XX, a.YY, b.XX, b.YY, c.XX, c.YY, d.XX, d.YY};
    //Helping vector this vector is at torque_angle from cd.
    //Therefore in resting state angle difference between ab and ce is zero. This vector is used to compute the strength of torque
    Vector ce;
    ce.XX =  cd.XX*ang.XX + cd.YY*ang.YY;
    ce.YY = -cd.XX*ang.YY + cd.YY*ang.XX;
    //normalize
    const real abn = ab.norm();
    const real abnS= ab.normSqr();
    const real cdn = cd.norm();
    const real cdnS= cd.normSqr();
    if (abn < REAL_EPSILON || cdn < REAL_EPSILON ) return;
    
    //normalize the vectors
    ab /= abn; cd /= cdn; ce /= cdn;
    
    //Coordinates of normal vectors yielding the direction of the force
    //fa = torque_weight*dangle*(h[0], h[1]) = torque_weight*dangle*na/la
    const real h[]={ ab.YY/abn, -ab.XX/abn, -ab.YY/abn, ab.XX/abn, -cd.YY/cdn, cd.XX/cdn, cd.YY/cdn, -cd.XX/cdn };
    
    //dangle = angle - torque_angle
    //real dangle = std::atan2( cross(ab, ce), dot(ab, ce) );
    real dangle = std::atan2( ab.XX*ce.YY - ab.YY*ce.XX, dot(ab, ce) );
    //Computation of the jacobian for the linearization
    //M = d_x f = M1 + M2
    //M1 = w1/l normal d_x dangle
    //M2 = w2 * dangle  d_x normal/l
    real w1 = weight;
    real w2 = weight*dangle;
    
    //Matrix M1 with k*hxh (outer product) this yieald a matrix stored with its lower triangular part in m. The -w1 is because ab = b-a
    real m[36] = { 0 };
    //blas::xspr('U', 8, -w1, h, 1, m);
    blas::xspr('L', 8, -w1, h, 1, m);
    
    //Matrix M2
    real Da = w2*( -2*ab.XX*ab.YY )/abnS;
    real da = w2*( ab.XX*ab.XX-ab.YY*ab.YY )/abnS;
    real Dc = w2*( -2*cd.XX*cd.YY )/cdnS;
    real dc = w2*(  cd.XX*cd.XX-cd.YY*cd.YY )/cdnS;
    real entrya[] = {-Da, -da, Da,  da}; //={d(na_x/la)/dxa, d(na_x/la)/dya, d(na_x/l)/dxb, ...}
    real entryc[] = { Dc,  dc,  -Dc,  -dc};//={ d(nc_x/lc)/dxc, d(nc_x/lc)/dyc, d(nc_x/l)/dxd, ...}
    int shifta = 0;
    int shiftc= 26;
    int mm;
    
    //Add second part of matrix.
    //The pos(-1, jj) accounts for the different signs of the matrix
    for ( int jj=0; jj <  4; ++jj) {
        for ( int ii=jj ; ii < 4; ++ii ) {
            m[ii + shifta] += std::pow(-1,jj)*entrya[ii-jj];
            m[ii + shiftc] += std::pow(-1,jj)*entryc[ii-jj];
        }
        shifta += 7 - jj;
        shiftc += 3 - jj;
    }
    
    //very Cumbersome!!!
    //vBAS = fa - M*P0
    for ( int ii = 0; ii < 8; ++ii )
    {
        vBAS[idx[ii]] += w2*h[ii];
        for (int jj = 0; jj < 8; ++jj) {
            if (jj < ii)
                mm = int(jj*(7.5-0.5*jj)+ii);
            else {
                mm = int(ii*(7.5-0.5*ii)+jj);
                mFUL(idx[ii], idx[jj]) += m[mm];
            }
            vBAS[idx[ii]] -= m[mm]*coord[jj];
        }
    }
}
#endif

/**
 Add torque between segments AB and CD defined by `pt1` and `pt2`.
 Opposite forces are applied at the end of the segments resulting in pure torque.
 
     force_A = weight * sin( angle - equilibrium_angle ) / |AB|
     force_B = - force_A
     force_C = weight * sin( angle - equilibrium_angle ) / |CD|
     force_D = - force_C
 
 These force vectors are orthogonal to the segments on which they are applied.
 The equilibrium angle, and the axis of rotation is specified by the rotation matrix R
 The interpolation coefficients of `pt1` and `pt2` are ignored.

 3D implicit torque implementation
 Serge Dmitrieff and FJN, 21.01.2019 -- 11.02.2019
 www.biophysics.fr and Cambridge University
*/
#if ( DIM > 1 )
void Meca::addTorque(Interpolation const& pt1,
                     Interpolation const& pt2,
                     MatrixBlock const& R,
                     const real weight)
{
    assert_true( weight >= 0 );

    const Vector AB = pt1.diff();
    const Vector CD = pt2.diff();
    const real iU = AB.inv_norm();
    const real iV = CD.inv_norm();
    const real wU = weight * iU;
    const real wV = weight * iV;
    const Vector u = AB * iU;
    const Vector v = CD * iV;

    //const MatrixBlock Id(0,1);  // identity matrix

    Vector Ru = R.vecmul(u);
    //Vector Tv = T.vecmul(v);
    Vector Tv = R.trans_vecmul(v);
    
#if ( 1 )
    real Tvu = dot(Tv, u);
    real Ruv = dot(Ru, v);
    assert_small(Tvu-Ruv); // Tvu and Ruv should be equal
    // current forces, exact formula
    Vector Fu = ( Tv - Tvu * u ) * wU;
    Vector Fv = ( Ru - Ruv * v ) * wV;
#else
    // approximate formula, if close to equilibrium: Tvu = Ruv = 1
    Vector Fu = ( Tv - u ) * wU;
    Vector Fv = ( Ru - v ) * wV;
#endif

    //std::clog << std::fixed;
    //std::clog << std::setw(9) << Tvu << " " << std::setw(9) << Ruv << std::setw(9) << sine << "\n";
    //std::clog << "u " << std::setw(12) << u << " Tv " << std::setw(12) << Tv << "\n";
    //std::clog << "v " << std::setw(12) << v << " Ru " << std::setw(12) << Ru << "\n";

#if ( 0 )
    // EXPLICIT
    sub_base(iiA, Fu);        // F(A)=-Fu
    add_base(iiB, Fu);        // F(B)=+Fu
    sub_base(iiC, Fv);        // F(C)=-Fv
    add_base(iiD, Fv);        // F(D)=+Fv
    return;
#endif

    const real wUU = wU * iU;
    const real wUV = wU * iV;
    const real wVV = wV * iV;
    
#if ( 0 )
    /// matrices used during development for testing different formula
    const MatrixBlock uxu = MatrixBlock::outerProduct(u);
    //const MatrixBlock uxv = MatrixBlock::outerProduct(u,v);
    const MatrixBlock vxu = MatrixBlock::outerProduct(v,u);
    const MatrixBlock vxv = MatrixBlock::outerProduct(v);
    //const MatrixBlock uxRu = MatrixBlock::outerProduct(u,Ru);
    const MatrixBlock vxRu = MatrixBlock::outerProduct(v,Ru);
    const MatrixBlock uxTv = MatrixBlock::outerProduct(u,Tv);
    const MatrixBlock vxTv = MatrixBlock::outerProduct(v,Tv);
    const MatrixBlock Ruxu = MatrixBlock::outerProduct(Ru,u);
    const MatrixBlock Ruxv = MatrixBlock::outerProduct(Ru,v);
    const MatrixBlock Tvxu = MatrixBlock::outerProduct(Tv,u);
    //const MatrixBlock Tvxv = MatrixBlock::outerProduct(Tv,v);

    // EXACT 1: formula obtained by derivation:
    const MatrixBlock duFu = (( uxu * 3 - Id ) * Tvu - Tvxu - uxTv ) * wUU;
    //const MatrixBlock duFu = MatrixBlock::offsetOuterProduct(-Tvu*wUU, u, 3*Tvu*wUU) - ( Tvxu + uxTv ) * wUU;

    //const MatrixBlock dvFu = ( T + uxv * Tvu - uxRu - Tvxv ) * wUV;
    const MatrixBlock duFv = ( R + vxu * Ruv - Ruxu - vxTv ) * wUV;
    const MatrixBlock dvFv = (( vxv * 3 - Id ) * Ruv - Ruxv - vxRu ) * wVV;
    //const MatrixBlock dvFv = MatrixBlock::offsetOuterProduct(-Ruv*wVV, v, 3*Ruv*wVV) - ( Ruxv + vxRu ) * wUU;
#else
    /*
    const MatrixBlock TvxTv = MatrixBlock::outerProduct(Tv);
    const MatrixBlock RuxRu = MatrixBlock::outerProduct(Ru);
    const MatrixBlock TvxRu = MatrixBlock::outerProduct(Tv,Ru);
    const MatrixBlock RuxTv = MatrixBlock::outerProduct(Ru,Tv);
     */
 /*
    // derivatives taken at u->Tv and v:
    const MatrixBlock duFu = ( TvxTv * ( 3 * Tvu - 2 ) - Id * Tvu ) * wUU;
    const MatrixBlock dvFu = ( T + Tvxv * ( Tvu - 2 )) * wUV;
    const MatrixBlock duFv = ( R + vxTv * ( Ruv - 2 )) * wUV;
    const MatrixBlock dvFv = ( vxv   * ( 3 * Ruv - 2 ) - Id * Ruv ) * wVV;
*/
/*
    // derivatives at u and v->Ru:
    const MatrixBlock duFu = ( uxu   * ( 3 * Tvu - 2 ) - Id * Ruv ) * wUU;
    const MatrixBlock dvFu = ( T + uxRu * ( Tvu - 2 )) * wUV;
    const MatrixBlock duFv = ( R + Ruxu * ( Ruv - 2 )) * wUV;
    const MatrixBlock dvFv = ( RuxRu * ( 3 * Ruv - 2 ) - Id * Tvu ) * wVV;
 */
/*
    // SIMPLIFIED 2: combining the two formula and assuming that ( Tvu = Ruv = 1 ):
    // This works great!
    const MatrixBlock duFu = ( uxu + TvxTv - Id * 2 ) * ( wUU * 0.5 );
    const MatrixBlock dvFu = ( T * 2 - uxRu - Tvxv ) * ( wUV * 0.5 );
    const MatrixBlock duFv = ( R * 2 - Ruxu - vxTv ) * ( wUV * 0.5 );
    const MatrixBlock dvFv = ( RuxRu + vxv - Id * 2 ) * ( wVV * 0.5 );

    // near equilibrium, further assuming that ( Tv = u ) and ( Ru = v )
    const MatrixBlock duFu = ( uxu - Id ) * wUU;
    //const MatrixBlock dvFu = ( T - uxv ) * wUV;
    const MatrixBlock duFv = ( R - vxu ) * wUV;
    const MatrixBlock dvFv = ( vxv - Id ) * wVV;
*/

    // SIMPLIFIED 3: same as above, but directly setting the matrices
    const MatrixBlock duFu = MatrixBlock::offsetOuterProduct(-wUU, u, wUU);
    //const MatrixBlock dvFu = ( T - MatrixBlock::outerProduct(u,v) ) * wUV;
    const MatrixBlock duFv = ( R - MatrixBlock::outerProduct(v,u) ) * wUV;
    const MatrixBlock dvFv = MatrixBlock::offsetOuterProduct(-wVV, v, wVV);
/*
    // SIMPLIFIED 4: why not just put the formula in?
    // Fu = ( Tv - u * Tvu ) * wU
    // Fv = ( Ru - v * Ruv ) * wV
    const MatrixBlock duFu = Id * ( -wUU * Tvu );
    const MatrixBlock dvFu = T * wUV;
    const MatrixBlock duFv = R * wUV;
    const MatrixBlock dvFv = Id * ( -wVV * Ruv );
 */

#endif
    
    // indices:
    const index_t iiA = pt1.matIndex1();
    const index_t iiB = pt1.matIndex2();
    const index_t iiC = pt2.matIndex1();
    const index_t iiD = pt2.matIndex2();

    add_block_diag(iiA, duFu);
    sub_block(iiB, iiA, duFu);
    add_block_diag(iiB, duFu);
    if ( iiA < iiC )
    {
        add_block(iiC, iiA, duFv);
        sub_block(iiD, iiA, duFv);
        sub_block(iiC, iiB, duFv);
        add_block(iiD, iiB, duFv);
    }
    else
    {
        const MatrixBlock dvFu = duFv.transposed();
        add_block(iiA, iiC, dvFu);
        sub_block(iiA, iiD, dvFu);
        sub_block(iiB, iiC, dvFu);
        add_block(iiB, iiD, dvFu);
    }
    add_block_diag(iiC, dvFv);
    sub_block(iiD, iiC, dvFv);
    add_block_diag(iiD, dvFv);
    
    // remaining part of the force (should be small near equilibrium)
    Vector Fu0 = Fu - duFu.vecmul(AB) - duFv.trans_vecmul(CD);
    Vector Fv0 = Fv - duFv.vecmul(AB) - dvFv.vecmul(CD);

    // add constant terms into vBAS[]:
    sub_base(iiA, Fu0);        // F(A) = -Fu
    add_base(iiB, Fu0);        // F(B) = +Fu
    sub_base(iiC, Fv0);        // F(C) = -Fv
    add_base(iiD, Fv0);        // F(D) = +Fv

#if ( 0 )
    const MatrixBlock dvFu = duFv.transposed();
    //std::clog << iiA << " " << iiB << " " << iiC << " " << iiD << " " << sine << "\n";
    std::clog << "duFu " << std::setw(12) << duFu << " duFv " << std::setw(12) << duFv << "\n";
    std::clog << "dvFu " << std::setw(12) << dvFu << " dvFv " << std::setw(12) << dvFv << "\n";
    std::clog << "Fu0 " << std::fixed << std::setw(12) << Fu0 << "    Fv0 " << std::setw(12) << Fv0 << "\n";
#endif
}


void Meca::addTorque(Interpolation const& pt1,
                     Interpolation const& pt2,
                     Vector2 const& ang, const real weight)
{
    assert_true( weight >= 0 );
    assert_small( ang.normSqr() - 1.0 );
    
#if ( DIM >= 3 )
    const Vector AB = pt1.diff();
    const Vector CD = pt2.diff();
    Vector axis = cross(AB, CD);
    real n = axis.norm();
    if ( n > REAL_EPSILON )
        axis /= n;
    else
        axis = Vector::randU();
    MatrixBlock rot = MatrixBlock::rotationAroundAxis(axis, ang.XX, ang.YY);
    addTorque(pt1, pt2, rot, weight);
#elif ( DIM == 2 )
    MatrixBlock rot(ang.XX, ang.YY, -ang.YY, ang.XX);
    addTorque(pt1, pt2, rot, weight);
#endif
}


/// derived from above, with R = Identity, 22.06.2024
void Meca::addTorque(Interpolation const& pt1,
                     Interpolation const& pt2,
                     const real weight)
{
    assert_true( weight >= 0 );

    const Vector AB = pt1.diff();
    const Vector CD = pt2.diff();
    const real iU = AB.inv_norm();
    const real iV = CD.inv_norm();
    const real wU = weight * iU;
    const real wV = weight * iV;
    const Vector u = AB * iU;
    const Vector v = CD * iV;

    real uv = dot(v, u);
    // current forces, exact formula
    Vector Fu = ( v - uv * u ) * wU;
    Vector Fv = ( u - uv * v ) * wV;

    const real wUU = wU * iU;
    const real wUV = wU * iV;
    const real wVV = wV * iV;
    
    const MatrixBlock Id(0, 1);

    // directly setting the matrices
    const MatrixBlock duFu = MatrixBlock::offsetOuterProduct(-wUU, u, wUU);
    const MatrixBlock duFv = ( Id - MatrixBlock::outerProduct(v,u) ) * wUV;
    const MatrixBlock dvFv = MatrixBlock::offsetOuterProduct(-wVV, v, wVV);
    
    // indices:
    const index_t iiA = pt1.matIndex1();
    const index_t iiB = pt1.matIndex2();
    const index_t iiC = pt2.matIndex1();
    const index_t iiD = pt2.matIndex2();

    add_block_diag(iiA, duFu);
    sub_block(iiB, iiA, duFu);
    add_block_diag(iiB, duFu);
    if ( iiC > iiA )
    {
        add_block(iiC, iiA, duFv);
        sub_block(iiD, iiA, duFv);
        sub_block(iiC, iiB, duFv);
        add_block(iiD, iiB, duFv);
    }
    else
    {
        const MatrixBlock dvFu = duFv.transposed();
        add_block(iiA, iiC, dvFu);
        sub_block(iiA, iiD, dvFu);
        sub_block(iiB, iiC, dvFu);
        add_block(iiB, iiD, dvFu);
    }
    add_block_diag(iiC, dvFv);
    sub_block(iiD, iiC, dvFv);
    add_block_diag(iiD, dvFv);
    
    // remaining part of the force (should be small near equilibrium)
    Vector Fu0 = Fu - duFu.vecmul(AB) - duFv.trans_vecmul(CD);
    Vector Fv0 = Fv - duFv.vecmul(AB) - dvFv.vecmul(CD);

    // add constant terms into vBAS[]:
    sub_base(iiA, Fu0); // F(A) = -Fu
    add_base(iiB, Fu0); // F(B) = +Fu
    sub_base(iiC, Fv0); // F(C) = -Fv
    add_base(iiD, Fv0); // F(D) = +Fv
}


#endif



MatrixBlock Meca::torqueMatrix(real weight, Torque const& axi, Vector2 const& ang)
{
#if ( DIM >= 3 )
    return (-weight) * MatrixBlock::rotationAroundAxis(axi, ang.XX, ang.YY);
#elif ( DIM == 2 )
    return (-weight) * MatrixBlock(ang.XX, axi*ang.YY, -axi*ang.YY, ang.XX);
#else
    return MatrixBlock(0, -weight);  //should not be used!
#endif
}


/**
 Add Torque between 3 points to align AB with BC:
 
 F = k * ( -A + 2*B - C )
 
 F_A = F
 F_B = -2 * F
 F_C = F
 
 FJN, 11.05.2021
 */
void Meca::addTorque3(Mecapoint const& ptA,
                      Mecapoint const& ptB,
                      Mecapoint const& ptC,
                      const real scale, const real weight)
{
    assert_true( weight >= 0 );
    const MatrixBlock W(0, weight);

    const index_t iiA = ptA.matIndex0();
    const index_t iiB = ptB.matIndex0();
    const index_t iiC = ptC.matIndex0();
    
    sub_block_diag(iiA, W);
    if ( iiA < iiB )
        add_block(iiB, iiA, 2, W);
    else
        add_block(iiA, iiB, 2, W);
    if ( iiA < iiC )
        sub_block(iiC, iiA, W);
    else
        sub_block(iiA, iiC, W);
    
    add_block_diag(iiB, -4, W);
    
    if ( iiB < iiC )
        add_block(iiC, iiB, 2, W);
    else
        add_block(iiB, iiC, 2, W);
    sub_block_diag(iiC, W);
}


/**
 Add Torque between 3 points.
 This version does not impose any particular distance between the points,
 and just move them to enforce the angle described by ABC
 */
void Meca::addTorque3(Mecapoint const& ptA,
                      Mecapoint const& ptB,
                      Mecapoint const& ptC,
                      const MatrixBlock & R, //already multiplied by -weight
                      const real weight)
{
    assert_true( weight >= 0 );
    const MatrixBlock W(0, -weight);
    const MatrixBlock T = R.transposed();

    // indices:
    const index_t iiA = ptA.matIndex0();
    const index_t iiB = ptB.matIndex0();
    const index_t iiC = ptC.matIndex0();

    /*
    Vector CD = position_diff(ptC, ptB) + R * AB;
    Vector fA = T * CD;
    Vector fB = ( W + T ) * ( -CD );
    Vector fC = weight * CD;
     */
    
    add_block_diag(iiA, W);
    if ( iiA < iiB )
        sub_block(iiB, iiA, W+R);
    else
        sub_block(iiA, iiB, W+T);
    if ( iiA < iiC )
        add_block(iiC, iiA, R);
    else
        add_block(iiA, iiC, T);
    
#if ( DIM == 2 )
    // small optimization in 2D, as the term (R+W)+(T+W) is diagonal
    real dd = R.trace() - 2 * weight;
    add_block_diag(iiB, MatrixBlock(0, dd));
#else
    add_block_diag(iiB, (R+W)+(T+W));
#endif
    
    if ( iiB < iiC )
        sub_block(iiC, iiB, W+R);
    else
        sub_block(iiB, iiC, W+T);
    add_block_diag(iiC, W);
}


/**
 Add Torque between 3 points, projecting the force in the plane of rotation
 This version does not impose any particular distance between the points,
 and just move them to enforce the angle described by ABC
 */
void Meca::addTorque3Plane(Mecapoint const& ptA,
                           Mecapoint const& ptB,
                           Mecapoint const& ptC,
                           const Torque & axi, Vector2 const& ang, const real weight)
{
    assert_true( weight >= 0 );
    assert_small( ang.normSqr() - 1.0 );

#if ( DIM >= 3 )
    /*
    const Vector3 AB = position_diff(ptA, ptB);
    const Vector3 BC = position_diff(ptB, ptC);
    Vector3 axi = normalize(cross(AB, BC));
    if ( axi != axi )
        return;
     */
    const MatrixBlock X = MatrixBlock::outerProduct(axi);
    const MatrixBlock R = MatrixBlock::rotationAroundAxis(axi, ang.XX, ang.YY) - X;
    const MatrixBlock P = MatrixBlock(0,1) - X; // symmetric
#elif ( DIM == 2 )
    const MatrixBlock R = MatrixBlock(ang.XX, ang.YY,-ang.YY, ang.XX);
    const MatrixBlock P(0,1);
#else
    const MatrixBlock R(0,1);  //should not be used!
    const MatrixBlock P(0,1);
#endif
    
    // indices:
    const index_t iiA = ptA.matIndex0();
    const index_t iiB = ptB.matIndex0();
    const index_t iiC = ptC.matIndex0();

    const MatrixBlock wP = -weight * P;
    const MatrixBlock wR = -weight * R;
    const MatrixBlock wT = wR.transposed();

    add_block_diag(iiA, wP);
    if ( iiA < iiB )
        sub_block(iiB, iiA, wR+wP);
    else
        sub_block(iiA, iiB, wT+wP);
    if ( iiA < iiC )
        add_block(iiC, iiA, wR);
    else
        add_block(iiA, iiC, wT);

    add_block_diag(iiB, 2*wP+(wR+wT));
    if ( iiB < iiC )
        sub_block(iiC, iiB, wP+wR);
    else
        sub_block(iiB, iiC, wP+wT);

    add_block_diag(iiC, wP);
}



/** This is variation 3, 20.08.2019
 It combines addTorque() without length with a LongLink(ptA, ptB);
 */
void Meca::addTorque3Long(Mecapoint const& ptA,
                          Mecapoint const& ptB,
                          Mecapoint const& ptC,
                          const MatrixBlock & R,
                          const real weight,
                          const real len, const real weightL)
{
    assert_true( weight >= 0 );
    assert_true( weightL >= 0 );
    const MatrixBlock W(0, -weight);
    const Vector AB = position_diff(ptA, ptB);
    const MatrixBlock T = R.transposed();
    
    // indices:
    const index_t iiA = ptA.matIndex0();
    const index_t iiB = ptB.matIndex0();
    const index_t iiC = ptC.matIndex0();
    
    // this is a LongLink(A, B):
    MatrixBlock wL(0,0);
    const real ab2 = AB.normSqr();
    if ( ab2 > REAL_EPSILON )
    {
        real ab = std::sqrt(ab2);
        const real wla = weightL * len / ab;
        sub_base(iiA, AB, wla);
        add_base(iiB, AB, wla);

        // regularize interaction to avoid creating negative eigen values
        if ( ab < len )
            wL = MatrixBlock::outerProduct(AB, -weightL/ab2);
        else
            wL = MatrixBlock::offsetOuterProduct(wla-weightL, AB, -wla/ab2);
    }
    
    add_block_diag(iiA, W+wL);
    if ( iiA < iiB )
        sub_block(iiB, iiA, W+R+wL);
    else
        sub_block(iiA, iiB, W+T+wL);
    if ( iiA < iiC )
        add_block(iiC, iiA, R);
    else
        add_block(iiA, iiC, T);

#if ( DIM == 2 )
    // in 2D, the term (R+W)+(T+W) is diagonal
    real dd = R.trace() - 2 * weight;
    add_block_diag(iiB, wL+MatrixBlock(0,dd));
#else
    add_block_diag(iiB, (R+W)+(T+W)+wL);
#endif
    if ( iiB < iiC )
        sub_block(iiC, iiB, W+R);
    else
        sub_block(iiB, iiC, W+T);

    add_block_diag(iiC, W);
}


/**
 Add Torque between 4 points to align AB with CD, assuming |AB| = |CD|:
     
     F = weight * [ ( B - A ) - ( D - C ) ]
     
     F_A =  F
     F_B = -F
     F_C = -F
     F_D =  F
 
 with ptA = { A, B } and ptC = { C, D }
 FJN, 11.05.2021
 */
#if USE_ISO_MATRIX
void Meca::addTorque4(index_t iiA, index_t iiB, index_t iiC, index_t iiD, const real weight)
{
    assert_true( weight >= 0 );
    sub_iso_diag(iiA, weight);
    sub_iso_diag(iiB, weight);
    sub_iso_diag(iiC, weight);
    sub_iso_diag(iiD, weight);
    
    add_iso(iiB, iiA, weight);
    if ( iiC < iiD )
        add_iso(iiD, iiC, weight);
    else
        add_iso(iiC, iiD, weight);
    
    if ( iiA < iiC )
    {
        add_iso(iiC, iiA, weight);
        sub_iso(iiD, iiA, weight);
        sub_iso(iiC, iiB, weight);
        add_iso(iiD, iiB, weight);
    }
    else
    {
        add_iso(iiA, iiC, weight);
        sub_iso(iiB, iiC, weight);
        sub_iso(iiA, iiD, weight);
        add_iso(iiB, iiD, weight);
    }
}
#else
void Meca::addTorque4(index_t iiA, index_t iiB, index_t iiC, index_t iiD, const real weight)
{
    assert_true( weight >= 0 );
    
    const MatrixBlock W(0, weight);
    
    sub_block_diag(iiA, W);
    add_block(iiB, iiA, W);
    sub_block_diag(iiB, W);
    sub_block_diag(iiC, W);
    if ( iiC < iiD )
        add_block(iiD, iiC, W);
    else
        add_block(iiC, iiD, W);
    sub_block_diag(iiD, W);
    
    if ( iiA < iiC )
    {
        add_block(iiC, iiA, W);
        sub_block(iiD, iiA, W);
        sub_block(iiC, iiB, W);
        add_block(iiD, iiB, W);
    }
    else
    {
        add_block(iiA, iiC, W);
        sub_block(iiA, iiD, W);
        sub_block(iiB, iiC, W);
        add_block(iiB, iiD, W);
    }
}
#endif


void Meca::addTorque4(Mecapoint const& ptA, Mecapoint const& ptC, const real weight)
{
    const index_t iiA = ptA.matIndex0();
    const index_t iiC = ptC.matIndex0();
    addTorque4(iiA, iiA+1, iiC, iiC+1, weight);
}

//------------------------------------------------------------------------------
#pragma mark - Isotropic links between Mecables
//------------------------------------------------------------------------------

/**
 Link `ptA` (A) and `ptB` (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 In practice, Meca::addLink() will update the matrix,
 adding `weight` at the indices corresponding to `A` and `B`.
 
 Note: with modulo, the position of the fibers may be shifted in space,
 and a correction is necessary to make the force calculation correct:
 
     force_A = weight * ( B - A - offset )
     force_B = weight * ( A - B + offset )

 Here 'offset' is a multiple of the space periodicity, corresponding to B-A:
 offset = modulo->offset( A - B )

 In practice, Meca::addLink() will update the vector vBAS[]:
 
     vBAS[A] += weight * offset;
     vBAS[B] -= weight * offset;
 
 In principle, what goes to vBAS[] with modulo can be derived
 simply by multiplying the matrix block by 'offset'.
 */
void Meca::addLink(Mecapoint const& ptA,
                   Mecapoint const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    const index_t ii0 = ptA.matIndex0();
    const index_t ii1 = ptB.matIndex0();

    if ( ii0 == ii1 )
        return;

    sub_iso_diag(ii0, weight);
    sub_iso_diag(ii1, weight);
    add_iso(std::max(ii0, ii1), std::min(ii0, ii1), weight);

    if ( modulo )
    {
#if USE_GLOBAL_POSITION
        Vector off = modulo_offset(ptA, ptB);
#else
        Vector off = interpolate1(vPTS, ii1) - interpolate1(vPTS, ii0);
        off = modulo->offset(off);
#endif
        if ( off.is_not_zero() )
        {
            sub_base(ii0, off, weight);
            add_base(ii1, off, weight);
        }
    }
    DRAW_LINK(ptA, ptB.pos());
}


/**
 Link interpolated vertex (A) and vertex (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 */
void Meca::addLink(Interpolation const& ptA,
                   Mecapoint const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();
    
    if ( any_equal(ii2, ii0, ii1) )
        return;
    
    //coefficients to form B-A:
    const real cc0 = -ptA.coef0();
    const real cc1 = -ptA.coef1();
    //const real cc2 = 1.0;
    
    const real ww0 = weight * cc0;
    const real ww1 = weight * cc1;
    const real ww2 = weight;
    
    sub_iso_diag(ii0, ww0 * cc0);
    sub_iso_diag(ii1, ww1 * cc1);
    add_iso_diag(ii2, ww2);         // since cc2 == 1

    sub_iso(ii1, ii0, ww1 * cc0);
    if ( ii0 < ii2 )
    {
        add_iso(ii2, ii0, ww0);         // since cc2 == 1
        add_iso(ii2, ii1, ww1);         // since cc2 == 1
    }
    else
    {
        add_iso(ii0, ii2, ww0);         // since cc2 == 1
        add_iso(ii1, ii2, ww1);         // since cc2 == 1
    }
    if ( modulo )
    {
#if USE_GLOBAL_POSITION
        Vector off = modulo_offset(ptA, ptB);
#else
        Vector off = interpolate2(vPTS, ii0, cc0, cc1) + interpolate1(vPTS, ii2);
        off = modulo->offset(off);
#endif
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
        }
    }
    DRAW_LINK(ptA, ptB.pos());
}


/**
 Link vertex (A) and interpolated vertex (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )

 */
void Meca::addLink(Mecapoint const& ptA,
                   Interpolation const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex0();
    const index_t ii1 = ptB.matIndex1();
    const index_t ii2 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    //const real cc0 = -1.0;
    const real cc1 = ptB.coef0();
    const real cc2 = ptB.coef1();
    
    const real ww0 = -weight;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    
    add_iso_diag(ii0, ww0);  // since cc0 == -1.0
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);

    if ( ii0 < ii1 )
    {
        add_iso(ii1, ii0, ww1);  // since cc0 == -1.0
        add_iso(ii2, ii0, ww2);  // since cc0 == -1.0
    }
    else
    {
        add_iso(ii0, ii1, ww1);  // since cc0 == -1.0
        add_iso(ii0, ii2, ww2);  // since cc0 == -1.0
    }
    sub_iso(ii2, ii1, ww2 * cc1);
  
    if ( modulo )
    {
#if USE_GLOBAL_POSITION
        Vector off = modulo_offset(ptA, ptB);
#else
        Vector off = interpolate2(vPTS, ii1, cc1, cc2) - interpolate1(vPTS, ii0);
        off = modulo->offset(off);
#endif
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
        }
    }
    DRAW_LINK(ptB, ptA.pos());
}


/**
 Link `ptA` (A) and `ptB` (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )

 */
void Meca::addLink(Interpolation const& ptA,
                   Interpolation const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    const real cc0 = -ptA.coef0();
    const real cc1 = -ptA.coef1();
    const real cc2 = ptB.coef0();
    const real cc3 = ptB.coef1();
    
    assert_small(cc0+cc1+cc2+cc3);

    const real ww0 = weight * cc0;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    const real ww3 = weight * cc3;
    
    sub_iso_diag(ii0, ww0 * cc0);
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);

    sub_iso(ii1, ii0, ww1 * cc0);
    if ( ii0 < ii2 )
    {
        sub_iso(ii2, ii0, ww2 * cc0);
        sub_iso(ii3, ii0, ww3 * cc0);
        sub_iso(ii2, ii1, ww2 * cc1);
        sub_iso(ii3, ii1, ww3 * cc1);
    }
    else
    {
        sub_iso(ii0, ii2, ww2 * cc0);
        sub_iso(ii1, ii2, ww2 * cc1);
        sub_iso(ii0, ii3, ww3 * cc0);
        sub_iso(ii1, ii3, ww3 * cc1);
    }
    sub_iso(ii3, ii2, ww3 * cc2);

    if ( modulo )
    {
#if USE_GLOBAL_POSITION
        Vector off = modulo_offset(ptA, ptB);
#else
        Vector off = interpolate2(vPTS, ii0, cc0, cc1) + interpolate2(vPTS, ii2, cc2, cc3);
        off = modulo->offset(off);
#endif
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
        }
    }
    DRAW_LINK(ptA, ptB.pos());
}

//------------------------------------------------------------------------------
#pragma mark - Isotropic links between Mecables (higher order interpolation)
//------------------------------------------------------------------------------

/**
 Link interpolation (A) and vertex at index 'pts' (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in the vertex of a Mecable, at index 'pts'.
 */
void Meca::addLink1(Interpolation const& pti,
                    const index_t ii0,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii1 = pti.matIndex1();
    const index_t ii2 = pti.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    //const real cc0 = -1.0;
    const real cc1 = pti.coef0();
    const real cc2 = pti.coef1();
    
    const real ww0 = -weight;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    
    add_iso_diag(ii0, ww0); // since cc0 == -1.0
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    
    if ( ii0 < ii1 )
    {
        add_iso(ii1, ii0, ww1);
        add_iso(ii2, ii0, ww2);
    }
    else
    {
        add_iso(ii0, ii1, ww1);
        add_iso(ii0, ii2, ww2);
    }
    sub_iso(ii2, ii1, ww2 * cc1);
    if ( modulo )
    {
        Vector off = interpolate2(vPTS, ii1, cc1, cc2) - interpolate1(vPTS, ii0);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
        }
    }
}


/**
 Link vertex (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is a Mecapoint (ptA).
 Point B in interpolated over 2 vertices, specified by index 'pts',
 using the coefficients `cc1`, `cc2`.
 */
void Meca::addLink2(Mecapoint const& ptA, const index_t ii1,
                    const real cc1, const real cc2,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex0();
    const index_t ii2 = ii1+1;
    
    if ( any_equal(ii0, ii1, ii2) )
        return;

    //const real cc0 = -1.0;
    
    assert_small(cc1+cc2-1.0);

    const real ww0 = -weight;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    
    add_iso_diag(ii0, ww0); // since cc0 == -1.0
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);

    if ( ii0 < ii1 )
    {
        add_iso(ii1, ii0, ww1); // since cc0 == -1.0
        add_iso(ii2, ii0, ww2); // since cc0 == -1.0
    }
    else
    {
        add_iso(ii0, ii1, ww1); // since cc0 == -1.0
        add_iso(ii0, ii2, ww2); // since cc0 == -1.0
    }
    sub_iso(ii2, ii1, ww2 * cc1);
        
    if ( modulo )
    {
        Vector off = interpolate2(vPTS, ii1, cc1, cc2) - interpolate1(vPTS, ii0);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
        }
    }
}

/**
 Link `pti` (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in interpolated over 2 vertices, specified by index 'pts',
 using the coefficients `cc2`, `cc3`.
 */
void Meca::addLink2(Interpolation const& pti, const index_t ii2,
                    const real cc2, const real cc3,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    const index_t ii3 = ii2+1;
    
    const real cc0 = -pti.coef0();
    const real cc1 = -pti.coef1();
    
    assert_small(cc2+cc3-1.0);

    const real ww0 = weight * cc0;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    const real ww3 = weight * cc3;
    
    sub_iso_diag(ii0, ww0 * cc0);
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);

    sub_iso(ii1, ii0, ww1 * cc0);
    if ( ii0 < ii2 )
    {
        sub_iso(ii2, ii0, ww2 * cc0);
        sub_iso(ii3, ii0, ww3 * cc0);
        sub_iso(ii2, ii1, ww2 * cc1);
        sub_iso(ii3, ii1, ww3 * cc1);
    }
    else
    {
        sub_iso(ii0, ii2, ww2 * cc0);
        sub_iso(ii1, ii2, ww2 * cc1);
        sub_iso(ii0, ii3, ww3 * cc0);
        sub_iso(ii1, ii3, ww3 * cc1);
    }
    sub_iso(ii3, ii2, ww3 * cc2);
    
    if ( modulo )
    {
        Vector off = interpolate2(vPTS, ii0, cc0, cc1) + interpolate2(vPTS, ii2, cc2, cc3);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
        }
    }
}


/**
 Link vertex (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is a Mecapoint.
 Point B in interpolated over 3 vertices, specified by index in 'pts',
 using the coefficients `cc1`, `cc2`, `cc3`.
 */
void Meca::addLink3(Mecapoint const& ptA, const index_t ii1,
                    const real cc1, const real cc2, const real cc3,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex0();
    const index_t ii2 = ii1+1;
    const index_t ii3 = ii1+2;
    
    //const real cc0 = -1.0;
    
    assert_small(cc1+cc2+cc3-1.0);

    const real ww0 = -weight; // since cc0 = -1
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    const real ww3 = weight * cc3;
    
    add_iso_diag(ii0, ww0);  // since cc0 = -1
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);

    if ( ii0 < ii1 )
    {
        add_iso(ii1, ii0, ww1);
        add_iso(ii2, ii0, ww2);
        add_iso(ii3, ii0, ww3);
    }
    else
    {
        add_iso(ii0, ii1, ww1);
        add_iso(ii0, ii2, ww2);
        add_iso(ii0, ii3, ww3);
    }
    sub_iso(ii2, ii1, ww2 * cc1);
    sub_iso(ii3, ii1, ww3 * cc1);
    sub_iso(ii3, ii2, ww3 * cc2);
        
    if ( modulo )
    {
        Vector off = interpolate3(vPTS, ii1, cc1, cc2, cc3) - interpolate1(vPTS, ii0);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
        }
    }
}


/**
 Link `pti` (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in interpolated over 4 vertices, specified by index in 'pts',
 using the coefficients `cc2`, `cc3`, `cc4`.
*/
void Meca::addLink3(Interpolation const& pti, const index_t ii2,
                    const real cc2, const real cc3, const real cc4,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    const index_t ii3 = ii2+1;
    const index_t ii4 = ii2+2;
    
    const real cc0 = -pti.coef0();
    const real cc1 = -pti.coef1();
    
    assert_small(cc2+cc3+cc4-1.0);

    const real ww0 = weight * cc0;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    const real ww3 = weight * cc3;
    const real ww4 = weight * cc4;
    
    sub_iso_diag(ii0, ww0 * cc0);
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);
    sub_iso_diag(ii4, ww4 * cc4);

    sub_iso(ii1, ii0, ww1 * cc0);
    if ( ii0 < ii2 )
    {
        sub_iso(ii2, ii0, ww2 * cc0);
        sub_iso(ii3, ii0, ww3 * cc0);
        sub_iso(ii4, ii0, ww4 * cc0);
        sub_iso(ii2, ii1, ww2 * cc1);
        sub_iso(ii3, ii1, ww3 * cc1);
        sub_iso(ii4, ii1, ww4 * cc1);
    }
    else
    {
        sub_iso(ii0, ii2, ww2 * cc0);
        sub_iso(ii1, ii2, ww2 * cc1);
        sub_iso(ii0, ii3, ww3 * cc0);
        sub_iso(ii1, ii3, ww3 * cc1);
        sub_iso(ii0, ii4, ww4 * cc0);
        sub_iso(ii1, ii4, ww4 * cc1);
    }
    sub_iso(ii3, ii2, ww3 * cc2);
    sub_iso(ii4, ii2, ww4 * cc2);
    sub_iso(ii4, ii3, ww4 * cc3);
        
    if ( modulo )
    {
        Vector off = interpolate2(vPTS, ii0, cc0, cc1) + interpolate3(vPTS, ii2, cc2, cc3, cc4);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
            add_base(ii4, off, ww4);
        }
    }
}


/**
 Link vertex (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is a Mecapoint.
 Point B in interpolated over 4 vertices, specified by index in 'pts',
 using the coefficients `cc1`, `cc2`, `cc3`, `cc4`.
 */
void Meca::addLink4(Mecapoint const& ptA, const index_t ii1,
                    const real cc1, const real cc2, const real cc3, const real cc4,
                    const real weight)
{
    assert_true( weight >= 0 );

    // indices:
    const index_t ii0 = ptA.matIndex0();
    const index_t ii2 = ii1+1;
    const index_t ii3 = ii1+2;
    const index_t ii4 = ii1+3;

    //const real cc0 = -1.0;
    assert_small(cc1+cc2+cc3+cc4-1.0);

    const real ww0 = -weight;  // since cc0 = -1
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    const real ww3 = weight * cc3;
    const real ww4 = weight * cc4;
    
    add_iso_diag(ii0, ww0);  // since cc0 = -1
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);
    sub_iso_diag(ii4, ww4 * cc4);

    if ( ii0 < ii1 )
    {
        add_iso(ii1, ii0, ww1);
        add_iso(ii2, ii0, ww2);
        add_iso(ii3, ii0, ww3);
        add_iso(ii4, ii0, ww4);
    }
    else
    {
        add_iso(ii0, ii1, ww1);
        add_iso(ii0, ii2, ww2);
        add_iso(ii0, ii3, ww3);
        add_iso(ii0, ii4, ww4);
    }
    
    sub_iso(ii2, ii1, ww2 * cc1);
    sub_iso(ii3, ii1, ww3 * cc1);
    sub_iso(ii4, ii1, ww4 * cc1);
    sub_iso(ii3, ii2, ww3 * cc2);
    sub_iso(ii4, ii2, ww4 * cc2);
    sub_iso(ii4, ii3, ww4 * cc3);
        
    if ( modulo )
    {
        Vector off = interpolate4(vPTS, ii1, cc1, cc2, cc3, cc4) - interpolate1(vPTS, ii0);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
            add_base(ii4, off, ww4);
        }
    }
}


/**
 Link `pti` (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in interpolated over 4 vertices, specified by index in 'pts',
 using the coefficients `cc1`, `cc2`, `cc3`, `cc4`, `cc5`.
*/
void Meca::addLink4(Interpolation const& pti, const index_t ii2,
                    const real cc2, const real cc3, const real cc4, const real cc5,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    const index_t ii3 = ii2+1;
    const index_t ii4 = ii2+2;
    const index_t ii5 = ii2+3;

    const real cc0 = -pti.coef0();
    const real cc1 = -pti.coef1();
    
    assert_small(cc2+cc3+cc4+cc5-1.0);

    const real ww0 = weight * cc0;
    const real ww1 = weight * cc1;
    const real ww2 = weight * cc2;
    const real ww3 = weight * cc3;
    const real ww4 = weight * cc4;
    const real ww5 = weight * cc5;
    
    sub_iso_diag(ii0, ww0 * cc0);
    sub_iso_diag(ii1, ww1 * cc1);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);
    sub_iso_diag(ii4, ww4 * cc4);
    sub_iso_diag(ii5, ww5 * cc5);

    sub_iso(ii1, ii0, ww1 * cc0);
    if ( ii0 < ii2 )
    {
        sub_iso(ii2, ii0, ww2 * cc0);
        sub_iso(ii3, ii0, ww3 * cc0);
        sub_iso(ii4, ii0, ww4 * cc0);
        sub_iso(ii5, ii0, ww5 * cc0);
        
        sub_iso(ii2, ii1, ww2 * cc1);
        sub_iso(ii3, ii1, ww3 * cc1);
        sub_iso(ii4, ii1, ww4 * cc1);
        sub_iso(ii5, ii1, ww5 * cc1);
    }
    else
    {
        sub_iso(ii0, ii2, ww2 * cc0);
        sub_iso(ii1, ii2, ww2 * cc1);
        sub_iso(ii0, ii3, ww3 * cc0);
        sub_iso(ii1, ii3, ww3 * cc1);
        sub_iso(ii0, ii4, ww4 * cc0);
        sub_iso(ii1, ii4, ww4 * cc1);
        sub_iso(ii0, ii5, ww5 * cc0);
        sub_iso(ii1, ii5, ww5 * cc1);
    }
    sub_iso(ii3, ii2, ww3 * cc2);
    sub_iso(ii4, ii2, ww4 * cc2);
    sub_iso(ii5, ii2, ww5 * cc2);
    sub_iso(ii4, ii3, ww4 * cc3);
    sub_iso(ii5, ii3, ww5 * cc3);
    sub_iso(ii5, ii4, ww5 * cc4);

    if ( modulo )
    {
        Vector off = interpolate2(vPTS, ii0, cc0, cc1) + interpolate4(vPTS, ii2, cc2, cc3, cc4, cc5);
        off = modulo->offset(off);
        if ( off.is_not_zero() )
        {
            add_base(ii0, off, ww0);
            add_base(ii1, off, ww1);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
            add_base(ii4, off, ww4);
            add_base(ii5, off, ww5);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Links with resting length
//------------------------------------------------------------------------------

/**
 Link `ptA` (A) and `ptB` (B),
 The force is affine with non-zero resting length:
 
     force_A = weight * ( B - A ) * ( length / |AB| - 1 )
     force_B = weight * ( A - B ) * ( length / |AB| - 1 )
 
 This streamlined version of addLongLink() is used for Steric interaction, with:
 - axi = position(ptB) - position(ptA)
 - ab2 = normSqr(axi)
 - periodic boundary conditions have been applied to `axi` if necessary
 - 'ab2 < len * len', thus leading to a repulsive force only.
 
 */

void Meca::addLongLink1(Mecapoint const& ptA,
                        Mecapoint const& ptB,
                        Vector const& axi,
                        const real ab2,
                        const real len,
                        const real weight)
{
    assert_true( ab2 > REAL_EPSILON );
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    const index_t aa = ptA.matIndex0();  // coef is +weight
    const index_t bb = ptB.matIndex0();  // coef is -weight
    assert_false( aa == bb );

    DRAW_LINK(ptA, axi, len);

    const real abn = std::sqrt(ab2);
    const real wab = -weight / ab2;
    MatrixBlock wT = MatrixBlock::outerProduct(axi, wab);
    
    Vector tmp = axi * (( wab * len ) * abn);
    add_base(aa, tmp);
    sub_base(bb, tmp);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(aa, off);
            sub_base(bb, off);
        }
    }

    add_block_diag(aa, wT);
    add_block_diag(bb, wT);
    sub_block(std::max(aa, bb), std::min(aa, bb), wT);
}


/**
Link `ptA` (A) and `ptB` (B),
The force is affine with non-zero resting length:

    force_A = weight * ( B - A ) * ( length / |AB| - 1 )
    force_B = weight * ( A - B ) * ( length / |AB| - 1 )

This streamlined version of addLongLink() is used for Steric interaction, with:
 - axi = position(ptB) - position(ptA)
 - ab2 = normSqr(axi)
 - periodic boundary conditions have been applied to `axi` if necessary
*/
void Meca::addLongLink2(Mecapoint const& ptA,
                        Mecapoint const& ptB,
                        Vector const& axi,
                        const real ab2,
                        const real len,
                        const real weight)
{
    assert_true( ab2 > REAL_EPSILON );
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    const index_t aa = ptA.matIndex0();  // coef is +weight
    const index_t bb = ptB.matIndex0();  // coef is -weight
    assert_false( aa == bb );

    DRAW_LINK(ptA, axi, len);
    
    const real iab = 1.0 / ab2;
    const real abn = std::sqrt(ab2);
    const real wla = weight * len * abn * iab; // weight * len / abn
    
    MatrixBlock wT;
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight*iab);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla*iab);
    
    Vector tmp = axi * wla;
    sub_base(aa, tmp);
    add_base(bb, tmp);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(aa, off);
            sub_base(bb, off);
        }
    }

    add_block_diag(aa, wT);
    add_block_diag(bb, wT);
    sub_block(std::max(aa, bb), std::min(aa, bb), wT);
}


/**
 Link `ptA` (A) and `ptB` (B),
 The force is affine with non-zero resting length: 
 
     force_A = weight * ( B - A ) * ( length / |AB| - 1 )
     force_B = weight * ( A - B ) * ( length / |AB| - 1 )
 
 */

void Meca::addLongLink(Mecapoint const& ptA,
                       Mecapoint const& ptB,
                       const real len,
                       const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    const index_t aa = ptA.matIndex0();  // coef is +weight
    const index_t bb = ptB.matIndex0();  // coef is -weight

    if ( aa == bb )
        return;
    
    Vector off, axi = position_diff(ptA, ptB);

    if ( modulo )
        off = modulo->foldOffset(axi);

    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = std::sqrt(ab2);
    const real wla = weight * len / abn;
    
    DRAW_LINK(ptA, axi, len);
    
    MatrixBlock wT;
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target: real wla = weight * min_real(len/abn, 1); */
    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);
    
    axi *= wla;
    sub_base(aa, axi);
    add_base(bb, axi);
    
    if ( modulo && off.is_not_zero() )
    {
        off = wT * off;
        sub_base(aa, off);
        add_base(bb, off);
    }

    add_block_diag(aa, wT);
    add_block_diag(bb, wT);
    sub_block(std::max(aa, bb), std::min(aa, bb), wT);
}


/**
 Link vertex (A) and interpolation (B),
 The force is affine with non-zero resting length: 
 
     force_A = weight * ( B - A ) * ( len / |AB| - 1 )
     force_B = weight * ( A - B ) * ( len / |AB| - 1 )
 
 */

void Meca::addLongLink(Mecapoint const& ptA,
                       Interpolation const& ptB,
                       const real len,
                       const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    // indices:
    const index_t ii0 = ptB.matIndex1();
    const index_t ii1 = ptB.matIndex2();
    const index_t ii2 = ptA.matIndex0();

    if ( any_equal(ii2, ii0, ii1) )
        return;

    Vector off, axi = position_diff(ptA, ptB);

    if ( modulo )
        off = modulo->foldOffset(axi);

    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = std::sqrt(ab2);
    const real wla = weight * len / abn;

    // interpolation coefficients:
    const real cc0 = ptB.coef0();
    const real cc1 = ptB.coef1();

    add_base(ii0, axi, cc0 * wla);
    add_base(ii1, axi, cc1 * wla);
    sub_base(ii2, axi, wla);

    MatrixBlock wT;
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);
    
    if ( modulo && off.is_not_zero() )
    {
        off = -( wT * off );
        add_base(ii0, off, cc0);
        add_base(ii1, off, cc1);
        add_base(ii2, off);
    }
    DRAW_LINK(ptB, axi, len);

    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block_diag(ii2, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    if ( ii1 < ii2 )
    {
        sub_block(ii2, ii0, cc0, wT);
        sub_block(ii2, ii1, cc1, wT);
    }
    else
    {
        sub_block(ii0, ii2, cc0, wT);
        sub_block(ii1, ii2, cc1, wT);
    }
}


/**
 Link `ptA` (A) and `ptB` (B),
 The force is affine with non-zero resting length: 
 
     force_A = weight * ( B - A ) * ( len / |AB| - 1 )
     force_B = weight * ( A - B ) * ( len / |AB| - 1 )

 */

void Meca::addLongLink(Interpolation const& ptA,
                       Interpolation const& ptB,
                       const real len,
                       const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();

    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    // coefficients to form B-A:
    const real cc0 = -ptA.coef0();
    const real cc1 = -ptA.coef1();
    const real cc2 =  ptB.coef0();
    const real cc3 =  ptB.coef1();

    Vector off, axi = position_diff(ptA, ptB);

    if ( modulo )
        off = modulo->foldOffset(axi);

    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = std::sqrt(ab2);
    const real wla = weight * len / abn;

    add_base(ii0, axi, cc0*wla);
    add_base(ii1, axi, cc1*wla);
    add_base(ii2, axi, cc2*wla);
    add_base(ii3, axi, cc3*wla);
    
    MatrixBlock wT;
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);
    
    if ( modulo && off.is_not_zero() )
    {
        off = -( wT * off );
        add_base(ii0, off, cc0);
        add_base(ii1, off, cc1);
        add_base(ii2, off, cc2);
        add_base(ii3, off, cc3);
    }
    DRAW_LINK(ptA, axi, len);

    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block_diag(ii2, cc2*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block(ii3, ii2, cc3*cc2, wT);
}


/**
Link `ptA` (A) and the point described by `inx` and `cc2`, `cc3`, `cc4`, `cc5` (B),
The force is affine with non-zero resting length:

    force_A = weight * ( B - A ) * ( len / |AB| - 1 )
    force_B = weight * ( A - B ) * ( len / |AB| - 1 )

*/
void Meca::addLongLink4(Interpolation const& ptA, const index_t ii2,
                        const real cc2, const real cc3, const real cc4, const real cc5,
                        const real len,
                        const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii3 = ii2+1;
    const index_t ii4 = ii2+2;
    const index_t ii5 = ii2+3;

    if ( any_equal(ii0, ii2, ii3, ii4, ii5) || any_equal(ii1, ii2, ii3, ii4, ii5) )
        return;
    
    // coefficients to form -A:
    const real cc0 = -ptA.coef0();
    const real cc1 = -ptA.coef1();
    
    assert_small(cc2+cc3+cc4+cc5-1.0);

    Vector axi = interpolate4(vPTS, ii2, cc2, cc3, cc4, cc5) + interpolate2(vPTS, ii0, cc0, cc1);
    Vector off;

    if ( modulo )
        off = modulo->foldOffset(axi);

    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = std::sqrt(ab2);
    const real wla = weight * len / abn;

    add_base(ii0, axi, cc0*wla);
    add_base(ii1, axi, cc1*wla);
    add_base(ii2, axi, cc2*wla);
    add_base(ii3, axi, cc3*wla);
    add_base(ii4, axi, cc4*wla);
    add_base(ii5, axi, cc5*wla);

    MatrixBlock wT;
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);
    
    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block_diag(ii2, cc2*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);
    add_block_diag(ii4, cc4*cc4, wT);
    add_block_diag(ii5, cc5*cc5, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii4, ii0, cc4*cc0, wT);
        add_block(ii5, ii0, cc5*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
        add_block(ii4, ii1, cc4*cc1, wT);
        add_block(ii5, ii1, cc5*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
        add_block(ii0, ii4, cc4*cc0, wT);
        add_block(ii1, ii4, cc4*cc1, wT);
        add_block(ii0, ii5, cc5*cc0, wT);
        add_block(ii1, ii5, cc5*cc1, wT);
    }
    add_block(ii3, ii2, cc3*cc2, wT);
    add_block(ii4, ii2, cc4*cc2, wT);
    add_block(ii5, ii2, cc5*cc2, wT);
    add_block(ii4, ii3, cc4*cc3, wT);
    add_block(ii5, ii3, cc5*cc3, wT);
    add_block(ii5, ii4, cc5*cc4, wT);

    if ( modulo && off.is_not_zero() )
    {
        off = -( wT * off );
        add_base(ii0, off, cc0);
        add_base(ii1, off, cc1);
        add_base(ii2, off, cc2);
        add_base(ii3, off, cc3);
        add_base(ii4, off, cc4);
        add_base(ii5, off, cc5);
    }
    DRAW_LINK(ptA, axi, len);
}


//------------------------------------------------------------------------------
#pragma mark - Off-axis links between Mecables
//------------------------------------------------------------------------------

/**
 Link `ptA` (A) and `ptB` (B),
 through an intermediate point S located on the side of the segment supporting A:
     S = A + len * N,
 where N is a unit vector that is orthogonal to the A-fiber in A.
 S is linearly related to the two vertices located on each side of A.
 The resulting force is linear of zero resting length, between B and S:
 
     force_S = weight * ( S - B )
     force_B = weight * ( B - S )
 
 */


#if ( DIM == 2 )

void Meca::addSideLink2D(Interpolation const& ptA,
                         Mecapoint const& ptB,
                         const real arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();

    if ( any_equal(ii2, ii0, ii1) )
        return;
    const real ee = arm * ptA.lenInv(), we = weight * ee;

    // interpolation coefficients and weights:
    const real cc0 = ptA.coef0(),  ww0 = weight * cc0;
    const real cc1 = ptA.coef1(),  ww1 = weight * cc1;
    
    sub_iso_diag(ii0, ww0 * cc0 + we * ee);
    sub_iso_diag(ii1, ww1 * cc1 + we * ee);
    sub_iso_diag(ii2, weight);
    
    real wd = ww0 * cc1 - we * ee;
    sub_block(ii1, ii0, Matrix22(wd, -we, we, wd));
    
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, Matrix22(ww0, -we,  we, ww0));
        add_block(ii2, ii1, Matrix22(ww1,  we, -we, ww1));
    }
    else
    {
        add_block(ii0, ii2, Matrix22(ww0,  we, -we, ww0));
        add_block(ii1, ii2, Matrix22(ww1, -we,  we, ww1));
    }
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            Matrix22 wAt(ww0,  we, -we, ww0);
            Matrix22 wBt(ww1, -we,  we, ww1);
            sub_base(ii0, wAt*off);
            sub_base(ii1, wBt*off);
            add_base(ii2, off, weight);
        }
    }
    DRAW_LINK(ptA, cross(arm, ptA.dir()), ptB.pos());
}

#endif


/**
 Link `B` to an interpolated point `S` on the side of `A`:
 
     S = A.pos() + cross( arm, A.dir() )
 
 Where A.dir() is the direction of the Fiber supporting `A`, in `A`.
 The vector `arm` should ideally be perpendicular to A.dir(), and in this case,
 `A` and `S` are separated by norm(arm). For best results, `arm` should also be
 perpendicular to the link axis `B-A` to separate B from A by |arm|.
 The forces are linear and sum up to zero:

     force(B) = weight * ( S - B )
     force(S) = -force(B)

 The `force_S` is redistributed on the vertices on each side of `A`,
 according to the interpolation coefficients, as a function of `arm`.
 
 This code is valid in any dimension and works in 2 and 3D
 */
void Meca::addSideLink3D(Interpolation const& ptA,
                         Mecapoint const& ptB,
                         Torque const& arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();
    
    if ( any_equal(ii2, ii0, ii1) )
        return;
    
    // interpolation coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    const Torque leg = arm * ptA.lenInv();

    MatrixBlock aR = MatrixBlock::vectorProduct(cc0, -leg);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1,  leg);
    
#if 0
    std::cerr.precision(3);
    // check that image D calculated from (AB) is near its target C
    Vector D = aR * ptA.pos1() + bR * ptA.pos2();
    std::cerr << "D: " << std::setw(9) << D << '\n';
    std::cerr << "C: " << std::setw(9) << ptB.pos() << '\n';
#endif
    
    MatrixBlock wAt = aR.transposed(-weight);
    MatrixBlock wBt = bR.transposed(-weight);
    
    //the following 3 terms are symmetric but not diagonal
    add_block_diag(ii0, wAt.mul(aR));
    add_block(ii1, ii0, wBt.mul(aR));
    add_block_diag(ii1, wBt.mul(bR));

    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, weight, aR);
        add_block(ii2, ii1, weight, bR);
    }
    else
    {
        sub_block(ii0, ii2, wAt);
        sub_block(ii1, ii2, wBt);
    }
    add_block_diag(ii2, MatrixBlock(0, -weight));
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            add_base(ii0, wAt*off);
            add_base(ii1, wBt*off);
            add_base(ii2, off, weight);
        }
    }
    DRAW_LINK(ptA, cross(arm, ptA.dir()), ptB.pos());
}


void Meca::addSideLink(Interpolation const& ptA,
                       Mecapoint const& ptB,
                       const real len,
                       const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real arm = std::copysign(len, cross(ptA.diff(), ptB.pos()-ptA.pos1()));
    addSideLink2D(ptA, ptB, arm, weight);

#else
    
    // set 'arm' perpendicular to Fiber and link:
    //Vector arm = cross(ptA.diff(), ptB.pos()-ptA.pos1());
    Vector arm = ptB.pos()-ptA.pos1();
    if ( modulo )
        modulo->fold(arm);
    arm = cross(ptA.diff(), arm);
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSideLink3D(ptA, ptB, arm*(len/n), weight);
    
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Off-axis links between Mecable (Interpolation)

#if ( DIM == 2 )

void Meca::addSideLink2D(Interpolation const& ptA,
                         Interpolation const& ptB,
                         const real arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    const real ee = -arm * ptA.lenInv();
    const real we = weight * ee;

    // coefficients to form B-A:
    const real cc0 = -ptA.coef0(),  ww0 = weight * cc0;
    const real cc1 = -ptA.coef1(),  ww1 = weight * cc1;
    const real cc2 =  ptB.coef0(),  ww2 = weight * cc2;
    const real cc3 =  ptB.coef1(),  ww3 = weight * cc3;

    sub_iso_diag(ii0, ww0 * cc0 + we * ee);
    sub_iso_diag(ii1, ww1 * cc1 + we * ee);
    sub_iso_diag(ii2, ww2 * cc2);
    sub_iso_diag(ii3, ww3 * cc3);
    sub_iso(ii3, ii2, ww3 * cc2);
    
#if ( 1 )
    // equivalent shortcut, since cc0+cc1 = -1
    real wd = ww0 * cc1 - we * ee;
    sub_block(ii1, ii0, Matrix22(wd, we, -we, wd));
#else
    Matrix22 A(cc0, -ee,  ee, cc0);
    Matrix22 B(cc1,  ee, -ee, cc1);
    sub_block(ii1, ii0, weight, B.trans_mul(A));
#endif
    
    if ( ii0 < ii2 )
    {
        Matrix22 wA(ww0, -we,  we, ww0);
        Matrix22 wB(ww1,  we, -we, ww1);
        sub_block(ii2, ii0, cc2, wA);
        sub_block(ii3, ii0, cc3, wA);
        sub_block(ii2, ii1, cc2, wB);
        sub_block(ii3, ii1, cc3, wB);
    }
    else
    {
        Matrix22 wAt(ww0,  we, -we, ww0);
        Matrix22 wBt(ww1, -we,  we, ww1);
        sub_block(ii0, ii2, cc2, wAt);
        sub_block(ii1, ii2, cc2, wBt);
        sub_block(ii0, ii3, cc3, wAt);
        sub_block(ii1, ii3, cc3, wBt);
    }
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            Matrix22 wAt(ww0,  we, -we, ww0);
            Matrix22 wBt(ww1, -we,  we, ww1);
            add_base(ii0, wAt*off);
            add_base(ii1, wBt*off);
            add_base(ii2, off, ww2);
            add_base(ii3, off, ww3);
        }
    }
    DRAW_LINK(ptA, cross(arm, ptA.dir()), ptB.pos());
}

#endif


/**
 Link `B` to an interpolated point `S` on the side of `A`:
 
     S = position(A) + cross( arm, A.dir() )
 
 Where A.dir() is the direction of the Fiber supporting `A`, in `A`.
 The vector `arm` should ideally be perpendicular to A.dir(), and in this case,
 `A` and `S` are separated by norm(arm). For best results, `arm` should also be
 perpendicular to the link axis `B-A` to separate B from A by |arm|.
 The forces are linear and sum up to zero:
 
     force(B) = weight * ( S - B )
     force(S) = -force(B)
 
 The `force_S` is redistributed on the vertices on each side of `A`,
 according to the interpolation coefficients, as a function of `arm`.
 
 This code is valid in any dimension and works in 2 and 3D
 */
void Meca::addSideLink3D(Interpolation const& ptA,
                         Interpolation const& ptB,
                         Torque const& arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();

    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    //coefficients to form A-B:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    const Torque leg = arm * ptA.lenInv();

    MatrixBlock aR = MatrixBlock::vectorProduct(cc0, -leg);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1,  leg);

    const real wcc2 = -weight * cc2;
    const real wcc3 = -weight * cc3;

    MatrixBlock wAt = aR.transposed(-weight);
    MatrixBlock wBt = bR.transposed(-weight);
    
    //the following 3 terms are symmetric but not diagonal
    add_block_diag(ii0, wAt.mul(aR));
    add_block(ii1, ii0, wBt.mul(aR));
    add_block_diag(ii1, wBt.mul(bR));
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, wcc2, aR);
        add_block(ii3, ii0, wcc3, aR);
        add_block(ii2, ii1, wcc2, bR);
        add_block(ii3, ii1, wcc3, bR);
    }
    else
    {
        add_block(ii0, ii2, cc2, wAt);
        add_block(ii1, ii2, cc2, wBt);
        add_block(ii0, ii3, cc3, wAt);
        add_block(ii1, ii3, cc3, wBt);
    }
    add_block_diag(ii2, MatrixBlock(0, wcc2*cc2));
    add_block(ii3, ii2, MatrixBlock(0, wcc3*cc2));
    add_block_diag(ii3, MatrixBlock(0, wcc3*cc3));
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            add_base(ii0, wAt*off);
            add_base(ii1, wBt*off);
            add_base(ii2, off, wcc2);
            add_base(ii3, off, wcc3);
        }
    }
    DRAW_LINK(ptA, cross(arm, ptA.dir()), ptB.pos());
}

 
 
/**
 Debug code to compare interactions in conditions where they should be identical
 you should run this and call 'opendiff x y'
*/
void Meca::testSideLink(Interpolation const& ptA,
                        Mecapoint const& ptB,
                        Torque const& arm,
                        const real weight)
{
    mFUL.reset();
    addSideLink3D(ptA, ptB, arm, weight);
    {
        std::ofstream o("x");
        o << "testSideLink " << ptA.matIndex1() << " " << ptB.matIndex0() << "\n";
        mFUL.printSparse(o, REAL_EPSILON);
    }
    
    real alpha = 0;
    index_t P = ptB.point();
    if ( P > 0 ) {
        alpha = 1;
        P -= 1;
    }
    
    mFUL.reset();
    addSideLink3D(ptA, Interpolation(ptB.mecable(), alpha, P, P+1), arm, weight);
    {
        std::ofstream o("y");
        o << "testSideLink " << ptA.matIndex1() << " " << ptB.matIndex0() << "\n";
        mFUL.printSparse(o, REAL_EPSILON);
    }
}



/**
 Link `ptA` (A) and `ptB` (B),
 Which is taken between B and a point S located on the side of A:
 
     S = A + len * N
 
 where N is a normalized vector orthogonal to the fiber in A.
 S is linearly related to the two vertices on the sides of A, P1 and P2
 In 3D, S is choosen in the plane of P1, P2 and B.
 The force is linear of zero resting length:
 
     force_S = weight * ( S - B )
     force_B = weight * ( B - S )
 
 */

void Meca::addSideLink(Interpolation const& ptA,
                       Interpolation const& ptB,
                       const real len,
                       const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideLink is meaningless in 1D");

#elif ( DIM == 2 )
    
    real arm = std::copysign(len, cross(ptA.diff(), ptB.pos()-ptA.pos1()));
    addSideLink2D(ptA, ptB, arm, weight);
    
#else

    // set 'arm' perpendicular to Fiber and link:
    //Vector arm = cross(ptA.diff(), ptB.pos()-ptA.pos1());
    Vector arm = ptB.pos()-ptA.pos1();
    if ( modulo )
        modulo->fold(arm);
    arm = cross(ptA.diff(), arm);
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSideLink3D(ptA, ptB, arm*(len/n), weight);
    
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Symmetric off-axis links and tilted links
//------------------------------------------------------------------------------

#if ( DIM == 2 )

/// this is old style code, addressing mFUL() directly, but it works!
void Meca::addSideSideLink2D(Interpolation const& ptA, const real armA,
                             Interpolation const& ptB, const real armB,
                             const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    index_t ia1 = ptA.matIndex1(), ia2 = ptA.matIndex2();
    index_t ib1 = ptB.matIndex1(), ib2 = ptB.matIndex2();
    
    if ( any_equal(ia1, ib1, ib2) || any_equal(ia2, ib1, ib2) )
        return;

    const real ca1 =  ptA.coef0(), ca2 =  ptA.coef1();
    const real cb1 = -ptB.coef0(), cb2 = -ptB.coef1();
    
    const real ee1 = armA * ptA.lenInv();
    const real ee2 = armB * ptB.lenInv();
    
    const real W = -weight;
    const real ca1w = ca1 * W, ca2w = ca2 * W;
    const real cb1w = cb1 * W, cb2w = cb2 * W;
   
    const real ee1w = ee1 * W, ee1ee1w = ee1 * ee1w;
    const real ee2w = ee2 * W, ee2ee2w = ee2 * ee2w;
    const real ee1ee2w = ee1 * ee2w;
    
    //isotropic terms
    add_iso_diag(ia1, ca1w * ca1 + ee1ee1w);
    add_iso_diag(ia2, ca2w * ca2 + ee1ee1w);
    add_iso_diag(ib1, cb1w * cb1 + ee2ee2w);
    add_iso_diag(ib2, cb2w * cb2 + ee2ee2w);

    add_iso(ia2, ia1, ca1w * ca2 - ee1ee1w);
    add_iso(ib2, ib1, cb1w * cb2 - ee2ee2w);
    
    if ( ia1 > ib1 )
    {
        add_iso(ia1, ib1, ca1w * cb1 - ee1ee2w);
        add_iso(ia1, ib2, ca1w * cb2 + ee1ee2w);
        add_iso(ia2, ib1, ca2w * cb1 + ee1ee2w);
        add_iso(ia2, ib2, ca2w * cb2 - ee1ee2w);
    }
    else
    {
        add_iso(ib1, ia1, ca1w * cb1 - ee1ee2w);
        add_iso(ib2, ia1, ca1w * cb2 + ee1ee2w);
        add_iso(ib1, ia2, ca2w * cb1 + ee1ee2w);
        add_iso(ib2, ia2, ca2w * cb2 - ee1ee2w);
    }
    // convert to full indices to directly address mFUL:
    ia1 *= DIM;
    ia2 *= DIM;
    ib1 *= DIM;
    ib2 *= DIM;
    
    mFUL(ia1  , ia2+1) -= ee1w;
    mFUL(ia1+1, ia2  ) += ee1w;
    
    mFUL(ib1  , ib2+1) -= ee2w;
    mFUL(ib1+1, ib2  ) += ee2w;
    
    const real ee1cb1w = ee1w * cb1;
    const real ee1cb2w = ee1w * cb2;
    const real ee2ca1w = ee2w * ca1;
    const real ee2ca2w = ee2w * ca2;
    
    mFUL(ia1, ib1+1) -=  ee2ca1w + ee1cb1w;
    mFUL(ia1, ib2+1) +=  ee2ca1w - ee1cb2w;
    
    mFUL(ia1+1, ib1) +=  ee2ca1w + ee1cb1w;
    mFUL(ia1+1, ib2) -=  ee2ca1w - ee1cb2w;
    
    mFUL(ia2, ib1+1) -=  ee2ca2w - ee1cb1w;
    mFUL(ia2, ib2+1) +=  ee2ca2w + ee1cb2w;
    
    mFUL(ia2+1, ib1) +=  ee2ca2w - ee1cb1w;
    mFUL(ia2+1, ib2) -=  ee2ca2w + ee1cb2w;
    
    DRAW_LINK(ptA, cross(armA, ptA.dir()), cross(armB, ptB.dir()), ptB.pos());
    
    if ( modulo )
        throw Exception("addSideSideLink2D is not usable with periodic boundary conditions");
}

#endif

// new style code, validated on 19.01.2020
void Meca::addSideSideLink(Interpolation const& ptA, Torque const& armA,
                           Interpolation const& ptB, Torque const& armB,
                           const real weight)
{
    assert_true( weight >= 0 );

    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
 
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;

    //coefficients to form A-B:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    
    const Torque legA = armA * ptA.lenInv();
    const Torque legB = armB * ptB.lenInv();

    MatrixBlock A = MatrixBlock::vectorProduct(cc0, -legA);
    MatrixBlock B = MatrixBlock::vectorProduct(cc1,  legA);
    MatrixBlock C = MatrixBlock::vectorProduct(cc2,  legB);
    MatrixBlock D = MatrixBlock::vectorProduct(cc3, -legB);

    MatrixBlock wAt = A.transposed(-weight);
    MatrixBlock wBt = B.transposed(-weight);
    MatrixBlock wCt = C.transposed(-weight);
    MatrixBlock wDt = D.transposed(-weight);
        
    /*
     We use block operations to set the matrix lower blocks:
            ii0  ii1  ii2  ii3
           ---------------------
    ii0   | A'A                |
    ii1   | B'A  B'B           |
    ii2   | C'A  C'B  C'C      |
    ii3   | D'A  D'B  D'C  D'D |
     */

    add_block_diag(ii0, wAt.mul(A));
    add_block_diag(ii1, wBt.mul(B));
    add_block_diag(ii2, wCt.mul(C));
    add_block_diag(ii3, wDt.mul(D));
 
    add_block(ii1, ii0, wBt.mul(A));
    add_block(ii3, ii2, wDt.mul(C));

    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, wCt.mul(A));
        add_block(ii3, ii0, wDt.mul(A));
        add_block(ii2, ii1, wCt.mul(B));
        add_block(ii3, ii1, wDt.mul(B));
    }
    else
    {
        add_block(ii0, ii2, wAt.mul(C));
        add_block(ii1, ii2, wBt.mul(C));
        add_block(ii0, ii3, wAt.mul(D));
        add_block(ii1, ii3, wBt.mul(D));
    }
 
    if ( modulo )
    {
        //this was not tested!
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            add_base(ii0, wAt*off);
            add_base(ii1, wBt*off);
            add_base(ii2, wCt*off);
            add_base(ii3, wDt*off);
        }
    }
    DRAW_LINK(ptA, cross(armA, ptA.dir()), cross(armB, ptB.dir()), ptB.pos());
}


/**
 Link `ptA` (A) and `ptB` (B),
 but the links are made between SA and SB which are located
 on the side of A and B, respectively:
 
     SA = A + len * N_A,
     SB = B + len * N_B,
 
 N_{A/B} is a normalized vector orthogonal to the fiber carrying each point:
 The force is linear of zero resting length,
 
     force_SA = weight * ( SA - SB )
     force_SB = weight * ( SB - SA )
 
 */

void Meca::addSideSideLink(Interpolation const& ptA,
                           Interpolation const& ptB,
                           const real len,
                           const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideSideLink meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector dir = position_diff(ptA, ptB);
    real armA = std::copysign(0.5*len, cross(ptA.diff(), dir));
    real armB = std::copysign(0.5*len, cross(dir, ptB.diff()));
    addSideSideLink(ptA, armA, ptB, armB, weight);

#else
    
    Vector dir = position_diff(ptA, ptB);
    Vector armA = cross(ptA.diff(), dir).normalized(0.5*len);
    Vector armB = cross(dir, ptB.diff()).normalized(0.5*len);
    addSideSideLink(ptA, armA, ptB, armB, weight);

#endif
}

//------------------------------------------------------------------------------
#pragma mark - Frictionless Links
//------------------------------------------------------------------------------

/**
 Link `ptA` (A) and `ptB` (B),
 The force is linear of zero resting length, but is anisotropic:
 The component of the force parallel to the fiber in A is removed

 If T is the normalized direction of the fiber in A:
 
     force_A = weight * ( 1 - T T' ) ( A - B )
     force_B = weight * ( 1 - T T' ) ( B - A )
 
 */

void Meca::addSlidingLink(Interpolation const& ptA,
                          Mecapoint const& ptB,
                          const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();
    
    if ( any_equal(ii2, ii0, ii1) )
        return;
    
    // interpolation coefficients:
    const real A = ptA.coef0();
    const real B = ptA.coef1();
    
    Vector dir = ptA.dir();

    /*
     Points are (a, b, e) with (ab) the Interpolation, and e the Mecapoint,
     P is the projection on the plane perpendicular to (ab)
         P = I - dir (x) dir / normSqr(dir)
     the interaction is  -weigth * transpose(bb, aa, -1) * P * ( bb, aa, -1 )
     we set only the upper part of this symmetric matrix:
     */

    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, A);
            add_base(ii1, off, B);
            sub_base(ii2, off);
        }
    }
    DRAW_LINK(ptA, ptB.pos());

    add_block_diag(ii0, A*A, wT);
    add_block_diag(ii1, B*B, wT);
    add_block_diag(ii2, wT);

    add_block(ii1, ii0, A*B, wT);
    if ( ii0 < ii2 )
    {
        sub_block(ii2, ii0, A, wT);
        sub_block(ii2, ii1, B, wT);
    }
    else
    {
        sub_block(ii0, ii2, A, wT);
        sub_block(ii1, ii2, B, wT);
    }
}


/**
Link `ptA` (A) and `ptB` (B),
 The force is linear of zero resting length, but is anisotropic:
 The component of the force parallel to the fiber in A is removed
 
 If T is the normalized direction of the fiber in A:
 
     force_A = weight * ( 1 - T T' ) ( A - B )
     force_B = weight * ( 1 - T T' ) ( B - A )
 
 */

void Meca::addSlidingLink(Interpolation const& ptA,
                          Interpolation const& ptB,
                          const real weight)
{
    assert_true( weight >= 0 );
 
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    //coefficients to form A-B:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    
    // on points (a, b, e), (ab) being the Interpolation, and e the Mecapoint,
    // P is the projection on the plane perpendicular to (ab): P.v= (v - (T.v)T/normSqr(T))
    // the interaction is  -wh' * P * h
    // we set only the upper part of this symmetric matrix:
    
    Vector dir = ptA.dir();

    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            add_base(ii2, off, cc2);
            add_base(ii3, off, cc3);
        }
    }
    DRAW_LINK(ptA, ptB.pos());

    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block_diag(ii2, cc2*cc2, wT);
    add_block(ii3, ii2, cc3*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);
}

//------------------------------------------------------------------------------
#pragma mark - Off-axis frictionless links (used for steric interactions)


/**
 Alternative 2D method in which we add an offset to vBAS
 Vector 'arm' must be parallel to the link and orthogonal to 'ptA'

  Without tangential force, a 'long link' is in the perpendicular direction.
  In the local reference frame, the matrix of interaction coefficients would be:
  real T[9] = { 0, 0, 0, 0, -weight, 0, 0, 0, 0 };
  we could transform it with a change-of-coordinates matrix R:
  Vector a = ptA.dir();
  Vector b = dir;
  Vector c = cross(a, b);
  real R[9] = { a.XX, a.YY, a.ZZ, b.XX, b.YY, b.ZZ, c.XX, c.YY, c.ZZ };
  real TR[3*3];
  blas::xgemm('N','T', 3, 3, 3, 1.0, T, 3, R, 3, 0.0, TR, 3);
  blas::xgemm('N','N', 3, 3, 3, 1.0, R, 3, TR, 3, 0.0, T, 3);
  equivalently, we can set directly the interaction coefficient matrix:
  */
void Meca::addSideSlidingLinkS(Interpolation const& ptA,
                               Mecapoint const& ptB,
                               Torque const& arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    // indices
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();
    
    if ( any_equal(ii2, ii0, ii1) )
        return;

    // interpolation coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();

    // wT = -weight * [ I - dir (x) dir ]
#if ( DIM == 3 )
    MatrixBlock wT = MatrixBlock::outerProduct(arm, -weight/arm.normSqr());
    Vector warm = -weight * arm;
#elif ( DIM == 2 )
    // set vector 'axi' perpendicular to Fiber
    real iseg = static_cast<Fiber const*>(ptA.mecable())->segmentationInv();
    Vector axi = cross(iseg, ptA.diff());
    MatrixBlock wT = MatrixBlock::outerProduct(axi, -weight);
    Vector warm = axi * ( -weight / arm );
#else
    MatrixBlock wT(0,0);
    Vector warm(0,0,0);
#endif

    add_base(ii0, warm, cc0);
    add_base(ii1, warm, cc1);
    sub_base(ii2, warm);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            sub_base(ii2, off);
        }
    }
    DRAW_LINK(ptA, cross(arm, ptA.dir()), ptB.pos());

    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block_diag(ii2, wT);

    add_block(ii1, ii0, cc1*cc0, wT);
    if ( ii0 < ii2 )
    {
        sub_block(ii2, ii0, cc0, wT);
        sub_block(ii2, ii1, cc1, wT);
    }
    else
    {
        sub_block(ii0, ii2, cc0, wT);
        sub_block(ii1, ii2, cc1, wT);
    }
}


#if ( DIM == 2 )
/**
 It is assumed that: norm(leg) = len / ptA.segmentation()
*/
void Meca::addSideSlidingLink2D(Interpolation const& ptA,
                                const real leg,
                                Mecapoint const& ptB,
                                Vector const& dir,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();
    
    if ( any_equal(ii2, ii0, ii1) )
        return;
    
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    
    // the projection matrix: wP = -weight * [ I - dir (x) dir ]
    MatrixBlock wP = MatrixBlock::offsetOuterProduct(-weight, dir, weight);
    
    // anti-symmetric matrix blocks:
    const Matrix22 A(cc0, -leg,  leg, cc0);
    const Matrix22 B(cc1,  leg, -leg, cc1);

    /*
     We use block operations to set the matrix block by block:
     | A'PA  A'PB  A'P |
     | B'PA  B'PB  B'P |
     |   PA    PB    P |
     This matrix has symmetric and anti-symmetric blocks,
     since P' = P whereas A and B are anti-symmetric
     */
    assert_true( ii0 < ii1 );
    
    if ( ii0 < ii2 )
    {
        const Matrix22 wPA = wP.mul(A);
        const Matrix22 wPB = wP.mul(B);
        add_block_diag(ii0, A.trans_mul(wPA));
        add_block(ii1, ii0, B.trans_mul(wPA));
        sub_block(ii2, ii0, wPA);
        add_block_diag(ii1, B.trans_mul(wPB));
        sub_block(ii2, ii1, wPB);
        add_block_diag(ii2, wP);
    }
    else
    {
        // in this case, swap indices to address lower triangle
        const Matrix22 AtwP = A.trans_mul(wP);
        const Matrix22 BtwP = B.trans_mul(wP);
        add_block_diag(ii2, wP);
        sub_block(ii0, ii2, AtwP);
        sub_block(ii1, ii2, BtwP);
        add_block_diag(ii0, AtwP.mul(A));
        add_block(ii1, ii0, BtwP.mul(A));
        add_block_diag(ii1, BtwP.mul(B));
    }
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            //off = -weight * ( off - dot(off, dir) * dir );
            off = wP * off;
            add_base(ii0, A.trans_vecmul(off));
            add_base(ii1, B.trans_vecmul(off));
            sub_base(ii2, off);
        }
    }
    DRAW_LINK(ptA, cross(leg, ptA.diff()), ptB.pos());
}

#endif


/**
 It is assumed that: norm(leg) = len / ptA.segmentation()
 */
void Meca::addSideSlidingLink3D(Interpolation const& ptA,
                                Torque const& leg,
                                Mecapoint const& ptB,
                                Vector const& dir,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex0();
    
    if ( any_equal(ii2, ii0, ii1) )
        return;
    
    // interpolation coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    
    MatrixBlock aR = MatrixBlock::vectorProduct(cc0, -leg);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1,  leg);

    // the projection matrix: P = -weight * [ I - dir (x) dir ]
    MatrixBlock wP = MatrixBlock::offsetOuterProduct(-weight, dir, weight);
    
    MatrixBlock aTwP = aR.trans_mul(wP);
    MatrixBlock bTwP = bR.trans_mul(wP);
    
    add_block_diag(ii0, aTwP.mul(aR));
    add_block(ii1, ii0, bTwP.mul(aR));
    add_block_diag(ii1, bTwP.mul(bR));
    if ( ii0 < ii2 )
    {
        sub_block(ii2, ii0, aTwP.transposed());
        sub_block(ii2, ii1, bTwP.transposed());
    }
    else
    {
        sub_block(ii0, ii2, aTwP);
        sub_block(ii1, ii2, bTwP);
    }
    add_block_diag(ii2, wP);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            add_base(ii0, aTwP*off);
            add_base(ii1, bTwP*off);
            sub_base(ii2, wP*off);
        }
    }
    DRAW_LINK(ptA, cross(leg, ptA.diff()), ptB.pos());
}


/// return a vector of norm 1.0, perpendicular to 'diff' and aligned with `off`:
Vector calculateArm(Vector off, Vector const& diff, real len)
{
    if ( modulo )
        modulo->fold(off);
    // remove component parallel to diff:
    off -= ( dot(off, diff) / diff.normSqr() ) * diff;
    real n = off.norm();
    if ( n > REAL_EPSILON )
        return off * ( len / n );
    else
        return diff.orthogonal(len);
}

/**
 Link `segA` (A) and `ptB` (B),
 This is a combination of a SideLink with a Sliding Link:
 The force is linear of zero resting length, but it is taken between B,
 and another point S located on the side of A:
 
     S = A + len * N,
 
 where N is a normalized vector orthogonal to the fiber in A, in the direction of B,
 which is obtained by cross product with the segment direction.
 In addition, the tangential part of the force is removed.
 
 If T is the normalized direction of the fiber in A:
 
     force_S = weight * ( 1 - T T' ) ( S - B )
     force_B = weight * ( 1 - T T' ) ( B - S )

 */
void Meca::addSideSlidingLink(FiberSegment const& segA, real abs,
                              Mecapoint const& ptB,
                              const real len,
                              const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideSlidingLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector dir = segA.dir();
    Vector ab = ptB.pos() - segA.pos1();
    if ( modulo )
        modulo->fold(ab);
    real leg = std::copysign(len*segA.lenInv(), cross(dir, ab));
    addSideSlidingLink2D(Interpolation(segA, abs), leg, ptB, dir, weight);
    
#else
    
    Vector dir = segA.dir();
    Vector leg = ptB.pos() - segA.pos1();
    if ( modulo )
        modulo->fold(leg);
    leg = cross(dir, leg);
    real n = leg.norm();
    if ( n > REAL_EPSILON )
    {
        leg *= ( len * segA.lenInv() ) / n;
        addSideSlidingLink3D(Interpolation(segA, abs), leg, ptB, dir, weight);
    }
    
#endif
}


#pragma mark - More off-axis frictionless links


/**
 Older code
 Vector 'arm' must be parallel to the link and orthogonal to 'ptA'
 */
void Meca::addSideSlidingLinkS(Interpolation const& ptA,
                               Interpolation const& ptB,
                               Torque const& arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    assert_false( ptA.invalid() );
    assert_false( ptB.invalid() );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    // coefficients to form A-B
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
   
    // wT = -weight * [ I - dir (x) dir ]
#if ( DIM == 3 )
    MatrixBlock wT = MatrixBlock::outerProduct(arm, -weight/arm.normSqr());
    Vector warm = (-weight) * arm;
#elif ( DIM == 2 )
    // set vector 'axi' perpendicular to Fiber:
    real iseg = static_cast<Fiber const*>(ptA.mecable())->segmentationInv();
    Vector axi = cross(iseg, ptA.diff());
    MatrixBlock wT = MatrixBlock::outerProduct(axi, -weight);
    Vector warm = axi * ( -weight / arm );
#else
    MatrixBlock wT(0,0);
    Vector warm(0,0,0);
#endif
    
    add_base(ii0, warm, cc0);
    add_base(ii1, warm, cc1);
    add_base(ii2, warm, cc2);
    add_base(ii3, warm, cc3);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            add_base(ii2, off, cc2);
            add_base(ii3, off, cc3);
        }
    }
    DRAW_LINK(ptA, cross(arm, ptA.dir()), ptB.pos());

    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block_diag(ii2, cc2*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);

    add_block(ii1, ii0, cc1*cc0, wT);
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block(ii3, ii2, cc3*cc2, wT);
}


#if ( DIM == 2 )

/**
 Link `ptA` (A) and `ptB` (B),
 This is a combination of a SideLink with a Sliding Link:
 The force is linear of zero resting length, but it is taken between B,
 and another point S located on the side of A:
 
     S = A + len * N,
 
 where N is a normalized vector orthogonal to the fiber in A, in the direction of B.
 In addition, the tangential part of the force is removed.
 
 If T is the normalized direction of the fiber in A:
 
     force_S = weight * ( 1 - T T' ) ( S - B )
     force_B = weight * ( 1 - T T' ) ( B - S )

 It is assumed that: norm(leg) = len / ptA.segmentation()
 */
void Meca::addSideSlidingLink2D(Interpolation const& ptA,
                                const real leg,
                                Interpolation const& ptB,
                                Vector const& dir,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    // coefficients to form A-B
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
        
    Matrix22 aR(cc0,-leg, leg, cc0);  // aR = alpha - len * R
    Matrix22 bR(cc1, leg,-leg, cc1);  // bR = beta + len * R
    
    // the projection matrix: P = -weight * [ I - dir (x) dir ]
    Matrix22 wP = Matrix22::offsetOuterProduct(-weight, dir, weight);
    
    Matrix22 aTwP = aR.trans_mul(wP);
    Matrix22 bTwP = bR.trans_mul(wP);

    add_block_diag(ii0, aTwP*aR);
    add_block(ii1, ii0, bTwP*aR);
    add_block_diag(ii1, bTwP*bR);
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, cc2, aTwP.transposed());
        add_block(ii3, ii0, cc3, aTwP.transposed());
        add_block(ii2, ii1, cc2, bTwP.transposed());
        add_block(ii3, ii1, cc3, bTwP.transposed());
    }
    else
    {
        add_block(ii0, ii2, cc2, aTwP);
        add_block(ii1, ii2, cc2, bTwP);
        add_block(ii0, ii3, cc3, aTwP);
        add_block(ii1, ii3, cc3, bTwP);
    }
    add_block_diag(ii2, cc2*cc2, wP);
    add_block(ii3, ii2, cc3*cc2, wP);
    add_block_diag(ii3, cc3*cc3, wP);

    if ( modulo )
    {
        throw Exception("addSideSlidingLink2D untested with periodic boundary conditions");
        
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            add_base(ii0, aTwP*off);
            add_base(ii1, bTwP*off);
            add_base(ii2, wP*off, cc2);
            add_base(ii3, wP*off, cc3);
        }
    }
    DRAW_LINK(ptA, cross(leg, ptA.dir()), ptB.pos());
}

#endif


/**
 Link `ptA` (A) and `ptB` (B),
 This is a combination of a SideLink with a Sliding Link:
 The force is linear of zero resting length, but it is taken between B,
 and another point S located on the side of A:
 
     S = A + cross( leg, T )
 
 where T is the unit vector tangent to the fiber in A, towards the plus-end.
 In addition, the tangential part of the force is removed.
 
 If T is the normalized direction of the fiber in A:
 
     force_S = weight * ( 1 - T T' ) ( S - B )
     force_B = weight * ( 1 - T T' ) ( B - S )

 It is assumed that: norm(leg) = len / ptA.segmentation()
*/
void Meca::addSideSlidingLink3D(Interpolation const& ptA,
                                Torque const& leg,
                                Interpolation const& ptB,
                                Vector const& dir,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii2, ii3) || any_equal(ii1, ii2, ii3) )
        return;
    
    // coefficients to form A-B
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    
    MatrixBlock aR = MatrixBlock::vectorProduct(cc0, -leg);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1,  leg);

    // the projection matrix: P = -weight * [ I - dir (x) dir ]
    MatrixBlock wP = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    MatrixBlock aTwP = aR.trans_mul(wP);
    MatrixBlock bTwP = bR.trans_mul(wP);
    
    add_block_diag(ii0, aTwP.mul(aR));
    add_block(ii1, ii0, bTwP.mul(aR));
    add_block_diag(ii1, bTwP.mul(bR));
    if ( ii0 < ii2 )
    {
        add_block(ii2, ii0, cc2, aTwP.transposed());
        add_block(ii3, ii0, cc3, aTwP.transposed());
        add_block(ii2, ii1, cc2, bTwP.transposed());
        add_block(ii3, ii1, cc3, bTwP.transposed());
    }
    else
    {
        add_block(ii0, ii2, cc2, aTwP);
        add_block(ii1, ii2, cc2, bTwP);
        add_block(ii0, ii3, cc3, aTwP);
        add_block(ii1, ii3, cc3, bTwP);
    }
    add_block_diag(ii2, cc2*cc2, wP);
    add_block(ii3, ii2, cc3*cc2, wP);
    add_block_diag(ii3, cc3*cc3, wP);
    
    if ( modulo )
    {
        Vector off = modulo_offset(ptA, ptB);
        if ( off.is_not_zero() )
        {
            add_base(ii0, aTwP*off);
            add_base(ii1, bTwP*off);
            add_base(ii2, wP*off, cc2);
            add_base(ii3, wP*off, cc3);
        }
    }
    DRAW_LINK(ptA, cross(leg, ptA.diff()), ptB.pos());
}


/**
 Link `segA` (A) and `ptB` (B),
 This is a combination of Side- and Sliding Links:
 The force is linear of zero resting length, but it is taken between B
 and another point S which is located on the side of A:
 S = A + len * N,
 where N is a normalized vector orthogonal to the fiber in A, in the direction of B.
 In addition, the part of the force tangential to A is removed.
 
 If T is the normalized direction of the fiber in A:
 
     force_S = weight * ( 1 - T T' ) ( S - B )
     force_B = weight * ( 1 - T T' ) ( B - S )
 
 */
void Meca::addSideSlidingLink(FiberSegment const& segA, real abs,
                              Interpolation const& ptB,
                              const real len,
                              const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideSlidingLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector dir = segA.dir();
    Vector ab = ptB.pos() - segA.pos1();
    if ( modulo )
        modulo->fold(ab);
    real leg = std::copysign(len*segA.lenInv(), cross(dir, ab));
    addSideSlidingLink2D(Interpolation(segA, abs), leg, ptB, dir, weight);
    
#else

    Vector dir = segA.dir();
    Vector leg = ptB.pos() - segA.pos1();
    if ( modulo )
        modulo->fold(leg);
    leg = cross(dir, leg);
    real n = leg.norm();
    if ( n > REAL_EPSILON )
    {
        leg *= ( len * segA.lenInv() ) / n;
        addSideSlidingLink3D(Interpolation(segA, abs), leg, ptB, dir, weight);
    }
    
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Links to fixed positions
//------------------------------------------------------------------------------

/**
 Link `ptA` (A) and a fixed position `pos` (G)
 The force is linear:
 
     force_A = weight * ( G - A );
 
 There is no counter-force in G, since G is immobile.
 */

void Meca::addPointClamp(Mecapoint const& ptA,
                         Vector pos,
                         const real weight)
{
    assert_true( weight >= 0 );
    const index_t ii0 = ptA.matIndex0();
    
    sub_iso_diag(ii0, weight);
    
    if ( modulo )
        modulo->fold(pos, ptA.pos());

    add_base(ii0, pos, weight);
    
    DRAW_LINK(ptA, pos);
}


/**
 Link `pti` (A) and a fixed position `pos` (G)
 The force is linear:
 
     force_A = weight * ( G - A );
 
 The point G is not associated to a Mecable, and there is no counter-force in G.
 */

void Meca::addPointClamp(Interpolation const& pti,
                         Vector pos,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    
    const real cc0 = pti.coef0();
    const real cc1 = pti.coef1();
    
    const real ww0 = weight * cc0;
    const real ww1 = weight * cc1;

    sub_iso_diag(ii0, ww0*cc0);
    sub_iso_diag(ii1, ww1*cc1);
    sub_iso(ii1, ii0, ww1*cc0);
    
    if ( modulo )
        modulo->fold(pos, pti.pos());
    
    add_base(ii0, pos, ww0);
    add_base(ii1, pos, ww1);
    
    DRAW_LINK(pti, pos);
}

/**
 This creates a Hookean link only in the X dimension
 */
void Meca::addPointClampX(Mecapoint const& ptA,
                          real x_pos,
                          const real weight)
{
    const index_t ii0 = ptA.matIndex0();

    MatrixBlock B(0, 0);
    B(0, 0) = weight;
    sub_block_diag(ii0, B);

    add_base(ii0, Vector(x_pos, 0, 0), weight);
}

/**
 This creates Hookean links only in the X and Y dimension
 in 2D this is equivalent to addPointClamp()
 in 3D there is no force in the Z, if pos.ZZ = 0
 */
void Meca::addPointClampXY(Mecapoint const& ptA,
                           Vector pos,
                           const real weight)
{
    const index_t ii0 = ptA.matIndex0();

#if ( DIM == 2 )
    sub_iso_diag(ii0, weight);
#elif ( DIM > 2 )
    MatrixBlock B(0, 0);
    B(0, 0) = weight;
    B(1, 1) = weight;
    sub_block_diag(ii0, B);
    assert_true( pos.ZZ == 0 );
#endif
    add_base(ii0, pos, weight);
}


void Meca::addPointClampToAll(Vector const& pos, const real weight)
{
    ABORT_NOW("Unfinished");
    /*
    Vector vec = weight * pos;
    for ( index_t p = 0; p < nbVertices(); ++p )
    {
        sub_iso_diag(p, weight);
        add_base(p, vec);
    }
     */
}

//------------------------------------------------------------------------------
#pragma mark - Links to fixed sphere
//------------------------------------------------------------------------------

/**
 Link `pte` (P) and a fixed sphere of radius `rad` and center `center` (C)
 The force is affine with non-zero resting length:

      force = weight * ( C - P ) * ( 1 - rad / |PC| )

 The constant part is:
 
      weight * ( C - P ) * ( 1 - rad / |PC| )

 for a point inside, this is directed outward.
 There is no force on the center C, which is an immobile position.
 */

void Meca::addSphereClamp(Vector const& off,
                          Mecapoint const& ptA,
                          Vector const& center,
                          real rad,
                          const real weight)
{
    assert_true( rad >= 0 );
    assert_true( weight >= 0 );
    
    const index_t ii0 = ptA.matIndex0();
    
    real len = off.norm();
    
    if ( len > REAL_EPSILON )
    {
        real wla = weight * rad / len;
        MatrixBlock wT;
        /* To stabilize the matrix with compression, we remove negative eigenvalues
         This is done by using len = 1 in the formula for links that are shorter
         than the desired target. */
        if ( rad < len )
            wT = MatrixBlock::offsetOuterProduct(wla-weight, off/len, -wla);
        else
            wT = MatrixBlock::outerProduct(off/len, -weight);
        
        add_block_diag(ii0, wT);
        add_base(ii0, wla*off-wT*center);
    }
}


void Meca::addSphereClamp(Vector const& off,
                          Interpolation const& ptA,
                          Vector const& center,
                          real rad,
                          const real weight)
{
    assert_true( rad >= 0 );
    assert_true( weight >= 0 );
    
    real len = off.norm();
    
    if ( len > REAL_EPSILON )
    {
        real wla = weight * rad / len;
        
        const index_t ii0 = ptA.matIndex1();
        const index_t ii1 = ptA.matIndex2();
        
        MatrixBlock wT;
        /* To stabilize the matrix with compression, we remove negative eigenvalues
         This is done by using len = 1 in the formula for links that are shorter
         than the desired target. */
        if ( rad < len )
            wT = MatrixBlock::offsetOuterProduct(wla-weight, off/len, -wla);
        else
            wT = MatrixBlock::outerProduct(off/len, -weight);

        // interpolation coefficients:
        const real cc0 = ptA.coef0();
        const real cc1 = ptA.coef1();

        add_block_diag(ii0, cc0*cc0, wT);
        add_block(ii1, ii0, cc1*cc0, wT);
        add_block_diag(ii1, cc1*cc1, wT);
        
        Vector vec = wla*off-wT*center;
        add_base(ii0, vec, cc0);
        add_base(ii1, vec, cc1);
    }
}


void Meca::addSphereClamp(Mecapoint const& pte,
                          Vector center,
                          real rad,
                          const real weight)
{
    Vector pos = pte.pos();
    if ( modulo )
        modulo->fold(center, pos);
    addSphereClamp(pos-center, pte, center, rad, weight);
}


void Meca::addSphereClamp(Interpolation const& pti,
                          Vector center,
                          real rad,
                          const real weight)
{
    Vector pos = pti.pos();
    if ( modulo )
        modulo->fold(center, pos);
    addSphereClamp(pos-center, pti, center, rad, weight);
}


//------------------------------------------------------------------------------
#pragma mark - Links to fixed cylinder
//------------------------------------------------------------------------------

/**
 Link `pte` (P) to a cylinder of axis X and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(P.XX, 0, 0)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the YZ plane.
 */
void Meca::addCylinderClampX(Mecapoint const& pte, const real rad, const real weight)
{
    assert_true( weight >= 0 );
    const index_t ii0 = pte.matIndex0();
    
#if ( DIM == 2 )
    
    //mFUL(inx+1, inx+1) -= weight;
    mFUL.diag_block(ii0)(1, 1) -= weight;
    vBAS[DIM*ii0+1] += weight * std::copysign(rad, pte.pos().YY);
    
#elif ( DIM >= 3 )

    Vector pos = pte.pos();
    real len = pos.normYZ();
    real dY, dZ;
        
    // we use: real wla = weight * min_real(rad/len, 1);
    MatrixBlock B(0, 0);
    if ( rad < len )
    {
        dY = pos.YY / len;
        dZ = pos.ZZ / len;
        // attractive
        real wla = weight * rad / len;
        B(1, 1) = wla * ( dY * dY - 1.0 ) + weight;
        B(2, 1) = wla * dZ * dY;
        B(2, 2) = wla * ( dZ * dZ - 1.0 ) + weight;
    }
    else if ( len > REAL_EPSILON )
    {
        dY = pos.YY / len;
        dZ = pos.ZZ / len;
        // repulsive
        B(1, 1) = weight * dY * dY;
        B(2, 1) = weight * dZ * dY;
        B(2, 2) = weight * dZ * dZ;
    }
    else
        return;
    
    B(1,2) = B(2,1);
    sub_block_diag(ii0, B);

#if ( 0 )
    real ab2 = pos.normYZSqr();
    /// offsetOuterProduct(dia, dir, len) = dia * Id + [ dir (x) dir ] * len
    if ( len < rad )
        B = MatrixBlock::outerProduct(pos, weight/ab2);
    else
    {
        real wla = weight * rad / len;
        B = MatrixBlock::offsetOuterProduct(weight-wla, pos, wla/ab2);
    }
    std::clog << B;
#endif

    real fac = weight * rad;
    // there should be no X component here!
    unsigned inx = ii0;
    vBAS[DIM*inx+1] += fac * dY;
    vBAS[DIM*inx+2] += fac * dZ;
    
#endif
}


/**
 Link `pte` (P) to a cylinder of axis Y and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, P.YY, 0)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the XZ plane.
 */
void Meca::addCylinderClampY(Mecapoint const& pte, const real rad, const real weight)
{
    assert_true( weight >= 0 );
    const index_t ii0 = pte.matIndex0();
    
#if ( DIM == 2 )
    
    //mFUL(inx, inx) -= weight;
    mFUL.diag_block(ii0)(0, 0) -= weight;
    vBAS[DIM*ii0] += weight * std::copysign(rad, pte.pos().XX);
    
#elif ( DIM >= 3 )

    Vector pos = pte.pos();
    real len = pos.normXZ();
    real dX, dZ;
    
    // we use: real wla = weight * min_real(rad/len, 1);
    MatrixBlock B(0, 0);
    if ( rad < len )
    {
        dX = pos.XX / len;
        dZ = pos.ZZ / len;
        // attractive
        real wla = weight * rad / len;
        B(0, 0) = wla * ( dX * dX - 1.0 ) + weight;
        B(2, 0) = wla * dZ * dX;
        B(2, 2) = wla * ( dZ * dZ - 1.0 ) + weight;
    }
    else if ( len > REAL_EPSILON )
    {
        dX = pos.XX / len;
        dZ = pos.ZZ / len;
        // repulsive
        B(0, 0) = weight * dX * dX;
        B(2, 0) = weight * dZ * dX;
        B(2, 2) = weight * dZ * dZ;
    }
    else
        return;
    
    B(0,2) = B(2,0);
    sub_block_diag(ii0, B);

    real fac = weight * rad;
    unsigned inx = ii0;
    vBAS[DIM*inx  ] += fac * dX;
    // there should be no Y component here!
    vBAS[DIM*inx+2] += fac * dZ;
    
#endif
}


/**
 Link `pte` (P) to a cylinder of axis Z and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, 0, P.ZZ)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the XY plane.
 */
void Meca::addCylinderClampZ(Mecapoint const& pte, const real rad, const real weight)
{
    assert_true( weight >= 0 );
    
#if ( DIM >= 3 )

    const index_t ii0 = pte.matIndex0();
    Vector pos = pte.pos();
    real len = pos.normXY();
    real dX, dY;
    
    // we use: real wla = weight * min_real(rad/len, 1);
    MatrixBlock B(0, 0);
    if ( rad < len )
    {
        dX = pos.XX / len;
        dY = pos.YY / len;
        // attractive
        real wla = weight * rad / len;
        B(0, 0) = wla * ( dX * dX - 1.0 ) + weight;
        B(1, 0) = wla * dY * dX;
        B(1, 1) = wla * ( dY * dY - 1.0 ) + weight;
    }
    else if ( len > REAL_EPSILON )
    {
        dX = pos.XX / len;
        dY = pos.YY / len;
        // repulsive
        B(0, 0) = weight * dX * dX;
        B(1, 0) = weight * dY * dX;
        B(1, 1) = weight * dY * dY;
    }
    else
        return;

    B(0,1) = B(1,0);
    sub_block_diag(ii0, B);

    real fac = weight * rad;
    unsigned inx = ii0;
    vBAS[DIM*inx  ] += fac * dX;
    vBAS[DIM*inx+1] += fac * dY;
    // there should be no Z component here!

#endif
}

/**
 Link `pte` (P) to a cylinder of given `axis` and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, 0, P.ZZ)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the plane orthogonal to `axis` anchored at `center`.
 */
void Meca::addCylinderClamp(Mecapoint const& pte,
                            Vector const& axis, Vector const& center,
                            const real rad, const real weight)
{
#if ( DIM > 2 )
    
    assert_true( weight >= 0 );
    const index_t ii0 = pte.matIndex0();

    //Projection removing components parallel to the cylinder axis:  P = [ I - axis (x) axis ]
    MatrixBlock P = MatrixBlock::offsetOuterProduct(1.0, axis, -1.0/axis.normSqr());
    
    Vector dir = P.vecmul( pte.pos() - center );
    real len = dir.norm();
    
    if ( len > REAL_EPSILON )
    {
        real wla = weight * rad / len;
        
        MatrixBlock W;
        if ( rad < len )
            W = MatrixBlock::offsetOuterProduct(wla-weight, dir/len, -wla);
        else
            W = MatrixBlock::outerProduct(dir/len, -weight);
        P = P.mul(W);
        
        add_block_diag(ii0, P);
        add_base(ii0, wla*dir - P*center);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Off-axis links to fixed positions
//------------------------------------------------------------------------------


#if ( DIM == 2 )

/**
 This assumes that Modulo has been applied, such that 0 == modulo->offset( ptA.pos() - pos );
 */
void Meca::addSidePointClamp2D(Interpolation const& ptA,
                               Vector pos,
                               const real arm,
                               const real weight)
{
    // interpolation coefficients:
    const real A = ptA.coef0(),  wA = weight * A;
    const real B = ptA.coef1(),  wB = weight * B;
    
    const real E = arm * ptA.lenInv();
    const real wE = weight * E;
    const real wEE = weight * E * E;
    
    // indices:
    index_t ii0 = ptA.matIndex1();
    index_t ii1 = ptA.matIndex2();
    
    //isotropic terms:
    sub_iso_diag(ii0,  wA * A + wEE);
    sub_iso_diag(ii1,  wB * B + wEE);
    sub_iso(ii1, ii0,  wA * B - wEE);
    
    // non-isotropic:
    mFUL(DIM*ii0  , DIM*ii1+1) += wE;
    mFUL(DIM*ii0+1, DIM*ii1  ) -= wE;

    //it seems to works also fine without the term in eew* below:
    Vector off = wE * Vector(-pos.YY, pos.XX);

    add_base(ii0, wA*pos+off);
    add_base(ii1, wB*pos-off);
    
    DRAW_LINK(ptA, cross(arm, ptA.dir()), pos);
}

#endif

/**
 A link of stiffness `weight`, between `offset_point` on the side of `ptA`,
 and a fixed position `pos` (G).

 This uses the vector product to offset the point on which the link is attached:
 
     offset_point = fiber_point + cross(arm, fiber_dir)
 
 with fiber_point = position(ptA) and fiber_dir = ptA.diff().normalized.
 `arm` must be perpendicular to link ( G - position(ptA) )

 F. Nedelec, March 2011, April 2024

 This assumes that Modulo has been applied, such that 0 == modulo->offset( ptA.pos() - pos );
 
 */
void Meca::addSidePointClamp3D(Interpolation const& ptA,
                               Vector pos,
                               Torque const& arm,
                               real const weight)
{
    assert_true( weight >= 0 );
    // indices:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    
    // interpolation coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    const Torque leg = arm * ptA.lenInv();
    
#if ( DIM == 1 )
    // this implementation works in any dimension, but is slower than the one below
    MatrixBlock aR = MatrixBlock::vectorProduct(cc0, -leg);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1,  leg);

    MatrixBlock At = aR.transposed();
    MatrixBlock Bt = bR.transposed();

    // the diagonal blocs are symmetric but not diagonal
    add_block_diag(ii0, -weight, At.mul(aR));
    add_block(ii1, ii0, -weight, Bt.mul(aR)); // not symmetric
    add_block_diag(ii1, -weight, Bt.mul(bR));
#else
    /* Since aR and bR only differ on their diagonals, the calculation above
     can be simplified: aR = cc0 - leg (x) Id; bR = cc1 + leg (x) Id 
     We used:
         cc0 + cc1 = 1
         cross(L,Id)^2 = L(x)L - normSqr(L)*Id
     and:
     
     */
#if ( DIM >= 3 )
    MatrixBlock LL = MatrixBlock::offsetOuterProduct(leg);
#else
    MatrixBlock LL = MatrixBlock(0, -leg * leg);
#endif
    MatrixBlock MM = MatrixBlock::vectorProduct(-cc0*cc1, leg);

    //std::clog<<-At.mul(aR)<<" "<<-Bt.mul(bR)<<" "<<-Bt.mul(aR)<<"\n";
    //std::clog<<LL+MatrixBlock(0,-cc0*cc0)<<" "<<LL+MatrixBlock(0,-cc1*cc1)<<" "<<MM-LL<<"/\n\n";
    //std::clog << ii0 << "  " << ii1 << "\n";
    
    // the diagonal blocs are symmetric but not diagonal
    add_block_diag(ii0, weight, LL, -cc0*cc0);
    add_block(ii1, ii0, weight, MM-LL);
    add_block_diag(ii1, weight, LL, -cc1*cc1);
    
#endif
 
    //add_base(ii0, At*pos, weight); // At * pos = cc0 * pos + cross(len, pos)
    //add_base(ii1, Bt*pos, weight); // Bt * pos = cc1 * pos - cross(len, pos)

    Vector vec = cross(leg, pos);
    add_base(ii0, vec + cc0*pos, weight);
    sub_base(ii1, vec - cc1*pos, weight);

    DRAW_LINK(ptA, cross(arm, ptA.dir()), pos);
}


/**
 Update Meca to include a connection between `ptA` (A) and a fixed position `pos` (G).
 The force is of zero resting length, but it is taken between G
 and another point S which is located on the side of the segment supporting A:
 
     S = A + len * N,
     force_S = weight * ( G - S )
 
 There is no counter-force in G, since G is immobile.
 */

void Meca::addSidePointClamp(Interpolation const& ptA,
                             Vector pos,
                             const real len,
                             const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSidePointClamp is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    // 'arm' is a pseudo-vector in the Z direction
    real arm = std::copysign(len, cross(ptA.diff(), pos-ptA.pos1()));
    addSidePointClamp2D(ptA, pos, arm, weight);
   
#else
    
    // 'arm' perpendicular to link and fiber is obtained by vector product:
    Vector arm = ptA.pos1() - pos;
    if ( modulo )
        pos += modulo->foldOffset(arm);
    arm = cross(arm, ptA.diff());
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSidePointClamp3D(ptA, pos, arm*(len/n), weight);

#endif  
}

//------------------------------------------------------------------------------
#pragma mark - Links to fixed lines and planes
//------------------------------------------------------------------------------

/**
 Link `ptA` to the X-axis with a stiffness `weight`.
 The component parallel to `X` are removed corresponding to a frictionless line
 
     force_Y = -weight * Y
     force_Z = -weight * Z

 */
void Meca::addLineClampX(Mecapoint const& pte, const real weight)
{
    assert_true( weight >= 0 );
    const index_t ii0 = pte.matIndex0();
    
#if ( DIM == 3 )
    MatrixBlock wT(0, 0, 0, 0, -weight, 0, 0, 0, -weight);
    add_block_diag(ii0, wT);
#elif ( DIM == 2 )
    MatrixBlock wT(0, 0, 0, -weight);
    add_block_diag(ii0, wT);
#endif
}


/**
 Link `ptA` (X) to the line defined by `G` and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 component parallel to `dir` are removed corresponding to a frictionless line:
 
     matrix M = 1 - dir (x) dir'
     force_X = weight * M * ( G - X )
 
 The vector `dir` should be of norm = 1.
 */

void Meca::addLineClamp(Mecapoint const& pte,
                        Vector const& pos,
                        Vector const& dir,
                        const real weight)
{
    assert_true( weight >= 0 );
    
    const index_t ii0 = pte.matIndex0();
    
    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_block_diag(ii0, wT);
    sub_base(ii0, wT*pos);
}


/**
 Link `ptA` and the line defined by `pos` (C) and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 components parallel to `dir` are removed corresponding to a frictionless line:
 
     M = I - dir (x) dir'
     force = weight * M * ( C - P )
 
 The vector `dir` should be of norm = 1.
 */

void Meca::addLineClamp(Interpolation const& pti,
                        Vector const& pos,
                        Vector const& dir,
                        const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();

    // interpolation coefficients:
    const real cc0 = pti.coef0();
    const real cc1 = pti.coef1();

    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii1, ii0, cc0*cc1, wT);
    
    //add the constant term:
    sub_base(ii0, wT*pos, cc0);
    sub_base(ii1, wT*pos, cc1);
}


/**
 This constrains a single degree of freedom indicated by 'xyz = { 0, 1, 2 }',
 to a plane in 3D. Hence 'off' is the X, Y or Z component of the plane.
*/
void Meca::addPlaneClampXYZ(Mecapoint const& P, int xyz, real off, real weight)
{
    assert_true( weight >= 0 );
    assert_true( xyz < DIM );
#if USE_MATRIX_BLOCK
    index_t inx = P.matIndex0();
    mFUL.diag_block(inx)(xyz, xyz) -= weight;
    vBAS[DIM*inx+xyz] += weight * off;
#else
    index_t inx = DIM * P.matIndex0() + xyz;
    mFUL(inx, inx) -= weight;
    vBAS[inx] += weight * off;
#endif
}


/**
 Link `ptA` (X) and the plane defined one of its point `pos` and the normal `dir`.
 The force is linear and the components parallel to the plane are removed,
 corresponding to an interaction with a frictionless plane:
 
     projection matrix M = dir (x) dir'
     force_X = weight * M * ( pos - X )
 
 The vector `dir` should be of norm = 1, or alternatively
 `weight` can be the true weigth divided by |dir|^2.
 */

void Meca::addPlaneClamp(Mecapoint const& pte,
                         Vector const& pos,
                         Vector const& dir,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    index_t ii0 = pte.matIndex0();
    
    // vBAS[DIM*inx] += dir * ( weigth * dot(pos,dir) );
    add_base(ii0, dir, weight*dot(pos, dir));
    
#if ( DIM == 1 ) && USE_ISO_MATRIX
    mISO.diagonal(ii0) -= weight;
#else
    // wT = -weight * [ dir (x) dir ]
    MatrixBlock wT = MatrixBlock::outerProduct(dir, -weight);
    add_block_diag(ii0, wT);
#endif
}


/**
 Link `ptA` (X) and the plane defined by `pos` and the normal `dir`.
 The force is linear and the perpendicular forces are removed, to create a frictionless plane:
 
     projection matrix M = dir (x) dir'
     force = weight * M * ( pos - X )
 
 The vector `dir` should be of norm = 1, or alternatively 
 `weight` can be the true weigth divided by |dir|^2.
 */

void Meca::addPlaneClamp(Interpolation const& pti,
                         Vector const& pos,
                         Vector const& dir,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();

    // interpolation coefficients:
    const real cc0 = pti.coef0();
    const real cc1 = pti.coef1();
    
    // wT = -weight * [ dir (x) dir ]
    MatrixBlock wT = MatrixBlock::outerProduct(dir, -weight);
    
    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    
    //add the constant term:
    sub_base(ii0, wT*pos, cc0);
    sub_base(ii1, wT*pos, cc1);
}


//------------------------------------------------------------------------------
#pragma mark - Experimental interactions
//------------------------------------------------------------------------------

/**
Links { pt1, pt2, pt3 } to a dragless junction `X` with stiffness { w1, w2, w3 }.

                           pt2
                          /
                  pt1 -- X
                          \
                           pt3

The position of the virtual point `X` is always determined by force balance:

    0 = w1 * ( pt1 - X ) + w2 * ( pt2 - X ) + w3 * ( pt3 - X )

The force on `pt1` is then Hookean:

    w1 * ( X - pt1 )

and similarly for the other points.
 
We first derive:

    X = ( w1 * pt1 + w2 * pt2 + w3 * pt3 ) / sum

and for the first point:

    f1 = ( w1 / sum ) * [ w2 * ( pt2 - pt1 ) + w3 * ( pt3 - pt1 ) ]
*/
void Meca::addTriLink(Interpolation const& pt1, const real w1,
                      Interpolation const& pt2, const real w2,
                      Interpolation const& pt3, const real w3)
{
    const real sum = w1 + w2 + w3;
    assert_true( sum > REAL_EPSILON );
    addLink(pt1, pt2, w1*w2/sum);
    addLink(pt1, pt3, w1*w3/sum);
    addLink(pt2, pt3, w2*w3/sum);
}


