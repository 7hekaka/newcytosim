// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "space.h"
#include "simul.h"
#include "modulo.h"
#include "mecapoint.h"
#include "interpolation.h"
#include "matrix11.h"
#include "matrix22.h"
#include "matrix33.h"
#include "matrix44.h"

//#include "vecprint.h"

extern Modulo const* modulo;

/// set TRUE to update matrix mC using block directives
/** This is significantly faster on machine with the AVX instruction set */
#define USE_MATRIX_BLOCK 1

//------------------------------------------------------------------------------
#pragma mark - Accessory functions

#if DRAW_MECA_LINKS

#  include "gle.h"
#  include "gle_color_list.h"

/// this performs the modulo on `c`
void drawLink(Vector const& a, Vector const& ab, Vector c)
{
    if ( modulo ) modulo->fold(c, a);
    gle::drawLink(a, ab, c);
}

#endif


/// true if any two values are equal
bool any_equal(const index_t a, const index_t b,
               const index_t c)
{
    //if ( a == b ) return true;
    if ( a == c ) return true;
    if ( b == c ) return true;
    return false;
}


/// true if any two values are equal
bool any_equal(const index_t a, const index_t b,
               const index_t c, const index_t d)
{
    //if ( a == b ) return true;
    if ( a == c ) return true;
    if ( a == d ) return true;
    if ( b == c ) return true;
    if ( b == d ) return true;
    //if ( c == d ) return true;
    return false;
}


//------------------------------------------------------------------------------
#pragma mark - Matrix access functions
//------------------------------------------------------------------------------

inline void Meca::add_block(index_t i, index_t j, MatrixBlock const& T)
{
#if ( DIM == 1 )
    mB(i,j) += T.value();
#elif USE_MATRIX_BLOCK
    mC.block(i, j).add_full(T);
#else
    assert_true( i != j );
    for ( int x = 0; x < DIM; ++x )
    for ( int y = 0; y < DIM; ++y )
        mC(i+y, j+x) += T(y,x);
#endif
}

inline void Meca::add_block(index_t i, index_t j, real alpha, MatrixBlock const& T)
{
#if ( DIM == 1 )
    mB(i,j) += alpha * T.value();
#elif USE_MATRIX_BLOCK
    mC.block(i, j).add_full(alpha, T);
#else
    assert_true( i != j );
    for ( int x = 0; x < DIM; ++x )
    for ( int y = 0; y < DIM; ++y )
        mC(i+y, j+x) += alpha * T(y,x);
#endif
}

inline void Meca::sub_block(index_t i, index_t j, MatrixBlock const& T)
{
#if ( DIM == 1 )
    mB(i,j) -= T.value();
#elif USE_MATRIX_BLOCK
    mC.block(i, j).sub_full(T);
#else
    assert_true( i != j );
    for ( int x = 0; x < DIM; ++x )
    for ( int y = 0; y < DIM; ++y )
        mC(i+y, j+x) -= T(y,x);
#endif
}


inline void Meca::add_diag_block(index_t i, MatrixBlock const& T)
{
#if ( DIM == 1 )
    mB(i,i) += T.value();
#elif USE_MATRIX_BLOCK
    mC.diag_block(i).add_half(T);
#else
    // add lower part of block
    for ( int x = 0; x < DIM; ++x )
    for ( int y = x; y < DIM; ++y )
        mC(i+y, i+x) += T(y,x);
#endif
}

inline void Meca::add_diag_block(index_t i, real alpha, MatrixBlock const& T)
{
#if ( DIM == 1 )
    mB(i,i) += alpha * T.value();
#elif USE_MATRIX_BLOCK
    mC.diag_block(i).add_half(alpha, T);
#else
    // add lower part of block
    for ( int x = 0; x < DIM; ++x )
    for ( int y = x; y < DIM; ++y )
        mC(i+y, i+x) += alpha * T(y,x);
#endif
}

inline void Meca::sub_diag_block(index_t i, MatrixBlock const& T)
{
#if ( DIM == 1 )
    mB(i,i) -= T.value();
#elif USE_MATRIX_BLOCK
    mC.diag_block(i).sub_half(T);
#else
    // add lower part of block
    for ( int x = 0; x < DIM; ++x )
    for ( int y = x; y < DIM; ++y )
        mC(i+y, i+x) -= T(y,x);
#endif
}


inline void Meca::add_iso(index_t i, index_t j, real val)
{
#if ( 1 )
    mB(i,j) += val;
#else
    mC.block(DIM*i, DIM*j).add_diag(val);
#endif
}

inline void Meca::sub_iso(index_t i, index_t j, real val)
{
#if ( 1 )
    mB(i,j) -= val;
#else
    mC.block(DIM*i, DIM*j).sub_diag(val);
#endif
}


inline void Meca::add_base(index_t i, Vector const& vec)
{
    assert_true( i % DIM == 0 );
    vec.add_to(vBAS+i);
}

inline void Meca::add_base(index_t i, Vector const& vec, real alpha)
{
    assert_true( i % DIM == 0 );
    vec.add_to(alpha, vBAS+i);
}

inline void Meca::sub_base(index_t i, Vector const& vec)
{
    assert_true( i % DIM == 0 );
    vec.sub_to(vBAS+i);
}

//------------------------------------------------------------------------------
#pragma mark - Explicit (constant) Forces
//------------------------------------------------------------------------------


/**
 Add constant force to a vertex
 */
void Meca::addForce(const Mecapoint & pte, Vector const& force)
{
    const index_t inx = DIM * pte.matIndex();
    add_base(inx, force);
}


/**
Add constant force to an interpolated position
 */
void Meca::addForce(const Interpolation & pti, Vector const& force)
{
    const index_t ii0 = DIM * pti.matIndex1();
    const index_t ii1 = DIM * pti.matIndex2();
    
    add_base(ii0, force, pti.coef2());
    add_base(ii1, force, pti.coef1());
}


//------------------------------------------------------------------------------
#pragma mark - Explicit (constant) Torque
//------------------------------------------------------------------------------

/**
 Add constant torque in `pti`:
 
     force = cross(torque, position)
 
 This is explicit and all contributions go in the force vector vBAS[]
*/
void Meca::addTorque(const Interpolation & pti, const Torque & torque)
{
    const index_t ii0 = DIM * pti.matIndex1();
    const index_t ii1 = DIM * pti.matIndex2();
    
    Vector d = pti.diff();
    Vector f = cross(torque/d.normSqr(), d);
    
    sub_base(ii0, f);
    add_base(ii1, f);
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
void Meca::addTorqueClamp(const Interpolation & pti,
                          Vector const& dir,
                          const real weight)
{
    assert_true( weight >= 0 );
    const index_t ii0 = DIM * pti.matIndex1();
    const index_t ii1 = DIM * pti.matIndex2();
    
    Vector d = pti.diff();
    real n = d.normSqr();

    Torque Tq = cross(d, dir);

#if ( DIM >= 3 )
    const real Tn = Tq.norm();
#else
    const real Tn = fabs(Tq);
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = atan2(Tn, dot(d, dir));
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     knowing that Tn = Tq.norm() = sin(angle) * sqrt(n)
     
     To have a Torque proportional to sin(angle), use:
     real nn = weight / ( n * sqrt(n) );
     */
    real nn = weight * angle / ( n * Tn );
    
    Vector f = cross(Tq * nn, d);
    
    sub_base(ii0, f);
    add_base(ii1, f);
}


/**
 Add an explicit torque to bring two segments parallel to each other,
 with a given weight:
 
     torque = weigth * cross(dirA, dirB)
     forceA =  cross(torque, dirA)
     forceB = -cross(torque, dirB)
 
 This is explicit and all contributions go in the force vector vBAS[]
 */
void Meca::addTorqueExplicit(const Interpolation & ptA,
                             const Interpolation & ptB,
                             const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    Vector da = ptA.diff();
    Vector db = ptB.diff();

    real na = da.normSqr();
    real nb = db.normSqr();

    Torque Tq = cross(da, db);
    
#if ( DIM >= 3 )
    const real Tn = Tq.norm();
#else
    const real Tn = fabs(Tq);
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = atan2(Tn, dot(da, db));
    
    /**
     Scale torque to make it proportional to angle:
     we multiply the vector Tq by angle / sin(angle),
     knowing that Tn = Tq.norm() = sin(angle) * sqrt( na * nb )
     
     To have a Torque proportional to sin(angle), use:
     real nn = sqrt( na * nb ); or nn = Tn / angle
     na = weight / ( na * nn );
     nb = weight / ( nb * nn );
     */
    na = weight * angle / ( na * Tn );
    nb = weight * angle / ( nb * Tn );
    
    Vector fa = cross(Tq * na, da);
    Vector fb = cross(db, Tq * nb);

    sub_base(ii0, fa);
    add_base(ii1, fa);
    sub_base(ii2, fb);
    add_base(ii3, fb);
}


/**
 Add an explicit torque to induce two segments to make an angle
 defined by (cosinus, sinus) relative to each other:
 
     torque = weigth * cross( dirA , dirB.rotated(angle) )
     forceA =  cross(torque, dirA)
     forceB = -cross(torque, dirB)
 
 The direction of `ptB` is rotated around `axis` defined as cross( dirA, dirB ).
 The calculation is explicit and all contributions go in the force vector vBAS[]
 It is assumed that `cosinus^2 + sinus^2 = 1`
 Note that if ( sinus == 0 ), you can use addTorque(ptA, ptB, weight)
 */
void Meca::addTorqueExplicit(const Interpolation & ptA,
                             const Interpolation & ptB,
                             const real cosinus, const real sinus,
                             const real weight)
{
    assert_true( weight >= 0 );

    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    Vector da = ptA.diff();
    Vector db = ptB.diff();

    real na = da.normSqr();
    real nb = db.normSqr();
    
#if ( DIM >= 3 )
    
    /*
     in 3D the axis of torque is perpendicular to both `da` and `db`,
     and the angle is only defined between 0 and PI,
     */
    Vector axis = cross(db, da).normalized(std::copysign(1.0, sinus));
    
    // rotate vector `db` around `arm` by angle specified as (cosinus, sinus):
    Vector rot = cosinus * db + sinus * cross(axis, db);
    
#elif ( DIM == 2 )

    // this correspond to the Z-direction, up or down:
    real dir = std::copysign(1, cross(da, db));

    // rotate vector `db` by angle defined by (cosinus, sinus) around Z
    Vector rot( db.XX*cosinus + db.YY*sinus*dir, db.YY*cosinus - db.XX*sinus*dir );
    
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
    const real Tn = fabs(Tq);
#endif
    
    // calculate current angle between segments in [0, pi]:
    const real angle = atan2(Tn, dot(da, rot));
    
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
    sub_base(ii0, fa);
    add_base(ii1, fa);
    sub_base(ii2, fb);
    add_base(ii3, fb);
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
void Meca::addTorquePoliti(const Interpolation & pt1,
                           const Interpolation & pt2,
                           const real cosinus, const real sinus,
                           const real weight)
{
    assert_true( weight >= 0 );
    if ( pt1.overlapping(pt2) )
        return;
    
    //index in the matrix mC:
    const index_t index[] = { DIM*pt1.matIndex1(), DIM*pt1.matIndex1()+1,
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
    ce.XX =  cd.XX*cosinus + cd.YY*sinus;
    ce.YY = -cd.XX*sinus   + cd.YY*cosinus;
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
    //real dangle = atan2( cross(ab, ce), dot(ab, ce) );
    real dangle = atan2( ab.XX*ce.YY - ab.YY*ce.XX, dot(ab, ce) );
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
            m[ii + shifta] += pow(-1,jj)*entrya[ii-jj];
            m[ii + shiftc] += pow(-1,jj)*entryc[ii-jj];
        }
        shifta += 7 - jj;
        shiftc += 3 - jj;
    }
    
    
    //very Cumbersome!!!
    //Entries for Matrix mC  and vector vBAS
    //vBAS = fa - M*P0
    for ( int ii = 0; ii < 8; ++ii )
    {
        vBAS[index[ii]] += w2*h[ii];
        for (int jj = 0; jj < 8; ++jj) {
            if (jj < ii)
                mm = int(jj*(7.5-0.5*jj)+ii);
            else {
                mm = int(ii*(7.5-0.5*ii)+jj);
                mC( index[ii], index[jj] ) += m[mm];
            }
            vBAS[index[ii]] -= m[mm]*coord[jj];
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
 The parameter `equilibrium_angle` is the resting angle specified as ( cosinus, sinus )
 The interpolation coefficients of `pt1` and `pt2` are ignored.

 3D implicit torque implementation
 Serge Dmitrieff and FJN, 21.01.2019 -- 11.02.2019
 www.biophysics.fr and Cambridge University
*/
#if ( DIM > 1 )
void Meca::addTorque(const Interpolation & pt1,
                     const Interpolation & pt2,
                     const real cosinus, const real sinus,
                     const real weight)
{
    const Vector AB = pt1.diff();
    const Vector CD = pt2.diff();
    const real iU = AB.inv_norm();
    const real iV = CD.inv_norm();
    const real wU = weight * iU;
    const real wV = weight * iV;
    const Vector u = AB * iU;
    const Vector v = CD * iV;

#if ( DIM == 3 )
    Vector axis = normalize(cross(AB,CD));
    const Matrix33 R = Matrix33::rotationAroundAxis(axis, cosinus, sinus);
    //const Matrix33 T = R.transposed();
    const Matrix33 Id(1,0,0,0,1,0,0,0,1);  // identity matrix
#else
    const Matrix22 R(cosinus, sinus, -sinus, cosinus);
    //const Matrix22 T(cosinus, -sinus, sinus, cosinus);
    const Matrix22 Id(1,0,0,1);  // identity matrix
#endif

    Vector Ru = R.vecmul(u);
    //Vector Tv = T.vecmul(v);
    Vector Tv = R.trans_vecmul(v);

    // indices in matrix mC:
    const index_t iiA = DIM * pt1.matIndex1();
    const index_t iiB = DIM * pt1.matIndex2();
    const index_t iiC = DIM * pt2.matIndex1();
    const index_t iiD = DIM * pt2.matIndex2();
    
#if ( 1 )
    real Tvu = dot(Tv, u);
    real Ruv = dot(Ru, v);
    assert_small(Tvu-Ruv); // Tvu and Ruv should be equal
    // exact formula
    Vector Fu = ( Tv - Tvu * u ) * wU;
    Vector Fv = ( Ru - Ruv * v ) * wV;
#else
    // approximate formula, if close to equilibrium: Tvu = Ruv = 1
    Vector Fu = ( Tv - u ) * wU;
    Vector Fv = ( Ru - v ) * wV;
#endif

    //std::clog << std::fixed;
    //std::clog << std::setw(9) << Tvu << " " << std::setw(9) << Ruv << std::setw(9) << sinus << "\n";
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

    add_diag_block(iiA, duFu);
    sub_block(iiB, iiA, duFu);
    add_diag_block(iiB, duFu);
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
    add_diag_block(iiC, dvFv);
    sub_block(iiD, iiC, dvFv);
    add_diag_block(iiD, dvFv);
    
    // remaining part of the force (these vectors should be small near equilibrium)
    Vector Fu0 = Fu - duFu.vecmul(AB) - duFv.trans_vecmul(CD);
    Vector Fv0 = Fv - duFv.vecmul(AB) - dvFv.vecmul(CD);

    // add constant terms into vBAS[]:
    sub_base(iiA, Fu0);        // F(A) = -Fu
    add_base(iiB, Fu0);        // F(B) = +Fu
    sub_base(iiC, Fv0);        // F(C) = -Fv
    add_base(iiD, Fv0);        // F(D) = +Fv

#if ( 0 )
    const MatrixBlock dvFu = duFv.transposed();
    //std::clog << iiA << " " << iiB << " " << iiC << " " << iiD << " " << sinus << "\n";
    std::clog << "duFu " << std::setw(12) << duFu << " duFv " << std::setw(12) << duFv << "\n";
    std::clog << "dvFu " << std::setw(12) << dvFu << " dvFv " << std::setw(12) << dvFv << "\n";
    std::clog << "Fu0 " << std::fixed << std::setw(12) << Fu0 << "    Fv0 " << std::setw(12) << Fv0 << "\n";
#endif
}
#endif


/**
 Add Torque between 3 points.
 This version does not impose any particular distance between the points,
 and just move them to enforce the angle described by ABC
 */
void Meca::addTorque(const Mecapoint & ptA,
                     const Mecapoint & ptB,
                     const Mecapoint & ptC,
                     const real cosinus, const real sinus,
                     const real weight)
{
    assert_true( weight >= 0 );
    const MatrixBlock W(0, weight);

#if ( DIM == 3 )
    const Vector AB = ptB.pos() - ptA.pos();
    const Vector BC = ptC.pos() - ptB.pos();
    Vector axis = normalize(cross(AB,BC));
    const Matrix33 R = Matrix33::rotationAroundAxis(axis, cosinus, sinus) * weight;
#elif ( DIM == 2 )
    const Matrix22 R = weight * Matrix22(cosinus, sinus,-sinus, cosinus);
#else
    const Matrix11 R(1);  //should not be used!
#endif

    const MatrixBlock T = R.transposed();

    // indices in matrix mC:
    const index_t iiA = DIM * ptA.matIndex();
    const index_t iiB = DIM * ptB.matIndex();
    const index_t iiC = DIM * ptC.matIndex();
    
    // the term (R+W)+(T+W) is diagonal
    MatrixBlock D(0, 2 * ( weight * cosinus + weight ));
    
    sub_diag_block(iiA, W);
    sub_diag_block(iiB, D); //(R+W)+(T+W)
    sub_diag_block(iiC, W);
    if ( iiB > iiA )
        add_block(iiB, iiA, W+R);
    else
        add_block(iiA, iiB, W+T);
    if ( iiC > iiA )
        sub_block(iiC, iiA, R);
    else
        sub_block(iiA, iiC, T);
    if ( iiC > iiB )
        add_block(iiC, iiB, W+R);
    else
        add_block(iiB, iiC, W+T);
}



/** This is variation 3, 20.08.2019
 It combines addTorque() without length with a LongLink(ptA, ptB);
 */
void Meca::addTorque(const Mecapoint & ptA,
                     const Mecapoint & ptB,
                     const Mecapoint & ptC,
                     const real cosinus, const real sinus,
                     const real len, const real weight)
{
    assert_true( weight >= 0 );
    const MatrixBlock W(0, weight);
    const Vector AB = ptB.pos() - ptA.pos();
    
#if ( DIM == 3 )
    const Vector BC = ptC.pos() - ptB.pos();
    Vector axis = normalize(cross(AB,BC));
    const Matrix33 R = Matrix33::rotationAroundAxis(axis, cosinus, sinus) * weight;
#elif ( DIM == 2 )
    const Matrix22 R = weight * Matrix22(cosinus, sinus,-sinus, cosinus);
#else
    const Matrix11 R(1);  //should not be used!
#endif
    
    const MatrixBlock T = R.transposed();
    
    // indices in matrix mC:
    const index_t iiA = DIM * ptA.matIndex();
    const index_t iiB = DIM * ptB.matIndex();
    const index_t iiC = DIM * ptC.matIndex();
    
    // this is a LongLink(A, B):
    MatrixBlock P(0,0);
    const real ab2 = AB.normSqr();
    if ( ab2 > REAL_EPSILON )
    {
        real ab = sqrt(ab2);
        const real wla = weight * len / ab;
        add_base(iiA, AB, -wla);
        add_base(iiB, AB,  wla);

        // regularize interaction to avoid creating negative eigen values
        if ( ab < len )
            P = MatrixBlock::outerProduct(AB, weight/ab2);
        else
            P = MatrixBlock::offsetOuterProduct(weight-wla, AB, wla/ab2);
    }    
    
    // the term (R+W)+(T+W) is diagonal
    MatrixBlock D(0, 2 * ( weight * cosinus + weight ));
    
    sub_diag_block(iiA, W+P);
    sub_diag_block(iiB, D+P); //(R+W)+(T+W)+P
    sub_diag_block(iiC, W);
    if ( iiB > iiA )
        add_block(iiB, iiA, W+R+P);
    else
        add_block(iiA, iiB, W+T+P);
    if ( iiC > iiA )
        sub_block(iiC, iiA, R);
    else
        sub_block(iiA, iiC, T);
    if ( iiC > iiB )
        add_block(iiC, iiB, W+R);
    else
        add_block(iiB, iiC, W+T);
    
    if ( modulo )
        throw Exception("addTorque(A,B,C) is not usable with periodic boundary conditions");
}


//------------------------------------------------------------------------------
#pragma mark - Interpolation over position vector
//------------------------------------------------------------------------------


Vector Meca::position2(const index_t inx[2], const real coef[2]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    return coef[0] * P0 + coef[1] * P1;
}


Vector Meca::position3(const index_t inx[3], const real coef[3]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    return ( coef[0] * P0 + coef[1] * P1 ) + coef[2] * P2;
}


Vector Meca::position4(const index_t inx[4], const real coef[4]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    Vector P3(vPTS+inx[3]);
    return ( coef[0] * P0 + coef[1] * P1 ) + ( coef[2] * P2 + coef[3] * P3 );
}


Vector Meca::position5(const index_t inx[5], const real coef[5]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    Vector P3(vPTS+inx[3]);
    Vector P4(vPTS+inx[4]);
    return ( coef[0] * P0 + coef[1] * P1 ) + ( coef[2] * P2 + coef[3] * P3 ) + coef[4] * P4;
}


Vector Meca::position6(const index_t inx[6], const real coef[6]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    Vector P3(vPTS+inx[3]);
    Vector P4(vPTS+inx[4]);
    Vector P5(vPTS+inx[5]);
    return ( coef[0] * P0 + coef[1] * P1 ) + ( coef[2] * P2 + coef[3] * P3 ) + ( coef[4] * P4 + coef[5] * P5 );
}

//------------------------------------------------------------------------------
#pragma mark - Links between Mecables
//------------------------------------------------------------------------------

/**
 Link `ptA` (A) and `ptB` (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 In practice, Meca::addLink() will update the matrix mB,
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
void Meca::addLink(const Mecapoint & ptA,
                   const Mecapoint & ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    const index_t ii0 = ptA.matIndex();
    const index_t ii1 = ptB.matIndex();

    if ( ii0 == ii1 )
        return;

    sub_iso(ii0, ii0, weight);
    sub_iso(ii1, ii1, weight);
    add_iso(ii0, ii1, weight);

    if ( modulo )
    {
        const real ww[] = { weight, -weight };
#if ( 1 )
        Vector off = modulo->offset(ptA.pos() - ptB.pos());
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
        }
#else
        const index_t inx[] = { DIM*ii0, DIM*ii1 };
        Vector off = modulo->offset(position2(inx, ww));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off);
            sub_base(DIM*ii1, off);
        }
#endif
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), ptB.pos());
    }
#endif
}


/**
 Link interpolated vertex (A) and vertex (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 */
void Meca::addLink(const Interpolation & ptA,
                   const Mecapoint & ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;

    //coefficients on the points:
    const real cc[] = { ptA.coef2(),   ptA.coef1(),    -1.0 };
    const real ww[] = { weight*cc[0], weight*cc[1], -weight };
    
    sub_iso(ii0, ii0, ww[0] * cc[0]);
    sub_iso(ii1, ii0, ww[1] * cc[0]);
    add_iso(ii2, ii0, ww[0]);         // since cc[2] == -1
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    add_iso(ii2, ii1, ww[1]);         // since cc[2] == -1
    add_iso(ii2, ii2, ww[2]);         // since cc[2] == -1
 
    if ( modulo )
    {
#if ( 1 )
        Vector off = modulo->offset(ptA.pos() - ptB.pos());
#else
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
#endif
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), ptB.pos());
    }
#endif
}


/**
 Link vertex (A) and interpolated vertex (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )

 */
void Meca::addLink(const Mecapoint & ptA,
                   const Interpolation & ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = ptA.matIndex();
    const index_t ii1 = ptB.matIndex1();
    const index_t ii2 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    //coefficients on the points:
    const real cc[] = {    1.0, -ptB.coef2(), -ptB.coef1() };
    const real ww[] = { weight, weight*cc[1], weight*cc[2] };
    
    sub_iso(ii0, ii0, ww[0]);  // since cc[0]==1
    sub_iso(ii1, ii0, ww[1]);  // since cc[0]==1
    sub_iso(ii2, ii0, ww[2]);  // since cc[0]==1
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
  
    if ( modulo )
    {
#if ( 1 )
        Vector off = modulo->offset(ptA.pos() - ptB.pos());
#else
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
#endif
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptB.mecable()->signature()).load();
        gle::drawLink(ptB.pos(), ptA.pos());
    }
#endif
}


/**
 Link `ptA` (A) and `ptB` (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )

 */
void Meca::addLink(const Interpolation & ptA,
                   const Interpolation & ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = ptA.matIndex1();
    const index_t ii1 = ptA.matIndex2();
    const index_t ii2 = ptB.matIndex1();
    const index_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    //interpolation coefficients:
    const real cc[] = {  ptA.coef2(),  ptA.coef1(), -ptB.coef2(), -ptB.coef1() };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*cc[2], weight*cc[3] };
    
    sub_iso(ii0, ii0, ww[0] * cc[0]);
    sub_iso(ii1, ii0, ww[1] * cc[0]);
    sub_iso(ii2, ii0, ww[2] * cc[0]);
    sub_iso(ii3, ii0, ww[3] * cc[0]);
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    sub_iso(ii3, ii1, ww[3] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    sub_iso(ii3, ii2, ww[3] * cc[2]);
    
    sub_iso(ii3, ii3, ww[3] * cc[3]);
    
    if ( modulo )
    {
#if ( 1 )
        Vector off = modulo->offset(ptA.pos() - ptB.pos());
#else
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3 };
        Vector off = modulo->offset(position4(inx, cc));
#endif
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
            add_base(DIM*ii3, off, ww[3]);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), ptB.pos());
    }
#endif
    
}

//------------------------------------------------------------------------------
#pragma mark - Links between Mecable (higher order interpolation)
//------------------------------------------------------------------------------

/**
 Link interpolation (A) and vertex at index 'pts' (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in the vertex of a Mecable, at index 'pts'.
 Diagonal and lower elements of mB are set.
 */
void Meca::addLink1(const Interpolation & pti,
                    const index_t pts,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = pts;
    const index_t ii1 = pti.matIndex1();
    const index_t ii2 = pti.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    const real cc[] = {    1.0, -pti.coef2(), -pti.coef1() };
    const real ww[] = { weight, weight*cc[1], weight*cc[2] };
    
    sub_iso(ii0, ii0, ww[0]); // since cc[0] == 1.0
    sub_iso(ii1, ii0, ww[1]); // since cc[0] == 1.0
    sub_iso(ii2, ii0, ww[2]); // since cc[0] == 1.0
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    
    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
        }
    }
}


/**
 Link vertex (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is a Mecapoint (ptA).
 Point B in interpolated over 2 vertices, specified by index in 'pts[]'
 using the coefficients given in `coef[]`.
 Diagonal and lower elements of mB are set.
 */
void Meca::addLink2(const Mecapoint & ptA,
                    const index_t pts[2], const real coef[2],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = ptA.matIndex();
    const index_t ii1 = pts[0];
    const index_t ii2 = pts[1];
    
    if ( any_equal(ii0, ii1, ii2) )
        return;

    const real cc[] = {   -1.0,      coef[0],      coef[1] };
    const real ww[] = {-weight, weight*cc[1], weight*cc[2] };
    
    assert_small(coef[0]+coef[1]-1.0);
    
    add_iso(ii0, ii0, ww[0]); // since cc[0] == -1.0
    add_iso(ii1, ii0, ww[1]); // since cc[0] == -1.0
    add_iso(ii2, ii0, ww[2]); // since cc[0] == -1.0
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    
    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
        }
    }
}

/**
 Link `pti` (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in interpolated over 2 vertices, specified by index in 'pts[]'
 using the coefficients given in `coef[]`.
 Diagonal and lower elements of mB are set.
 */
void Meca::addLink2(const Interpolation & pti,
                    const index_t pts[2], const real coef[2],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    const index_t ii2 = pts[0];
    const index_t ii3 = pts[1];
    
    const real cc[] = { -pti.coef2(), -pti.coef1(),      coef[0],      coef[1] };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*cc[2], weight*cc[3] };

    assert_small(coef[0]+coef[1]-1.0);
    
    sub_iso(ii0, ii0, ww[0] * cc[0]);
    sub_iso(ii1, ii0, ww[1] * cc[0]);
    sub_iso(ii2, ii0, ww[2] * cc[0]);
    sub_iso(ii3, ii0, ww[3] * cc[0]);

    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    sub_iso(ii3, ii1, ww[3] * cc[1]);

    sub_iso(ii2, ii2, ww[2] * cc[2]);
    sub_iso(ii3, ii2, ww[3] * cc[2]);

    sub_iso(ii3, ii3, ww[3] * cc[3]);
    
    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3 };
        Vector off = modulo->offset(position4(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
            add_base(DIM*ii3, off, ww[3]);
        }
    }
}


/**
 Link vertex (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is a Mecapoint.
 Point B in interpolated over 3 vertices, specified by index in 'pts[]'
 using the coefficients given in `coef[]`.
 Diagonal and lower elements of mB are set.
 */
void Meca::addLink3(const Mecapoint & ptA,
                    const index_t pts[3], const real coef[3],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = ptA.matIndex();
    const index_t ii1 = pts[0];
    const index_t ii2 = pts[1];
    const index_t ii3 = pts[2];

    const real cc[] = {   -1.0,      coef[0],      coef[1],      coef[2] };
    const real ww[] = {-weight, weight*cc[1], weight*cc[2], weight*cc[3] };

    assert_small(coef[0]+coef[1]+coef[2]-1.0);
    
    add_iso(ii0, ii0, ww[0]);  // since cc[0] = -1
    add_iso(ii1, ii0, ww[1]);
    add_iso(ii2, ii0, ww[2]);
    add_iso(ii3, ii0, ww[3]);
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    sub_iso(ii3, ii1, ww[3] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    sub_iso(ii3, ii2, ww[3] * cc[2]);
    
    sub_iso(ii3, ii3, ww[3] * cc[3]);
    
    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3 };
        Vector off = modulo->offset(position4(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
            add_base(DIM*ii3, off, ww[3]);
        }
    }
}


/**
 Link `pti` (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in interpolated over 4 vertices, specified by index in 'pts[]'
 using the coefficients given in `coef[]`.
 Diagonal and lower elements of mB are set.
*/
void Meca::addLink3(const Interpolation & pti,
                    const index_t pts[3], const real coef[3],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    const index_t ii2 = pts[0];
    const index_t ii3 = pts[1];
    const index_t ii4 = pts[2];

    const real cc[] = { -pti.coef2(), -pti.coef1(),      coef[0],      coef[1],      coef[2] };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*cc[2], weight*cc[3], weight*cc[4] };
    
    assert_small(coef[0]+coef[1]+coef[2]-1.0);
    
    sub_iso(ii0, ii0, ww[0] * cc[0]);
    sub_iso(ii1, ii0, ww[1] * cc[0]);
    sub_iso(ii2, ii0, ww[2] * cc[0]);
    sub_iso(ii3, ii0, ww[3] * cc[0]);
    sub_iso(ii4, ii0, ww[4] * cc[0]);
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    sub_iso(ii3, ii1, ww[3] * cc[1]);
    sub_iso(ii4, ii1, ww[4] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    sub_iso(ii3, ii2, ww[3] * cc[2]);
    sub_iso(ii4, ii2, ww[4] * cc[2]);
    
    sub_iso(ii3, ii3, ww[3] * cc[3]);
    sub_iso(ii4, ii3, ww[4] * cc[3]);
    
    sub_iso(ii4, ii4, ww[4] * cc[4]);
    
    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3, DIM*ii4 };
        Vector off = modulo->offset(position5(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
            add_base(DIM*ii3, off, ww[3]);
            add_base(DIM*ii4, off, ww[4]);
        }
    }
}

/**
 Link vertex (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is a Mecapoint.
 Point B in interpolated over 4 vertices, specified by index in 'pts[]'
 using the coefficients given in `coef[]`.
 Diagonal and lower elements of mB are set.
 */
void Meca::addLink4(const Mecapoint & ptA,
                    const index_t pts[4], const real coef[4],
                    const real weight)
{
    assert_true( weight >= 0 );

    //index in the matrix mB:
    const index_t ii0 = ptA.matIndex();
    const index_t ii1 = pts[0];
    const index_t ii2 = pts[1];
    const index_t ii3 = pts[2];
    const index_t ii4 = pts[3];

    const real cc[] = {   -1.0,      coef[0],      coef[1],      coef[2],      coef[3] };
    const real ww[] = {-weight, weight*cc[1], weight*cc[2], weight*cc[3], weight*cc[4] };
    
    assert_small(coef[0]+coef[1]+coef[2]+coef[3]-1.0);
    
    add_iso(ii0, ii0, ww[0]);
    add_iso(ii1, ii0, ww[1]);
    add_iso(ii2, ii0, ww[2]);
    add_iso(ii3, ii0, ww[3]);
    add_iso(ii4, ii0, ww[4]);
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    sub_iso(ii3, ii1, ww[3] * cc[1]);
    sub_iso(ii4, ii1, ww[4] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    sub_iso(ii3, ii2, ww[3] * cc[2]);
    sub_iso(ii4, ii2, ww[4] * cc[2]);
    
    sub_iso(ii3, ii3, ww[3] * cc[3]);
    sub_iso(ii4, ii3, ww[4] * cc[3]);
    
    sub_iso(ii4, ii4, ww[4] * cc[4]);
    
    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3, DIM*ii4 };
        Vector off = modulo->offset(position5(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
            add_base(DIM*ii3, off, ww[3]);
            add_base(DIM*ii4, off, ww[4]);
        }
    }
}


/**
 Link `pti` (A) and interpolated point (B)
 The force is linear with a zero resting length:
 
     force_A = weight * ( B - A )
     force_B = weight * ( A - B )
 
 Point A is an Interpolation.
 Point B in interpolated over 4 vertices, specified by index in 'pts[]'
 using the coefficients given in `coef[]`.
 Diagonal and lower elements of mB are set.
*/
void Meca::addLink4(const Interpolation & pti,
                    const index_t pts[4], const real coef[4],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const index_t ii0 = pti.matIndex1();
    const index_t ii1 = pti.matIndex2();
    const index_t ii2 = pts[0];
    const index_t ii3 = pts[1];
    const index_t ii4 = pts[2];
    const index_t ii5 = pts[3];

    const real cc[] = { -pti.coef2(), -pti.coef1(),      coef[0],      coef[1],      coef[2],      coef[3] };
    const real ww[] = { weight*cc[0], weight*cc[1], weight*cc[2], weight*cc[3], weight*cc[4], weight*cc[5] };
    
    assert_small(coef[0]+coef[1]+coef[2]+coef[3]-1.0);
    
    sub_iso(ii0, ii0, ww[0] * cc[0]);
    sub_iso(ii1, ii0, ww[1] * cc[0]);
    sub_iso(ii2, ii0, ww[2] * cc[0]);
    sub_iso(ii3, ii0, ww[3] * cc[0]);
    sub_iso(ii4, ii0, ww[4] * cc[0]);
    sub_iso(ii5, ii0, ww[5] * cc[0]);
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    sub_iso(ii3, ii1, ww[3] * cc[1]);
    sub_iso(ii4, ii1, ww[4] * cc[1]);
    sub_iso(ii5, ii1, ww[5] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    sub_iso(ii3, ii2, ww[3] * cc[2]);
    sub_iso(ii4, ii2, ww[4] * cc[2]);
    sub_iso(ii5, ii2, ww[5] * cc[2]);
    
    sub_iso(ii3, ii3, ww[3] * cc[3]);
    sub_iso(ii4, ii3, ww[4] * cc[3]);
    sub_iso(ii5, ii3, ww[5] * cc[3]);
    
    sub_iso(ii4, ii4, ww[4] * cc[4]);
    sub_iso(ii5, ii4, ww[5] * cc[4]);

    sub_iso(ii5, ii5, ww[5] * cc[5]);

    if ( modulo )
    {
        const index_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3, DIM*ii4, DIM*ii5 };
        Vector off = modulo->offset(position6(inx, cc));
        if ( !off.null() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
            add_base(DIM*ii2, off, ww[2]);
            add_base(DIM*ii3, off, ww[3]);
            add_base(DIM*ii4, off, ww[4]);
            add_base(DIM*ii5, off, ww[5]);
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
 
 */

void Meca::addLongLink(const Mecapoint & ptA,
                       const Mecapoint & ptB,
                       const real len,
                       const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    const index_t ia = DIM * ptA.matIndex();  // coef is +weight
    const index_t ib = DIM * ptB.matIndex();  // coef is -weight

    if ( ia == ib )
        return;
    
    Vector off, axi = ptB.pos() - ptA.pos();

    if ( modulo )
        modulo->foldOffset(axi, off);
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), axi, len);
    }
#endif

    const real abn = axi.norm();
    if ( abn < REAL_EPSILON )
        return;
    
    const real wla = weight * len / abn;

    add_base(ia, axi,-wla);
    add_base(ib, axi, wla);
    
    axi /= abn;

    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    const bool cooked = ( len > abn );
    
    MatrixBlock T;

    if ( cooked )
        T = MatrixBlock::outerProduct(axi, weight);
    else
        T = MatrixBlock::offsetOuterProduct(weight-wla, axi, wla);
    
    sub_diag_block(ia, T);
    sub_diag_block(ib, T);
    if ( ia > ib )
        add_block(ia, ib, T);
    else
        add_block(ib, ia, T);
    
    if ( modulo && !off.null() )
    {
        if ( cooked )
            off = ( weight * dot(off, axi) ) * axi;
        else
            off = ( wla * dot(off, axi) ) * axi + ( weight - wla ) * off;
        sub_base(ia, off);
        add_base(ib, off);
    }
}


/**
 Link vertex (A) and interpolation (B),
 The force is affine with non-zero resting length: 
 
     force_A = weight * ( B - A ) * ( len / |AB| - 1 )
     force_B = weight * ( A - B ) * ( len / |AB| - 1 )
 
 */

void Meca::addLongLink(const Mecapoint & ptA,
                       const Interpolation & pti,
                       const real len, const real weight )
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex();
    const index_t ii1 = DIM * pti.matIndex1();
    const index_t ii2 = DIM * pti.matIndex2();

    if ( any_equal(ii0, ii1, ii2) )
        return;

    //force coefficients on the points:
    const real cc[] = {   -1.0,   pti.coef2(),   pti.coef1() };
    const real ww[] = { weight, -weight*cc[1], -weight*cc[2] };

    Vector off, axi = ptA.pos() - pti.pos();

    if ( modulo )
        modulo->foldOffset(axi, off);
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(pti.mecable()->signature()).load();
        gle::drawLink(pti.pos(), axi, len);
    }
#endif
    
    const real abn = axi.norm();
    if ( abn < REAL_EPSILON ) return;
    
    axi /= abn;

    add_base(ii0, axi, ww[0]);
    add_base(ii1, axi, ww[1]);
    add_base(ii2, axi, ww[2]);

    real lab = len / abn;

    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    const bool cooked = ( len > abn );
    
    MatrixBlock T;

    if ( cooked )
        T = MatrixBlock::outerProduct(axi);
    else
        T = MatrixBlock::offsetOuterProduct(1.0-lab, axi, lab);

    Matrix33 W = Matrix33::outerProduct(cc, ww);
    
    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_diag_block(ii2, W(2,2), T);

    if ( modulo && !off.null() )
    {
        if ( cooked )
            off = dot(off, axi) * axi;
        else
            off = ( lab * dot(off, axi) ) * axi + ( 1.0 - lab ) * off;
        
        add_base(ii0, off, ww[0]);
        add_base(ii1, off, ww[1]);
        add_base(ii2, off, ww[2]);
    }
}


/**
 Link `ptA` (A) and `ptB` (B),
 The force is affine with non-zero resting length: 
 
     force_A = weight * ( B - A ) * ( len / |AB| - 1 )
     force_B = weight * ( A - B ) * ( len / |AB| - 1 )

 */

void Meca::addLongLink(const Interpolation & ptA,
                       const Interpolation & ptB,
                       const real len,
                       const real weight )
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();

    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    //force coefficients on the points:
    const real cc[] = { ptA.coef2(),   ptA.coef1(), -ptB.coef2(), -ptB.coef1() };
    const real ww[] = {-weight*cc[0],-weight*cc[1],-weight*cc[2],-weight*cc[3] };
    
    Vector off, axi = ptB.pos() - ptA.pos();

    if ( modulo )
        modulo->foldOffset(axi, off);
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), axi, len);
    }
#endif

    const real abn = axi.norm();
    if ( abn < REAL_EPSILON ) return;
    
    axi /= abn;
    
    add_base(ii0, axi, len*ww[0]);
    add_base(ii1, axi, len*ww[1]);
    add_base(ii2, axi, len*ww[2]);
    add_base(ii3, axi, len*ww[3]);

    real lab = len / abn;
    
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    const bool cooked = ( len > abn );
    
    MatrixBlock T;
    
    if ( cooked )
        T = MatrixBlock::outerProduct(axi);
    else
        T = MatrixBlock::offsetOuterProduct(1.0-lab, axi, lab);
    
    Matrix44 W = Matrix44::outerProduct(cc, ww);
    
    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_block(ii3, ii0, W(3,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_block(ii3, ii1, W(3,1), T);
    add_diag_block(ii2, W(2,2), T);
    add_block(ii3, ii2, W(3,2), T);
    add_diag_block(ii3, W(3,3), T);

    if ( modulo && !off.null() )
    {
        if ( cooked )
            off = dot(off, axi) * axi;
        else
            off = ( lab * dot(off, axi) ) * axi + ( 1.0 - lab ) * off;

        add_base(ii0, off, ww[0]);
        add_base(ii1, off, ww[1]);
        add_base(ii2, off, ww[2]);
        add_base(ii3, off, ww[3]);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Off-axis links between Mecable
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
 
 @todo interSideLink2D should use block operations
 */


#if ( DIM == 2 )

void Meca::addSideLink2D(const Interpolation & ptA,
                         const Mecapoint & ptB,
                         const real arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    index_t ii0 = ptA.matIndex1();
    index_t ii1 = ptA.matIndex2();
    index_t ii2 = ptB.matIndex();

    if ( any_equal(ii0, ii1, ii2) )
        return;

    //force coefficients on the points:
    const real ca1 = ptA.coef2();
    const real ca2 = ptA.coef1();
    const real eps = arm / ptA.len();
    
    const real ca1w = weight * ca1;
    const real ca2w = weight * ca2;
    const real epsw = weight * eps;
    const real e2sw = eps * epsw;
    
    //\todo put all terms in blocks (mC)
    //we put the isotropic terms in mB
    sub_iso(ii0, ii0, ca1w * ca1 + e2sw);
    sub_iso(ii0, ii1, ca1w * ca2 - e2sw);
    sub_iso(ii1, ii1, ca2w * ca2 + e2sw);
    
    sub_iso(ii2, ii2, weight);
    add_iso(ii0, ii2, ca1w);
    add_iso(ii1, ii2, ca2w);
    
    //index in the matrix mC:
    ii0 *= DIM;
    ii1 *= DIM;
    ii2 *= DIM;
    
    mC(ii0  , ii1+1) += epsw;
    mC(ii0+1, ii1  ) -= epsw;
    
    mC(ii0  , ii2+1) -= epsw;
    mC(ii0+1, ii2  ) += epsw;
    mC(ii1,   ii2+1) += epsw;
    mC(ii1+1, ii2  ) -= epsw;
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        
        real offx = off.XX;
        if ( offx != 0 )
        {
            vBAS[ii0  ] += ca1w * offx;
            vBAS[ii0+1] += epsw * offx;
            vBAS[ii1  ] += ca2w * offx;
            vBAS[ii1+1] -= epsw * offx;
            vBAS[ii2  ] += offx;
        }
        real offy = off.YY;
        if ( offy != 0 )
        {
            vBAS[ii0  ] -= epsw * offy;
            vBAS[ii0+1] += ca1w * offy;
            vBAS[ii1  ] += epsw * offy;
            vBAS[ii1+1] += ca2w * offy;
            vBAS[ii2+1] += offy;
        }
    }
}

#elif ( DIM >= 3 )


/**
 This is experimental and should not be used.
 
 Link `ptA` (A) and `ptB` (B)
 and a point S which is on the side of `ptA`.
 
     S = pos_a + cross( arm, dir_a)

 arm must be perpendicular to link
 
 @todo interSideLink3D should use block operations
 */
void Meca::addSideLink3D(const Interpolation & ptA,
                         const Mecapoint & ptB,
                         Vector const& arm,
                         const real weight)
{
    assert_true( weight >= 0 );

    // indices to mC:
    const index_t ia1 = ptA.matIndex1();
    const index_t ia2 = ptA.matIndex2();
    const index_t ib  = ptB.matIndex();
    
    if ( any_equal(ia1, ia2, ib) )
        return;
    
    const index_t inx[6] = { DIM*ia1, DIM*ia1+1, DIM*ia1+2, DIM*ia2, DIM*ia2+1, DIM*ia2+2 };

    real a = ptA.coef2();
    real b = ptA.coef1();
    real s = 1.0 / ptA.len();
    
    real ex = s * arm.XX;
    real ey = s * arm.YY;
    real ez = s * arm.ZZ;
    
    /* The transfer matrix transforms the two Mecapoint in ptA,
     to the side point S:
     S = aa * pt1 + bb * pt2 + cross(arm, normalize( pt2 - pt1 ))
     
     It was generated in Maxima:
     MVP: matrix([0, -ez, ey], [ez, 0, -ex], [-ey, ex, 0]);
     MD: addcol(-ident(3), ident(3));
     MC: addcol(aa*ident(3), bb*ident(3));
     T: MC+MVP.MD;
     */
    const real T[18] = {
        a,  ez, -ey,   b, -ez,  ey,
      -ez,   a,  ex,  ez,   b, -ex,
       ey, -ex,   a, -ey,  ex,   b
    };
    
    real a2 = a * a, ab = a * b;
    real b2 = b * b;
    
    real exx = ex * ex, exy = ex*ey, exz = ex*ez;
    real eyy = ey * ey, eyz = ey*ez;
    real ezz = ez * ez;
    
    // TT = transpose(T) * T is symmetric, and thus we only set half of it:
    /* Maxima code:
     TT: expand(transpose(T) . T);
     */
    real TT[36] = {
        eyy+ezz+a2,  0,           0,           0,           0,           0,
        -exy,        exx+ezz+a2,  0,           0,           0,           0,
        -exz,       -eyz,         exx+eyy+a2,  0,           0,           0,
        -ezz-eyy+ab, ez+exy,      exz-ey,      eyy+ezz+b2,  0,           0,
        -ez+exy,    -ezz-exx+ab,  eyz+ex,     -exy,         exx+ezz+b2,  0,
        exz+ey,      eyz-ex,     -eyy-exx+ab, -exz,        -eyz,         exx+eyy+b2
    };
    
    // we project to bring all forces in the plane perpendicular to 'arm'
    //real sca = arm.inv_norm();
    //real an = a * sca;
    //real bn = b * sca;
    // Maxima code: matrix([ex, ey, ez]) . T;
    //real TP[9] = { an*ex, an*ey, an*ez, bn*ex, bn*ey, bn*ez, -ex, -ey, -ez };
    //blas::xgemm('N','N', 6, 1, 3, sca, T, 6, arm, 3, 0.0, TP, 6);
    
    //blas::xsyrk('U','N', 6, 1, weight, TP, 6, -weight, TT, 6);
    blas::xscal(36, -weight, TT, 1);
    
    for ( int ii=0 ; ii<6; ++ii )
    for ( int jj=ii; jj<6; ++jj )
        mC(inx[ii], inx[jj]) += TT[ii+6*jj];
    
    //sub_iso(ia1, ib, a * weight);
    //sub_iso(ia2, ib, b * weight);
    sub_iso(ib, ib, weight);
    
    for ( int ii=0; ii<6; ++ii )
    for ( int jj=0; jj<3; ++jj )
        mC(inx[ii], DIM*ib+jj) += weight * T[ii+6*jj];
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm,ptA.dir()), ptB.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("addSideLink3D is not usable with periodic boundary conditions");
}


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

void Meca::addSideLinkS(const Interpolation & ptA,
                        const Mecapoint & ptB,
                        Vector const& arm,
                        const real len,
                        const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );

    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex();

    if ( any_equal(ii0, ii1, ii2) )
        return;

    MatrixBlock T;
    // vector 'a' is parallel to first Fiber
    {
        Vector a = ptA.diff();
        Vector v = a / a.normSqr();
        Vector b = arm / len;
        
        // we can set directly the interaction coefficient matrix:
        T = MatrixBlock::outerProduct(a, v) + MatrixBlock::outerProduct(b);
    }
    
    // weights and indices:
    const real cc[3] = {   ptA.coef2(),   ptA.coef1(),   -1.0 };
    const real ww[3] = { -weight*cc[0], -weight*cc[1], weight };

    Matrix33 W = Matrix33::outerProduct(cc, ww);
    
    add_base(ii0, arm, ww[0]);
    add_base(ii1, arm, ww[1]);
    add_base(ii2, arm, ww[2]);

    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_diag_block(ii2, W(2,2), T);

#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), arm, ptB.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("addSideLinkS is not usable with periodic boundary conditions");
}
#endif


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


void Meca::addSideLink(const Interpolation & ptA,
                       const Mecapoint & ptB,
                       const real len,
                       const real weight )
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real arm = len * RNG.sign_exc(cross(ptA.diff(), ptB.pos()-ptA.pos()));
    addSideLink2D(ptA, ptB, arm, weight);

#else
    
    // set 'arm' perpendicular to direction of the Fiber associated with `ptA`:
    Vector arm = calculateArm(ptB.pos()-ptA.pos(), ptA.diff(), len);
    addSideLinkS(ptA, ptB, arm, len, weight);
    
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Off-axis links between Mecable (Interpolation)

#if ( DIM == 2 )

void Meca::addSideLink2D(const Interpolation & ptA,
                         const Interpolation & ptB,
                         const real arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB and mC:
    const index_t ia0 = ptA.matIndex1(),  ii0 = DIM * ia0;
    const index_t ia1 = ptA.matIndex2(),  ii1 = DIM * ia1;
    const index_t ib2 = ptB.matIndex1(),  ii2 = DIM * ib2;
    const index_t ib3 = ptB.matIndex2(),  ii3 = DIM * ib3;
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    // weights and indices:
    const real w = -weight;
    const real cc0 =  ptA.coef2(),  ww0 = w * cc0;
    const real cc1 =  ptA.coef1(),  ww1 = w * cc1;
    const real cc2 = -ptB.coef2(),  ww2 = w * cc2;
    const real cc3 = -ptB.coef1(),  ww3 = w * cc3;
    
    const real ee = arm / ptA.len(),  we = w * ee;

    Matrix22 A(cc0, -ee,  ee, cc0);
    Matrix22 B(cc1,  ee, -ee, cc1);

    add_iso(ia0, ia0, ww0 * cc0 + we * ee);
    add_iso(ia1, ia1, ww1 * cc1 + we * ee);
    add_iso(ib2, ib2, ww2 * cc2);
    add_iso(ib3, ib3, ww3 * cc3);
    add_iso(ib3, ib2, ww3 * cc2);

    add_block(ii1, ii0, w, B.trans_mul(A));
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, ww2, A);
        add_block(ii3, ii0, ww3, A);
        add_block(ii2, ii1, ww2, B);
        add_block(ii3, ii1, ww3, B);
    }
    else
    {
        Matrix22 At = A.transposed();
        Matrix22 Bt = B.transposed();
        add_block(ii0, ii2, ww2, At);
        add_block(ii1, ii2, ww2, Bt);
        add_block(ii0, ii3, ww3, At);
        add_block(ii1, ii3, ww3, Bt);
    }
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off *= -weight;
            add_base(ii0, A.trans_vecmul(off));
            add_base(ii1, B.trans_vecmul(off));
            add_base(ii2, off, cc2);
            add_base(ii3, off, cc3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), cross(arm,ptA.diff()), ptB.pos());
    }
#endif
}

#elif ( DIM >= 3 )

void Meca::addSideLinkS(const Interpolation & ptA,
                        const Interpolation & ptB,
                        Vector const& arm,
                        const real len,
                        const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );

    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();

    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    MatrixBlock T;
    {
        Vector a = ptA.diff();
        Vector v = a / a.normSqr();
        Vector b = arm / len;
        // Vector c = cross(a, b);
        
        // we can set directly the interaction coefficient matrix:
        T = MatrixBlock::outerProduct(a, v) + MatrixBlock::outerProduct(b);
    }
    
    // weights and indices:
    const real cc[4] = {   ptA.coef2(),   ptA.coef1(),  -ptB.coef2(),  -ptB.coef1() };
    const real ww[4] = { -weight*cc[0], -weight*cc[1], -weight*cc[2], -weight*cc[3] };
    
    Matrix44 W = Matrix44::outerProduct(cc, ww);
    
    add_base(ii0, arm, ww[0]);
    add_base(ii1, arm, ww[1]);
    add_base(ii2, arm, ww[2]);
    add_base(ii3, arm, ww[3]);
 
    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_block(ii3, ii0, W(3,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_block(ii3, ii1, W(3,1), T);
    add_diag_block(ii2, W(2,2), T);
    add_block(ii3, ii2, W(3,2), T);
    add_diag_block(ii3, W(3,3), T);

#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), arm, ptB.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("addSideLinkS is not usable with periodic boundary conditions");
}

#endif


/**
 Link `ptA` (A) and `ptB` (B),
 Which is taken between B and a point S located on the side of A:
 S = A + len * N,
 where N is a normalized vector orthogonal to the fiber in an.
 S is linearly related to the two vertices on the sides of A, P1 and P2
 In 3D S is choosen in the plane of P1, P2 and B.
 The force is linear of zero resting length:
 
     force_S = weight * ( S - B )
     force_B = weight * ( B - S )
 
 */

void Meca::addSideLink(const Interpolation & ptA,
                       const Interpolation & ptB,
                       const real len,
                       const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideLink is meaningless in 1D");

#elif ( DIM == 2 )
    
    real arm = len * RNG.sign_exc( cross(ptA.diff(), ptB.pos()-ptA.pos()) );
    addSideLink2D(ptA, ptB, arm, weight);
    
#else

    // set 'arm' perpendicular to direction of the Fiber associated with `ptA`:
    Vector arm = calculateArm(ptB.pos()-ptA.pos(), ptA.diff(), len);
    addSideLinkS(ptA, ptB, arm, len, weight);
    
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Symmetric off-axis links
//------------------------------------------------------------------------------

#if ( DIM == 2 )

/*
void Meca::addSideSideLink2D(const Interpolation & ptA,
                             const Interpolation & ptB,
                             const real len,
                             const real weight,
                             real side1, real side2 )
{
    assert_true( weight >= 0 );
 
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();
 
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    // weights and indices:
    const real w = -weight;
    const real cc0 =  ptA.coef2(),  ww0 = w * cc0;
    const real cc1 =  ptA.coef1(),  ww1 = w * cc1;
    const real cc2 = -ptB.coef2(),  ww2 = w * cc2;
    const real cc3 = -ptB.coef1(),  ww3 = w * cc3;

    const real ee1 = side1 / ( 2 * ptA.len() ), we1 = w * ee1;
    const real ee2 = side2 / ( 2 * ptB.len() ), we2 = w * ee2;

    Matrix22 A(cc0, -ee1,  ee1, cc0);
    Matrix22 B(cc1,  ee1, -ee1, cc1);
    Matrix22 C(cc2, -ee2,  ee2, cc2);
    Matrix22 D(cc3,  ee2, -ee2, cc3);
    
    Matrix22 Aw(ww0, -ew1,  ew1, ww0);
    Matrix22 Bw(ww1,  ew1, -ew1, ww1);
    Matrix22 Cw(ww2, -ew2,  ew2, ww2);
    Matrix22 Dw(ww3,  ew2, -ew2, ww3);
    
    add_diag_block(ii0, A.trans_mul(Aw));
    add_diag_block(ii1, B.trans_mul(Bw));
    add_diag_block(ii2, C.trans_mul(Cw));
    add_diag_block(ii3, D.trans_mul(Dw));
 
    add_block(ii1, ii0, B.trans_mul(Aw));
    add_block(ii3, ii2, D.trans_mul(Cw));

    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, C.trans_mul(Aw));
        add_block(ii3, ii0, D.trans_mul(Aw));
        add_block(ii2, ii1, C.trans_mul(Bw));
        add_block(ii3, ii1, D.trans_mul(Bw));
    }
    else
    {
        add_block(ii0, ii2, A.trans_mul(Cw));
        add_block(ii1, ii2, B.trans_mul(Cw));
        add_block(ii0, ii3, A.trans_mul(Dw));
        add_block(ii1, ii3, B.trans_mul(Dw));
    }
 
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            add_base(ii0, Aw.trans_vecmul(off));
            add_base(ii1, Bw.trans_vecmul(off));
            add_base(ii2, Cw.trans_vecmul(off));
            add_base(ii3, Dw.trans_vecmul(off));
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(ee1,ptA.diff()), cross(ee2,ptB.diff()), ptB.pos());
    }
#endif
}
*/

void Meca::addSideSideLink2D(const Interpolation & ptA,
                             const Interpolation & ptB,
                             const real len,
                             const real weight,
                             real side1, real side2 )
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );
    
    //index in the matrix mB:
    index_t ia1 = ptA.matIndex1(), ia2 = ptA.matIndex2();
    index_t ib1 = ptB.matIndex1(), ib2 = ptB.matIndex2();
    
    if ( any_equal(ia1, ia2, ib1, ib2) )
        return;

    const real ca1 =  ptA.coef2(), ca2 =  ptA.coef1();
    const real cb1 = -ptB.coef2(), cb2 = -ptB.coef1();
    
    const real ee1 = side1 * len / ( 2 * ptA.len() );
    const real ee2 = side2 * len / ( 2 * ptB.len() );
    
    const real w = -weight;
    const real ca1w = ca1 * w, ca2w = ca2 * w;
    const real cb1w = cb1 * w, cb2w = cb2 * w;
   
    const real ee1w = ee1 * w, ee1ee1w = ee1 * ee1w;
    const real ee2w = ee2 * w, ee2ee2w = ee2 * ee2w;
    const real ee1ee2w = ee1 * ee2w;
    
    //we put the isotropic terms in mB
    add_iso(ia1, ia1, ca1w * ca1 + ee1ee1w);
    add_iso(ia1, ia2, ca1w * ca2 - ee1ee1w);
    add_iso(ia2, ia2, ca2w * ca2 + ee1ee1w);
    
    add_iso(ib1, ib1, cb1w * cb1 + ee2ee2w);
    add_iso(ib1, ib2, cb1w * cb2 - ee2ee2w);
    add_iso(ib2, ib2, cb2w * cb2 + ee2ee2w);
    
    add_iso(ia1, ib1, ca1w * cb1 - ee1ee2w);
    add_iso(ia1, ib2, ca1w * cb2 + ee1ee2w);
    add_iso(ia2, ib1, ca2w * cb1 + ee1ee2w);
    add_iso(ia2, ib2, ca2w * cb2 - ee1ee2w);
    
    //index in the matrix mC:
    ia1 *= DIM;
    ia2 *= DIM;
    ib1 *= DIM;
    ib2 *= DIM;
    
    mC(ia1  , ia2+1) -= ee1w;
    mC(ia1+1, ia2  ) += ee1w;
    
    mC(ib1  , ib2+1) -= ee2w;
    mC(ib1+1, ib2  ) += ee2w;
    
    const real ee1cb1w = ee1w * cb1;
    const real ee1cb2w = ee1w * cb2;
    const real ee2ca1w = ee2w * ca1;
    const real ee2ca2w = ee2w * ca2;
    
    mC(ia1, ib1+1) -=  ee2ca1w + ee1cb1w;
    mC(ia1, ib2+1) +=  ee2ca1w - ee1cb2w;
    
    mC(ia1+1, ib1) +=  ee2ca1w + ee1cb1w;
    mC(ia1+1, ib2) -=  ee2ca1w - ee1cb2w;
    
    mC(ia2, ib1+1) -=  ee2ca2w - ee1cb1w;
    mC(ia2, ib2+1) +=  ee2ca2w + ee1cb2w;
    
    mC(ia2+1, ib1) +=  ee2ca2w - ee1cb1w;
    mC(ia2+1, ib2) -=  ee2ca2w + ee1cb2w;
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(ee1,ptA.diff()), cross(ee2,ptB.diff()), ptB.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("addSideSideLink2D is not usable with periodic boundary conditions");
}

#endif


/**
 Link `ptA` (A) and `ptB` (B),
 but the links are maded between SA and SB which are located
 on the side of A and B, respectively:
 
     SA = A + len * N_A,
     SB = B + len * N_B,
 
 N_X is a normalized vector orthogonal to the fiber carrying X, in X:
 The force is linear of zero resting length,
 
     force_SA = weight * ( SA - SB )
     force_SB = weight * ( SB - SA )
 
 */

void Meca::addSideSideLink(const Interpolation & ptA,
                           const Interpolation & ptB,
                           const real len,
                           const real weight )
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideSideLink meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector dir = ptB.pos() - ptA.pos();
    real side1 = RNG.sign_exc( cross(ptA.diff(), dir) );
    real side2 = RNG.sign_exc( cross(dir, ptB.diff()) );
    addSideSideLink2D(ptA, ptB, len, weight, side1, side2);
    
#else
    
    throw Exception("Meca::addSideSideLink was not implemented in 3D");
    
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

void Meca::addSlidingLink(const Interpolation & ptA,
                          const Mecapoint & ptB,
                          const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    //force coefficients on the points:
    const real A = ptA.coef2();
    const real B = ptA.coef1();
    const real AA = A * A, AB = A * B, BB = B * B;
    
    Vector dir = ptA.dir();

    /*
     Points are (a, b, e) with (ab) the Interpolation, and e the Mecapoint,
     P is the projection on the plane perpendicular to (ab)
         P = I - dir (x) dir / normSqr(dir)
     the interaction is  -weigth * transpose(bb, aa, -1) * P * ( bb, aa, -1 )
     we set only the upper part of this symmetric matrix:
     */

    // T = -weight * [ I - dir (x) dir ]
    MatrixBlock T = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_diag_block(ii0, AA, T);
    add_diag_block(ii1, BB, T);
    add_diag_block(ii2, T);

    add_block(ii0, ii1, AB, T);
    add_block(ii0, ii2, -A, T);
    add_block(ii1, ii2, -B, T);
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off = weight * ( off - dot(dir, off) * dir );
            sub_base(ii0, A*off);
            sub_base(ii1, B*off);
            add_base(ii2, off);
        }
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

void Meca::addSlidingLink(const Interpolation & ptA,
                          const Interpolation & ptB,
                          const real weight)
{
    assert_true( weight >= 0 );
 
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    //interpolation coefficients
    const real ww[4] = { ptA.coef2(), ptA.coef1(), -ptB.coef2(), -ptB.coef1() };
    
    // on points (a, b, e), (ab) being the Interpolation, and e the Mecapoint,
    // P is the projection on the plane perpendicular to (ab): P.v= (v - (T.v)T/normSqr(T))
    // the interaction is  -wh' * P * h
    // we set only the upper part of this symmetric matrix:
    
    Vector dir = ptA.dir();

    // T = -weight * [ I - dir (x) dir ]
    MatrixBlock T = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    Matrix44 W = Matrix44::outerProduct(ww);
    
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_block(ii3, ii0, W(3,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_block(ii3, ii1, W(3,1), T);
    add_diag_block(ii2, W(2,2), T);
    add_block(ii3, ii2, W(3,2), T);
    add_diag_block(ii3, W(3,3), T);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off = weight * ( off - dot(dir, off) * dir );
            add_base(ii0, off, ww[0]);
            add_base(ii1, off, ww[1]);
            add_base(ii2, off, ww[2]);
            add_base(ii3, off, ww[3]);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Off-axis frictionless links (used for steric interactions)

#if ( DIM == 2 )

void Meca::addSideSlidingLink2D(const Interpolation & ptA,
                                const Mecapoint & ptB,
                                const real arm,
                                const real weight)
{
    assert_true( weight >= 0 );

    // indices
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    Vector dir = ptA.dir();
    const real aa = ptA.coef2();
    const real bb = ptA.coef1();
    const real ee = arm / ptA.len();
    
    // the (symmetric) projection matrix:
    // P = -weight * [ I - dir (x) dir ]
    MatrixBlock P = MatrixBlock::offsetOuterProduct(-weight, dir, weight);
    
    // anti-symmetric matrix blocks:
    const Matrix22 A( -aa,  ee, -ee, -aa );
    const Matrix22 B( -bb, -ee,  ee, -bb );

    /*
     We use block operations to set the matrix block by block:
     | A'PA  A'PB  A'P |
     | B'PA  B'PB  B'P |
     |   PA    PB    P |
     This matrix has symmetric and anti-symmetric blocks,
     since P' = P but A and B are not symmetric
     */
    assert_true( ii0 < ii1 );
    
    if ( ii2 > ii1 )
    {
        const Matrix22 PA = P.mul(A);
        const Matrix22 PB = P.mul(B);
        add_diag_block(ii0, A.trans_mul(PA));
        add_block(ii1, ii0, B.trans_mul(PA));
        add_block(ii2, ii0, PA);
        add_diag_block(ii1, B.trans_mul(PB));
        add_block(ii2, ii1, PB);
        add_diag_block(ii2, P);
    }
    else
    {
        // in this case, swap indices to address lower triangle
        const Matrix22 AtP = A.trans_mul(P);
        const Matrix22 BtP = B.trans_mul(P);
        add_diag_block(ii2, P);
        add_block(ii0, ii2, AtP);
        add_block(ii1, ii2, BtP);
        add_diag_block(ii0, AtP.mul(A));
        add_block(ii1, ii0, BtP.mul(A));
        add_diag_block(ii1, BtP.mul(B));
    }
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off = weight * ( off - dot(off, dir) * dir );
            add_base(ii0, A.trans_vecmul(off));
            add_base(ii1, B.trans_vecmul(off));
            add_base(ii2, off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}


/**
 Alternative 2D method in which we add an offset to vBAS
 */
void Meca::addSideSlidingLinkS(const Interpolation & ptA,
                               const Mecapoint & ptB,
                               const real arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    // indices
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    // set vector 'axi' perpendicular to Fiber:
    Vector axi = cross(1.0, ptA.dir());
    
    MatrixBlock T = MatrixBlock::outerProduct(axi);
    
    // we set directly the transformed offset vector:
    axi *= arm;
    
    // weights and indices:
    const real cc[3] = {   ptA.coef2(),   ptA.coef1(),   -1.0 };
    const real ww[3] = { -weight*cc[0], -weight*cc[1], weight };
    
    Matrix33 W = Matrix33::outerProduct(cc, ww);
    
    add_base(ii0, axi, ww[0]);
    add_base(ii1, axi, ww[1]);
    add_base(ii2, axi, ww[2]);

    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_diag_block(ii2, W(2,2), T);

    if ( modulo )
        throw Exception("addSideSlidingLinkS is not usable with periodic boundary conditions");
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), axi, ptB.pos());
    }
#endif
}


#elif ( DIM >= 3 )

/**
 Vector 'arm' must be parallel to the link and orthogonal to 'ptA'
 */

void Meca::addSideSlidingLinkS(const Interpolation & ptA,
                               const Mecapoint & ptB,
                               Vector const& arm,
                               const real len,
                               const real weight)
{    
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );
    
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex();

    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    /*
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
    
    Vector axi = arm / len;
    
    MatrixBlock T = MatrixBlock::outerProduct(axi);
    
    // weights and indices:
    const real w = -weight;
    const real cc[3] = { ptA.coef2(), ptA.coef1(),  -1.0 };
    const real ww[3] = { w * cc[0], w * cc[1], w * cc[2] };
    
    Matrix33 W = Matrix33::outerProduct(cc, ww);
    
    add_base(ii0, arm, ww[0]);
    add_base(ii1, arm, ww[1]);
    add_base(ii2, arm, ww[2]);

    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_diag_block(ii2, W(2,2), T);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off = dot(axi, off) * axi;
            add_base(ii0, off, ww[0]);
            add_base(ii1, off, ww[1]);
            add_base(ii2, off, ww[2]);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), arm, ptB.pos());
    }
#endif
    
}
#endif

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

 */
void Meca::addSideSlidingLink(const Interpolation & ptA,
                              const Mecapoint & ptB,
                              const real len,
                              const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector as = ptB.pos()-ptA.pos();
    if ( modulo )
        modulo->fold(as);
    real arm  = len * RNG.sign_exc( cross(ptA.diff(), as) );
    addSideSlidingLink2D(ptA, ptB, arm, weight);
    
#else
    
    // set 'arm' perpendicular to direction of the Fiber associated with `ptA`:
    Vector arm = calculateArm(ptB.pos()-ptA.pos(), ptA.diff(), len);
    addSideSlidingLinkS(ptA, ptB, arm, len, weight);
    
#endif
}

#pragma mark -


#if ( DIM == 2 )

// @todo interSideSlidingLink2D should use block operations

void Meca::addSideSlidingLink2D(const Interpolation & ptA,
                                const Interpolation & ptB,
                                const real arm,
                                const real weight)
{
    assert_true( weight >= 0 );

    const index_t inx[] = { DIM*ptA.matIndex1(),  DIM*ptA.matIndex1()+1,
                            DIM*ptA.matIndex2(),  DIM*ptA.matIndex2()+1,
                            DIM*ptB.matIndex1(),  DIM*ptB.matIndex1()+1,
                            DIM*ptB.matIndex2(),  DIM*ptB.matIndex2()+1 };

    if ( any_equal(inx[0], inx[2], inx[4], inx[6]) )
        return;

    Vector dir = ptA.dir();
    const real A1 =  ptA.coef2(), A2 =  ptA.coef1();
    const real B1 = -ptB.coef2(), B2 = -ptB.coef1();
    
    const real ee = arm / ptA.len();

    //this is done the 'hard' way by multiplying all matrices
    //coefficient matrix:
    real T[2*8] = { A1, -ee, ee, A1, A2, ee, -ee,  A2,
                    B1,   0,  0, B1, B2,  0,   0,  B2 };
    
    //the projection matrix:
    const real P[4] = { 1-dir.XX*dir.XX, -dir.XX*dir.YY, -dir.XX*dir.YY, 1-dir.YY*dir.YY };
    
    real PT[2*8], TPT[8*8];
    blas::xgemm('N','N', 2, 8, 2, -weight, P, 2, T, 2, 0.0, PT, 2);
    blas::xgemm('T','N', 8, 8, 2, 1.0, T, 2, PT, 2, 0.0, TPT, 8);
    
    //printf("\n");  VecPrint::print(8,8, TPT);
    
    for ( int ii=0; ii<8; ++ii )
    for ( int jj=ii; jj<8; ++jj )
        mC(inx[ii], inx[jj]) += TPT[ii+8*jj];
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            for ( int ii=0; ii<8; ++ii )
            {
                vBAS[inx[ii]] -= TPT[ii+8*4] * off.XX;
                vBAS[inx[ii]] -= TPT[ii+8*5] * off.YY;
                vBAS[inx[ii]] -= TPT[ii+8*6] * off.XX;
                vBAS[inx[ii]] -= TPT[ii+8*7] * off.YY;
            }
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}

/**
 Alternative 2D method in which we add an offset to vBAS
 */
void Meca::addSideSlidingLinkS(const Interpolation & ptA,
                               const Interpolation & ptB,
                               const real arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    // set vector 'axi' perpendicular to Fiber:
    Vector axi = cross(1.0, ptA.dir());
    
    MatrixBlock T = MatrixBlock::outerProduct(axi);
    
    // we set directly the transformed offset vector:
    axi *= arm;

    // weights and indices:
    const real cc[4] = {   ptA.coef2(),   ptA.coef1(),  -ptB.coef2(),  -ptB.coef1() };
    const real ww[4] = { -weight*cc[0], -weight*cc[1], -weight*cc[2], -weight*cc[3] };
    
    Matrix44 W = Matrix44::outerProduct(cc, ww);
    
    add_base(ii0, axi, ww[0]);
    add_base(ii1, axi, ww[1]);
    add_base(ii2, axi, ww[2]);
    add_base(ii3, axi, ww[3]);
    
    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_block(ii3, ii0, W(3,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_block(ii3, ii1, W(3,1), T);
    add_diag_block(ii2, W(2,2), T);
    add_block(ii3, ii2, W(3,2), T);
    add_diag_block(ii3, W(3,3), T);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off = off - dot(axi, off) * axi;
            add_base(ii0, off, ww[0]);
            add_base(ii1, off, ww[1]);
            add_base(ii2, off, ww[2]);
            add_base(ii3, off, ww[3]);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), axi, ptB.pos());
    }
#endif
}


#elif ( DIM >= 3 )

    /**
     Vector 'arm' must be parallel to the link and orthogonal to 'ptA'
     */
void Meca::addSideSlidingLinkS(const Interpolation & ptA,
                               const Interpolation & ptB,
                               Vector const& arm,
                               const real len,
                               const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len > REAL_EPSILON );
    
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t ii2 = DIM * ptB.matIndex1();
    const index_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    /*
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
    
    Vector axi = arm / len;
    
    MatrixBlock T = MatrixBlock::outerProduct(axi);
    
    // weights and indices:
    const real cc[4] = {   ptA.coef2(),   ptA.coef1(),   -ptB.coef2(), -ptB.coef1() };
    const real ww[4] = { -weight*cc[0], -weight*cc[1], -weight*cc[2], -weight*cc[3] };
    
    Matrix44 W = Matrix44::outerProduct(cc, ww);

    add_base(ii0, arm, ww[0]);
    add_base(ii1, arm, ww[1]);
    add_base(ii2, arm, ww[2]);
    add_base(ii3, arm, ww[3]);

    // fill the matrix mC
    add_diag_block(ii0, W(0,0), T);
    add_block(ii1, ii0, W(1,0), T);
    add_block(ii2, ii0, W(2,0), T);
    add_block(ii3, ii0, W(3,0), T);
    add_diag_block(ii1, W(1,1), T);
    add_block(ii2, ii1, W(2,1), T);
    add_block(ii3, ii1, W(3,1), T);
    add_diag_block(ii2, W(2,2), T);
    add_block(ii3, ii2, W(3,2), T);
    add_diag_block(ii3, W(3,3), T);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( !off.null() )
        {
            off = dot(axi, off) * axi;
            add_base(ii0, off, ww[0]);
            add_base(ii1, off, ww[1]);
            add_base(ii2, off, ww[2]);
            add_base(ii3, off, ww[3]);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLink(ptA.pos(), arm, ptB.pos());
    }
#endif
    
}

#endif

/**
 Link `ptA` (A) and `ptB` (B),
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

void Meca::addSideSlidingLink(const Interpolation & ptA,
                              const Interpolation & ptB,
                              const real len,
                              const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideSlidingLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real arm = len * RNG.sign_exc( cross(ptA.diff(), ptB.pos()-ptA.pos()) );
    addSideSlidingLink2D(ptA, ptB, arm, weight);
    
#else
    
    // set 'arm' perpendicular to direction of the Fiber associated with `ptA`:
    Vector arm = calculateArm(ptB.pos()-ptA.pos(), ptA.diff(), len);
    addSideSlidingLinkS(ptA, ptB, arm, len, weight);
    
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
    const index_t inx = ptA.matIndex();
    
    sub_iso(inx, inx, weight);
    
    if ( modulo )
        modulo->fold(pos, ptA.pos());

    add_base(DIM*inx, pos, weight);
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
    
    const real c1 = pti.coef2();
    const real c2 = pti.coef1();
    
    assert_true( 0 <= c1  &&  c1 <= 1 );
    assert_true( 0 <= c2  &&  c2 <= 1 );

    const real c2w = weight * c2;
    const real c1w = weight * c1;
    
    sub_iso(ii0, ii0, c1w * c1);
    sub_iso(ii0, ii1, c2w * c1);
    sub_iso(ii1, ii1, c2w * c2);
    
    if ( modulo )
        modulo->fold(pos, pti.pos());
    
    add_base(DIM*ii0, pos, c1w);
    add_base(DIM*ii1, pos, c2w);
}


//------------------------------------------------------------------------------
#pragma mark - Links to fixed sphere and cylinder
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

void Meca::addSphereClamp(Vector const& pos,
                          Mecapoint const& pte,
                          Vector const& center,
                          real rad,
                          const real weight)
{
    assert_true( rad >= 0 );
    assert_true( weight >= 0 );
    const index_t inx = DIM * pte.matIndex();
    
    Vector dir = pos - center;
    
    real len = dir.norm();
    
    if ( len > REAL_EPSILON )
        dir /= len;
    else
    {
        dir = Vector::randU();
        len = rad;
    }

    MatrixBlock T;
    if ( rad < len )
    {
        // point is outside sphere
        real R = weight * rad / len;
        // T = dia * I - len * [ I - dir (x) dir ]
        T = MatrixBlock::offsetOuterProduct(weight-R, dir, R);
        
        real facX = weight * rad + R * dot(dir, center);
        real facC = weight - R;
        add_base(inx, facX * dir + facC * center);
    }
    else
    {
        // point is inside sphere
        T = MatrixBlock::outerProduct(dir, weight);
        real facX = weight * ( rad + dot(dir, center) );
        add_base(inx, dir, facX);
    }
    
    sub_diag_block(inx, T);
    
    if ( modulo )
        throw Exception("addSphereClamp is not usable with periodic boundary conditions");
}


void Meca::addSphereClamp(Mecapoint const& pte,
                          Vector const& center,
                          real rad,
                          const real weight)
{
    addSphereClamp(pte.pos(), pte, center, rad, weight);
}


void Meca::addSphereClamp(Interpolation const& pti,
                          Vector const& center,
                          real rad,
                          const real weight)
{
    // interpolate on the two flanking vertices using coefficients:
    addSphereClamp(pti.pos(), pti.exact1(), center, rad, weight*pti.coef2());
    addSphereClamp(pti.pos(), pti.exact2(), center, rad, weight*pti.coef1());
}


/**
 Link `pte` (P) to a cylinder of axis X and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(P.XX, 0, 0)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force is only in the YZ plane.
 */

void Meca::addCylinderClampX(const Mecapoint & pte,
                             real  rad,
                             const real weight)
{
    assert_true( weight >= 0 );
    const index_t inx = DIM * pte.matIndex();
    
#if ( DIM == 2 )
    
    mC(inx+1, inx+1) -= weight;
    vBAS[inx+1]      += weight * std::copysign(rad, pte.pos().YY);
    
#elif ( DIM >= 3 )

    Vector pos = pte.pos();
    real dir_n = pos.normYZ();
    if ( dir_n < REAL_EPSILON )
        return;
    
    Vector dir(0, pos.YY/dir_n, pos.ZZ/dir_n);
    
    real facX;
    
    if ( rad < dir_n )
    {
        rad /= dir_n;
        facX = weight * rad * dir_n;
        
        mC(inx+1, inx+1) -= weight * ( 1.0 - rad * ( 1.0 - dir.YY * dir.YY ) );
        mC(inx+1, inx+2) -= weight * rad * dir.YY * dir.ZZ;
        mC(inx+2, inx+2) -= weight * ( 1.0 - rad * ( 1.0 - dir.ZZ * dir.ZZ ) );
    }
    else
    {
        facX = weight * rad;

        mC(inx+1, inx+1) -= weight * dir.YY * dir.YY;
        mC(inx+1, inx+2) -= weight * dir.YY * dir.ZZ;
        mC(inx+2, inx+2) -= weight * dir.ZZ * dir.ZZ;
    }
    
    // there should be no XX component here!
    vBAS[inx+1] += facX * dir.YY;
    vBAS[inx+2] += facX * dir.ZZ;
    
#endif
}


/**
 Link `pte` (P) to a cylinder of axis Z and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, 0, P.ZZ)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force is constrained in the XY plane.
 */

void Meca::addCylinderClampZ(const Mecapoint & pte,
                             real  rad,
                             const real weight)
{
    assert_true( weight >= 0 );
    
#if ( DIM > 1 )

    const index_t inx = DIM * pte.matIndex();
    Vector pos = pte.pos();
    real dir_n = pos.normXY();
    if ( dir_n < REAL_EPSILON )
        return;
    
    Vector dir(pos.XX/dir_n, pos.YY/dir_n, 0);

    real facX;
    
    if ( rad < dir_n )
    {
        rad /= dir_n;
        facX = weight * rad * dir_n;
        
        mC(inx  , inx  ) -= weight * ( 1.0 - rad * ( 1.0 - dir.XX * dir.XX ) );
        mC(inx  , inx+1) -= weight * rad * dir.XX * dir.YY;
        mC(inx+1, inx+1) -= weight * ( 1.0 - rad * ( 1.0 - dir.YY * dir.YY ) );
    }
    else
    {
        facX = weight * rad;
        
        mC(inx  , inx  ) -= weight * dir.XX * dir.XX;
        mC(inx  , inx+1) -= weight * dir.XX * dir.YY;
        mC(inx+1, inx+1) -= weight * dir.YY * dir.YY;
    }
    
    vBAS[inx  ] += facX * dir.XX;
    vBAS[inx+1] += facX * dir.YY;
    // there should be no ZZ component here!

#endif
}


//------------------------------------------------------------------------------
#pragma mark - Off-axis links to fixed positions
//------------------------------------------------------------------------------


#if ( DIM == 2 )

void Meca::addSidePointClamp2D(Interpolation const& ptA,
                               Vector const& pos,
                               const real arm,
                               const real weight)
{
    //force coefficients on the points:
    const real A = ptA.coef2(),  wA = weight * A;
    const real B = ptA.coef1(),  wB = weight * B;
    
    const real E = arm / ptA.len();
    const real wE = weight * E;
    const real wEE = weight * E * E;
    
    //index in the matrix mB:
    index_t ii0 = ptA.matIndex1();
    index_t ii1 = ptA.matIndex2();
    
    //we put the isotropic terms in mB
    sub_iso(ii0, ii0,  wA * A + wEE);
    sub_iso(ii0, ii1,  wA * B - wEE);
    sub_iso(ii1, ii1,  wB * B + wEE);
    
    //index in the matrix mC:
    ii0 *= DIM;
    ii1 *= DIM;
    
    mC(ii0  , ii1+1) += wE;
    mC(ii0+1, ii1  ) -= wE;
    
    //it seems to works also fine without the term in eew* below:
    vBAS[ii0  ] += wA * pos.XX - wE * pos.YY;
    vBAS[ii0+1] += wA * pos.YY + wE * pos.XX;
    vBAS[ii1  ] += wB * pos.XX + wE * pos.YY;
    vBAS[ii1+1] += wB * pos.YY - wE * pos.XX;

#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), pos);
    }
#endif
    
    if ( modulo )
        throw Exception("addSidePointClamp2D is not usable with periodic boundary conditions");
}

#elif ( DIM >= 3 )

/**
 A link of stiffness `weight`, between offset_point on the side of `ptA`, and the fixed position `pos` (G).

 This uses the vector product x -> cross(arm, x) to offset the point on which the link is attached:
 
     offset_point = fiber_point + cross(arm, fiber_dir)
 
 with fiber_point = ptA.pos() and fiber_dir = ptA.diff().normalized.
 `arm` must be perpendicular to link ( G - ptA.pos() )

 F. Nedelec, March 2011
 
 @todo addSidePointClamp3D should use block operations
 */
void Meca::addSidePointClamp3D(Interpolation const& ptA,
                               Vector const& pos,
                               Vector const& arm,
                               real const weight)
{
    real aa = ptA.coef2();
    real bb = ptA.coef1();
    
    real s = 1.0 / ptA.len();

    real ex = s * arm.XX;
    real ey = s * arm.YY;
    real ez = s * arm.ZZ;
    
    // indices to mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();
    const index_t inx[6] = { ii0, ii0+1, ii0+2, ii1, ii1+1, ii1+2 };
    
    /* The transfer matrix transforms the two Mecapoint in ptA,
     to the side point S:
     S = aa * pt1 + bb * pt2 + cross(arm, normalize( pt2 - pt1 ))
     
     It was generated in Maxima:
     MVP: matrix([0, -ez, ey], [ez, 0, -ex], [-ey, ex, 0]);
     MD: addcol(-ident(3), ident(3));
     MC: addcol(aa*ident(3), bb*ident(3));
     T: MC+MVP.MD;
     */
    const real T[18] = {
         aa,  ez, -ey,  bb, -ez,  ey,
        -ez,  aa,  ex,  ez,  bb, -ex,
         ey, -ex,  aa, -ey,  ex,  bb
    };
    
#if ( 0 )
    
    real TT[36];
    // TT = transpose(T) * T
    blas::xsyrk('U','N', 6, 3, 1.0, T, 6, 0.0, TT, 6);
    
#else
    
    real a2 = aa * aa;
    real b2 = bb * bb;
    real ab = aa * bb;
    
    real exx = ex * ex, exy = ex*ey, exz = ex*ez;
    real eyy = ey * ey, eyz = ey*ez;
    real ezz = ez * ez;
    
    // TT = transpose(T) * T is symmetric, and thus we only set half of it:
    /* Maxima code:
    TT: expand(transpose(T) . T);
     */
    real TT[36] = {
        eyy+ezz+a2,  0,           0,           0,           0,           0,
        -exy,        exx+ezz+a2,  0,           0,           0,           0,
        -exz,       -eyz,         exx+eyy+a2,  0,           0,           0,
        -ezz-eyy+ab, ez+exy,      exz-ey,      eyy+ezz+b2,  0,           0,
        -ez+exy,    -ezz-exx+ab,  eyz+ex,     -exy,         exx+ezz+b2,  0,
        exz+ey,      eyz-ex,     -eyy-exx+ab, -exz,        -eyz,         exx+eyy+b2
    };
    
#endif
    
    // we project to bring all forces in the plane perpendicular to 'arm'
    real sca = arm.inv_norm();
    real aan = aa * sca;
    real bbn = bb * sca;
    real TP[6] = { aan*ex, aan*ey, aan*ez, bbn*ex, bbn*ey, bbn*ez };
    
    //blas::xgemm('N','N', 6, 1, 3, sca, T, 6, arm, 3, 0.0, TP, 6);
    blas::xsyrk('U','N', 6, 1, weight, TP, 6, -weight, TT, 6);
    
    for ( int ii=0; ii<6; ++ii )
    for ( int jj=ii; jj<6; ++jj )
        mC(inx[ii], inx[jj]) += TT[ii+6*jj];
    
    // { gx, gy, gz } is the projection of `pos` in the plane perpendicular to 'arm'
    real ws = dot(arm, pos) * sca * sca;
    real gx = weight * ( pos.XX - ws * arm.XX );
    real gy = weight * ( pos.YY - ws * arm.YY );
    real gz = weight * ( pos.ZZ - ws * arm.ZZ );
    
    for ( int ii=0; ii<6; ++ii )
        vBAS[inx[ii]] += T[ii] * gx + T[ii+6] * gy + T[ii+12] * gz;
                 
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), pos);
    }
#endif
    
    if ( modulo )
        throw Exception("addSidePointClamp3D is not usable with periodic boundary conditions");
}

#endif  

/**
 Update Meca to include a connection between `ptA` (A) and a fixed position `pos` (G).
 The force is of zero resting length, but it is taken between G
 and another point S which is located on the side of the segment supporting A:
 
     S = A + len * N,
     force_S = weight * ( G - S )
 
 There is no counter-force in G, since G is immobile.
 */

void Meca::addSidePointClamp(Interpolation const& ptA,
                             Vector const& pos,
                             const real len,
                             const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSidePointClamp is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    // 'arm' is a vector in the Z direction
    real arm = len * RNG.sign_exc( cross(ptA.diff(), pos-ptA.pos()));
    addSidePointClamp2D(ptA, pos, arm, weight);
   
#else
    
    // 'arm' perpendicular to link and fiber is obtained by vector product:
    Vector arm = cross( ptA.pos()-pos, ptA.diff() );
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSidePointClamp3D(ptA, pos, arm * ( len / n ), weight);

#endif  
}

//------------------------------------------------------------------------------
#pragma mark - Links to lines and planes
//------------------------------------------------------------------------------

/**
 Link `ptA` (X) to the line defined by `G` and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 components parallel to `dir` are removed corresponding to a frictionless line:
 
     matrix M = 1 - dir (x) dir'
     force_X = weight * M * ( G - X )
 
 The vector `dir` should be of norm = 1.
 */

void Meca::addLineClamp(const Mecapoint & ptA,
                        Vector const& pos,
                        Vector const& dir,
                        const real weight )
{
    assert_true( weight >= 0 );
    
    const index_t inx = DIM * ptA.matIndex();
    
    
    // T = -weight * [ I - dir (x) dir ]
    MatrixBlock T = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_diag_block(inx, T);
    
    Vector off = weight * ( pos - dot(pos, dir) * dir );
    add_base(inx, off);
}


/**
 Link `ptA` and the line defined by `pos` (C) and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 components parallel to `dir` are removed corresponding to a frictionless line:
 
     M = I - dir (x) dir'
     force = weight * M * ( C - P )
 
 The vector `dir` should be of norm = 1.
 */

void Meca::addLineClamp(const Interpolation & ptA,
                        Vector const& pos,
                        Vector const& dir,
                        const real weight )
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();

    //force coefficients on the points:
    const real A = ptA.coef2();
    const real B = ptA.coef1();

    // T = -weight * [ I - dir (x) dir ]
    MatrixBlock T = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_diag_block(ii0, A*A, T);
    add_diag_block(ii1, B*B, T);
    add_block(ii0, ii1, A*B, T);
    
    //add the constant term:
    Vector off = weight * ( pos - dot(pos, dir) * dir );
    add_base(ii0, off, A);
    add_base(ii1, off, B);
}


/**
 Link `ptA` (X) and the plane defined one of its point `G` and the normal `dir`.
 The force is linear and the components parallel to the plane are removed,
 corresponding to an interaction with a frictionless plane:
 
     matrix M = dir (x) dir'
     force_X = weight * M * ( G - X )
 
 The vector `dir` should be of norm = 1, or alternatively
 `weight` can be the true weigth divided by |dir|^2.
 */

void Meca::addPlaneClamp(const Mecapoint & ptA,
                         Vector const& pos,
                         Vector const& dir,
                         const real weight )
{
    assert_true( weight >= 0 );
    
    const index_t inx = DIM * ptA.matIndex();
    
    // vBAS[inx] += dir * ( weigth * dot(pos,dir) );
    add_base(inx, dir, weight*dot(pos, dir));
    
#if ( DIM == 1 )
    mC(inx, inx) -= weight;
#else
    MatrixBlock T = MatrixBlock::outerProduct(dir, -weight);
    add_diag_block(inx, T);
#endif
}


/**
 Link `ptA` (X) and the plane defined by `pos` (G) and the normal `dir`.
 The force is linear and the perpendicular forces are removed, to create a frictionless plane:
 
     M = dir (x) dir'
     force = weight * M * ( G - X )
 
 The vector `dir` should be of norm = 1, or alternatively 
 `weight` can be the true weigth divided by |dir|^2.
 */

void Meca::addPlaneClamp(const Interpolation & ptA,
                         Vector const& pos,
                         Vector const& dir,
                         const real weight )
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const index_t ii0 = DIM * ptA.matIndex1();
    const index_t ii1 = DIM * ptA.matIndex2();

    //force coefficients on the points:
    const real A = ptA.coef2();
    const real B = ptA.coef1();
    
    //add the constant term:
    Vector off = ( weight * dot(pos, dir)) * dir;
    add_base(ii0, off, A);
    add_base(ii1, off, B);
    
    MatrixBlock T = MatrixBlock::outerProduct(dir, -weight);
    
    add_diag_block(ii0, A*A, T);
    add_diag_block(ii1, B*B, T);
    add_block(ii0, ii1, A*B, T);
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
void  Meca::addTriLink(Interpolation const& pt1, const real w1,
                       Interpolation const& pt2, const real w2,
                       Interpolation const& pt3, const real w3)
{
    const real sum = w1 + w2 + w3;
    assert_true( sum > REAL_EPSILON );
    addLink(pt1, pt2, w1*w2/sum);
    addLink(pt1, pt3, w1*w3/sum);
    addLink(pt2, pt3, w2*w3/sum);
}


/**
 Do not use this function!
 
 If `weigth > 0`, this creates an attractive force that decreases like 1/R^3
 */
void Meca::addCoulomb( const Mecapoint & ptA, const Mecapoint & ptB, real weight )
{
    Vector ab = ptB.pos() - ptA.pos();
    real abnSqr = ab.normSqr(), abn=sqrt(abnSqr);
    
    const index_t inxA = DIM * ptA.matIndex();
    const index_t inxB = DIM * ptB.matIndex();
    
    if ( abn < REAL_EPSILON ) return;
    ab /= abn;
    
    real abn3 = weight / abnSqr;
    real abn5 = weight / ( abnSqr * abn );
    
    add_base(inxA, ab,-3*abn3);
    add_base(inxB, ab, 3*abn3);
    
    for ( int ii = 0; ii < DIM; ++ii )
    {
        for ( int jj = ii; jj < DIM; ++jj )
        {
            real m = abn5 * ( (ii==jj) - 3 * ab[ii] * ab[jj] );
            
            mC(inxA+ii, inxA+jj) -= m;
            mC(inxB+ii, inxB+jj) -= m;
            
            mC(inxA+ii, inxB+jj) += m;
            if ( ii != jj )
                mC(inxA+jj, inxB+ii) += m;
        }
    }
}

