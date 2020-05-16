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

/**
 Note that you can disable all the modulo code with this definition.
 The code under 'if ( modulo )' will never be executed, and should even
 be discarded during compilation if optimizations are enabled.
 */
//constexpr Modulo * modulo = nullptr;


/// set TRUE to update matrix mC using block directives
/** This is significantly faster on machine with the AVX instruction set */
#define USE_MATRIX_BLOCK 1


//------------------------------------------------------------------------------
#pragma mark - Accessory functions

#if DRAW_MECA_LINKS

#  include "gle.h"
#  include "gle_color_list.h"

/// this performs the modulo on `c`
void drawLinkM(Vector const& a, Vector const& ab, Vector c)
{
    if ( modulo )
        c = modulo->image(c, a);
    gle::drawLink(a, ab, c);
}

#endif


/// true if any two values are equal
inline bool any_equal(const size_t a, const size_t b,
                      const size_t c)
{
    return ( a == c ) | ( b == c );
}


/// true if any two values are equal
inline bool any_equal(const size_t a, const size_t b,
                      const size_t c, const size_t d)
{
    return ( a == c ) | ( a == d ) | ( b == c ) | ( b == d );
}


//------------------------------------------------------------------------------
#pragma mark - Functions to set matrix elements
//------------------------------------------------------------------------------

#define CHECK_INDICES(I,J,C) ((void) 0)
//#define CHECK_INDICES(I,J,C) { if (J>I) printf(" wrong-sided %s %lu %lu\n",C,I,J); }

#define PRINT_BLOCK(I,J,B) ((void) 0)
//#define PRINT_BLOCK(I,J,B) { std::clog<<std::setw(3)<<I<<" "<<J<<" "<<std::setw(10)<<B<<'\n'; }


// add alpha * T to mC.
inline void Meca::add_block(size_t i, size_t j, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"add");
#if USE_MATRIX_BLOCK
    assert_true( i > j );
    mC.block(i, j).add_full(T);
    PRINT_BLOCK(i,j,T);
#elif ( DIM == 1 )
    mB(i,j) += T.value();
#else
    assert_true( i > j );
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = 0; y < DIM; ++y )
        mC(i+y, j+x) += T(y,x);
#endif
}

// add T to mC.
inline void Meca::add_block(size_t i, size_t j, real alpha, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"add_alpha");
#if USE_MATRIX_BLOCK
    assert_true( i > j );
    mC.block(i, j).add_full(alpha, T);
    PRINT_BLOCK(i,j,alpha*T);
#elif ( DIM == 1 )
    mB(i,j) += alpha * T.value();
#else
    assert_true( i > j );
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = 0; y < DIM; ++y )
        mC(i+y, j+x) += alpha * T(y,x);
#endif
}

// subtract T to mC.
inline void Meca::sub_block(size_t i, size_t j, MatrixBlock const& T)
{
    CHECK_INDICES(i,j,"sub");
#if USE_MATRIX_BLOCK
    assert_true( i > j );
    mC.block(i, j).sub_full(T);
    PRINT_BLOCK(i,j,-T);
#elif ( DIM == 1 )
    mB(i,j) -= T.value();
#else
    assert_true( i > j );
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = 0; y < DIM; ++y )
        mC(i+y, j+x) -= T(y,x);
#endif
}

// add T to the diagonal of mC. `T` should be symmetric
inline void Meca::add_block_diag(size_t i, MatrixBlock const& T)
{
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mC.diag_block(i).add_half(T);
    PRINT_BLOCK(i,i,T);
#elif ( DIM == 1 )
    mB(i,i) += T.value();
#else
    // add lower part of block
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = x; y < DIM; ++y )
        mC(i+y, i+x) += T(y,x);
#endif
}

// add alpha * T to the diagonal of mC. `T` should be symmetric
inline void Meca::add_block_diag(size_t i, real alpha, MatrixBlock const& T)
{
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mC.diag_block(i).add_half(alpha, T);
    PRINT_BLOCK(i,i,alpha*T);
#elif ( DIM == 1 )
    mB(i,i) += alpha * T.value();
#else
    // add lower part of block
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = x; y < DIM; ++y )
        mC(i+y, i+x) += alpha * T(y,x);
#endif
}

// add -T to the diagonal of mC. `T` should be symmetric
inline void Meca::sub_block_diag(size_t i, MatrixBlock const& T)
{
#if USE_MATRIX_BLOCK
    assert_small(T.asymmetry());
    mC.diag_block(i).sub_half(T);
    PRINT_BLOCK(i,i,-T);
#elif ( DIM == 1 )
    mB(i,i) -= T.value();
#else
    // add lower part of block
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = x; y < DIM; ++y )
        mC(i+y, i+x) -= T(y,x);
#endif
}


// add val to mB, the XYZ-isometric component
inline void Meca::add_iso(size_t i, size_t j, real val)
{
    CHECK_INDICES(i,j,"add_iso");
#if USE_ISO_MATRIX
    mB(i,j) += val;
#else
    size_t ii = DIM * std::max(i, j);
    size_t jj = DIM * std::min(i, j);
    mC.block(ii, jj).add_diag(val);
#endif
}

// add -val to mB, the XYZ-isometric component
inline void Meca::sub_iso(size_t i, size_t j, real val)
{
    CHECK_INDICES(i,j,"sub_iso");
#if USE_ISO_MATRIX
    mB(i,j) -= val;
#else
    size_t ii = DIM * std::max(i, j);
    size_t jj = DIM * std::min(i, j);
    mC.block(ii, jj).sub_diag(val);
#endif
}


inline void Meca::add_base(size_t i, Vector const& vec)
{
    assert_true( i % DIM == 0 );
    vec.add_to(vBAS+i);
}

inline void Meca::add_base(size_t i, Vector const& vec, real alpha)
{
    assert_true( i % DIM == 0 );
    vec.add_to(alpha, vBAS+i);
}

inline void Meca::sub_base(size_t i, Vector const& vec)
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
void Meca::addForce(Mecapoint const& pte, Vector const& force)
{
    const size_t inx = DIM * pte.matIndex();
    add_base(inx, force);
}


/**
Add constant force to an interpolated position
 */
void Meca::addForce(Interpolation const& pti, Vector const& force)
{
    const size_t ii0 = DIM * pti.matIndex1();
    const size_t ii1 = DIM * pti.matIndex2();
    
    add_base(ii0, force, pti.coef0());
    add_base(ii1, force, pti.coef1());
}


void Meca::addForceToAll(Vector const& force)
{
    for ( size_t p = 0; p < nb_points(); ++p )
        add_base(DIM*p, force);
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
    const size_t ii0 = DIM * pti.matIndex1();
    const size_t ii1 = DIM * pti.matIndex2();
    
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
void Meca::addTorqueClamp(Interpolation const& pti,
                          Vector const& dir,
                          const real weight)
{
    assert_true( weight >= 0 );
    const size_t ii0 = DIM * pti.matIndex1();
    const size_t ii1 = DIM * pti.matIndex2();
    
    Vector d = pti.diff();
    real n = d.normSqr();

    Torque Tq = cross(d, dir);

#if ( DIM >= 3 )
    const real Tn = Tq.norm();
#else
    const real Tn = abs_real(Tq);
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
void Meca::addTorqueExplicit(Interpolation const& ptA,
                             Interpolation const& ptB,
                             const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
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
    const real Tn = abs_real(Tq);
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
void Meca::addTorqueExplicit(Interpolation const& ptA,
                             Interpolation const& ptB,
                             const real cosinus, const real sinus,
                             const real weight)
{
    assert_true( weight >= 0 );
    assert_small( cosinus*cosinus + sinus*sinus - 1.0 );

    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
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
    Vector axis = cross(db, da).normalized(sign_real(sinus));
    
    // rotate vector `db` around `arm` by angle specified as (cosinus, sinus):
    Vector rot = cosinus * db + sinus * cross(axis, db);
    
#elif ( DIM == 2 )

    // this correspond to the Z-direction, up or down:
    real dir = sign_real(cross(da, db));

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
    const real Tn = abs_real(Tq);
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
void Meca::addTorquePoliti(Interpolation const& pt1,
                           Interpolation const& pt2,
                           const real cosinus, const real sinus,
                           const real weight)
{
    assert_true( weight >= 0 );
    assert_small( cosinus*cosinus + sinus*sinus - 1.0 );

    if ( pt1.overlapping(pt2) )
        return;
    
    //index in the matrix mC:
    const size_t index[] = { DIM*pt1.matIndex1(), DIM*pt1.matIndex1()+1,
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
                mC(index[ii], index[jj]) += m[mm];
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

    // indices in matrix mC:
    const size_t iiA = DIM * pt1.matIndex1();
    const size_t iiB = DIM * pt1.matIndex2();
    const size_t iiC = DIM * pt2.matIndex1();
    const size_t iiD = DIM * pt2.matIndex2();
    
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


void Meca::addTorque(Interpolation const& pt1,
                     Interpolation const& pt2,
                     const real cosinus, const real sinus,
                     const real weight)
{
    assert_true( weight >= 0 );
    assert_small( cosinus*cosinus + sinus*sinus - 1.0 );
    
#if ( DIM == 3 )
    const Vector AB = pt1.diff();
    const Vector CD = pt2.diff();
    Vector axis = cross(AB, CD);
    real n = axis.norm();
    if ( n > REAL_EPSILON )
        axis /= n;
    else
        axis = Vector::randU();
    addTorque(pt1, pt2, MatrixBlock::rotationAroundAxis(axis, cosinus, sinus), weight);
#elif ( DIM == 2 )
    addTorque(pt1, pt2, MatrixBlock(cosinus, sinus, -sinus, cosinus), weight);
#endif
}

#endif



MatrixBlock Meca::torqueMatrix(real weight, Torque const& axi, real cosinus, real sinus)
{
#if ( DIM == 3 )
    return (-weight) * MatrixBlock::rotationAroundAxis(axi, cosinus, sinus);
#elif ( DIM == 2 )
    return (-weight) * MatrixBlock(cosinus, axi*sinus, -axi*sinus, cosinus);
#else
    return MatrixBlock(-weight);  //should not be used!
#endif
}



/**
 Add Torque between 3 points.
 This version does not impose any particular distance between the points,
 and just move them to enforce the angle described by ABC
 */
void Meca::addTorque(Mecapoint const& ptA,
                     Mecapoint const& ptB,
                     Mecapoint const& ptC,
                     const MatrixBlock & R, //already multiplied by -weight
                     const real weight)
{
    assert_true( weight >= 0 );
    const MatrixBlock W(0, -weight);

/*
#if ( DIM == 3 )
    const Vector3 AB = ptB.pos() - ptA.pos();
    const Vector3 BC = ptC.pos() - ptB.pos();
    Vector3 axi = normalize(cross(AB, BC));
    if ( axi != axi )
        return;
    const MatrixBlock R = -weight * MatrixBlock::rotationAroundAxis(axi, cosinus, sinus);
#elif ( DIM == 2 )
    const Matrix22 R = -weight * Matrix22(cosinus, sinus,-sinus, cosinus);
#else
    const Matrix11 R(-weight);  //should not be used!
#endif
*/
    const MatrixBlock T = R.transposed();

    // indices in matrix mC:
    const size_t iiA = DIM * ptA.matIndex();
    const size_t iiB = DIM * ptB.matIndex();
    const size_t iiC = DIM * ptC.matIndex();

    /*
    Vector CD = ptB.pos() + R * AB - ptC.pos();
    Vector fA = T * CD;
    Vector fB = ( W + T ) * ( -CD );
    Vector fC = weight * CD;
     */
    
    add_block_diag(iiA, W);
    if ( iiB > iiA )
        sub_block(iiB, iiA, W+R);
    else
        sub_block(iiA, iiB, W+T);
    if ( iiC > iiA )
        add_block(iiC, iiA, R);
    else
        add_block(iiA, iiC, T);
    
#if ( DIM == 2 )
    // small optimization in 2D, as the term (R+W)+(T+W) is diagonal
    real dd = R.trace() - 2.0 * weight;
    add_block_diag(iiB, MatrixBlock(0, dd));
#else
    add_block_diag(iiB, (R+W)+(T+W));
#endif
    
    if ( iiC > iiB )
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
void Meca::addTorquePlane(Mecapoint const& ptA,
                          Mecapoint const& ptB,
                          Mecapoint const& ptC,
                          const Torque & axi,
                          const real cosinus, const real sinus,
                          const real weight)
{
    assert_true( weight >= 0 );
    assert_small( cosinus*cosinus + sinus*sinus - 1.0 );

#if ( DIM == 3 )
    /*
    const Vector3 AB = ptB.pos() - ptA.pos();
    const Vector3 BC = ptC.pos() - ptB.pos();
    Vector3 axi = normalize(cross(AB, BC));
    if ( axi != axi )
        return;
     */
    const MatrixBlock X = MatrixBlock::outerProduct(axi);
    const MatrixBlock R = MatrixBlock::rotationAroundAxis(axi, cosinus, sinus) - X;
    const MatrixBlock P = MatrixBlock(0,1) - X;
#elif ( DIM == 2 )
    const MatrixBlock R = MatrixBlock(cosinus, sinus,-sinus, cosinus);
    const MatrixBlock P(0,1);
#else
    const MatrixBlock R(1.0);  //should not be used!
    const MatrixBlock P(0,1);
#endif
    
    // indices in matrix mC:
    const size_t iiA = DIM * ptA.matIndex();
    const size_t iiB = DIM * ptB.matIndex();
    const size_t iiC = DIM * ptC.matIndex();

    const MatrixBlock wP = -weight * P;
    const MatrixBlock wR = -weight * R;
    const MatrixBlock wT = wR.transposed();

    add_block_diag(iiA, wP);
    if ( iiB > iiA )
        sub_block(iiB, iiA, wR+wP);
    else
        sub_block(iiA, iiB, (wR+wP).transposed());
    if ( iiC > iiA )
        add_block(iiC, iiA, wR);
    else
        add_block(iiA, iiC, wT);

    add_block_diag(iiB, (wP+wR)+(wP+wT));
    if ( iiC > iiB )
        sub_block(iiC, iiB, wP+wR);
    else
        sub_block(iiB, iiC, wP+wT);

    add_block_diag(iiC, wP);
}



/** This is variation 3, 20.08.2019
 It combines addTorque() without length with a LongLink(ptA, ptB);
 */
void Meca::addTorqueLong(Mecapoint const& ptA,
                         Mecapoint const& ptB,
                         Mecapoint const& ptC,
                         const MatrixBlock & R,
                         const real weight,
                         const real len, const real weightL)
{
    assert_true( weight >= 0 );
    assert_true( weightL >= 0 );
    const MatrixBlock W(0, -weight);
    const Vector AB = ptB.pos() - ptA.pos();
    
/*
#if ( DIM == 3 )
    const Vector3 BC = ptC.pos() - ptB.pos();
    Vector3 axi = normalize(cross(AB, BC));
    if ( axi != axi )
        return;
    const MatrixBlock R = -weight * MatrixBlock::rotationAroundAxis(axi, cosinus, sinus);
#elif ( DIM == 2 )
    const Matrix22 R = -weight * Matrix22(cosinus, sinus,-sinus, cosinus);
#else
    const Matrix11 R(-weight);  //should not be used!
#endif
*/

    const MatrixBlock T = R.transposed();
    
    // indices in matrix mC:
    const size_t iiA = DIM * ptA.matIndex();
    const size_t iiB = DIM * ptB.matIndex();
    const size_t iiC = DIM * ptC.matIndex();
    
    // this is a LongLink(A, B):
    MatrixBlock wL(0,0);
    const real ab2 = AB.normSqr();
    if ( ab2 > REAL_EPSILON )
    {
        real ab = sqrt(ab2);
        const real wla = weightL * len / ab;
        add_base(iiA, AB, -wla);
        add_base(iiB, AB,  wla);

        // regularize interaction to avoid creating negative eigen values
        if ( ab < len )
            wL = MatrixBlock::outerProduct(AB, -weightL/ab2);
        else
            wL = MatrixBlock::offsetOuterProduct(wla-weightL, AB, -wla/ab2);
    }
    
    add_block_diag(iiA, W+wL);
    if ( iiB > iiA )
        sub_block(iiB, iiA, W+R+wL);
    else
        sub_block(iiA, iiB, W+T+wL);
    if ( iiC > iiA )
        add_block(iiC, iiA, R);
    else
        add_block(iiA, iiC, T);

#if ( DIM == 2 )
    // in 2D, the term (R+W)+(T+W) is diagonal
    real dd = R.trace() - 2.0 * weight;
    add_block_diag(iiB, wL+MatrixBlock(0,dd));
#else
    add_block_diag(iiB, (R+W)+(T+W)+wL);
#endif
    if ( iiC > iiB )
        sub_block(iiC, iiB, W+R);
    else
        sub_block(iiB, iiC, W+T);

    add_block_diag(iiC, W);
}


//------------------------------------------------------------------------------
#pragma mark - Interpolation of positions
//------------------------------------------------------------------------------


Vector Meca::position1(const size_t inx) const
{
    return Vector(vPTS+inx);
}

Vector Meca::position2(const size_t inx[2], const real coef[2]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    return coef[0] * P0 + coef[1] * P1;
}


Vector Meca::position3(const size_t inx[3], const real coef[3]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    return ( coef[0] * P0 + coef[1] * P1 ) + coef[2] * P2;
}


Vector Meca::position4(const size_t inx[4], const real coef[4]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    Vector P3(vPTS+inx[3]);
    return ( coef[0] * P0 + coef[1] * P1 ) + ( coef[2] * P2 + coef[3] * P3 );
}


Vector Meca::position5(const size_t inx[5], const real coef[5]) const
{
    Vector P0(vPTS+inx[0]);
    Vector P1(vPTS+inx[1]);
    Vector P2(vPTS+inx[2]);
    Vector P3(vPTS+inx[3]);
    Vector P4(vPTS+inx[4]);
    return ( coef[0] * P0 + coef[1] * P1 ) + ( coef[2] * P2 + coef[3] * P3 ) + coef[4] * P4;
}


Vector Meca::position6(const size_t inx[6], const real coef[6]) const
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
void Meca::addLink(Mecapoint const& ptA,
                   Mecapoint const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    const size_t ii0 = ptA.matIndex();
    const size_t ii1 = ptB.matIndex();

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
        if ( off.is_not_zero() )
        {
            add_base(DIM*ii0, off, ww[0]);
            add_base(DIM*ii1, off, ww[1]);
        }
#else
        const size_t inx[] = { DIM*ii0, DIM*ii1 };
        Vector off = modulo->offset(position2(inx, ww));
        if ( off.is_not_zero() )
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
void Meca::addLink(Interpolation const& ptA,
                   Mecapoint const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = ptA.matIndex1();
    const size_t ii1 = ptA.matIndex2();
    const size_t ii2 = ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;

    //coefficients on the points:
    const real cc[] = { ptA.coef0(),   ptA.coef1(),    -1.0 };
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
#endif
        if ( off.is_not_zero() )
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
void Meca::addLink(Mecapoint const& ptA,
                   Interpolation const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = ptA.matIndex();
    const size_t ii1 = ptB.matIndex1();
    const size_t ii2 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    //coefficients on the points:
    const real cc[] = {    1.0, -ptB.coef0(), -ptB.coef1() };
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
#endif
        if ( off.is_not_zero() )
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
void Meca::addLink(Interpolation const& ptA,
                   Interpolation const& ptB,
                   const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = ptA.matIndex1();
    const size_t ii1 = ptA.matIndex2();
    const size_t ii2 = ptB.matIndex1();
    const size_t ii3 = ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    //interpolation coefficients:
    const real cc[] = {  ptA.coef0(),  ptA.coef1(), -ptB.coef0(), -ptB.coef1() };
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3 };
        Vector off = modulo->offset(position4(inx, cc));
#endif
        if ( off.is_not_zero() )
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
#pragma mark - Links between Mecables (higher order interpolation)
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
void Meca::addLink1(Interpolation const& pti,
                    const size_t pts,
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = pts;
    const size_t ii1 = pti.matIndex1();
    const size_t ii2 = pti.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    const real cc[] = {    1.0, -pti.coef0(), -pti.coef1() };
    const real ww[] = { weight, weight*cc[1], weight*cc[2] };
    
    sub_iso(ii0, ii0, ww[0]); // since cc[0] == 1.0
    sub_iso(ii1, ii0, ww[1]); // since cc[0] == 1.0
    sub_iso(ii2, ii0, ww[2]); // since cc[0] == 1.0
    
    sub_iso(ii1, ii1, ww[1] * cc[1]);
    sub_iso(ii2, ii1, ww[2] * cc[1]);
    
    sub_iso(ii2, ii2, ww[2] * cc[2]);
    
    if ( modulo )
    {
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
        if ( off.is_not_zero() )
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
void Meca::addLink2(Mecapoint const& ptA,
                    const size_t pts[2], const real coef[2],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = ptA.matIndex();
    const size_t ii1 = pts[0];
    const size_t ii2 = pts[1];
    
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2 };
        Vector off = modulo->offset(position3(inx, cc));
        if ( off.is_not_zero() )
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
void Meca::addLink2(Interpolation const& pti,
                    const size_t pts[2], const real coef[2],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = pti.matIndex1();
    const size_t ii1 = pti.matIndex2();
    const size_t ii2 = pts[0];
    const size_t ii3 = pts[1];
    
    const real cc[] = { -pti.coef0(), -pti.coef1(),      coef[0],      coef[1] };
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3 };
        Vector off = modulo->offset(position4(inx, cc));
        if ( off.is_not_zero() )
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
void Meca::addLink3(Mecapoint const& ptA,
                    const size_t pts[3], const real coef[3],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = ptA.matIndex();
    const size_t ii1 = pts[0];
    const size_t ii2 = pts[1];
    const size_t ii3 = pts[2];

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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3 };
        Vector off = modulo->offset(position4(inx, cc));
        if ( off.is_not_zero() )
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
void Meca::addLink3(Interpolation const& pti,
                    const size_t pts[3], const real coef[3],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = pti.matIndex1();
    const size_t ii1 = pti.matIndex2();
    const size_t ii2 = pts[0];
    const size_t ii3 = pts[1];
    const size_t ii4 = pts[2];

    const real cc[] = { -pti.coef0(), -pti.coef1(),      coef[0],      coef[1],      coef[2] };
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3, DIM*ii4 };
        Vector off = modulo->offset(position5(inx, cc));
        if ( off.is_not_zero() )
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
void Meca::addLink4(Mecapoint const& ptA,
                    const size_t pts[4], const real coef[4],
                    const real weight)
{
    assert_true( weight >= 0 );

    //index in the matrix mB:
    const size_t ii0 = ptA.matIndex();
    const size_t ii1 = pts[0];
    const size_t ii2 = pts[1];
    const size_t ii3 = pts[2];
    const size_t ii4 = pts[3];

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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3, DIM*ii4 };
        Vector off = modulo->offset(position5(inx, cc));
        if ( off.is_not_zero() )
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
void Meca::addLink4(Interpolation const& pti,
                    const size_t pts[4], const real coef[4],
                    const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    const size_t ii0 = pti.matIndex1();
    const size_t ii1 = pti.matIndex2();
    const size_t ii2 = pts[0];
    const size_t ii3 = pts[1];
    const size_t ii4 = pts[2];
    const size_t ii5 = pts[3];

    const real cc[] = { -pti.coef0(), -pti.coef1(),      coef[0],      coef[1],      coef[2],      coef[3] };
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
        const size_t inx[] = { DIM*ii0, DIM*ii1, DIM*ii2, DIM*ii3, DIM*ii4, DIM*ii5 };
        Vector off = modulo->offset(position6(inx, cc));
        if ( off.is_not_zero() )
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

void Meca::addLongLink(Mecapoint const& ptA,
                       Mecapoint const& ptB,
                       const real len,
                       const real weight)
{
    assert_true( weight >= 0 );
    assert_true( len >= 0 );

    const size_t ia = DIM * ptA.matIndex();  // coef is +weight
    const size_t ib = DIM * ptB.matIndex();  // coef is -weight

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

    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = sqrt(ab2);

    const real wla = weight * len / abn;

    add_base(ia, axi,-wla);
    add_base(ib, axi, wla);
    
    MatrixBlock wT;

    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */

    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);
    
    add_block_diag(ia, wT);
    add_block_diag(ib, wT);
    if ( ia > ib )
        sub_block(ia, ib, wT);
    else
        sub_block(ib, ia, wT);
    
    if ( modulo && off.is_not_zero() )
    {
        off = wT * off;
        add_base(ia, off);
        sub_base(ib, off);
    }
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

    //index in the matrix mC:
    const size_t ii0 = DIM * ptB.matIndex1();
    const size_t ii1 = DIM * ptB.matIndex2();
    const size_t ii2 = DIM * ptA.matIndex();

    if ( any_equal(ii0, ii1, ii2) )
        return;

    Vector off, axi = ptA.pos() - ptB.pos();

    if ( modulo )
        modulo->foldOffset(axi, off);
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptB.mecable()->signature()).load();
        gle::drawLink(ptB.pos(), axi, len);
    }
#endif
    
    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = sqrt(ab2);

    const real wla = weight * len / abn;

    // coefficients:
    const real cc0 = ptB.coef0();
    const real cc1 = ptB.coef1();

    add_base(ii0, axi, -cc0 * wla);
    add_base(ii1, axi, -cc1 * wla);
    add_base(ii2, axi, wla);

    MatrixBlock wT;

    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */

    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);

    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block(ii2, ii0,    -cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii2, ii1,    -cc1, wT);
    add_block_diag(ii2, wT);

    if ( modulo && off.is_not_zero() )
    {
        off = wT * off;
        add_base(ii0, off, cc0);
        add_base(ii1, off, cc1);
        add_base(ii2, off);
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

    //indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();

    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();

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

    const real ab2 = axi.normSqr();
    if ( ab2 < REAL_EPSILON ) return;
    const real abn = sqrt(ab2);

    const real wla = weight * len / abn;

    add_base(ii0, axi, -cc0 * wla);
    add_base(ii1, axi, -cc1 * wla);
    add_base(ii2, axi, -cc2 * wla);
    add_base(ii3, axi, -cc3 * wla);
    
    /* To stabilize the matrix with compression, we remove negative eigenvalues
     This is done by using len = 1 in the formula for links that are shorter
     than the desired target. */
    
    MatrixBlock wT;
    
    if ( len > abn )
        wT = MatrixBlock::outerProduct(axi, -weight/ab2);
    else
        wT = MatrixBlock::offsetOuterProduct(wla-weight, axi, -wla/ab2);
    
    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block_diag(ii2, cc2*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block(ii3, ii2, cc3*cc2, wT);

    if ( modulo && off.is_not_zero() )
    {
        off = wT * off;
        add_base(ii0, off, cc0);
        add_base(ii1, off, cc1);
        add_base(ii2, off, cc2);
        add_base(ii3, off, cc3);
   }
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
    
    // indices in the matrix mB and mC:
    const size_t ia0 = ptA.matIndex1(),  ii0 = DIM * ia0;
    const size_t ia1 = ptA.matIndex2(),  ii1 = DIM * ia1;
    const size_t ib2 = ptB.matIndex(),   ii2 = DIM * ib2;

    if ( any_equal(ii0, ii1, ii2) )
        return;
    const real ee = arm / ptA.len(), we = weight * ee;

    // coefficients and weights:
    const real cc0 = ptA.coef0(),  ww0 = weight * cc0;
    const real cc1 = ptA.coef1(),  ww1 = weight * cc1;
    
    sub_iso(ia0, ia0, ww0 * cc0 + we * ee);
    sub_iso(ia1, ia1, ww1 * cc1 + we * ee);
    sub_iso(ib2, ib2, weight);
    
    real wd = ww0 * cc1 - we * ee;
    sub_block(ii1, ii0, Matrix22(wd, -we, we, wd));
    
    if ( ii2 > ii0 )
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
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            Matrix22 waT(ww0,  we, -we, ww0);
            Matrix22 wbT(ww1, -we,  we, ww1);
            sub_base(ii0, waT*off);
            sub_base(ii1, wbT*off);
            add_base(ii2, off, weight);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}

#endif


/**
 Link `B` to an interpolated point `S` on the side of `A`:
 
     S = pos(A) + cross( arm, A.dir() )
 
 Where A.dir() is the direction of the Fiber supporting `A`, in `A`.
 The vector `arm` should ideally be perpendicular to A.dir(), and in this case,
 `A` and `B` are separated by norm(arm). The force is linear:

     force_B = weight * ( S - B )
     force_S = -force_B
 
 The `force_S` is redistributed on the vertices on each side of `A`,
 according to the interpolation coefficients, as usual.
 
 This code is valid in any dimension and works in 2 and 3D
 */
void Meca::addSideLink3D(Interpolation const& ptA,
                         Mecapoint const& ptB,
                         Torque const& arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    // coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    const real eps = -1.0 / ptA.len();

    MatrixBlock aR = MatrixBlock::vectorProduct(cc0,  eps*arm);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1, -eps*arm);
    
#if 0
    std::cerr.precision(3);
    // check that image D calculated from (AB) is near its target C
    Vector D = aR * ptA.pos1() + bR * ptA.pos2();
    std::cerr << "D: " << std::setw(9) << D << std::endl;
    std::cerr << "C: " << std::setw(9) << ptB.pos() << std::endl;
#endif
    
    MatrixBlock waT = aR.transposed(-weight);
    MatrixBlock wbT = bR.transposed(-weight);
    
    // fill the matrix mC
    add_block_diag(ii0, waT*aR); //this term is symmetric but not diagonal
    add_block(ii1, ii0, wbT*aR);
    add_block_diag(ii1, wbT*bR); //this term is symmetric but not diagonal

    if ( ii2 > ii0 )
    {
        sub_block(ii2, ii0, waT.transposed());
        sub_block(ii2, ii1, wbT.transposed());
    }
    else
    {
        sub_block(ii0, ii2, waT);
        sub_block(ii1, ii2, wbT);
    }
    add_block_diag(ii2, MatrixBlock(0, -weight));
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, waT*off);
            add_base(ii1, wbT*off);
            add_base(ii2, off, weight);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
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
    
    // indices in the matrix mB and mC:
    const size_t ia0 = ptA.matIndex1(),  ii0 = DIM * ia0;
    const size_t ia1 = ptA.matIndex2(),  ii1 = DIM * ia1;
    const size_t ib2 = ptB.matIndex1(),  ii2 = DIM * ib2;
    const size_t ib3 = ptB.matIndex2(),  ii3 = DIM * ib3;
    
    if ( any_equal(ia0, ia1, ib2, ib3) )
        return;
    
    const real ee = arm / ptA.len(), we = weight * ee;

    // coefficients and weights:
    const real cc0 =  ptA.coef0(),  ww0 = weight * cc0;
    const real cc1 =  ptA.coef1(),  ww1 = weight * cc1;
    const real cc2 = -ptB.coef0(),  ww2 = weight * cc2;
    const real cc3 = -ptB.coef1(),  ww3 = weight * cc3;

    sub_iso(ia0, ia0, ww0 * cc0 + we * ee);
    sub_iso(ia1, ia1, ww1 * cc1 + we * ee);
    sub_iso(ib2, ib2, ww2 * cc2);
    sub_iso(ib3, ib2, ww3 * cc2);
    sub_iso(ib3, ib3, ww3 * cc3);

#if ( 1 )
    real wd = ww0 * cc1 - we * ee;
    sub_block(ii1, ii0, Matrix22(wd, -we, we, wd));
#else
    Matrix22 A(cc0, -ee,  ee, cc0);
    Matrix22 B(cc1,  ee, -ee, cc1);
    add_block(ii1, ii0, -weight, B.trans_mul(A));
#endif
    
    if ( ii2 > ii0 )
    {
        Matrix22 A(ww0, -we,  we, ww0);
        Matrix22 B(ww1,  we, -we, ww1);
        add_block(ii2, ii0, -cc2, A);
        add_block(ii3, ii0, -cc3, A);
        add_block(ii2, ii1, -cc2, B);
        add_block(ii3, ii1, -cc3, B);
    }
    else
    {
        Matrix22 At(ww0,  we, -we, ww0);
        Matrix22 Bt(ww1, -we,  we, ww1);
        add_block(ii0, ii2, -cc2, At);
        add_block(ii1, ii2, -cc2, Bt);
        add_block(ii0, ii3, -cc3, At);
        add_block(ii1, ii3, -cc3, Bt);
    }
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            Matrix22 wA(ww0,  we, -we, ww0);
            Matrix22 wB(ww1, -we,  we, ww1);
            sub_base(ii0, wA*off);
            sub_base(ii1, wB*off);
            add_base(ii2, off, -ww2);
            add_base(ii3, off, -ww3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}

#endif


/**
Link `B` to an interpolated point `S` on the side of `A`:

    S = pos(A) + cross( arm, A.dir() )

Where A.dir() is the direction of the Fiber supporting `A`, in `A`.
The vector `arm` should ideally be perpendicular to A.dir(), and in this case,
`A` and `B` are separated by norm(arm). The force is linear:

    force_B = weight * ( S - B )
    force_S = -force_B

The `force_S` is redistributed on the vertices on each side of `A`,
according to the interpolation coefficients, as usual.

This code is valid in any dimension and works in 2 and 3D
*/
void Meca::addSideLink3D(Interpolation const& ptA,
                         Interpolation const& ptB,
                         Torque const& arm,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();

    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    const real eps = -1.0 / ptA.len();

    MatrixBlock aR = MatrixBlock::vectorProduct(cc0,  eps*arm);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1, -eps*arm);

    const real wcc2 = -weight * cc2;
    const real wcc3 = -weight * cc3;

    MatrixBlock waT = aR.transposed(-weight);
    MatrixBlock wbT = bR.transposed(-weight);
    
    // fill the matrix mC
    add_block_diag(ii0, waT*aR);  //this term is symmetric but not diagonal
    add_block(ii1, ii0, wbT*aR);
    add_block_diag(ii1, wbT*bR);  //this term is symmetric but not diagonal
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, wcc2, aR);
        add_block(ii3, ii0, wcc3, aR);
        add_block(ii2, ii1, wcc2, bR);
        add_block(ii3, ii1, wcc3, bR);
    }
    else
    {
        add_block(ii0, ii2, cc2, waT);
        add_block(ii0, ii3, cc3, waT);
        add_block(ii1, ii2, cc2, wbT);
        add_block(ii1, ii3, cc3, wbT);
    }
    add_block_diag(ii2, MatrixBlock(0, wcc2*cc2));
    add_block(ii3, ii2, MatrixBlock(0, wcc3*cc2));
    add_block_diag(ii3, MatrixBlock(0, wcc3*cc3));
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, waT*off);
            add_base(ii1, wbT*off);
            add_base(ii2, off, -wcc2);
            add_base(ii3, off, -wcc3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
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
    mC.reset();
    addSideLink3D(ptA, ptB, arm, weight);
    
    {
        std::ofstream o("x");
        o << "testSideLink " << ptA.matIndex1() << " " << ptB.matIndex() << "\n";
        mC.printSparse(o, REAL_EPSILON);
        o.close();
    }
    
    mC.reset();
    size_t P = ptB.point();
    if ( P > 0 )
        addSideLink3D(ptA, Interpolation(ptB.mecable(), P-1, P, 1), arm, weight);
    else
        addSideLink3D(ptA, Interpolation(ptB.mecable(), P, P+1, 0), arm, weight);
    
    {
        std::ofstream o("y");
        o << "testSideLink " << ptA.matIndex1() << " " << ptB.matIndex() << "\n";
        mC.printSparse(o, REAL_EPSILON);
        o.close();
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
#pragma mark - Symmetric off-axis links
//------------------------------------------------------------------------------

#if ( DIM == 2 )

// this is old style code, addressing mC() directly, but it works!
void Meca::addSideSideLink2Dalt(Interpolation const& ptA,
                                Interpolation const& ptB,
                                const real armA, const real armB,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mB:
    size_t ia1 = ptA.matIndex1(), ia2 = ptA.matIndex2();
    size_t ib1 = ptB.matIndex1(), ib2 = ptB.matIndex2();
    
    if ( any_equal(ia1, ia2, ib1, ib2) )
        return;

    const real ca1 =  ptA.coef0(), ca2 =  ptA.coef1();
    const real cb1 = -ptB.coef0(), cb2 = -ptB.coef1();
    
    const real ee1 = armA / ptA.len();
    const real ee2 = armB / ptB.len();
    
    const real W = -weight;
    const real ca1w = ca1 * W, ca2w = ca2 * W;
    const real cb1w = cb1 * W, cb2w = cb2 * W;
   
    const real ee1w = ee1 * W, ee1ee1w = ee1 * ee1w;
    const real ee2w = ee2 * W, ee2ee2w = ee2 * ee2w;
    const real ee1ee2w = ee1 * ee2w;
    
    //we put the isotropic terms in mB
    add_iso(ia1, ia1, ca1w * ca1 + ee1ee1w);
    add_iso(ia2, ia1, ca1w * ca2 - ee1ee1w);
    add_iso(ia2, ia2, ca2w * ca2 + ee1ee1w);
    
    add_iso(ib1, ib1, cb1w * cb1 + ee2ee2w);
    add_iso(ib2, ib1, cb1w * cb2 - ee2ee2w);
    add_iso(ib2, ib2, cb2w * cb2 + ee2ee2w);
    
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
        gle::drawLink(ptA.pos(), cross(armA, ptA.dir()), cross(armB, ptB.dir()), ptB.pos());
    }
#endif
    
    if ( modulo )
        throw Exception("addSideSideLink2D is not usable with periodic boundary conditions");
}


// new style code, validated on 19.01.2020
void Meca::addSideSideLink2D(Interpolation const& ptA,
                             Interpolation const& ptB,
                             const real armA, const real armB,
                             const real weight)
{
    assert_true( weight >= 0 );

    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
 
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    // weights and coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    const real ee1 =  armA / ptA.len();
    const real ee2 = -armB / ptB.len();

    Matrix22 A(cc0, -ee1,  ee1, cc0);
    Matrix22 B(cc1,  ee1, -ee1, cc1);
    Matrix22 C(cc2, -ee2,  ee2, cc2);
    Matrix22 D(cc3,  ee2, -ee2, cc3);

    Matrix22 wAt = A.transposed(-weight); //(ww0,  we1, -we1, ww0);
    Matrix22 wBt = B.transposed(-weight); //(ww1, -we1,  we1, ww1);
    Matrix22 wCt = C.transposed(-weight); //(ww2,  we2, -we2, ww2);
    Matrix22 wDt = D.transposed(-weight); //(ww3, -we2,  we2, ww3);
        
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

    if ( ii2 > ii0 )
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
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, wAt*off);
            add_base(ii1, wBt*off);
            add_base(ii2, wCt*off);
            add_base(ii3, wDt*off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(armA, ptA.dif()), cross(armB, ptB.dir()), ptB.pos());
    }
#endif
}


#endif

// new style code, validated on 19.01.2020
void Meca::addSideSideLink3D(Interpolation const& ptA,
                             Interpolation const& ptB,
                             Torque const& armA, Torque const& armB,
                             const real weight)
{
    assert_true( weight >= 0 );

    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
 
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    
    const real epsA = -1.0 / ptA.len();
    const real epsB =  1.0 / ptA.len();

    MatrixBlock A = MatrixBlock::vectorProduct(cc0,  epsA*armA);
    MatrixBlock B = MatrixBlock::vectorProduct(cc1, -epsA*armA);
    MatrixBlock C = MatrixBlock::vectorProduct(cc2,  epsB*armB);
    MatrixBlock D = MatrixBlock::vectorProduct(cc3, -epsB*armB);

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

    if ( ii2 > ii0 )
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
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, wAt*off);
            add_base(ii1, wBt*off);
            add_base(ii2, wCt*off);
            add_base(ii3, wDt*off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(armA, ptA.dir()), cross(armB, ptB.dir()), ptB.pos());
    }
#endif
}


/**
 Link `ptA` (A) and `ptB` (B),
 but the links are maded between SA and SB which are located
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
    
    Vector dir = ptB.pos() - ptA.pos();
    real armA = std::copysign(0.5*len, cross(ptA.diff(), dir));
    real armB = std::copysign(0.5*len, cross(dir, ptB.diff()));
    addSideSideLink2D(ptA, ptB, armA, armB, weight);

#else
    
    Vector dir = ptB.pos() - ptA.pos();
    Vector armA = cross(ptA.diff(), dir).normalized(0.5*len);
    Vector armB = cross(dir, ptB.diff()).normalized(0.5*len);
    addSideSideLink3D(ptA, ptB, armA, armB, weight);

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
    
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    //force coefficients on the points:
    const real A = ptA.coef0();
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

    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_block_diag(ii0, AA, wT);
    add_block_diag(ii1, BB, wT);
    add_block_diag(ii2, wT);

    add_block(ii0, ii1, AB, wT);
    add_block(ii0, ii2, -A, wT);
    add_block(ii1, ii2, -B, wT);
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, A*off);
            add_base(ii1, B*off);
            sub_base(ii2, off);
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
 
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    // coefficients:
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
    
    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block_diag(ii2, cc2*cc2, wT);
    add_block(ii3, ii2, cc3*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            add_base(ii2, off, cc2);
            add_base(ii3, off, cc3);
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
#pragma mark - Off-axis frictionless links (used for steric interactions)

#if ( DIM == 2 )

void Meca::addSideSlidingLink2D(Interpolation const& ptA,
                                Mecapoint const& ptB,
                                const real arm,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    /*
     The length of fiber's segment is known, so we could spare the calculation
     of `seg` below, which involves a square root
     */
    const real seg = ptA.len();
    const real ee = arm / seg;
    Vector dir = ptA.diff() / seg;
    
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    
    // the projection matrix: wP = -weight * [ I - dir (x) dir ]
    MatrixBlock wP = MatrixBlock::offsetOuterProduct(-weight, dir, weight);
    
    // anti-symmetric matrix blocks:
    const Matrix22 A(cc0, -ee,  ee, cc0);
    const Matrix22 B(cc1,  ee, -ee, cc1);

    /*
     We use block operations to set the matrix block by block:
     | A'PA  A'PB  A'P |
     | B'PA  B'PB  B'P |
     |   PA    PB    P |
     This matrix has symmetric and anti-symmetric blocks,
     since P' = P whereas A and B are anti-symmetric
     */
    assert_true( ii0 < ii1 );
    
    if ( ii2 > ii0 )
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
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            //off = -weight * ( off - dot(off, dir) * dir );
            off = wP * off;
            add_base(ii0, A.trans_vecmul(off));
            add_base(ii1, B.trans_vecmul(off));
            sub_base(ii2, off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}


/**
 Alternative 2D method in which we add an offset to vBAS
 */
void Meca::addSideSlidingLinkS(Interpolation const& ptA,
                               Mecapoint const& ptB,
                               const real arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    // indices
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    // set vector 'axi' perpendicular to Fiber:
    const Vector2 axi = cross(1.0, ptA.dir());
    
    // P = -weight * [ I - dir (x) dir ]
    Matrix22 wT = Matrix22::outerProduct(axi, -weight);
    
    // coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    
    real wa = -weight * arm;
    add_base(ii0, axi, wa*cc0);
    add_base(ii1, axi, wa*cc1);
    add_base(ii2, axi, wa);
    
    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block(ii2, ii0,    -cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii2, ii1,    -cc1, wT);
    add_block_diag(ii2, wT);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            sub_base(ii2, off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), axi, ptB.pos());
    }
#endif
}

#elif ( DIM == 3 )
/**
 Older code
 Vector 'arm' must be parallel to the link and orthogonal to 'ptA'
 */
void Meca::addSideSlidingLinkS(Interpolation const& ptA,
                               Mecapoint const& ptB,
                               Vector3 const& arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex();

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
    
    Vector3 warm = -weight * arm;
    MatrixBlock wT = MatrixBlock::outerProduct(arm, -weight/arm.normSqr());
    
    // coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();

    add_base(ii0, warm, cc0);
    add_base(ii1, warm, cc1);
    sub_base(ii2, warm);

    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block(ii2, ii0,    -cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii2, ii1,    -cc1, wT);
    add_block_diag(ii2, wT);

    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            sub_base(ii2, off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), arm, ptB.pos());
    }
#endif
}

#endif



void Meca::addSideSlidingLink3D(Interpolation const& ptA,
                                Mecapoint const& ptB,
                                Torque const& arm,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex();
    
    if ( any_equal(ii0, ii1, ii2) )
        return;
    
    // coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    const real eps = -1.0 / ptA.len();
    
    MatrixBlock aR = MatrixBlock::vectorProduct(cc0,  eps*arm);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1, -eps*arm);

    // the projection matrix: P = -weight * [ I - dir (x) dir ]
    MatrixBlock wP = MatrixBlock::offsetOuterProduct(-weight, ptA.diff()*eps, weight);
    
    MatrixBlock aTwP = aR.trans_mul(wP);
    MatrixBlock bTwP = bR.trans_mul(wP);
    
    // fill the matrix mC
    add_block_diag(ii0, aTwP*aR);
    add_block(ii1, ii0, bTwP*aR);
    add_block_diag(ii1, bTwP*bR);
    if ( ii2 > ii0 )
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
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, aTwP*off);
            add_base(ii1, bTwP*off);
            sub_base(ii2, wP*off);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
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
void Meca::addSideSlidingLink(Interpolation const& ptA,
                              Mecapoint const& ptB,
                              const real len,
                              const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    Vector ab = ptB.pos()-ptA.pos();
    if ( modulo )
        modulo->fold(ab);
    real arm = std::copysign(len, cross(ptA.diff(), ab));
    addSideSlidingLink2D(ptA, ptB, arm, weight);
    
#else
    /*
    // old version used before 2020 for steric interactions
    Vector arm = calculateArm(ptB.pos()-ptA.pos1(), ptA.diff(), len);
    addSideSlidingLinkS(ptA, ptB, arm, weight);
    */
    //Vector arm = cross(ptA.diff(), ptB.pos()-ptA.pos1());
    Vector arm = ptB.pos()-ptA.pos1();
    if ( modulo )
        modulo->fold(arm);
    arm = cross(ptA.diff(), arm);
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSideSlidingLink3D(ptA, ptB, arm*(len/n), weight);

#endif
}


#pragma mark - More off-axis frictionless links


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
*/
void Meca::addSideSlidingLink2D(Interpolation const& ptA,
                                Interpolation const& ptB,
                                const real arm,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    const real eps = -1.0 / ptA.len();
    
    real ee = arm * eps;
    
    Matrix22 aR(cc0, ee,-ee, cc0);  // aR = alpha - len * R
    Matrix22 bR(cc1,-ee, ee, cc1);  // bR = beta + len * R
    
    // the projection matrix: P = -weight * [ I - dir (x) dir ]
    Matrix22 wP = Matrix22::offsetOuterProduct(-weight, ptA.diff()*eps, weight);
    
    Matrix22 aTwP = aR.trans_mul(wP);
    Matrix22 bTwP = bR.trans_mul(wP);

    // fill the matrix mC
    add_block_diag(ii0, aTwP*aR);
    add_block(ii1, ii0, bTwP*aR);
    add_block_diag(ii1, bTwP*bR);
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, cc2, aTwP.transposed());
        add_block(ii3, ii0, cc3, aTwP.transposed());
        add_block(ii2, ii1, cc2, bTwP.transposed());
        add_block(ii3, ii1, cc3, bTwP.transposed());
    }
    else
    {
        add_block(ii0, ii2, cc2, aTwP);
        add_block(ii0, ii3, cc3, aTwP);
        add_block(ii1, ii2, cc2, bTwP);
        add_block(ii1, ii3, cc3, bTwP);
    }
    add_block_diag(ii2, cc2*cc2, wP);
    add_block(ii3, ii2, cc3*cc2, wP);
    add_block_diag(ii3, cc3*cc3, wP);

    if ( modulo )
    {
        throw Exception("addSideSlidingLink2D untested with periodic boundary conditions");
        
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, aTwP*off);
            add_base(ii1, bTwP*off);
            add_base(ii2, wP*off, cc2);
            add_base(ii3, wP*off, cc3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}


/**
 Alternative old 2D method in which we add an offset to vBAS
 */
void Meca::addSideSlidingLinkS(Interpolation const& ptA,
                               Interpolation const& ptB,
                               const real arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;

    // set vector 'axi' perpendicular to Fiber:
    Vector2 axi = cross(1.0, ptA.dir());
    
    Matrix22 wT = Matrix22::outerProduct(axi, -weight);

    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    
    Vector warm = axi * ( -weight / arm );
    
    add_base(ii0, warm, cc0);
    add_base(ii1, warm, cc1);
    add_base(ii2, warm, cc2);
    add_base(ii3, warm, cc3);
    
    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block_diag(ii2, cc2*cc2, wT);
    add_block(ii3, ii2, cc3*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            add_base(ii2, off, cc2);
            add_base(ii3, off, cc3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), axi, ptB.pos());
    }
#endif
}


#elif ( DIM >= 3 )

/**
 Older code
 Vector 'arm' must be parallel to the link and orthogonal to 'ptA'
 */
void Meca::addSideSlidingLinkS(Interpolation const& ptA,
                               Interpolation const& ptB,
                               Vector3 const& arm,
                               const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
   
    Vector3 warm = -weight * arm;
    MatrixBlock wT = MatrixBlock::outerProduct(arm, -weight/arm.normSqr());

    add_base(ii0, warm, cc0);
    add_base(ii1, warm, cc1);
    add_base(ii2, warm, cc2);
    add_base(ii3, warm, cc3);
    
    // fill the matrix mC
    add_block_diag(ii0, cc0*cc0, wT);
    add_block(ii1, ii0, cc1*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
    }
    else
    {
        add_block(ii0, ii2, cc2*cc0, wT);
        add_block(ii0, ii3, cc3*cc0, wT);
        add_block(ii1, ii2, cc2*cc1, wT);
        add_block(ii1, ii3, cc3*cc1, wT);
    }
    add_block_diag(ii2, cc2*cc2, wT);
    add_block(ii3, ii2, cc3*cc2, wT);
    add_block_diag(ii3, cc3*cc3, wT);
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            off = wT * off;
            add_base(ii0, off, cc0);
            add_base(ii1, off, cc1);
            add_base(ii2, off, cc2);
            add_base(ii3, off, cc3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
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

void Meca::addSideSlidingLink3D(Interpolation const& ptA,
                                Interpolation const& ptB,
                                Torque const& arm,
                                const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    const size_t ii2 = DIM * ptB.matIndex1();
    const size_t ii3 = DIM * ptB.matIndex2();
    
    if ( any_equal(ii0, ii1, ii2, ii3) )
        return;
    
    // coefficients:
    const real cc0 =  ptA.coef0();
    const real cc1 =  ptA.coef1();
    const real cc2 = -ptB.coef0();
    const real cc3 = -ptB.coef1();
    const real eps = -1.0 / ptA.len();
    
    MatrixBlock aR = MatrixBlock::vectorProduct(cc0,  eps*arm);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1, -eps*arm);

    // the projection matrix: P = -weight * [ I - dir (x) dir ]
    MatrixBlock wP = MatrixBlock::offsetOuterProduct(-weight, ptA.diff()*eps, weight);

    MatrixBlock aTwP = aR.trans_mul(wP);
    MatrixBlock bTwP = bR.trans_mul(wP);
    
    // fill the matrix mC
    add_block_diag(ii0, aTwP*aR);
    add_block(ii1, ii0, bTwP*aR);
    add_block_diag(ii1, bTwP*bR);
    if ( ii2 > ii0 )
    {
        add_block(ii2, ii0, cc2, aTwP.transposed());
        add_block(ii3, ii0, cc3, aTwP.transposed());
        add_block(ii2, ii1, cc2, bTwP.transposed());
        add_block(ii3, ii1, cc3, bTwP.transposed());
    }
    else
    {
        add_block(ii0, ii2, cc2, aTwP);
        add_block(ii0, ii3, cc3, aTwP);
        add_block(ii1, ii2, cc2, bTwP);
        add_block(ii1, ii3, cc3, bTwP);
    }
    add_block_diag(ii2, cc2*cc2, wP);
    add_block(ii3, ii2, cc3*cc2, wP);
    add_block_diag(ii3, cc3*cc3, wP);
    
    if ( modulo )
    {
        Vector off = modulo->offset( ptB.pos() - ptA.pos() );
        if ( off.is_not_zero() )
        {
            add_base(ii0, aTwP*off);
            add_base(ii1, bTwP*off);
            add_base(ii2, wP*off, cc2);
            add_base(ii3, wP*off, cc3);
        }
    }
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        drawLinkM(ptA.pos(), cross(arm, ptA.dir()), ptB.pos());
    }
#endif
}


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

void Meca::addSideSlidingLink(Interpolation const& ptA,
                              Interpolation const& ptB,
                              const real len,
                              const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSideSlidingLink is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    real arm = std::copysign(len, cross(ptA.diff(), ptB.pos()-ptA.pos1()));
    addSideSlidingLink2D(ptA, ptB, arm, weight);
    
#else
    
    // old version for steric interactions
    //Vector arm = calculateArm(ptB.pos()-ptA.pos1(), ptA.diff(), len);
    //addSideSlidingLinkS(ptA, ptB, arm, weight);
    //Vector arm = cross(ptA.diff(), ptB.pos()-ptA.pos1());
    Vector arm = ptB.pos()-ptA.pos1();
    if ( modulo )
        modulo->fold(arm);
    arm = cross(ptA.diff(), arm);
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSideSlidingLink3D(ptA, ptB, arm*(len/n), weight);

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
    const size_t inx = ptA.matIndex();
    
    sub_iso(inx, inx, weight);
    
    if ( modulo )
        pos = modulo->image(pos, ptA.pos());

    add_base(DIM*inx, pos, weight);
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), pos);
    }
#endif
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
    
    const size_t ii0 = pti.matIndex1();
    const size_t ii1 = pti.matIndex2();
    
    const real c1 = pti.coef0();
    const real c2 = pti.coef1();
    
    assert_true( 0 <= c1  &&  c1 <= 1 );
    assert_true( 0 <= c2  &&  c2 <= 1 );

    const real c2w = weight * c2;
    const real c1w = weight * c1;
    
    sub_iso(ii0, ii0, c1w * c1);
    sub_iso(ii0, ii1, c2w * c1);
    sub_iso(ii1, ii1, c2w * c2);
    
    if ( modulo )
        pos = modulo->image(pos, pti.pos());
    
    add_base(DIM*ii0, pos, c1w);
    add_base(DIM*ii1, pos, c2w);
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(pti.mecable()->signature()).load();
        gle::drawLink(pti.pos(), pos);
    }
#endif
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
    const size_t inx = DIM * ptA.matIndex();

#if ( DIM == 2 )
    sub_iso(ptA.matIndex(), ptA.matIndex(), weight);
#elif ( DIM > 2 )
    mC(inx,   inx  ) -= weight;
    mC(inx+1, inx+1) -= weight;
    assert_true( pos.ZZ == 0 );
#endif
    add_base(inx, pos, weight);
}


void Meca::addPointClampToAll(Vector const& pos, const real weight)
{
    Vector vec = weight * pos;
    for ( size_t p = 0; p < nb_points(); ++p )
    {
        sub_iso(p, p, weight);
        add_base(p, vec);
    }
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
    
    const size_t inx = DIM * ptA.matIndex();
    
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
        
        add_block_diag(inx, wT);
        add_base(inx, wla*off-wT*center);
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
        
        const size_t ii0 = DIM * ptA.matIndex1();
        const size_t ii1 = DIM * ptA.matIndex2();
        
        MatrixBlock wT;
        /* To stabilize the matrix with compression, we remove negative eigenvalues
         This is done by using len = 1 in the formula for links that are shorter
         than the desired target. */
        if ( rad < len )
            wT = MatrixBlock::offsetOuterProduct(wla-weight, off/len, -wla);
        else
            wT = MatrixBlock::outerProduct(off/len, -weight);

        // coefficients:
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
        center = modulo->image(center, pos);
    addSphereClamp(pos-center, pte, center, rad, weight);
}


void Meca::addSphereClamp(Interpolation const& pti,
                          Vector center,
                          real rad,
                          const real weight)
{
    Vector pos = pti.pos();
    if ( modulo )
        center = modulo->image(center, pos);
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
void Meca::addCylinderClampX(Mecapoint const& pte,
                             real rad, const real weight)
{
    assert_true( weight >= 0 );
    const size_t inx = DIM * pte.matIndex();
    
#if ( DIM == 2 )
    
    mC(inx+1, inx+1) -= weight;
    vBAS[inx+1]      += weight * std::copysign(rad, pte.pos().YY);
    
#elif ( DIM >= 3 )

    Vector pos = pte.pos();
    real fac = weight * rad;
    real len = pos.normYZ();
    if ( len < REAL_EPSILON )
        return;
    
    Vector dir(0, pos.YY/len, pos.ZZ/len);
    
    if ( rad < len )
    {
        real wla = weight * rad / len;
        mC(inx+1, inx+1) -= wla * dir.YY * dir.YY + weight - wla;
        mC(inx+2, inx+1) -= wla * dir.YY * dir.ZZ;
        mC(inx+2, inx+2) -= wla * dir.ZZ * dir.ZZ + weight - wla;
    }
    else
    {
        fac = weight * rad;

        mC(inx+1, inx+1) -= weight * dir.YY * dir.YY;
        mC(inx+2, inx+1) -= weight * dir.YY * dir.ZZ;
        mC(inx+2, inx+2) -= weight * dir.ZZ * dir.ZZ;
    }
    
    // there should be no XX component here!
    vBAS[inx+1] += fac * dir.YY;
    vBAS[inx+2] += fac * dir.ZZ;
    
#endif
}


/**
 Link `pte` (P) to a cylinder of axis Y and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, P.YY, 0)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the XZ plane.
 */
void Meca::addCylinderClampY(Mecapoint const& pte,
                             real rad, const real weight)
{
    assert_true( weight >= 0 );
    const size_t inx = DIM * pte.matIndex();
    
#if ( DIM == 2 )
    
    mC(inx, inx) -= weight;
    vBAS[inx]    += weight * std::copysign(rad, pte.pos().XX);
    
#elif ( DIM >= 3 )

    Vector pos = pte.pos();
    real fac = weight * rad;
    real len = pos.normXZ();
    if ( len < REAL_EPSILON )
        return;
    
    Vector dir(pos.XX/len, 0, pos.ZZ/len);
    
    if ( rad < len )
    {
        real wla = weight * rad / len;
        mC(inx  , inx  ) -= wla * dir.XX * dir.XX + weight - wla;
        mC(inx+2, inx  ) -= wla * dir.XX * dir.ZZ;
        mC(inx+2, inx+2) -= wla * dir.ZZ * dir.ZZ + weight - wla;
    }
    else
    {
        mC(inx  , inx  ) -= weight * dir.XX * dir.XX;
        mC(inx+2, inx  ) -= weight * dir.XX * dir.ZZ;
        mC(inx+2, inx+2) -= weight * dir.ZZ * dir.ZZ;
    }
    
    vBAS[inx  ] += fac * dir.XX;
    // there should be no YY component here!
    vBAS[inx+2] += fac * dir.ZZ;
    
#endif
}


/**
 Link `pte` (P) to a cylinder of axis Z and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, 0, P.ZZ)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the XY plane.
 */
void Meca::addCylinderClampZ(Mecapoint const& pte,
                             const real rad, const real weight)
{
    assert_true( weight >= 0 );
    
#if ( DIM > 1 )

    const size_t inx = DIM * pte.matIndex();
    Vector pos = pte.pos();
    real fac = weight * rad;
    real len = pos.normXY();
    if ( len < REAL_EPSILON )
        return;
    
    Vector dir(pos.XX/len, pos.YY/len, 0);
    
    if ( rad < len )
    {
        real wla = weight * rad / len;
        mC(inx  , inx  ) -= wla * dir.XX * dir.XX + weight - wla;
        mC(inx+1, inx  ) -= wla * dir.XX * dir.YY;
        mC(inx+1, inx+1) -= wla * dir.YY * dir.YY + weight - wla;
    }
    else
    {
        mC(inx  , inx  ) -= weight * dir.XX * dir.XX;
        mC(inx+1, inx  ) -= weight * dir.XX * dir.YY;
        mC(inx+1, inx+1) -= weight * dir.YY * dir.YY;
    }
    
    vBAS[inx  ] += fac * dir.XX;
    vBAS[inx+1] += fac * dir.YY;
    // there should be no ZZ component here!

#endif
}

/**
 Link `pte` (P) to a cylinder of axis Z and radius `rad`
 The force is affine with non-zero resting length:
 
     G = Vector(0, 0, P.ZZ)
     force = weight * ( G - P ) * ( rad / |PG| - 1 )
 
 The force resides in the XY plane.
 */
void Meca::addCylinderClamp(Mecapoint const& pte,
                            Vector const& axis, Vector const& center,
                            const real rad, const real weight)
{
#if ( DIM > 2 )
    
    assert_true( weight >= 0 );
    const size_t inx = DIM * pte.matIndex();

    //Projection matrix along the cylinder axis:  P = [ I - axis (x) axis ]
    MatrixBlock P = MatrixBlock::offsetOuterProduct(1.0, axis, -1.0/axis.normSqr());
    
    Vector dir = P * ( pte.pos() - center );
    real len = dir.norm();
    
    if ( len > REAL_EPSILON )
    {
        real wla = weight * rad / len;
        
        if ( rad < len )
            P = P * MatrixBlock::offsetOuterProduct(wla-weight, dir/len, -wla);
        else
            P = P * MatrixBlock::outerProduct(dir/len, -weight);
        
        add_block_diag(inx, P);
        add_base(inx, wla*dir - P*center);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Off-axis links to fixed positions
//------------------------------------------------------------------------------


#if ( DIM == 2 )

void Meca::addSidePointClamp2D(Interpolation const& ptA,
                               Vector pos,
                               const real arm,
                               const real weight)
{
    //force coefficients on the points:
    const real A = ptA.coef0(),  wA = weight * A;
    const real B = ptA.coef1(),  wB = weight * B;
    
    const real E = arm / ptA.len();
    const real wE = weight * E;
    const real wEE = weight * E * E;
    
    //index in the matrix mB:
    size_t ii0 = ptA.matIndex1();
    size_t ii1 = ptA.matIndex2();
    
    //we put the isotropic terms in mB
    sub_iso(ii0, ii0,  wA * A + wEE);
    sub_iso(ii0, ii1,  wA * B - wEE);
    sub_iso(ii1, ii1,  wB * B + wEE);
    
    //index in the matrix mC:
    ii0 *= DIM;
    ii1 *= DIM;
    
    mC(ii0  , ii1+1) += wE;
    mC(ii0+1, ii1  ) -= wE;

    if ( modulo )
        pos += modulo->offset( ptA.pos() - pos );
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), pos);
    }
#endif

    //it seems to works also fine without the term in eew* below:
    Vector off = wE * Vector(-pos.YY, pos.XX);

    add_base(ii0, wA*pos+off);
    add_base(ii1, wB*pos-off);
}

#endif

/**
 A link of stiffness `weight`, between offset_point on the side of `ptA`, and the fixed position `pos` (G).

 This uses the vector product x -> cross(arm, x) to offset the point on which the link is attached:
 
     offset_point = fiber_point + cross(arm, fiber_dir)
 
 with fiber_point = ptA.pos() and fiber_dir = ptA.diff().normalized.
 `arm` must be perpendicular to link ( G - ptA.pos() )

 F. Nedelec, March 2011
 
 */
void Meca::addSidePointClamp3D(Interpolation const& ptA,
                               Vector pos,
                               Torque const& arm,
                               real const weight)
{
    assert_true( weight >= 0 );
    // indices in the matrix mC:
    const size_t ii0 = DIM * ptA.matIndex1();
    const size_t ii1 = DIM * ptA.matIndex2();
    
    // coefficients:
    const real cc0 = ptA.coef0();
    const real cc1 = ptA.coef1();
    const real eps = -1.0 / ptA.len();

    MatrixBlock aR = MatrixBlock::vectorProduct(cc0,  eps*arm);
    MatrixBlock bR = MatrixBlock::vectorProduct(cc1, -eps*arm);

    MatrixBlock waT = aR.transposed(-weight);
    MatrixBlock wbT = bR.transposed(-weight);
    
    // fill the matrix mC
    add_block_diag(ii0, waT*aR); //this term is symmetric but not diagonal
    add_block(ii1, ii0, wbT*aR);
    add_block_diag(ii1, wbT*bR); //this term is symmetric but not diagonal

    if ( modulo )
        pos += modulo->offset( ptA.pos() - pos );
    
#if DRAW_MECA_LINKS
    if ( drawLinks )
    {
        gle::bright_color(ptA.mecable()->signature()).load();
        gle::drawLink(ptA.pos(), cross(arm, ptA.dir()), pos);
    }
#endif

    sub_base(ii0, waT*pos);
    sub_base(ii1, wbT*pos);
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
                             Vector const& pos,
                             const real len,
                             const real weight)
{
#if ( DIM == 1 )
    
    throw Exception("Meca::addSidePointClamp is meaningless in 1D");
    
#elif ( DIM == 2 )
    
    // 'arm' is a vector in the Z direction
    real arm = std::copysign(len, cross(ptA.diff(), pos-ptA.pos1()));
    addSidePointClamp2D(ptA, pos, arm, weight);
   
#else
    
    // 'arm' perpendicular to link and fiber is obtained by vector product:
    //Vector arm = cross(ptA.diff(), pos-ptA.pos1());
    Vector arm = pos-ptA.pos1();
    if ( modulo )
        modulo->fold(arm);
    arm = cross(ptA.diff(), arm);
    real n = arm.norm();
    if ( n > REAL_EPSILON )
        addSidePointClamp3D(ptA, pos, arm*(len/n), weight);

#endif  
}

//------------------------------------------------------------------------------
#pragma mark - Links to fixed lines and planes
//------------------------------------------------------------------------------

/**
 This constrains a single degree of freedom indicated by 'inx', representing a
 planar constraint in 3D. Hence 'pos' is the X, Y or Z component of the plane.
*/
void Meca::addPlaneClamp(size_t inx,
                         real off,
                         const real weight)
{
    mC(inx, inx) -= weight;
    vBAS[inx] += weight * off;
}

void Meca::addPlaneClampX(Mecapoint const& P, real off, real weight)
{
    size_t inx = DIM * P.matIndex();
    mC(inx, inx) -= weight;
    vBAS[inx] += weight * off;
}

void Meca::addPlaneClampY(Mecapoint const& P, real off, real weight)
{
    size_t inx = DIM * P.matIndex() + 1;
    mC(inx, inx) -= weight;
    vBAS[inx] += weight * off;
}

void Meca::addPlaneClampZ(Mecapoint const& P, real off, real weight)
{
    size_t inx = DIM * P.matIndex() + 2;
    mC(inx, inx) -= weight;
    vBAS[inx] += weight * off;
}

/**
 Link `ptA` (X) to the line defined by `G` and tangent vector `dir`.
 The force is linear with position with a stiffness `weight`, and its
 components parallel to `dir` are removed corresponding to a frictionless line:
 
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
    
    const size_t inx = DIM * pte.matIndex();
    
    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_block_diag(inx, wT);
    sub_base(inx, wT*pos);
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
    
    //index in the matrix mC:
    const size_t ii0 = DIM * pti.matIndex1();
    const size_t ii1 = DIM * pti.matIndex2();

    //force coefficients on the points:
    const real cc0 = pti.coef0();
    const real cc1 = pti.coef1();

    // wT = -weight * [ I - dir (x) dir ]
    MatrixBlock wT = MatrixBlock::offsetOuterProduct(-weight, dir, weight);

    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii0, ii1, cc0*cc1, wT);
    
    //add the constant term:
    add_base(ii0, wT*pos, -cc0);
    add_base(ii1, wT*pos, -cc1);
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

void Meca::addPlaneClamp(Mecapoint const& pte,
                         Vector const& pos,
                         Vector const& dir,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    const size_t inx = DIM * pte.matIndex();
    
    // vBAS[inx] += dir * ( weigth * dot(pos,dir) );
    add_base(inx, dir, weight*dot(pos, dir));
    
#if ( DIM == 1 ) && USE_ISO_MATRIX
    mB(inx, inx) -= weight;
#else
    // wT = -weight * [ dir (x) dir ]
    MatrixBlock wT = MatrixBlock::outerProduct(dir, -weight);
    add_block_diag(inx, wT);
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

void Meca::addPlaneClamp(Interpolation const& pti,
                         Vector const& pos,
                         Vector const& dir,
                         const real weight)
{
    assert_true( weight >= 0 );
    
    //index in the matrix mC:
    const size_t ii0 = DIM * pti.matIndex1();
    const size_t ii1 = DIM * pti.matIndex2();

    //force coefficients on the points:
    const real cc0 = pti.coef0();
    const real cc1 = pti.coef1();
    
    // wT = -weight * [ dir (x) dir ]
    MatrixBlock wT = MatrixBlock::outerProduct(dir, -weight);
    
    add_block_diag(ii0, cc0*cc0, wT);
    add_block_diag(ii1, cc1*cc1, wT);
    add_block(ii0, ii1, cc0*cc1, wT);
    
    //add the constant term:
    add_base(ii0, wT*pos, -cc0);
    add_base(ii1, wT*pos, -cc1);
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
void Meca::addCoulomb(Mecapoint const& ptA, Mecapoint const& ptB, real weight)
{
    Vector ab = ptB.pos() - ptA.pos();
    real abnSqr = ab.normSqr(), abn=sqrt(abnSqr);
    
    const size_t inxA = DIM * ptA.matIndex();
    const size_t inxB = DIM * ptB.matIndex();
    
    if ( abn < REAL_EPSILON ) return;
    ab /= abn;
    
    real abn3 = weight / abnSqr;
    real abn5 = weight / ( abnSqr * abn );
    
    add_base(inxA, ab,-3*abn3);
    add_base(inxB, ab, 3*abn3);
    
    for ( unsigned ii = 0; ii < DIM; ++ii )
    {
        for ( unsigned jj = ii; jj < DIM; ++jj )
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

