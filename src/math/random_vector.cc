// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "messages.h"
#include "random_vector.h"
#include "vector1.h"
#include "vector4.h"
#include "random.h"

/**
 Random vectors are generated using the global Random Generator `RNG`
 */


//------------------------------------------------------------------------------
#pragma mark - 1D Vectors

Vector1 Vector1::randS()       { return Vector1(  RNG.sreal()); }
Vector1 Vector1::randH()       { return Vector1(  RNG.shalf()); }
Vector1 Vector1::randS(real n) { return Vector1(n*RNG.sreal()); }
Vector1 Vector1::randP()       { return Vector1(  RNG.preal()); }
Vector1 Vector1::randP(real n) { return Vector1(n*RNG.preal()); }
Vector1 Vector1::randU()       { return Vector1(  RNG.sflip()); }
Vector1 Vector1::randU(real n) { return Vector1(n*RNG.sflip()); }
void Vector1::addRand(real n) { XX += n*RNG.sreal(); }

Vector1 Vector1::randB()       { return Vector1(  RNG.sreal()); }
Vector1 Vector1::randB(real n) { return Vector1(n*RNG.sreal()); }
Vector1 Vector1::randG(real n) { return Vector1(n*RNG.gauss()); }

Vector1 Vector1::randOrthoU(real) const { return Vector1(0.0); }
Vector1 Vector1::randOrthoB(real) const { return Vector1(0.0); }


//------------------------------------------------------------------------------
#pragma mark - 2D Vectors

Vector2 Vector2::randS()       { return Vector2(  RNG.sreal(),   RNG.sreal()); }
Vector2 Vector2::randH()       { return Vector2(  RNG.shalf(),   RNG.shalf()); }
Vector2 Vector2::randS(real n) { return Vector2(n*RNG.sreal(), n*RNG.sreal()); }
Vector2 Vector2::randP()       { return Vector2(  RNG.preal(),   RNG.preal()); }
Vector2 Vector2::randP(real n) { return Vector2(n*RNG.preal(), n*RNG.preal()); }
Vector2 Vector2::randG(real n) { return Vector2(n*RNG.gauss(), n*RNG.gauss()); }
void Vector2::addRand(real n) { XX += n*RNG.sreal(); YY += n*RNG.sreal(); }

#if ( 0 )

Vector2 Vector2::randU()
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while ( d > 1.0  ||  d < 0.01 );
    return Vector2(x, y) / std::sqrt(d);
}

Vector2 Vector2::randU(const real n)
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while ( d > 1.0  ||  d < 0.01 );
    return Vector2(x, y) * (n/std::sqrt(d));
}

#else

Vector2 Vector2::randU()
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while (( d > 1.0 )|( d == 0 ));
    // avoiding the square-root!
    return Vector2(x*x-y*y, 2*x*y) / d;
}

Vector2 Vector2::randU(const real n)
{
    real d, x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = x*x + y*y;
    } while (( d > 1.0 )|( d == 0 ));
    return Vector2(x*x-y*y, 2*x*y) * (n/d);
}

#endif


Vector2 Vector2::randB()
{
    real x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
    } while ( x*x + y*y > 1.0 );
    return Vector2(x, y);
}


Vector2 Vector2::randB(const real n)
{
    real x, y;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
    } while ( x*x + y*y > 1.0 );
    return Vector2(x, y) * n;
}


Vector2 Vector2::randOrthoU(const real len) const
{
    return Vector2(-YY, XX) * RNG.sflip(len);
}

/** this assumes norm(*this) == 1 **/
Vector2 Vector2::randOrthoB(const real len) const
{
    real s = RNG.sreal() * len;
    return Vector2(-YY, XX) * s;
}

//------------------------------------------------------------------------------
#pragma mark - 3D Vectors

Vector3 Vector3::randS()       { return Vector3(RNG.sreal(),     RNG.sreal(),   RNG.sreal()); }
Vector3 Vector3::randH()       { return Vector3(RNG.shalf(),     RNG.shalf(),   RNG.shalf()); }
Vector3 Vector3::randS(real n) { return Vector3(n*RNG.sreal(), n*RNG.sreal(), n*RNG.sreal()); }
Vector3 Vector3::randP()       { return Vector3(RNG.preal(),     RNG.preal(),   RNG.preal()); }
Vector3 Vector3::randP(real n) { return Vector3(n*RNG.preal(), n*RNG.preal(), n*RNG.preal()); }
Vector3 Vector3::randG(real n) { return Vector3(n*RNG.gauss(), n*RNG.gauss(), n*RNG.gauss()); }
void Vector3::addRand(real n) { XX += n*RNG.sreal(); YY += n*RNG.sreal(); ZZ += n*RNG.sreal(); }


#if ( 0 )

/// hypercube rejection method
Vector3 Vector3::randU()
{
    real x, y, z, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
        d = x*x + y*y + z*z;
    } while ( d > 1.0  ||  d < 0.01 );
    return Vector3(x, y, z) / std::sqrt(d);

}

/// hypercube rejection method
Vector3 Vector3::randU(real n)
{
    real x, y, z, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
        d = x*x + y*y + z*z;
    } while ( d > 1.0  ||  d < 0.01 );
    return Vector3(x, y, z) * (n/std::sqrt(d));
}

#elif ( 1 )

/**
 Derived from Marsaglia (1972)
 Allen & Tildesley "Computer Simulation of Liquids" Clarendon Pres, Oxford 1987
 http://mathworld.wolfram.com/SpherePointPicking.html
 This uses only 2 random-numbers!
*/
Vector3 Vector3::randU()
{
    real x, y, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = 1.0 - x*x - y*y;
    } while ( d <= 0 );
    real h = 2.0 * std::sqrt(d);
    return Vector3(x*h, y*h, 1.0-2.0*d);
}

Vector3 Vector3::randU(const real n)
{
    real x, y, d;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        d = 1.0 - x*x - y*y;
    } while ( d <= 0 );
    real h = ( n + n ) * std::sqrt(d);
    return Vector3(x*h, y*h, n*(1.0-2.0*d));
}

#endif


Vector3 Vector3::randB()
{
    real x, y, z;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
    } while ( x*x + y*y + z*z > 1.0 );
    return Vector3(x, y, z);
}


Vector3 Vector3::randB(const real n)
{
    real x, y, z;
    do {
        x = RNG.sreal();
        y = RNG.sreal();
        z = RNG.sreal();
    } while ( x*x + y*y + z*z > 1.0 );
    return Vector3(x, y, z) * n;
}


/** this assumes norm(*this) == 1 **/
Vector3 Vector3::randOrthoU(const real len) const
{
    real C, S;
    RNG.urand2(C, S, len);
    return orthogonalNCS(1.0, C, S);
}

/** this assumes norm(*this) == 1 **/
Vector3 Vector3::randOrthoB(const real len) const
{
    const Vector2 V = Vector2::randB();
    return orthogonalNCS(1.0, len * V.XX, len * V.YY);
}

//------------------------------------------------------------------------------
#pragma mark - 4-component Vectors treated as 3D

Vector4 Vector4::randS()       { return Vector4(RNG.sreal(),     RNG.sreal(),   RNG.sreal()); }
Vector4 Vector4::randH()       { return Vector4(RNG.shalf(),     RNG.shalf(),   RNG.shalf()); }
Vector4 Vector4::randS(real n) { return Vector4(n*RNG.sreal(), n*RNG.sreal(), n*RNG.sreal()); }
Vector4 Vector4::randP()       { return Vector4(RNG.preal(),     RNG.preal(),   RNG.preal()); }
Vector4 Vector4::randP(real n) { return Vector4(n*RNG.preal(), n*RNG.preal(), n*RNG.preal()); }
Vector4 Vector4::randG(real n) { return Vector4(n*RNG.gauss(), n*RNG.gauss(), n*RNG.gauss()); }
void  Vector4::addRand(real n) { XX += n*RNG.sreal(); YY += n*RNG.sreal(); ZZ += n*RNG.sreal(); }


//------------------------------------------------------------------------------
#pragma mark - Functions to distribute multiple points


/// Distribute points on the ball ( norm <= 1 ).
size_t tossPointsBall(std::vector<Vector3>& pts, real sep, size_t max_trials)
{
    const real ss = sep * sep;
    size_t ouf = 0;
    size_t n = 0;
    
    for ( Vector3& vec : pts )
    {
    toss:
        if ( ++ouf > max_trials )
            break;
        
        const Vector3 V = Vector3::randB();
        
        for ( size_t i = 0; i < n; ++i )
            if ( distanceSqr(V, pts[i]) < ss )
                goto toss;
        
        vec = V;
        ouf = 0;
        ++n;
    }
    return n;
}


/**
 Generate a random distribution of points on the unit disc,
 with the distance between two points never below `sep`.
 @return number of points stored in 'pts[]'
 */
size_t tossPointsDisc(std::vector<Vector2>& pts, real sep, size_t max_trials)
{
    const real ss = sep * sep;
    size_t ouf = 0;
    size_t n = 0;
    
    for ( Vector2& vec : pts )
    {
    toss:
        if ( ++ouf > max_trials )
            break;
        
        const Vector2 V = Vector2::randB();
        
        for ( size_t i = 0; i < n; ++i )
            if ( distanceSqr(V, pts[i]) < ss )
                goto toss;
        
        vec = V;
        ouf = 0;
        ++n;
    }
    return n;
}


size_t tossPointsCap(std::vector<Vector2>& pts, real cap, real sep, size_t max_trials)
{
    assert_true( cap < 2.0 );
    const real ss = sep * sep;
    size_t ouf = 0;
    size_t n = 0;
    real angle = std::acos(1.0-cap);
    for ( Vector2& vec : pts )
    {
    toss:
        if ( ++ouf > max_trials )
            break;
        
        real a = angle * RNG.sreal();

        for ( size_t i = 0; i < n; ++i )
            if ( abs_real(a-pts[i].XX) < ss )
                goto toss;
        
        vec.XX = a;
        ouf = 0;
        ++n;
    }
    for ( Vector2& vec : pts )
        vec.set(std::cos(vec.XX), std::sin(vec.XX));
    return n;
}

/**
 Generate a random distribution of points over a portion of the unit sphere,
 near (1,0,0) and defined by a thickness `cap`, that is symmetric around the X-axis.
 No two points should be closer than `sep`.
 @return number of points stored in `pts[]`
 */
size_t tossPointsCap(std::vector<Vector3>& pts, real cap, real sep, size_t max_trials)
{
    assert_true( cap < 2.0 );
    const real ss = sep * sep;
    size_t ouf = 0;
    size_t n = 0;
    real C, S;

    for ( Vector3& vec : pts )
    {
    toss:
        if ( ++ouf > max_trials )
            break;
        
        real u = max_real(1.0 - cap * RNG.preal(), -1.0);
        real v = std::sqrt(1.0 - u*u);
        RNG.urand2(C, S, v);
        Vector3 pos(u, C, S);
        
        for ( size_t i = 0; i < n; ++i )
            if ( distanceSqr(pos, pts[i]) < ss )
                goto toss;
        
        vec = pos;
        ouf = 0;
        ++n;
    }
    return n;
}


size_t distributePointsSphere(std::vector<Vector3>& pts, real sep, size_t max_trials)
{
    /*
     Estimate max separation by dividing area, given highest packing possible:
     Each hexagon covers 3 discs of diameter D has surface area 3/4*PI*D*D,
     Since sphere surface is 4 * PI, the surface of hexagon = 3*4*PI/nb_points
     Hence 1/4 * D*D = 4 / nb_points
     */
    real sup = std::sqrt(16.0/pts.size());

    real dis = std::min(sep, sup);
    size_t ouf = 0;
    size_t res = 0;
    while ( res < pts.size() )
    {
        res = tossPointsSphere(pts, dis, max_trials);
        if ( ++ouf > 64 )
        {
            ouf = 0;
            dis /= 1.0905044;
        }
    }
    if ( dis < sep )
        Cytosim::log("distributePointsSphere separation: ", sep, " --> ", dis, "\n");
    return res;
}
