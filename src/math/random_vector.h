// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef RANDOM_VECTOR_H
#define RANDOM_VECTOR_H

#include <vector>
#include "vector2.h"
#include "vector3.h"

/**
 Functions to generate single random vectors are members of Vector1, Vector2, etc.
 and defined in random_vector.cc
 */


/// Distribute points on the ball ( norm <= 1 ).
size_t tossPointsBall(std::vector<Vector3>& pts, real sep, size_t max_trials);


/// Distribute points on the unit disc, at distance never below `sep`.
size_t tossPointsDisc(std::vector<Vector2>& pts, real sep, size_t limit_trials);

/// Distribute points on the unit circle, over an arc of thickness `cap`.
size_t tossPointsCap(std::vector<Vector2>& pts, real cap, real sep, size_t limit_trials);

/// Distribute points on the unit disc, over a spherical cap of thickness `cap`.
size_t tossPointsCap(std::vector<Vector3>& pts, real cap, real sep, size_t limit_trials);


/// Distribute points on the sphere (3D).
size_t distributePointsSphere(std::vector<Vector3>& pts, real sep, size_t max_trials);

/// Distribute points on the unit sphere, at distance never below `sep`.
/**
 Generate a random distribution of points on the unit circle,
 with the distance between two points never below `sep`.
 @return number of points stored in 'pts[]'
 */
template <typename VECTOR>
size_t tossPointsSphere(std::vector<VECTOR>& pts, real sep, size_t max_trials)
{
    const real ss = sep * sep;
    size_t ouf = 0;
    size_t n = 0;
    
    VECTOR pos;
    for ( VECTOR& vec : pts )
    {
    toss:
        if ( ++ouf > max_trials )
            break;
        
        pos = VECTOR::randU();
        
        // check distance will all the other points:
        for ( size_t i = 0; i < n; ++i )
            if ( distanceSqr(pos, pts[i]) < ss )
                goto toss;
        
        vec = pos;
        ouf = 0;
        ++n;
    }
    return n;
}

template <typename VECTOR>
size_t tossPointsCap(std::vector<VECTOR>& pts, real cap, size_t max_trials)
{
    size_t num = pts.size();
    size_t cnt = 0, ouf = 0;
    real sep, sup = std::sqrt( 2 * M_PI * cap / num );
    if ( VECTOR::dimensionality() == 2 )
        sup = M_SQRT2 * std::acos(1-cap) / num;
    do {
        // we decrease gradually the separation, to reach a good solution...
        sep = 512 * sup / real(ouf+512);
        cnt = tossPointsCap(pts, cap, sep, 1024);
        //std::clog << "tossCap(" << num << ") placed " << cnt << " with sep = " << sep << "\n";
        if ( ++ouf >= max_trials )
            break;
    } while ( cnt < num );
    return cnt;
}

#endif

