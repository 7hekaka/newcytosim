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


/// Distribute points on the unit disc, at distance never below `sep`.
size_t tossPointsDisc(std::vector<Vector2>& pts, real sep, size_t limit_trials);


/// Distribute points on the unit disc, over a spherical cap of thickness `cap`.
size_t tossPointsCap(std::vector<Vector3>& pts, real cap, real sep, size_t limit_trials);


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
size_t distributePointsSphere(std::vector<VECTOR>& pts, real sep, size_t max_trials)
{
    // estimate separation by dividing area:
    if ( sep <= 0 )
        sep = std::sqrt(4*M_PI/pts.size());

    real dis = sep;
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
        std::clog << "distributePointsSphere " << sep << " --> " << dis << "\n";
    return res;
}

#endif

