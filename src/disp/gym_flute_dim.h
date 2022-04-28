// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef GYM_FLUTE_DIM_H
#define GYM_FLUTE_DIM_H

#include "gym_flute.h"
#include "dim.h"

/*
 Attention:
 Since the types depend on the dimensionality, this unit must be compiled
 separately for each different value of DIM.
 */

#if ( DIM >= 3 )
typedef flute3 fluteD;
typedef flute8 fluteD4;
#else
typedef flute2 fluteD;
typedef flute6 fluteD4;
#endif

namespace gym
{
    
    /// load points of dimension 'dim' into mapped buffer
    void loadPoints(size_t, const float pts[]);
    
    /// load points of dimension 'dim' into mapped buffer
    void loadPoints(size_t, const double pts[]);

    /// map / unmap GPU buffer for vertex data only
    inline fluteD* mapBufferVD(size_t n) { return (fluteD*)mapFloatBuffer((DIM>2?3:2)*n); }
    inline void  unmapBufferVD() { unmap(); setBufferV((DIM>2?3:2)); }
    inline void rebindBufferVD(size_t gap) { rebind(); setBufferV((DIM>2?3:2), gap); }

    /// map / unmap GPU buffer for color data + vertex
    inline fluteD4* mapBufferC4VD(size_t n) { return (fluteD4*)mapFloatBuffer((DIM>2?8:6)*n); }
    inline void   unmapBufferC4VD() { unmap(); setBufferCV(4, (DIM>2?4:2)); }
    inline void  rebindBufferC4VD(size_t gap) { rebind(); setBufferCV(4, (DIM>2?4:2), gap); }
};

#endif
