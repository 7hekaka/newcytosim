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
typedef flute8 flute4D;
#else
typedef flute2 fluteD;
typedef flute6 flute4D;
#endif

namespace gym
{
    /// load points of dimension 'dim' into mapped buffer
    void loadPoints(size_t, const float pts[]);
    
    /// load points of dimension 'dim' into mapped buffer
    void loadPoints(size_t, const double pts[]);
    
    /// reinterpret buffer with different layout
    inline void rebindBufferVD(size_t gap) { rebind(); setBufferV((DIM>2?3:2), gap); }
    
    /// reinterpret buffer with different layout
    inline void rebindBufferC4VD(size_t gap) { rebind(); setBufferCV(4, (DIM>2?4:2), gap); }
    
    /// these variables are needed in 1D to reset the Y-coordinates
    extern float * last_map;
    /// these variables are needed in 1D to reset the Y-coordinates
    extern size_t last_cnt;

#if DIM > 1
    /// map / unmap GPU buffer for vertex data only
    inline fluteD* mapBufferVD(size_t n) { return (fluteD*)mapFloatBuffer((DIM>2?3:2)*n); }
    inline void  unmapBufferVD(bool = 1) { unmap(); setBufferV((DIM>2?3:2)); }
    
    /// map / unmap GPU buffer for color + vertex data
    inline flute4D* mapBufferC4VD(size_t n) { return (flute4D*)mapFloatBuffer((DIM>2?8:6)*n); }
    inline void   unmapBufferC4VD(bool = 1) { unmap(); setBufferCV(4, (DIM>2?4:2)); }
#else
    inline fluteD* mapBufferVD(size_t n)
    {
        last_cnt = n;
        last_map = mapFloatBuffer(2*n);
        return (fluteD*)last_map;
    }

    inline void unmapBufferVD(bool reset = true)
    {
        if ( reset )
        for ( size_t i = 0; i < last_cnt; ++i )
            last_map[2*i+1] = 0.f;
        unmap();
        setBufferV(2);
    }

    inline flute4D* mapBufferC4VD(size_t n)
    {
        last_cnt = n;
        last_map = mapFloatBuffer(6*n);
        return (flute4D*)last_map;
    }

    inline void unmapBufferC4VD(bool reset = true)
    {
        if ( reset )
        for ( size_t i = 0; i < last_cnt; ++i )
            last_map[6*i+5] = 0.f;
        unmap();
        setBufferCV(4, 2);
    }
#endif
};

#endif
