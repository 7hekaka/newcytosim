// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef GYM_FLUTE_DIM_H
#define GYM_FLUTE_DIM_H

#include "dim.h"

/*
 Attention:
 Since the types depend on the dimensionality, this unit must be compiled
 separately for each different value of DIM.
 */

#if ( DIM >= 3 )
typedef flute3 fluteD;
typedef flute8 flute4D;
#elif ( DIM == 2 )
typedef flute2 fluteD;
typedef flute6 flute4D;
#else
typedef flute1 fluteD;
typedef flute6 flute4D;
#endif

namespace gym
{
    /// reinterpret buffer with different layout
    inline void rebindBufferVD(size_t stride) { rebind(); setBufferV((DIM>2?3:2), stride); }
    
    /// reinterpret buffer with different layout
    inline void rebindBufferC4VD(size_t stride) { rebind(); setBufferCV(4, (DIM>2?4:2), stride); }

    inline fluteD* mapBufferVD(size_t n)
    {
        float* res = mapFloatBuffer((DIM>2?3:2)*n);
        return (fluteD*)res;
    }

    inline void unmapBufferVD()
    {
        unmap();
        setBufferV((DIM>2?3:2));
    }

    inline flute4D* mapBufferC4VD(size_t n)
    {
        float* res = mapFloatBuffer((DIM>2?8:6)*n);
        return (flute4D*)res;
    }

    inline void unmapBufferC4VD()
    {
        unmap();
        setBufferCV(4, (DIM>2?4:2));
    }
    
    inline void loadPoints(size_t cnt, const float pts[])
    {
#if ( DIM == 1 )
        // OpenGL cannot handle only 1 coordinate per vertex
        flute2 * flt = mapBufferV2(2*cnt);
        for ( size_t i = 0; i < cnt; ++i )
            flt[i].set(pts[i], 0.f);
#else
        float * flt = mapFloatBuffer(DIM*cnt);
        for ( size_t i = 0; i < DIM*cnt; ++i )
            flt[i] = pts[i];
#endif
        unmap();
        setBufferV((DIM>2?3:2));
    }

    inline void loadPoints(size_t cnt, const double pts[])
    {
#if ( DIM == 1 )
        // OpenGL cannot handle only 1 coordinate per vertex
        flute2 * flt = mapBufferV2(2*cnt);
        for ( size_t i = 0; i < cnt; ++i )
            flt[i].set(pts[i], 0.f);
#else
        float * flt = mapFloatBuffer(DIM*cnt);
        for ( size_t i = 0; i < DIM*cnt; ++i )
            flt[i] = float(pts[i]);
#endif
        unmap();
        setBufferV((DIM>2?3:2));
    }

};

#endif
