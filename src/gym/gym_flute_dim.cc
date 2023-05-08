// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "gym_flute_dim.h"

float * gym::last_map = nullptr;
size_t gym::last_cnt = 0;

void gym::loadPoints(size_t cnt, const float pts[])
{
#if ( DIM == 1 )
    // OpenGL cannot have only 1 coordinate per vertex
    flute2 * flt = mapBufferV2(2*cnt);
    for ( size_t i = 0; i < cnt; ++i )
        flt[i].set(pts[i], 0.f);
#else
    float * flt = mapFloatBuffer(DIM*cnt);
    for ( size_t i = 0; i < DIM*cnt; ++i )
        flt[i] = pts[i];
#endif
    unmapBufferVD();
}

void gym::loadPoints(size_t cnt, const double pts[])
{
#if ( DIM == 1 )
    // OpenGL cannot have only 1 coordinate per vertex
    flute2 * flt = mapBufferV2(2*cnt);
    for ( size_t i = 0; i < cnt; ++i )
        flt[i].set(pts[i], 0.f);
#else
    float * flt = mapFloatBuffer(DIM*cnt);
    for ( size_t i = 0; i < DIM*cnt; ++i )
        flt[i] = float(pts[i]);
#endif
    unmapBufferVD();
}
