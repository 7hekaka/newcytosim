// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "gym_flute_dim.h"

#if DIM == 1

static flute4D * last_map = nullptr;
static size_t last_cnt = 0;


/// in 1D we need to reset the Y-coordinate
flute4D* gym::mapBufferC4VD(size_t n)
{
    last_cnt = n;
    last_map = (flute4D*)mapFloatBuffer(6*n);
    return last_map;
}

/// in 1D we need to reset the Y-coordinate
void gym::unmapBufferC4VD(bool reset)
{
    if ( reset )
    for ( size_t i = 0; i < last_cnt; ++i )
        last_map[i].setY(0.f);
    unmap();
    setBufferCV(4, (DIM>2?4:2));
}

#endif


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
