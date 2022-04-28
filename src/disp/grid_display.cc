// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "grid_display.h"
#include "gym_draw.h"

/**
 This uses the current OpenGL color and line width.
 */
void drawBoundaries(Map<1> const& map, float width)
{
    size_t sup = 1 + map.breadth(0);
    flute2 * flt = gym::mapBufferV2(2*sup);
    for ( size_t n = 0; n < sup; ++n )
    {
        float x = map.position(0, n);
        flt[2*n  ] = { x, -0.5f };
        flt[2*n+1] = { x,  0.5f };
    }
    gym::unmapBufferV2();
    gym::drawLines(width, 0, 2*sup);
}


/**
 This uses the current OpenGL color and line width.
 */
void drawBoundaries(Map<2> const& map, float width)
{
    size_t supX = 1 + map.breadth(0);
    size_t supY = 1 + map.breadth(1);
    flute2 * flt = gym::mapBufferV2(2*std::max(supX, supY));

    float i = map.inf(0);
    float s = map.sup(0);
    for ( size_t n = 0; n < supY; ++n )
    {
        float y = map.position(1, n);
        flt[2*n  ] = { i, y };
        flt[2*n+1] = { s, y };
    }
    gym::unmapBufferV2();
    gym::drawLines(width, 0, 2*supY);
    
    i = map.inf(1);
    s = map.sup(1);
    for ( size_t n = 0; n < supX; ++n )
    {
        float x = map.position(0, n);
        flt[2*n  ] = { x, i };
        flt[2*n+1] = { x, s };
    }
    gym::unmapBufferV2();
    gym::drawLines(width, 0, 2*supX);
}


/**
 This uses the current OpenGL color and line width.
 */
void drawBoundaries(Map<3> const& map, float width)
{
    size_t supX = 1 + map.breadth(0);
    size_t supY = 1 + map.breadth(1);
    size_t supZ = 1 + map.breadth(2);

    float i = map.inf(0);
    float s = map.sup(0);
    for ( size_t iz = 0; iz < supZ; ++iz )
    {
        float z = map.position(2, iz);
        flute3 * flt = gym::mapBufferV3(2*supY);
        for ( size_t n = 0; n < supY; ++n )
        {
            float y = map.position(1, n);
            flt[2*n  ] = { i, y, z };
            flt[2*n+1] = { s, y, z };
        }
        gym::unmapBufferV3();
        gym::drawLines(width, 0, 2*supY);
    }
    
    i = map.inf(1);
    s = map.sup(1);
    for ( size_t iz = 0; iz < supZ; ++iz )
    {
        float z = map.position(2, iz);
        flute3 * flt = gym::mapBufferV3(2*supX);
        for ( size_t n = 0; n < supX; ++n )
        {
            float x = map.position(0, n);
            flt[2*n  ] = { x, i, z };
            flt[2*n+1] = { x, s, z };
        }
        gym::unmapBufferV3();
        gym::drawLines(width, 0, 2*supX);
    }

    i = map.inf(2);
    s = map.sup(2);
    for ( size_t ix = 0; ix < supX; ++ix )
    {
        float x = map.position(0, ix);
        flute3 * flt = gym::mapBufferV3(2*supY);
        for ( size_t n = 0; n < supY; ++n )
        {
            float y = map.position(1, n);
            flt[2*n  ] = { x, y, i };
            flt[2*n+1] = { x, y, s };
        }
        gym::unmapBufferV3();
        gym::drawLines(width, 0, 2*supY);
    }
}

