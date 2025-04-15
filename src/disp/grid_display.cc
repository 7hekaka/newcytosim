// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "grid_display.h"
#include "gym_draw.h"

/**
 This uses the current OpenGL color.
 */
void drawBoundaries(Map<1> const& map, float width)
{
    index_t sup = 1 + map.breadth(0);
    flute2 * flt = gym::mapBufferV2(2*sup);
    for ( index_t n = 0; n < sup; ++n )
    {
        float x = map.position(0, n);
        flt[0] = { x, -0.5f };
        flt[1] = { x,  0.5f };
        flt += 2;
    }
    gym::unmapBufferV2();
    gym::drawLines(width, 0, 2*sup);
}


/**
 This uses the current OpenGL color.
 */
void drawBoundaries(Map<2> const& map, float width)
{
    const index_t supX = 1 + map.breadth(0);
    const index_t supY = 1 + map.breadth(1);
    const index_t sup = 2 * ( supX + supY );
    flute2 * flt = gym::mapBufferV2(sup);

    float iX = map.inf(0);
    float sX = map.sup(0);
    for ( index_t n = 0; n < supY; ++n )
    {
        float y = map.position(1, n);
        flt[0] = { iX, y };
        flt[1] = { sX, y };
        flt += 2;
    }

    float iY = map.inf(1);
    float sY = map.sup(1);
    for ( index_t n = 0; n < supX; ++n )
    {
        float x = map.position(0, n);
        flt[0] = { x, iY };
        flt[1] = { x, sY };
        flt += 2;
    }
    gym::unmapBufferV2();
    gym::drawLines(width, 0, sup);
}


/**
 This uses the current OpenGL color.
 */
void drawBoundaries(Map<3> const& map, float width)
{
    const index_t supX = 1 + map.breadth(0);
    const index_t supY = 1 + map.breadth(1);
    const index_t supZ = 1 + map.breadth(2);

    float iX = map.inf(0);
    float sX = map.sup(0);
    flute3 * flt = gym::mapBufferV3(2*supY*supZ);
    for ( index_t iz = 0; iz < supZ; ++iz )
    {
        float z = map.position(2, iz);
        for ( index_t n = 0; n < supY; ++n )
        {
            float y = map.position(1, n);
            flt[0] = { iX, y, z };
            flt[1] = { sX, y, z };
            flt += 2;
        }
    }
    gym::unmapBufferV3();
    gym::drawLines(width, 0, 2*supY*supZ);

    float iY = map.inf(1);
    float sY = map.sup(1);
    flt = gym::mapBufferV3(2*supX*supZ);
    for ( index_t iz = 0; iz < supZ; ++iz )
    {
        float z = map.position(2, iz);
        for ( index_t n = 0; n < supX; ++n )
        {
            float x = map.position(0, n);
            flt[0] = { x, iY, z };
            flt[1] = { x, sY, z };
            flt += 2;
        }
    }
    gym::unmapBufferV3();
    gym::drawLines(width, 0, 2*supX*supZ);

    float iZ = map.inf(2);
    float sZ = map.sup(2);
    flt = gym::mapBufferV3(2*supX*supY);
    for ( index_t ix = 0; ix < supX; ++ix )
    {
        float x = map.position(0, ix);
        for ( index_t n = 0; n < supY; ++n )
        {
            float y = map.position(1, n);
            flt[0] = { x, y, iZ };
            flt[1] = { x, y, sZ };
            flt += 2;
        }
    }
    gym::unmapBufferV3();
    gym::drawLines(width, 0, 2*supX*supY);
}

