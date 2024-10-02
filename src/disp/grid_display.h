// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Created 09/03/2015 by Francois Nedelec

#ifndef GRID_DISPLAY_H
#define GRID_DISPLAY_H

#include "vector1.h"
#include "vector2.h"
#include "vector3.h"
#include "map.h"
#include "grid.h"
#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"

/// display the edges of a 1D grid using OpenGL
template<int ORD> void drawBoundaries(Map<ORD> const&, float) {};

void drawBoundaries(Map<1> const&, float);
void drawBoundaries(Map<2> const&, float);
void drawBoundaries(Map<3> const&, float);

//------------------------------------------------------------------------------
#pragma mark - paint cells, with a color-scale defined by an external function


/// display the values stored in the cells of a 1D grid using OpenGL
/**
 OpenGL color is to be specified by the argument function pointer:
 gym_color color(TYPE, CELL const&, Vector1 const&)
 */
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 1> const& grid,
                gym_color color(TYPE, CELL const&, Vector1 const&),
                TYPE arg)
{
    float dx = grid.cellWidth(0), cx = 0.5f * dx;
    float dy = 1;

    flute6 * flu = gym::mapBufferC4V2(4*grid.breadth(0)+4);
    flute6 * ptr = flu;
    for ( index_t ix = 0; ix < grid.breadth(0); ++ix )
    {
        float x = grid.position(0, ix);
        gym_color col = color(arg, grid.icell1D(ix), Vector1(x+cx));
        // use a full quad to achieve uniform color within
        ptr[0] = { col, x+dx, dy };
        ptr[1] = { col, x, 0 };
        ptr[2] = { col, x, dy };
        ptr[3] = { col, x+dx, 0 };
        ptr += 4;
    }
    gym::unmapBufferC4V2();
    gym::drawTriangleStrip(0, ptr-flu);
    gym::cleanupCV();
}


/// display the values stored in the cells of a 2D grid using OpenGL
/**
 OpenGL color is to be specified by the argument function pointer:
 gym_color color(TYPE, CELL const&, Vector1 const&)
*/
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 2> const& grid,
                gym_color color(TYPE, CELL const&, Vector2 const&),
                TYPE arg)
{
    float dx = grid.cellWidth(0), cx = 0.5f * dx;
    float dy = grid.cellWidth(1), cy = 0.5f * dy;
 
    for ( index_t iy = 0; iy < grid.breadth(1); ++iy )
    {
        flute6 * flu = gym::mapBufferC4V2(4*grid.breadth(0)+4);
        flute6 * ptr = flu;
        for ( index_t ix = 0; ix < grid.breadth(0); ++ix )
        {
            float x = grid.position(0, ix);
            float y = grid.position(1, iy);
            gym_color col = color(arg, grid.icell2D(ix, iy), Vector2(x+cx, y+cy));
            // use a full quad to achieve uniform color within
            ptr[0] = { col, x, y+dy };
            ptr[1] = { col, x, y };
            ptr[2] = { col, x+dx, y+dy };
            ptr[3] = { col, x+dx, y };
            ptr += 4;
        }
        gym::unmapBufferC4V2();
        gym::drawTriangleStrip(0, ptr-flu);
    }
    gym::cleanupCV();
}


/// display the slice of a 3D grid in a plane parallel to XY at `Z = zzz`
/**
 OpenGL color is to be specified by the argument function pointer:
 gym_color color(TYPE, CELL const&, Vector1 const&)
 */
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 3> const& grid,
                gym_color color(TYPE, CELL const&, Vector3 const&),
                TYPE arg,
                real zzz = 0)
{
    assert_true(grid.hasCells());
    float dx = grid.cellWidth(0), cx = 0.5f * dx;
    float dy = grid.cellWidth(1), cy = 0.5f * dy;
    float z = (float)zzz;
    
    index_t iz = grid.index(2, zzz);
    for ( index_t iy = 0; iy < grid.breadth(1); ++iy )
    {
        flute8 * flu = gym::mapBufferC4V4(4*grid.breadth(0)+4);
        flute8 * ptr = flu;
        for ( index_t ix = 0; ix < grid.breadth(0); ++ix )
        {
            float x = grid.position(0, ix);
            float y = grid.position(1, iy);
            gym_color col = color(arg, grid.icell3D(ix, iy, iz), Vector3(x+cx, y+cy, zzz));
            // use a full quad to achieve uniform color within
            ptr[0] = { col, x, y+dy, z };
            ptr[1] = { col, x, y, z };
            ptr[2] = { col, x+dx, y+dy, z };
            ptr[3] = { col, x+dx, y, z };
            ptr += 4;
        }
        gym::unmapBufferC4V4();
        gym::drawTriangleStrip(0, ptr-flu);
    }
    gym::cleanupCV();
}


/// display the slice of a 3D grid in a plane parallel to Y: `Y=pos`
/**
 OpenGL color is to be specified by the argument function pointer:
 gym_color color(TYPE, CELL const&, Vector1 const&)
*/
template <typename CELL, typename TYPE>
void drawValuesXZ(Grid<CELL, 3> const& grid,
                  gym_color color(TYPE, CELL const&, Vector3 const&),
                  TYPE arg,
                  real yyy)
{
    assert_true(grid.hasCells());
    float dx = grid.cellWidth(0), cx = 0.5f * dx;
    float dz = grid.cellWidth(1), cz = 0.5f * dz;
    float y = (float)yyy;

    index_t iy = grid.index(1, yyy);
    for ( index_t iz = 0; iz < grid.breadth(2); ++iz )
    {
        flute8 * flu = gym::mapBufferC4V4(4*grid.breadth(0)+4);
        flute8 * ptr = flu;
        for ( index_t ix = 0; ix < grid.breadth(0); ++ix )
        {
            float x = grid.position(0, ix);
            float z = grid.position(2, iz);
            gym_color col = color(arg, grid.icell3D(ix, iy, iz), Vector3(x+cx, yyy, z+cz));
            ptr[0] = { col, x, y, z+dz };
            ptr[1] = { col, x, y, z };
            ptr[2] = { col, x+dx, y, z+dz };
            ptr[3] = { col, x+dx, y, z };
            ptr += 4;
        }
        gym::unmapBufferC4V4();
        gym::drawTriangleStrip(0, ptr-flu);
    }
    gym::cleanupCV();
}


// display the slice of a 3D grid in a plane parallel to X: `X=pos`
/**
 OpenGL color is to be specified by the argument function pointer:
 gym_color color(TYPE, CELL const&, Vector1 const&)
 */
template <typename CELL, typename TYPE>
void drawValuesYZ(Grid<CELL, 3> const& grid,
                  gym_color color(TYPE, CELL const&, Vector3 const&),
                  TYPE arg,
                  real xxx)
{
    assert_true(grid.hasCells());
    float dy = grid.cellWidth(0), cy = 0.5f * dy;
    float dz = grid.cellWidth(1), cz = 0.5f * dz;
    float x = (float)xxx;

    index_t ix = grid.index(0, xxx);
    for ( index_t iz = 0; iz < grid.breadth(2); ++iz )
    {
        flute8 * flu = gym::mapBufferC4V4(4*grid.breadth(1)+4);
        flute8 * ptr = flu;
        for ( index_t iy = 0; iy < grid.breadth(1); ++iy )
        {
            float y = grid.position(1, iy);
            float z = grid.position(2, iz);
            gym_color col = color(arg, grid.icell3D(ix, iy, iz), Vector3(xxx, y+cy, z+cz));
            ptr[0] = { col, x, y, z+dz };
            ptr[1] = { col, x, y, z };
            ptr[2] = { col, x, y+dy, z+dz };
            ptr[3] = { col, x, y+dy, z };
            ptr += 4;
        }
        gym::unmapBufferC4V4();
        gym::drawTriangleStrip(0, ptr-flu);
    }
    gym::cleanupCV();
}


/// display a slice of the 3D field in a plane perpendicular to 'dir'
/**
 a Cell color is specified by `bool set_color(TYPE, CELL const&, Vector2 const&)`
 The return value of this function is ignored.
 */
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 3> const& grid,
                gym_color color(TYPE, CELL const&, Vector3 const&),
                TYPE arg,
                Vector3 const& dir,
                real alpha)
{
    assert_true(grid.hasCells());
    
    // this defines the finesse of the triangular mesh:
    real cel = 0.2 * grid.minimumWidth(1);
    int R = (int)ceil( grid.radius() / cel );

    Vector3 dx, dy;
    dir.orthonormal(dx, dy, cel);
    dy *= 0.5 * M_SQRT3;
    
    for ( int y = -R; y <= R; y+=2 )
    {
        flute8 * flu = gym::mapBufferC4V4(4*R+2);
        Vector3 A = y * dy + alpha * dir;
        Vector3 B = A + dy + 0.5 * dx;
        flute8 * ptr = flu;
        for ( int n = -R; n <= R; ++n )
        {
            Vector3 V = A + n * dx;
            Vector3 W = B + n * dx;
            ptr[0] = { color(arg, grid.interpolate3D(V), V), V };
            ptr[1] = { color(arg, grid.interpolate3D(W), W), W };
            ptr += 2;
        }
        gym::unmapBufferC4V4();
        gym::drawTriangleStrip(0, ptr-flu);
        flu = gym::mapBufferC4V4(4*R+2);
        ptr = flu;
        A += 2 * dy;
        for ( int n = -R; n <= R; ++n )
        {
            Vector3 V = A + n * dx;
            Vector3 W = B + n * dx;
            ptr[0] = { color(arg, grid.interpolate3D(V), V), V };
            ptr[1] = { color(arg, grid.interpolate3D(W), W), W };
            ptr += 2;
       }
        gym::unmapBufferC4V4();
        gym::drawTriangleStrip(0, ptr-flu);
    }
    gym::cleanupCV();
}


#endif

