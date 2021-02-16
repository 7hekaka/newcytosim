// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Created 09/03/2015 by Francois Nedelec

#ifndef GRID_DISPLAY_H
#define GRID_DISPLAY_H

#include "vector1.h"
#include "vector2.h"
#include "vector3.h"
#include "map.h"
#include "grid.h"
#include "opengl.h"
#include "gle.h"
#include "flute.h"

/// display the edges of a 1D grid using OpenGL
void drawEdges(Map<1> const&);

/// display the edges of a 2D grid using OpenGL
void drawEdges(Map<2> const&);

/// display the edges of a 3D grid using OpenGL
void drawEdges(Map<3> const&);

//------------------------------------------------------------------------------
#pragma mark -


/// display the values stored in the cells of a 1D grid using OpenGL
/**
 OpenGL color is to be specified by the provided function:
 bool set_color(void*, CELL const&, Vector2 const&);
 Each particular cell is displayed only if `set_color' returns true.
 */
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 1> const& grid,
                gle_color color(TYPE, CELL const&, Vector1 const&),
                TYPE arg)
{
    float d = grid.cellWidth(0);
    float e = 2;
    
    // paint all cells one by one
    for ( size_t c = 0; c < grid.breadth(0); ++c )
    {
        float x = grid.position(0, c);
        gle_color col = color(arg, grid[c], Vector1(x));
        if ( col.visible() )
        {
            col.load();
            GLfloat pts[8] = {x, -e, x+d, -e, x, e, x+d, e};
            glVertexPointer(2, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }
    }
}


/// display the values stored in the cells of a 2D grid using OpenGL
/**
 OpenGL color is to be specified by the provided function:
 bool set_color(void*, CELL const&, Vector2 const&);
 Each particular cell is displayed only if `set_color' returns true.
*/
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 2> const& grid,
                gle_color color(TYPE, CELL const&, Vector2 const&),
                TYPE arg)
{
    float d = 0.5 * grid.cellWidth(0);
    float e = 0.5 * grid.cellWidth(1);

    // paint all cells one by one
    for ( size_t c = 0; c < grid.nbCells(); ++c )
    {
        Vector2 w;
        grid.setPositionFromIndex(w, c, 0.5);
        gle_color col = color(arg, grid[c], w);
        if ( col.visible() )
        {
            col.load();
            GLfloat X(w.XX), Y(w.YY);
            GLfloat pts[8] = {X-d, Y-e, X+d, Y-e, X-d, Y+e, X+d, Y+e};
            glVertexPointer(2, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }
    }
}


/// display the slice of a 3D grid in a plane parallel to XY at `Z = zzz`
/**
 a Cell color is specified by `bool set_color(TYPE, CELL const&, Vector2 const&)`
 Each particular cell is displayed only if `set_color' returns true.
 */
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 3> const& grid,
                gle_color color(TYPE, CELL const&, Vector3 const&),
                TYPE arg,
                real zzz = 0)
{
    assert_true(grid.hasCells());
    GLfloat d(0.5 * grid.cellWidth(0));
    GLfloat e(0.5 * grid.cellWidth(1));
    
    size_t z = grid.index(2, zzz);
    for ( size_t y = 0; y < grid.breadth(1); ++y )
    for ( size_t x = 0; x < grid.breadth(0); ++x )
    {
        Vector3 w(grid.position(0, x+0.5), grid.position(1, y+0.5), zzz);
        gle_color col = color(arg, grid.icell3D(x,y,z), w);
        if ( col.visible() )
        {
            col.load();
            GLfloat X(w.XX), Y(w.YY), Z(zzz);
            GLfloat pts[12] = {X-d, Y-e, Z, X+d, Y-e, Z, X-d, Y+e, Z, X+d, Y+e, Z};
            glVertexPointer(3, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }
    }
}


/// display the slice of a 3D grid in a plane parallel to Y: `Y=pos`
/**
 OpenGL color is to be specified by the provided function:
 bool set_color(void*, CELL const&, Vector2 const&);
 Each particular cell is displayed only if `set_color' returns true.
*/
template <typename CELL, typename TYPE>
void drawValuesXZ(Grid<CELL, 3> const& grid,
                  gle_color color(TYPE, CELL const&, Vector3 const&),
                  TYPE arg,
                  real yyy)
{
    assert_true(grid.hasCells());
    GLfloat d(0.5 * grid.cellWidth(0));
    GLfloat e(0.5 * grid.cellWidth(1));
    
    size_t y = grid.index(1, yyy);
    for ( size_t z = 0; z < grid.breadth(2); ++z )
    for ( size_t x = 0; x < grid.breadth(0); ++x )
    {
        Vector3 w(grid.position(0, x+0.5), yyy, grid.position(2, y+0.5));
        gle_color col = color(arg, grid.icell3D(x,y,z), w);
        if ( col.visible() )
        {
            col.load();
            GLfloat X(w.XX), Y(yyy), Z(w.ZZ);
            GLfloat pts[12] = {X-d, Y, Z-e, X+d, Y, Z-e, X-d, Y, Z+e, X+d, Y, Z+e};
            glVertexPointer(3, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }
    }
}


// display the slice of a 3D grid in a plane parallel to X: `X=pos`
/**
 OpenGL color is to be specified by the provided function:
 bool set_color(void*, CELL const&, Vector2 const&);
 Each particular cell is displayed only if `set_color' returns true.
 */
template <typename CELL, typename TYPE>
void drawValuesYZ(Grid<CELL, 3> const& grid,
                  gle_color color(TYPE, CELL const&, Vector3 const&),
                  TYPE arg,
                  real xxx)
{
    assert_true(grid.hasCells());
    GLfloat d(0.5 * grid.cellWidth(0));
    GLfloat e(0.5 * grid.cellWidth(2));
    
    size_t x = grid.index(0, xxx);
    for ( size_t z = 0; z < grid.breadth(2); ++z )
    for ( size_t y = 0; y < grid.breadth(1); ++y )
    {
        Vector3 w(xxx, grid.position(1, y+0.5), grid.position(2, z+0.5));
        gle_color col = color(arg, grid.icell3D(x,y,z), w);
        if ( col.visible() )
        {
            col.load();
            GLfloat X(xxx), Y(w.YY), Z(w.ZZ);
            GLfloat pts[12] = {X, Y-d, Z-e, X, Y+d, Z-e, X, Y-d, Z+e, X, Y+d, Z+e};
            glVertexPointer(3, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }
    }
}


/// display a slice of the field in a plane perpendicular to 'dir'
/**
 a Cell color is specified by `bool set_color(TYPE, CELL const&, Vector2 const&)`
 The return value of this function is ignored.
 */
template <typename CELL, typename TYPE>
void drawValues(Grid<CELL, 3> const& grid,
                gle_color color(TYPE, CELL const&, Vector3 const&),
                TYPE arg,
                Vector3 const& dir,
                real pos)
{
    assert_true(grid.hasCells());
    
    // this defines the finesse of the triangular mesh:
    real cel = 0.2 * grid.minimumWidth(1);
    int R = (int)ceil( grid.radius() / cel );

    Vector3 dx, dy;
    dir.orthonormal(dx, dy, cel);
    dy *= M_SQRT3_2;
    
    flute3 * pts = new flute3[4*R+2];
    flute4 * col = new flute4[4*R+2];
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_FLOAT, 0, col);
    glVertexPointer(3, GL_FLOAT, 0, pts);

    for ( int y = -R; y <= R; y+=2 )
    {
        Vector3 A = y * dy + pos * dir;
        Vector3 B = A + dy + 0.5 * dx;
        size_t i = 0;
        for ( int n = -R; n <= R; ++n )
        {
            Vector3 V = A + n * dx;
            Vector3 W = B + n * dx;
            col[i] = color(arg, grid.interpolate3D(V), V);
            pts[i++] = V;
            col[i] = color(arg, grid.interpolate3D(W), W);
            pts[i++] = W;
        }
        glDrawArrays(GL_TRIANGLE_STRIP, 0, i);
        i = 0;
        A += 2 * dy;
        for ( int n = -R; n <= R; ++n )
        {
            Vector3 V = A + n * dx;
            Vector3 W = B + n * dx;
            col[i] = color(arg, grid.interpolate3D(V), V);
            pts[i++] = V;
            col[i] = color(arg, grid.interpolate3D(W), W);
            pts[i++] = W;
        }
        glDrawArrays(GL_TRIANGLE_STRIP, 0, i);
    }
    glDisableClientState(GL_COLOR_ARRAY);
    delete[] pts;
    delete[] col;
}


#endif


