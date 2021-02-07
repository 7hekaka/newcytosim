// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
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
                bool set_color(TYPE, CELL const&, Vector1 const&),
                TYPE arg)
{
    float d = grid.cellWidth(0);
    float e = 2;
    
    // paint all cells one by one
    for ( size_t c = 0; c < grid.breadth(0); ++c )
    {
        float x = grid.position(0, c);
        if ( set_color(arg, grid[c], Vector1(x)) )
        {
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
                bool set_color(TYPE, CELL const&, Vector2 const&),
                TYPE arg)
{
    float d = 0.5 * grid.cellWidth(0);
    float e = 0.5 * grid.cellWidth(1);

    // paint all cells one by one
    for ( size_t c = 0; c < grid.nbCells(); ++c )
    {
        Vector2 w;
        grid.setPositionFromIndex(w, c, 0.5);
        if ( set_color(arg, grid[c], w) )
        {
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
                bool set_color(TYPE, CELL const&, Vector3 const&),
                TYPE arg,
                real zzz = 0)
{
    assert_true(grid.hasCells());
    
    real d = 0.5 * grid.cellWidth(0);
    real e = 0.5 * grid.cellWidth(1);
    
    size_t z = grid.index(2, zzz);
    
    for ( size_t y = 0; y < grid.breadth(1); ++y )
    for ( size_t x = 0; x < grid.breadth(0); ++x )
    {
        Vector3 w(grid.position(0, x+0.5), grid.position(1, y+0.5), zzz);
        if ( set_color(arg, grid.icell3D(x,y,z), w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex3f(w.XX-d, w.YY-e, zzz);
            glVertex3f(w.XX+d, w.YY-e, zzz);
            glVertex3f(w.XX-d, w.YY+e, zzz);
            glVertex3f(w.XX+d, w.YY+e, zzz);
            glEnd();
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
                  bool set_color(TYPE, CELL const&, Vector3 const&),
                  TYPE arg,
                  real yyy)
{
    assert_true(grid.hasCells());
    
    real d = 0.5 * grid.cellWidth(0);
    real e = 0.5 * grid.cellWidth(2);
    
    size_t y = grid.index(1, yyy);
    
    for ( size_t z = 0; z < grid.breadth(2); ++z )
    for ( size_t x = 0; x < grid.breadth(0); ++x )
    {
        Vector3 w(grid.position(0, x+0.5), yyy, grid.position(2, y+0.5));
        if ( set_color(arg, grid.icell3D(x,y,z), w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex3f(w.XX-d, yyy, w.ZZ-e);
            glVertex3f(w.XX+d, yyy, w.ZZ-e);
            glVertex3f(w.XX-d, yyy, w.ZZ+e);
            glVertex3f(w.XX+d, yyy, w.ZZ+e);
            glEnd();
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
                  bool set_color(TYPE, CELL const&, Vector3 const&),
                  TYPE arg,
                  real xxx)
{
    assert_true(grid.hasCells());
    
    real d = 0.5 * grid.cellWidth(0);
    real e = 0.5 * grid.cellWidth(2);
    
    size_t x = grid.index(0, xxx);
    
    for ( size_t z = 0; z < grid.breadth(2); ++z )
    for ( size_t y = 0; y < grid.breadth(1); ++y )
    {
        Vector3 w(xxx, grid.position(1, y+0.5), grid.position(2, z+0.5));
        if ( set_color(arg, grid.icell3D(x,y,z), w) )
        {
            glBegin(GL_TRIANGLE_STRIP);
            glVertex3f(xxx, w.YY-d, w.ZZ-e);
            glVertex3f(xxx, w.YY+d, w.ZZ-e);
            glVertex3f(xxx, w.YY-d, w.ZZ+e);
            glVertex3f(xxx, w.YY+d, w.ZZ+e);
            glEnd();
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
                bool set_color(TYPE, CELL const&, Vector3 const&),
                TYPE arg,
                Vector3 const& dir,
                real pos)
{
    assert_true(grid.hasCells());
    
    // this defines the finesse of the triangular mesh:
    real n = 0.2 * grid.minimumWidth(1);
    int m = (int)( grid.radius() / n );
    
    Vector3 dx, dy;
    dir.orthonormal(dx, dy, n);
    Vector3 dh = dy * std::cos(M_PI/6);
    
    for ( int y = -m; y <= m; y+=2 )
    {
        const Vector3 a = y * dh + pos * dir;
        const Vector3 b = a + dh;
        glBegin(GL_TRIANGLE_STRIP);
        for ( int x = -m; x <= m; ++x )
        {
            const Vector3 v = a + x * dx;
            set_color(arg, grid.interpolate3D(v), v);
            gle::gleVertex(v);
            
            const Vector3 w = b + ( x + 0.5 ) * dx;
            set_color(arg, grid.interpolate3D(w), w);
            gle::gleVertex(w);
        }
        glEnd();
        glBegin(GL_TRIANGLE_STRIP);
        for ( int x = -m; x <= m; ++x )
        {
            const Vector3 v = b + x * dx + dh;
            set_color(arg, grid.interpolate3D(v), v);
            gle::gleVertex(v);
            
            const Vector3 w = b + ( x + 0.5 ) * dx;
            set_color(arg, grid.interpolate3D(w), w);
            gle::gleVertex(w);
        }
        glEnd();
   }
}


#endif


