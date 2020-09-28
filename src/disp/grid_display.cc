// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "grid_display.h"


/**
 This uses the current OpenGL color and line width.
 */
void drawEdges(Map<1> const& map)
{
    const real u =  0.5;
    const real d = -0.5;
    glBegin(GL_LINES);
    for ( real ix = 0; ix <= map.breadth(0); ++ix )
    {
        real x = map.position(0, ix);
        gle::gleVertex(x, d);
        gle::gleVertex(x, u);
    }
    glEnd();
}


/**
 This uses the current OpenGL color and line width.
 */
void drawEdges(Map<2> const& map)
{
    real i = map.inf(0);
    real s = map.sup(0);
    glBegin(GL_LINES);
    for ( float iy = 0; iy <= map.breadth(1); ++iy )
    {
        real y = map.position(1, iy);
        gle::gleVertex(i, y);
        gle::gleVertex(s, y);
    }
    glEnd();
    
    i = map.inf(1);
    s = map.sup(1);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= map.breadth(0); ++ix )
    {
        real x = map.position(0, ix);
        gle::gleVertex(x, i);
        gle::gleVertex(x, s);
    }
    glEnd();
}


/**
 This uses the current OpenGL color and line width.
 */
void drawEdges(Map<3> const& map)
{
    real i = map.inf(0);
    real s = map.sup(0);
    glBegin(GL_LINES);
    for ( float iy = 0; iy <= map.breadth(1); ++iy )
    for ( float iz = 0; iz <= map.breadth(2); ++iz )
    {
        real y = map.position(1, iy);
        real z = map.position(2, iz);
        gle::gleVertex(i, y, z);
        gle::gleVertex(s, y, z);
    }
    glEnd();
    
    i = map.inf(1);
    s = map.sup(1);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= map.breadth(0); ++ix )
    for ( float iz = 0; iz <= map.breadth(2); ++iz )
    {
        real x = map.position(0, ix);
        real z = map.position(2, iz);
        gle::gleVertex(x, i, z);
        gle::gleVertex(x, s, z);
    }
    glEnd();
    
    i = map.inf(2);
    s = map.sup(2);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= map.breadth(0); ++ix )
    for ( float iy = 0; iy <= map.breadth(1); ++iy )
    {
        real x = map.position(0, ix);
        real y = map.position(1, iy);
        gle::gleVertex(x, y, i);
        gle::gleVertex(x, y, s);
    }
    glEnd();
}

