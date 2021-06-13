// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "grid_display.h"

/**
 This uses the current OpenGL color and line width.
 */
void drawBoundaries(Map<1> const& map)
{
    size_t sup = 1 + map.breadth(0);
    flute4 * pts = new flute4[sup];
    
    for ( size_t n = 0; n < sup; ++n )
    {
        float x = map.position(0, n);
        pts[n] = {x, -0.5f, x, 0.5f};
    }
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_LINES, 0, 2*sup);
    delete[] pts;
}


/**
 This uses the current OpenGL color and line width.
 */
void drawBoundaries(Map<2> const& map)
{
    size_t sup0 = 1 + map.breadth(0);
    size_t sup1 = 1 + map.breadth(1);
    flute4 * pts = new flute4[std::max(sup0, sup1)];

    float i = map.inf(0);
    float s = map.sup(0);
    for ( size_t n = 0; n < sup1; ++n )
    {
        float y = map.position(1, n);
        pts[n] = {i, y, s, y};
    }
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_LINES, 0, 2*sup1);
    
    i = map.inf(1);
    s = map.sup(1);
    for ( size_t n = 0; n < sup0; ++n )
    {
        float x = map.position(0, n);
        pts[n] = {x, i, x, s};
    }
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_LINES, 0, 2*sup0);
    delete[] pts;
}


/**
 This uses the current OpenGL color and line width.
 */
void drawBoundaries(Map<3> const& map)
{
    size_t sup0 = 1 + map.breadth(0);
    size_t sup1 = 1 + map.breadth(1);
    size_t sup2 = 1 + map.breadth(2);
    flute6 * pts = new flute6[std::max(sup1, sup2)];
    glVertexPointer(3, GL_FLOAT, 0, pts);

    GLfloat i = map.inf(0);
    GLfloat s = map.sup(0);

    for ( size_t iy = 0; iy < sup1; ++iy )
    {
        GLfloat y = map.position(1, iy);
        for ( size_t n = 0; n < sup2; ++n )
        {
            GLfloat z = map.position(2, n);
            pts[n] = {i, y, z, s, y, z};
        }
        glDrawArrays(GL_LINES, 0, 2*sup2);
    }
    
    i = map.inf(1);
    s = map.sup(1);
    GLfloat b = map.inf(2);
    GLfloat t = map.sup(2);
    for ( size_t ix = 0; ix < sup0; ++ix )
    {
        GLfloat x = map.position(0, ix);
        for ( size_t n = 0; n < sup2; ++n )
        {
            GLfloat z = map.position(2, n);
            pts[n] = {x, i, z, x, s, z};
        }
        glDrawArrays(GL_LINES, 0, 2*sup2);
        for ( size_t n = 0; n < sup1; ++n )
        {
            GLfloat y = map.position(1, n);
            pts[n] = {x, y, b, x, y, t};
        }
        glDrawArrays(GL_LINES, 0, 2*sup1);
    }
    delete[] pts;
}

