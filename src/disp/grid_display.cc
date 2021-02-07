// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "grid_display.h"
using namespace gle;

/**
 This uses the current OpenGL color and line width.
 */
void drawEdges(Map<1> const& map)
{
    const float u =  0.5f;
    const float d = -0.5f;
    glBegin(GL_LINES);
    for ( real ix = 0; ix <= map.breadth(0); ++ix )
    {
        float x = map.position(0, ix);
        glVertex2f(x, d);
        glVertex2f(x, u);
    }
    glEnd();
}


/**
 This uses the current OpenGL color and line width.
 */
void drawEdges(Map<2> const& map)
{
    float i = map.inf(0);
    float s = map.sup(0);
    glBegin(GL_LINES);
    for ( float iy = 0; iy <= map.breadth(1); ++iy )
    {
        float y = map.position(1, iy);
        glVertex2f(i, y);
        glVertex2f(s, y);
    }
    glEnd();
    
    i = map.inf(1);
    s = map.sup(1);
    glBegin(GL_LINES);
    for ( float ix = 0; ix <= map.breadth(0); ++ix )
    {
        float x = map.position(0, ix);
        glVertex2f(x, i);
        glVertex2f(x, s);
    }
    glEnd();
}


/**
 This uses the current OpenGL color and line width.
 */
void drawEdges(Map<3> const& map)
{
    GLfloat i = map.inf(0);
    GLfloat s = map.sup(0);
    glBegin(GL_LINES);
    for ( real iy = 0; iy <= map.breadth(1); ++iy )
    for ( real iz = 0; iz <= map.breadth(2); ++iz )
    {
        GLfloat y = map.position(1, iy);
        GLfloat z = map.position(2, iz);
        glVertex3f(i, y, z);
        glVertex3f(s, y, z);
    }
    glEnd();
    
    i = map.inf(1);
    s = map.sup(1);
    glBegin(GL_LINES);
    for ( real ix = 0; ix <= map.breadth(0); ++ix )
    for ( real iz = 0; iz <= map.breadth(2); ++iz )
    {
        GLfloat x = map.position(0, ix);
        GLfloat z = map.position(2, iz);
        glVertex3f(x, i, z);
        glVertex3f(x, s, z);
    }
    glEnd();
    
    i = map.inf(2);
    s = map.sup(2);
    glBegin(GL_LINES);
    for ( real ix = 0; ix <= map.breadth(0); ++ix )
    for ( real iy = 0; iy <= map.breadth(1); ++iy )
    {
        GLfloat x = map.position(0, ix);
        GLfloat y = map.position(1, iy);
        glVertex3f(x, y, i);
        glVertex3f(x, y, s);
    }
    glEnd();
}

