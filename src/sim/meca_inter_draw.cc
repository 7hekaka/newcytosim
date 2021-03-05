// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University


#include "../disp/gle.h"
#include "../disp/gle_color_list.h"

#define DRAW_LINK(PTA, ...)\
{ if ( drawLinks ) drawLink(gle::bright_color(PTA.mecable()->signature()), __VA_ARGS__); }


/// Display link between 2 positions
void drawLink(gle_color const& col, Vector const& a, Vector const& b)
{
    col.load();
    glLineStipple(1, 0xFFFF);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(b);
    glEnd();
}

/// Display link between 2 positions, with resting length `len`
void drawLink(gle_color const& col, Vector const& a, Vector const& ab, real len)
{
    col.load();
    Vector b = a + ab;
    Vector dx = ab * (( 1 - len / ab.norm() ) / 2);
    glLineStipple(1, 0x3333);
    glBegin(GL_LINES);
    gle::gleVertex(a+dx);
    gle::gleVertex(b-dx);
    glEnd();
    glLineStipple(1, 0xFFFF);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(a+dx);
    gle::gleVertex(b-dx);
    gle::gleVertex(b);
    glEnd();
    glBegin(GL_POINTS);
    gle::gleVertex(a);
    gle::gleVertex(b);
    glEnd();
}

/// Display link between 3 positions
void drawLink(gle_color const& col, Vector const& a, Vector const& ab, Vector c)
{
    col.load();
    if ( modulo )
        modulo->fold(c, a);
    Vector b = a + ab;
    glLineStipple(1, 0x7310);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(b);
    glEnd();
    glLineStipple(1, 0x5555);
    glBegin(GL_LINES);
    gle::gleVertex(b);
    gle::gleVertex(c);
    glEnd();
    glBegin(GL_POINTS);
    gle::gleVertex(b);
    glEnd();
}

/// Display link between 4 positions
void drawLink(gle_color const& col, Vector const& a, Vector const& ab, Vector const& dc, Vector const& d)
{
    col.load();
    Vector b = a + ab;
    Vector c = d + dc;
    glLineStipple(1, 0x7171);
    glBegin(GL_LINES);
    gle::gleVertex(a);
    gle::gleVertex(b);
    gle::gleVertex(c);
    gle::gleVertex(d);
    glEnd();
    glLineStipple(1, 0xFFFF);
    glBegin(GL_LINES);
    gle::gleVertex(b);
    gle::gleVertex(c);
    glEnd();
    glBegin(GL_POINTS);
    gle::gleVertex(b);
    gle::gleVertex(c);
    glEnd();
}
