// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University


#include "../disp/gle.h"
#include "../disp/gle_color_list.h"
#include "../disp/gym_flute.h"

#define DRAW_LINK(PT, ...)\
{ if ( drawLinks ) drawLink(gle::bright_color(PT.mecable()->signature()), PT.pos(), __VA_ARGS__); }


/// Display link between 2 positions
void drawLink(gle_color const& col, Vector const& a, Vector const& b)
{
    fluteD4* flu = gym::mapBufferC4VD(2);
    flu[0] = { col, a };
    flu[1] = { col, b };
    gym::unmapBufferC4VD();
    glLineStipple(1, 0xFFFF);
    glDrawArrays(GL_LINES, 0, 2);
}

/// Display link between 2 positions, with resting length `len`
void drawLink(gle_color const& col, Vector const& a, Vector const& ab, real len)
{
    Vector b = a + ab;
    Vector dx = ab * (( 1 - len / ab.norm() ) / 2);
    fluteD4* flu = gym::mapBufferC4VD(4);
    flu[0] = { col, a };
    flu[1] = { col, a+dx };
    flu[2] = { col, b-dx };
    flu[3] = { col, b };
    gym::unmapBufferC4VD();
    glLineStipple(1, 0x3333);
    glDrawArrays(GL_LINES, 1, 2);
    glLineStipple(1, 0xFFFF);
    glDrawArrays(GL_LINES, 0, 2);
    glDrawArrays(GL_LINES, 2, 2);
    glDrawArrays(GL_POINTS, 0, 1);
    glDrawArrays(GL_POINTS, 3, 1);
}

/// Display link between 3 positions
void drawLink(gle_color const& col, Vector const& a, Vector const& ab, Vector c)
{
    if ( modulo )
        modulo->fold(c, a);
    Vector b = a + ab;
    fluteD4* flu = gym::mapBufferC4VD(4);
    flu[0] = { col, a };
    flu[1] = { col, b };
    flu[2] = { col, c };
    gym::unmapBufferC4VD();
    glLineStipple(1, 0x7310);
    glDrawArrays(GL_LINES, 0, 2);
    glLineStipple(1, 0x5555);
    glDrawArrays(GL_LINES, 1, 2);
    glDrawArrays(GL_POINTS, 1, 1);
}

/// Display link between 4 positions
void drawLink(gle_color const& col, Vector const& a, Vector const& ab, Vector const& dc, Vector const& d)
{
    Vector b = a + ab;
    Vector c = d + dc;
    fluteD4* flu = gym::mapBufferC4VD(4);
    flu[0] = { col, a };
    flu[1] = { col, b };
    flu[2] = { col, c };
    flu[3] = { col, d };
    gym::unmapBufferC4VD();
    glLineStipple(1, 0x7171);
    glDrawArrays(GL_LINES, 0, 4);
    glLineStipple(1, 0xFFFF);
    glDrawArrays(GL_LINES, 1, 2);
    glDrawArrays(GL_POINTS, 1, 2);
}
