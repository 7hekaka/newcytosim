// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University


#include "../gym/gym_color.h"
#include "../gym/gym_color_list.h"
#include "../gym/gym_flute.h"
#include "../gym/gym_flute_dim.h"
#include "../gym/gym_draw.h"
#include "../gym/gym_cap.h"

#define DRAW_LINK(PT, ...)\
{ if ( drawLinks ) drawLink(gym::bright_color(PT.mecable()->signature()), PT.pos(), __VA_ARGS__); }

constexpr float LINE_WIDTH = 2;
constexpr float POINT_SIZE = 5;


/// Display link between 2 positions
void drawLink(gym_color const& col, Vector const& a, Vector const& b)
{
    flute4D* flu = gym::mapBufferC4VD(2);
    flu[0] = { col, a };
    flu[1] = { col, b };
    gym::unmapBufferC4VD();
    gym::enableLineStipple(0xFFFF);
    gym::drawLines(LINE_WIDTH, 0, 2);
}

/// Display link between 2 positions, with resting length `len`
void drawLink(gym_color const& col, Vector const& a, Vector const& ab, real len)
{
    Vector b = a + ab;
    Vector dx = ab * (( 1 - len / ab.norm() ) / 2);
    flute4D* flu = gym::mapBufferC4VD(4);
    flu[0] = { col, a };
    flu[1] = { col, a+dx };
    flu[2] = { col, b-dx };
    flu[3] = { col, b };
    gym::unmapBufferC4VD();
    gym::enableLineStipple(0x3333);
    gym::drawLines(LINE_WIDTH, 1, 2);
    gym::enableLineStipple(0xFFFF);
    gym::drawLines(LINE_WIDTH, 0, 2);
    gym::drawLines(LINE_WIDTH, 2, 2);
    gym::drawPoints(POINT_SIZE, 0, 1);
    gym::drawPoints(POINT_SIZE, 3, 1);
}

/// Display link between 3 positions
void drawLink(gym_color const& col, Vector const& a, Vector const& ab, Vector c)
{
    if ( modulo )
        modulo->fold(c, a);
    Vector b = a + ab;
    flute4D* flu = gym::mapBufferC4VD(4);
    flu[0] = { col, a };
    flu[1] = { col, b };
    flu[2] = { col, c };
    gym::unmapBufferC4VD();
    gym::enableLineStipple(0x7310);
    gym::drawLines(LINE_WIDTH, 0, 2);
    gym::enableLineStipple(0x5555);
    gym::drawLines(LINE_WIDTH, 1, 2);
    gym::drawPoints(POINT_SIZE, 1, 1);
}

/// Display link between 4 positions
void drawLink(gym_color const& col, Vector const& a, Vector const& ab, Vector const& dc, Vector const& d)
{
    Vector b = a + ab;
    Vector c = d + dc;
    flute4D* flu = gym::mapBufferC4VD(4);
    flu[0] = { col, a };
    flu[1] = { col, b };
    flu[2] = { col, c };
    flu[3] = { col, d };
    gym::unmapBufferC4VD();
    gym::enableLineStipple(0x7171);
    gym::drawLines(LINE_WIDTH, 0, 4);
    gym::enableLineStipple(0xFFFF);
    gym::drawLines(LINE_WIDTH, 1, 2);
    gym::drawPoints(POINT_SIZE, 1, 2);
}
