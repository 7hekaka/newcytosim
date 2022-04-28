// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef GYM_TEXT_H
#define GYM_TEXT_H

class gle_color;

/// Fonts inherited from GLUT
enum FontType
{
    BITMAP_9_BY_15 = 2,
    BITMAP_8_BY_13 = 3,
    BITMAP_TIMES_ROMAN_10 = 4,
    BITMAP_TIMES_ROMAN_24 = 5,
    BITMAP_HELVETICA_10 = 6,
    BITMAP_HELVETICA_12 = 7,
    BITMAP_HELVETICA_18 = 8
};

namespace gym
{
    
    /// set position of text in window
    float textPosition(float&, float&, int width, int height, int line, int W, int H, const int position);

    /// display text on a rectangle of color `bcol`, in a corner of the center of the display window
    void placeText(int position, FontType, const float color[4], const char text[], const float back[4], int width, int height);
    
    /// stroke character
    void strokeCharacter(float X, float Y, float S, int mono, const float color[4], unsigned char character, float width);

#ifdef VECTOR3_H
    /// draw `text` at position `pos`
    void drawText(Vector3 const& pos, FontType font, const float color[4], const char text[], float offset=0);
#endif
#ifdef VECTOR2_H
    /// draw `text` at position `pos`
    inline void drawText(Vector2 const& pos, FontType font, const float color[4], const char text[], float offset=0)
    {
        drawText(Vector3(pos.XX, pos.YY, 0), font, color, text, offset);
    }
#endif
#ifdef VECTOR1_H
    /// draw `text` at position `pos`
    inline void drawText(Vector1 const& pos, FontType font, const float color[4], const char text[], float offset=0)
    {
        drawText(Vector3(pos.XX, 0, 0), font, color, text, offset);
    }
#endif

}

#endif
