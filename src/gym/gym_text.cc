// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "gym_view.h"
#include "gym_draw.h"
#include "gle_color.h"
#include "gym_flute.h"
#include "gym_flat.h"
#include "vector1.h"
#include "vector2.h"
#include "vector3.h"
#include "fg_font.h"
#include "fg_stroke.h"
#include "gym_text.h"

/**
 draw text at position `vec`
 */
void gym::drawText(Vector3 const& vec, FontType font, const float color[4], const char text[], const float offset)
{
    gym::disableDepthTest();
    gym::disableAlphaTest();
    gym::disableLighting();
#if 1
    int H = fgFontHeight(font);
    int L, W = fgTextWidth(font, text, L);
    gym::translate(vec.XX, vec.YY, vec.ZZ);
    fgBitmapString(-W*offset, -H/3, 0.0625, font, color, text, H);
    gym::translate(-vec.XX, -vec.YY, -vec.ZZ);
#else
    gym::color(color);
    fgStrokeString(vec.XX, vec.YY, 0.001, 1, text, 1, 1, 0);
#endif
    gym::restoreDepthTest();
    gym::restoreLighting();
    gym::restoreAlphaTest();
}

/// width = text_width; height = text_heigth, (W, H) = window_size
float gym::textPosition(float& px, float& py, int width, int height, int lines,
                        int W, int H, const int position)
{
    assert_true( W > 0 );
    assert_true( H > 0 );
    
    switch( position )
    {
        case 0:
            //bottom-left, text going up
            px = height/2;
            py = height/2;
            return height;
        case 1:
            //bottom-right, text going up
            px = W - width - height/2;
            if ( px < 0 ) px = 0;
            py = height/2;
            return height;
        case 2:
            //top-right, text going down
            px = W - width - height/2;
            if ( px < 0 ) px = 0;
            py = H - height;
            return -height;
        default:
        case 3:
            //top-left, text going down
            px = height/2;
            py = H - height;
            return -height;
        case 4:
            //center, text going down
            px = ( W - width ) / 2;
            if ( px < 0 ) px = 0;
            py = ( H + lines*height ) / 2;
            return -height;
    }
    return height;
}

    
/**
 The text is displayed in the current color.
 A background rectangle is displayed only if `bcol` is visible.
 
 Possible values for `position`:
 - 0: bottom-left, text going up
 - 1: bottom-right, text going up
 - 2: top-right, text going down
 - 3: top-left, text going down
 - 4: center, text going down
 .
 
 Note: width and height are the current size of the viewport (window)
 */
void gym::placeText(int position, FontType font, const float color[4],
                    const char text[], const float back[4], int W, int H)
{
    int lines = 1;
    int height = fgFontHeight(font);
    int width = fgTextWidth(font, text, lines);
    
    float X, Y;
    float vshift = gym::textPosition(X, Y, width, height, lines, W, H, position);
    
    if ( back && back[3] > 0 )
    {
        float E = height;
        float T = Y + lines * vshift;
        float B = std::min(Y, T) - E/4;
        T = std::max(Y, T) + E + E/4;
        float R = X + width + E;
        gym::paintOctagon(X-E, B, R, T, back, 5);
        if ( position == 4 )
            gym::drawOctagon(X-E, B, R, T, color, 5, 1);
    }
    
    fgBitmapString(X, Y, 1.f, font, color, text, vshift);
}


/// stroke character
void gym::strokeCharacter(float X, float Y, float S, int mono, const float color[4], unsigned char character, float width)
{
    gym::color(color);
    fgStrokeCharacter(X, Y, S, mono, character, width, width);
}

