// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

class gle_color;
class Vector1;
class Vector2;
class Vector3;

/// Fonts inherited from GLUT
enum FontType
{
    STROKE_ROMAN = 0,
    STROKE_MONO_ROMAN = 1,
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

    /// return height in pixel of font
    int fontHeight(FontType font);
    
    /// compute horizontal size of text
    int maxTextWidth(FontType, const char text[], int& lines);
    
    /// draw text at the current OpenGL raster position and raster color
    void bitmapString(FontType, const char text[], float vshift);

    /// display text on a rectangle of color `bcol`, in a corner of the center of the display window
    void drawText(FontType, const char text[], gle_color bcol, int position, int width, int height);
    
    /// draw `text` at position `pos`
    void drawText(Vector1 const& pos, FontType, const char text[], float dx=0);
    
    /// draw `text` at position `pos`
    void drawText(Vector2 const& pos, FontType, const char text[], float dx=0);
    
    /// draw `text` at position `pos`
    void drawText(Vector3 const& pos, FontType, const char text[], float dx=0);

    /// stroke character
    void strokeCharacter(FontType, int character);

}
