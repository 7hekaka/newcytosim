// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

class gle_color;
class Vector1;
class Vector2;
class Vector3;

namespace gym
{
    /// return height in pixel of GLUT font
    int fontHeight(void* font);
    
    /// compute size of text
    int maxTextWidth(const char text[], void* font, int& lines);
    
    /// draw text at the current OpenGL raster position and raster color
    void bitmapString(const char text[], void* font = nullptr, GLfloat vshift=0);

    /// display text on a rectangle of color `bcol`, in a corner of the center of the display window
    void drawText(const char text[], void* font, gle_color bcol, int position, int width, int height);
    
    /// draw `text` at position `pos`
    void drawText(Vector1 const& pos, const char text[], void* font, float dx=0);
    
    /// draw `text` at position `pos`
    void drawText(Vector2 const& pos, const char text[], void* font, float dx=0);
    
    /// draw `text` at position `pos`
    void drawText(Vector3 const& pos, const char text[], void* font, float dx=0);

}
