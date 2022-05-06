// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University


namespace gym
{
    /// draw Bitmap using a texture over a rectangle
    void drawPixels(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits);

    /// draw Bitmap using a texture over a rectangle
    void drawBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits);
    
    /// draw all pixels of the bitmap
    void paintBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits);
    
    /// draw all pixels of the bitmap
    void paintSequence(float X, float Y, float W, float H, const char str[]);

    /// draw rectangle ( rect = [ left, bottom, right, top ] )
    void drawRectangle(const int rect[4], float width);

    /// draw rectangle ( left, bottom, right, top )
    void drawRectangle(float L, float B, float R, float T, float Z, float width);
    
    /// paint rectangle ( left, bottom, right, top )
    void fillRectangle(float L, float B, float R, float T, float Z, const float col[4]);
    
    /// draw plane with squares of alternating colors
    void drawTiledFloor(int R, float T, float Z, const float col1[], const float col2[]);
    
    /// draw Octagon
    void paintOctagon(float L, float B, float R, float T, const float col[4], float D);

    /// paint octagon
    void drawOctagon(float L, float B, float R, float T, const float col[4], float D, float W);

}
