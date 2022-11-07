// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_image.h"
#include <string.h>

/// sign-extension: copying the sign bit of `arg` to all bits
inline static int8_t sex(const int8_t& arg)
{
    union {
        int16_t i;
        struct { int8_t l, h; };
    } u { arg };
    return u.h;
}


/// promote 0/1 bits values to 1 byte per bit
void gym::unpackBitmap(unsigned char * input, unsigned W, unsigned H, const unsigned char* bits)
{
    const unsigned Wb = ( ( W + 7 ) & ~7 ) / 8;
    // sign-extension is only possible on signed types:
    int8_t * data = (int8_t*)input;
    // each row is independent:
    for ( unsigned i = 0; i < H; ++i )
    {
        int8_t * ptr = data + i * W;
        int8_t const* row = (int8_t*)bits + i * Wb;
        // the operation should be vectorized, processing 16 bytes at a time
        for ( unsigned j = 0; j < Wb; ++j )
        {
            int8_t b = row[j];
            ptr[0] = sex( b << 7 );
            ptr[1] = sex( b << 6 );
            ptr[2] = sex( b << 5 );
            ptr[3] = sex( b << 4 );
            ptr[4] = sex( b << 3 );
            ptr[5] = sex( b << 2 );
            ptr[6] = sex( b << 1 );
            ptr[7] = sex( b << 0 );
            ptr += 8;
        }
    }
}


/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `dst` should be of size `4*W*H` with 4 bytes per pixels: R, G, B and A,
 while `src` should be `bin*bin` times larger. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 */
void gym::downsampleRGBA(uint8_t dst[], unsigned W, unsigned H,
                         uint8_t const src[], unsigned bin)
{
    const size_t BB = bin * bin;

#if ( 0 )
    //reset destination:
    for ( size_t u = 0; u < W*H; ++u )
    {
        dst[4*u  ] = 0xFF;
        dst[4*u+1] = 0xFF;
        dst[4*u+2] = 0xFF;
        dst[4*u+3] = 0xFF;
    }
#endif
    
    for ( unsigned x = 0; x < W; ++x )
    for ( unsigned y = 0; y < H; ++y )
    {
        uint8_t const* ptr = src + 4 * bin * ( x + bin*W*y );
        size_t r = 0, g = 0, b = 0, a = 0;
        for ( unsigned dx = 0; dx < bin; ++dx )
        for ( unsigned dy = 0; dy < bin; ++dy )
        {
            uint8_t const* p = ptr + 4 * ( dx + bin*W*dy );
            r += p[0];
            g += p[1];
            b += p[2];
            a += p[3];
        }
            
        dst[4*(x+W*y)  ] = (uint8_t)( r / BB );
        dst[4*(x+W*y)+1] = (uint8_t)( g / BB );
        dst[4*(x+W*y)+2] = (uint8_t)( b / BB );
        dst[4*(x+W*y)+3] = (uint8_t)( a / BB );
    }
}


/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `src` should be of size `3*W*H` with 3 bytes per pixels: R, G, B,
 while `src` will be `bin*bin` times smaller. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 */
void gym::downsampleRGB(uint8_t dst[], unsigned W, unsigned H,
                        const uint8_t src[], unsigned bin)
{
    const size_t BB = bin * bin;
#if ( 0 )
    //reset destination:
    for ( size_t u = 0; u < W*H; ++u )
    {
        dst[3*u  ] = 0xFF;
        dst[3*u+1] = 0xFF;
        dst[3*u+2] = 0xFF;
    }
#endif
    
    for ( size_t x = 0; x < W; ++x )
    for ( size_t y = 0; y < H; ++y )
    {
        uint8_t const* ptr = src + 3 * bin * ( x + bin*W*y );
        size_t r = 0, g = 0, b = 0;
        for ( size_t dx = 0; dx < bin; ++dx )
        for ( size_t dy = 0; dy < bin; ++dy )
        {
            uint8_t const* p = ptr + 3 * ( dx + bin*W*dy );
            r += p[0];
            g += p[1];
            b += p[2];
        }
        
        dst[3*(x+W*y)  ] = (uint8_t)( r / BB );
        dst[3*(x+W*y)+1] = (uint8_t)( g / BB );
        dst[3*(x+W*y)+2] = (uint8_t)( b / BB );
    }
}


/// print pixel map in ASCII
void gym::printPixels(FILE* f, uint8_t const* pix, unsigned W, unsigned H)
{
    for ( size_t y = 0; y < H; ++y )
    {
        for ( size_t x = 0; x < W; ++x )
        {
            uint8_t const* p = pix + 4 * (x+W*y);
            uint16_t lum = ( p[0] + p[1] + p[2] ) * p[3] / ( 3 * 255 );
            fprintf(f, "%02X", lum);
        }
        fprintf(f, "\n");
    }
}

