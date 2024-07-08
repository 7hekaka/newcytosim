// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_image.h"


/// promote 0/1 bits values to 1 byte per bit, either full 1 or full 0
void gym::unpackBitmap(unsigned char * bytes, unsigned W, unsigned H, const unsigned char* bits, unsigned lda)
{
    const unsigned char ONE = 0xFF;
    // number of bytes in a row of 'input'
    const unsigned Wb = ( W + 7 ) >> 3;
    // each line is independent:
    for ( unsigned i = 0; i < H; ++i )
    {
        unsigned char const* row = bits + i * Wb;
        unsigned char * dst = bytes + ( H-1 - i ) * lda;
        // the operation could be vectorized, processing 16 bytes at a time
        for ( unsigned j = 0; j < Wb; ++j )
        {
            unsigned char b = row[j];
            dst[0] = ( b & 128 ) ? ONE : 0;
            dst[1] = ( b & 64  ) ? ONE : 0;
            dst[2] = ( b & 32  ) ? ONE : 0;
            dst[3] = ( b & 16  ) ? ONE : 0;
            dst[4] = ( b & 8   ) ? ONE : 0;
            dst[5] = ( b & 4   ) ? ONE : 0;
            dst[6] = ( b & 2   ) ? ONE : 0;
            dst[7] = ( b & 1   ) ? ONE : 0;
            dst += 8;
        }
    }
}


#if 1
/** Transform bitmap into a triangle strip */
size_t gym::unpackBitmap(flute2* flu, unsigned W, unsigned H, float X0, float Y0, float S, const unsigned char* bits)
{
    const unsigned Wb = ( W + 7 ) >> 3;
    flute2* ptr = flu;
    for ( unsigned i = 0; i < H; ++i )
    {
        float X = X0;
        float Y = Y0 + S * i, T = Y + S;
        const unsigned char* row = bits + i * Wb;
        unsigned byte = row[0];
        for ( unsigned b = 1; b < Wb; ++b )
            byte = ( byte << 8 ) | row[b];
#if 1
        // join next rows if identical to current one:
        while ( i+1 < H  &&  bits[i*Wb+Wb] == *row )
        {
            unsigned next = row[Wb];
            for ( unsigned b = 1; b < Wb; ++b )
                next = ( next << 8 ) | row[Wb+b];
            if ( byte != next ) break;
            T += S;
            ++i;
        }
#endif
        byte <<= 8 * ( sizeof(unsigned) - Wb );
        unsigned p, k = 0;
        while ( byte )
        {
            // find next '1':
            p = __builtin_clz(byte);
            k += p;
            byte <<= p;
            X = X0 + k * S;
            ptr[0] = { X, Y };
            ptr[1] = { X, Y };
            ptr[2] = { X, T };
            ptr += 3;
            // find next '0':
            p = __builtin_clz(~byte);
            k += p;
            byte <<= p;
            X = X0 + k * S;
            ptr[0] = { X, Y };
            ptr[1] = { X, T };
            ptr[2] = { X, T };
            ptr += 3;
        }
    }
    return ptr-flu;
}

#else
/** Transform bitmap into a triangle strip */
size_t gym::unpackBitmap(flute2* flu, unsigned W, unsigned H, float X0, float Y0, float S, const unsigned char* bits)
{
    const unsigned Wb = ( W + 7 ) >> 3;
    flute2* ptr = flu;
    for ( unsigned j = 0; j < Wb; ++j )
    {
        float X = X0, R = X + S;
        unsigned char old = 0;
        for ( int k = 7; k >= 0; --k )
        {
            float Y = Y0;
            for ( unsigned i = 0; i < H; ++i )
            {
                unsigned char bit = ( bits[i*Wb+j] >> k ) & 1;
                if ( bit != old )
                {
                    old = bit;
                    ptr[0] = { X, Y };
                    ptr[1] = { (bit?X:R), Y };
                    ptr[2] = { R, Y };
                    ptr += 3;
                }
                Y += S;
            }
            if ( old )
            {
                ptr[0] = { X, Y };
                ptr[1] = { R, Y };
                ptr[2] = { R, Y };
                ptr += 3;
            }
            X += S;
            R += S;
        }
    }
    return ptr-flu;
}
#endif

/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `dst` should be of size `4*W*H` with 4 bytes per pixels: R, G, B and A,
 while `src` should be `bin*bin` times larger. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 Note that 'dst' may be equal to 'src'.
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
    
    for ( unsigned y = 0; y < H; ++y )
    for ( unsigned x = 0; x < W; ++x )
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
 Note that 'dst' may be equal to 'src'.
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
    
    for ( size_t y = 0; y < H; ++y )
    for ( size_t x = 0; x < W; ++x )
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

