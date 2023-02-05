// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include <cstdio>
#include <stdint.h>

namespace gym
{

    /// reduce image size by 'bin'
    void downsampleRGB(uint8_t dst[], unsigned W, unsigned H, uint8_t const src[], unsigned bin);

    /// reduce image size by 'bin'
    void downsampleRGBA(uint8_t dst[], unsigned W, unsigned H, uint8_t const src[], unsigned bin);

    /// print pixel map in ASCII
    void printPixels(FILE*, uint8_t const* pix, unsigned W, unsigned H);
    
}

/// convert binary image into one byte per pixel
void unpackBitmap(unsigned char data[], unsigned W, unsigned H, const unsigned char bits[], unsigned);
