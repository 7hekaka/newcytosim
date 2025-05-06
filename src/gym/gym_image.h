// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include <cstdio>
#include <cstdint>
#include "gym_flute.h"

namespace gym
{

    /// reduce image size by 'bin'
    void downsampleRGB(uint8_t dst[], unsigned W, unsigned H, uint8_t const src[], unsigned bin);

    /// reduce image size by 'bin'
    void downsampleRGBA(uint8_t dst[], unsigned W, unsigned H, uint8_t const src[], unsigned bin);

    /// print pixel map in ASCII
    void printPixels(FILE*, uint8_t const* pix, unsigned W, unsigned H);

    /// convert binary image into one byte per pixel
    void unpackBitmap(unsigned char data[], unsigned W, unsigned H, const unsigned char bits[], unsigned);

    /// unpack binary image into vertex data to be rendered as a triangle strip
    size_t unpackBitmap(flute2*, unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits);

}
