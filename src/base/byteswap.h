// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <stdint.h>

/**
 Could use the Intel SIMD function _bswap() and _bswap64() here
 */

template <class T>
void byteswap(T) = delete;

/// reverse byte order
inline static uint16_t byteswap(uint16_t i)
{
    //return __builtin_bswap16(i);
    return ( i & 0xff ) << 8 | ( i & 0xff00 ) >> 8;
}

/// reverse byte order
inline static int16_t byteswap(int16_t i)
{
    return (int16_t)byteswap((uint16_t)i);
}

/// reverse byte order
inline static uint32_t byteswap(uint32_t i)
{
    //return __builtin_bswap32(i);
    return (i & 0xff) << 24 |
    (i & 0x0000ff00) << 8 |
    (i & 0x00ff0000) >> 8 |
    (i & 0xff000000) >> 24;
}

/// reverse byte order
inline static int32_t byteswap(int32_t i)
{
    return (int32_t)byteswap((uint32_t)i);
}

/// reverse byte order
inline static uint64_t byteswap(uint64_t i)
{
    //return __builtin_bswap64(i);
    return (i & 0xff) << 56 |
    (i & 0x000000000000ff00) << 40 |
    (i & 0x0000000000ff0000) << 24 |
    (i & 0x00000000ff000000) <<  8 |
    (i & 0x000000ff00000000) >>  8 |
    (i & 0x0000ff0000000000) >> 24 |
    (i & 0x00ff000000000000) >> 40 |
    (i & 0xff00000000000000) >> 56;
}

/// reverse byte order of float
inline static float byteswap(float& x)
{
    uint32_t i = byteswap(reinterpret_cast<uint32_t&>(x));
    return reinterpret_cast<float&>(i);
}

/// reverse byte order of double
inline static double byteswap(double& x)
{
    uint64_t i = byteswap(reinterpret_cast<uint64_t&>(x));
    return reinterpret_cast<double&>(i);
}
