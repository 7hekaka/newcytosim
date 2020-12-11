
#include <stdint.h>

/**
 Can use the Intel SIMD function _bswap() and _bswap64()
 */

template <class T>
void byteswap(T) = delete;

/// reverse byte order
inline static uint16_t byteswap(uint16_t i)
{
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
    return (i & 0xff) << 24 |
    (i & 0xff00) << 8 |
    (i & 0xff0000) >> 8 |
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
    return (i & 0xff) << 56
    | (i & 0xff00) << 40 |
    (i & 0xff0000) << 24 |
    (i & 0xff000000) << 8 |
    (i & 0xff00000000) >> 8 |
    (i & 0xff0000000000) >> 24 |
    (i & 0xff000000000000) >> 40 |
    (i & 0xff00000000000000) >> 56;
}

/// reverse byte order of float
inline static float byteswap(float i)
{
    union { uint32_t i; float f; } tmp;
    tmp.f = i;
    tmp.i = byteswap(tmp.i);
    return tmp.f;
}

/// reverse byte order of double
inline static double byteswap(double i)
{
    union { uint64_t i; double f; } tmp;
    tmp.f = i;
    tmp.i = byteswap(tmp.i);
    return tmp.f;
}
