// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <cstdint>

template <typename T> void byteswap16(T) = delete;
template <typename T> void byteswap32(T) = delete;
template <typename T> void byteswap64(T) = delete;

/// reverse byte order
static inline uint16_t byteswap16(uint16_t i)
{
    return ( i & 0xff ) << 8 | ( i & 0xff00 ) >> 8;
}

/// reverse byte order
static inline int16_t byteswap16(int16_t i)
{
    return (int16_t)byteswap16((uint16_t)i);
}

/// reverse byte order
static inline uint32_t byteswap32(uint32_t i)
{
    return (i & 0xff) << 24 |
    (i & 0x0000ff00) << 8 |
    (i & 0x00ff0000) >> 8 |
    (i & 0xff000000) >> 24;
}

/// reverse byte order
static inline int32_t byteswap32(int32_t i)
{
    return (int32_t)byteswap32((uint32_t)i);
}

/// reverse byte order
static inline uint64_t byteswap64(uint64_t i)
{
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
static inline float byteswap32(float& x)
{
    uint32_t i = byteswap32(reinterpret_cast<uint32_t&>(x));
    return reinterpret_cast<float&>(i);
}

/// reverse byte order of double
static inline double byteswap64(double& x)
{
    uint64_t i = byteswap64(reinterpret_cast<uint64_t&>(x));
    return reinterpret_cast<double&>(i);
}
