// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector1.h"
#include "vector2.h"
#include "vector3.h"

// construct from Vector1
Vector3::Vector3(const Vector1& arg) : XX(arg.XX), YY(0.0), ZZ(0.0)
{
#if VECTOR3_USES_AVX
    vec[3] = 0;
#endif
}

// construct from Vector2
Vector3::Vector3(const Vector2& arg) : XX(arg.XX), YY(arg.YY), ZZ(0.0)
{
#if VECTOR3_USES_AVX
    vec[3] = 0;
#endif
}

/**
 This accepts 'X Y Z' but also 'X' and 'X Y'.
 At least one scalar must be read to be valid
 */
std::istream& operator >> (std::istream& is, Vector3& v)
{
#if VECTOR3_USES_AVX
    v.vec[3] = 0;
#endif
    if ( is >> v.XX )
    {
        if ( is >> v.YY )
        {
            if ( is >> v.ZZ )
                ;
            else
            {
                v.ZZ = 0;
                is.clear();
            }
        }
        else
        {
            v.YY = 0;
            v.ZZ = 0;
            is.clear();
        }
    }
    return is;
}


std::ostream& operator << (std::ostream& os, Vector3 const& v)
{
    int w = (int)os.width();
    os << std::setw(w) << std::showpos << v.XX << " ";
    os << std::setw(w) << v.YY << " ";
    os << std::setw(w) << v.ZZ << std::noshowpos;
    return os;
}

