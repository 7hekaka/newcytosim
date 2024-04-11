// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector3.h"
#include "vector2.h"
#include "vector1.h"

// construct from Vector1
Vector3::Vector3(const Vector1& arg) : XX(arg.XX), YY(0.0), ZZ(0.0)
{
#if VECTOR3_USES_AVX
    xyz[3] = 0;
#endif
}

// construct from Vector2
Vector3::Vector3(const Vector2& arg) : XX(arg.XX), YY(arg.YY), ZZ(0.0)
{
#if VECTOR3_USES_AVX
    xyz[3] = 0;
#endif
}

/**
 This accepts 'X Y Z' but also 'X' and 'X Y'.
 If no scalar is read the stream state is set to 'fail'
 */
std::istream& operator >> (std::istream& is, Vector3& v)
{
#if VECTOR3_USES_AVX
    v.xyz[3] = 0;
#endif
    if ( is >> v.XX )
    {
        if ( is.eof() )
        {
            v.YY = 0;
            v.ZZ = 0;
            return is;
        }
        std::streampos isp = is.tellg();
        if ( is >> v.YY )
        {
            if ( is.eof() )
            {
                v.ZZ = 0;
                return is;
            }
            isp = is.tellg();
            if ( is >> v.ZZ )
                ;
            else
            {
                v.ZZ = 0;
                is.clear();
                is.seekg(isp);
            }
        }
        else
        {
            v.YY = 0;
            v.ZZ = 0;
            is.clear();
            is.seekg(isp);
        }
    }
    return is;
}

