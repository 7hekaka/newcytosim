// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector4.h"
#include "vector3.h"
#include "vector2.h"
#include "vector1.h"


// construct from Vector1
Vector4::Vector4(const Vector1& arg) : XX(arg.XX), YY(0.0), ZZ(0.0), TT(0.0)
{
}

// construct from Vector2
Vector4::Vector4(const Vector2& arg) : XX(arg.XX), YY(arg.YY), ZZ(0.0), TT(0.0)
{
}

// construct from Vector3
Vector4::Vector4(const Vector3& arg) : XX(arg.XX), YY(arg.YY), ZZ(arg.ZZ), TT(0.0)
{
}

/**
 This accepts 'X Y Z' but also 'X' and 'X Y'.
 At least one scalar must be read to be valid
 */
std::istream& operator >> (std::istream& is, Vector4& v)
{
    if ( is >> v.XX )
    {
        if ( is.eof() )
        {
            v.YY = 0;
            v.ZZ = 0;
            v.TT = 0;
            return is;
        }
        std::streampos isp = is.tellg();
        if ( is >> v.YY )
        {
            if ( is.eof() )
            {
                v.ZZ = 0;
                v.TT = 0;
                return is;
            }
            isp = is.tellg();
            if ( is >> v.ZZ )
            {
                if ( is.eof() )
                {
                    v.TT = 0;
                    return is;
                }
                isp = is.tellg();
                if ( is >> v.TT )
                    ;
                else
                {
                    v.TT = 0;
                    is.clear();
                    is.seekg(isp);
                }
            }
            else
            {
                v.ZZ = 0;
                v.TT = 0;
                is.clear();
                is.seekg(isp);
            }
        }
        else
        {
            v.YY = 0;
            v.ZZ = 0;
            v.TT = 0;
            is.clear();
            is.seekg(isp);
        }
    }
    return is;
}

