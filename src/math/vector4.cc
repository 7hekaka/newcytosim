// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector4.h"


/**
 This accepts 'X Y Z' but also 'X' and 'X Y'.
 At least one scalar must be read to be valid
 */
std::istream& operator >> (std::istream& is, Vector4& v)
{
    if ( is >> v.XX )
    {
        if ( is >> v.YY )
        {
            if ( is >> v.ZZ )
            {
                if ( is >> v.TT )
                    ;
                else
                {
                    v.TT = 0;
                    is.clear();
                }
            }
            else
            {
                v.ZZ = 0;
                v.TT = 0;
                is.clear();
            }
        }
        else
        {
            v.YY = 0;
            v.ZZ = 0;
            v.TT = 0;
            is.clear();
        }
    }
    return is;
}

