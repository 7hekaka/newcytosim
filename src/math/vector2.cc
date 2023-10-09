// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector2.h"

/**
 Eat one scalar value if it is equal to zero
 */
static void eatOneZero(std::istream& is)
{
    std::streampos isp = is.tellg();
    while ( isspace(is.peek()) )
        is.get();
    if ( is.fail() )
    {
        // restore initial state:
        is.seekg(isp);
        is.clear();
    }
    else if ( is.peek() == '0' )
    {
        real Z = 0;
        is >> Z;
        if ( is.fail() || Z != 0 )
        {
            // restore initial state:
            is.seekg(isp);
            is.clear();
        }
    }
}


/**
 This accepts 'X Y 0' but also 'X' and 'X Y'.
 At least one scalar must be read to be valid
 */
std::istream& operator >> (std::istream& is, Vector2& v)
{
    if ( is >> v.XX )
    {
        if ( is.eof() )
        {
            v.YY = 0;
            return is;
        }
        std::streampos isp = is.tellg();
        if ( is >> v.YY )
            eatOneZero(is);
        else
        {
            v.YY = 0;
            is.clear();
            is.seekg(isp);
        }
    }
    return is;
}

