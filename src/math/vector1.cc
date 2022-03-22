// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector1.h"

/**
 Eat up to two scalar values if they are equal to zero
 */
static void eatTwoZeros(std::istream& is)
{
    if ( is.good() )
    {
        std::streampos isp = is.tellg();
        real Z = 0;
        is >> Z;
        if ( is.fail() || Z != 0 )
        {
            is.seekg(isp);
            is.clear();
        }
        isp = is.tellg();
        is >> Z;
        if ( is.fail() || Z != 0 )
        {
            is.seekg(isp);
            is.clear();
        }
    }
}


/**
 This accepts 'X 0 0' but also 'X' and 'X 0'.
 At least one scalar must be read to be valid
 */
std::istream& operator >> (std::istream& is, Vector1& v)
{
    if ( is >> v.XX )
        eatTwoZeros(is);
    return is;
}
