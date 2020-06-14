// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <string>
#include <iostream>
#include "random.h"



/**
 Returns a bitwise representation of the argument 'val'
 */
template < typename T >
std::string to_bits(const T& val, char spc = 0)
{
    unsigned char * ptr = (unsigned char*) & val;
    char res[10*sizeof(T)];
    char * out = res;
    
    for ( int i = sizeof(T)-1; i >= 0; --i)
    {
        unsigned char byte = ptr[i];
        for ( int i = 0; i < 8; ++i )
        {
            *out++ = (byte&0x80?'1':'0');
            byte <<= 1;
        }
        if ( spc )
            *out++ = spc;
    }
    *out = 0;
    return std::string(res);
}



int main(int argc, char* argv[])
{
    RNG.seed();
    
    for ( int i = 0; i < 16; ++i )
        std::cout << to_bits(i, ' ') << '\n';
    
    printf("\ndone\n");
    return 0;
}
