// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// branchless code with sign extension;  FJN 14.06.2020
// https://en.wikipedia.org/wiki/Sign_extension

#include <string>
#include <iostream>
#include <climits>
#include <cstdint>
#include "random.h"


/**
 Returns a bitwise representation of the argument 'val'.
 */
template < typename T >
std::string to_bits(const T& val, char spc = 0)
{
    unsigned char * ptr = (unsigned char*) & val;
    char res[1+sizeof(T)*(CHAR_BIT+1)] = {0};
    char * out = res + CHAR_BIT * sizeof(T) + (spc?sizeof(T):0) - 1;
    
    for ( size_t i = 0; i < sizeof(T); ++i )
    {
        unsigned char byte = ptr[i];
        for ( size_t j = 0; j < CHAR_BIT; ++j )
        {
            *out-- = '0' + ( 1 & byte );
            byte >>= 1;
        }
        if ( spc )
            *out-- = spc;
    }
    return std::string(res);
}



template < typename T >
inline T sox(const T& arg)
{
    return arg >> (CHAR_BIT*sizeof(T)-1);
}

/// implements sign-extension, copying the sign bit of `arg` to all bits
inline int sex(const int32_t& arg)
{
    union
    {
        int64_t i;
        struct { int32_t l, h; };
    } u { arg };
    return u.h;
}

/// implements sign-extension, copying the sign bit of `arg` to all bits
inline int sex(const int16_t& arg)
{
    union
    {
        int32_t i;
        struct { int16_t l, h; };
    } u { arg };
    return u.h;
}

/// implements sign-extension, copying the sign bit of `arg` to all bits
inline int8_t sex(const int8_t& arg)
{
    union
    {
        int16_t i;
        struct { int8_t l, h; };
    } u { arg };
    return u.h;
}

/// implements sign-extension, copying the sign bit of `arg` to all bits
inline uint32_t sex(const float& arg)
{
    union { float d; int32_t i; } udi { arg };
    return udi.i >> (CHAR_BIT*sizeof(float)-1);
}

/// implements sign-extension, copying the sign bit of `arg` to all bits
inline uint64_t sex(const double& arg)
{
    union { double d; int64_t i; } udi { arg };
    return udi.i >> (CHAR_BIT*sizeof(double)-1);
}


/// returns 'neg' if ( arg < 0 ) and 'pos' otherwise
inline double sign_select_hack(const double& val, const double& neg, const double& pos)
{
    union { double d; int64_t i; } v { val };
    union { double d; int64_t i; } n { neg };
    union { double d; int64_t i; } p { pos };
    // using (a & ~mask) | (b & mask)  = a ^ ((a ^ b) & mask);
    p.i ^= ( p.i ^ n.i ) & ( v.i >> (CHAR_BIT*sizeof(double)-1) );
    return p.d;
}


int main(int argc, char* argv[])
{
    RNG.seed();
    
    std::cout << "true  " << to_bits(true, ' ') << '\n';
    std::cout << "false " << to_bits(false, ' ') << '\n';
    
    std::cout << " 0        " << to_bits(0, ' ') << '\n';
    std::cout << "-1        " << to_bits(-1, ' ') << '\n';
    std::cout << " 7        " << to_bits( 7, ' ') << '\n';
    std::cout << "-7        " << to_bits(-7, ' ') << '\n';
    std::cout << "-255      " << to_bits(-255, ' ') << '\n';
    
    std::cout << "sox(0)    " << to_bits(sox(0), ' ') << '\n';
    std::cout << "sox(7)    " << to_bits(sox(7), ' ') << '\n';
    std::cout << "sox(-7)   " << to_bits(sox(-7), ' ') << '\n';
    std::cout << "sox(-255) " << to_bits(sox(-255), ' ') << '\n';

    std::cout << "-1.0f     " << to_bits(-1.f, ' ') << '\n';
    std::cout << "-0.0f     " << to_bits(-0.f, ' ') << '\n';
    std::cout << " 0.0f     " << to_bits( 0.f, ' ') << '\n';
    std::cout << " 1.0f     " << to_bits( 1.f, ' ') << '\n';
    
    std::cout << "sex(-1.f) " << to_bits(sex(-1.f), ' ') << '\n';
    std::cout << "sex(-0.f) " << to_bits(sex(-0.f), ' ') << '\n';
    std::cout << "sex(+0.f) " << to_bits(sex( 0.f), ' ') << '\n';
    std::cout << "sex(+1.f) " << to_bits(sex( 1.f), ' ') << '\n';

    std::cout << "-1.0      " << to_bits(-1.0, ' ') << '\n';
    std::cout << "-0.0      " << to_bits(-0.0, ' ') << '\n';
    std::cout << " 0.0      " << to_bits( 0.0, ' ') << '\n';
    std::cout << " 1.0      " << to_bits( 1.0, ' ') << '\n';
    
    std::cout << "sex(-1.0) " << to_bits(sex(-1.0), ' ') << '\n';
    std::cout << "sex(-0.0) " << to_bits(sex(-0.0), ' ') << '\n';
    std::cout << "sex(+0.0) " << to_bits(sex( 0.0), ' ') << '\n';
    std::cout << "sex(+1.0) " << to_bits(sex( 1.0), ' ') << '\n';

    /*
    std::cout << " 0.0 " << to_bits( 0.0, ' ') << '\n';
    std::cout << " 1.0 " << to_bits( 1.0, ' ') << '\n';
    std::cout << "-1.0 " << to_bits(-1.0, ' ') << '\n';
    std::cout << "-2.0 " << to_bits(-2.0, ' ') << '\n';

    std::cout << "sex(0.0)  " << to_bits(sex( 0.0), ' ') << '\n';
    std::cout << "sex(1.0)  " << to_bits(sex( 1.0), ' ') << '\n';
    std::cout << "sex(-1.0) " << to_bits(sex(-1.0), ' ') << '\n';
    std::cout << "sex(-2.0) " << to_bits(sex(-2.0), ' ') << '\n';
     */
    
    std::cout << "sign_select(-1.0, 1.0, 2.0)  " << sign_select(-1.0, 1.0, 2.0) << '\n';
    std::cout << "sign_select(-0.0, 1.0, 2.0)  " << sign_select(-0.0, 1.0, 2.0) << '\n';
    std::cout << "sign_select( 0.0, 1.0, 2.0)  " << sign_select( 0.0, 1.0, 2.0) << '\n';
    std::cout << "sign_select( 1.0, 1.0, 2.0)  " << sign_select( 1.0, 1.0, 2.0) << '\n';

    printf("\ndone\n");
}
