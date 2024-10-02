// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cstdio>
#include <cstdint>


/// print 'cnt' components of 'ptr' on a line
static inline void print_real(FILE* out, size_t cnt, real const* ptr, const char prefix[])
{
    if ( !ptr || cnt == 0 )
        fprintf(out, "void");
    else
    {
        fprintf(out, "%s", prefix);
        for ( size_t i = 0; i < cnt; ++i )
            fprintf(out, " %9.4f", ptr[i]);
    }
}


/// check memory alignement of a pointer for AVX load/store
static inline void check_alignment(void const* ptr)
{
    uintptr_t a = ((uintptr_t)ptr & 31);
    if ( a )
        fprintf(stderr, "missaligned pointer %p (%lu)\n", ptr, a);
}


#if ( 0 )

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "tokenizer.h"
/**
 Redefining the standard extraction operator is normally not necessary,
 but this can be useful to track certain bugs
 */
static inline std::istream& operator >> (std::istream& is, real& x)
{
    std::string str = Tokenizer::get_real(is);
    if ( str.empty() )
        return is;
    try {
        x = std::stod(str);
        //std::clog << " custom operator >> |" << str << "| -> " << x << '\n';
    }
    catch ( std::invalid_argument & e ) {
        std::cerr << " error in operator >> |" << str << "| " << '\n';
    }
    return is;
}

#endif
#if ( 0 )

/// return the usual base-10 representation of a number
/** This is equivalent to std::to_string from the C++11 standard */
template <typename T>
std::string to_string(T const& x)
{
    std::ostringstream oss;
    oss << x;
    return oss.str();
}


template <typename T>
std::string to_string(T const& x, int width, unsigned precision)
{
    std::ostringstream oss;
    oss.precision(precision);
    oss << std::setw(width) << std::fixed << x;
    return oss.str();
}


#endif
