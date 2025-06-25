// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#ifndef VECPRINT_H
#define VECPRINT_H

#include <cstdio>
#include <cmath>
#include <string>
#include <cfloat>
#include <string.h>


/// Templated functions to print Vectors and Matrices with minimal formatting
/** This is mostly used for debugging */
namespace VecPrint
{
    constexpr size_t SUP = 12UL;
    /// default destination:
    #define OUT stderr
    
    /// print 'len' components of 'vec[]' on a line
    template < typename T >
    void print(FILE* file, size_t len, const T* vec, int digits = 2, size_t dim = 0)
    {
        if ( !vec )
            fprintf(file, " null");
        else if ( len == 0 )
            fprintf(file, " void");
        else
        {
            char fmt[32];
            snprintf(fmt, sizeof(fmt), " %%%i.%if", digits+3, digits);
            for ( size_t i = 0; i < len; ++i )
            {
                fprintf(file, fmt, vec[i]);
                if ( dim && ( 0 == (i+1) % dim ))
                    putc(39, file);
            }
        }
    }

    /// print vector to OUT
    template < typename T >
    void print(size_t len, const T* vec, int digits = 2, size_t dim = 0)
    {
        print(OUT, len, vec, digits, dim);
    }
    
    /// print vector to OUT
    template < typename T >
    void print(char const* msg, size_t len, const T* vec, int digits = 2, size_t dim = 0)
    {
        fprintf(OUT, "%6s ", msg);
        print(OUT, len, vec, digits, dim);
        putc('\n', OUT);
    }

    /// print up to 16 scalars from given vector, from the start
    template < typename T >
    void head(size_t len, const T* vec, int digits)
    {
        if ( len <= SUP )
            print(len, vec, digits);
        else
        {
            print(std::min(SUP, len), vec, digits);
            printf("...");
        }
    }
    
    template < typename T >
    void head(char const* msg, size_t len, const T* vec, int digits = 2)
    {
        fprintf(OUT, "%6s ", msg);
        head(len, vec, digits);
        putc('\n', OUT);
    }

    /// print up to 16 scalars from given vector, taken from the edges
    template < typename T >
    void edges(size_t len, const T* vec, int digits = 2)
    {
        if ( len <= SUP )
            print(len, vec, digits);
        else
        {
            print(SUP/2, vec, digits);
            fprintf(OUT, "...");
            print(SUP/2, vec+len-SUP/2, digits);
        }
    }
    
    template < typename T >
    void edges(char const* msg, size_t len, const T* vec, int digits = 2)
    {
        fprintf(OUT, "%6s ", msg);
        edges(len, vec, digits);
        putc('\n', OUT);
    }

    template < typename T >
    void print(std::string const& msg, size_t len, const T* vec, int digits = 2)
    {
        fprintf(OUT, "%s %6zu ", msg.c_str(), len);
        print(OUT, std::min(len, SUP), vec, digits);
        putc('\n', OUT);
    }
    
    /// print 'len' components of 'alpha * vec[]' on a line
    template < typename T >
    void print(FILE* file, size_t len, const T* vec, int digits, T alpha)
    {
        if ( !vec )
            fprintf(file, " null");
        else if ( len == 0 )
            fprintf(file, " void");
        else
        {
            char fmt[32];
            snprintf(fmt, sizeof(fmt), " %%%i.%if", digits+3, digits);
            for ( size_t i = 0; i < len; ++i )
            {
                fprintf(file, fmt, alpha*vec[i]);
                if ( 0 == ( i % 3 ) )
                    putc(' ', file);
            }
        }
        fflush(file);
    }


    /// print 'len' components of 'vec[]' on separate lines
    template < typename T >
    void dump(FILE* file, size_t len, const T* vec, int digits = 8)
    {
        if ( !vec )
            fprintf(file, " null");
        else if ( len == 0  )
            fprintf(file, " void");
        else
        {
            char fmt[32];
            snprintf(fmt, sizeof(fmt), " %%%i.%ie\n", 9, digits);
            for ( size_t i = 0; i < len; ++i )
                fprintf(file, fmt, vec[i]);
        }
        fflush(file);
    }
    
    
    /// print matrix `mat[]` of size 'lin*col', and leading dimension `ldd` with precision 'digits'
    template < typename T >
    void full(FILE* file, size_t lin, size_t col, const T* mat, size_t ldd, int digits = 2)
    {
        if ( !mat )
            fprintf(file, " null");
        else if ( lin+col == 0  )
            fprintf(file, " void");
        {
            const T threshold = std::pow(0.1, digits);
            char zer[32] = { 0 }, eps[32] = { 0 };
            char fmt[32] = " %4.0f";
            
            { // build format strings:
                snprintf(fmt, sizeof(fmt), " %%%i.%if", digits+4, digits);
                snprintf(zer, sizeof(zer), fmt, 0.0);
                bool dot = false; char * d = zer;
                for ( char * c = zer; *c; ++c )
                {
                    if ( *c == '0' ) { *c = ' '; d = c; }
                    dot |= ( *c == '.' );
                }
                if ( !dot ) *d = '.';
                memcpy(eps, zer, sizeof(eps));
                for ( char * c = eps; *c; ++c )
                    if ( *c == '.' ) *c = '*';
            }
            
            for ( size_t ii = 0; ii < lin; ++ii )
            {
                for ( size_t jj = 0; jj < col; ++jj )
                {
                    T val = mat[ii+ldd*jj];
                    if ( std::fabs(val) < FLT_EPSILON )
                        fputs(zer, file);
                    else if ( std::fabs(val) < threshold )
                        fputs(eps, file);
                    else
                        fprintf(file, fmt, mat[ii+ldd*jj]);
                }
                putc('\n', file);
            }
        }
        fflush(file);
    }
    
    
    template < typename T >
    void full(size_t lin, size_t col, const T* mat, size_t ldd, int digits = 2)
    {
        full(OUT, lin, col, mat, ldd, digits);
    }
    
    template < typename T >
    void full(std::string const& msg, size_t lin, size_t col, const T* mat, size_t ldd, int digits = 2)
    {
        fprintf(OUT, "%s %lux%lu :\n", msg.c_str(), lin, col);
        full(OUT, std::min(lin, SUP), std::min(col, SUP), mat, ldd, digits);
    }

    /// print matrix in sparse format: line_index, column_index, value
    template < typename T >
    void sparse(FILE * file, size_t lin, size_t col, const T* mat, size_t ldd, int digits = 8, T threshold = 0)
    {
        if ( !mat )
            fprintf(file, " null");
        else if ( lin+col == 0  )
            fprintf(file, " void");
        else
        {
            char fmt[64];
            snprintf(fmt, sizeof(fmt), " %%3i %%3i %%9.%if\n", digits);
            for (size_t ii = 0; ii < lin; ++ii )
                for (size_t jj = 0; jj < col; ++jj )
                {
                    T val = mat[ii+ldd*jj];
                    if ( std::fabs(val) > threshold )
                        fprintf(file, fmt, ii, jj, val);
                }
        }
        putc('\n', file);
    }
    
    
    /// print a matrix in sparse format, but adding `off` to all line and column indices
    template < typename T >
    void sparse_off(FILE * file, size_t lin, size_t col, const T* mat, size_t ldd, size_t off, int digits = 8)
    {
        if ( !mat )
            fprintf(file, " null");
        else if ( lin+col == 0  )
            fprintf(file, " void");
        else
        {
            char fmt[32];
            snprintf(fmt, sizeof(fmt), " %%3i %%3i  %%9.%if\n", digits);
            for (size_t ii = 0; ii < lin; ++ii )
                for (size_t jj = 0; jj < col; ++jj )
                    fprintf(file, fmt, ii+off, jj+off, mat[ii+ldd*jj]);
        }
        putc('\n', file);
    }
    
    /// print matrix `mat[]` of size lin*col, and leading dimension `ldd` in ASCII art...
    template < typename T >
    void image(FILE * file, size_t lin, size_t col, const T* mat, size_t ldd, T scale)
    {
        if ( !mat )
            fprintf(file, " null");
        else if ( lin+col == 0  )
            fprintf(file, " void");
        else
        {
            char map[] = ".:+*O%#$";
            
            const T threshold = 0.01 * scale;
            for ( size_t ii = 0; ii < lin; ++ii )
            {
                putc('|', file);
                for ( size_t jj = 0; jj < col; ++jj )
                {
                    T val = mat[ii+ldd*jj];
                    if ( val != val )
                        putc('@', file);
                    else if ( val < threshold )
                        putc(' ', file);
                    else
                    {
                        int x = std::max(7, 2+std::log10(std::fabs(val)/scale));
                        putc(map[x], file);
                    }
                }
                putc('|', file);
                putc('\n', file);
            }
        }
        putc('\n', file);
    }
}

#endif
