/// FJ Nedelec, La Foret Fouesnant, 10.08.2020

#include <stdio.h>
#include "arm_neon.h"

typedef float real;
typedef float32x4_t vec4f;

const size_t CNT = 1024;

void zero(real * vec)
{
    for ( int i = 0; i < CNT; ++i )
        vec[i] = 0;
}

//------------------------- PLAIN SCALAR CODE ----------------------------------

void sum(real const* a, real const* b, real* c)
{
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT; ++i )
        c[i] = a[i] + b[i];
}

real dot(real const* a, real const* b)
{
    real sum = 0;
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT; ++i )
        sum += a[i] * b[i];
    return sum;
}

void fmadd(real const* a, real const* b, real* c)
{
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT; ++i )
        c[i] = a[i] * b[i] + c[i];
}

real sum(real const* vec)
{
    real sum = 0;
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT; ++i )
        sum = sum + vec[i];
    return sum;
}

//--------------------------- ARM's SIMD NEON ----------------------------------

#ifdef __ARM_NEON__

void sum(vec4f const* a, vec4f const* b, vec4f* c)
{
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT/4; ++i )
        c[i] = vaddq_f32(a[i], b[i]);
}

real dot(vec4f const* a, vec4f const* b)
{
    vec4f sum = { 0 };
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT/4; ++i )
        sum = vfmaq_f32(sum, a[i], b[i]);
    return sum[0] + sum[1] + sum[2] + sum[3];
}

void fmadd(vec4f const* a, vec4f const* b, vec4f* c)
{
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT/4; ++i )
        c[i] = vfmaq_f32(c[i], a[i], b[i]);
}

real sum(vec4f const* vec)
{
    vec4f sum = { 0 };
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT/4; ++i )
        sum = vaddq_f32(sum, vec[i]);
    return sum[0] + sum[1] + sum[2] + sum[3];
}

#endif

//--------------------------------- MAIN ---------------------------------------

int main(int argc, char* argv[])
{
    real a[CNT] = { 1 };
    real b[CNT] = { 0.5 };
    real c[CNT] = { 0 };
    
    fmadd(a, b, c);
    real s = sum(c);
    printf("scalar %f  %f\n", c[0], s);

#ifdef __ARM_NEON__
    zero(c);
    fmadd((vec4f*)a, (vec4f*)b, (vec4f*)c);
    s = sum((vec4f*)c);
    printf("neon   %f  %f\n", c[0], s);
#endif
}
