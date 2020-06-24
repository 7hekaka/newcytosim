// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// Wednesday 24 June 2020 was a very nice day in Strasbourg

#ifndef SIMD_FLOAT_H
#define SIMD_FLOAT_H

//--------------------------- SSE Single Precision -----------------------------

#ifdef __SSE3__

#include <pmmintrin.h>

/// Vector of 4 floats
typedef __m128 vec4f;

inline vec4f load4f(float const* a)      { return _mm_load_ps(a); }
inline void store4f(float* a, vec4f b)   { return _mm_store_ps(a, b); }
inline vec4f max4f(vec4f a, vec4f b)     { return _mm_max_ps(a,b); }
inline vec4f min4f(vec4f a, vec4f b)     { return _mm_min_ps(a,b); }
inline vec4f and4f(vec4f a, vec4f b)     { return _mm_and_ps(a,b); }
inline vec4f andnot4f(vec4f a, vec4f b)  { return _mm_andnot_ps(a,b); }
inline vec4f abs4f(vec4f a)              { return _mm_andnot_ps(_mm_set1_ps(-0.0), a); }
#define permute4f(a,b)       _mm_permute_ps(a,b)       // same as shuffle2(a,a,b)

#endif

//-------------------------- AVX Single Precision-------------------------------

#ifdef __AVX__

/// Vector of 8 floats
typedef __m256 vec8f;

/// Vector of doubles
typedef __m256d vec4;

inline vec8f load8f(float const* a)     { return _mm256_load_ps(a); }
inline void store8f(float* a, vec8f b)  { return _mm256_store_ps(a, b); }
inline void store4f(float* a, vec4 b)   { return _mm_store_ps(a, _mm256_cvtpd_ps(b)); }

inline vec8f setzero8f()                { return _mm256_setzero_ps(); }
inline vec8f set8f(float const& a)      { return _mm256_set1_ps(a); }
inline vec8f mulf(vec8f a, vec8f b)     { return _mm256_mul_ps(a,b); }
inline vec8f divf(vec8f a, vec8f b)     { return _mm256_div_ps(a,b); }
inline vec8f addf(vec8f a, vec8f b)     { return _mm256_add_ps(a,b); }
inline vec8f subf(vec8f a, vec8f b)     { return _mm256_sub_ps(a,b); }

inline vec8f max8f(vec8f a, vec8f b)    { return _mm256_max_ps(a,b); }
inline vec8f min8f(vec8f a, vec8f b)    { return _mm256_min_ps(a,b); }
inline vec8f and8f(vec8f a, vec8f b)    { return _mm256_and_ps(a,b); }
inline vec8f andnot8f(vec8f a, vec8f b) { return _mm256_andnot_ps(a,b); }
inline vec8f abs8f(vec8f a)             { return _mm256_andnot_ps(_mm256_set1_ps(-0.0), a); }

#define permute8f128(a,b,c)  _mm256_permute4f128_ps(a,b,c)

/// approximate inverse
inline vec8f rcpf(vec8f a)              { return _mm256_rcp_ps(a); }
/// approximate inverse square root
inline vec8f rsqrtf(vec8f a)            { return _mm256_rsqrt_ps(a); }

inline vec4  cvt4sd(vec4f a)            { return _mm256_cvtps_pd(a); }
inline vec4f cvt4ds(vec4 a)             { return _mm256_cvtpd_ps(a); }
inline vec4f getlof(vec8f a)            { return _mm256_castps256_ps128(a); }
inline vec4f gethif(vec8f a)            { return _mm256_extractf128_ps(a,1); }

inline vec8f cvt8i(__m256i a)           { return _mm256_cvtepi32_ps(a); }

#define load8si(a)           _mm256_load_si256(a)
#define cmpf(a,b,c)          _mm256_cmp_ps(a,b,c)
#define permute2f128f(a,b,c) _mm256_permute2f128_ps(a,b,c)

#endif  // AVX

#endif
