// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// Wednesday 24 June 2020 was a very nice day in Strasbourg

#ifndef SIMD_FLOAT_H
#define SIMD_FLOAT_H

#include <immintrin.h>

//--------------------------- SSE Single Precision -----------------------------

#if defined(__SSE3__)

/// Vector of 4 floats
typedef __m128 vec4f;

inline static vec4f setzero4f()                   { return _mm_setzero_ps(); }
inline static vec4f set4f(float a)                { return _mm_set1_ps(a); }
inline static vec4f load2f(float const* a)        { return (vec4f)_mm_loadl_epi64((__m128i*)a); }
inline static vec4f load4f(float const* a)        { return _mm_load_ps(a); }
inline static vec4f loadu4f(float const* a)       { return _mm_loadu_ps(a); }
inline static void  store3f(float* a, vec4f b)    { a[0]=b[0]; a[1]=b[1]; a[2]=b[2]; }
inline static void  store4f(float* a, vec4f b)    { _mm_store_ps(a, b); }
inline static void  storeu4f(float* a, vec4f b)   { _mm_storeu_ps(a, b); }
inline static vec4f add4f(vec4f a, vec4f b)       { return _mm_add_ps(a,b); }
inline static vec4f mul4f(vec4f a, vec4f b)       { return _mm_mul_ps(a,b); }
inline static vec4f max4f(vec4f a, vec4f b)       { return _mm_max_ps(a,b); }
inline static vec4f min4f(vec4f a, vec4f b)       { return _mm_min_ps(a,b); }
inline static vec4f and4f(vec4f a, vec4f b)       { return _mm_and_ps(a,b); }
inline static vec4f andnot4f(vec4f a, vec4f b)    { return _mm_andnot_ps(a,b); }
inline static vec4f abs4f(vec4f a)                { return _mm_andnot_ps(_mm_set1_ps(-0.0f), a); }
inline static vec4f unpacklo4f(vec4f a, vec4f b)  { return _mm_unpacklo_ps(a,b); }
inline static vec4f unpackhi4f(vec4f a, vec4f b)  { return _mm_unpackhi_ps(a,b); }
inline static vec4f duplo4f(vec4f a)              { return _mm_unpacklo_ps(a,a); }
inline static vec4f duphi4f(vec4f a)              { return _mm_unpackhi_ps(a,a); }

#define shuffle4f(a,b,k) _mm_shuffle_ps(a,b,k)

inline static vec4f load3f(float const* a) { return _mm_blend_ps(_mm_loadu_ps(a), _mm_setzero_ps(), 0b1000); }


#endif  // __SSE3__


#if defined(__SSE4_1__)

#  define blend4f(a,b,k) _mm_blend_ps(a,b,k)

inline static vec4f sign_select(vec4f const& val, vec4f const& neg, vec4f const& pos)
{
#if defined(__AVX512VL__)
    return _mm_mask_mov_ps(pos, val, neg);
#else
    return _mm_blendv_ps(pos, neg, val);
#endif
}

#endif


#if defined(__AVX__)
// copy a[0] into all elements of destination
inline static vec4f broadcast1f(vec4f a)          { return _mm_permute_ps(a,0x00); }
inline static vec4f broadcast1f(float const* a)   { return _mm_broadcast_ss(a); }
inline static vec4f streamload4f(float const* a)  { return (vec4f)_mm_stream_load_si128((__m128i*)a); }
#define permute4f(a,k)    _mm_permute_ps(a,k)
// Convert between vector types
inline static vec4f cvt4ds(__m256d a)             { return _mm256_cvtpd_ps(a); }
inline static __m256d  cvt4sd(vec4f a)            { return _mm256_cvtps_pd(a); }
inline static void store4f(float* a, __m256d b)   { _mm_store_ps(a, _mm256_cvtpd_ps(b)); }
#elif defined(__SSE3__)
inline static vec4f broadcast1f(vec4f a)          { return _mm_shuffle_ps(a,a,0x00); }
inline static vec4f broadcast1f(float const* a)   { return _mm_load1_ps(a); }
inline static vec4f streamload4f(float const* a)  { return _mm_load_ps(a); }
#define permute4f(a,k)    _mm_shuffle_ps(a,a,k)
#endif

//-------------------------- FMA Single Precision-------------------------------

#if defined(__FMA__)
inline static vec4f fmadd4f (vec4f a, vec4f b, vec4f c) { return _mm_fmadd_ps(a,b,c); }
inline static vec4f fmsub4f (vec4f a, vec4f b, vec4f c) { return _mm_fmsub_ps(a,b,c); }
inline static vec4f fnmadd4f(vec4f a, vec4f b, vec4f c) { return _mm_fnmadd_ps(a,b,c); }
#elif defined(__SSE3__)
// erzatz functions
inline static vec4f fmadd4f (vec4f a, vec4f b, vec4f c) { return _mm_add_ps(_mm_mul_ps(a,b), c); }  // a * b + c
inline static vec4f fmsub4f (vec4f a, vec4f b, vec4f c) { return _mm_sub_ps(_mm_mul_ps(a,b), c); }
inline static vec4f fnmadd4f(vec4f a, vec4f b, vec4f c) { return _mm_sub_ps(c, _mm_mul_ps(a,b)); }
#endif

//-------------------------- AVX Single Precision-------------------------------

#if defined(__AVX__)

/// Vector of 8 floats
typedef __m256 vec8f;

inline static vec8f setzero8f()                   { return _mm256_setzero_ps(); }
inline static vec8f set8f(float const& a)         { return _mm256_set1_ps(a); }
inline static vec8f load8f(float const* a)        { return _mm256_load_ps(a); }
inline static vec8f loadu8f(float const* a)       { return _mm256_loadu_ps(a); }
inline static void  store8f(float* a, vec8f b)    { _mm256_store_ps(a, b); }
inline static void  storeu8f(float* a, vec8f b)   { _mm256_storeu_ps(a, b); }
inline static vec8f mulf(vec8f a, vec8f b)        { return _mm256_mul_ps(a,b); }
inline static vec8f divf(vec8f a, vec8f b)        { return _mm256_div_ps(a,b); }
inline static vec8f addf(vec8f a, vec8f b)        { return _mm256_add_ps(a,b); }
inline static vec8f subf(vec8f a, vec8f b)        { return _mm256_sub_ps(a,b); }

inline static vec8f max8f(vec8f a, vec8f b)       { return _mm256_max_ps(a,b); }
inline static vec8f min8f(vec8f a, vec8f b)       { return _mm256_min_ps(a,b); }
inline static vec8f and8f(vec8f a, vec8f b)       { return _mm256_and_ps(a,b); }
inline static vec8f andnot8f(vec8f a, vec8f b)    { return _mm256_andnot_ps(a,b); }
inline static vec8f abs8f(vec8f a)                { return _mm256_andnot_ps(_mm256_set1_ps(-0.0), a); }

#define permute8f128(a,b,c)  _mm256_permute4f128_ps(a,b,c)

/// approximate inverse
inline static vec8f rcpf(vec8f a)                 { return _mm256_rcp_ps(a); }
/// approximate inverse square root
inline static vec8f rsqrtf(vec8f a)               { return _mm256_rsqrt_ps(a); }

inline static vec4f getlo4f(vec8f a)              { return _mm256_castps256_ps128(a); }
inline static vec4f gethi4f(vec8f a)              { return _mm256_extractf128_ps(a,1); }

inline static vec8f cvt8i(__m256i a)              { return _mm256_cvtepi32_ps(a); }

#define load8si(a)           _mm256_load_si256(a)
#define cmp8f(a,b,c)         _mm256_cmp_ps(a,b,c)
#define permute2f128f(a,b,c) _mm256_permute2f128_ps(a,b,c)

inline static vec8f sign_select(vec8f const& val, vec8f const& neg, vec8f const& pos)
{
#if defined(__AVX512VL__)
    return _mm256_mask_mov_ps(pos, val, neg);
#else
    return _mm256_blendv_ps(pos, neg, val);
#endif
}

#endif  // AVX

#endif // SIMD_FLOAT_H
