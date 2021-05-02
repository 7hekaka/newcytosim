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
inline static vec4f set4fi(uint32_t a)            { return _mm_castsi128_ps(_mm_set1_epi32(a)); }
inline static vec4f load1f(float const* a)        { return _mm_load_ss(a); }
inline static vec4f load2f(float const* a)        { return (vec4f)_mm_loadl_epi64((__m128i*)a); }
inline static vec4f load4f(float const* a)        { return _mm_load_ps(a); }
inline static vec4f loadu4f(float const* a)       { return _mm_loadu_ps(a); }

inline static void store1f(float* a, vec4f b)     { _mm_store_ss(a,b); }
//inline static void store2f(float* a, vec4f b)     { _mm_storeu_si64((void*)a,_mm_castps_si128(b)); } //_mm_storel_pi((__m64*)a, b); }
inline static void store3f(float* a, vec4f b)     { a[0]=b[0]; a[1]=b[1]; a[2]=b[2]; }
inline static void store4f(float* a, vec4f b)     { _mm_store_ps(a,b); }
inline static void storeu4f(float* a, vec4f b)    { _mm_storeu_ps(a,b); }

inline static vec4f add4f(vec4f a, vec4f b)       { return _mm_add_ps(a,b); }
inline static vec4f sub4f(vec4f a, vec4f b)       { return _mm_sub_ps(a,b); }
inline static vec4f mul4f(vec4f a, vec4f b)       { return _mm_mul_ps(a,b); }

inline static vec4f max4f(vec4f a, vec4f b)       { return _mm_max_ps(a,b); }
inline static vec4f min4f(vec4f a, vec4f b)       { return _mm_min_ps(a,b); }

inline static vec4f or4f(vec4f a, vec4f b)        { return _mm_or_ps(a,b); }
inline static vec4f and4f(vec4f a, vec4f b)       { return _mm_and_ps(a,b); }
inline static vec4f andnot4f(vec4f a, vec4f b)    { return _mm_andnot_ps(a,b); }
inline static vec4f abs4f(vec4f a)                { return _mm_andnot_ps(_mm_set1_ps(-0.0f), a); }

inline static vec4f unpacklo4f(vec4f a, vec4f b)  { return _mm_unpacklo_ps(a,b); }
inline static vec4f unpackhi4f(vec4f a, vec4f b)  { return _mm_unpackhi_ps(a,b); }

inline static vec4f duplo4f(vec4f a)              { return _mm_unpacklo_ps(a,a); }
inline static vec4f duphi4f(vec4f a)              { return _mm_unpackhi_ps(a,a); }

// return { A1, B0 } from a = { A0, A1 } and b = { B0, B1 }
inline static vec4f twine2f64(vec4f a, vec4f b) { return _mm_shuffle_ps(a, b, 0x4E); }

inline static vec4f cmplt4f(vec4f a, vec4f b) { return _mm_cmplt_ps(a, b); }
inline static vec4f cmpgt4f(vec4f a, vec4f b) { return _mm_cmpgt_ps(a, b); }
inline static vec4f cmple4f(vec4f a, vec4f b) { return _mm_cmple_ps(a, b); }
inline static vec4f cmpge4f(vec4f a, vec4f b) { return _mm_cmpge_ps(a, b); }


/// convert integer to float:
inline static vec4f cvt4if(__m128i a) { return _mm_cvtepi32_ps(a); }
inline static vec4f cast4f(__m128i a) { return _mm_castsi128_ps(a); }

/// approximate reciprocal square root: 1 / sqrt(a)
inline static vec4f rsqrt4f(vec4f a)  { return _mm_rsqrt_ps(a); }

#define shuffle4f(a,b,k) _mm_shuffle_ps(a,b,k)

#endif  // __SSE3__


#if defined(__SSE4_1__)

inline static vec4f blend31f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1000); }
inline static vec4f blend22f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1100); }
inline static vec4f blend13f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1110); }

#  define blend4f(a,b,k) _mm_blend_ps(a,b,k)

inline static vec4f sign_select4f(vec4f val, vec4f neg, vec4f pos)
{
#if defined(__AVX512VL__)
    return _mm_mask_mov_ps(pos, val, neg);
#else
    return _mm_blendv_ps(pos, neg, val);
#endif
}

// loading 4 and clearing one
inline static vec4f load3f(float const* a) { return _mm_blend_ps(_mm_loadu_ps(a), _mm_setzero_ps(), 0b1000); }

#elif defined(__SSE3__)

inline static vec4f blend31f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, _mm_shuffle_ps(a,b,0xEE), 0xC4); }

// loading 3 elements
inline static vec4f load3f(float const* a) { vec4f b; b[0]=a[0]; b[1]=a[2]; b[1]=a[2]; return b; }

#endif


#if defined(__AVX__)
// comparison
#define cmp4f(a,b,c)     _mm_cmp_ps(a,b,c)
// copy a[0] into all elements of destination
inline static vec4f broadcastlof(vec4f a)         { return _mm_permute_ps(a,0x00); }
inline static vec4f broadcast1f(float const* a)   { return _mm_broadcast_ss(a); }
inline static vec4f streamload4f(float const* a)  { return (vec4f)_mm_stream_load_si128((__m128i*)a); }
#define permute4f(a,k)    _mm_permute_ps(a,k)
// Convert between single and double types
inline static vec4f cvt4ds(__m256d a)             { return _mm256_cvtpd_ps(a); }
inline static __m256d cvt4sd(vec4f a)             { return _mm256_cvtps_pd(a); }
inline static void store4d(float* a, __m256d b)   { _mm_storeu_ps(a, _mm256_cvtpd_ps(b)); }
#elif defined(__SSE3__)
inline static vec4f broadcastlof(vec4f a)         { return _mm_shuffle_ps(a,a,0x00); }
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

inline static vec8f setzero8f()                  { return _mm256_setzero_ps(); }
inline static vec8f set8f(float a)               { return _mm256_set1_ps(a); }
inline static vec8f set8fi(uint32_t a)           { return _mm256_castsi256_ps(_mm256_set1_epi32(a)); }
inline static vec8f load8f(float const* a)       { return _mm256_load_ps(a); }
inline static vec8f loadu8f(float const* a)      { return _mm256_loadu_ps(a); }

inline static void storelof(float* a, vec8f b)   { _mm_store_ss(a,_mm256_castps256_ps128(b)); }
inline static void store8f(float* a, vec8f b)    { _mm256_store_ps(a, b); }
inline static void storeu8f(float* a, vec8f b)   { _mm256_storeu_ps(a, b); }

inline static vec8f mul8f(vec8f a, vec8f b)      { return _mm256_mul_ps(a,b); }
inline static vec8f div8f(vec8f a, vec8f b)      { return _mm256_div_ps(a,b); }
inline static vec8f add8f(vec8f a, vec8f b)      { return _mm256_add_ps(a,b); }
inline static vec8f sub8f(vec8f a, vec8f b)      { return _mm256_sub_ps(a,b); }

inline static vec8f max8f(vec8f a, vec8f b)      { return _mm256_max_ps(a,b); }
inline static vec8f min8f(vec8f a, vec8f b)      { return _mm256_min_ps(a,b); }

inline static vec8f or8f(vec8f a, vec8f b)       { return _mm256_or_ps(a,b); }
inline static vec8f and8f(vec8f a, vec8f b)      { return _mm256_and_ps(a,b); }
inline static vec8f andnot8f(vec8f a, vec8f b)   { return _mm256_andnot_ps(a,b); }
inline static vec8f abs8f(vec8f a)               { return _mm256_andnot_ps(_mm256_set1_ps(-0.0), a); }
inline static vec8f flipsign8f(vec8f a)          { return _mm256_xor_ps(a, _mm256_set1_ps(-0.0)); }

inline static vec8f unpacklo8f(vec8f a, vec8f b) { return _mm256_unpacklo_ps(a,b); }
inline static vec8f unpackhi8f(vec8f a, vec8f b) { return _mm256_unpackhi_ps(a,b); }

inline static vec8f isnan8f(vec8f a) { return _mm256_cmp_ps(a,a,3); }
// blend to select `b` if `k == true`, and `a` otherwise:
inline static vec8f blendv8f(vec8f a, vec8f b, vec8f k) { return _mm256_blendv_ps(a,b,k); }
inline static vec8f swap2f128(vec8f a) { return _mm256_permute2f128_ps(a, a, 0x01); }
// permute positions 1&2, 3&4, 5&6, 7&8
inline static vec8f permute8f(vec8f a) { return _mm256_permute_ps(a, 0xB1); }
// permute positions 12 and 34, 45 and 67
inline static vec8f permute44f(vec8f a) { return _mm256_permute_ps(a, 0x4E); }


#define permute8f128(a,b,c)  _mm256_permute4f128_ps(a,b,c)


/// approximate inverse: 1/a
inline static vec8f rcp8f(vec8f a)    { return _mm256_rcp_ps(a); }
/// square root
inline static vec8f sqrt8f(vec8f a)   { return _mm256_sqrt_ps(a); }
/// approximate reciprocal square root: 1 / sqrt(a)
inline static vec8f rsqrt8f(vec8f a)  { return _mm256_rsqrt_ps(a); }

inline static vec4f getlo4f(vec8f a)  { return _mm256_castps256_ps128(a); }
inline static vec4f gethi4f(vec8f a)  { return _mm256_extractf128_ps(a,1); }
/// concatenate two vec4f into a vec8f
inline static vec8f cat4f(vec4f h, vec4f l) { return _mm256_set_m128(h, l); }

inline static vec8f cvt8if(__m256i a) { return _mm256_cvtepi32_ps(a); }
inline static vec8f cast8f(__m256i a) { return _mm256_castsi256_ps(a); }

#define load8si(a)           _mm256_load_si256(a)
#define cmp8f(a,b,c)         _mm256_cmp_ps(a,b,c)
#define permute2f128f(a,b,c) _mm256_permute2f128_ps(a,b,c)

inline static vec8f sign_select8f(vec8f val, vec8f neg, vec8f pos)
{
#if defined(__AVX512VL__)
    return _mm256_mask_mov_ps(pos, val, neg);
#else
    return _mm256_blendv_ps(pos, neg, val);
#endif
}

#endif // AVX

#endif // SIMD_FLOAT_H
