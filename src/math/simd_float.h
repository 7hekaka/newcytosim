// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// Wednesday 24 June 2020 was a very nice day in Strasbourg

#ifndef SIMD_FLOAT_H
#define SIMD_FLOAT_H

//--------------------------- SSE Single Precision -----------------------------

#if defined(__SSE3__)

#include <immintrin.h>

/// Vector of 4 floats
typedef __m128 vec4f;

inline static vec4f setzero4f()                   { return _mm_setzero_ps(); }
inline static vec4f set4f(float a)                { return _mm_set1_ps(a); }
inline static vec4f set4fi(int a)                 { return _mm_castsi128_ps(_mm_set1_epi32(a)); }
inline static vec4f load1f(float const* a)        { return _mm_load_ss(a); }
inline static vec4f load2f(float const* a)        { return _mm_castsi128_ps(_mm_loadl_epi64((__m128i*)a)); }
inline static vec4f load4f(float const* a)        { return _mm_load_ps(a); }
inline static vec4f loadu4f(float const* a)       { return _mm_loadu_ps(a); }

inline static void store1f(float* a, vec4f b)     { _mm_store_ss(a,b); }
inline static void store2f(float* a, vec4f b)     { _mm_storel_pi((__m64*)a, b); }
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

// returns { a[0], b[0], a[2], b[2] }
inline static vec4f unpacklo4f(vec4f a, vec4f b)  { return _mm_unpacklo_ps(a,b); }
// returns { a[1], b[1], a[3], b[3] }
inline static vec4f unpackhi4f(vec4f a, vec4f b)  { return _mm_unpackhi_ps(a,b); }

// returns { a[0], a[0], a[2], a[2] }
inline static vec4f duplo4f(vec4f a)             { return _mm_shuffle_ps(a, a, 0xA0); }
// returns { a[1], a[1], a[3], a[3] }
inline static vec4f duphi4f(vec4f a)             { return _mm_movehdup_ps(a); } //_mm_shuffle_ps(a, a, 0xF5); }

// return { B1, A1 } from a = { A0, A1 } and b = { B0, B1 }
inline static vec4f movehl4f(vec4f a, vec4f b) { return _mm_movehl_ps(a, b); }
// return { A0, B0 } from a = { A0, A1 } and b = { B0, B1 }
inline static vec4f movelh4f(vec4f a, vec4f b) { return _mm_movelh_ps(a, b); }

// return { A0, A1, A0, A1 }
inline static vec4f getlof(vec4f a)            { return _mm_movelh_ps(a, a); }
// return { A2, A3, A2, A3 }
inline static vec4f gethif(vec4f a)            { return _mm_movehl_ps(a, a); }

// return { A1, A2, A3, B0 } from a = { A0, A1, A2, A3 } and b = { B0, B1, B2, B3 }
inline static vec4f catshift1f(vec4f a, vec4f b) { return _mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(b), _mm_castps_si128(a), 4)); }
// return { A2, A3, B0, B1 } from a = { A0, A1, A2, A3 } and b = { B0, B1, B2, B3 }
inline static vec4f catshift2f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, b, 0x4E); }
// return { A3, B0, B1, B2 } from a = { A0, A1, A2, A3 } and b = { B0, B1, B2, B3 }
inline static vec4f catshift3f(vec4f a, vec4f b) { return _mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(b), _mm_castps_si128(a), 12)); }

inline static vec4f cmplt4f(vec4f a, vec4f b) { return _mm_cmplt_ps(a, b); }
inline static vec4f cmpgt4f(vec4f a, vec4f b) { return _mm_cmpgt_ps(a, b); }
inline static vec4f cmple4f(vec4f a, vec4f b) { return _mm_cmple_ps(a, b); }
inline static vec4f cmpge4f(vec4f a, vec4f b) { return _mm_cmpge_ps(a, b); }


// set i-th bit in returned value if a[i] < b[i], for i = {0, 1, 2, 3}
inline static int lower_mask4f(vec4f a, vec4f b) { return _mm_movemask_ps(_mm_cmplt_ps(a,b)); }
inline static int any_true4f(vec4f a) { return !_mm_test_all_zeros((__m128i)a, (__m128i{-1l, -1l})); }

/// convert integer to float:
inline static vec4f cvt4if(__m128i a) { return _mm_cvtepi32_ps(a); }
inline static vec4f cast4f(__m128i a) { return _mm_castsi128_ps(a); }

/// approximate reciprocal square root: 1 / sqrt(a)
inline static vec4f rsqrt4f(vec4f a)  { return _mm_rsqrt_ps(a); }

#endif  // __SSE3__

#if defined(__SSE4_1__)

inline static vec4f blend31f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1000); }
inline static vec4f blend22f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1100); }
inline static vec4f blend13f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1110); }

inline static vec4f clear4th(vec4f a) { return _mm_blend_ps(a,_mm_setzero_ps(),0b1000); }

/// return { a[0], a[1], b[2], a[3] }
inline static vec4f blend0010f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b0100); }

/// return `neg` if `val < 0` and `pos` otherwise
inline static vec4f sign_select4f(vec4f val, vec4f neg, vec4f pos) { return _mm_blendv_ps(pos, neg, val); }

inline static vec4f load3f(float const* a) { return _mm_blend_ps(_mm_loadu_ps(a), _mm_setzero_ps(), 0b1000); }
// loading 4 and clearing one
inline static vec4f load3fZ(float const* a) { return _mm_blend_ps(_mm_loadu_ps(a), _mm_setzero_ps(), 0b1000); }

#elif defined(__SSE3__)

// emulating the blend function using two shuffles
inline static vec4f blend31f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, _mm_shuffle_ps(a,b,0xEE), 0xC4); }
inline static vec4f blend22f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, b, 0xE4); }
inline static vec4f blend13f(vec4f a, vec4f b) { return _mm_shuffle_ps(_mm_shuffle_ps(a, b, 0x44), b, 0xEC); }

inline static vec4f clear4th(vec4f a) { return _mm_shuffle_ps(a, _mm_shuffle_ps(a,_mm_setzero_ps(),0xEE), 0xC4); }

// loading 3 elements
inline static vec4f load3f(float const* a) { return vec4f{a[0], a[1], a[2], 0}; }
// loading 4 and clearing one
inline static vec4f load3fZ(float const* a) { return clear4th(loadu4f(a)); }

#endif


#if defined(__AVX__)
// comparison
#define cmp4f(a,b,c)     _mm_cmp_ps(a,b,c)
inline static vec4f broadcast1f(float const* a)   { return _mm_broadcast_ss(a); }

// copy a[0] into all elements of destination
inline static vec4f broadcastXf(vec4f a)          { return _mm_permute_ps(a,0x00); }
inline static vec4f broadcastYf(vec4f a)          { return _mm_permute_ps(a,0x55); }
inline static vec4f broadcastZf(vec4f a)          { return _mm_permute_ps(a,0xAA); }
inline static vec4f broadcastTf(vec4f a)          { return _mm_permute_ps(a,0xFF); }

// non-temporal load
inline static vec4f streamload4f(float const* a)  { return _mm_castsi128_ps(_mm_stream_load_si128((__m128i*)a)); }

#define permute4f(a,k)    _mm_permute_ps(a,k)
// Convert between single and double types
inline static vec4f cvt4ds(__m256d a)             { return _mm256_cvtpd_ps(a); }
inline static __m256d cvt4sd(vec4f a)             { return _mm256_cvtps_pd(a); }

/// convert double and store them in single precision
inline static void store4df(float* a, __m256d b)  { _mm_storeu_ps(a, _mm256_cvtpd_ps(b)); }

#elif defined(__SSE3__)

inline static vec4f broadcast1f(float const* a)   { return _mm_load1_ps(a); }
inline static vec4f streamload4f(float const* a)  { return _mm_load_ps(a); }

inline static vec4f broadcastXf(vec4f a)          { return _mm_shuffle_ps(a,a,0x00); }
inline static vec4f broadcastYf(vec4f a)          { return _mm_shuffle_ps(a,a,0x55); }
inline static vec4f broadcastZf(vec4f a)          { return _mm_shuffle_ps(a,a,0xAA); }
inline static vec4f broadcastTf(vec4f a)          { return _mm_shuffle_ps(a,a,0xFF); }

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
inline static vec8f set8fi(int a)                { return _mm256_castsi256_ps(_mm256_set1_epi32(a)); }
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
inline static vec8f swap4f128(vec8f a) { return _mm256_permute2f128_ps(a, a, 0x01); }
// return { 1 0 3 2  5 4 7 6 }, permutting positions 1&2, 3&4, 5&6, 7&8
inline static vec8f permute8f(vec8f a) { return _mm256_permute_ps(a, 0xB1); }
// return { 2 3 0 1  6 7 4 5 }, permutting positions 1+2 with 3+4, 4+5 with 6+7
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
inline static vec8f cat44f(vec4f l, vec4f h) { return _mm256_insertf128_ps(_mm256_castps128_ps256(l), h, 1); }

inline static vec8f cvt8if(__m256i a) { return _mm256_cvtepi32_ps(a); }
inline static vec8f cast8f(__m256i a) { return _mm256_castsi256_ps(a); }

#define load8si(a)           _mm256_load_si256(a)
#define cmp8f(a,b,c)         _mm256_cmp_ps(a,b,c)
#define permute2f128f(a,b,c) _mm256_permute2f128_ps(a,b,c)

/// return `neg` if `val < 0` and `pos` otherwise
inline static vec8f sign_select8f(vec8f val, vec8f neg, vec8f pos) { return _mm256_blendv_ps(pos, neg, val); }

#endif // AVX

#if 0

/// convert 1 half-float to float
inline static unsigned short cvt1sh(float a) { return _cvtss_sh(a, _MM_FROUND_NO_EXC); }
/// convert 1 half-floats to floats
inline static float cvt4hs(unsigned short a) { return _cvtsh_ss(a); }
/// convert 4 floats to half-floats
inline static __m128i cvt4sh(__m128 a) { return _mm_cvtps_ph(a, _MM_FROUND_NO_EXC); }
/// convert 4 half-floats to floats
inline static vec4f cvt4hs(__m128i a) { return _mm_cvtph_ps(a); }

#endif

#endif // SIMD_FLOAT_H
