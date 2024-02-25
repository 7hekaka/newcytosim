// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// Wednesday 24 June 2020 was a very nice day in Strasbourg

#ifndef SIMD_FLOAT_H
#define SIMD_FLOAT_H

// restrict functions to local scope
#define LOCAL static inline

//--------------------------- SSE Single Precision -----------------------------

#if defined(__SSE3__)

#include <immintrin.h>

// use the 4-vector since there is no penalty for using larger size:
typedef __m128 vec2f;

LOCAL vec2f setzero2f() { return _mm_setzero_ps(); }

LOCAL vec2f load1f(float const* a) { return _mm_load_ss(a); }
LOCAL vec2f load2f(float const* a) { return _mm_castsi128_ps(_mm_loadl_epi64((__m128i*)a)); }
LOCAL vec2f loaddupf(float const* a) { return _mm_load1_ps(a); }

LOCAL void store1f(float* a, vec2f b) { _mm_store_ss(a, b); }
LOCAL void store2f(float* a, vec2f b) { _mm_storel_pi((__m64*)a, b); }

LOCAL vec2f add2f(vec2f a, vec2f b) { return _mm_add_ps(a,b); }
LOCAL vec2f sub2f(vec2f a, vec2f b) { return _mm_sub_ps(a,b); }
LOCAL vec2f mul2f(vec2f a, vec2f b) { return _mm_mul_ps(a,b); }

/// return { a[0], b[1] }
//LOCAL vec2f blend11f(vec2f a, vec2f b) { return _mm_blend_ps(a,b,0b1010); }
LOCAL vec2f blend11f(vec2f a, vec2f b) { return _mm_shuffle_ps(_mm_shuffle_ps(a, b, 0x44), b, 0xEC); }

// returns { a[0], a[0] }
LOCAL vec2f duplo2f(vec2f a) { return _mm_moveldup_ps(a); }
// returns { a[1], a[1] }
LOCAL vec2f duphi2f(vec2f a) { return _mm_movehdup_ps(a); }

// returns { a[0], b[0] }
LOCAL vec2f unpacklo2f(vec2f a, vec2f b) { return _mm_movelh_ps(a, b); }
// returns { b[1], a[1] }
LOCAL vec2f unpackhi2f(vec2f a, vec2f b) { return _mm_movehl_ps(a, b); }

#if defined(__FMA__)
/// a * b + c
LOCAL vec2f fmadd2f(vec2f a, vec2f b, vec2f c)  { return _mm_fmadd_ps(c,a,b); }
/// c - a * b
LOCAL vec2f fnmadd2f(vec2f a, vec2f b, vec2f c) { return _mm_fnmadd_ps(c,a,b); }
#else
LOCAL vec2f fmadd2f(vec2f a, vec2f b, vec2f c)  { return _mm_add_ps(_mm_mul_ps(a,b), c); }
LOCAL vec2f fnmadd2f(vec2f a, vec2f b, vec2f c) { return _mm_sub_ps(c, _mm_mul_ps(a,b)); }
#endif

//----------------------------- 4-floats Vectors -------------------------------

/// Vector of 4 floats
typedef __m128 vec4f;
typedef __m128i vec4i;

LOCAL vec4i cast4fi(vec4f a) { return _mm_castps_si128(a); }
LOCAL vec4f cast4if(vec4i a) { return _mm_castsi128_ps(a); }

LOCAL vec4f setzero4f()    { return _mm_setzero_ps(); }
LOCAL vec4f set4f(float a) { return _mm_set1_ps(a); }
LOCAL vec4i set4i(int a)   { return _mm_set1_epi32(a); }
LOCAL vec4f set4fi(int a)  { return cast4if(_mm_set1_epi32(a)); }

/// _mm_load_ss loads a single element and zero the upper three:
LOCAL vec4f load4f(float const* a)     { return _mm_load_ps(a); }
LOCAL vec4f loadu4f(float const* a)    { return _mm_loadu_ps(a); }

LOCAL void store3f(float* a, vec4f b)  { a[0]=b[0]; a[1]=b[1]; a[2]=b[2]; }
LOCAL void store4f(float* a, vec4f b)  { _mm_store_ps(a,b); }
LOCAL void storeu4f(float* a, vec4f b) { _mm_storeu_ps(a,b); }

LOCAL void storeL2f(float* a, vec4f b) { _mm_storel_pi((__m64*)a, b); }
LOCAL void storeU2f(float* a, vec4f b) { _mm_storeh_pi((__m64*)a, b); }
/// convert 2 singles in lower positions, and store them in double precision
LOCAL void store2d(double* a, vec4f b) { _mm_store_pd(a, _mm_cvtps_pd(b)); }

LOCAL vec4f add4f(vec4f a, vec4f b)    { return _mm_add_ps(a,b); }
LOCAL vec4f sub4f(vec4f a, vec4f b)    { return _mm_sub_ps(a,b); }
LOCAL vec4f mul4f(vec4f a, vec4f b)    { return _mm_mul_ps(a,b); }
LOCAL vec4f div4f(vec4f a, vec4f b)    { return _mm_div_ps(a,b); }

LOCAL vec4f min4f(vec4f a, vec4f b)    { return _mm_min_ps(a,b); }
LOCAL vec4f max4f(vec4f a, vec4f b)    { return _mm_max_ps(a,b); }

LOCAL vec4f or4f(vec4f a, vec4f b)     { return _mm_or_ps(a,b); }
LOCAL vec4f and4f(vec4f a, vec4f b)    { return _mm_and_ps(a,b); }
LOCAL vec4f andnot4f(vec4f a, vec4f b) { return _mm_andnot_ps(a,b); }

LOCAL vec4f abs4f(vec4f a)             { return _mm_andnot_ps(_mm_set1_ps(-0.0f), a); }


// returns { a[0], b[0], a[1], b[1] }
LOCAL vec4f unpacklo4f(vec4f a, vec4f b) { return _mm_unpacklo_ps(a,b); }
// returns { a[2], b[2], a[3], b[3] }
LOCAL vec4f unpackhi4f(vec4f a, vec4f b) { return _mm_unpackhi_ps(a,b); }

// returns { a[0], a[0], a[2], a[2] }
LOCAL vec4f duplo4f(vec4f a) { return _mm_moveldup_ps(a); }
// returns { a[1], a[1], a[3], a[3] }
LOCAL vec4f duphi4f(vec4f a) { return _mm_movehdup_ps(a); }

// return { a[0], a[1], b[0], b[1] }
LOCAL vec4f movelh4f(vec4f a, vec4f b) { return _mm_movelh_ps(a, b); }
// return { b[2], b[3], a[2], a[2] }
LOCAL vec4f movehl4f(vec4f a, vec4f b) { return _mm_movehl_ps(a, b); }

// return { a[0], a[1] }
LOCAL vec2f getlo2f(vec4f a) { return (vec2f)a; }
// return { a[2], a[3] }
LOCAL vec2f gethi2f(vec4f a) { return (vec2f)_mm_movehl_ps(a, a); }

// return { A1, A2, A3, B0 } from a = { A0, A1, A2, A3 } and b = { B0, B1, B2, B3 }
LOCAL vec4f catshift1f(vec4f a, vec4f b) { return cast4if(_mm_alignr_epi8(cast4fi(b), cast4fi(a), 4)); }
// return { A2, A3, B0, B1 } from a = { A0, A1, A2, A3 } and b = { B0, B1, B2, B3 }
LOCAL vec4f catshift2f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, b, 0x4E); }
// return { A3, B0, B1, B2 } from a = { A0, A1, A2, A3 } and b = { B0, B1, B2, B3 }
LOCAL vec4f catshift3f(vec4f a, vec4f b) { return cast4if(_mm_alignr_epi8(cast4fi(b), cast4fi(a), 12)); }

LOCAL vec4f cmplt4f(vec4f a, vec4f b) { return _mm_cmplt_ps(a, b); }
LOCAL vec4f cmpgt4f(vec4f a, vec4f b) { return _mm_cmpgt_ps(a, b); }
LOCAL vec4f cmple4f(vec4f a, vec4f b) { return _mm_cmple_ps(a, b); }
LOCAL vec4f cmpge4f(vec4f a, vec4f b) { return _mm_cmpge_ps(a, b); }


/// set i-th bit in returned value to ( a[i] < b[i] ), for i = {0, 1, 2, 3}
LOCAL int lower_mask4f(vec4f a, vec4f b) { return _mm_movemask_ps(_mm_cmplt_ps(a,b)); }
/// true if any float component is non-zero
LOCAL int any_true4f(vec4f a) { return !_mm_test_all_zeros((__m128i)a, (__m128i{-1l, -1l})); }

LOCAL vec4i shiftbitsR4(vec4f a, int b) { return _mm_srli_epi32(cast4fi(a), b); }
LOCAL vec4i shiftbitsL4(vec4f a, int b) { return _mm_slli_epi32(cast4fi(a), b); }

/// extract component
//LOCAL int32_t getlane4i(vec4i a, int i) { return _mm_extract_epi32(a, i); }
#define getlane4i(a, i) _mm_extract_epi32(a, i);


/// load 4 32-bit integers
LOCAL vec4i load4i(int32_t const* a) { return _mm_load_si128((__m128i*)a); }
/// convert integers to floats:
LOCAL vec4f cvt4if(__m128i a) { return _mm_cvtepi32_ps(a); }
/// load integers and convert to float:
LOCAL vec4f load4if(int32_t const* a) { return _mm_cvtepi32_ps(_mm_load_si128((__m128i const*)a)); }

/// square root
LOCAL vec4f sqrt4f(vec4f a) { return _mm_sqrt_ps(a); }
/// approximate reciprocal square root: 1 / sqrt(a)
LOCAL vec4f rsqrt4f(vec4f a) { return _mm_rsqrt_ps(a); }

LOCAL vec4f negative4f(vec4f a) { return _mm_cmplt_ps(a, setzero4f()); }
LOCAL vec4f positive4f(vec4f a) { return _mm_cmpgt_ps(a, setzero4f()); }
LOCAL vec4f notpositive4f(vec4f a) { return _mm_cmpngt_ps(a, setzero4f()); }
LOCAL vec4f greaterequal4f(vec4f a, vec4f b) { return _mm_cmpge_ps(a, b); }
LOCAL vec4f lowerthan4f(vec4f a, vec4f b) { return _mm_cmplt_ps(a, b); }

#endif  // __SSE3__

#if defined(__SSE4_1__)

/// return { a[0], a[1], a[2], b[3] }
LOCAL vec4f blend31f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1000); }
/// return { a[0], a[1], b[2], b[3] }
LOCAL vec4f blend22f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1100); }
/// return { a[0], b[1], b[2], b[3] }
LOCAL vec4f blend13f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b1110); }

LOCAL vec4f clear4th(vec4f a) { return _mm_blend_ps(a,_mm_setzero_ps(),0b1000); }

/// return { a[0], a[1], b[2], a[3] }
LOCAL vec4f blend0010f(vec4f a, vec4f b) { return _mm_blend_ps(a,b,0b0100); }

/// blend `a` and `b` based on mask `k`: select 'b' if 'topmost bit of k == 1' and 'a' otherwise
LOCAL vec4f blendv4f(vec4f a, vec4f b, vec4f k) { return _mm_blendv_ps(a,b,k); }

/// return `neg` if `val < 0` and `pos` otherwise
LOCAL vec4f signselect4f(vec4f val, vec4f neg, vec4f pos) { return _mm_blendv_ps(pos, neg, val); }

LOCAL vec4f load3f(float const* a) { return _mm_blend_ps(_mm_loadu_ps(a), _mm_setzero_ps(), 0b1000); }
// loading 4 and clearing one
LOCAL vec4f load3fZ(float const* a) { return _mm_blend_ps(_mm_loadu_ps(a), _mm_setzero_ps(), 0b1000); }

#elif defined(__SSE3__)

// emulating the blend function using two shuffles
LOCAL vec4f blend31f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, _mm_shuffle_ps(a,b,0xEE), 0xC4); }
LOCAL vec4f blend22f(vec4f a, vec4f b) { return _mm_shuffle_ps(a, b, 0xE4); }
LOCAL vec4f blend13f(vec4f a, vec4f b) { return _mm_shuffle_ps(_mm_shuffle_ps(a, b, 0x44), b, 0xEC); }

LOCAL vec4f clear4th(vec4f a) { return _mm_shuffle_ps(a, _mm_shuffle_ps(a,_mm_setzero_ps(),0xEE), 0xC4); }

// loading 3 elements
LOCAL vec4f load3f(float const* a) { return vec4f{a[0], a[1], a[2], 0}; }
// loading 4 and clearing one
LOCAL vec4f load3fZ(float const* a) { return clear4th(loadu4f(a)); }

#endif

#if defined(__AVX__)
// comparison
#define cmp4f(a,b,c) _mm_cmp_ps(a,b,c)

LOCAL vec4f broadcast1f(float const* a) { return _mm_broadcast_ss(a); }

// copy a[0] into all elements of destination
LOCAL vec4f broadcastXf(vec4f a) { return _mm_permute_ps(a,0x00); }
// copy a[1] into all elements of destination
LOCAL vec4f broadcastYf(vec4f a) { return _mm_permute_ps(a,0x55); }
// copy a[2] into all elements of destination
LOCAL vec4f broadcastZf(vec4f a) { return _mm_permute_ps(a,0xAA); }
// copy a[3] into all elements of destination
LOCAL vec4f broadcastTf(vec4f a) { return _mm_permute_ps(a,0xFF); }

// non-temporal load
LOCAL vec4f streamload4f(float const* a) { return cast4if(_mm_stream_load_si128((__m128i*)a)); }

#define permute4f(a,k) _mm_permute_ps(a,k)
// Convert between single and double types
LOCAL vec4f cvt4ds(__m256d a) { return _mm256_cvtpd_ps(a); }
LOCAL __m256d cvt4sd(vec4f a) { return _mm256_cvtps_pd(a); }

/// convert double and store them in single precision
LOCAL void store4df(float* a, __m256d b) { _mm_storeu_ps(a, _mm256_cvtpd_ps(b)); }

#elif defined(__SSE3__)

LOCAL vec4f broadcast1f(float const* a) { return _mm_load1_ps(a); }
LOCAL vec4f streamload4f(float const* a) { return _mm_load_ps(a); }

LOCAL vec4f broadcastXf(vec4f a) { return _mm_shuffle_ps(a,a,0x00); }
LOCAL vec4f broadcastYf(vec4f a) { return _mm_shuffle_ps(a,a,0x55); }
LOCAL vec4f broadcastZf(vec4f a) { return _mm_shuffle_ps(a,a,0xAA); }
LOCAL vec4f broadcastTf(vec4f a) { return _mm_shuffle_ps(a,a,0xFF); }

#define permute4f(a,k) _mm_shuffle_ps(a,a,k)

#endif

//-------------------------- FMA Single Precision-------------------------------

#if defined(__FMA__)
LOCAL vec4f fmadd4f (vec4f a, vec4f b, vec4f c) { return _mm_fmadd_ps(a,b,c); }
LOCAL vec4f fnmadd4f(vec4f a, vec4f b, vec4f c) { return _mm_fnmadd_ps(a,b,c); }
#elif defined(__SSE3__)
// erzatz functions
LOCAL vec4f fmadd4f (vec4f a, vec4f b, vec4f c) { return _mm_add_ps(_mm_mul_ps(a,b), c); }  // a * b + c
LOCAL vec4f fnmadd4f(vec4f a, vec4f b, vec4f c) { return _mm_sub_ps(c, _mm_mul_ps(a,b)); }
#endif

//-------------------------- AVX Single Precision-------------------------------

#if defined(__AVX__)

/// Vector of 8 floats
typedef __m256 vec8f;
typedef __m256i vec8i;

LOCAL vec8f setzero8f()            { return _mm256_setzero_ps(); }
LOCAL vec8f set8f(float a)         { return _mm256_set1_ps(a); }
LOCAL vec8i set8i(int a)           { return _mm256_set1_epi32(a); }
LOCAL vec8f set8fi(int a)          { return _mm256_castsi256_ps(_mm256_set1_epi32(a)); }
LOCAL vec8f load8f(float const* a) { return _mm256_load_ps(a); }
LOCAL vec8f loadu8f(float const* a) { return _mm256_loadu_ps(a); }

LOCAL void storelof(float* a, vec8f b)   { _mm_store_ss(a,_mm256_castps256_ps128(b)); }
LOCAL void store8f(float* a, vec8f b)    { _mm256_store_ps(a, b); }
LOCAL void storeu8f(float* a, vec8f b)   { _mm256_storeu_ps(a, b); }

LOCAL vec8f mul8f(vec8f a, vec8f b)      { return _mm256_mul_ps(a,b); }
LOCAL vec8f div8f(vec8f a, vec8f b)      { return _mm256_div_ps(a,b); }
LOCAL vec8f add8f(vec8f a, vec8f b)      { return _mm256_add_ps(a,b); }
LOCAL vec8f sub8f(vec8f a, vec8f b)      { return _mm256_sub_ps(a,b); }

LOCAL vec8f max8f(vec8f a, vec8f b)      { return _mm256_max_ps(a,b); }
LOCAL vec8f min8f(vec8f a, vec8f b)      { return _mm256_min_ps(a,b); }

LOCAL vec8f or8f(vec8f a, vec8f b)       { return _mm256_or_ps(a,b); }
LOCAL vec8f and8f(vec8f a, vec8f b)      { return _mm256_and_ps(a,b); }
LOCAL vec8f andnot8f(vec8f a, vec8f b)   { return _mm256_andnot_ps(a,b); }

LOCAL vec8f abs8f(vec8f a)               { return _mm256_andnot_ps(_mm256_set1_ps(-0.0), a); }
LOCAL vec8f flipsign8f(vec8f a)          { return _mm256_xor_ps(a, _mm256_set1_ps(-0.0)); }
LOCAL vec8f lowerthan8f(vec8f a, vec8f b) { return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }

LOCAL vec8f unpacklo8f(vec8f a, vec8f b) { return _mm256_unpacklo_ps(a,b); }
LOCAL vec8f unpackhi8f(vec8f a, vec8f b) { return _mm256_unpackhi_ps(a,b); }

LOCAL vec8f isnan8f(vec8f a) { return _mm256_cmp_ps(a,a,3); }
// blend to select `b` if `k == true`, and `a` otherwise:
LOCAL vec8f blendv8f(vec8f a, vec8f b, vec8f k) { return _mm256_blendv_ps(a,b,k); }
LOCAL vec8f swap4f128(vec8f a) { return _mm256_permute2f128_ps(a, a, 0x01); }
// return { 1 0 3 2  5 4 7 6 }, permutting positions 1&2, 3&4, 5&6, 7&8
LOCAL vec8f permute8f(vec8f a) { return _mm256_permute_ps(a, 0xB1); }
// return { 2 3 0 1  6 7 4 5 }, permutting positions 1+2 with 3+4, 4+5 with 6+7
LOCAL vec8f permute44f(vec8f a) { return _mm256_permute_ps(a, 0x4E); }


#define permute8f128(a,b,c)  _mm256_permute4f128_ps(a,b,c)
#define getlane8i(a, b) _mm256_extract_epi32(a, b);


/// approximate inverse: 1/a
LOCAL vec8f rcp8f(vec8f a)    { return _mm256_rcp_ps(a); }
/// square root
LOCAL vec8f sqrt8f(vec8f a)   { return _mm256_sqrt_ps(a); }
/// approximate reciprocal square root: 1 / sqrt(a)
LOCAL vec8f rsqrt8f(vec8f a)  { return _mm256_rsqrt_ps(a); }

LOCAL vec4f getlo4f(vec8f a)  { return _mm256_castps256_ps128(a); }
LOCAL vec4f gethi4f(vec8f a)  { return _mm256_extractf128_ps(a,1); }
/// concatenate two vec4f into a vec8f
LOCAL vec8f cat44f(vec4f l, vec4f h) { return _mm256_insertf128_ps(_mm256_castps128_ps256(l), h, 1); }

LOCAL vec8f cvt8if(vec8i a) { return _mm256_cvtepi32_ps(a); }
LOCAL vec8f cast8f(vec8i a) { return _mm256_castsi256_ps(a); }

LOCAL vec8i load8i(vec8i const* a) { return _mm256_load_si256(a); }
/// load integers and convert to float:
LOCAL vec8f load8if(vec8i const* a) { return _mm256_cvtepi32_ps(_mm256_load_si256(a)); }

#define cmp8f(a,b,c)         _mm256_cmp_ps(a,b,c)
#define permute2f128f(a,b,c) _mm256_permute2f128_ps(a,b,c)

/// return `neg` if `val < 0` and `pos` otherwise
LOCAL vec8f signselect8f(vec8f val, vec8f neg, vec8f pos) { return _mm256_blendv_ps(pos, neg, val); }

#endif // AVX

#if defined(__AVX2__)

LOCAL vec4i and4i(vec4i a, vec4i b) { return (vec4i)_mm_and_ps((vec4f)a,(vec4f)b); }
LOCAL vec8i shiftbitsR8(vec8f a, int b) { return _mm256_srli_epi32(_mm256_castps_si256(a), b); }
LOCAL vec8i shiftbitsL8(vec8f a, int b) { return _mm256_slli_epi32(_mm256_castps_si256(a), b); }

#endif

#if 0

/// convert 1 half-float to float
LOCAL unsigned short cvt1sh(float a) { return _cvtss_sh(a, _MM_FROUND_NO_EXC); }
/// convert 1 half-floats to floats
LOCAL float cvt4hs(unsigned short a) { return _cvtsh_ss(a); }
/// convert 4 floats to half-floats
LOCAL __m128i cvt4sh(__m128 a) { return _mm_cvtps_ph(a, _MM_FROUND_NO_EXC); }
/// convert 4 half-floats to floats
LOCAL vec4f cvt4hs(__m128i a) { return _mm_cvtph_ps(a); }

#endif

#undef LOCAL

#endif // SIMD_FLOAT_H
