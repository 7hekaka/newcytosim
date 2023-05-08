// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Started on Monday 5 June 2018, which was a very nice day in Strasbourg


/// Vector holding 4 double precision floats
typedef __m256d vec4;

#define set64x(a,b,c,d) _mm256_setr_epi64x(a,b,c,d)

constexpr vec4 sgn1111 = {-0.0, -0.0, -0.0, -0.0};


LOCAL vec4 setr4(double a, double b, double c, double d) { return _mm256_setr_pd(a,b,c,d); }
LOCAL vec4 set4(double a, double b, double c, double d)  { return _mm256_set_pd(a,b,c,d); }

LOCAL vec4 set4(double a)               { return _mm256_set1_pd(a); }
LOCAL vec4 setzero4()                   { return _mm256_setzero_pd(); }

LOCAL vec4 cast4(vec2 a)                { return _mm256_castpd128_pd256(a); }
LOCAL vec2 cast2(vec4 a)                { return _mm256_castpd256_pd128(a); }

LOCAL vec4 load4(double const* a)       { return _mm256_load_pd(a); }
LOCAL vec4 loadu4(double const* a)      { return _mm256_loadu_pd(a); }

/// unaligned load 2 values, and zeros out the upper two
LOCAL vec4 load2Z(double const* a)      { return _mm256_blend_pd(_mm256_castpd128_pd256(_mm_loadu_pd(a)), _mm256_setzero_pd(), 0b1100); }

/// unaligned load 3 values, and zeros out the upper one
LOCAL vec4 load3Z(double const* a)      { return _mm256_blend_pd(_mm256_loadu_pd(a), _mm256_setzero_pd(), 0b1000); }

/// unaligned load 2 values, allowing for undefined value in upper positions
LOCAL vec4 load2crap(double const* a)   { return _mm256_castpd128_pd256(_mm_load_pd(a)); }

//LOCAL void store1(double* a, vec4 b)    { _mm_store_sd(a, cast2(b)); }
//LOCAL void store2(double* a, vec4 b)    { _mm_store_pd(a, cast2(b)); }
//LOCAL void storeu2(double* a, vec4 b)   { _mm_storeu_pd(a, cast2(b)); }
LOCAL void store4(double* a, vec4 b)    { _mm256_store_pd(a,b); }
LOCAL void storeu4(double* a, vec4 b)   { _mm256_storeu_pd(a,b); }

/// convert 4 singles and store them in double precision
LOCAL void store4d(double* a, __m128 b)  { _mm256_store_pd(a, _mm256_cvtps_pd(b)); }

LOCAL __m256i makemask(long i)
{
    constexpr __m256d ramp{0.5, 1.5, 2.5, 3.5};
    return _mm256_castpd_si256(_mm256_cmp_pd(ramp, _mm256_set1_pd((double)i), _CMP_LT_OQ));
}

LOCAL vec4 maskload4(double const* a, __m256i k)     { return _mm256_maskload_pd(a,k); }
LOCAL void maskstore4(double* a, __m256i k, vec4 b)  { _mm256_maskstore_pd(a,k,b); }

/// load 1 double into all 4 positions
LOCAL vec4 broadcast1(double const* a)  { return _mm256_broadcast_sd(a); }
/// load 2 doubles and duplicate: X, Y, X, Y
LOCAL vec4 broadcast2(double const* a)  { return _mm256_broadcast_pd((__m128d const*)a); }

/// return { A0, A1 }
LOCAL vec2 getlo(vec4 a)                { return _mm256_castpd256_pd128(a); }
/// return { A2, A3 }
LOCAL vec2 gethi(vec4 a)                { return _mm256_extractf128_pd(a,1); }
/// concatenate two vec2 into a vec4: { l[0], l[0], h[0], h[1] }
LOCAL vec4 cat22(vec2 h, vec2 l) { return _mm256_insertf128_pd(_mm256_castpd128_pd256(l), h, 1); }

LOCAL vec4 add4(vec4 a, vec4 b)         { return _mm256_add_pd(a,b); }
LOCAL vec4 sub4(vec4 a, vec4 b)         { return _mm256_sub_pd(a,b); }
LOCAL vec4 mul4(vec4 a, vec4 b)         { return _mm256_mul_pd(a,b); }
LOCAL vec4 div4(vec4 a, vec4 b)         { return _mm256_div_pd(a,b); }
LOCAL vec4 sqrt4(vec4 a)                { return _mm256_sqrt_pd(a); }

LOCAL vec4 hadd4(vec4 a, vec4 b)        { return _mm256_hadd_pd(a,b); }

LOCAL vec4 max4(vec4 a, vec4 b)         { return _mm256_max_pd(a,b); }
LOCAL vec4 min4(vec4 a, vec4 b)         { return _mm256_min_pd(a,b); }
LOCAL vec4 and4(vec4 a, vec4 b)         { return _mm256_and_pd(a,b); }
LOCAL vec4 andnot4(vec4 a, vec4 b)      { return _mm256_andnot_pd(a,b); }
LOCAL vec4 abs4(vec4 a)                 { return _mm256_andnot_pd(sgn1111, a); }
LOCAL vec4 flipsign4(vec4 a)            { return _mm256_xor_pd(a, sgn1111); }

/// returns { a[0], b[0], a[2], b[2] }
LOCAL vec4 unpacklo4(vec4 a, vec4 b)    { return _mm256_unpacklo_pd(a,b); }
/// returns { a[1], b[1], a[3], b[3] }
LOCAL vec4 unpackhi4(vec4 a, vec4 b)    { return _mm256_unpackhi_pd(a,b); }

/// returns { a[0], a[0], a[2], a[2] }
LOCAL vec4 duplo4(vec4 a)               { return _mm256_movedup_pd(a); } //_mm256_unpacklo_pd(a,a)
/// returns { a[1], a[1], a[3], a[3] }
LOCAL vec4 duphi4(vec4 a)               { return _mm256_permute_pd(a,15); } //_mm256_unpackhi_pd(a,a)

/// copy a[0] into all elements of dst.
LOCAL vec4 broadcastX(vec4 a) { return _mm256_movedup_pd(_mm256_permute2f128_pd(a, a, 0x00)); }

/* Unused functions:
 LOCAL vec4 loadu22(double const* a, double const* b) { return _mm256_loadu2_m128d(a,b); }
 LOCAL void store22(double* a, double* b, vec4 c) { return _mm256_storeu2_m128d(a,b,c); }
 */

#define permute2f128(a,b,k) _mm256_permute2f128_pd(a,b,k)

/// swap the two 128 bit lanes, return { a[2], a[3], a[0], a[1] }
LOCAL vec4 swap2f128(vec4 a) { return _mm256_permute2f128_pd(a, a, 0x03); }

/// return { a[0], a[1], a[0], a[1] }
LOCAL vec4 duplo2f128(vec4 a) { return _mm256_permute2f128_pd(a, a, 0x20); }
/// return { a[2], a[3], a[2], a[3] }
LOCAL vec4 duphi2f128(vec4 a) { return _mm256_permute2f128_pd(a, a, 0x31); }

/// extract the lower 128 bit lanes, return { a[0], a[1], b[0], b[1] }
LOCAL vec4 unpacklo2f128(vec4 a, vec4 b) { return _mm256_permute2f128_pd(a, b, 0x20); }
/// extract the higher 128 bit lanes, return { a[2], a[3], b[2], b[3] }
LOCAL vec4 unpackhi2f128(vec4 a, vec4 b) { return _mm256_permute2f128_pd(a, b, 0x31); }

LOCAL vec2 permute2(vec2 a) { return _mm_permute_pd(a,1); }

/// return { a[0], a[0], a[3], b[3] }
LOCAL vec4 duplohi4(vec4 a) { return _mm256_permute_pd(a,0b1100); }

#define insertf128(a,b,k)   _mm256_insertf128_pd(a,b,k)
#define permute4(a,k)       _mm256_permute_pd(a,k)
#define shuffle4(a,b,k)     _mm256_shuffle_pd(a,b,k)
#define cmp4(a,b,k)         _mm256_cmp_pd(a,b,k)

/// return { a[0], a[1], a[2], b[3] }
LOCAL vec4 blend31(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b1000); }
/// return { a[0], a[1], b[2], b[3] }
LOCAL vec4 blend22(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b1100); }
/// return { a[0], b[1], b[2], b[3] }
LOCAL vec4 blend13(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b1110); }

/// return { a[0], b[1], a[2], b[3] }
LOCAL vec4 blend0101(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b1010); }
/// return { a[0], b[1], b[2], a[3] }
LOCAL vec4 blend0110(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b0110); }
/// return { a[0], a[1], b[2], a[3] }
LOCAL vec4 blend0010(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b0100); }
/// return { b[0], b[1], a[2], b[3] }
LOCAL vec4 blend1101(vec4 a, vec4 b) { return _mm256_blend_pd(a,b,0b1011); }

/// concatenate, making { ABCD } from a={ AB } b={ CD }
LOCAL vec4 concatenate22(vec2 a, vec2 b) { return _mm256_insertf128_pd(_mm256_castpd128_pd256(a),b,1); }

/// concatenate and shift left by 1 steps, making { BCDE } from a={ ABCD } b={ EFGH }
LOCAL vec4 catshift1(vec4 a, vec4 b) { return _mm256_shuffle_pd(a, _mm256_permute2f128_pd(a, b, 0x21), 0b0101); } // { CDEF }
/// concatenate and shift left by 2 steps, making { CDEF } from a={ ABCD } b={ EFGH }
LOCAL vec4 catshift2(vec4 a, vec4 b) { return _mm256_permute2f128_pd(a, b, 0x21); }
/// concatenate and shift left by 3 steps, making { DEFG } from a={ ABCD } b={ EFGH }
LOCAL vec4 catshift3(vec4 a, vec4 b) { return _mm256_shuffle_pd(_mm256_permute2f128_pd(a, b, 0x21), b, 0b0101); }

/// zero out last scalar
LOCAL vec4 clear4th(vec4 a) { return _mm256_blend_pd(a,_mm256_setzero_pd(),0b1000); }

/// load 4 single precision and convert to double precision
LOCAL vec4 load4d(float const* a) { return _mm256_cvtps_pd(_mm_loadu_ps(a)); }

/// return `neg` if `val < 0` and `pos` otherwise
LOCAL vec4 signselect4(vec4 val, vec4 neg, vec4 pos) { return _mm256_blendv_pd(pos, neg, val); }

#if 0
  LOCAL vec4  load3(double const* a)    { return blend0010(load2Z(a)), broadcast1(a+2)); }
  LOCAL void store3(double* a, vec4 b)  { storeu2(a, getlo(b)); store1(a+2, gethi(b)); }
#else
  //LOCAL vec4  load3(double const* a)  { return _mm256_loadu_pd(a); }
  constexpr __m256i msk1110 = {-1, -1, -1, 0};  // -1 = all bits to 1
  LOCAL vec4  load3(double const* a)    { return _mm256_maskload_pd(a, msk1110); }
  LOCAL void store3(double* a, vec4 b)  { _mm256_maskstore_pd(a, msk1110, b); }
#endif


/// returns the sum of the elements, broadcasted
LOCAL vec4 esum4(vec4 v)
{
    vec4 s = add4(v, swap2f128(v));
    return add4(s, permute4(s, 0b0101));
}

/// returns the dot product of two vectors, broadcasted
LOCAL vec4 dot4(vec4 a, vec4 b)
{
    vec4 m = mul4(a, b);
    vec4 s = add4(m, swap2f128(m));
    return add4(s, permute4(s, 0b0101));
}

/// square of vector norm, broadcasted
LOCAL vec4 normsqr4(vec4 vec)
{
    vec4 m = mul4(vec, vec);
    vec4 s = add4(m, swap2f128(m));
    return add4(s, permute4(s, 0b0101));
}

/// normalize vector
LOCAL vec4 normalize4(vec4 vec)
{
    vec4 m = mul4(vec, vec);
    vec4 s = add4(m, swap2f128(m));
    m = add4(s, permute4(s, 0b0101));
    return div4(vec, sqrt4(m));
}

/// normalize vector to 'n'
LOCAL vec4 normalize4(vec4 vec, double n)
{
    vec4 m = mul4(vec, vec);
    vec4 s = add4(m, swap2f128(m));
    m = add4(s, permute4(s, 0b0101));
    return mul4(vec, div4(set4(n), sqrt4(m)));
}


//----------------------------- Half Precision -------------------------------

/// convert doubles to half-precision floats
LOCAL void convert_to_halfs(size_t cnt, double const* src, uint16_t* dst)
{
    for ( size_t i = 0; i < cnt; i += 4 )
    {
        __m256d d = _mm256_load_pd(src+i);
        __m128 f = _mm256_cvtpd_ps(d);
        __m128i h = _mm_cvtps_ph(f, _MM_FROUND_TO_NEAREST_INT);
        _mm_storel_epi64((__m128i*)(dst+i), h);
    }
}

//---------------------------------- AVX2 --------------------------------------

#if defined(__AVX2__)

#define permute4x64(a,k)    _mm256_permute4x64_pd(a,k)

/// return { X X Y Y } from { X Y }
LOCAL vec4 interleave2(vec2 a) { return _mm256_permute4x64_pd(cast4(a), 0x50); }
/// return { X X Y Y } from { X Y - - }
LOCAL vec4 interleave4(vec4 a) { return _mm256_permute4x64_pd(a, 0x50); }

/// copy a[0] into all elements of dst.
LOCAL vec4 broadcastX(vec2 a)  { return _mm256_broadcastsd_pd(a); }


/// cross product of two 3D vectors ( X Y Z T )
LOCAL vec4 cross4(vec4 a, vec4 b)
{
    // this will not fly because permute accross 128bit boudaries are slow
    vec4 at = permute4x64(a, 0xC9); // Y Z X T
    vec4 bt = permute4x64(b, 0xD2); // Z X Y T
    return sub4(mul4(at,bt), permute4x64(mul4(at,b), 0xC9));
}

// Non-temporal load: an aligned load intruction that bypasses the cache.
LOCAL vec4 streamload4(double const* a) { return _mm256_castsi256_pd(_mm256_stream_load_si256((__m256i const*)a)); }
//LOCAL vec4 streamload4(double const* a) { return _mm256_loadu_pd(a); }

#elif defined(__AVX__)

LOCAL vec4 streamload4(double const* a) { return _mm256_loadu_pd(a); }

/// return { X X Y Y } from { X Y }
LOCAL vec4 interleave2(vec2 a) { return permute4(permute2f128(cast4(a), cast4(a), 0x00), 0b1100); }
/// return { X X Y Y } from { X Y - - }
LOCAL vec4 interleave4(vec4 a) { return permute4(permute2f128(a, a, 0x00), 0b1100); }

#endif

//----------------------------------- FMA --------------------------------------

#if defined(__FMA__)
LOCAL vec4 fmadd4(vec4 a, vec4 b, vec4 c)  { return _mm256_fmadd_pd(a,b,c); }
LOCAL vec4 fnmadd4(vec4 a, vec4 b, vec4 c) { return _mm256_fnmadd_pd(a,b,c); }
#else
LOCAL vec4 fmadd4(vec4 a, vec4 b, vec4 c)  { return _mm256_add_pd(_mm256_mul_pd(a,b), c); }
LOCAL vec4 fnmadd4(vec4 a, vec4 b, vec4 c) { return _mm256_sub_pd(c, _mm256_mul_pd(a,b)); }
#endif
