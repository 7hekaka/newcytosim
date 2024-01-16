// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Started on Monday 5 June 2018, which was a very nice day in Strasbourg

/// Vector holding 2 double precision floats
typedef __m128d vec2;

constexpr vec2 sgn11{-0.0, -0.0};

// Attention: the second value returned by load1() is not set and will be garbage!
LOCAL vec2 load1(double const* a)           { return _mm_load_sd(a); }
LOCAL vec2 load1Z(double const* a)          { return _mm_loadl_pd(_mm_setzero_pd(), a); }

// Attention: load2() does not initialize the upper AVX registers
LOCAL vec2 load2(double const* a)           { return _mm_load_pd(a); }

// unaligned load
LOCAL vec2 loadu2(double const* a)          { return _mm_loadu_pd(a); }

// load 1 double and duplicate into both positions
LOCAL vec2 loaddup2(double const* a)        { return _mm_loaddup_pd(a); }

LOCAL vec2 loadhi2(vec2 a, double const* b) { return _mm_loadh_pd(a,b); }
LOCAL vec2 loadlo2(vec2 a, double const* b) { return _mm_loadl_pd(a,b); }

// convert two single-precision floats in lower registers, to double precision
LOCAL vec2 cvtsd2(__m128 a) { return _mm_cvtps_pd(a); }

// load 1 float and convert to double and zero
LOCAL vec2 load1d(float const* a) { return _mm_cvtps_pd(_mm_load_ss(a)); }
// load 2 floats and convert to double
LOCAL vec2 load2d(float const* a) { return _mm_cvtps_pd(_mm_castsi128_ps(_mm_loadl_epi64((__m128i*)a))); }

LOCAL void store1(double* a, vec2 b)   { _mm_store_sd(a, b); }
LOCAL void store2(double* a, vec2 b)   { _mm_store_pd(a, b); }
LOCAL void storedup(double* a, vec2 b) { _mm_store1_pd(a, b); }
LOCAL void storelo(double* a, vec2 b)  { _mm_store_sd(a, b); }
LOCAL void storeu2(double* a, vec2 b)  { _mm_storeu_pd(a, b); }

LOCAL vec2 duplo2(vec2 a)            { return _mm_movedup_pd(a); }
LOCAL vec2 duphi2(vec2 a)            { return _mm_shuffle_pd(a, a, 0b11); }

LOCAL vec2 add1(vec2 a, vec2 b)      { return _mm_add_sd(a,b); }
LOCAL vec2 sub1(vec2 a, vec2 b)      { return _mm_sub_sd(a,b); }
LOCAL vec2 mul1(vec2 a, vec2 b)      { return _mm_mul_sd(a,b); }
LOCAL vec2 div1(vec2 a, vec2 b)      { return _mm_div_sd(a,b); }

LOCAL vec2 add2(vec2 a, vec2 b)      { return _mm_add_pd(a,b); }
LOCAL vec2 sub2(vec2 a, vec2 b)      { return _mm_sub_pd(a,b); }
LOCAL vec2 mul2(vec2 a, vec2 b)      { return _mm_mul_pd(a,b); }
LOCAL vec2 div2(vec2 a, vec2 b)      { return _mm_div_pd(a,b); }

LOCAL vec2 sqrt2(vec2 a)             { return _mm_sqrt_pd(a); }
LOCAL vec2 hadd2(vec2 a, vec2 b)     { return _mm_hadd_pd(a,b); }

LOCAL vec2 min2(vec2 a, vec2 b)      { return _mm_min_pd(a,b); }
LOCAL vec2 max2(vec2 a, vec2 b)      { return _mm_max_pd(a,b); }
LOCAL vec2 and2(vec2 a, vec2 b)      { return _mm_and_pd(a,b); }
LOCAL vec2 andnot2(vec2 a, vec2 b)   { return _mm_andnot_pd(a,b); }
LOCAL vec2 abs2(vec2 a)              { return _mm_andnot_pd(sgn11, a); }
LOCAL vec2 flipsign2(vec2 a)         { return _mm_xor_pd(a, sgn11); }

LOCAL vec2 set2(double a, double b)  { return _mm_set_pd(a, b); }
LOCAL vec2 set2(double a)            { return _mm_set1_pd(a); }
LOCAL vec2 setr2(double a, double b) { return _mm_setr_pd(a,b); }
LOCAL vec2 setzero2()                { return _mm_setzero_pd(); }

/// return { a[0], b[0] }
LOCAL vec2 unpacklo2(vec2 a, vec2 b) { return _mm_unpacklo_pd(a,b); }
/// return { a[1], b[1] }
LOCAL vec2 unpackhi2(vec2 a, vec2 b) { return _mm_unpackhi_pd(a,b); }
/// return { a[1], a[0] }
LOCAL vec2 swap2(vec2 a)             { return _mm_shuffle_pd(a, a, 0b01); }

/// concatenate and shift left, returning { BC } from a={ AB } b={ CD }
LOCAL vec2 catshift(vec2 a, vec2 b) { return _mm_shuffle_pd(a, b, 0b01); }

/// blend to return { low = a[0], high = b[1] }
//LOCAL vec2 blend11(vec2 a, vec2 b) { return _mm_shuffle_pd(a, b, 0b10); }
LOCAL vec2 blend11(vec2 a, vec2 b) { return _mm_move_sd(b, a); }

#define cmp2(a,b,k) _mm_cmp_pd(a,b,k)

/// returns the sum of the elements, broadcasted
LOCAL vec2 esum2(vec2 v)
{
    return add2(v, swap2(v));
}

/// returns the dot product of two vectors, broadcasted
LOCAL vec2 dot2(vec2 a, vec2 b)
{
    vec2 p = mul2(a, b);
    return add2(p, swap2(p));
}

/// square of vector norm, broadcasted
LOCAL vec2 normsqr2(vec2 vec)
{
    vec2 p = mul2(vec, vec);
    return add2(p, swap2(p));
}

/// normalize vector
LOCAL vec2 normalize2(vec2 vec)
{
    vec2 p = mul2(vec, vec);
    vec2 s = add2(p, swap2(p));
    return div2(vec, sqrt2(s));
}

/// normalize vector to 'n'
LOCAL vec2 normalize2(vec2 vec, double n)
{
    vec2 p = mul2(vec, vec);
    vec2 s = add2(p, swap2(p));
    return mul2(vec, div2(set2(n), sqrt2(s)));
}


#if defined(__SSE4_1__)

/// return `neg` if `val < 0` and `pos` otherwise
LOCAL vec2 signselect2(vec2 val, vec2 neg, vec2 pos) { return _mm_blendv_pd(pos, neg, val); }

#endif

#if defined(__FMA__)
LOCAL vec2 fmadd1(vec2 a, vec2 b, vec2 c)  { return _mm_fmadd_sd(a,b,c); }  // a * b + c
LOCAL vec2 fnmadd1(vec2 a, vec2 b, vec2 c) { return _mm_fnmadd_sd(a,b,c); } // c - a * b

LOCAL vec2 fmadd2(vec2 a, vec2 b, vec2 c)  { return _mm_fmadd_pd(a,b,c); }
LOCAL vec2 fnmadd2(vec2 a, vec2 b, vec2 c) { return _mm_fnmadd_pd(a,b,c); }
#else
/// Erzatz Fused Multiply-Add functions
LOCAL vec2 fmadd1(vec2 a, vec2 b, vec2 c)  { return _mm_add_sd(_mm_mul_sd(a,b), c); }
LOCAL vec2 fnmadd1(vec2 a, vec2 b, vec2 c) { return _mm_sub_sd(c, _mm_mul_sd(a,b)); }

LOCAL vec2 fmadd2(vec2 a, vec2 b, vec2 c)  { return _mm_add_pd(_mm_mul_pd(a,b), c); }
LOCAL vec2 fnmadd2(vec2 a, vec2 b, vec2 c) { return _mm_sub_pd(c, _mm_mul_pd(a,b)); }
#endif
