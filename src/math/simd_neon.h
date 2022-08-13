// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Started on Monday 5 June 2018, which was a very nice day in Strasbourg

/* Unfinished work: needs to be converted to ARM's NEON intrinsics */

/// Vector holding 2 double precision floats
typedef float64x1_t vec1;
typedef float64x2_t vec2;

constexpr vec2 sgn11 = {-0.0, -0.0};

// Attention: the second value returned by load1() is not set and will be garbage!
LOCAL vec2 load1(double const* a)           { return vsetq_lane_f64(*a, vdupq_n_f64(0), 0); }
LOCAL vec2 load1Z(double const* a)          { return vsetq_lane_f64(*a, vdupq_n_f64(0), 0); }

// Attention: load2() does not initialize the upper AVX registers
LOCAL vec2 load2(double const* a)           { return vld1q_f64(a); }

// unaligned load
LOCAL vec2 loadu2(double const* a)          { return vld1q_f64(a); }

// load 1 double and duplicate
LOCAL vec2 loaddup2(double const* a)        { return vld1q_dup_f64(a); }

//LOCAL vec2 loadhi2(vec2 a, double const* b) { return _mm_loadh_pd(a,b); }
//LOCAL vec2 loadlo2(vec2 a, double const* b) { return _mm_loadl_pd(a,b); }

LOCAL vec2 cvtsd2(float32x2_t a) { return vcvt_f64_f32(a); }
// convert two single-precision floats in lower registers, to double precision
LOCAL vec2 cvtsd2(float32x4_t a) { return vcvt_f64_f32(vget_low_f32(a)); }

// load 1 float and convert to double and zero
LOCAL vec2 load1d(float const* a) { return vcvt_f64_f32(vset_lane_f32(*a, vdup_n_f32(0), 0)); }
// load 2 floats and convert to double
LOCAL vec2 load2d(float const* a) { return vcvt_f64_f32(vld1_f32(a)); }

LOCAL void store1(double* a, vec1 b)   { vst1_f64(a, b); }
LOCAL void store1(double* a, vec2 b)   { vst1_f64(a, vget_low_f64(b)); }
LOCAL void store2(double* a, vec2 b)   { vst1q_f64(a, b); }
//LOCAL void storedup(double* a, vec2 b) { _mm_store1_pd(a, b); }
//LOCAL void storelo(double* a, vec2 b)  { _mm_store_sd(a, b); }
LOCAL void storeu2(double* a, vec2 b)  { vst1q_f64(a, b); }

LOCAL vec2 duplo2(vec2 a)            { return vdupq_laneq_f64(a, 0); }
LOCAL vec2 duphi2(vec2 a)            { return vdupq_laneq_f64(a, 1); }

LOCAL vec1 mul1(vec1 a, vec1 b)      { return vmul_f64(a,b); }
LOCAL vec1 div1(vec1 a, vec1 b)      { return vdiv_f64(a,b); }
LOCAL vec1 add1(vec1 a, vec1 b)      { return vadd_f64(a,b); }
LOCAL vec1 sub1(vec1 a, vec1 b)      { return vsub_f64(a,b); }

LOCAL vec2 mul2(vec2 a, vec2 b)      { return vmulq_f64(a, b); }
LOCAL vec2 div2(vec2 a, vec2 b)      { return vdivq_f64(a, b); }
LOCAL vec2 add2(vec2 a, vec2 b)      { return vaddq_f64(a, b); }
LOCAL vec2 sub2(vec2 a, vec2 b)      { return vsubq_f64(a, b); }
LOCAL vec2 hadd2(vec2 a, vec2 b)     { return vpaddq_f64(a, b); }

LOCAL vec2 sqrt2(vec2 a)             { return vsqrtq_f64(a); }
LOCAL vec2 max2(vec2 a, vec2 b)      { return vmaxq_f64(a, b); }
LOCAL vec2 min2(vec2 a, vec2 b)      { return vminq_f64(a, b); }
LOCAL vec2 and2(vec2 a, vec2 b)      { return vandq_s64(a, b); }
LOCAL vec2 andnot2(vec2 a, vec2 b)   { return vbicq_s64(a, b); }
LOCAL vec2 abs2(vec2 a)              { return vbicq_s64(sgn11, a); }
LOCAL vec2 flipsign2(vec2 a)         { return veorq_s64(a, sgn11); }

LOCAL vec2 setr2(double a, double b) { double u[2] = { a, b }; return vld1q_f64(u); }
LOCAL vec2 set2(double a)            { return vdupq_n_f64(a); }
LOCAL vec2 setzero2()                { return vdupq_n_f64(0); }

/// return { a[0], b[0] }
LOCAL vec2 unpacklo2(vec2 a, vec2 b) { return vzip1q_f64(a, b); }
/// return { a[1], b[1] }
LOCAL vec2 unpackhi2(vec2 a, vec2 b) { return vzip2q_f64(a, b); }
LOCAL vec2 swap2(vec2 a)             { return vextq_f64(a, a, 1); }

/// concatenate and shift left, returning { BC } from a={ AB } b={ CD }
LOCAL vec2 catshift(vec2 a, vec2 b) { return vextq_f64(b, a, 1); } //OKAY?

/// blend to return { low = a[0], high = b[1] }
LOCAL vec2 blend11(vec2 a, vec2 b) { return vcombine_f64(vget_low_f64(a), vget_high_f64(a)); }

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


//---------------------------- Multiply-Accumulate -----------------------------

/// a * b + c
LOCAL vec1 fmadd1(vec1 a, vec1 b, vec1 c)  { return vfma_f64(c,a,b); }
/// c - a * b
LOCAL vec1 fnmadd1(vec1 a, vec1 b, vec1 c) { return vfms_f64(c,a,b); }


/// a * b + c
LOCAL vec2 fmadd2(vec2 a, vec2 b, vec2 c)  { return vfmaq_f64(c,a,b); }
/// c - a * b
LOCAL vec2 fnmadd2(vec2 a, vec2 b, vec2 c) { return vfmsq_f64(c,a,b); }

