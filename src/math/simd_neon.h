// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// SIMD for the ARM architecture; Split from simd.h on 12.08.2022


typedef float64x1_t vec1;

LOCAL vec1 load1S(double const* a) { return vld1_f64(a); }
LOCAL void store1(double* a, vec1 b)   { vst1_f64(a, b); }

LOCAL vec1 add1(vec1 a, vec1 b)      { return vadd_f64(a,b); }
LOCAL vec1 sub1(vec1 a, vec1 b)      { return vsub_f64(a,b); }
LOCAL vec1 mul1(vec1 a, vec1 b)      { return vmul_f64(a,b); }
LOCAL vec1 div1(vec1 a, vec1 b)      { return vdiv_f64(a,b); }
/// a * b + c
LOCAL vec1 fmadd1(vec1 a, vec1 b, vec1 c)  { return vfma_f64(c,a,b); }
/// c - a * b
LOCAL vec1 fnmadd1(vec1 a, vec1 b, vec1 c) { return vfms_f64(c,a,b); }



/// Vector holding 2 double precision floats
typedef float64x2_t vec2;

LOCAL void _mm_empty() {};

LOCAL vec2 setr2(double a, double b) { return vec2{a, b}; }
LOCAL vec2 set2(double a)  { return vdupq_n_f64(a); }
LOCAL vec2 setzero2()      { return vdupq_n_f64(0); }
LOCAL vec2 negate2(vec2 a) { return vnegq_f64(a); }

LOCAL vec2 load1(double const* a) { return vec2{a[0], 0}; }
LOCAL vec2 load1Z(double const* a) { return vec2{a[0], 0}; }
LOCAL vec2 load1upper(double const* b) { return vsetq_lane_f64(*b, setzero2(), 1); }
LOCAL vec2 load1lower(vec2 a, double const* b) { return vsetq_lane_f64(*b, a, 0); }
LOCAL vec2 load1upper(vec2 a, double const* b) { return vsetq_lane_f64(*b, a, 1); }

LOCAL vec2 load2(double const* a) { return vld1q_f64(a); }
LOCAL vec2 loadu2(double const* a) { return vld1q_f64(a); }

// load 1 double and duplicate into both positions
LOCAL vec2 loaddup2(double const* a) { return vld1q_dup_f64(a); }

//LOCAL vec2 loadhi2(vec2 a, double const* b) { return _mm_loadh_pd(a,b); }
//LOCAL vec2 loadlo2(vec2 a, double const* b) { return _mm_loadl_pd(a,b); }

LOCAL vec2 cvtsd2(float32x2_t a) { return vcvt_f64_f32(a); }
// convert two single-precision floats in lower registers, to double precision
LOCAL vec2 cvtsd2(float32x4_t a) { return vcvt_f64_f32(vget_low_f32(a)); }

//LOCAL vec2 load1d(float const* a) { return vec2{a[0], 0}; }
LOCAL vec2 load1d(float const* a) { return vcvt_f64_f32(float32x2_t{a[0], 0}); }
// load 2 floats and convert to double
LOCAL vec2 load2d(float const* a) { return vcvt_f64_f32(vld1_f32(a)); }

//------ half-precision, float16
LOCAL float load_half(__fp16 const* a) { float16x4_t h = vld1_f16(a); float32x4_t f = vcvt_f32_f16(h); return f[0]; }

// load 1 float on 2 bytes and convert to double and zero
LOCAL float32x2_t load1_half(__fp16 const* a) { return vget_low_f32(vcvt_f32_f16(float16x4_t{a[0], 0, 0, 0})); }
// load 1 float on 2 bytes and convert to double and zero
LOCAL float32x2_t load2_half(__fp16 const* a) { return vget_low_f32(vcvt_f32_f16(float16x4_t{a[0], a[1], 0, 0})); }

LOCAL vec2 load1h(__fp16 const* a) { return vcvt_f64_f32(load1_half(a)); }
// load 2 floats on 2x2 bytes and convert to double
LOCAL vec2 load2h(__fp16 const* a) { return vcvt_f64_f32(load2_half(a)); }

//------ bfloat16
typedef uint16_t bfloat16;

LOCAL float load_bf16(bfloat16 const* a) { union { float f; uint16_t i[2]; } u; u.i[1]=a[0]; u.i[0]=0; return u.f; }

// load 1 float on 2 bytes and convert to double and zero
LOCAL float32x2_t load1_bf16(bfloat16 const* a) { return uint16x4_t{0, a[0], 0, 0}; }
// load 1 float on 2 bytes and convert to double and zero
LOCAL float32x2_t load2_bf16(bfloat16 const* a) { return uint16x4_t{0, a[0], 0, a[1]}; }

LOCAL vec2 load1b(bfloat16 const* a) { return vcvt_f64_f32(load1_bf16(a)); }
// load 2 floats on 2x2 bytes and convert to double
LOCAL vec2 load2b(bfloat16 const* a) { return vcvt_f64_f32(load2_bf16(a)); }


LOCAL void store1(double* a, vec2 b)   { vst1_f64(a, vget_low_f64(b)); }
LOCAL void store2(double* a, vec2 b)   { vst1q_f64(a, b); }
LOCAL void store1upper(double* a, vec2 b) { vst1_f64(a, vget_high_f64(b)); }

//LOCAL void storedup(double* a, vec2 b) { _mm_store1_pd(a, b); }
//LOCAL void storelo(double* a, vec2 b)  { _mm_store_sd(a, b); }
LOCAL void storeu2(double* a, vec2 b)  { vst1q_f64(a, b); }

LOCAL vec2 duplo2(vec1 a)            { return vcombine_f64(a, a); }
LOCAL vec2 duplo2(vec2 a)            { return vdupq_laneq_f64(a, 0); }
LOCAL vec2 duphi2(vec2 a)            { return vdupq_laneq_f64(a, 1); }

LOCAL vec2 add1(vec2 a, vec2 b)      { return vaddq_f64(a, b); }
LOCAL vec2 sub1(vec2 a, vec2 b)      { return vsubq_f64(a, b); }
LOCAL vec2 mul1(vec2 a, vec2 b)      { return vmulq_f64(a, b); }
LOCAL vec2 div1(vec2 a, vec2 b)      { return vdivq_f64(a, b); }

LOCAL vec2 add2(vec2 a, vec2 b)      { return vaddq_f64(a, b); }
LOCAL vec2 sub2(vec2 a, vec2 b)      { return vsubq_f64(a, b); }
LOCAL vec2 mul2(vec2 a, vec2 b)      { return vmulq_f64(a, b); }
LOCAL vec2 div2(vec2 a, vec2 b)      { return vdivq_f64(a, b); }
LOCAL vec2 hadd2(vec2 a, vec2 b)     { return vpaddq_f64(a, b); }
LOCAL vec2 sqrt2(vec2 a)             { return vsqrtq_f64(a); }

LOCAL vec2 min2(vec2 a, vec2 b)      { return vminq_f64(a, b); }
LOCAL vec2 max2(vec2 a, vec2 b)      { return vmaxq_f64(a, b); }
LOCAL vec2 and2(vec2 a, vec2 b)      { return vandq_s64(a, b); }
LOCAL vec2 andnot2(vec2 a, vec2 b)   { return vbicq_s64(a, b); }
LOCAL vec2 abs2(vec2 a)              { return vbicq_s64(a, float64x2_t{-0.,-0.}); }
LOCAL vec2 flipsign2(vec2 a)         { return veorq_s64(a, float64x2_t{-0.,-0.}); }

/// return { a[0], b[0] }
LOCAL vec2 unpacklo2(vec2 a, vec2 b) { return vzip1q_f64(a, b); }
/// return { a[1], b[1] }
LOCAL vec2 unpackhi2(vec2 a, vec2 b) { return vzip2q_f64(a, b); }
/// return { a[1], a[0] }
LOCAL vec2 swap2(vec2 a)             { return vextq_f64(a, a, 1); }

/// change arguments to a = { a[0] b[0] } and b = { b[0] a[0] }
LOCAL void swap2upper(vec2& a, vec2& b) { a=vzip1q_f64(a, b); b=vzip1q_f64(b, a); }
/// change arguments to a = { b[0] a[0] } and b = { a[0] b[0] }
LOCAL void swap2lower(vec2& a, vec2& b) { vec2 t=b; b=vzip1q_f64(a, b); a=vzip1q_f64(t, a); }

/// concatenate and shift left, returning { BC } from a={ AB } b={ CD }
LOCAL vec2 catshift(vec2 a, vec2 b) { return vextq_f64(a, b, 1); }
/// concatenate and shift left, returning { BC } from a={ AB } b={ CD }
LOCAL vec2 catshift(vec2 a, vec1 b) { return vcombine_f64(vget_high_f64(a), b); }

/// blend to return { low = a[0], high = b[1] }
//LOCAL vec2 blend11(vec2 a, vec2 b) { return vcombine_f64(vget_low_f64(a), vget_high_f64(b)); }
LOCAL vec2 blend11(vec2 a, vec2 b) { return vbslq_f64(int64x2_t{~0,0},a,b); }

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
LOCAL vec2 fmadd1(vec2 a, vec2 b, vec2 c)  { return vfmaq_f64(c,a,b); }
/// c - a * b
LOCAL vec2 fnmadd1(vec2 a, vec2 b, vec2 c) { return vfmsq_f64(c,a,b); }

/// a * b + c
LOCAL vec2 fmadd2(vec2 a, vec2 b, vec2 c)  { return vfmaq_f64(c,a,b); } //vmlaq_f64 will issue mul & add
/// c - a * b
LOCAL vec2 fnmadd2(vec2 a, vec2 b, vec2 c) { return vfmsq_f64(c,a,b); } //vmlsq_f64 will issue mul & sub

//----------------------------- Single Precision -------------------------------

typedef float32x2_t vec2f;

LOCAL vec2f setzero2f() { return vdup_n_f32(0); }

LOCAL vec2f load1f(float const* a) { return vset_lane_f32(a[0], setzero2f(), 0); }
LOCAL vec2f load2f(float const* a) { return vld1_f32(a); }
LOCAL vec2f loaddupf(float const* a) { return vld1_dup_f32(a); }

LOCAL void store1f(float* a, vec2f b) { *a = b[0]; }
LOCAL void store2f(float* a, vec2f b) { vst1_f32(a, b); }

LOCAL vec2f add2f(vec2f a, vec2f b) { return vadd_f32(a,b); }
LOCAL vec2f sub2f(vec2f a, vec2f b) { return vsub_f32(a,b); }
LOCAL vec2f mul2f(vec2f a, vec2f b) { return vmul_f32(a,b); }

/// return { a[0], b[1] }
LOCAL vec2f blend11f(vec2f a, vec2f b) { return vbsl_f32(int32x2_t{~0,0},a,b); }

// returns { a[0], a[0] }
LOCAL vec2f duplo2f(vec2f a) { return vzip1_f32(a, a); }
// returns { a[1], a[1] }
LOCAL vec2f duphi2f(vec2f a) { return vzip2_f32(a, a); }

// returns { a[0], b[0] }
LOCAL vec2f unpacklo2f(vec2f a, vec2f b) { return vzip1_f32(a, b); }
// returns { b[1], a[1] }
LOCAL vec2f unpackhi2f(vec2f a, vec2f b) { return vzip2_f32(a, b); }

/// a * b + c
LOCAL vec2f fmadd2f(vec2f a, vec2f b, vec2f c)  { return vfma_f32(c,a,b); }
/// c - a * b
LOCAL vec2f fnmadd2f(vec2f a, vec2f b, vec2f c) { return vfms_f32(c,a,b); }


//----------------------------- Single Precision -------------------------------
typedef float32x4_t vec4f;
typedef int32x4_t vec4i;

LOCAL vec4i cast4fi(vec4f a) { return vreinterpretq_s32_f32(a); }
LOCAL vec4f cast4if(vec4i a) { return vreinterpretq_f32_s32(a); }

LOCAL vec4f setzero4f() { return vdupq_n_f32(0); }
LOCAL vec4f negate4f(vec4f a) { return vnegq_f32(a); }

LOCAL vec4f broadcast1f(float const* a) { return vld1q_dup_f32(a); }
//LOCAL vec4f load1f(float const* a) { return vec4f{a[0], 0, 0, 0}; }
//LOCAL vec4f load2f(float const* a) { return vec4f{a[0], a[1], 0, 0}; }

LOCAL void store1f(float* a, vec4f b) { *a = b[0]; }
LOCAL void store2f(float* a, vec4f b) { vst1_f32(a, vget_low_f32(b)); }
LOCAL void store3f(float* a, vec4f b) { a[0]=b[0]; a[1]=b[1]; a[2]=b[2]; }
LOCAL void store4f(float* a, vec4f b) { vst1q_f32(a, b); }
LOCAL void storeu4f(float* a, vec4f b) { vst1q_f32(a, b); }

/// convert single and store them in double precision
LOCAL void store2d(double* a, vec2f b) { vst1q_f64(a, vcvt_f64_f32(b)); }

LOCAL vec4f set4f(float a) { return vdupq_n_f32(a); }
LOCAL vec4i set4i(int a)   { return vdupq_n_u32(a); }
LOCAL vec4f set4fi(int a) { return vdupq_n_u32(a); }

LOCAL vec4f load4f(float const* a) { return vld1q_f32(a); }
LOCAL vec4f loadu4f(float const* a) { return vld1q_f32(a); }

//LOCAL vec4f shiftbitsR4(vec4f a, int i) { return vshrq_n_u32((int32x4_t)a, i); }
//LOCAL vec4f shiftbitsL4(vec4f a, int i) { return vshlq_n_u32((int32x4_t)a, i); }
//LOCAL int32_t getlane4i(vec4i a, int i) { return vgetq_lane_u32(a, i); }

#define shiftbitsR4(a,i) vshrq_n_u32((int32x4_t)a,i)
#define shiftbitsL4(a,i) vshlq_n_u32((int32x4_t)a,i)
#define getlane4i(a, i) vgetq_lane_u32(a, i);

/// load integers
LOCAL int32x4_t load4i(int32_t const* a) { return vld1q_s32(a); }
/// convert integer to float:
LOCAL vec4f cvt4if(int32x4_t a) { return vcvtq_f32_s32(a); }

/// load integers and convert to float:
LOCAL vec4f load4if(int32_t const* a) { return vcvtq_f32_s32(vld1q_s32(a)); }
LOCAL vec4f load4uf(uint32_t const* a) { return vcvtq_f32_u32(vld1q_u32(a)); }

LOCAL vec4i negative4f(vec4f a) { return vcltq_f32(a, setzero4f()); }
LOCAL vec4i positive4f(vec4f a) { return vcgtq_f32(a, setzero4f()); }
LOCAL vec4i notpositive4f(vec4f a) { return vcleq_f32(a, setzero4f()); }
LOCAL vec4i greaterequal4f(vec4f a, vec4f b) { return vcgeq_f32(a,b); }
LOCAL vec4i lowerthan4f(vec4f a, vec4f b) { return vcltq_f32(a,b); }

LOCAL vec4f add4f(vec4f a, vec4f b) { return vaddq_f32(a,b); }
LOCAL vec4f sub4f(vec4f a, vec4f b) { return vsubq_f32(a,b); }
LOCAL vec4f mul4f(vec4f a, vec4f b) { return vmulq_f32(a,b); }
LOCAL vec4f div4f(vec4f a, vec4f b) { return vdivq_f32(a,b); }
LOCAL vec4f sqrt4f(vec2 a)          { return vsqrtq_f32(a); }

LOCAL vec4f min4f(vec4f a, vec4f b) { return vminq_f32(a,b); }
LOCAL vec4f max4f(vec4f a, vec4f b) { return vmaxq_f32(a,b); }

LOCAL vec4f hadd4f(vec4f a, vec4f b) { return vpaddq_f32(a, b); }

LOCAL vec4i or4f(vec4i a, vec4i b) { return vorrq_s32(a, b); }
LOCAL vec4i and4f(vec4i a, vec4i b) { return vandq_s32(a, b); }
LOCAL vec4i andnot4f(vec4i a, vec4i b) { return vbicq_s32(a, b); }
LOCAL vec4i and4i(vec4i a, vec4i b) { return vandq_u32(a,b); }

LOCAL vec4f abs4f(vec4f a)          { return vabsq_f32(a); }
LOCAL vec4f flipsign4f(vec4f a)     { return veorq_s32(a, vec4f{-0.,-0.,-0.,-0.}); }

/// select 'b' if 'k==1' and 'a' otherwise
LOCAL vec4f blendv4f(vec4f a, vec4f b, vec4f k) { return vbslq_f32(k,b,a); }
/// return `neg` if `val < 0` and `pos` otherwise
LOCAL vec4f signselect4f(vec4f val, vec4f neg, vec4f pos) { return vbslq_f32(vcleq_f32(val, setzero4f()), neg, pos); }

/// return { a[0], a[1], b[2], b[3] }
LOCAL vec4f blend31f(vec4f a, vec4f b) { return vbslq_f32(int32x4_t{~0,0,0,0},a,b); }
/// return { a[0], b[1], b[2], b[3] }
LOCAL vec4f blend22f(vec4f a, vec4f b) { return vbslq_f32(int32x4_t{~0,~0,0,0},a,b); }
/// return { a[0], a[1], a[2], b[3] }
LOCAL vec4f blend13f(vec4f a, vec4f b) { return vbslq_f32(int32x4_t{~0,~0,~0,0},a,b); }

/// concatenate and shift left
LOCAL vec4f catshift1f(vec4f a, vec4f b) { return vextq_f32(a, b, 1); }
/// concatenate and shift left
LOCAL vec4f catshift2f(vec4f a, vec4f b) { return vextq_f32(a, b, 2); }
/// concatenate and shift left
LOCAL vec4f catshift3f(vec4f a, vec4f b) { return vextq_f32(a, b, 3); }

// returns { a[0], a[0], a[2], a[2] }
LOCAL vec4f duplo4f(vec4f a) { return vzip1q_f32(a, a); }
// returns { a[1], a[1], a[3], a[3] }
LOCAL vec4f duphi4f(vec4f a) { return vzip2q_f32(a, a); }

LOCAL vec2f getlo2f(vec4f a) { return vget_low_f32(a); }
LOCAL vec2f gethi2f(vec4f a) { return vget_high_f32(a); }

// returns { a[0], b[0], a[1], b[1] }
LOCAL vec4f unpacklo4f(vec4f a, vec4f b) { return vzip1q_f32(a, b); }
// returns { a[2], b[2], a[3], b[3] }
LOCAL vec4f unpackhi4f(vec4f a, vec4f b) { return vzip2q_f32(a, b); }

// return { a[0], a[1], b[0], b[1] }
//LOCAL vec4f movelh4f(vec4f a, vec4f b) { return vcombine_f32(vget_low_f32(a), vget_low_f32(b)); }
// return { b[2], b[3], a[2], a[2] }
//LOCAL vec4f movehl4f(vec4f a, vec4f b) { return vcombine_f32(vget_high_f32(b), vget_high_f32(a)); }

// return { a[0], a[1], b[0], b[1] }
LOCAL vec4f movelh4f(vec4f a, vec4f b) { return vtrn1q_f64(a, b); }
// return { b[2], b[3], a[2], a[2] }
LOCAL vec4f movehl4f(vec4f a, vec4f b) { return vtrn2q_f64(b, a); }

// copy a[0] into all elements of destination
LOCAL vec4f broadcastXf(vec4f a) { return vdupq_laneq_f32(a,0); }
// copy a[1] into all elements of destination
LOCAL vec4f broadcastYf(vec4f a) { return vdupq_laneq_f32(a,1); }
// copy a[2] into all elements of destination
LOCAL vec4f broadcastZf(vec4f a) { return vdupq_laneq_f32(a,2); }
// copy a[3] into all elements of destination
LOCAL vec4f broadcastTf(vec4f a) { return vdupq_laneq_f32(a,3); }

// returns { a[0], a[1], a[2], 0 }
LOCAL vec4f clear4th(vec4f a) { return vsetq_lane_f32(0, a, 3); }
// loading 3 elements
LOCAL vec4f load3f(float const* a) { return vec4f{a[0], a[1], a[2], 0}; }
// loading 4 and clearing one
LOCAL vec4f load3fZ(float const* a) { return clear4th(loadu4f(a)); }

/// a * b + c
LOCAL vec4f fmadd4f(vec4f a, vec4f b, vec4f c)  { return vfmaq_f32(c,a,b); }
/// c - a * b
LOCAL vec4f fnmadd4f(vec4f a, vec4f b, vec4f c) { return vfmsq_f32(c,a,b); }

/// set i-th bit in returned value to ( a[i] < b[i] ), for i = {0, 1, 2, 3}
LOCAL int lower_mask4f(vec4f a, vec4f b)
{
    uint32x4_t lowerthan = vcltq_f32(a, b); // all bits set in each component
    constexpr int32x4_t mask = { 1, 2, 4, 8 };
    int res = vaddvq_u32(vandq_u32(lowerthan, mask));
    return res;
}

//----------------------------- Half Precision -------------------------------

/// convert doubles to half-precision floats
LOCAL void convert_to_halfs(size_t cnt, double const* src, __fp16* dst)
{
    #pragma omp simd
    for ( size_t i = 0; i < cnt; i += 4 )
    {
        float32x2_t f = vcvt_f32_f64(vld1q_f64(src+i));
        float32x2_t g = vcvt_f32_f64(vld1q_f64(src+i+2));
        float16x4_t h = vcvt_f16_f32(vcombine_f32(f, g));
        vst1_f16(dst+i, h);
    }
}
