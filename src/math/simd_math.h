
#if defined(__AVX__)


/// approximate reciprocal square root : 1 / sqrt(x) for 4 doubles
inline static vec4 rsqrt4(vec4 x)
{
    vec4 t = _mm256_set1_pd(1.5);
    vec4 h = _mm256_mul_pd(x, _mm256_set1_pd(0.5));
    x = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)));
    // refinement: u  <-  1.5 * u - 0.5 * x * u^3
    vec4 a = _mm256_mul_pd(x, x);
    vec4 b = _mm256_mul_pd(h, x);
    vec4 c = _mm256_mul_pd(t, x);
#if defined(__FMA__)
    x = _mm256_fnmadd_pd(a, b, c);
#else
    x = _mm256_sub_pd(c, _mm256_mul_pd(a, b));
#endif
    // refinement: u  <-  1.5 * u - 0.5 * x * u^3
    a = _mm256_mul_pd(x, x);
    b = _mm256_mul_pd(h, x);
    c = _mm256_mul_pd(t, x);
#if defined(__FMA__)
    return _mm256_fnmadd_pd(a, b, c);
#else
    return _mm256_sub_pd(c, _mm256_mul_pd(a, b));
#endif
}


/// approximate reciprocal square root : 1 / sqrt(x) for 8 floats
inline static vec8f rsqrt8fi(vec8f x)
{
    vec8f h = _mm256_mul_ps(x, _mm256_set1_ps(0.5f));
    x = _mm256_rsqrt_ps(x);
    // refinement: u  <-  1.5 * u - 0.5 * x * u^3
    vec8f a = _mm256_mul_ps(x, x);
    vec8f b = _mm256_mul_ps(h, x);
    vec8f c = _mm256_mul_ps(_mm256_set1_ps(1.5f), x);
#if defined(__FMA__)
    return _mm256_fnmadd_ps(a, b, c);
#else
    return _mm256_sub_ps(c, _mm256_mul_ps(a, b));
#endif
}


// natural logarithm:
#if defined(__INTEL_COMPILER)

inline static vec8f log8f(vec8f const x) { return _mm256_log_ps(x); }

#else


inline vec4f log4f(vec4f xxx)
{
    // masks:
    const vec4f mant = _mm_castsi128_ps(_mm_set1_epi32(0x007fffff));
    const vec4f expo = _mm_castsi128_ps(_mm_set1_epi32(0x3f800000));
    // polynomial coefficients
    const vec4f a = _mm_set1_ps(+3.529304993f);
    const vec4f b = _mm_set1_ps(-2.461222105f);
    const vec4f c = _mm_set1_ps(+1.130626167f);
    const vec4f d = _mm_set1_ps(-0.288739945f);
    const vec4f e = _mm_set1_ps(+3.110401639e-2f);
    const vec4f f = _mm_set1_ps(-89.970756366f);
    const vec4f g = _mm_set1_ps(0.6931471805f);
    // used to clear negative / NaN arguments:
    vec4f invalid = _mm_cmp_ps(xxx, setzero4f(), _CMP_NGT_US);
    // extract exponent:
    vec4f cst = cvt4i(_mm_srli_epi32(_mm_castps_si128(xxx), 23));
    cst = add4f(mul4f(cst, g), f);
    // clear exponents:
    xxx = or4f(expo, and4f(mant, xxx));
    // evaluate polynom:
    vec4f tmp = add4f(mul4f(xxx, e), d);
    tmp = add4f(mul4f(xxx, tmp), c);
    tmp = add4f(mul4f(xxx, tmp), b);
    tmp = add4f(mul4f(xxx, tmp), a);
    tmp = add4f(mul4f(xxx, tmp), cst);
    // clear negative arguments:
    return or4f(tmp, invalid);
}


/*
 Absolute error bounded by 1e-5 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error.
 By Jacques-Henri Jourdan, SIMD by Francois Nedelec 12.01.2021
 */
inline vec8f log8f(vec8f xxx)
{
    // masks:
    const vec8f mant = _mm256_castsi256_ps(_mm256_set1_epi32(0x007fffff));
    const vec8f expo = _mm256_castsi256_ps(_mm256_set1_epi32(0x3f800000));
    // polynomial coefficients
    const vec8f a = _mm256_set1_ps(+3.529304993f);
    const vec8f b = _mm256_set1_ps(-2.461222105f);
    const vec8f c = _mm256_set1_ps(+1.130626167f);
    const vec8f d = _mm256_set1_ps(-0.288739945f);
    const vec8f e = _mm256_set1_ps(+3.110401639e-2f);
    const vec8f f = _mm256_set1_ps(-89.970756366f);
    const vec8f g = _mm256_set1_ps(0.6931471805f);
    // used to clear negative / NaN arguments:
    vec8f invalid = _mm256_cmp_ps(xxx, setzero8f(), _CMP_NGT_US);
    // extract exponent:
#ifdef __AVX2__
    vec8f cst = cvt8i(_mm256_srli_epi32(_mm256_castps_si256(xxx), 23));
#else
    __m128 hi = _mm256_extractf128_ps(xxx, 1);
    __m128 lo = _mm256_extractf128_ps(xxx, 0);
    hi = _mm_cvtepi32_ps(_mm_srli_epi32(_mm_castps_si128(hi), 23));
    lo = _mm_cvtepi32_ps(_mm_srli_epi32(_mm_castps_si128(lo), 23));
    vec8f cst = _mm256_set_m128(hi, lo);
#endif
    cst = add8f(mul8f(cst, g), f);
    // clear exponents:
    xxx = or8f(expo, and8f(mant, xxx));
    // evaluate polynom:
    vec8f tmp = add8f(mul8f(xxx, e), d);
    tmp = add8f(mul8f(xxx, tmp), c);
    tmp = add8f(mul8f(xxx, tmp), b);
    tmp = add8f(mul8f(xxx, tmp), a);
    tmp = add8f(mul8f(xxx, tmp), cst);
    // clear negative arguments:
    return or8f(tmp, invalid);
}

#endif

#endif
