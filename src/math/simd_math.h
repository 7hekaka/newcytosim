// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

/** Approximate versions of log() and 1/sqrt() for SIMD vectors */

#if defined(__AVX__)

/// approximate reciprocal square root : 1 / sqrt(x) for 4 doubles
inline static vec4 rsqrt4(vec4 x)
{
    vec4 t = set4(1.5);
    vec4 h = mul4(x, set4(0.5));
    x = cvt4sd(rsqrt4f(cvt4ds(x)));
    // refinement: u  <-  1.5 * u - 0.5 * x * u^3
    vec4 a = mul4(x, x);
    vec4 b = mul4(h, x);
    vec4 c = mul4(t, x);
    x = sub4(c, mul4(a, b));
    // refinement: u  <-  1.5 * u - 0.5 * x * u^3
    a = mul4(x, x);
    b = mul4(h, x);
    c = mul4(t, x);
    return sub4(c, mul4(a, b));
}


/// approximate reciprocal square root : 1 / sqrt(x) for 8 floats
inline static vec8f rsqrt8fi(vec8f x)
{
    vec8f h = mul8f(x, set8f(0.5f));
    x = rsqrt8f(x);
    // refinement: u  <-  1.5 * u - 0.5 * x * u^3
    vec8f a = mul8f(x, x);
    vec8f b = mul8f(h, x);
    vec8f c = mul8f(set8f(1.5f), x);
    return sub8f(c, mul8f(a, b));
}


#if defined(__INTEL_COMPILER)

/// natural logarithm, part of Intel's SVML library
inline static vec4f log4f(vec4f const x) { return _mm_log_ps(x); }
/// natural logarithm, part of Intel's SVML library
inline static vec8f log8f(vec8f const x) { return _mm256_log_ps(x); }

#endif


/// Approximate natural logarithm by Jacques-Henri Jourdan
/**
 Absolute error bounded by 1e-5 for normalized inputs
 Returns a finite number for +inf input
 Returns -inf for nan and <= 0 inputs.
 Continuous error.
 SIMD by Francois Nedelec 12.01.2021
 */
inline vec4f logapprox4f(vec4f x)
{
    // masks:
    const vec4f mant = set4fi(0x007fffff);
    const vec4f expo = set4fi(0x3f800000);
    // polynomial coefficients
    const vec4f a = set4f(+3.529304993f);
    const vec4f b = set4f(-2.461222105f);
    const vec4f c = set4f(+1.130626167f);
    const vec4f d = set4f(-0.288739945f);
    const vec4f e = set4f(+3.110401639e-2f);
    const vec4f f = set4f(-89.970756366f);
    const vec4f g = set4f(0.6931471805f);
    // used to clear negative / NaN arguments:
    vec4f invalid = cmp4f(x, setzero4f(), _CMP_NGT_UQ);
    // extract exponent:
    vec4f cst = cvt4if(_mm_srli_epi32(_mm_castps_si128(x), 23));
    cst = add4f(mul4f(cst, g), f);
    // clear exponents:
    x = or4f(expo, and4f(mant, x));
    // evaluate polynom:
    vec4f tmp = add4f(mul4f(x, e), d);
    tmp = add4f(mul4f(x, tmp), c);
    tmp = add4f(mul4f(x, tmp), b);
    tmp = add4f(mul4f(x, tmp), a);
    tmp = add4f(mul4f(x, tmp), cst);
    // clear negative arguments:
    return or4f(tmp, invalid);
}


/// Approximate natural logarithm by Jacques-Henri Jourdan
/**
 Absolute error bounded by 1e-5 for normalized inputs
 Returns a finite number for +inf input
 Returns -inf for nan and <= 0 inputs.
 Continuous error.
 SIMD by Francois Nedelec 12.01.2021
 */
inline vec8f logapprox8f(vec8f x)
{
    // masks:
    const vec8f mant = set8fi(0x007fffff);
    const vec8f expo = set8fi(0x3f800000);
    // polynomial coefficients
    const vec8f a = set8f(+3.529304993f);
    const vec8f b = set8f(-2.461222105f);
    const vec8f c = set8f(+1.130626167f);
    const vec8f d = set8f(-0.288739945f);
    const vec8f e = set8f(+3.110401639e-2f);
    const vec8f f = set8f(-89.970756366f);
    const vec8f g = set8f(0.6931471805f);
    // used to clear negative / NaN arguments:
    vec8f invalid = cmp8f(x, setzero8f(), _CMP_NGT_UQ);
    // extract exponent:
#ifdef __AVX2__
    vec8f cst = cvt8if(_mm256_srli_epi32(_mm256_castps_si256(x), 23));
#else
    __m128 hi = gethi4f(x);
    __m128 lo = getlo4f(x);
    hi = cvt4if(_mm_srli_epi32(_mm_castps_si128(hi), 23));
    lo = cvt4if(_mm_srli_epi32(_mm_castps_si128(lo), 23));
    vec8f cst = cat4f(hi, lo);
#endif
    cst = add8f(mul8f(cst, g), f);
    // clear exponents:
    x = or8f(expo, and8f(mant, x));
    // evaluate polynom:
    vec8f tmp = add8f(mul8f(x, e), d);
    tmp = add8f(mul8f(x, tmp), c);
    tmp = add8f(mul8f(x, tmp), b);
    tmp = add8f(mul8f(x, tmp), a);
    tmp = add8f(mul8f(x, tmp), cst);
    // clear negative arguments:
    return or8f(tmp, invalid);
}

#endif
