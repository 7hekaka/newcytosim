// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

/** Approximate versions of log() and 1/sqrt() for SIMD vectors */

#if defined(__AVX__)

/// approximate reciprocal square root : 1 / sqrt(x) for 4 doubles
inline vec4 rsqrt4(vec4 x)
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
inline vec8f rsqrt8fi(vec8f x)
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
inline vec4f log4f(vec4f const x) { return _mm_log_ps(x); }
/// natural logarithm, part of Intel's SVML library
inline vec8f log8f(vec8f const x) { return _mm256_log_ps(x); }

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


#define USE_HORNER_RULE 1

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
    const vec8f a1 = set8f(+3.529304993f);
    const vec8f a2 = set8f(-2.461222105f);
    const vec8f a3 = set8f(+1.130626167f);
    const vec8f a4 = set8f(-0.288739945f);
    const vec8f a5 = set8f(+3.110401639e-2f);
    const vec8f F = set8f(-89.970756366f);
    const vec8f G = set8f(0.6931471805f);
    // used to clear negative / NaN arguments:
    vec8f invalid = cmp8f(x, setzero8f(), _CMP_NGT_UQ);
    // extract exponent:
#ifdef __AVX2__
    vec8f a0 = cvt8if(_mm256_srli_epi32(_mm256_castps_si256(x), 23));
#else
    __m128 hi = gethi4f(x);
    __m128 lo = getlo4f(x);
    hi = cvt4if(_mm_srli_epi32(_mm_castps_si128(hi), 23));
    lo = cvt4if(_mm_srli_epi32(_mm_castps_si128(lo), 23));
    vec8f a0 = cat4f(hi, lo);
#endif
    a0 = add8f(mul8f(a0, G), F);
    // clear exponents:
    x = or8f(expo, and8f(mant, x));
    /* Mathematically equivalent polynomial evaluations:
     a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5
     a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5))))
     [a0 + a1*x] + xx*([a2 + a3*x] + xx*[a4 + a5*x]))
     */
#if USE_HORNER_RULE
    vec8f tmp = add8f(mul8f(x, a5), a4);
    tmp = add8f(mul8f(x, tmp), a3);
    tmp = add8f(mul8f(x, tmp), a2);
    tmp = add8f(mul8f(x, tmp), a1);
    tmp = add8f(mul8f(x, tmp), a0);
#else
    // could reduce latency, at the cost of one additional multiplication
    vec8f xx = mul8f(x, x);
    vec8f t01 = add8f(mul8f(x, a1), a0);
    vec8f t23 = add8f(mul8f(x, a3), a2);
    vec8f tmp = add8f(mul8f(x, a5), a4);
    tmp = add8f(mul8f(xx, tmp), t23);
    tmp = add8f(mul8f(xx, tmp), t01);
#endif
    // clear negative arguments:
    return or8f(tmp, invalid);
}

/// Approximate cos+sin by Jacques-Henri Jourdan
/* Correct only in [-pi, pi]
   Absolute error bounded by 5e-5
   Continuous error */
inline void sincosapprox8f(vec8f& S, vec8f& C, const vec8f x)
{
    vec8f xx = mul8f(x, x);
    
    const vec8f c4 = set8f(1.8929864824e-5f);
    const vec8f c3 = set8f(-1.3422947025e-3f);
    const vec8f c2 = set8f(4.1518035216e-2f);
    const vec8f c1 = set8f(-0.4998515820f);
    const vec8f c0 = set8f(1.f);
    
    /* Mathematically equivalent polynomial evaluations:
     a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4
     a0 + x*(a1 + x*(a2 + x*(a3 + x*a4)))
     ([a0 + a1*x] + xx*([a2 + a3*x] + a4*xx)

     a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 + a6*x^6
     a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + a6*x)))))
     ([a0 + a1*x] + xx*[a2 + a3*x]) + [xx*xx]*([a4 + a5*x] + a6*x)
     */
#if USE_HORNER_RULE
    // Horner's rule for 4th order polynom
    C = add8f(mul8f(xx, c4), c3);
    C = add8f(mul8f(xx, C), c2);
    C = add8f(mul8f(xx, C), c1);
    C = add8f(mul8f(xx, C), c0);
#else
    vec8f x4 = mul8f(xx, xx);
    vec8f c01 = add8f(mul8f(xx, c1), c0);
    vec8f c23 = add8f(mul8f(xx, c3), c2);
    C = add8f(mul8f(x4, c4), c23);
    C = add8f(mul8f(x4, C), c01);
#endif
    
    const vec8f s4 = set8f(2.1478401777e-6f);
    const vec8f s3 = set8f(-1.9264918228e-4f);
    const vec8f s2 = set8f(8.3089787513e-3f);
    const vec8f s1 = set8f(-0.1666243672f);
    const vec8f s0 = set8f(0.9999793767f);

#if USE_HORNER_RULE
    S = add8f(mul8f(xx, s4), s3);
    S = add8f(mul8f(xx, S), s2);
    S = add8f(mul8f(xx, S), s1);
    S = add8f(mul8f(xx, S), s0);
#else
    vec8f s01 = add8f(mul8f(xx, s1), s0);
    vec8f s23 = add8f(mul8f(xx, s3), s2);
    S = add8f(mul8f(x4, s4), s23);
    S = add8f(mul8f(x4, S), s01);
#endif
    S = mul8f(x, S);
}


/**
 Initialize CS[] to a circle:
 delta = 2 * PI / cnt;
     CS[2*i  ] = rad * cos(start+i*delta)
     CS[2*i+1] = rad * sin(start+i*delta)
 */
inline void circleAVX(size_t cnt, float CS[], double rad, double start)
{
    const double theta = 2.0 * M_PI / (double)cnt;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const double c2 = c * c - s * s;
    const double s2 = c * s + c * s;

    vec4 cs{ c2, s2,  c2, s2};
    vec4 sc{-s2, c2, -s2, c2};
    
    const double x0 = rad * std::cos(start);
    const double y0 = rad * std::sin(start);
    vec4 pp{x0, y0, c*x0-s*y0, s*x0+c*y0};
    
    float * ptr = CS;
    float * const end = CS + 2 * cnt;
    while ( ptr < end )
    {
        store4d(ptr, pp);
        ptr += 4;
        // apply the rotation matrix
        // x = c * x - s * y;
        // y = s * y + c * y;
        pp = add4(mul4(cs, duplo4(pp)), mul4(sc, duphi4(pp)));
    }
    end[0] = x0;
    end[1] = y0;
}

#endif
