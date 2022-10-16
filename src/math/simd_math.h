// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#if defined(__INTEL_COMPILER)

/// natural logarithm, part of Intel's SVML library
inline vec4f log4f(vec4f const x) { return _mm_log_ps(x); }
/// natural logarithm, part of Intel's SVML library
inline vec8f log8f(vec8f const x) { return _mm256_log_ps(x); }

#endif


#if USE_SIMD
/// Approximate natural logarithm by Jacques-Henri Jourdan
/**
 Absolute error bounded by 1e-5 for normalized inputs
 Returns a finite number for +inf input
 Returns -inf for nan and <= 0 inputs.
 Continuous error.
 SIMD by FJN 12.01.2021 derived from:
 http://gallium.inria.fr/blog/fast-vectorizable-math-approx/
 */
inline vec4f logapprox4f(vec4f x)
{
    // masks:
    const vec4f mant = set4fi(0x007FFFFF);
    const vec4f expo = set4fi(0x3F800000);
    // polynomial coefficients
    const vec4f a = set4f(+3.529304993f);
    const vec4f b = set4f(-2.461222105f);
    const vec4f c = set4f(+1.130626167f);
    const vec4f d = set4f(-0.288739945f);
    const vec4f e = set4f(+3.110401639e-2f);
    const vec4f f = set4f(-89.970756366f);
    const vec4f g = set4f(0.6931471805f);
    // used to clear negative / NaN arguments:
    vec4f invalid = notpositive4f(x);
    // extract exponent:
    vec4f cst = cvt4if(shiftbitsR4(x, 23));
    cst = fmadd4f(cst, g, f);
    // reset exponents to '127':
    x = or4f(expo, and4f(mant, x));
    // evaluate polynom:
    vec4f tmp = fmadd4f(x, e, d);
    tmp = fmadd4f(x, tmp, c);
    tmp = fmadd4f(x, tmp, b);
    tmp = fmadd4f(x, tmp, a);
    tmp = fmadd4f(x, tmp, cst);
    // set invalid arguments to all 1s which is not-a-number:
    return or4f(tmp, invalid);
    return tmp;
}

#endif


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


/// Approximate natural logarithm by Jacques-Henri Jourdan
/**
 Absolute error bounded by 1e-5 for normalized inputs
 Returns a finite number for +inf input
 Returns -inf for nan and <= 0 inputs.
 Continuous error.
 SIMD by FJN 12.01.2021 derived from:
 */
inline vec8f logapprox8f(vec8f x)
{
    // masks:
    const vec8f mant = set8fi(0x007FFFFF);
    const vec8f expo = set8fi(0x3F800000);
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
#if defined(__AVX2__)
    vec8f a0 = cvt8if(shiftbitsR8(x, 23));
#else
    vec4f h = cvt4if(shiftbitsR4(gethi4f(x), 23));
    vec4f l = cvt4if(shiftbitsR4(getlo4f(x), 23));
    vec8f a0 = cat44f(l, h);
#endif
    a0 = add8f(mul8f(a0, G), F);
    // clear exponents:
    x = or8f(expo, and8f(mant, x));
    /* Mathematically equivalent polynomial evaluations:
     a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5
     a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5))))
     [a0 + a1*x] + xx*([a2 + a3*x] + xx*[a4 + a5*x]))
     */
    vec8f tmp = add8f(mul8f(x, a5), a4);
    tmp = add8f(mul8f(x, tmp), a3);
    tmp = add8f(mul8f(x, tmp), a2);
    tmp = add8f(mul8f(x, tmp), a1);
    tmp = add8f(mul8f(x, tmp), a0);
    // set invalid arguments to all 1s which is not-a-number:
    return or8f(tmp, invalid);
}

/// Approximate cos+sin by Jacques-Henri Jourdan
/**
   This is correct only for x in [-pi, pi]
   Absolute error bounded is by 5e-5
   Continuous error
 SIMD by FJN 12.01.2021 derived from:
 */
inline void sincosapprox8f(vec8f& S, vec8f& C, const vec8f x)
{
    const vec8f c4 = set8f(1.8929864824e-5f);
    const vec8f c3 = set8f(-1.3422947025e-3f);
    const vec8f c2 = set8f(4.1518035216e-2f);
    const vec8f c1 = set8f(-0.4998515820f);
    const vec8f c0 = set8f(1.f);
    
    const vec8f s4 = set8f(2.1478401777e-6f);
    const vec8f s3 = set8f(-1.9264918228e-4f);
    const vec8f s2 = set8f(8.3089787513e-3f);
    const vec8f s1 = set8f(-0.1666243672f);
    const vec8f s0 = set8f(0.9999793767f);

    vec8f xx = mul8f(x, x);

    /* Mathematically equivalent polynomial evaluations:
     a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4
     a0 + x*(a1 + x*(a2 + x*(a3 + x*a4)))
     ([a0 + a1*x] + xx*([a2 + a3*x] + a4*xx)

     a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 + a6*x^6
     a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + a6*x)))))
     ([a0 + a1*x] + xx*[a2 + a3*x]) + [xx*xx]*([a4 + a5*x] + a6*x)
     */
    // Horner's rule for 4th order polynoms
    S = add8f(mul8f(xx, s4), s3);
    C = add8f(mul8f(xx, c4), c3);
    S = add8f(mul8f(xx, S), s2);
    C = add8f(mul8f(xx, C), c2);
    S = add8f(mul8f(xx, S), s1);
    C = add8f(mul8f(xx, C), c1);
    S = add8f(mul8f(xx, S), s0);
    C = add8f(mul8f(xx, C), c0);
    S = mul8f(x, S);
}

#endif


#if defined(__SSE3__)

/**
 Initialize ptr[] to a circle:
 delta = 2 * PI / cnt;
     ptr[2*i  ] = rad * cos(start+i*theta) + cX
     ptr[2*i+1] = rad * sin(start+i*theta) + cY
 */
void set_arc_SEE(size_t cnt, float ptr[], float rad, float start,
                 float delta, float cX, float cY)
{
    const float c0 = cosf(start);
    const float s0 = sinf(start);
    const float c1 = cosf(delta);
    const float s1 = sinf(delta);
    const float C = c1*c1 - s1*s1;
    const float S = c1*s1 + c1*s1;
    
    vec4f CS { C, S,  C, S};
    vec4f SC {-S, C, -S, C};

    vec4f cc { cX, cY, cX, cY };
    vec4f xx { c0, s0, c0*c1-s0*s1, s0*c1+s1*c0 };
    xx = mul4f(set4f(rad), xx);

    float * const end = ptr + 2 * cnt - 2;
    while ( ptr < end )
    {
        storeu4f(ptr, add4f(xx, cc));
        // apply the rotation matrix
        // x = c * x - s * y;
        // y = s * x + c * y;
        xx = add4f(mul4f(CS, duplo4f(xx)), mul4f(SC, duphi4f(xx)));
        ptr += 4;
    }
    storeu4f(ptr, add4f(xx, cc));
}

#endif
