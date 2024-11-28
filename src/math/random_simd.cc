// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.


#if USE_SIMD

/**
 This folds the corners of [-1, 1] x [-1, 1] that are outside the unit circle,
 to map these points back onto the original square. For randomly distributed points,
 this increases the number of points within the unit circle by a factor 3-2*sqrt(2)
 without changing the property of being equidistributed within the unit circle.
 */
static inline void fold_corners4f(vec4f& x, vec4f& y)
{
    // test if point is close to a corner, crossing lines |x|+|y| > sqrt(2)
    vec4f mut = greaterequal4f(add4f(abs4f(x), abs4f(y)), set4f(M_SQRT2));
    // coordinates of nearest corner, scaled by S: copysign(S, x)
    constexpr float S = float(M_SQRT2 + 1.0);
    constexpr vec4f pos{S, S, S, S};
    constexpr vec4f neg{-S, -S, -S, -S};
    vec4f cx = signselect4f(x, neg, pos);
    vec4f cy = signselect4f(y, neg, pos);
    // subtract corner to recover a square of side 2, covering the circle
    cx = sub4f(mul4f(pos, x), cx);
    cy = sub4f(mul4f(pos, y), cy);
    // apply mutation, if selected:
    x = blendv4f(x, cx, mut);
    y = blendv4f(y, cy, mut);
}


static real * makeGaussians_SIMD(real dst[], size_t cnt, const uint32_t arg[])
{
    const int32_t * src = (int32_t*)arg;
    const int32_t * end = src + cnt;
    const vec4f fac = set4f(0x1p-31);
    const vec4f two = set4fi(2);
    const vec4f half = set4f(-0.5f);

    while ( src < end )
    {
        vec4f x = mul4f(fac, load4if(src));
        vec4f y = mul4f(fac, load4if(src+4));
        src += 8;
        fold_corners4f(x, y); // increases from 490 to 574 expected!
        vec4f n = add4f(mul4f(x,x), mul4f(y,y));
        // the log() will be negative if `n < 1`, returning 0 if n == 0:
        vec4f l = log_approx4f(n);
        // set valid[i] to 2 whenever `l[i] < 0`, and 0 otherwise:
        vec4i valid = cast4f4i(and4f(negative4f(l), two));
        // n = sqrt( log(n) / ( -0.5 * n ) );
        n = sqrt4f(div4f(l, mul4f(half, n)));
        vec4f z;
        z = mul4f(n, x);
        y = mul4f(n, y);
        // place corresponding X and Y values next to each other:
        x = unpacklo4f(z, y);
        y = unpackhi4f(z, y);
#if REAL_IS_DOUBLE
        // convert 8 single-precision values
        store2d(dst, getlo2f(x)); dst += getlane4i(valid, 0);
        store2d(dst, gethi2f(x)); dst += getlane4i(valid, 1);
        store2d(dst, getlo2f(y)); dst += getlane4i(valid, 2);
        store2d(dst, gethi2f(y)); dst += getlane4i(valid, 3);
#else
        // convert 8 single-precision values
        store2f(dst, getlo2f(x)); dst += getlane4i(valid, 0);
        store2f(dst, gethi2f(x)); dst += getlane4i(valid, 1);
        store2f(dst, getlo2f(y)); dst += getlane4i(valid, 2);
        store2f(dst, gethi2f(y)); dst += getlane4i(valid, 3);
#endif
    }
    return dst;
}


/** This follows Box-Muller's method with approximate Log and Sin/Cos functions */
static real * makeGaussiansBM_SIMD(real dst[], size_t cnt, const uint32_t* arg)
{
    const uint32_t* src = arg;
    const uint32_t* end = src + cnt;

    constexpr float tmp(M_PI*0x1p-31);
    constexpr vec4f PI{tmp, tmp, tmp, tmp};

    //vec4f mm{1.f, 1.f, 1.f, 1.f};
    while ( src < end )
    {
        // generate angle in ]-PI, PI[:
        vec4f t = mul4f(PI, load4if((int32_t*)src)); //load 32-bit signed integer as float
        vec4f x, y;
        sincos_approx4f(x, y, t);

#ifdef __ARM_NEON__
        //load 32-bit unsigned integer as float
        vec4f n = minuslog_approx4f32(load4uf(src+4));
#else
        //load 32-bit signed integer as float
        vec4f n = minuslog_approx4f31(abs4f(load4if((int32_t*)(src+4))));
#endif
        //mm = min4f(mm, n);
        // this should be equivalent to: n = sqrt( -2 * log(R) )
        n = sqrt4f(n);
        x = mul4f(n, x);
        y = mul4f(n, y);
#if REAL_IS_DOUBLE
        // convert 8 single-precision values
        store2d(dst  , getlo2f(x));
        store2d(dst+2, getlo2f(y));
        store2d(dst+4, gethi2f(x));
        store2d(dst+6, gethi2f(y));
#else
        // store 8 single-precision values
        store4f(dst, x);
        store4f(dst+4, y);
#endif
        src += 8;
        dst += 8;
    }
    /*
    mm = min4f(mm, duphi4f(mm));
    float m = std::min(mm[0], mm[1]);
    static float om = 1.0;
    if ( m < om ) {
        printf("min: %.12f\n", m);
        om = m;
    }
     */
    return dst;
}

/// compute approximate exponential derivates
/** Using -log(1-R) where R is randomly distributed in [0, 1[ */
static float* makeExponentials_SIMD(float dst[], size_t cnt, const uint32_t* arg)
{
    const uint32_t * src = arg;
    const uint32_t * end = src + cnt;
    
    const vec4f half{0.5f, 0.5f, 0.5f, 0.5f};
    while ( src < end )
    {
#ifdef __ARM_NEON__
        //load 32-bit unsigned integer as float
        vec4f x = minuslog_approx4f32(load4uf(src));
#else
        //load 32-bit signed integer as float
        vec4f x = minuslog_approx4f31(abs4f(load4if((int32_t*)(src))));
#endif
        x = mul4f(x, half);
        src += 4;
        store4f(dst, x);
        dst += 4;
    }
    return dst;
}


#endif


#if defined(__AVX__)


/**
 This folds the corners of [-1, 1] x [-1, 1] that are outside the unit circle,
 to map these points back onto the original square. For randomly distributed points,
 this increases the number of points within the unit circle by a factor 3-2*sqrt(2)
 without changing the property of being equidistributed within the unit circle.
 */
static inline void fold_corners8f(vec8f& x, vec8f& y)
{
    // test if point is close to corner, separated by line |x|+|y| > sqrt(2)
    vec8f mut = cmp8f(add8f(abs8f(x), abs8f(y)), set8f(M_SQRT2), _CMP_GE_OQ);
    // coordinates of nearest corner, scaled: copysign(S, x)
    constexpr float S = float(M_SQRT2 + 1.0);
    const vec8f pos = set8f(S), neg = set8f(-S);
    vec8f cx = blendv8f(pos, neg, x); //use the sign of 'x'
    vec8f cy = blendv8f(pos, neg, y); //use the sign of 'y'
    // subtract corner to recover a square of side 2, covering the circle
    cx = sub8f(mul8f(pos, x), cx);
    cy = sub8f(mul8f(pos, y), cy);
    // apply mutation, if selected:
    x = blendv8f(x, cx, mut);
    y = blendv8f(y, cy, mut);
}


/** This follows Box-Muller's method with approximate Log and Sin/Cos functions */
static real * makeGaussiansBM_AVX(real dst[], size_t cnt, const uint32_t* arg)
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8;

    const vec8f eps = set8f(0x1p-31);
    const vec8f two = set8f(2.0f);
    const vec8f PI = set8f(M_PI*0x1p-31);

    while ( src < end )
    {
        // generate norm in ]0, 1]: eps * float(uint32)
        vec8f n = mul8f(eps, abs8f(load8if(src++)));
        // generate angle in ]-PI, PI[:
        vec8f t = mul8f(PI, load8if(src++));
        // transform norm:
        n = sqrt8f(mul8f(abs8f(log_approx8f(n)), two));
        vec8f x, y;
        sincos_approx8f(x, y, t);
        x = mul8f(n, x);
        y = mul8f(n, y);
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4d(dst   , getlo4f(x));
        store4d(dst+4 , gethi4f(x));
        store4d(dst+8 , getlo4f(y));
        store4d(dst+12, gethi4f(y));
#else
        // store 16 single-precision values
        store8f(dst, x);
        store8f(dst+8, y);
#endif
        dst += 16;
    }
    return dst;
}


/// pack array by removing all 'not-a-number's in arrray [s, e[
/** Attention: this may not work with option `-ffast-math` */
static inline real * remove_nan_pairs(real * s, real * e)
{
    //real * start = s;
    e -= 2;
    if ( s < e )
    {
        // find the next `nan` going upward:
        while ( s[0] == s[0] )
        {
            s += 2;
            if ( s >= e )
                break;
        }
        // skip `nan` values going downward:
        while ( e[0] != e[0] )
        {
            e -= 2;
            if ( e <= s )
                break;
        }
        //printf("- %4lu  %4lu\n", s-start, e-start);
        // swap numbers over, putting NaNs toward the end of array
        std::swap(s[0], e[0]);
        std::swap(s[1], e[1]);
        s += 2;
        e -= 2;
    }
    while ( s < e )
    {
        // find the next `nan` going upward:
        while ( s[0] == s[0] )
            s += 2;
        // skip `nan` values going downward:
        while ( e[0] != e[0] )
            e -= 2;
        if ( s >= e )
            break;
        //printf("  %4lu  %4lu\n", s-start, e-start);
        // copy numbers over:
        s[0] = e[0];
        s[1] = e[1];
        s += 2;
        e -= 2;
    }
    return s;
}

static real * makeGaussians_AVX(real dst[], size_t cnt, const uint32_t* arg)
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8;

    const vec8f one = set8f(1.0f);
    const vec8f two = set8fi(2);
    const vec8f fac = set8f(0x1p-31);
    const vec8f half = set8f(-0.5f);

    while ( src < end )
    {
        vec8f x = mul8f(fac, load8if(src++));
        vec8f y = mul8f(fac, load8if(src++));
        fold_corners8f(x, y); // increases from 490 to 574 expected!
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        // set valid[i] to 2 whenever 'n[i] < 1.0', and 0 otherwise:
        vec8i valid = _mm256_castps_si256(and8f(lowerthan8f(n, one), two));
        //w = std::sqrt( -2 * std::log(n) / n );
        n = sqrt8f(div8f(log_approx8f(n), mul8f(half, n)));
        y = mul8f(n, y);
        n = mul8f(n, x);
        // place corresponding X and Y values next to each other:
        x = unpacklo8f(n, y);
        y = unpackhi8f(n, y);
        vec4f a = getlo4f(x);
        vec4f b = getlo4f(y);
        vec4f c = gethi4f(x);
        vec4f d = gethi4f(y);
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store2d(dst, getlo2f(a)); dst += getlane8i(valid, 0);
        store2d(dst, gethi2f(a)); dst += getlane8i(valid, 1);
        store2d(dst, getlo2f(b)); dst += getlane8i(valid, 2);
        store2d(dst, gethi2f(b)); dst += getlane8i(valid, 3);
        store2d(dst, getlo2f(c)); dst += getlane8i(valid, 4);
        store2d(dst, gethi2f(c)); dst += getlane8i(valid, 5);
        store2d(dst, getlo2f(d)); dst += getlane8i(valid, 6);
        store2d(dst, gethi2f(d)); dst += getlane8i(valid, 7);
#else
        // store 16 single-precision values
        store2f(dst, getlo2f(a)); dst += getlane8i(valid, 0);
        store2f(dst, gethi2f(a)); dst += getlane8i(valid, 1);
        store2f(dst, getlo2f(b)); dst += getlane8i(valid, 2);
        store2f(dst, gethi2f(b)); dst += getlane8i(valid, 3);
        store2f(dst, getlo2f(c)); dst += getlane8i(valid, 4);
        store2f(dst, gethi2f(c)); dst += getlane8i(valid, 5);
        store2f(dst, getlo2f(d)); dst += getlane8i(valid, 6);
        store2f(dst, gethi2f(d)); dst += getlane8i(valid, 7);
#endif
    }
    return dst;
}


#endif

