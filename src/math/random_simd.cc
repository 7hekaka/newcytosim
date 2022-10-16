// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.


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


#if USE_SIMD

/**
 This folds the corners of [-1, 1] x [-1, 1] that are outside the unit circle,
 to map these points back onto the original square. For randomly distributed points,
 this increases the number of points within the unit circle by a factor 3-2*sqrt(2)
 without changing the property of being equidistributed within the unit circle.
 */
static inline void fold_corners4f(vec4f& x, vec4f& y)
{
    // test if point is close to corner, separated by line |x|+|y| > sqrt(2)
    vec4f mut = greaterequal4f(add4f(abs4f(x), abs4f(y)), set4f(M_SQRT2));
    // coordinates of nearest corner, scaled: copysign(S, x)
    constexpr float S = M_SQRT2 + 1.0f;
    const vec4f pos = set4f(S), neg = set4f(-S);
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
    const vec4i * src = reinterpret_cast<const vec4i*>(arg);
    const vec4f fac = set4f(0x1p-31);
    const vec4f one = set4f(1.0f);
    const vec4f half = set4f(-0.5f);

    vec4i const* end = src + cnt / 4;
    while ( src < end )
    {
        vec4f x = mul4f(fac, load4if(src));
        vec4f y = mul4f(fac, load4if(src+1));
        src += 2;
        fold_corners4f(x, y); // increases from 490 to 574 expected!
        vec4f n = add4f(mul4f(x,x), mul4f(y,y));
        // set valid[i] to 2 whenever 'n[i] < 1.0', and 0 otherwise:
        vec4i valid = shiftbitsL4(shiftbitsR4(lowerthan4f(n, one), 31), 1);
        n = sqrt4f(div4f(logapprox4f(n), mul4f(half, n)));
        x = mul4f(n, x);
        y = mul4f(n, y);
        // place corresponding X and Y values next to each other:
        n = unpacklo4f(x, y);
        y = unpackhi4f(x, y);
#if REAL_IS_DOUBLE
        // convert 8 single-precision values
        store2d(dst, getlo2f(n));
        dst += getlane4i(valid, 0);
        store2d(dst, gethi2f(n));
        dst += getlane4i(valid, 1);
        store2d(dst, getlo2f(y));
        dst += getlane4i(valid, 2);
        store2d(dst, gethi2f(y));
        dst += getlane4i(valid, 3);
#else
        // convert 8 single-precision values
        store2f(dst, getlo2f(n));
        dst += getlane4i(valid, 0);
        store2f(dst, gethi2f(n));
        dst += getlane4i(valid, 1);
        store2f(dst, getlo2f(y));
        dst += getlane4i(valid, 2);
        store2f(dst, gethi2f(y));
        dst += getlane4i(valid, 3);
#endif
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
    constexpr float S = M_SQRT2 + 1.0f;
    const vec8f ss = set8f(S), mm = set8f(-S);
    vec8f cx = blendv8f(ss, mm, x); //use the sign of 'x'
    vec8f cy = blendv8f(ss, mm, y); //use the sign of 'y'
    // subtract corner to recover a square of side 2, covering the circle
    cx = sub8f(mul8f(ss, x), cx);
    cy = sub8f(mul8f(ss, y), cy);
    // apply mutation, if selected:
    x = blendv8f(x, cx, mut);
    y = blendv8f(y, cy, mut);
}


/** This follows Box-Muller's method with approximate Log and Sin/Cos functions */
static real * makeGaussians_AVXBM(real dst[], size_t cnt, const uint32_t* arg)
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8;

    const vec8f eps = set8f(0x1p-31);
    const vec8f one = set8f(1.0f);
    const vec8f two = set8f(-2.0f);
    const vec8f PI = set8f(M_PI*0x1p-31);

    size_t i = 0;
    while ( src < end )
    {
        // generate norm in ]0, 1]: 1 - eps * float(uint32)
        vec8f n = sub8f(one, mul8f(eps, abs8f(load8if(src++))));
        // generate angle in ]-PI, PI[:
        vec8f t = mul8f(PI, load8if(src++));
        // transform norm:
        n = sqrt8f(mul8f(logapprox8f(n), two));
        vec8f x, y;
        sincosapprox8f(x, y, t);
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
        i += 16;
    }
    return dst;
}


static real * makeGaussians_AVX1(real dst[], size_t cnt, const uint32_t* arg)
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8;

    const vec8f fac = set8f(0x1p-31);
    const vec8f half = set8f(-0.5f);

    real * d = dst;
    while ( src < end )
    {
        vec8f x = mul8f(fac, load8if(src++));
        vec8f y = mul8f(fac, load8if(src++));
        fold_corners8f(x, y); // increases from 490 to 574 expected!
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        //w = std::sqrt( -2 * std::log(n) / n );
        n = sqrt8f(div8f(logapprox8f(n), mul8f(half, n)));
        x = mul8f(n, x);
        y = mul8f(n, y);
        // place corresponding X and Y values next to each other:
        n = unpacklo8f(x, y);
        y = unpackhi8f(x, y);
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4d(d   , getlo4f(n));
        store4d(d+4 , getlo4f(y));
        store4d(d+8 , gethi4f(n));
        store4d(d+12, gethi4f(y));
#else
        // store 16 single-precision values
        store8f(d  , n);
        store8f(d+8, y);
#endif
        d += 16;
    }
    return remove_nan_pairs(dst, d);
    return d;
}



static inline void sort_nans(vec8f& x, vec8f& y)
{
    vec8f u = x;
    vec8f k = isnan8f(x);
    x = blendv8f(u, y, k);
    y = blendv8f(y, u, k);
}

/**
 Calculates Gaussian-distributed, single precision random number,
 using SIMD AVX instructions
 Array `dst` should be able to hold as many 32-bit numbers as `src`.
 if 'real==float', for 256 bits of input, this produces ~64*PI bits of numbers.
 if 'real==double', this produces more output bits than input!

 The function used to calculate logarithm on SIMD data is part of the
 Intel SVML library, and is provided by the Intel compiler.

 F. Nedelec 02.01.2017
 */
static real * makeGaussians_AVX2(real dst[], size_t cnt, const uint32_t* arg)
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8 - 2;

    const vec8f eps = set8f(0x1p-31);
    const vec8f half = set8f(-0.5f);

    real * d = dst;  //cnt=78=39*2 // cnt*8 = 640
    real * e = dst+4*(cnt-2);
    real * f = dst+6*(cnt-2);
    while ( src < end )
    {
        // generate 16 random floats in [-1, 1]:
        vec8f x = mul8f(eps, load8if(src++));
        vec8f y = mul8f(eps, load8if(src++));
        fold_corners8f(x, y);
        // calculate norm:
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        // n = sqrt( -2 * log(n) / n )
        n = sqrt8f(div8f(logapprox8f(n), mul8f(half, n)));
        //n = sqrt8f(div8f(log8f(n), mul8f(half, n)));
        //n = rsqrt8f(div8f(mul8f(half, n), logapprox8f(n)));
        x = mul8f(n, x);
        y = mul8f(n, y);
        // place corresponding X and Y values next to each other:
        n = unpacklo8f(x, y);
        y = unpackhi8f(x, y);
        // swap the NaNs to move them to 'y':
        vec8f k = isnan8f(n);
        x = blendv8f(n, y, k);
        y = blendv8f(y, n, k);
#if 1
        // generate 16 random floats in [-1, 1]:
        vec8f z = mul8f(eps, load8if(src++));
        vec8f t = mul8f(eps, load8if(src++));
        fold_corners8f(z, t);
        vec8f p = add8f(mul8f(z,z), mul8f(t,t));
        p = sqrt8f(div8f(logapprox8f(p), mul8f(half, p)));
        z = mul8f(p, z);
        t = mul8f(p, t);
        p = unpacklo8f(z, t);
        t = unpackhi8f(z, t);
        // swap the NaNs to move them to 't':
        vec8f l = isnan8f(p);
        z = blendv8f(p, t, l);
        t = blendv8f(t, p, l);
        //sorting network for 4 inputs [[1-2][3-4][1 3][2-4][2 3]]
        sort_nans(x, z);
        sort_nans(y, t);
        sort_nans(y, z);
        //t = permute8f(t);
        //sort_nans(z, t);
#endif
#if 0   // can push more NaNs out by swapping within the lanes:
        vec8f k = isnan8f(u);
        x = blendv8f(u, y, k);
        y = blendv8f(y, u, k);
        // swap the NaNs to move them to 'y':
        u = permute44f(x);
        k = isnan8f(u);
        x = blendv8f(u, y, k);
        y = blendv8f(y, u, k);
        // swap the NaNs to move them to 'y':
        u = swap2f128(x);
        k = isnan8f(u);
        x = blendv8f(u, y, k);
        y = blendv8f(y, u, k);
        // swap the NaNs to move them to 'y':
        u = permute44f(x);
        k = isnan8f(u);
        x = blendv8f(u, y, k);
        y = blendv8f(y, u, k);
#endif
#if REAL_IS_DOUBLE
        // convert to single-precision values
        store4d(d   , getlo4f(x));
        store4d(d+4 , gethi4f(x));
        store4d(d+8 , getlo4f(y));
        store4d(d+12, gethi4f(y));
        store4d(e   , getlo4f(z));
        store4d(e+4 , gethi4f(z));
        store4d(f   , getlo4f(t));
        store4d(f+4 , gethi4f(t));
#else
        // store 16 single-precision values
        store8f(d, x);
        store8f(d+8, y);
        store8f(e, z);
        store8f(f, t);
#endif
        d += 16;
        e += 8;
        f += 8;
    }
    return remove_nan_pairs(dst, f);
    return dst+8*cnt;
}


/// compute exponential derivates
static real* makeExponentials_AVX(real dst[], size_t cnt, const uint32_t* arg)
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8;

    const vec8f eps = set8f(0x1p-31f);
    const vec8f one = set8f(1.0f);
    
    while ( src < end )
    {
        // generate random floats in ]0, 1]:
        vec8f x = sub8f(one, mul8f(eps, abs8f(cvt8if(load8i(src)))));
        x = abs8f(logapprox8f(x));
#if REAL_IS_DOUBLE
        // convert 8 single-precision values
        store4d(dst  , getlo4f(x));
        store4d(dst+4, gethi4f(x));
#else
        // store 8 single-precision values
        store8f(dst, x);
#endif
        dst += 8;
        ++src;
    }
    return dst;
}


#endif

