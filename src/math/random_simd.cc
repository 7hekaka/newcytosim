// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.


// pack array by removing all 'not-a-number's in arrray [s, e[
static inline real * remove_nans(real * s, real * e)
{
    while ( s < e )
    {
        --e;
        // find the next `nan` going upward:
        while ( !std::isnan(*s) )
        {
            if ( ++s > e )
                return s;
        }
        // skip `nan` values going downward:
        while ( std::isnan(*e) )
        {
            if ( --e <= s )
                return s;
        }
        // copy number over:
        *s++ = *e;
        //*e = 0.0;
    }
    return s;
}


#if defined(__AVX__)
/*
 // generate random floats in [1, 2]:
 const vec8f mant = _mm256_castsi256_ps(_mm256_set1_epi32(0x807fffff));
 const vec8f expo = _mm256_castsi256_ps(_mm256_set1_epi32(0x3f800000));
 const vec8f ninf = _mm256_castsi256_ps(_mm256_set1_epi32(0xff800000)); // -INFINITY
 const vec8f one = _mm256_set1_ps(1.0f);
 x = or8f(expo, and8f(mant, cast8f(load8si(src++))));
 y = or8f(expo, and8f(mant, cast8f(load8si(src++))));

inline static vec8f magic8f(vec8f xxx)
{
    // bad approximation of sqrt(x*ln(x)) for x in [1, 2]
    const vec8f one = set8f(1.0f);
    const vec8f alpha = set8f(-2*M_LN2);
    return sqrt8f(mul8f(alpha, mul8f(xxx, sub8f(xxx, one))));
}
 */


/// calculate  sqrt( -2 * log(n) / n )
static inline vec8f sqrtlogdiv8f(const vec8f n)
{
    const vec8f half = set8f(-0.5);
    return sqrt8f(div8f(log8f(n), mul8f(half, n)));
    //return rsqrt8f(div8f(mul8f(half, n), log8f(n)));
    //return rsqrt8f(mul8f(mul8f(half, n), rcp8f(log8f(n))));
}


static real * gauss_fill_AVX0(real dst[], size_t cnt, const __m256i src[])
{
    const vec8f fac = set8f(TWO_POWER_MINUS_31);

    real * d = dst;
    __m256i const* end = src + cnt;
    while ( src < end )
    {
        vec8f x = mul8f(fac, cvt8if(load8si(src++)));
        vec8f y = mul8f(fac, cvt8if(load8si(src++)));
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        //w = std::sqrt( -2 * std::log(n) / n );
        n = sqrtlogdiv8f(n);
        x = mul8f(n, x);
        y = mul8f(n, y);
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4(d   , cvt4sd(getlo4f(x)));
        store4(d+4 , cvt4sd(getlo4f(y)));
        store4(d+8 , cvt4sd(gethi4f(x)));
        store4(d+12, cvt4sd(gethi4f(y)));
#else
        // store 16 single-precision values
        store8f(d  , x);
        store8f(d+8, y);
#endif
        d += 16;
    }
    return remove_nans(dst, d);
    return d;
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
static real * gauss_fill_AVX(real dst[], size_t cnt, const __m256i src[])
{
    const vec8f eps = set8f(0x1p-31); //TWO_POWER_MINUS_31
    const vec8f half = set8f(-0.5);
    __m256i const* end = src + cnt;

    real * d = dst;
    real * e = dst+8*(cnt-2);
    while ( src < end )
    {
        // generate 16 random floats in [-1, 1]:
        vec8f x = mul8f(eps, cvt8if(load8si(src++)));
        vec8f y = mul8f(eps, cvt8if(load8si(src++)));
        // calculate norm:
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        // complex transformation here:
        //n = sqrt8f(div8f(logapprox8f(n), mul8f(half, n)));
        //n = rsqrt8f(div8f(mul8f(half, n), logapprox8f(n)));
        n = sqrtlogdiv8f(n);
        x = mul8f(n, x);
        y = mul8f(n, y);
        // place corresponding X and Y values next to each other:
        vec8f t = unpacklo8f(x, y);
        y = unpackhi8f(x, y);
        // swap the NaNs to move them to 'y':
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
#if 0   // can push more NaNs out by swapping within the lanes:
        // swap the NaNs to move them to 'y':
        t = permute44f(x);
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
        // swap the NaNs to move them to 'y':
        t = swap2f128(x);
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
        // swap the NaNs to move them to 'y':
        t = permute44f(x);
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
#endif
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4(d  , cvt4sd(getlo4f(x)));
        store4(d+4, cvt4sd(gethi4f(x)));
        store4(e  , cvt4sd(getlo4f(y)));
        store4(e+4, cvt4sd(gethi4f(y)));
#else
        // store 16 single-precision values
        store8f(d, x);
        store8f(e, y);
#endif
        d += 8;
        e -= 8;
    }
    //return remove_nans(dst, d) - dst;
    return dst+8*cnt;
}

#endif
