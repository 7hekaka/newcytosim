#pragma once
/**
 * @file  SFMT-avx2.h
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) for Intel AVX2
 *
 * @author Seizh 
 * @note this file is based from SFMT-sse2.h 
 *  (see original SFMT implementation)
 * Copyright (C) 2015 SeizhLab. All rights reserved.
 * 
 * The new BSD License is applied to this software, see LICENSE.txt
 */

#ifndef SFMT_AVX2_H
#define SFMT_AVX2_H

#if SFMT_SL1 < 16
  #error sorry, assumed SFMT_SL1 >= 16.
#endif
#if SFMT_N & 1
  #error sorry, assumed SFMT_N is even number. 
#endif
#if SFMT_POS1 & 1
  #error sorry, assumed SFMT_POS1 is even number. 
#endif

/**
 * parameters used by 256-bit AVX2.
 */
static const w256_t avx2_param_mask = {{SFMT_MSK1, SFMT_MSK2, SFMT_MSK3, SFMT_MSK4,
    SFMT_MSK1, SFMT_MSK2, SFMT_MSK3, SFMT_MSK4}};

static const w128_t avx_mask = {{SFMT_MSK1, SFMT_MSK2, SFMT_MSK3, SFMT_MSK4}};

#ifdef __AVX2__
/**
 * This function represents the recursion formula.
 * @param r an output
 * @param a double 128-bit part of the interal state array 
 * @param b double 128-bit part of the interal state array
 * @param c double 128-bit part of the interal state array
 */
inline static void mm_recursion2(__m256i* r, __m256i a, __m256i b, __m256i c)
{
    __m256i x, y, z;

    x = _mm256_slli_si256(a, SFMT_SL2);
    x = _mm256_xor_si256(x, a);
    y = _mm256_srli_epi32(b, SFMT_SR1);
    y = _mm256_and_si256(y, avx2_param_mask.ymm);
    x = _mm256_xor_si256(x, y);
    z = _mm256_srli_si256(c, SFMT_SR2);
    x = _mm256_xor_si256(x, z);

    /* assume SFMT_SL1 >= 16 */
    z = _mm256_permute2f128_si256(c, x, 0x21);     /* [c.upper, x.lower] */
    z = _mm256_slli_epi32(z, SFMT_SL1);
    x = _mm256_xor_si256(x, z);

    _mm256_store_si256(r, x);
}
#else
inline static void mm_recursion2(__m256i* r, __m256i a, __m256i b, __m256i c)
{
    __m128i la = _mm256_castsi256_si128(a), ua = _mm256_extractf128_si256(a, 1);
    __m128i lb = _mm256_castsi256_si128(b), ub = _mm256_extractf128_si256(b, 1);
    __m128i lc = _mm256_castsi256_si128(c), uc = _mm256_extractf128_si256(c, 1);
    __m128i lx, ly, lz;
    __m128i ux, uy, uz;

    lx = _mm_slli_si128(la, SFMT_SL2);
    ux = _mm_slli_si128(ua, SFMT_SL2);
    lx = _mm_xor_si128(lx, la);
    ux = _mm_xor_si128(ux, ua);
    ly = _mm_srli_epi32(lb, SFMT_SR1);
    uy = _mm_srli_epi32(ub, SFMT_SR1);
    ly = _mm_and_si128(ly, avx_mask.xmm);
    uy = _mm_and_si128(uy, avx_mask.xmm);
    lx = _mm_xor_si128(lx, ly);
    ux = _mm_xor_si128(ux, uy);
    lz = _mm_srli_si128(lc, SFMT_SR2);
    uz = _mm_srli_si128(uc, SFMT_SR2);
    lx = _mm_xor_si128(lx, lz);
    ux = _mm_xor_si128(ux, uz);

    /* assume SFMT_SL1 >= 16 */
    //z = _mm256_permute2f128_si256(c, x, 0x21);     /* [c.upper, x.lower] */
    //lz = uc; uz = lx;
    lz = _mm_slli_epi32(uc, SFMT_SL1);
    uz = _mm_slli_epi32(lx, SFMT_SL1);
    lx = _mm_xor_si128(lx, lz);
    ux = _mm_xor_si128(ux, uz);

    _mm256_store_si256(r, _mm256_set_m128i(ux, lx));
}
#endif

/**
 * This function fills the internal state array with pseudorandom
 * integers.
 * @param sfmt SFMT internal state
 */
void sfmt_gen_rand_all(sfmt_t * sfmt)
{
    //printf("AVX2 sfmt_gen_rand_all\n");
    int i;
    int pos = SFMT_POS1 / 2;
    __m256i r;
    __m256i * pstate = sfmt->state256;
    
    r = _mm256_load_si256(pstate+SFMT_N256-1);
    for (i = 0; i < SFMT_N256-pos; ++i) {
        mm_recursion2(
            pstate+i,
            _mm256_load_si256(pstate+i),
            _mm256_load_si256(pstate+i+pos),
            r
        );
        r = _mm256_load_si256(pstate+i);
    }
    for (; i < SFMT_N256; ++i) {
        mm_recursion2(
            pstate+i,
            _mm256_load_si256(pstate+i),
            _mm256_load_si256(pstate+i+pos-SFMT_N256),
            r
        );
        r = _mm256_load_si256(pstate+i);
    }
}

/**
 * This function fills the user-specified array with pseudorandom
 * integers.
 * @param sfmt SFMT internal state.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
#if ( 0 )
static void gen_rand_array(sfmt_t * sfmt, w128_t * array, unsigned size)
{
    int i, j;
    __m256i r;
    w128_t * pstate = sfmt->state;

    r = _mm256_loadu_si256((__m256i*)&pstate[SFMT_N - 2]);
    for (i = 0; i < SFMT_N - SFMT_POS1; i+=2) {
        mm_recursion2(
            (__m256i*)&array[i], 
            _mm256_load_si256((__m256i*)&pstate[i] ), 
            _mm256_load_si256((__m256i*)&pstate[i + SFMT_POS1]), 
            r
        );
        r = _mm256_load_si256((__m256i*)&array[i]);
    }
    for (; i < SFMT_N; i+=2) {
        mm_recursion2(
            (__m256i*)&array[i], 
            _mm256_load_si256((__m256i*)&pstate[i] ), 
            _mm256_load_si256((__m256i*)&array[i + SFMT_POS1 - SFMT_N]), 
            r
        );
        r = _mm256_loadu_si256((__m256i*)&array[i]);
    }
    for (; i < size - SFMT_N; i+=2) {
        mm_recursion2(
            (__m256i*)&array[i], 
            _mm256_load_si256((__m256i*)&array[i - SFMT_N] ), 
            _mm256_load_si256((__m256i*)&array[i + SFMT_POS1 - SFMT_N]), 
            r
        );
        r = _mm256_load_si256((__m256i*)&array[i]);
        
    }
    for (j = 0; j < 2 * SFMT_N - size; j += 2) {
        _mm256_store_si256((__m256i*)&pstate[j], _mm256_load_si256((__m256i*)&array[j + size - SFMT_N]));
    }
    for (; i < size; i+=2, j+=2) {
        mm_recursion2(
            (__m256i*)&array[i], 
            _mm256_load_si256((__m256i*)&array[i - SFMT_N] ), 
            _mm256_load_si256((__m256i*)&array[i + SFMT_POS1 - SFMT_N]), 
            r
        );
        r = _mm256_load_si256((__m256i*)&array[i]);
        _mm256_store_si256((__m256i*)&pstate[j], r);
    }
}
#else
static void gen_rand_array(sfmt_t * sfmt, w128_t * input, unsigned size128)
{
    unsigned size = size128 / 2;
    unsigned pos = SFMT_POS1 / 2;
    unsigned i, j;
    __m256i r;
    __m256i * pstate = sfmt->state256;
    __m256i * array = (__m256i*)input;

    r = _mm256_loadu_si256(&pstate[SFMT_N256 - 1]);
    for (i = 0; i < SFMT_N256-pos; ++i) {
        mm_recursion2(
            &array[i],
            _mm256_load_si256(&pstate[i]),
            _mm256_load_si256(&pstate[i + pos]),
            r
        );
        r = _mm256_load_si256(&array[i]);
    }
    for (; i < SFMT_N256; ++i) {
        mm_recursion2(
            &array[i],
            _mm256_load_si256(&pstate[i]),
            _mm256_load_si256(&array[i + pos - SFMT_N256]),
            r
        );
        r = _mm256_loadu_si256(&array[i]);
    }
    for (; i < size - SFMT_N256; ++i) {
        mm_recursion2(
            &array[i],
            _mm256_load_si256(&array[i - SFMT_N256] ),
            _mm256_load_si256(&array[i + pos - SFMT_N256]),
            r
        );
        r = _mm256_load_si256(&array[i]);
    }
    for (j = 0; j < SFMT_N - size; ++j) {
        _mm256_store_si256(&pstate[j], _mm256_load_si256(&array[j + size - SFMT_N256]));
    }
    for (; i < size; ++i, ++j) {
        mm_recursion2(
            &array[i],
            _mm256_load_si256(&array[i - SFMT_N256] ),
            _mm256_load_si256(&array[i + pos - SFMT_N256]),
            r
        );
        r = _mm256_load_si256(&array[i]);
        _mm256_store_si256(&pstate[j], r);
    }
}
#endif
#endif
