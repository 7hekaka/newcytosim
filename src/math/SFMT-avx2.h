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
 * https://github.com/seizh/sfmt_mod
 *
 * This has been rewritten by Francois Nedelec on 16.01.2021
 */

#ifndef SFMT_AVX2_H
#define SFMT_AVX2_H


#define PRINT_TO_CHECK_AVX_CODE 0

#if PRINT_TO_CHECK_AVX_CODE
void sfmt_gen_rand_all_generic(sfmt_t * sfmt)
{
    int i;
    w128_t *r1, *r2;

    r1 = &sfmt->state[SFMT_N - 2];
    r2 = &sfmt->state[SFMT_N - 1];
    for (i = 0; i < SFMT_N - SFMT_POS1; i++) {
        do_recursion(&sfmt->state[i], &sfmt->state[i],
                     &sfmt->state[i + SFMT_POS1], r1, r2);
        r1 = r2;
        r2 = &sfmt->state[i];
    }
    for (i = SFMT_N - SFMT_POS1; i < SFMT_N; i++) {
        do_recursion(&sfmt->state[i], &sfmt->state[i],
                     &sfmt->state[i + SFMT_POS1 - SFMT_N], r1, r2);
        r1 = r2;
        r2 = &sfmt->state[i];
    }
}

#include <stdio.h>
void print_state(w128_t const* ptr)
{
    for ( int u = 0; u < SFMT_N/2; u += 2 )
        printf("%016llx %016llx %016llx %016llx\n",
               ptr[2*u].u64[0], ptr[2*u].u64[1], ptr[2*u+1].u64[0], ptr[2*u+1].u64[1]);
    printf("\n");
}
#endif


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


/**
 * This function represents the recursion formula.
 * @param r an output
 * @param a double 128-bit part of the interal state array 
 * @param b double 128-bit part of the interal state array
 * @param c double 128-bit part of the interal state array
 */
inline static __m256i mm256_recursion(__m256i a, __m256i b, __m256i c)
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
    return _mm256_xor_si256(x, z);
}

/**
 * This function fills the internal state array with pseudorandom
 * integers.
 * @param sfmt SFMT internal state
 */
void sfmt_gen_rand_all(sfmt_t * sfmt)
{
    //printf("AVX2 sfmt_gen_rand_all\n");
    int i;
    __m256i r;
    __m256i * pstate = (__m256i*)sfmt->state;
#if PRINT_TO_CHECK_AVX_CODE
    uint64_t tmp[SFMT_N64];
    memcpy(tmp, sfmt->state, 8*SFMT_N64);
    sfmt_gen_rand_all_generic(sfmt);
    print_state(sfmt->state);
    memcpy(sfmt->state, tmp, 8*SFMT_N64);
#endif
    r = pstate[SFMT_N256-1];
    for (i = 0; i < (SFMT_N-SFMT_POS1)/2; ++i) {
        r = mm256_recursion(pstate[i], pstate[i+SFMT_POS1/2], r);
        pstate[i] = r;
    }
    for (; i < SFMT_N256; ++i) {
        r = mm256_recursion(pstate[i], pstate[i+(SFMT_POS1-SFMT_N)/2], r);
        pstate[i] = r;
    }
#if PRINT_TO_CHECK_AVX_CODE
    print_state(sfmt->state);
#endif
}

/**
 * This function fills the user-specified array with pseudorandom
 * integers.
 * @param sfmt SFMT internal state.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
static void gen_rand_array(sfmt_t * sfmt, w128_t * input, unsigned size128)
{
    int size = size128 / 2;
    unsigned i, j;
    __m256i r;
    __m256i * pstate = (__m256i*)sfmt->state;
    __m256i * array = (__m256i*)input;

    r = pstate[SFMT_N256 - 1];
    for (i = 0; i < (SFMT_N-SFMT_POS1)/2; ++i) {
        r = mm256_recursion(pstate[i], pstate[i + SFMT_POS1/2], r);
        pstate[i] = r;
    }
    for (i = (SFMT_N-SFMT_POS1)/2; i < SFMT_N256; ++i) {
        r = mm256_recursion(pstate[i], array[i + (SFMT_POS1-SFMT_N)/2], r);
        pstate[i] = r;
    }
    for (i = SFMT_N256; i < size - SFMT_N256; ++i) {
        r = mm256_recursion(array[i - SFMT_N/2], array[i + (SFMT_POS1-SFMT_N)/2], r);
        pstate[i] = r;
    }
    for (j = 0; j < SFMT_N - size; ++j) {
        pstate[j] = array[j + size - SFMT_N/2];
    }
    for (i = size - SFMT_N256; i < size; ++i, ++j) {
        r = mm256_recursion(array[i - SFMT_N/2], array[i + (SFMT_POS1-SFMT_N)/2], r);
        pstate[i] = r;
    }
}

#endif
