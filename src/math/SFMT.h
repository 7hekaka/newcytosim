#pragma once
/**
 * @file SFMT.h
 *
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) pseudorandom
 * number generator using C structure.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2006, 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University.
 * Copyright (C) 2012 Mutsuo Saito, Makoto Matsumoto, Hiroshima
 * University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t
 #define PRIu64 "llu"
 #define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#ifndef SFMTST_H
#define SFMTST_H

/* F. Nedelec */
#ifdef __SSE2__
#define HAVE_SSE2 1
#endif
#if defined(__AVX2__)
#define HAVE_AVX2 1
#endif
#if defined(__ARM_NEON__)
#define HAVE_NEON 1
#endif
/* F. Nedelec */


#if defined(__cplusplus)
extern "C" {
#endif

#include <assert.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
  #include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
  typedef unsigned int uint32_t;
  typedef unsigned __int64 uint64_t;
  #define inline __inline
#else
  #include <inttypes.h>
  #if defined(__GNUC__)
    #define inline __inline__
  #endif
#endif

#ifndef PRIu64
  #if defined(_MSC_VER) || defined(__BORLANDC__)
    #define PRIu64 "I64u"
    #define PRIx64 "I64x"
  #else
    #define PRIu64 "llu"
    #define PRIx64 "llx"
  #endif
#endif

#include "SFMT-params.h"

/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
  #if !defined(__APPLE__)
    #include <altivec.h>
  #endif
/** 128-bit data structure */
union W128_T {
    vector unsigned int s;
    uint32_t u[4];
    uint64_t u64[2];
};
#elif defined(HAVE_NEON)
#include <arm_neon.h>
    
/** 128-bit data structure */
union W128_T {
    uint32_t u[4];
    uint64_t u64[2];
    uint32x4_t si;
};
#elif defined(HAVE_SSE2)
  #include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
    uint32_t u[4];
    uint64_t u64[2];
     __m128i xmm;
};
#else
/** 128-bit data structure */
union W128_T {
    uint32_t u[4];
    uint64_t u64[2];
};
#endif

#ifdef HAVE_AVX2
  #include <immintrin.h>
/** 256-bit data structure */
union W256_T {
    uint32_t u[8];
    uint64_t u64[4];
     __m128i xmm[2];
     __m256i ymm;
};
typedef union W256_T w256_t;
#endif

/** 128-bit data type */
typedef union W128_T w128_t;

/**
 * SFMT internal state
 */
struct SFMT_T {
    /** the 128-bit internal state array */
    w128_t state[SFMT_N];
    /** index counter to the 32-bit internal state array */
    int idx;
};

typedef struct SFMT_T sfmt_t;

void sfmt_fill_array32(sfmt_t * sfmt, uint32_t * array, int size);
void sfmt_fill_array64(sfmt_t * sfmt, uint64_t * array, int size);
void sfmt_init_gen_rand(sfmt_t * sfmt, uint32_t seed);
void sfmt_init_by_array(sfmt_t * sfmt, uint32_t * init_key, int key_length);
const char * sfmt_get_idstring(sfmt_t * sfmt);
int sfmt_get_min_array_size32(sfmt_t * sfmt);
int sfmt_get_min_array_size64(sfmt_t * sfmt);
void sfmt_gen_rand_all(sfmt_t * sfmt);

#if defined(__cplusplus)
}
#endif

#endif
