// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Started on Monday 5 June 2018, which was a very nice day in Strasbourg

#ifndef SIMD_H
#define SIMD_H

/// disable SIMD code by default:
#define USE_SIMD 0

// restrict functions to local scope
#define LOCAL static inline

#if defined(__SSE3__)
#  include <immintrin.h>
#  include "simd_sse.h"
#  undef USE_SIMD
#  define USE_SIMD 1
#endif

#if defined(__AVX__)
#  include <immintrin.h>
#  include "simd_avx.h"
#  undef USE_SIMD
#  define USE_SIMD 3
#endif

#if defined(__ARM_NEON__)
#  include <arm_neon.h>
#  include "simd_neon.h"
#  undef USE_SIMD
#  define USE_SIMD 8
#endif

#undef LOCAL

#endif // SIMD_H
