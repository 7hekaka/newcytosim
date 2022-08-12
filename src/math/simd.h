// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// Started on Monday 5 June 2018, which was a very nice day in Strasbourg

#ifndef SIMD_H
#define SIMD_H

// restrict functions to local scope
#define LOCAL static inline

#if defined(__SSE3__)
#include <immintrin.h>
#include "simd_sse.h"
#endif

#if defined(__AVX__)
#include <immintrin.h>
#include "simd_avx.h"
#endif

#if defined(__ARM_NEON__)
#include <arm_neon.h>
//#include "simd_neon.h"
#endif

#undef LOCAL

#endif // SIMD_H
