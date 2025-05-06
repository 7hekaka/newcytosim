// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef RANDOM_H
#define RANDOM_H

#include <cstdint>

#ifndef REAL_H
#  include "real.h"
#endif

#define SFMT_MEXP 19937

#include "SFMT.h"

/// the maximum value of a signed 32-bit integer is 2^31-1
#define TWO_POWER_MINUS_31 0x1p-31
/// the maximum value of a unsigned 32-bit integer is 2^32-1
#define TWO_POWER_MINUS_32 0x1p-32
/// the maximum value of a unsigned 64-bit integer is 2^64-1
#define TWO_POWER_MINUS_64 0x1p-64


/// Random Number Generator
/**
 The generation of random bits is done with Mersenne Twister from U. of Hiroshima
 
 http://en.wikipedia.org/wiki/Mersenne_twister
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 
 This class provides a convenient interface, and functions to generate floating
 point values following various distributions (Gaussian, Exponential, Poisson),
 and other facilities.
 
 F. Nedelec, reclassified on 13.09.2018
*/
class alignas(32) Random
{
    /// reserve of random integers
    uint32_t integers_[SFMT_N32];

    /// reserve of normally-distributed numbers, with zero mean and unit variance
    real gaussians_[SFMT_N32];
    
    /// reserve of exponentially-distributed numbers, with unit mean
    float exponentials_[SFMT_N32];

    /// Mersenne Twister Generator
    sfmt_t twister_;

    /// pointer to integer reserve, at unused lower end
    uint32_t const* start_;

    /// pointer to integer reserve, at unused upper end + 1
    uint32_t const* finish_;
    
    /// pointer to access the next value in `gaussians_[]`
    real * next_gaussian_;
    
    /// pointer to access the next value in `exponentials_[]`
    float * next_exponential_;
    
protected:
    
    /// replenish state vector
    void refill()
    {
        memcpy(integers_, twister_.state, 4*SFMT_N32);
        start_ = integers_;
        finish_ = start_ + SFMT_N32;
        //printf("Random:%p:refill %p\n", this, pthread_self());
        sfmt_gen_rand_all(&twister_);
    }

    /// extract next 32 random bits
    uint32_t URAND32()
    {
        if ( finish_ <= start_ )
            refill();
        --finish_;
        return *finish_;
    }

    /// extract next 32 random bits as signed integer
    int32_t SRAND32()
    {
        if ( finish_ <= start_ )
            refill();
        --finish_;
        return *reinterpret_cast<int32_t const*>(finish_);
    }
    
    /// extract next 64 random bits
    uint64_t URAND64()
    {
        // consume data from the 'start_' to preserve alignment
        if ( finish_ < 2 + start_ )
            refill();
        uint64_t const* p = reinterpret_cast<uint64_t const*>(start_);
        start_ += 2;
        return *p;
    }
    
    /// extract next 64 random bits as signed integer
    int64_t SRAND64()
    {
        // consume data from the 'start_' to preserve alignment
        if ( finish_ < 2 + start_ )
            refill();
        int64_t const* p = reinterpret_cast<int64_t const*>(start_);
        start_ += 2;
        return *p;
    }
    
    /// a non-negative 'real' using 31 random bits, in [0, 1[
    real ZERO2ONE()
    {
        // the cast gives a number with (x <= 2^31-1)
        return std::fabs(static_cast<real>(SRAND32())) * TWO_POWER_MINUS_31;
    }

    /// a signed 'real' using 31 random bits + sign bit (-2^31+1 < x <= 2^31-1)
    real BIGREAL()
    {
        return static_cast<real>(SRAND32());
    }

public:
    
    /// Constructor sets the state vector to zero
    Random();
    
    /// destructor
    ~Random();

    /// true if state vector is not entirely zero
    bool seeded();

    /// seed with given 32 bit integer
    void seed(const uint32_t s);
    
    /// seed from std::random_device
    void seed();

    /// signed integer in [-2^31+1, 2^31-1];
    int32_t  sint32() { return SRAND32(); }

    /// unsigned integer in [0, 2^32-1]
    uint32_t pint32() { return URAND32(); }
    
    /// unsigned integer in [0, 2^64-1]
    int64_t  sint64() { return SRAND64(); }

    /// unsigned integer in [0, 2^64-1]
    uint64_t pint64() { return URAND64(); }

#if ( 0 )
    /// unsigned integer in [0,n-1] for n < 2^32
    uint32_t pint32(const uint32_t& n) { return uint32_t(ZERO2ONE()*n); }

    /// unsigned integer in [0,n-1] for n < 2^64
    uint64_t pint64(const uint64_t& n) { return uint64_t(ZERO2ONE()*n); }
#else
    /// unsigned integer in [0,n-1] for n < 2^32, Daniel Lemire's method
    /** Fast Random Integer Generation in an Interval, 2019 https://doi.org/10.1145/3230636 */
    uint32_t pint32(const uint32_t& n) { return (uint32_t)(((uint64_t)URAND32() * (uint64_t)n) >> 32); }
    
    /// unsigned integer in [0,n-1] for n < 2^64, Daniel Lemire's method
    uint64_t pint64(const uint64_t& p) {
#ifdef __SIZEOF_INT128__ // then we know we have 128-bit integers
        return (uint64_t)(((__uint128_t)URAND64() * (__uint128_t)p) >> 64);
#else
        return URAND64() % p; // fallback
#endif
    }
#endif
 
    /// unsigned integer in [0,n-1] for n < 2^32, Daniel Lemire's fair method
    uint32_t pint32_fair(const uint32_t& range)
    {
        uint64_t multiresult = (uint64_t)URAND32() * (uint64_t)range;
        uint32_t leftover = (uint32_t) multiresult;
        if ( leftover < range )
        {
            uint32_t threshold = -range % range;
            while ( leftover < threshold )
            {
                multiresult = (uint64_t)URAND32() * (uint64_t)range;
                leftover = (uint32_t) multiresult;
            }
        }
        return (uint32_t)(multiresult >> 32);
    }

    /// integer in [0,n] for n < 2^32, (slow) bitwise algorithm
    uint32_t pint32_slow(uint32_t n);
    
    /// integer in [a,b] for a, b < 2^32
    uint32_t pint32(uint32_t a, uint32_t b) { return a + pint32(1+b-a); }

    /// integer in [0 N], with probabilities given in ratio[] of size N, with sum(ratio)>0
    uint32_t pint32_ratio(uint32_t n, const uint32_t ratio[]);

    /// integer k of probability distribution p(k,E) = exp(-E) * pow(E,k) / factorial(k)
    uint32_t poisson(double E);
    
    /// integer k of probability distribution p(k,E) = EL * pow(E,k) / factorial(k)
    uint32_t poissonE(double EL);
    
    /// integer k of probability distribution p(k,E) = exp(-E) * pow(E,k) / factorial(k)
    uint32_t poisson_knuth(real E);

    /// number of successive unsuccessful trials, when success has probability p (result >= 0)
    uint32_t geometric(real p);

    /// number of sucesses among n trials of probability p
    uint32_t binomial(int n, real p);
    
    
    /// returns true with probability (p), and false with probability (1-p)
    bool test(real p)     { return ( ZERO2ONE() <  p ); }
    
    /// returns true with probability (1-p), and false with probability (p)
    bool test_not(real p) { return ( ZERO2ONE() >= p ); }
    
    /// 0  or  1  with equal chance
    int flip()            { return URAND32() >> 31; }
    
    /// returns -1  or  1 with equal chance
    int flipsign()        { return (SRAND32() > 0) ? 1 : -1; }
    
    /// returns 1 with probability P and -1 with probability 1-P
    int flipsign(real p)  { return 2*(int)test(p) - 1; }

    /// True with probability 1/8
    bool flip_8th()       { return URAND32() < 1<<29; }

    /// random float in [0,1[, requires IEEE Standard 754 
    float pfloat()        { return float(URAND32() >> 8) * float(0x1.fp-24); }
    
    /// random float in ]-1,1[, requires IEEE Standard 754
    float sfloat()        { return float(SRAND32() >> 8) * float(0x1.fp-23); }
    
    /// slow random double in [0,1[, using two uint32_t to set all the fraction bits, requires IEEE Standard 754
    double pdouble()       { return double(URAND64() >> 11) * 0x1.0p-53; }
    
    /// slow random double in ]-1,1[, using two uint32_t to set all the fraction bits, requires IEEE Standard 754
    double sdouble()       { return double(SRAND64() >> 11) * 0x1.0p-52; }
    
    /// positive real number in [0,1[, zero included
    real preal()           { return ZERO2ONE(); }
    
    /// positive real number in [0,n[ = n * preal() : deprecated, use preal() * n
    real preal(real n)     { return n * ( ZERO2ONE() ); }
    
    /// signed real number in ]-1,1[, boundaries excluded
    real sreal()           { return BIGREAL() * TWO_POWER_MINUS_31; }
    
    /// signed real number in ]-1/2, 1/2[, boundaries excluded
    real shalf()           { return BIGREAL() * TWO_POWER_MINUS_32; }
    
    /// returns -1.0 or 1.0 with equal chance
    real sflip()           { return (SRAND32() > 0) ? real(1.0) : real(-1.0); }
    
    /// returns -a or a with equal chance
    real sflip(real a)     { return (SRAND32() > 0) ? a : -a; }

    /// non-zero real number in ]0,1]
    real preal_exc()       { return 1 - ZERO2ONE(); }
    
    /// non-zero real number in ]0,n]
    real preal_exc(real n) { return preal_exc() * n; }
    
    /// set 2 random number in [-1, 1]
    void sreal2(real&, real&);
    
    /// set 4 random number in [-1, 1]
    void sreal4(real&, real&, real&, real&);

    /// set 2 random numbers such that C*C + S*S = 1
    void urand2(real&, real&);
 
    /// set 2 random numbers such that C*C + S*S = R
    void urand2(real&, real&, real R);

    /// real number uniformly distributed in [a,b[
    real real_uniform(real a, real b) { return a + preal() * ( b - a ); }
    
    /// add random number in ]-mag, mag[ to each component of 2D vec[]
    void add_srand1(real dst[1], const real vec[1], real mag);
    
    /// add random number in ]-mag, mag[ to each component of 2D vec[]
    void add_srand2(real dst[2], const real vec[2], real mag);
    
    /// add random number in ]-mag, mag[ to each component of 3D vec[]
    void add_srand3(real dst[3], const real vec[3], real mag);
    
    /// refill array `gaussians_[]`, resetting `next_gaussian_`
    void refill_gaussians();
  
    /// set two independent random numbers, both following a normal law N(0,1)
    void box_muller(real&, real&);
    
    /// set two independent random numbers, both following a normal law N(0,v*v)
    void gauss_set(real&, real&, real v);

    /// random Gaussian number, following a normal law N(0,1)
    real gauss()
    {
        while ( next_gaussian_ <= gaussians_ )
            refill_gaussians();
        --next_gaussian_;
        return *next_gaussian_;
    }

    /// fill array `vec` with independent random numbers following normal law N(0,1).
    void gauss_set(real vec[], size_t n);
    
    /// fill array `vec` with independent random numbers following normal law N(0,v*v).
    void gauss_set(real vec[], size_t n, real v);

    /// set 2 statistically independent signed real number, following the normal law N(0,1), slower algorithm
    void gauss_boxmuller(real &, real&);
    
    /// refill array `exponentials_[]`, resetting `next_exponential_`
    void refill_exponentials();

    /// random in [0, inf[, distributed as P(x>m) = exp(-m); mean = 1.0, variance = 1.0
    float exponential()
    {
        //return -std::log(1 - ZERO2ONE());
        if ( next_exponential_ <= exponentials_ )
            refill_exponentials();
        --next_exponential_;
        return *next_exponential_;
    }

    /// exponentially distributed positive real, with P(x) = exp(-x/E) / E,  parameter E is 1/Rate
    real exponential(const real E) { return -E * std::log(1 - ZERO2ONE());  }

    /// fair choice among two given values
    template<typename T>
    T choice(const T& x, const T& y)
    {
        if ( SRAND32() > 0 )
            return x;
        else
            return y;
    }
    
    /// fair choice within an array of `size` values
    template<typename T>
    T choice(const T val[], uint32_t size)
    {
        return val[ pint32(size) ];
    }
    
    /// uniform shuffling of array `T[]`.
    /** Algorithm from knuth's The Art of Programming, Vol 2 chp. 3.4.2 */
    template <typename T> 
    void shuffle(T val[], uint32_t size)
    {
        uint32_t jj = size, kk;
        while ( jj > 1 )
        {
            kk = pint32(jj);
            --jj;
            T tmp   = val[jj];
            val[jj] = val[kk];
            val[kk] = tmp;
        }
    }
};


/// The Random Number Generator is thread local to avoid data corruption
extern thread_local Random RNG;


#endif  //RANDOM_H
