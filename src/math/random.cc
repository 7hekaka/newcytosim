// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "random.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cstring>
#include <ctime>
#include <random>

/// static object
thread_local Random RNG;



/// the most significant bit in a 32-bits integer
[[maybe_unused]] constexpr uint32_t BIT31 = 1U << 31;

/// sign bit in double precision (double)
[[maybe_unused]] constexpr uint64_t BIT63 = 1ULL << 63;


/**
 The generator is initialized with a zero state vector,
 and seed() must be called before any random number can be produced.
 */
Random::Random()
{
    uintptr_t a = ((uintptr_t)this & 31);
    if ( a ) fprintf(stderr, "Random missaligned on 32 + %lu\n", a);

    //fprintf(stderr, "Random with SFMT_N32 = %i\n", SFMT_N32);
    
    // clear state (not necessary):
    memset(integers_, 0, 4*SFMT_N32);
    memset(gaussians_, 0, sizeof(real)*SFMT_N32);
    memset(twister_.state[0].u, 0, 4*SFMT_N32);

    // initialize pointers signalling an empty reserve:
    start_ = integers_;
    finish_ = start_;
    next_gaussian_ = gaussians_;
    next_exponential_ = exponentials_;
}


Random::~Random()
{
    //printf("Random Number Generator released\n");
}


void Random::seed(const uint32_t s)
{
    //printf("Random:%p:seed(%u)\n", this, s);
    sfmt_init_gen_rand(&twister_, s);
    refill();
}


void Random::seed()
{
    uint32_t s = 7;
    try {
        // read system source if available:
        std::random_device rd;
        s = rd();
        //std::cerr << "random device ---> seed " << s << '\n';
    }
    catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
    }
    seed(s);
}


bool Random::seeded()
{
    uint32_t * buf = twister_.state[0].u;
    for ( size_t n = 0; n < SFMT_N32; ++n )
        if ( buf[n] )
            return true;
    return false;
}


//------------------------------------------------------------------------------
#pragma mark - Vectors

void Random::add_srand1(real dst[1], const real ptr[1], real mag)
{
    mag *= TWO_POWER_MINUS_31;
    if ( finish_ <= start_ )
        refill();
    --finish_;
    dst[0] = ptr[0] + mag * static_cast<real>(finish_[0]);
}

void Random::add_srand2(real dst[2], const real ptr[2], real mag)
{
    mag *= TWO_POWER_MINUS_31;
    if ( finish_ < 2 + start_ )
        refill();
    finish_ -= 2;
    dst[0] = ptr[0] + mag * static_cast<real>(finish_[0]);
    dst[1] = ptr[1] + mag * static_cast<real>(finish_[1]);
}

void Random::add_srand3(real dst[3], const real ptr[3], real mag)
{
    mag *= TWO_POWER_MINUS_31;
    if ( finish_ < 3 + start_ )
        refill();
    finish_ -= 3;
    dst[0] = ptr[0] + mag * static_cast<real>(finish_[0]);
    dst[1] = ptr[1] + mag * static_cast<real>(finish_[1]);
    dst[2] = ptr[2] + mag * static_cast<real>(finish_[2]);
}

void Random::sreal2(real& a, real& b)
{
    if ( finish_ < 2 + start_ )
        refill();
    finish_ -= 2;
    int32_t const * ptr = reinterpret_cast<int32_t const *>(finish_);
    a = TWO_POWER_MINUS_31 * static_cast<real>(ptr[0]);
    b = TWO_POWER_MINUS_31 * static_cast<real>(ptr[1]);
}

void Random::sreal4(real& a, real& b, real& c, real& d)
{
    if ( finish_ < 4 + start_ )
        refill();
    finish_ -= 4;
    int32_t const * ptr = reinterpret_cast<int32_t const *>(finish_);
    a = TWO_POWER_MINUS_31 * static_cast<real>(ptr[0]);
    b = TWO_POWER_MINUS_31 * static_cast<real>(ptr[1]);
    c = TWO_POWER_MINUS_31 * static_cast<real>(ptr[2]);
    d = TWO_POWER_MINUS_31 * static_cast<real>(ptr[3]);
}

/** Set two signed random numbers in [-1, 1] such that `C*C + S*S == 1` */
void Random::urand2(real& C, real& S)
{
    real d, x, y;
    do {
        sreal2(x, y);
        d = x*x + y*y;
    } while (( d > 1.0 )|( d == 0 ));
    // we can avoid the square root here by using the half angle formula:
    // consider the imaginary number a = x + i*y = sqrt(d) * exp(i*angle), then
    // a*a = d * exp(i*2*angle), thus a*a/d is uniformly distributed in the unit disc:
    C = ( x*x - y*y ) / d;
    S = ( x*y + x*y ) / d;
}

/** Set two signed random numbers in [-R, R] such that `C*C + S*S == R*R` */
void Random::urand2(real& C, real& S, const real R)
{
    real d, x, y;
    do {
        sreal2(x, y);
        d = x*x + y*y;
    } while (( d > 1.0 )|( d == 0 ));
    d = R / d;
    C = ( x*x - y*y ) * d;
    S = ( x*y + x*y ) * d;
}

//------------------------------------------------------------------------------
#pragma mark - Gaussian derivates


/**
 Set two signed real number, following a normal law N(0,v*v)
 using Box-Muller method (George E. P. Box et Mervin E. Muller, 1958)
 */
void Random::box_muller(real& x, real& y)
{
    real w = std::sqrt( -2 * std::log(preal()) );
    real a = M_PI * sreal();
    x = w * std::cos(a);
    y = w * std::sin(a);
}

/**
 Set two signed real number, following a normal law N(0,v*v)
 using Marsaglia polar method (George Marsaglia, 1964)
 */
void Random::gauss_set(real& a, real& b, real v)
{
    real x, y, w;
    do {
        x = sreal();
        y = sreal();
        w = x * x + y * y;
    } while ( w >= 1.0 || w == 0 );
    /*
     formula below are only valid if ( w > 0 ),
     which may be false only with a minuscule probability
     */
    w = v * std::sqrt( -2 * std::log(w) / w );
    a = w * x;
    b = w * y;
}


/**
 Fill array `vec[]` with Gaussian values ~ N(0,1).
 the size of `vec` should be a multiple of 2, and sufficient to hold `end-src` values
 For each 4 input values, this produces ~PI values.
 @Return address past the last value stored in `dst`
 */
real * makeGaussians(real dst[], size_t cnt, const int32_t src[])
{
    const real alpha(TWO_POWER_MINUS_31);
    int32_t const*const end = src + cnt;
    while ( src < end )
    {
        real x = real(src[0]) * alpha;
        real y = real(src[1]) * alpha;
        src += 2;
#if 1
        /**
        This folds the corners of [-1, 1] x [-1, 1] that are outside the unit circle,
        to map these points back onto the original square. For randomly distributed points,
        this increases the number of points within the unit circle by a factor 3-2*sqrt(2)
        without changing the property of being equidistributed within the unit circle.
        */
        if ( abs_real(x) + abs_real(y) >= M_SQRT2 )
        {
            constexpr real S = M_SQRT1_2 + 1;
            // subtract corner and scale to recover a square of size sqrt(1/2)
            real cx = S * x - std::copysign(S, x);
            real cy = S * y - std::copysign(S, y);
            // apply rotation, scaling by sqrt(2): x' = y + x;  y' = y - x
            x = cy + cx;
            y = cy - cx;
        }
#endif
        real w = x * x + y * y;
        if (( 0 < w ) & ( w <= 1 ))
        {
            w = std::sqrt( -2 * std::log(w) / w );
            dst[0] = w * x;
            dst[1] = w * y;
            dst += 2;
        }
    }
    return dst;
}

/**
 Fill array `gaussians_` with approximately 500 Gaussian values ~ N(0,1).
 Set `next_gaussian` past the last position containing a valid number.
 The number of gaussian values set by this function is random,
 and it may even be zero.
 */
void Random::refill_gaussians()
{
    next_gaussian_ = makeGaussians(gaussians_, SFMT_N32, (int32_t*)twister_.state);
    sfmt_gen_rand_all(&twister_);
}


#if ( 0 )

/**
 Fill `n` Gaussian values ~ N(0,1) in array `vec[]`.
 */
void Random::gauss_set(real vec[], size_t cnt, real v = 1.0)
{
    unsigned u = cnt % 8;
    unsigned w = u % 2;
    
    if ( w )
        vec[0] = v * gauss();
    
    for ( ; w < u; w += 2 )
        gauss_set(vec[w], vec[w+1], v);
    
    for ( ; u < cnt; u += 8 )
    {
        gauss_set(vec[u  ], vec[u+1], v);
        gauss_set(vec[u+2], vec[u+3], v);
        gauss_set(vec[u+4], vec[u+5], v);
        gauss_set(vec[u+6], vec[u+7], v);
    }
}

#else

/**
 Fill `n` Gaussian values ~ N(0,1) in array `vec[]`.
 */
void Random::gauss_set(real vec[], size_t cnt)
{
    size_t n = (size_t)( next_gaussian_ - gaussians_ );
    // check if `vec` would consume all the buffer:
    while ( n <= cnt )
    {
        // use all values in buffer:
        copy_real(n, gaussians_, vec);
        vec += n;
        cnt -= n;
        refill_gaussians();
        n = (size_t)( next_gaussian_ - gaussians_ );
    };
    
    // use `cnt` values from buffer:
    next_gaussian_ -= cnt;
    copy_real(cnt, next_gaussian_, vec);
}

#endif


/**
 This is the Box & Muller method.
 A note on the generation of random normal deviates
 Box & Muller, 1958
    n = sqrt( -2 * log(R) ), where R is random in [0, 1]
 and with A = random in [-PI, PI]:
    x = n * cos(A)
    y = n * sin(A),
 */
void Random::gauss_boxmuller(real& x, real& y)
{
    real ang = real(SRAND32()) * ( TWO_POWER_MINUS_31 * M_PI );
    real nrm = std::sqrt( -2 * std::log( preal_exc() ));
    x = nrm * std::cos(ang);
    y = nrm * std::sin(ang);
}


//------------------------------------------------------------------------------
#pragma mark - Exponential derivates

/// fill array `dst` with random numbers following the Normal law.
/**
 Using idependent random numbers provided by 'src[]', or size `cnt`.
 The Normal law has Mean = 0 and Variance = 1.
 @return the address pass the last value set.
 The number of values is `return - dst` which should be 2*cnt/3.
 
 Using method from Luc Devroye: Non-uniform random variate generation. Page 394:
 From 3 independent identically distributed uniform [0,1] random variates U, V, W.
     Y = -log(U*V)
     RETURN W*Y, (l-W)*Y
 */
float * makeExponentials(float dst[], size_t cnt, const uint32_t src[])
{
    cnt -= 2;
    const float alpha(TWO_POWER_MINUS_32);
    for ( size_t i = 0; i < cnt; i += 3 )
    {
        float U = alpha * static_cast<float>(src[i+0]) + alpha;
        float V = alpha * static_cast<float>(src[i+1]) + alpha;
        float W = alpha * static_cast<float>(src[i+2]);
        float y = -logf(U*V);
        dst[0] = y * W;
        dst[1] = y * ( 1.f - W );
        dst += 2;
    }
    return dst;
}


void Random::refill_exponentials()
{
    float * ptr = makeExponentials(exponentials_, SFMT_N32, (uint32_t*)twister_.state);
    //assert_true( ptr == exponentials_ + SFMT_N32 );
    next_exponential_ = ptr;
    //printf("refill_exponentials\n");
    sfmt_gen_rand_all(&twister_);
}


//------------------------------------------------------------------------------
#pragma mark - Integers

/**
 integer in [0,n] for n < 2^32
 */
uint32_t Random::pint32_slow(const uint32_t n)
{
    // Find which bits are used in n
    uint32_t used = n | ( n >> 1 );
    used |= (used >> 2);
    used |= (used >> 4);
    used |= (used >> 8);
    used |= (used >> 16);
    
    // Draw numbers until one is found in [0,n]
    uint32_t i;
    do
        i = URAND32() & used;  // toss unused bits to shorten search
    while ( i > n );
    return i;
}


/**
 returns an integer in [0 n], with the probabilities given as arguments
 The sum of `ratio[]` is calculated
 */
uint32_t Random::pint32_ratio(const uint32_t n, const uint32_t ratio[])
{
    uint32_t ii = 0, sum = 0;
    for ( ii = 0; ii < n; ++ii )
        sum += ratio[ii];
    sum = (uint32_t) std::floor( preal() * sum );
    while ( sum >= ratio[ii] )
        sum -= ratio[ii++];
    return ii;
}


/**
 Return Poisson distributed integer, with expectation=E and variance=E
 http://en.wikipedia.org/wiki/Poisson_distribution
 
 This routine is slow for large values of E.
 If E > 256, this returns a Gaussian distribution of parameter (E, E),
 which is a good approximation of the Poisson distribution
 
 Knuth D.E. The art of computer programming, Vol II: Seminumerical algorithms.
 
 This method fails for E > 700, in double precision
 */
uint32_t Random::poisson_knuth(const real E)
{
    if ( E > 256 )
        return static_cast<uint32_t>( gauss() * std::sqrt(E) + E );
    if ( E < 0 )
        return 0;
    real L = std::exp(-E);
    real p = preal();
    uint32_t k = 0;
    while ( p > L )
    {
        ++k;
        p *= preal();
    }
    return k;
}


/**
 Return Poisson distributed integer, with expectation=E  variance=E
 http://en.wikipedia.org/wiki/Poisson_distribution
 
 This routine is slow for large values of E.
 If E > 512, this returs a Gaussian distribution of parameter (E, E),
 which is a good approximation of the Poisson distribution.
 
 This method fails for E > 700, in double precision
 */
uint32_t Random::poisson(const double E)
{
    if ( E > 256 )
        return static_cast<uint32_t>( gauss() * std::sqrt(E) + E );
    if ( E < 0 )
        return 0;
    double p = std::exp(-E);
    double s = p;
    uint32_t k = 0;
    double u = preal();
    while ( u > s )
    {
        ++k;
        p *= E / k;
        s += p;
    }
    return k;
}


/**
 This is equivalent to calling poisson(exp(-E))
 The argument is EL = exp(-E)
 expectation=E  variance=E (see wikipedia, Poisson Distribution)
 */
uint32_t Random::poissonE(const double EL)
{
    double p = preal();
    uint32_t k = 0;
    while ( p > EL )
    {
        ++k;
        p *= preal();
    }
    return k;
}


uint32_t Random::geometric(const real P)
{
    if ( P < 0 )
        return 0;
    const uint32_t pi = (uint32_t)( P * 0x1p32 );
    
    uint32_t s = 0;
    while ( URAND32() > pi )
        ++s;
    return s;
}


uint32_t Random::binomial(const int N, const real P)
{
    if ( P < 0 )
        return 0;
    const uint32_t pi = (uint32_t)( P * 0x1p32 );
    
    uint32_t s = 0;
    for ( int x = 0; x < N; ++x )
        if ( URAND32() < pi )
            ++s;
    return s;
}

