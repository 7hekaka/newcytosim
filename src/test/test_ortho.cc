// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include <algorithm>
#include <cstdlib>


typedef float real;

/// absolute value
inline real abs_real(const real x) { return std::abs(x); }

/// return `neg` if `arg < 0` and `pos` otherwise
inline real sign_select(real const& val, real const& neg, real const& pos)
{
    // this should be branchless, using Intel's VBLENDVPD instruction
    if ( val > 0 ) return pos;
    else return neg;
}


/// signed random real in [-1, 1]
real srand()
{
    int i = random() << 1;
    return static_cast<real>(i) * 0x1p-31;
}


/// keeping time using Intel's cycle counters
unsigned long long rdt = 0;
/// start timer
inline void tic() { rdt = __rdtsc(); }
/// stop timer and print time
inline double toc(double num) { return double(__rdtsc()-rdt)/num; }


class Vector3
{
public:
    
    real XX, YY, ZZ, TT;
    
    Vector3(real a, real b, real c) { XX=a; YY=b; ZZ=c; }
    void set(real a, real b, real c) { XX=a; YY=b; ZZ=c; }
    Vector3() { XX=srand(); YY=srand(); ZZ=srand(); }

public:
    
    Vector3 orthogonal() const
    {
        real ax = abs_real(XX);
        real ay = abs_real(YY);
        real az = abs_real(ZZ);
        
        real x = ax - std::min(ay, az); // negative if XX is smallest
        real y = ay - std::min(az, ax);
        real z = std::min(x, y); // use z if x and y are positive (could use bitwise OR)
        
        return Vector3(sign_select(z, 0, YY) - sign_select(y, ZZ, 0),
                       sign_select(x, ZZ, 0) - sign_select(z, 0, XX),
                       sign_select(y, XX, 0) - sign_select(x, YY, 0));
    }
    
    Vector3 orthogonalB() const
    {
        if ( abs_real(XX) < abs_real(YY) )
        {
            if ( abs_real(XX) < abs_real(ZZ) )
                return Vector3( 0, -ZZ, YY); //XX is the smallest
            else
                return Vector3(YY, -XX, 0); //ZZ is the smallest
        }
        else
        {
            if ( abs_real(YY) < abs_real(ZZ) )
                return Vector3(-ZZ,   0,  XX); //YY is the smallest
            else
                return Vector3( YY, -XX,   0); //ZZ is the smallest
        }
    }
};


void test(size_t cnt)
{
    const size_t MAX = 1 << 14;
    Vector3 vec[MAX];
    Vector3 vic[MAX];
    /*
    for ( size_t i = 0; i < MAX; ++i )
        printf("%9.3f %9.3f %9.3f\n", vec[i].XX, vec[i].YY, vec[i].ZZ);
     */
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        for ( size_t i = 0; i < MAX; ++i )
            vic[i] = vec[i].orthogonal();
        for ( size_t i = 0; i < MAX; ++i )
            vec[i] = vic[i].orthogonal();
    }
    printf("  orthogonal    %7.3f\n", toc(cnt*MAX));
    
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        for ( size_t i = 0; i < MAX; ++i )
            vic[i] = vec[i].orthogonalB();
        for ( size_t i = 0; i < MAX; ++i )
            vec[i] = vic[i].orthogonalB();
    }
    printf("  orthogonalB   %7.3f\n", toc(cnt*MAX));
}


int main ()
{
    test(1<<15);
    return 0;
}
