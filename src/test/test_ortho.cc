// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include <algorithm>
#include <cstdlib>


typedef double real;

/// absolute value
inline real abs_real(const real x) { return std::abs(x); }

/// return `neg` if `val < 0` and `pos` otherwise
inline float sign_select(float const val, float const neg, float const pos)
{
    return ( val < 0 ? neg : pos );
}

/// return `neg` if `val < 0` and `pos` otherwise
inline double sign_select(double const val, double const neg, double const pos)
{
    return ( val < 0 ? neg : pos );
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
    
    real XX, YY, ZZ;
    
    Vector3() { XX=0; YY=0; ZZ=0; }
    Vector3(real a, real b, real c) { XX=a; YY=b; ZZ=c; }
    void set(real a, real b, real c) { XX=a; YY=b; ZZ=c; }
    void set_rand() { XX=srand(); YY=srand(); ZZ=srand(); }
    
    real norm() { return sqrt(XX*XX + YY*YY + ZZ*ZZ); }

public:
    
    void set_unit(real N)
    {
        real nn;
        do {
            XX = srand();
            YY = srand();
            ZZ = srand();
            nn = XX*XX + YY*YY + ZZ*ZZ;
        } while ( nn > 1 );
        nn = N / sqrt(nn);
        XX *= nn;
        YY *= nn;
        ZZ *= nn;
    }

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


/**
 This works only if norm(z[]) == n!
 Derived from `Building an Orthonormal Basis, Revisited`,
 Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
 optimized by Marc B. Reynolds
*/
void orthonormal(const real z[3], real n, real x[3], real y[3])
{
    const real s = std::copysign(real(1), z[2]);
    const real a = z[1] / ( z[2] + s * n );
    const real b = z[1] * a;
    const real c = z[0] * a;
    x[0] = -z[2] - b;
    x[1] = c;
    x[2] = z[0];
    y[0] = s * c;
    y[1] = s * b - n;
    y[2] = s * z[1];
}

real dot(real A[3], real B[3]) { return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; }

void test_orthonormal()
{
    const real N = 2;
    const size_t MAX = 1 << 14;
    Vector3 vec[MAX];
    
    for ( size_t i = 0; i < MAX; ++i )
    {
        vec[i].set_unit(N);
        real Z[3] = { vec[i].XX, vec[i].YY, vec[i].ZZ };
        real Y[3] = { 0 };
        real X[3] = { 0 };

        orthonormal(Z, N, X, Y);
        printf("  %9.3f %9.3f %9.3f :", X[0], X[1], X[2]);
        printf("  %9.3f %9.3f %9.3f ", dot(X, X), dot(X,Y), dot(X,Z));
        printf("  %9.3f %9.3f %9.3f ", dot(Y, X), dot(Y,Y), dot(Y,Z));
        printf("  %9.3f %9.3f %9.3f\n", dot(Z, X), dot(Z,Y), dot(Z,Z));
    }
}

void test(size_t cnt)
{
    const size_t MAX = 1 << 14;
    Vector3 vec[MAX];
    Vector3 vic[MAX];
    
    for ( size_t i = 0; i < MAX; ++i )
        vec[i].set_rand();
    
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
    
#if 1
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        for ( size_t i = 0; i < MAX; ++i )
            vic[i] = vec[i].orthogonalB();
        for ( size_t i = 0; i < MAX; ++i )
            vec[i] = vic[i].orthogonalB();
    }
    printf("  orthogonalB   %7.3f\n", toc(cnt*MAX));
#endif
}


int main()
{
    printf("test_ortho --- %lu bytes real --- %s\n", sizeof(real), __VERSION__);
    //test(1<<15);
    test_orthonormal();
}
