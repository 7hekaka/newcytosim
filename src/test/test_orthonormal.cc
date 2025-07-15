// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "real.h"
#include "timer.h"
#include "vector3.h"
#include "random.h"


void randomize(Vector3 & vec)
{
    vec.XX = RNG.sreal();
    vec.YY = RNG.sreal();
    vec.ZZ = RNG.sreal();
}


void test_orthonormal()
{
    Vector3 Z, X, Y;
    randomize(Z);
    
    real N = Z.norm();
    real NN = 4.0;
    Z.orthonormal(N, X, Y, 2.0/N);
    printf("  %7.3f : %7.3f %7.3f %7.3f :", N, X[0]/N, X[1]/N, X[2]/N);
    printf("  %7.3f %7.3f %7.3f ", dot(X,X)/NN, dot(X,Y)/NN, dot(X,Z)/NN);
    printf("  %7.3f %7.3f %7.3f ", dot(Y,X)/NN, dot(Y,Y)/NN, dot(Y,Z)/NN);
    printf("  %7.3f %7.3f %7.3f ", dot(Z,X)/NN, dot(Z,Y)/NN, dot(Z,Z)/NN);
    
    real A = RNG.sreal() * M_PI;
    real C = std::cos(A), S = std::sin(A);
    Vector3 V = Z.orthogonalNCS(N, C, S);
    printf("  %7.3f : %7.3f\n", N, norm(V));
}


int main(int argc, char * argv[])
{
    printf("test_orthonormal --- %lu bytes real --- %s\n", sizeof(real), __VERSION__);

    size_t CNT = 1<<10;
    if ( argc > 1 )
        CNT = strtol(argv[1], nullptr, 10);
    RNG.seed();
    
    for ( size_t i = 0; i < CNT; ++i )
        test_orthonormal();
}
