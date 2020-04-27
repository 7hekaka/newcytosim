// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "random.h"

#include "vector2.h"
#include "vector3.h"
#include "vector4.h"

#include "matrix33.h"
#include "matrix34.h"
#include "matrix44.h"
#include "matfull.h"


Vector3 dir(0, 1, 0);
Vector3 off(0.1, 0.2, -0.3);
Vector3 sun(0, 0, 0);

real diff(size_t size, real const* a, real const* b)
{
    real s = 0;
    for ( size_t i=0; i<size; ++i )
        s += abs_real( a[i] - b[i] );
    return s;
}


template <typename MATRIX>
void checkMatrix(MATRIX & mat)
{
    mat = MATRIX::outerProduct(dir, off);
    Vector3 vec = mat.vecmul(sun);
    Vector3 vik = mat.trans_vecmul(sun);

    std::clog << vec << " | " << vik << "  ";
    std::clog << mat << "  " << mat.transposed() << "\n";

    //printf("  check %+16.6f  %+16.6f ", diff);
}


void checkMatrixFull(Matrix44 const& src)
{
    MatrixFull mat;
    mat.resize(3);
    
    for ( size_t i = 0; i < 3; ++i )
    for ( size_t j = 0; j < 3; ++j )
        mat(i,j) = src(i,j);
    
    Vector3 vec, vik(0,0,0), vok;
    mat.vecMul(sun, vec);
    mat.vecMulAdd(sun, vik);
    mat.vecMulAVX(sun, vok);

    std::clog << vec << " | " << vik << " | " << vok << "\n";
    std::clog << mat << "\n";

    //printf("  check %+16.6f  %+16.6f ", diff);
}

//------------------------------------------------------------------------------
#pragma mark -

int main( int argc, char* argv[] )
{
    printf("test_block --- %s\n", __VERSION__);
    RNG.seed();
    
    Matrix33 M33(1,0);
    Matrix34 M34(1,0);
    Matrix44 M44(1,0);

    sun = Vector3::randS();
    dir = Vector3::randU();
    off = Vector3::randU();
    
    std::clog.precision(3);
    std::clog.setf(std::ios::fixed);

    checkMatrix(M33);
    checkMatrix(M34);
    checkMatrix(M44);
    
    checkMatrixFull(M44);
    
    return EXIT_SUCCESS;
}


