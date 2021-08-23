// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
/// FJN 26.04.2020

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "vecprint.h"
#include "tictoc.h"
#include "random.h"
#include "blas.h"

#include "vector2.h"
#include "vector3.h"
#include "vector4.h"

#include "matrix33.h"
#include "matrix34.h"
#include "matrix44.h"
#include "matfull.h"

using namespace TicToc;
size_t NVAL = 15;

Vector4 dir(0, 1, 0);
Vector4 off(0.1, 0.2, -0.3);
Vector4 sun(0, 0, 0);

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
    Vector4 vec = mat.vecmul(sun);
    Vector4 vik = mat.trans_vecmul(sun);

    std::clog << vec << " | " << vik << "  ";
    std::clog << mat << "  " << mat.transposed() << "\n";

    //printf("  check %+16.6f  %+16.6f ", diff);
}


void checkMatrixFull(Matrix44 const& src)
{
    const size_t SUP = 4;
    MatrixFull mat;
    mat.resize(SUP);
    
    for ( size_t i = 0; i < SUP; ++i )
    for ( size_t j = 0; j < SUP; ++j )
        mat(i,j) = src(i,j);
    
    Vector3 vec, vik(0,0,0), vok;
    mat.vecMul0(sun, vec);
    mat.vecMulAdd(sun, vik);
    mat.vecMul(sun, vok);

    std::clog << vec << " | " << vik << " | " << vok << "\n";
    std::clog << mat << "\n";

    //printf("  check %+16.6f  %+16.6f ", diff);
}


template <typename MATRIX>
void speedBLAS(size_t cnt, MATRIX const& mat, real* src, real * x, real * y, real * z)
{
    size_t size = mat.size();
    real* mem = new_real(size*size);
    for ( size_t i = 0; i < size; ++i )
    for ( size_t j = 0; j < size; ++j )
        mem[i+size*j] = mat(i,j);
    
    const int N = (int)size;
    blas::xgemv('N', N, N, 1.0, mem, N, src, 1, 0.0, x, 1);
    VecPrint::print(std::min(size,NVAL), x);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        blas::xgemv('N', N, N, 1.0, mem, N, x, 1, 0.0, y, 1);
        blas::xgemv('N', N, N, 1.0, mem, N, y, 1, 0.0, z, 1);
        blas::xgemv('N', N, N, 1.0, mem, N, z, 1, 0.0, x, 1);
    }
    printf("  DGEMV %5.2f\n", toc(cnt));

    // transpose matrix
    for ( size_t i = 0  ; i < size; ++i )
    for ( size_t j = i+1; j < size; ++j )
    {
        real tmp = mem[i+size*j];
        mem[i+size*j] = mem[j+size*i];
        mem[j+size*i] = tmp;
    }
    
    blas::xgemv('T', N, N, 1.0, mem, N, src, 1, 0.0, x, 1);
    VecPrint::print(std::min(size,NVAL), x);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        blas::xgemv('T', N, N, 1.0, mem, N, x, 1, 0.0, y, 1);
        blas::xgemv('T', N, N, 1.0, mem, N, y, 1, 0.0, z, 1);
        blas::xgemv('T', N, N, 1.0, mem, N, z, 1, 0.0, x, 1);
    }
    printf(" tDGEMV %5.2f\n", toc(cnt));
    
    free_real(mem);
}
    
    
void speedMatrix(size_t size, size_t cnt)
{
    MatrixFull mat;
    mat.resize(size);
    
    real * x = new_real(size);
    real * y = new_real(size);
    real * z = new_real(size);
    real * s = new_real(size);

    for ( size_t i = 0; i < size; ++i )
        s[i] = RNG.sreal();
    
    for ( size_t i = 0; i < size; ++i )
    for ( size_t j = 0; j < size; ++j )
        mat(i,j) = RNG.preal();

    printf("Matrix %s size %lu\n", mat.what().c_str(), size);

    mat.vecMul0(s, x);
    VecPrint::print(std::min(size,NVAL), x);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        mat.vecMul0(x, y);
        mat.vecMul0(y, z);
        mat.vecMul0(z, x);
    }
    printf(" SCALAR %5.2f\n", toc(cnt));

    mat.vecMul(s, x);
    VecPrint::print(std::min(size,NVAL), x);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        mat.vecMul(x, y);
        mat.vecMul(y, z);
        mat.vecMul(z, x);
    }
    printf("   AVX? %5.2f\n", toc(cnt));

    speedBLAS(cnt, mat, s, x, y, z);
    
    mat.transpose();
    mat.transVecMul(s, x);
    VecPrint::print(std::min(size,NVAL), x);
    tic();
    for ( size_t n = 0; n < cnt; ++n )
    {
        mat.transVecMul(x, y);
        mat.transVecMul(y, z);
        mat.transVecMul(z, x);
    }
    printf(" TRANSP %5.2f\n", toc(cnt));
    
    free_real(x);
    free_real(y);
    free_real(z);
    free_real(s);
}


//------------------------------------------------------------------------------
#pragma mark -

int main( int argc, char* argv[] )
{
    printf("test_block --- %s\n", __VERSION__);
    std::clog.setf(std::ios::fixed);
    std::clog.precision(3);
    RNG.seed();
    
    if ( 0 )
    {
        Matrix33 M33(1,0);
        Matrix34 M34(1,0);
        Matrix44 M44(1,0);
        
        sun = Vector3::randS();
        dir = Vector3::randU();
        off = Vector3::randU();
        
        checkMatrix(M33);
        checkMatrix(M34);
        checkMatrix(M44);
        checkMatrixFull(M44);
    }
    
    speedMatrix(119, 1<<12);
}


