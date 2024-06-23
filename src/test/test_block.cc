// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
/// FJN 26.04.2020

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "vecprint.h"
#include "timer.h"
#include "random.h"
#include "blas.h"

#include "vector2.h"
#include "vector3.h"
#include "vector4.h"

#include "matrix33.h"
#include "matrix34.h"
#include "matrix44.h"
#include "matfull.h"


real diff(size_t size, real const* a, real const* b)
{
    real s = 0;
    for ( size_t i=0; i<size; ++i )
        s += abs_real( a[i] - b[i] );
    return s;
}


template <typename MATRIX, typename VECTOR>
void checkMatrix(MATRIX & mat, VECTOR dir, VECTOR off, VECTOR sun)
{
    mat = MATRIX::outerProduct(dir, off);
    VECTOR vec = mat.vecmul(sun);
    VECTOR vik = mat.trans_vecmul(sun);

    std::clog << "block " << std::setw(3) << mat.what() << ":  ";
    std::clog << vec << " | " << vik << "     ";
    std::clog << mat << "\n";
    //std::clog << mat.transposed() << "\n";

    //printf("  check %+16.6f  %+16.6f ", diff);
}

void checkMatrix44(Matrix44 const& src, Vector4 sun)
{
    const size_t SUP = 4;
    MatrixFull mat;
    mat.resize(SUP);
    
    for ( size_t i = 0; i < SUP; ++i )
    for ( size_t j = 0; j < SUP; ++j )
        mat(i,j) = src(i,j);
    
    Vector4 vec, vik(0,0,0), vok;
    mat.vecMul0(sun, vec);
    mat.vecMulAdd(sun, vik);
    mat.vecMul(sun, vok);

    std::clog << "m-full 16:  ";
    std::clog << vec << " = " << vik << " = " << vok << "\n";
    std::clog << mat << "\n";

    //printf("  check %+16.6f  %+16.6f ", diff);
}


template <typename MATRIX>
void speedBLAS(int cnt, MATRIX const& mat, real* src, real * x, real * y, real * z)
{
    const int N = (int)mat.size();
    real* mem = new_real(N*N);
    for ( int i = 0; i < N; ++i )
    for ( int j = 0; j < N; ++j )
        mem[i+N*j] = mat(i,j);
    
    blas::xgemv('N', N, N, 1.0, mem, N, src, 1, 0.0, x, 1);
    VecPrint::head(N, x, 2);
    tick();
    for ( int n = 0; n < cnt; ++n )
    {
        blas::xgemv('N', N, N, 1.0, mem, N, x, 1, 0.0, y, 1);
        blas::xgemv('N', N, N, 1.0, mem, N, y, 1, 0.0, z, 1);
        blas::xgemv('N', N, N, 1.0, mem, N, z, 1, 0.0, x, 1);
    }
    printf("  DGEMV %5.2f\n", tock(cnt));

    // transpose matrix
    for ( int i = 0  ; i < N; ++i )
    for ( int j = i+1; j < N; ++j )
    {
        real tmp = mem[i+N*j];
        mem[i+N*j] = mem[j+N*i];
        mem[j+N*i] = tmp;
    }
    
    blas::xgemv('T', N, N, 1.0, mem, N, src, 1, 0.0, x, 1);
    VecPrint::head(N, x, 2);
    tick();
    for ( int n = 0; n < cnt; ++n )
    {
        blas::xgemv('T', N, N, 1.0, mem, N, x, 1, 0.0, y, 1);
        blas::xgemv('T', N, N, 1.0, mem, N, y, 1, 0.0, z, 1);
        blas::xgemv('T', N, N, 1.0, mem, N, z, 1, 0.0, x, 1);
    }
    printf(" tDGEMV %5.2f\n", tock(cnt));
    
    free_real(mem);
}
    
    
void speedMatrix(int size, int cnt)
{
    MatrixFull mat;
    mat.resize(size);
    
    real * x = new_real(size);
    real * y = new_real(size);
    real * z = new_real(size);
    real * s = new_real(size);

    for ( int i = 0; i < size; ++i )
        s[i] = RNG.sreal();
    
    for ( int i = 0; i < size; ++i )
    for ( int j = 0; j < size; ++j )
        mat(i,j) = RNG.preal();

    printf("Matrix %s size %i\n", mat.what().c_str(), size);

    mat.vecMul0(s, x);
    VecPrint::head(size, x, 2);
    tick();
    for ( int n = 0; n < cnt; ++n )
    {
        mat.vecMul0(x, y);
        mat.vecMul0(y, z);
        mat.vecMul0(z, x);
    }
    printf(" SCALAR %5.2f\n", tock(cnt));

    mat.vecMul(s, x);
    VecPrint::head(size, x, 2);
    tick();
    for ( int n = 0; n < cnt; ++n )
    {
        mat.vecMul(x, y);
        mat.vecMul(y, z);
        mat.vecMul(z, x);
    }
    printf("   AVX? %5.2f\n", tock(cnt));

    speedBLAS(cnt, mat, s, x, y, z);
    
    mat.transpose();
    mat.transVecMul(s, x);
    VecPrint::head(size, x, 2);
    tick();
    for ( int n = 0; n < cnt; ++n )
    {
        mat.transVecMul(x, y);
        mat.transVecMul(y, z);
        mat.transVecMul(z, x);
    }
    printf(" TRANSP %5.2f\n", tock(cnt));
    
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
    
    if ( 1 )
    {
        Matrix33 M33(1,0);
        Matrix34 M34(1,0);
        Matrix44 M44(1,0);

        Vector3 dir(1, 0, 0);
        Vector3 off(0, 1, 0);
        
        dir = Vector3::randS();
        off = Vector3::randS();

        Vector3 vec = Vector3::randS();

        checkMatrix(M33, dir, off, vec);
        checkMatrix(M34, dir, off, vec);

        checkMatrix(M44, Vector4(dir), Vector4(off), Vector4(vec));
        checkMatrix44(M44, Vector4(vec));
    }
    
    speedMatrix(119, 1<<12);
}


