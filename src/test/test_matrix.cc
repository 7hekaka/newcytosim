// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "vecprint.h"
#include "random.h"

#include "dim.h"
#include "matrix33.h"
#include "sparmatsym.h"
#include "sparmatsym1.h"
#include "sparmatsym2.h"
#include "sparmatsymblk.h"
#include "sparmatsymblkdiag.h"
#include "sparmatblk.h"


typedef SparMatSym1       SparMat1;
typedef SparMatSym2       SparMat2;
typedef SparMatBlk        SparMatA;
typedef SparMatSymBlk     SparMatB;
typedef SparMatSymBlkDiag SparMatD;

typedef unsigned SIZE_T;

const size_t N_RUN = 16;
const size_t N_MUL = 99;

#define PAD 4

real alpha = 2.0;
real beta = -1.0;

#if 1
/// keeping time using Intel's cycle counters
unsigned long long rdt = 0;
/// start timer
inline void tic() { rdt = __rdtsc(); }
/// return time since last 'tic()'
inline double toc(double num) { return double(__rdtsc()-rdt) / num; }
#else
# include "tictoc.h"
using namespace TicToc;
#endif


//------------------------------------------------------------------------------
#pragma mark - Functions

real checksum(size_t size, real const* vec, real const* ptr)
{
    real s = 0;
    for ( size_t i=0; i<size; ++i )
        s += vec[i] * ptr[i];
    return s;
}


real diff(size_t size, real const* a, real const* b)
{
    real s = 0;
    for ( size_t i=0; i<size; ++i )
        s += abs_real( a[i] - b[i] );
    return s;
}

// fill lower triangle of matrix
void setIndices(size_t fill, size_t*& ii, size_t*& jj, unsigned sup, size_t dim)
{
    delete[] ii;
    delete[] jj;
    
    ii = new size_t[fill];
    jj = new size_t[fill];
    
    for ( size_t n = 0; n < fill; ++n )
    {
        size_t i, j;
        do {
            i = RNG.pint32(sup) / dim;
            j = RNG.pint32(sup) / dim;
        } while ( i == j );
        ii[n] = dim * std::max(i,j);
        jj[n] = dim * std::min(i,j);
    }
}

void setVectors(size_t size, real*& x, real*& y, real*& z)
{
    free_real(x);
    free_real(y);
    free_real(z);

    x = new_real(size+PAD);
    y = new_real(size+PAD);
    z = new_real(size+PAD);

    for ( size_t n = 0; n < size; ++n )
    {
        x[n] = RNG.sreal();
        y[n] = RNG.preal();
        z[n] = 0;
    }
    for ( size_t n = size; n < size+PAD; ++n )
    {
        x[n] = 0;
        y[n] = 0;
        z[n] = 0;
    }
}

template <typename MATRIX>
void fillMatrix(MATRIX& mat, const size_t i, const size_t j)
{
    assert_true( i+DIM-1 < mat.size() );
    assert_true( j < i );  // below diagonal
    
#if ( DIM == 3 )
    Matrix33 M(alpha, beta, -beta, beta, alpha, beta, -beta, beta, alpha);
#elif ( DIM == 2 )
    Matrix22 M(alpha, beta, -beta, alpha);
#else
    Matrix11 M(alpha);
#endif
    
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = x; y < DIM; ++y )
    {
        real v = M(y,x);
        mat(i+y, i+x) += v;
        mat(j+y, j+x) += v;
    }
    if ( 1 )
    {
        size_t I = std::max(i, j);
        size_t J = std::min(i, j);
        for ( size_t x = 0; x < DIM; ++x )
        for ( size_t y = 0; y < DIM; ++y )
            mat(I+y, J+x) += M(y,x);
    }
}


template <typename MATRIX>
void fillMatrixIso(MATRIX& mat, const size_t i, const size_t j)
{
    mat(i, i) += alpha;
    if ( i > j )
        mat(i, j) += beta;
    else
        mat(j, i) += beta;
    mat(j, j) += alpha;
}

//------------------------------------------------------------------------------
#pragma mark - Compare two matrices


template <typename MATRIX, typename MATROX>
void compareMatrix(size_t size,  MATRIX & mat1, MATROX& mat2, size_t fill)
{
    real * tmp1 = new_real(size*size);
    real * tmp2 = new_real(size*size);

    mat1.resize(size);
    mat2.resize(size);
    
    mat1.reset();
    mat2.reset();

    for ( size_t n = 0; n < fill; ++n )
    {
        real a = 0.1; //RNG.preal();
        size_t ii = RNG.pint32(size);
        size_t jj = RNG.pint32(size);
        if ( ii != jj )
        {
            size_t i = std::max(ii, jj);
            size_t j = std::min(ii, jj);
            
            mat1(i, i) += a;
            mat1(i, j) -= a;
            mat1(j, j) += a;

            mat2(i, i) += a;
            mat2(i, j) -= a;
            mat2(j, j) += a;
        }
    }
    
    bool error = false;
    for ( size_t cnt = DIM; cnt < size; cnt += DIM )
    {
        size_t inx = DIM * ( RNG.pint32(size-cnt) / DIM );
        
        zero_real(size*size, tmp1);
        zero_real(size*size, tmp2);

        mat1.addDiagonalBlock(tmp1, size, inx, cnt);
        mat2.addDiagonalBlock(tmp2, size, inx, cnt);
        
        for ( size_t i = 0; i < size; ++i )
        for ( size_t j = 0; j < size; ++j )
        {
            real e = abs_real(tmp1[i+size*j]-tmp2[i+size*j]);
            if ( e > 0.1 )
            {
                std::clog << "Error " << i << " " << j << "\n";
                error = true;
            }
        }
        
        if ( error )
        {
            std::clog << "Size " << size << " : " << mat1.what() << "  " << mat2.what();
            std::clog << " inx " << inx << " + " << cnt << " ";
            std::clog << ": error\n";
            std::clog << mat2.what() << ":\n";
            VecPrint::print(std::clog, cnt, cnt, tmp2, size);
            
            zero_real(size*size, tmp2);
            mat1.addDiagonalBlock(tmp2, size, inx, cnt);
            std::clog << mat1.what() << ":\n";
            VecPrint::print(std::clog, cnt, cnt, tmp2, size);
            break;
        }
    }
    if ( !error )
    {
        std::clog << mat1.what() << " and " << mat2.what();
        std::clog << " are identical\n";
        //VecPrint::print(std::clog, cnt, cnt, tmp2, size);
    }

    free_real(tmp1);
    free_real(tmp2);
}


//------------------------------------------------------------------------------
#pragma mark - Test Sparse Matrix

template <typename MATRIX>
void checkMatrix(MATRIX & mat, const size_t size,
                 real const* x, real const* y, real * z)
{
    assert_true(mat.size() == size);
    zero_real(size, z);
    mat.vecMulAdd(x, z);
    real sum1 = checksum(size, y, z);
    
    zero_real(size, z);
    mat.vecMulAdd_ALT(x, z);
    real sum2 = checksum(size, y, z);
    
    mat.vecMul(x, z);
    real sum3 = checksum(size, y, z);

    printf("  check %+16.6f %+16.6f %+16.6f ", sum1, sum2, sum3);
}


template <typename MATRIX>
void testMatrix(MATRIX & mat,
                const size_t size, real const* x, real const* y, real * z,
                const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);
    const size_t CNT = N_RUN * N_MUL;
    const size_t DIV = CNT << 18;

    tic();
    for ( size_t ii=0; ii<N_RUN; ++ii )
    {
        mat.reset();
        for ( size_t n=0; n<fill; ++n )
            fillMatrix(mat, iny[n], inx[n]);
        //            mat(iny[n], inx[n]) += alpha;
    }
    double ts = toc(DIV);
    mat.prepareForMultiply(1);

    tic();
    for ( size_t n=0; n<CNT; ++n )
    {
        mat.vecMulAdd(y, z);
        mat.vecMulAdd(x, z);
    }
    double t1 = toc(DIV);
    
    tic();
    for ( size_t n=0; n<CNT; ++n )
    {
        mat.vecMulAdd_ALT(x, z);
        mat.vecMulAdd_ALT(y, z);
    }
    double t2 = toc(DIV);

    tic();
    for ( size_t n=0; n<CNT; ++n )
    {
        mat.vecMul(y, z);
        mat.vecMul(x, z);
    }
    double t3 = toc(DIV);

    printf("\n%-32s ", mat.what().c_str());
    printf("set %9.2f  muladd %9.2f  alt %9.2f  mul %9.2f", ts, t1, t2, t3);
    checkMatrix(mat, size, x, y, z);
}


#ifdef _OPENMP
#include <omp.h>

constexpr size_t CHK = DIM*2;


template <typename MATRIX>
void checkMatrixParallel(MATRIX & mat, const size_t size,
                         real const* x, real const* y, real * z)
{
    assert_true(mat.size() == size);
    
    zero_real(size, z);
    for ( size_t i = 0; i < size; ++i )
        mat.vecMulAdd(x, z, i, i+1);
    real sum = checksum(size, y, z);
    printf(" check %+16.6f", sum);

    for ( int i = 0; i < 5; ++i )
    {
        zero_real(size, z);
        omp_set_num_threads(1<<i);
        #pragma omp parallel for
        for ( size_t u = 0; u < size; u += CHK )
            mat.vecMulAdd(x, z, u, u+CHK);
        real sum1 = checksum(size, y, z);
        if ( sum != sum1 )
            printf(" %luT %+16.6f", 1<<i, sum1);
    }
}


template <typename MATRIX>
void testMatrixParallel(MATRIX & mat,
                        const size_t size, real const* x, real const* y, real * z,
                        const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);
    mat.reset();
    for ( size_t n=0; n<fill; ++n )
        fillMatrix(mat, iny[n], inx[n]);
    mat.prepareForMultiply(1);

    const size_t sup = 7;
    double t[sup] = { 0 };
    
    for ( int i = 0; i < sup; ++i )
    {
        omp_set_num_threads(1<<i);
        tic();
        for ( size_t j=0; j < N_RUN*N_MUL; ++j )
        {
            #pragma omp parallel for
            for ( size_t i = 0; i < size; i += CHK )
            {
                mat.vecMulAdd(y, z, i, i+CHK);
                mat.vecMulAdd(x, z, i, i+CHK);
            }
        }
        t[i] = toc();
    }

    printf("\n%-28s", mat.what().c_str());
    for ( int i = 0; i < sup; ++i )
        printf("  %luT %6.3f", 1<<i, t[i]);
}

#endif

#if ( DIM == 1 )
#   define VECMULADDISO vecMulAdd
#elif ( DIM == 2 )
#   define VECMULADDISO vecMulAddIso2D
#elif ( DIM == 3 )
#   define VECMULADDISO vecMulAddIso3D
#endif

/// multidimensional isotropic multiplication
template <typename MATRIX>
void testIsoMatrix(MATRIX & mat,
                   const size_t size, real const* x, real const* y, real * z,
                   const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);
    const size_t DIV = ( N_RUN * N_MUL ) << 18;

    tic();
    for ( size_t i=0; i<N_RUN; ++i )
    {
        mat.reset();
        for ( size_t n=0; n<fill; ++n )
            fillMatrixIso(mat, inx[n], iny[n]);
    }
    double ts = toc(DIV);
    
    tic();
    for ( size_t i=0; i<N_RUN; ++i )
    {
        mat.prepareForMultiply(DIM);
        for ( size_t n=0; n<N_MUL; ++n )
        {
            mat.VECMULADDISO(x, z);
            mat.VECMULADDISO(y, z);
        }
    }
    double tm = toc(DIV);
    
    printf("\n%-29s ", mat.what().c_str());
    printf("isoset %9.2f  isomul %9.2f", ts, tm);

    zero_real(size, z);
    mat.VECMULADDISO(x, z);
    printf("  isocheck %+16.6f ", checksum(size, x, z));
    
#if ( 0 )
    // compute checksum with another matrix class
    SparMat1 ful;
    ful.resize(size*DIM);
    for ( size_t j = 0; j < size; ++j )
    for ( size_t i = j; i < size; ++i )
    {
        real a = mat(i, j);
        if ( abs_real(a) > 0 )
        {
            for ( int d = 0; d < DIM; ++d )
                ful(DIM*i+d, DIM*j+d) = a;
        }
    }
    ful.prepareForMultiply(1);
    zero_real(size, z);
    ful.vecMulAdd(x, z);
    printf("  %+16.6f ", checksum(size, x, z));
#endif
}


void testMatrices(const size_t size, const size_t fill)
{
    printf("------ %iD size %lu  filled %.1f %% :", DIM, size, fill*100.0/size/size);

    size_t * inx = nullptr;
    size_t * iny = nullptr;
    
    setIndices(fill, iny, inx, size, DIM);
    
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;

    setVectors(DIM*size, x, y, z);
    alpha = RNG.sreal();
    
    SparMat1 mat1; testMatrix(mat1, size, x, y, z, fill, inx, iny);
    SparMat2 mat2; testMatrix(mat2, size, x, y, z, fill, inx, iny);
#if 1
    SparMatB mat3; testMatrix(mat3, size, x, y, z, fill, inx, iny);
    SparMatD mat4; testMatrix(mat4, size, x, y, z, fill, inx, iny);
    //testMatrix(mat4, size, x, y, z, fill, inx, iny);
    SparMatA mat5; testMatrix(mat5, size, x, y, z, fill, inx, iny);
#endif
#ifdef _OPENMP
    testMatrixParallel(mat5, size, x, y, z, fill, inx, iny);
    checkMatrixParallel(mat5, size, x, y, z);
#endif
#if 0
    testIsoMatrix(mat1, size, x, y, z, fill, inx, iny);
    testIsoMatrix(mat2, size, x, y, z, fill, inx, iny);
#endif
    printf("\n");
    
#if ( 0 )
    std::ofstream os1("mat1.txt");
    std::ofstream os3("mat3.txt");
    mat1.printSparse(os1, 0);
    mat3.printSparse(os3, 0);
#endif
    
    free_real(x);
    free_real(y);
    free_real(z);
    delete[] inx;
    delete[] iny;
}


//------------------------------------------------------------------------------
#pragma mark - Block matrices

const real dir[4] = {  2, 1, -1, 3 };
const real vec[4] = { -1, 3,  1, 2 };

template < typename MATRIX >
void fillMatrixBlock(MATRIX& mat, const size_t fill, size_t inx[], size_t iny[])
{
    typename MATRIX::Block S = MATRIX::Block::outerProduct(dir);
    typename MATRIX::Block U = MATRIX::Block::outerProduct(dir, vec);
    
    for ( size_t n=0; n<fill; ++n )
    {
        size_t ii = inx[n] - inx[n] % MATRIX::Block::dimension();
        size_t jj = iny[n] - iny[n] % MATRIX::Block::dimension();
        mat.diag_block(ii).sub_half(S);
        mat.diag_block(jj).add_half(S);
        mat.block(ii, jj).add_full(U);
    }
}

template < typename MATRIX >
void fillBlockMatrix(MATRIX& mat, const size_t i, const size_t j)
{
    typename MATRIX::Block M(alpha, beta, -beta, beta, alpha, beta, -beta, beta, alpha);
    
    mat.diag_block(i).add_half(M);
    mat.diag_block(j).add_half(M);

    if ( i > j )
        mat.block(i,j).add_full(M);
    else
        mat.block(j,i).add_full(M.transposed());
}

/**
This compares the Scalar and SIMD implementations of one matrix
*/
 void testMatrixBlock(SparMatSymBlk & mat,
                     const size_t size, real const* x, real const* y, real * z,
                     const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);
    
    tic();
    for ( size_t r = 0; r < N_RUN; ++r )
    {
        mat.reset();
        fillMatrixBlock(mat, fill, inx, iny);
    }
    double ts = toc(N_RUN);

    mat.prepareForMultiply(1);
    
    // calculate and compare sum for two methods:
    copy_real(size+PAD, y, z);
    mat.vecMulAdd_ALT(x, z);
    real sum = checksum(size, z, x);
    
    copy_real(size+PAD, y, z);
    mat.vecMulAdd(x, z);
    real res = checksum(size, z, x);

    tic();
    for ( size_t n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( size_t m=0; m<N_MUL; ++m )
            mat.vecMulAdd_ALT(x, z);
    }
    double nop = N_MUL * N_RUN * mat.nbElements();
    double t1 = toc(nop);

    tic();
    for ( size_t n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( size_t m = 0; m < N_MUL; ++m )
            mat.vecMulAdd(x, z);
    }
    double t2 = toc(nop);
    
    printf("%6lu %26s ", size, mat.what().c_str());
    printf("set %8.1f mul %8.1f  alt %8.1f", ts, t1, t2);
    printf(" :  checksum  %+24.16f %+24.16f", sum, res);
    if ( sum != res )
        printf("  failed!\n");
    else
        printf("\n");
}


void testMatrixBlock(const size_t size, const size_t fill)
{
    size_t * inx = nullptr;
    size_t * iny = nullptr;
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;

    setIndices(fill, iny, inx, size, DIM);
    setVectors(size, x, y, z);
    alpha = RNG.sreal();
    
    SparMatSymBlk mat;
    testMatrixBlock(mat, size, x, y, z, fill, inx, iny);
    
    free_real(x);
    free_real(y);
    free_real(z);
    delete[] inx;
    delete[] iny;
}

//------------------------------------------------------------------------------
#pragma mark - Main

int compareInt(const void* A, const void* B)
{
    size_t a = *static_cast<size_t const*>(A);
    size_t b = *static_cast<size_t const*>(B);
    return ( a > b ) - ( b > a );
}

int main( int argc, char* argv[] )
{
    printf("Matrix test and timing code --- real %lu bytes --- %s\n", sizeof(real), __VERSION__);

    RNG.seed();
#if ( 0 )
    SparMat1 mat1;
    SparMat2 mat2;
    SparMatA matA; //asymetric!
    SparMatB mat3;
    SparMatD mat4;
    
    // check correctness:
    compareMatrix(4*3, mat1, mat2, 1<<4);
    compareMatrix(4*5, mat1, mat3, 1<<4);
    compareMatrix(4*7, mat2, mat3, 1<<5);
#endif
#if ( 0 )
    compareMatrix(4*11, mat1, mat3, 1<<6);
    compareMatrix(4*33, mat1, mat3, 1<<16);
    compareMatrix(4*3, mat1, mat4, 1<<16);
    compareMatrix(4*7, mat1, mat4, 1<<16);
    compareMatrix(4*11, mat1, mat4, 1<<16);
    compareMatrix(4*33, mat1, mat4, 1<<16);
#endif
#if ( 0 )
    testMatrices(DIM*2, 1);
    testMatrices(DIM*3, 3);
    testMatrices(DIM*4, 4);
    testMatrices(DIM*7, 5);
    testMatrices(DIM*17, 7);
    testMatrices(DIM*33, 1111);
#endif
#if ( 0 )
    printf("\ntest_matrix BLOCK_SIZE %i (%s)", DIM, SparMatSymB::Block::what().c_str());
    size_t siz = DIM;
    for ( int i = 0; i < 14; ++i )
    {
        siz = 1 + ( 2 << i ) * RNG.preal();
        size_t fill = 1 + siz * siz * ( 0.01 * RNG.preal() );
        testMatrixBlock(siz, fill);
    }
#endif
#if ( 0 )
    testMatrixBlock(DIM*17, 2);
    testMatrixBlock(DIM*347, 1019);
    testMatrixBlock(DIM*753, 43039);
#endif
#if ( 0 )
    testMatrixBlock(2253, 1<<14);
    testMatrixBlock(2253, 1<<14);
    testMatrixBlock(2253, 1<<14);
    testMatrixBlock(2253, 1<<14);
#endif
#if ( 0 )
    testMatrixBlock(DIM*1251, 25821);
    testMatrixBlock(DIM*1785, 153034);
    testMatrixBlock(DIM*2311, 231111);
    //testMatrixBlock(DIM*3217, 671234);
#endif
#if ( 0 )
    //testMatrices(DIM*17, 23);
    //testMatrices(DIM*91, 1<<12);
    testMatrices(DIM*196, 1<<13);
    //testMatrices(DIM*436, 1<<15);
    testMatrices(DIM*8*94, 1<<16);
    //testMatrices(DIM*8*169, 1<<17);
    testMatrices(DIM*8*331, 1<<18);
#endif
    testMatrices(DIM*8*111, 1<<18);
#if ( 0 )
    //testMatrices(DIM*17, 23);
    size_t dim[5] = { 0 };
    for ( int i = 0; i < 5; ++i ) dim[i] = RNG.pint32(1<<(i+7));
    qsort(dim, 5, sizeof(size_t), compareInt);
    for ( int i = 0; i < 5; ++i )
        testMatrices(DIM*dim[i], RNG.pint32(dim[i]*dim[i]));
#endif
}

