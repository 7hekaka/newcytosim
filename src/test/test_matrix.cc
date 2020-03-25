// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "random.h"
#include "tictoc.h"
#include "vecprint.h"

#include "dim.h"
#include "matsparsesym.h"
#include "matsparsesym1.h"
#include "matsparsesym2.h"
#include "matsparsesymblk.h"
#include "matsparseblk.h"

using namespace TicToc;

typedef MatrixSparseSymmetricBlock MatrixSparseSymmetricB;

const int N_RUN = 16;
const int N_MUL = 64;

#define PAD 4

real alpha = 2.0;
real beta = 1.0;


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
        s += std::abs( a[i] - b[i] );
    return s;
}

// fill lower triangle of matrix
void setIndices(size_t fill, size_t*& ii, size_t*& jj, size_t mx, size_t bs)
{
    delete[] ii;
    delete[] jj;
    
    ii = new size_t[fill];
    jj = new size_t[fill];
    
    for ( size_t n = 0; n < fill; ++n )
    {
        size_t i, j;
        do {
            i = RNG.pint(mx) / bs;
            j = RNG.pint(mx) / bs;
        } while ( i == j );
        ii[n] = bs * std::max(i,j);
        jj[n] = bs * std::min(i,j);
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
        z[n] = 0.0;
    }
    for ( size_t n = size; n < size+PAD; ++n )
    {
        x[n] = 0.0;
        y[n] = 0.0;
        z[n] = 0.0;
    }
}

//------------------------------------------------------------------------------
#pragma mark -


template <typename MATRIXA, typename MATRIXB>
void compare(size_t size,  MATRIXA & mat1, MATRIXB& mat2, size_t fill)
{
    real * tmp1 = new_real(size*size);
    real * tmp2 = new_real(size*size);
    
    mat1.reset();
    mat2.reset();
    
    mat1.resize(size);
    mat2.resize(size);
    
    for ( size_t n = 0; n < fill; ++n )
    {
        real a = 10.0 * RNG.preal();
        size_t ii = RNG.pint(size);
        size_t jj = RNG.pint(size);
        mat1(ii, jj) += a;
        mat2(ii, jj) += a;
    }
    
    for ( size_t nbc = DIM; nbc < size; nbc+=DIM )
    {
        size_t inx = DIM * ( RNG.pint(size-nbc) / DIM );
        std::clog << "Comparing matrices: size " << size << " inx " << inx << " nbc " << nbc << " ";
        
        zero_real(size*size, tmp1);
        
        mat1.addDiagonalBlock(tmp1, size, inx, nbc);
        
        //std::clog<<"mat1:\n";
        //VecPrint::print(std::clog, nbc, nbc, tmp1, size);
        
        zero_real(size*size, tmp2);
        mat2.addDiagonalBlock(tmp2, size, inx, nbc);
        
        //std::clog<<"mat2:\n";
        //VecPrint::print(std::clog, nbc, nbc, tmp2, size);
        
        real nrm = 0;
        for ( size_t n=0; n<size*size; ++n )
        {
            tmp1[n] -= tmp2[n];
            nrm += tmp1[n] * tmp1[n];
        }
        
        if ( nrm > 0 )
        {
            std::clog<<": diff:\n";
            VecPrint::print(std::clog, nbc, nbc, tmp1, size);
        }
        else
        {
            std::clog<<": identical\n";
            VecPrint::print(std::clog, nbc, nbc, tmp2, size);
        }
    }
    
    free_real(tmp1);
    free_real(tmp2);
}

#if ( DIM == 3 )
void fillMatrix(MatrixSparseSymmetricBlock& mat, const size_t i, const size_t j)
{
    Matrix33 M(alpha, -beta, beta, -beta, alpha, -beta, beta, -beta, alpha);
    
    mat.diag_block(i).add_half(M);
    mat.diag_block(j).add_half(M);

    if ( i > j )
        mat.block(i,j).add_full(M);
    else
        mat.block(j,i).add_full(M.transposed());
}


void fillMatrix(MatrixSparseBlock& mat, const size_t i, const size_t j)
{
    Matrix34 M(alpha, -beta, beta, -beta, alpha, -beta, beta, -beta, alpha);
    
    mat.diag_block(i).add_half(M);
    mat.diag_block(j).add_half(M);
    
    if ( i > j )
        mat.block(i,j).add_full(M);
    else
        mat.block(j,i).add_full(M.transposed());
}
#endif


template <typename MATRIX>
void fillMatrix(MATRIX& mat, const size_t i, const size_t j)
{
#if ( DIM == 3 )
    Matrix33 M(alpha, -beta, beta, -beta, alpha, -beta, beta, -beta, alpha);
#elif ( DIM == 2 )
    Matrix22 M(alpha, -beta, beta, alpha);
#else
    Matrix11 M(alpha);
#endif
    
    for ( int x = 0; x < DIM; ++x )
    for ( int y = x; y < DIM; ++y )
        mat(i+y, i+x) += M(y,x);
    
    if ( i > j )
    {
        for ( int x = 0; x < DIM; ++x )
        for ( int y = 0; y < DIM; ++y )
            mat(i+y, j+x) += M(y,x);
    }
    else
    {
        for ( int x = 0; x < DIM; ++x )
        for ( int y = 0; y < DIM; ++y )
            mat(j+y, i+x) += M(x,y);
    }
    
    for ( int x = 0; x < DIM; ++x )
    for ( int y = x; y < DIM; ++y )
        mat(j+y, j+x) += M(y,x);
}

template <typename MATRIX>
void fillMatrixIso(MATRIX& mat, const size_t i, const size_t j)
{
    mat(i, i) += alpha;
    if ( i > j )
        mat(i, j) -= beta;
    else
        mat(j, i) -= beta;
    mat(j, j) += alpha;
}

//------------------------------------------------------------------------------
#pragma mark -

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

    printf("  check %+16.6f  %+16.6f ", sum1, sum2);
}


template <typename MATRIX>
void testMatrix(MATRIX & mat,
                const size_t size, real const* x, real const* y, real * z,
                const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);

    tic();
    for ( int ii=0; ii<N_RUN; ++ii )
    {
        mat.reset();
        for ( size_t n=0; n<fill; ++n )
            fillMatrix(mat, iny[n], inx[n]);
    }
    double ts = toc();
    mat.prepareForMultiply(1);

    tic();
    for ( int n=0; n<N_RUN*N_MUL; ++n )
    {
        mat.vecMulAdd(y, z);
        mat.vecMulAdd(x, z);
    }
    double t1 = toc();
    
    tic();
    for ( int n=0; n<N_RUN*N_MUL; ++n )
    {
        mat.vecMulAdd_ALT(x, z);
        mat.vecMulAdd_ALT(y, z);
    }
    double t2 = toc();

    printf("\n %20s : ", mat.what().c_str());
    printf("set %8.3f  mul %8.3f  alt %8.3f", ts, t1, t2);
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
    #pragma omp parallel for num_threads(4)
    for ( int i = 0; i < size; i += CHK )
        mat.vecMulAdd(x, z, i, i+CHK);
    real sum1 = checksum(size, y, z);

    zero_real(size, z);
    #pragma omp parallel for num_threads(8)
    for ( int i = 0; i < size; i += CHK )
        mat.vecMulAdd(x, z, i, i+CHK);
    real sum2 = checksum(size, y, z);

    printf("  check %+16.6f  %+16.6f", sum1, sum2);
}


template <typename MATRIX>
void testMatrixParallel(MATRIX & mat,
                        const size_t size, real const* x, real const* y, real * z,
                        const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);

    tic();
    for ( int ii=0; ii<N_RUN; ++ii )
    {
        mat.reset();
        for ( size_t n=0; n<fill; ++n )
            fillMatrix(mat, iny[n], inx[n]);
    }
    double ts = toc();
    mat.prepareForMultiply(1);

    tic();
    for ( int n=0; n<N_RUN*N_MUL; ++n )
    {
        #pragma omp parallel for num_threads(2)
        for ( int i = 0; i < size; i += CHK )
        {
            mat.vecMulAdd(y, z, i, i+CHK);
            mat.vecMulAdd(x, z, i, i+CHK);
        }
    }
    double t2 = toc();

    tic();
    for ( int n=0; n<N_RUN*N_MUL; ++n )
    {
        #pragma omp parallel for num_threads(4)
        for ( int i = 0; i < size; i += CHK )
        {
            mat.vecMulAdd(y, z, i, i+CHK);
            mat.vecMulAdd(x, z, i, i+CHK);
        }
    }
    double t4 = toc();

    printf("\n %20s threaded mul :  x2  %8.3f  x4  %8.3f", mat.what().c_str(), t2, t4);
    checkMatrixParallel(mat, size, x, y, z);
}

#endif


/// multidimensional isotropic multiplication
template <typename MATRIX>
void testMatrixIso(MATRIX & mat,
                   const size_t size, real const* x, real const* y, real * z,
                   const size_t fill, size_t inx[], size_t iny[])
{
    tic();
    for ( size_t ii=0; ii<N_RUN; ++ii )
    {
        mat.reset();
        for ( size_t n=0; n<fill; ++n )
            fillMatrixIso(mat, inx[n], iny[n]);
    }
    double t2 = toc();
    
    tic();
    for ( size_t ii=0; ii<N_RUN; ++ii )
    {
        mat.prepareForMultiply(DIM);
        for ( int n=0; n<N_MUL; ++n )
#if ( DIM >= 3 )
            mat.vecMulAddIso3D(x, z);
#else
            mat.vecMulAddIso2D(x, z);
#endif
    }
    double t3 = toc();
    
    printf("  isoset %8.3f  isomul %8.3f", t2, t3);
}


void testMatrices(const size_t size, const size_t fill)
{
    printf("------ %iD size %lu  filled %.1f %% :", DIM, size, fill*100.0/size/size);
    //MatrixSparseSymmetric  mat0;
    MatrixSparseSymmetric1 mat1;
    //MatrixSparseSymmetric2 mat2;
    MatrixSparseSymmetricB mat3;
    MatrixSparseBlock      mat4;

    size_t * inx = nullptr;
    size_t * iny = nullptr;
    
    setIndices(fill, iny, inx, size, DIM);
    
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;

    setVectors(DIM*size, x, y, z);
    alpha = RNG.sreal();
    
    //testMatrix(mat0, size, x, y, z, fill, inx, iny);
    testMatrix(mat1, size, x, y, z, fill, inx, iny);
    //testMatrix(mat2, size, x, y, z, fill, inx, iny);
    testMatrix(mat3, size, x, y, z, fill, inx, iny);
    testMatrix(mat4, size, x, y, z, fill, inx, iny);
#ifdef _OPENMP
    testMatrixParallel(mat4, size, x, y, z, fill, inx, iny);
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
#pragma mark -

const real dir[4] = {  2, 1, -1, 3 };
const real vec[4] = { -1, 3,  1, 2 };

void fillMatrixBlock(MatrixSparseSymmetricBlock& mat, const size_t fill, size_t inx[], size_t iny[])
{
    SquareBlock S = SquareBlock::outerProduct(dir);
    SquareBlock U = SquareBlock::outerProduct(dir, vec);
    
    for ( size_t n=0; n<fill; ++n )
    {
        size_t ii = inx[n] - inx[n] % SquareBlock::dimension();
        size_t jj = iny[n] - iny[n] % SquareBlock::dimension();
        mat.diag_block(ii).sub_half(S);
        mat.diag_block(jj).add_half(S);
        mat.block(ii, jj).add_full(U);
    }
}

/**
This compares the Scalar and SIMD implementations of one matrix
*/
 void testMatrixBlock(MatrixSparseSymmetricBlock & mat,
                     const size_t size, real const* x, real const* y, real * z,
                     const size_t fill, size_t inx[], size_t iny[])
{
    mat.resize(size);
    
    tic();
    for ( int r = 0; r < N_RUN; ++r )
    {
        mat.reset();
        fillMatrixBlock(mat, fill, inx, iny);
    }
    double ts = toc();

    mat.prepareForMultiply(1);
    
    // calculate and compare sum for two methods:
    copy_real(size+PAD, y, z);
    mat.vecMulAdd_ALT(x, z);
    real sum = checksum(size, z, x);
    
    copy_real(size+PAD, y, z);
    mat.vecMulAdd(x, z);
    real res = checksum(size, z, x);

    unsigned long long time = __rdtsc();
    for ( int n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( int m=0; m<N_MUL; ++m )
            mat.vecMulAdd_ALT(x, z);
    }
    double nop = N_MUL * N_RUN * mat.nbElements();
    double t1 = ( __rdtsc() - time ) / nop;

    time = __rdtsc();
    for ( int n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( int m = 0; m < N_MUL; ++m )
            mat.vecMulAdd(x, z);
    }
    double t2 = ( __rdtsc() - time ) / nop;
    
    printf("%6lu %18s ", size, mat.what().c_str());
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
    
    MatrixSparseSymmetricBlock mat;
    testMatrixBlock(mat, size, x, y, z, fill, inx, iny);
    
    free_real(x);
    free_real(y);
    free_real(z);
    delete[] inx;
    delete[] iny;
}

int compareInt(const void* p, const void* q)
{
    size_t i = *(size_t const*)p;
    size_t j = *(size_t const*)q;
    return ( i > j ) - ( j > i );
}

int main( int argc, char* argv[] )
{
#ifdef COMPILER_VERSION
    printf("Matrix test and timing code, compiled with %s\n", COMPILER_VERSION);
#endif

    RNG.seed();
#if ( 0 )
        // small tests to check correctness:
        MatrixSparseSymmetric1 mat1;
        MatrixSparseSymmetricB mat3;
        MatrixSparseBlock mat4;
        
        compare(4*3, mat1, mat3, 1<<4);
        compare(4*7, mat1, mat3, 1<<5);
        compare(4*11, mat1, mat3, 1<<6);
        compare(4*33, mat1, mat3, 1<<16);
        compare(4*3, mat1, mat4, 1<<16);
        compare(4*7, mat1, mat4, 1<<16);
        compare(4*11, mat1, mat4, 1<<16);
        compare(4*33, mat1, mat4, 1<<16);
#endif
#if ( 0 )
        testMatrices(6, 1);
        testMatrices(12, 1);
        testMatrices(DIM*7, 1);
        testMatrices(DIM*17, 2);
        testMatrices(DIM*33, 1111);
#endif
#if ( 0 )
        printf("\ntest_matrix BLOCK_SIZE %i (%s)", DIM, SquareBlock::what().c_str());
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
#if ( 1 )
        //testMatrices(DIM*17, 23);
        //testMatrices(DIM*91, 1<<12);
        testMatrices(DIM*196, 1<<11);
        testMatrices(DIM*436, 1<<13);
        testMatrices(DIM*714, 1<<14);
        testMatrices(DIM*1358, 1<<15);
        testMatrices(DIM*2130, 1<<14);
        testMatrices(DIM*2130, 1<<15);
#endif
#if ( 0 )
        //testMatrices(DIM*17, 23);
        size_t dim[5] = { 0 };
        for ( int i = 0; i < 5; ++i ) dim[i] = RNG.pint(1<<(i+7));
        qsort(dim, 5, sizeof(size_t), compareInt);
        for ( int i = 0; i < 5; ++i )
            testMatrices(DIM*dim[i], RNG.pint(dim[i]*dim[i]));
#endif
return EXIT_SUCCESS;
}


