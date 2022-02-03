// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "vecprint.h"
#include "random.h"
#include "timer.h"
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

// number of multiplication in sequence
constexpr size_t N_MUL = 47;

// number of repeat of ( 1 prepare + N_MUL multiplications)
constexpr size_t N_RUN = 4;

constexpr size_t CNT = N_RUN * N_MUL;

#define PAD 4

size_t icnt_ = 0;
size_t * inx_ = nullptr;
size_t * iny_ = nullptr;

real alpha = 2.0;
real beta = -1.0;


//------------------------------------------------------------------------------
#pragma mark - Functions

int compareInt(const void* A, const void* B)
{
    size_t a = *static_cast<size_t const*>(A);
    size_t b = *static_cast<size_t const*>(B);
    return ( a > b ) - ( b > a );
}

real checksum(size_t size, real const* vec)
{
    real s = 0;
    for ( size_t i=0; i<size; ++i )
        s += vec[i];
    return s;
}


real diff(size_t size, real const* a, real const* b)
{
    real s = 0;
    for ( size_t i=0; i<size; ++i )
        s += abs_real( a[i] - b[i] );
    return s;
}

///set indices within [0, sup] that are multiples of 'dim'
void setIndices(size_t cnt, size_t sup, size_t dim)
{
    delete[] inx_;
    delete[] iny_;
    inx_ = nullptr;
    iny_ = nullptr;
    icnt_ = cnt;
    if ( cnt > 0 )
    {
        inx_ = new size_t[icnt_];
        iny_ = new size_t[icnt_];
        sup /= dim;
        for ( size_t n = 0; n < icnt_; ++n )
        {
            size_t i = RNG.pint32(sup);
            size_t j = RNG.pint32(sup);
            inx_[n] = dim * std::max(i,j);
            iny_[n] = dim * std::min(i,j);
        }
    }
}

void setVectors(size_t size, real*& x, real*& y, real*& z)
{
    free_real(x);
    free_real(y);
    free_real(z);
    x = nullptr;
    y = nullptr;
    z = nullptr;
    if ( size > 0 )
    {
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
}

//------------------------------------------------------------------------------
#pragma mark - Compare two matrices


template <typename MATRIX, typename MATROX>
void compareMatrix(size_t S, MATRIX & mat1, MATROX& mat2, size_t fill)
{
    size_t nume = S * S;
    real * tmp1 = new_real(nume);
    real * tmp2 = new_real(nume);

    mat1.resize(S);
    mat2.resize(S);
    
    mat1.reset();
    mat2.reset();

    for ( size_t n = 0; n < fill; ++n )
    {
        real a = 0.1; //RNG.preal();
        size_t ii = RNG.pint32(S);
        size_t jj = RNG.pint32(S);
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
    for ( size_t cnt = DIM; cnt < S; cnt += DIM )
    {
        size_t inx = DIM * ( RNG.pint32(S-cnt) / DIM );
        
        zero_real(nume, tmp1);
        zero_real(nume, tmp2);

        mat1.addDiagonalBlock(tmp1, S, inx, cnt);
        mat2.addDiagonalBlock(tmp2, S, inx, cnt);
        
        for ( size_t i = 0; i < S; ++i )
        for ( size_t j = 0; j < S; ++j )
        {
            real e = abs_real(tmp1[i+S*j]-tmp2[i+S*j]);
            if ( e > 0.1 )
            {
                std::clog << "Error " << i << " " << j << "\n";
                error = true;
            }
        }
        
        if ( error )
        {
            std::clog << "Size " << S << " : " << mat1.what() << "  " << mat2.what();
            std::clog << " inx " << inx << " + " << cnt << " ";
            std::clog << ": error\n";
            VecPrint::full(mat2.what(), cnt, cnt, tmp2, S);
            
            zero_real(nume, tmp2);
            mat1.addDiagonalBlock(tmp2, S, inx, cnt);
            VecPrint::full(mat1.what(), cnt, cnt, tmp2, S);
            break;
        }
    }
    if ( !error )
    {
        std::clog << mat1.what() << " and " << mat2.what();
        std::clog << " are identical\n";
        //VecPrint::full(cnt, cnt, tmp2, size);
    }

    free_real(tmp1);
    free_real(tmp2);
}

//------------------------------------------------------------------------------
#pragma mark - Test Sparse Matrix


template <typename MATRIX>
void fillMatrix(MATRIX& mat, const size_t i, const size_t j)
{
    assert_true( i+DIM <= mat.size() );
    assert_true( j+DIM <= mat.size() );
    
#if ( DIM == 3 )
    Matrix33 M(alpha, beta, -beta, beta, alpha, beta, -beta, beta, alpha);
#elif ( DIM == 2 )
    Matrix22 M(alpha, beta, -beta, alpha);
#else
    Matrix11 M(alpha);
#endif
    // this is a block on the diagonal:
    for ( size_t x = 0; x < DIM; ++x )
    for ( size_t y = x; y < DIM; ++y )
        mat(i+y, i+x) += M(y,x);
    
    if ( i != j )
    {
        // this is a block on the diagonal:
        for ( size_t x = 0; x < DIM; ++x )
        for ( size_t y = x; y < DIM; ++y )
            mat(j+y, j+x) += M(y,x);
        // off-diagonal
        size_t I = std::max(i, j);
        size_t J = std::min(i, j);
        for ( size_t x = 0; x < DIM; ++x )
        for ( size_t y = 0; y < DIM; ++y )
            mat(I+y, J+x) += M(y,x);
    }
}


template <typename MATRIX>
void checkMatrix(MATRIX & mat, real const* x, real * z)
{
    size_t S = mat.size();
    zero_real(S, z);
    mat.vecMulAdd(x, z);
    real sum1 = checksum(S, z);
    
    zero_real(S, z);
    mat.vecMulAdd_ALT(x, z);
    real sum2 = checksum(S, z);
    
    mat.vecMul(x, z);
    real sum3 = checksum(S, z);
    
    printf("  check %+16.6f %+16.6f %+16.6f ", sum1, sum2, sum3);
}


template <typename MATRIX>
void testMatrix(MATRIX & mat, real const* x, real const* y, real * z)
{
    tick();
    for ( size_t n=0; n<N_RUN; ++n )
    {
        mat.reset();
        for ( size_t i=0; i<icnt_; ++i )
            fillMatrix(mat, iny_[i], inx_[i]);
    }
    double ts = tock(N_RUN);
    mat.prepareForMultiply(1);

    tick();
    for ( size_t n=0; n<CNT; ++n )
    {
        mat.vecMulAdd(y, z);
        mat.vecMulAdd(x, z);
    }
    double t1 = tock(N_RUN);
    
    tick();
    for ( size_t n=0; n<CNT; ++n )
    {
        mat.vecMulAdd_ALT(x, z);
        mat.vecMulAdd_ALT(y, z);
    }
    double t2 = tock(N_RUN);

    tick();
    for ( size_t n=0; n<CNT; ++n )
    {
        mat.vecMul(y, z);
        mat.vecMul(x, z);
    }
    double t3 = tock(N_RUN);

    printf("\n%-32s ", mat.what().c_str());
    printf("set %9.0f  muladd %9.0f  alt %9.0f  mul %9.0f", ts, t1, t2, t3);
    checkMatrix(mat, x, z);
    fflush(stdout);
}

//------------------------------------------------------------------------------
#pragma mark - Test Parallel Matrix

#ifdef _OPENMP
#include <omp.h>

/// chunck for parallel execution
constexpr size_t PARALLEL_CHUNK = DIM*4;


template <typename MATRIX>
void checkMatrixParallel(MATRIX & mat, real const* x, real const* y, real * z)
{
    size_t S = mat.size();
    zero_real(S, z);
    // check processing columns one-by-one:
    for ( size_t i = 0; i < S; ++i )
        mat.vecMulAdd(x, z, i, i+1);
    real sum = checksum(S, z);
    printf(" check %+16.6f", sum);

    for ( int i = 0; i < 5; ++i )
    {
        zero_real(S, z);
        omp_set_num_threads(1<<i);
        #pragma omp parallel for
        for ( size_t u = 0; u < S; u += PARALLEL_CHUNK )
            mat.vecMulAdd(x, z, u, u+PARALLEL_CHUNK);
        real sum1 = checksum(S, z);
        if ( sum != sum1 )
            printf(" %luT %+16.6f", 1<<i, sum1);
    }
}


template <typename MATRIX>
void testMatrixParallel(MATRIX & mat, real const* x, real const* y, real * z)
{
    size_t S = mat.size();
    mat.reset();
    for ( size_t n=0; n<icnt_; ++n )
        fillMatrix(mat, iny_[n], inx_[n]);
    mat.prepareForMultiply(1);

    const size_t sup = 7;
    double t[sup] = { 0 };
    
    for ( int i = 0; i < sup; ++i )
    {
        omp_set_num_threads(1<<i);
        tick();
        for ( size_t j=0; j < CNT; ++j )
        {
            #pragma omp parallel for
            for ( size_t i = 0; i < S; i += PARALLEL_CHUNK )
                mat.vecMulAdd(y, z, i, i+PARALLEL_CHUNK);
            #pragma omp parallel for
            for ( size_t i = 0; i < S; i += PARALLEL_CHUNK )
                mat.vecMulAdd(x, z, i, i+PARALLEL_CHUNK);
        }
        t[i] = tock(N_RUN);
    }

    printf("\n%-32s", mat.what().c_str());
    for ( int i = 0; i < sup; ++i )
        printf("  %luT %6.3f", 1<<i, t[i]);
    fflush(stdout);
}

#endif

//------------------------------------------------------------------------------
#pragma mark - Test Iso Matrix

template <typename MATRIX>
void fillMatrixIso(MATRIX& mat, const size_t i, const size_t j)
{
    //printf("fillMatrixIso %lu %lu <---- %f\n", i, j, alpha);
    mat.diagonal(i) += alpha;
    if ( i > j )
    {
        mat.diagonal(j) += alpha;
        mat(i, j) += beta;
    }
    else if ( j > i )
    {
        mat.diagonal(j) += alpha;
        mat(j, i) += beta;
    }
}

#if ( DIM == 1 )
#   define VECMULADDISO vecMulAdd
#elif ( DIM == 2 )
#   define VECMULADDISO vecMulAddIso2D
#elif ( DIM == 3 )
#   define VECMULADDISO vecMulAddIso3D
#endif

template <typename MATRIX>
void checkIsoMatrix(MATRIX & mat, real const* x, real const* y, real * z)
{
    size_t SD = DIM * mat.size();
    zero_real(SD, z);
    mat.VECMULADDISO(x, z);
    mat.VECMULADDISO(y, z);
    printf("   check %+16.6f ", checksum(SD, z));
    
    if ( 1 )
    {
        size_t S = mat.size();
        // compute checksum after manually copying to every subspace:
        MATRIX ful;
        ful.resize(SD);
        for ( size_t j = 0; j < S; ++j )
        {
            real e = mat.diagonal(j);
            for ( int d = 0; d < DIM; ++d )
                ful.diagonal(DIM*j+d) = e;
            for ( size_t i = j+1; i < S; ++i )
            {
                real * p = mat.addr(i, j);
                if ( p )
                {
                    for ( int d = 0; d < DIM; ++d )
                        ful(DIM*i+d, DIM*j+d) = *p;
                }
            }
        }
        ful.prepareForMultiply(1);
        zero_real(SD, z);
        ful.vecMulAdd(x, z);
        ful.vecMulAdd(y, z);
        printf("  %+16.6f ", checksum(SD, z));
    }
}


/// multidimensional isotropic multiplication
template <typename MATRIX>
void testIsoMatrix(MATRIX & mat, real const* x, real const* y, real * z)
{
    tick();
    for ( size_t i=0; i<N_RUN; ++i )
    {
        mat.reset();
        for ( size_t n=0; n<icnt_; ++n )
            fillMatrixIso(mat, inx_[n], iny_[n]);
    }
    double ts = tock(N_RUN);
    
    tick();
    for ( size_t i=0; i<N_RUN; ++i )
    {
        mat.prepareForMultiply(DIM);
        for ( size_t n=0; n<N_MUL; ++n )
        {
            mat.VECMULADDISO(x, z);
            mat.VECMULADDISO(y, z);
        }
    }
    double tm = tock(N_RUN);
    
    printf("\n%-29s ", mat.what().c_str());
    printf("isoset %9.0f  isomul %9.0f", ts, tm);
    checkIsoMatrix(mat, x, y, z);
    fflush(stdout);
}


//------------------------------------------------------------------------------
#pragma mark - Block matrices

const real dir[4] = {  2, 1, -1, 3 };
const real vec[4] = { -1, 3,  1, 2 };

template < typename MATRIX >
void fillBlockMatrix(MATRIX& mat)
{
    typename MATRIX::Block S = MATRIX::Block::outerProduct(dir);
    typename MATRIX::Block U = MATRIX::Block::outerProduct(dir, vec);
    
    for ( size_t n=0; n<icnt_; ++n )
    {
        size_t ii = inx_[n] - inx_[n] % MATRIX::Block::dimension();
        size_t jj = iny_[n] - iny_[n] % MATRIX::Block::dimension();
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
 void testBlockMatrix(SparMatSymBlk & mat, real const* x, real const* y, real * z)
{
    const size_t S = mat.size();
    
    tick();
    for ( size_t r = 0; r < N_RUN; ++r )
    {
        mat.reset();
        fillBlockMatrix(mat);
    }
    double ts = tock(N_RUN);

    mat.prepareForMultiply(1);
    
    // calculate and compare sum for two methods:
    copy_real(S+PAD, y, z);
    mat.vecMulAdd_ALT(x, z);
    real sum = checksum(S, z);
    
    copy_real(S+PAD, y, z);
    mat.vecMulAdd(x, z);
    real res = checksum(S, z);

    tick();
    for ( size_t n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( size_t m=0; m<N_MUL; ++m )
            mat.vecMulAdd_ALT(x, z);
    }
    double nop = N_RUN * mat.nbElements();
    double t1 = tock(nop);

    tick();
    for ( size_t n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( size_t m = 0; m < N_MUL; ++m )
            mat.vecMulAdd(x, z);
    }
    double t2 = tock(nop);
    
    printf("%6lu %26s ", S, mat.what().c_str());
    printf("set %9.1f mul %9.1f  alt %9.1f", ts, t1, t2);
    printf(" :  checksum  %+24.16f %+24.16f", sum, res);
    if ( sum != res )
        printf("  failed!\n");
    else
        printf("\n");
}


//------------------------------------------------------------------------------
#pragma mark - Test Sparse Matrix

void testIsoMatrices(const size_t S, real const* x, real const* y, real * z)
{
    SparMat1 mat1; mat1.resize(S);
    SparMat2 mat2; mat2.resize(S);
    testIsoMatrix(mat1, x, y, z);
    testIsoMatrix(mat2, x, y, z);
    testIsoMatrix(mat1, y, x, z);
    testIsoMatrix(mat2, y, x, z);
}


void testMatrices(const size_t S, real const* x, real const* y, real * z)
{
    SparMat1 mat1; mat1.resize(S); testMatrix(mat1, x, y, z);
    SparMat2 mat2; mat2.resize(S); testMatrix(mat2, x, y, z);

    SparMatB mat3; mat3.resize(S); testMatrix(mat3, x, y, z);
    SparMatD mat4; mat4.resize(S); testMatrix(mat4, x, y, z);
    SparMatA mat5; mat5.resize(S); testMatrix(mat5, x, y, z);
#ifdef _OPENMP
    testMatrixParallel(mat5, x, y, z);
    checkMatrixParallel(mat5, x, y, z);
#endif
    
#if ( 0 )
    std::ofstream os1("mat1.txt");
    std::ofstream os3("mat3.txt");
    mat1.printSparse(os1, 0);
    mat3.printSparse(os3, 0);
#endif
}


void testMatrices(const size_t S, const size_t F)
{
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;
    setVectors(DIM*S, x, y, z);
    //beta = -RNG.preal();

    printf("------ iso %iD x %lu  filled %.1f %% :", DIM, S, F*100.0/S/S);
    setIndices(F, S, 1);
    testIsoMatrices(S, x, y, z);
    printf("\n");

    printf("------ %iD x %lu  filled %.1f %% :", DIM, S, F*100.0/S/S);
    setIndices(F, DIM*S, DIM);
    testMatrices(DIM*S, x, y, z);
    printf("\n");
    
    setIndices(0, S, DIM);
    free_real(x);
    free_real(y);
    free_real(z);
}


void testBlockMatrix(const size_t S, const size_t F)
{
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;

    setIndices(F, S, DIM);
    setVectors(S, x, y, z);
    //beta = -RNG.preal();
    
    SparMatSymBlk mat;
    mat.resize(S);
    testBlockMatrix(mat, x, y, z);
    
    setIndices(0, S, DIM);
    free_real(x);
    free_real(y);
    free_real(z);
}

//------------------------------------------------------------------------------
#pragma mark - Main

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
    testMatrices(DIM, 1);
    testMatrices(2, 1);
    testMatrices(3, 3);
    testMatrices(4, 4);
    testMatrices(7, 5);
    testMatrices(17, 7);
    testMatrices(33, 1111);
#endif
#if ( 0 )
    printf("\ntest_matrix BLOCK_SIZE %i (%s)", DIM, SparMatSymB::Block::what().c_str());
    size_t siz = DIM;
    for ( int i = 0; i < 14; ++i )
    {
        siz = 1 + ( 2 << i ) * RNG.preal();
        size_t fill = 1 + siz * siz * ( 0.01 * RNG.preal() );
        testBlockMatrix(siz, fill);
    }
#endif
#if ( 0 )
    testBlockMatrix(DIM*17, 2);
    testBlockMatrix(DIM*347, 1019);
    testBlockMatrix(DIM*753, 43039);
#endif
#if ( 0 )
    testBlockMatrix(2253, 1<<14);
    testBlockMatrix(2253, 1<<14);
    testBlockMatrix(2253, 1<<14);
    testBlockMatrix(2253, 1<<14);
#endif
#if ( 0 )
    testBlockMatrix(DIM*1251, 25821);
    testBlockMatrix(DIM*1785, 153034);
    testBlockMatrix(DIM*2311, 231111);
    //testBlockMatrix(DIM*3217, 671234);
#endif
#if ( 0 )
    //testMatrices(17, 23);
    //testMatrices(91, 1<<12);
    testMatrices(196, 1<<13);
    //testMatrices(436, 1<<15);
    testMatrices(8*94, 1<<16);
    //testMatrices(8*169, 1<<17);
    testMatrices(8*331, 1<<18);
#endif
    testMatrices(8*111, 1<<17);
#if ( 0 )
    //testMatrices(17, 23);
    size_t dim[5] = { 0 };
    for ( int i = 0; i < 5; ++i ) dim[i] = RNG.pint32(1<<(i+7));
    qsort(dim, 5, sizeof(size_t), compareInt);
    for ( int i = 0; i < 5; ++i )
        testMatrices(dim[i], RNG.pint32(dim[i]*dim[i]));
#endif
}

