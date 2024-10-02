// Cytosim was created by Francois Nedelec.  Copyright 2022 Cambridge University.

#include <sys/time.h>
#include <fstream>
#include <map>

#include "assert_macro.h"
#include "exceptions.h"
#include "vecprint.h"
#include "random.h"
#include "timer.h"
#include "real.h"
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
constexpr size_t N_MUL = 97;

// number of repeat of ( 1 prepare + N_MUL multiplications )
constexpr size_t N_RUN = 8;

// extra for allocation
#define PAD 4

size_t icnt_ = 0;
unsigned * inx_ = nullptr;
unsigned * iny_ = nullptr;

real alpha = 2.0;
real kappa = 0.7;
real beta = -1.0;
real iota = -0.3;

constexpr size_t REP = N_RUN * N_MUL;

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

///set indices within [0, sup] that are below the diagonal: Y >= X
void setIndices(unsigned sup, size_t cnt)
{
    delete[] inx_;
    delete[] iny_;
    inx_ = nullptr;
    iny_ = nullptr;
    icnt_ = cnt;
    if ( cnt > 0 )
    {
        inx_ = new unsigned[icnt_];
        iny_ = new unsigned[icnt_];
        for ( size_t n = 0; n < cnt; ++n )
        {
            unsigned i = RNG.pint32(sup);
            unsigned j = RNG.pint32(sup);
            inx_[n] = std::min(i,j);
            iny_[n] = std::max(i,j);
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
            z[n] = RNG.preal();
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
void compareMatrices(unsigned S, MATRIX & mat1, MATROX& mat2, size_t fill)
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
        unsigned ii = RNG.pint32(S);
        unsigned jj = RNG.pint32(S);
        if ( ii != jj )
        {
            unsigned i = std::max(ii, jj);
            unsigned j = std::min(ii, jj);
            
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

        mat1.addDiagonalBlock(tmp1, S, inx, cnt, 1, 1);
        mat2.addDiagonalBlock(tmp2, S, inx, cnt, 1, 1);
        
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
            mat1.addDiagonalBlock(tmp2, S, inx, cnt, 1, 1);
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
void fillMatrix(MATRIX& mat)
{
#if ( DIM == 3 )
    Matrix33 S(alpha, beta, -beta, beta, -alpha, beta, -beta, beta, alpha);
    Matrix33 U(kappa, iota, iota, -iota, kappa, iota, -iota, -iota, kappa);
    //Matrix33 U(kappa, iota, -iota, iota, kappa, iota, -iota, iota, -kappa);
#elif ( DIM == 2 )
    Matrix22 S(alpha, beta, -beta, alpha);
    Matrix22 U(kappa, iota,  iota, kappa);
#else
    Matrix11 S(alpha);
    Matrix11 U(kappa);
#endif
    for ( size_t n = 0; n < icnt_; ++n )
    {
        size_t i = DIM * iny_[n];
        size_t j = DIM * inx_[n];
        //printf("\n  ( %lu %lu ) ", i, j); S.print_smart(std::cout);
        // this is a block on the diagonal:
        for ( size_t x = 0; x < 1; ++x )
        for ( size_t y = x; y < DIM; ++y )
            mat(i+y, i+x) += S(y, x);
        
        if ( i != j )
        {
            // this is a block on the diagonal:
            for ( size_t x = 0; x < DIM; ++x )
            for ( size_t y = x; y < DIM; ++y )
                mat(j+y, j+x) += S(y,x);
            // off-diagonal non-symmetric block!
            for ( size_t x = 0; x < DIM; ++x )
            for ( size_t y = x; y < DIM; ++y )
                mat(i+y, j+x) += U(y,x);
        }
    }
}

// check the results of Z = Mat * X + Y with different methods
template <typename MATRIX>
real checkMatrix(MATRIX & mat, real const* x, real const* y, real * z, size_t NP)
{
    size_t N = mat.size();
    if ( NP > 1 ) VecPrint::print("xxxx", NP, x);
    if ( NP > 1 ) VecPrint::print("yyyy", NP, y);
    copy_real(N, y, z);
    mat.vecMulAdd(x, z);
    real sum1 = checksum(N, z);
    if ( NP > 1 ) VecPrint::print("vec1", NP, z);

    copy_real(N, y, z);
    mat.vecMulAdd_ALT(x, z);
    real sum2 = checksum(N, z);
    if ( NP > 1 ) VecPrint::print("vec2", NP, z);

    mat.vecMul(x, z);
    for ( size_t i = 0; i < N; ++i ) z[i] += y[i];
    real sum3 = checksum(N, z);
    if ( NP > 1 ) VecPrint::print("vec3", NP, z);
    
    real n = std::max(abs_real(sum1-sum2), abs_real(sum2-sum3));
    if ( NP > 0 )
    {
        printf(" check %+16.6f %+16.6f %+16.6f ", sum1, sum2, sum3);
        if ( n > 1e-6 )
            printf(" FAILED");
        else
            printf(" : pass");
    }
    return n;
}


template <typename MATRIX, void (*fill)(MATRIX& mat)>
void testMatrix(MATRIX & mat, real const* x, real const* y, real * z)
{
    tick();
    for ( size_t n=0; n<N_RUN; ++n )
    {
        mat.reset();
        fill(mat);
        mat.prepareForMultiply(1);
    }
    double ts = tock(N_RUN);

    tick();
    for ( size_t n=0; n<REP; ++n )
    {
        mat.vecMulAdd(y, z);
        mat.vecMulAdd(x, z);
    }
    double t1 = tock(N_RUN);
    
    tick();
    for ( size_t n=0; n<REP; ++n )
    {
        mat.vecMulAdd_ALT(x, z);
        mat.vecMulAdd_ALT(y, z);
    }
    double t2 = tock(N_RUN);

    tick();
    for ( size_t n=0; n<REP; ++n )
    {
        mat.vecMul(y, z);
        mat.vecMul(x, z);
    }
    double t3 = tock(N_RUN);

    printf("\n%-32s ", mat.what().c_str());
    printf("set %9.0f  muladd %9.0f  alt %9.0f  mul %9.0f", ts, t1, t2, t3);
    checkMatrix(mat, x, y, z, 1);
    fflush(stdout);
}

//------------------------------------------------------------------------------
#pragma mark - Multithreaded test with libdispatch

#if defined(__APPLE__)

#include <dispatch/dispatch.h>

struct Job
{
    SparMatBlk const* mat;
    real const* src;
    real * dst;
    size_t mul;
};

void worker(void* context, size_t i)
{
    Job * job = (Job*)context;
    //fprintf(stderr, "%lX ", i&15);
    job->mat->vecMulAdd(job->src, job->dst, job->mul*i, job->mul*(i+1));
}

void testMatrixDispatch(SparMatBlk& mat, real const* x, real const* y, real * z, std::map<int, double>& rec)
{
    size_t S = mat.size() / DIM;

    dispatch_queue_t queue = dispatch_queue_create("QUEUE", DISPATCH_QUEUE_CONCURRENT);
    for ( unsigned u : { 1, 2, 4, 8, 10, 12, 16, 32 } )
    {
        size_t cnt = (S+u-1)/u;
        Job job1 { &mat, x, z, u };
        Job job2 { &mat, y, z, u };
        zero_real(DIM*S, z);
        dispatch_apply_f(cnt, queue, &job1, worker);
        dispatch_apply_f(cnt, queue, &job2, worker);
        printf("\n%2u libdispatch %+16.6f", u, checksum(DIM*S, z));
        tick();
        for ( size_t j = 1; j < REP; ++j )
        {
            dispatch_apply_f(cnt, queue, &job1, worker);
            dispatch_apply_f(cnt, queue, &job2, worker);
        }
        rec[u] = tock(N_RUN);
    }
    dispatch_release(queue);
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Multithreaded test with OpenMP

#ifdef _OPENMP
#include <omp.h>

template <typename MATRIX>
void checkMatrixOMP(MATRIX & mat, real const* x, real const* y, real * z, size_t chunk)
{
    size_t S = mat.size() / DIM;
    zero_real(DIM*S, z);
    // check processing columns one-by-one:
    for ( size_t i = 0; i < S; ++i )
    {
        mat.vecMulAdd(x, z, i, i+1);
        mat.vecMulAdd(y, z, i, i+1);
    }
    printf("\n%s OpenMP check %+16.6f", mat.what().c_str(), checksum(DIM*S, z));

    // change number of parallel threads:
    for ( int i = 0; i < 4; ++i )
    {
        zero_real(DIM*S, z);
        omp_set_num_threads(1<<i);
        #pragma omp parallel for
        for ( size_t u = 0; u < S; u += chunk )
        {
            mat.vecMulAdd(x, z, u, u+chunk);
            mat.vecMulAdd(y, z, u, u+chunk);
        }
        printf(" %+16.6f", checksum(DIM*S, z));
    }
}


template <typename MATRIX>
void testMatrixOMP(MATRIX & mat, real const* x, real const* y, real * z, size_t chunk, std::map<int, double>& rec)
{
    size_t S = mat.size() / DIM;
    
    for ( int u : { 1, 2, 3, 4, 6, 8, 10, 16 } )
    {
        zero_real(DIM*S, z);
        omp_set_num_threads(u);
        tick();
        for ( size_t j=0; j < REP; ++j )
        {
            #pragma omp parallel for
            for ( size_t i = 0; i < S; i += chunk )
            {
                mat.vecMulAdd(x, z, i, i+chunk);
                mat.vecMulAdd(y, z, i, i+chunk);
            }
        }
        rec[u] = tock(N_RUN);
    }
}

#endif


void testParallelVecmul(const unsigned S, const size_t F)
{
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;
    setVectors(DIM*S, x, y, z);
    setIndices(S, F);

    SparMatA mat;
    mat.resize(DIM*S);
    mat.reset();
    fillMatrix(mat);
    mat.prepareForMultiply(1);

    printf("------ %i x %i  filled %.2f %%: %s", DIM, S, F*100.0/S/S, mat.what().c_str());
    
    zero_real(DIM*S, z);
    mat.vecMulAdd(x, z);
    mat.vecMulAdd(y, z);
    printf("\nsequential     %+16.6f ", checksum(DIM*S, z));
    tick();
    for ( size_t j = 1; j < REP; ++j )
    {
        mat.vecMulAdd(x, z);
        mat.vecMulAdd(y, z);
    }
    
    std::map<int, double> rec;
    rec[0] = tock(N_RUN);
    printf(" %-6.1f", rec[0]);

#ifdef _OPENMP
    size_t chunk = S / 16;
    checkMatrixOMP(mat, x, y, z, chunk);
    testMatrixOMP(mat, x, y, z, chunk, rec);
    printf("\n-->  OpenMP  ");
    for ( auto u : rec )
        printf(" %uT %-6.1f", u.first, u.second);
#endif
#if defined(__APPLE__)
    testMatrixDispatch(mat, x, y, z, rec);
    printf("\n--> dispatch ");
    for ( auto u : rec )
        printf(" %uT %-6.1f", u.first, u.second);
#endif
    fflush(stdout);
    printf("\n");

    setIndices(0, 0);
    setVectors(0, x, y, z);
}


//------------------------------------------------------------------------------
#pragma mark - Test Iso Matrix

template <typename MATRIX>
void fillMatrixIso(MATRIX& mat)
{
    mat.reset();
    for ( size_t n = 0; n < icnt_; ++n )
    {
        size_t i = iny_[n];
        size_t j = inx_[n];
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
        for ( index_t j = 0; j < S; ++j )
        {
            real e = mat.diagonal(j);
            for ( int d = 0; d < DIM; ++d )
                ful.diagonal(DIM*j+d) = e;
            for ( index_t i = j+1; i < S; ++i )
            {
                real * p = mat.address(i, j);
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
        fillMatrixIso(mat);
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
#if ( 0 )
    SparMat1 mat1; mat1.resize(S); testMatrix<SparMat1, fillMatrix>(mat1, x, y, z);
    SparMat2 mat2; mat2.resize(S); testMatrix<SparMat2, fillMatrix>(mat2, x, y, z);
#endif
    SparMatB mat3; mat3.resize(S); testMatrix<SparMatB, fillMatrix>(mat3, x, y, z);
    SparMatD mat4; mat4.resize(S); testMatrix<SparMatD, fillMatrix>(mat4, x, y, z);
    //SparMatA mat5; mat5.resize(S); testMatrix<SparMatA, fillMatrix>(mat5, x, y, z);
#if ( 0 )
    std::ofstream os1("mat1.txt");
    std::ofstream os3("mat3.txt");
    mat1.printSparse(os1, 0);
    mat3.printSparse(os3, 0);
#endif
}


template < typename MATRIX >
void testMatricesMulti(const unsigned S, const size_t F, const size_t rep)
{
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;
    setVectors(DIM*S, x, y, z);
    
    MATRIX mat;
    printf("------ %i x %i  multi-test %s:\n", DIM, S, mat.what().c_str());
    mat.resize(DIM*S);
    real err = 0;
    for ( size_t i = 0; i < rep; ++i )
    {
        setIndices(S, F);
        mat.reset();
        fillMatrix(mat);
        mat.prepareForMultiply(1);
        err = checkMatrix(mat, x, y, z, 0);
        if ( err > REAL_EPSILON ) break;
    }
    checkMatrix(mat, x, y, z, 1+15*(err>1e-6));
    printf(" ---> error %e\n", err);
    //mat.printSparse(std::cout, 0, 0, 1000);
    setIndices(0, 0);
    setVectors(0, x, y, z);
}


void testMatrices(const unsigned S, const size_t F)
{
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;
    setVectors(DIM*S, x, y, z);
    //beta = -RNG.preal();

    if ( 0 )
    {
        printf("------ iso %iD x %i  filled %.2f %%", DIM, S, F*100.0/S/S);
        setIndices(S, F);
        testIsoMatrices(S, x, y, z);
        printf("\n");
    }
    printf("------ %i x %i  filled %.2f %%:", DIM, S, F*100.0/S/S);
    setIndices(S, F);
    testMatrices(DIM*S, x, y, z);
    printf("\n");
    setIndices(0, 0);
    setVectors(0, x, y, z);
}


//------------------------------------------------------------------------------
#pragma mark - Block matrices

const real dir[4] = {  2, 1, -1, 3 };
const real vec[4] = { -1, 3,  1, 2 };

template < typename MATRIX >
void fillBlockMatrix(MATRIX& mat)
{
    typename MATRIX::Block D = MATRIX::Block::outerProduct(dir);
    typename MATRIX::Block B = MATRIX::Block::outerProduct(dir);
    typename MATRIX::Block U = MATRIX::Block::outerProduct(vec, dir);

    for ( size_t n = 0; n < icnt_; ++n )
    {
        size_t i = iny_[n];
        size_t j = inx_[n];
        mat.diag_block(i).add_half(D);
        if ( i > j )
        {
            mat.diag_block(j).add_half(D);
            mat.block(j+1, j).add_full(B);
            mat.block(i, j).sub_full(U);
        }
    }
}


void testBlockMatrix(const unsigned S, const size_t F)
{
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;
    setVectors(DIM*S, x, y, z);
    setIndices(S, F);
    
    printf("------ %i x %u with %lu blocks:", DIM, S, F);
    SparMatA A; A.resize(DIM*S); testMatrix<SparMatA, fillBlockMatrix>(A, x, y, z);
    SparMatB B; B.resize(DIM*S); testMatrix<SparMatB, fillBlockMatrix>(B, x, y, z);
    SparMatD D; D.resize(DIM*S); testMatrix<SparMatD, fillBlockMatrix>(D, x, y, z);
    printf("\n");
    
    //D.printSummary(std::cout, 0, DIM*S);

    setIndices(0, 0);
    setVectors(0, x, y, z);
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
    compareMatrices(4*3, mat1, mat2, 1<<4);
    compareMatrices(4*5, mat1, mat3, 1<<4);
    compareMatrices(4*7, mat2, mat3, 1<<5);
#endif
#if ( 0 )
    compareMatrices(4*11, mat1, mat3, 1<<6);
    compareMatrices(4*33, mat1, mat3, 1<<16);
    compareMatrices(4*3, mat1, mat4, 1<<16);
    compareMatrices(4*7, mat1, mat4, 1<<16);
    compareMatrices(4*11, mat1, mat4, 1<<16);
    compareMatrices(4*33, mat1, mat4, 1<<16);
#endif
#if ( 0 )
    testMatrices(1, 1);
    testMatrices(2, 1);
    testMatrices(3, 3);
    testMatrices(4, 4);
    testMatrices(7, 5);
    testMatrices(17, 7);
    testMatrices(33, 21);
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
    testBlockMatrix(17, 2);
    testBlockMatrix(347, 1019);
    testBlockMatrix(753, 43039);
#endif
#if ( 0 )
    testBlockMatrix(2253, 1<<14);
    testBlockMatrix(2253, 1<<14);
    testBlockMatrix(2253, 1<<14);
    testBlockMatrix(2253, 1<<14);
#endif
#if ( 0 )
    testBlockMatrix(1251, 25821);
    testBlockMatrix(1785, 153034);
    testBlockMatrix(2311, 231111);
    //testBlockMatrix(DIM*3217, 671234);
#endif
#if ( 0 )
    testMatrices(3, 1);
    //testMatrices(7, 23);
    //testMatrices(17, 23);
    //testMatrices(29, 1<<12);
    testMatrices(65, 1<<13);
    //testMatrices(145, 1<<15);
    testMatrices(8*31, 1<<16);
    //testMatrices(8*56, 1<<17);
    testMatrices(8*110, 1<<18);
#else
    testMatricesMulti<SparMatD>(37, 1<<17, 1024);
    testMatrices(5494, 131836);
    testBlockMatrix(5494, 131836);
    testBlockMatrix(5494, 431836);
    //testParallelVecmul(15464, 131836);
#endif
#if ( 0 )
    size_t dim[5] = { 0 };
    for ( int i = 0; i < 5; ++i ) dim[i] = RNG.pint32(1<<(i+7));
    qsort(dim, 5, sizeof(size_t), compareInt);
    for ( int i = 0; i < 5; ++i )
        testMatrices(dim[i], RNG.pint32(dim[i]*dim[i]));
#endif
}

