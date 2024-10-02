// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include <sstream>

#include "sparmatsymblk.h"
#include "assert_macro.h"
#include "vector2.h"
#include "vector3.h"
#include "simd.h"

// Flags to enable SIMD implementation
#if defined(__AVX__)
#  include "simd_float.h"
#  define SMSB_USES_AVX 1
#  define SMSB_USES_SSE 1
#elif USE_SIMD
#  include "simd_float.h"
#  define SMSB_USES_AVX 0
#  define SMSB_USES_SSE 1
#else
#  define SMSB_USES_AVX 0
#  define SMSB_USES_SSE 0
#endif


SparMatSymBlk::SparMatSymBlk()
: rsize_(0), alloc_(0), column_(nullptr)
{
    colidx_ = new unsigned[2]();
}


void SparMatSymBlk::allocate(index_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be increased to gain performance:
          more memory will be used, but reallocation will be less frequent
        */
        constexpr index_t chunk = 32;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "SMSB allocates %u\n", alc);
        Column * ptr = new Column[alc];
       
        if ( column_ )
        {
            for (index_t n = 0; n < alloc_; ++n )
                ptr[n] = column_[n];
            delete[] column_;
        }
        
        column_ = ptr;
        alloc_  = alc;
        
        delete[] colidx_;
        colidx_ = new unsigned[alc+2];
        for ( unsigned n = 0; n <= alc; ++n )
            colidx_[n] = n;
    }
}


void SparMatSymBlk::deallocate()
{
    delete[] column_;
    delete[] colidx_;
    column_ = nullptr;
    colidx_ = nullptr;
    alloc_ = 0;
}

//------------------------------------------------------------------------------
#pragma mark - Column

SparMatSymBlk::Column::Column()
{
    nbb_ = 0;
    alo_ = 0;
    inx_ = nullptr;
    blk_ = nullptr;
}


/*
 This may require some smart allocation scheme.
 */
void SparMatSymBlk::Column::allocate(index_t alc)
{
    if ( alc > alo_ )
    {
        //if ( inx_ ) fprintf(stderr, "SMSB reallocates column %lu for %lu\n", inx_[0], alc);
        //else fprintf(stderr, "SMSB allocates column for %u\n", alc);
        /*
         'chunk' can be increased, to possibly gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr index_t chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk - 1 );
        
        // use aligned memory:
        void * ptr = new_real(alc*sizeof(Block)/sizeof(real)+4);
        Block * blk_new  = new(ptr) Block[alc];

        if ( posix_memalign(&ptr, 32, alc*sizeof(unsigned)) )
            throw std::bad_alloc();
        auto * inx_new = (unsigned*)ptr;

        if ( inx_ )
        {
            for ( index_t n = 0; n < nbb_; ++n )
                inx_new[n] = inx_[n];
            free(inx_);
        }

        if ( blk_ )
        {
            for ( index_t n = 0; n < nbb_; ++n )
                blk_new[n] = blk_[n];
            free_real(blk_);
        }
        inx_ = inx_new;
        blk_ = blk_new;
        alo_ = alc;
        
        //std::clog << "Column " << this << "  " << alc << ": ";
        //std::clog << " alignment " << ((uintptr_t)elem_ & 63) << "\n";
    }
}


void SparMatSymBlk::Column::deallocate()
{
    //if ( inx_ ) fprintf(stderr, "SMSB deallocates column %lu of size %lu\n", inx_[0], alo_);
    free(inx_);
    free_real(blk_);
    inx_ = nullptr;
    blk_ = nullptr;
    alo_ = 0;
    nbb_ = 0;
}


void SparMatSymBlk::Column::operator = (SparMatSymBlk::Column & col)
{
    //if ( inx_ ) fprintf(stderr, "SMSB transfers column %u\n", inx_[0]);
    free(inx_);
    free_real(blk_);

    nbb_ = col.nbb_;
    alo_ = col.alo_;
    inx_ = col.inx_;
    blk_ = col.blk_;
    
    col.nbb_ = 0;
    col.alo_ = 0;
    col.inx_ = nullptr;
    col.blk_ = nullptr;
}




/* This is a silly search that could be optimized */
SparMatSymBlk::Block* SparMatSymBlk::Column::find_block(index_t ii) const
{
    for ( index_t n = 0; n < nbb_; ++n )
        if ( inx_[n] == ii )
            return blk_ + n;
    return nullptr;
}

/**
 This allocates to be able to hold the matrix element if necessary
 */
SparMatSymBlk::Block& SparMatSymBlk::Column::block(index_t ii, index_t jj)
{
    assert_true( ii >= jj );
    SparMatSymBlk::Block * B = find_block(ii);
    if ( !B )
    {
        allocate(nbb_+2);
        if ( nbb_ == 0 )
        {
            // put diagonal term always first:
            inx_[0] = jj;
            blk_[0].reset();
            nbb_ = 1;
            if ( ii == jj )
                return blk_[0];
        }
        assert_true( ii > jj );
        //add the requested term last:
        B = blk_ + nbb_;
        inx_[nbb_] = ii;
        B->reset();
        ++nbb_;
    }
    //printColumn(jj);
    return *B;
}


void SparMatSymBlk::Column::reset()
{
    nbb_ = 0;
}

SparMatSymBlk::Block& SparMatSymBlk::diag_block(index_t ii)
{
    assert_true( ii < rsize_ );
    Column & col = column_[ii];
    if ( col.nbb_ == 0 )
    {
        //fprintf(stderr, "new diagonal element for column %i\n", ii);
        col.allocate(1);
        col.nbb_ = 1;
        // put diagonal term always first:
        col.inx_[0] = ii;
        col.blk_[0].reset();
    }
    assert_true(col.inx_[0] == ii);
    return col.blk_[0];
}


real& SparMatSymBlk::element(index_t iii, index_t jjj)
{
    // branchless code to address lower triangle
    index_t ii = std::max(iii, jjj);
    index_t jj = std::min(iii, jjj);
#if ( S_BLOCK_SIZE == 1 )
    return column_[jj].block(ii, jj).value();
#else
    index_t i = ii / S_BLOCK_SIZE;
    index_t j = jj / S_BLOCK_SIZE;
    return column_[j].block(i, j)(ii%S_BLOCK_SIZE, jj%S_BLOCK_SIZE);
#endif
}


real* SparMatSymBlk::address(index_t iii, index_t jjj) const
{
    // branchless code to address lower triangle
    index_t ii = std::max(iii, jjj);
    index_t jj = std::min(iii, jjj);
#if ( S_BLOCK_SIZE == 1 )
    return column_[jj].block(ii, jj).data();
#else
    index_t i = ii / S_BLOCK_SIZE;
    index_t j = jj / S_BLOCK_SIZE;
    Block * B = column_[j].find_block(i);
    if ( B )
        return B->addr(ii%S_BLOCK_SIZE, jj%S_BLOCK_SIZE);
    return nullptr;
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

void SparMatSymBlk::reset()
{
    for ( index_t n = 0; n < rsize_; ++n )
        column_[n].reset();
}


bool SparMatSymBlk::notZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Column & col = column_[jj];
        for ( index_t n = 0 ; n < col.nbb_ ; ++n )
            if ( col[n] != 0.0 )
                return true;
    }
    //if here, the matrix is empty
    return false;
}


void SparMatSymBlk::scale(const real alpha)
{
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Column & col = column_[jj];
        for ( index_t n = 0 ; n < col.nbb_ ; ++n )
            col[n].scale(alpha);
    }
}


void SparMatSymBlk::addDiagonalBlock(real* mat, index_t ldd, const index_t start, const index_t cnt, index_t mul) const
{
    assert_true( mul == S_BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );
    
    for ( index_t jj = 0; jj < cnt; ++jj )
    {
        Column & col = column_[jj+start];
        if ( col.nbb_ > 0 )
        {
            assert_true(col.inx_[0] == jj + start);
            real * dst = mat + ( jj + ldd * jj ) * S_BLOCK_SIZE;
            col[0].addto_symm(dst, ldd);
            for ( index_t n = 1; n < col.nbb_; ++n )
            {
                // assuming lower triangle is stored:
                assert_true( col.inx_[n] > jj + start );
                auto ij = col.inx_[n] - ( jj + start );
                if ( jj+ij < cnt )
                {
                    //fprintf(stderr, "SMSBD %4lu %4lu\n", ii, jj); col[n].print(stderr);
                    col[n].addto(dst+ij*S_BLOCK_SIZE, ldd);
                    col[n].addto_trans(dst+ij*ldd*S_BLOCK_SIZE, ldd);
                }
            }
        }
    }
}


void SparMatSymBlk::addLowerBand(real alpha, real* mat, index_t ldd, index_t start, index_t cnt,
                                 const index_t mul, const index_t rank) const
{
    assert_true( mul == S_BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );

    for ( index_t jj = 0; jj < cnt; ++jj )
    {
        Column & col = column_[jj+start];
        if ( col.nbb_ > 0 )
        {
            assert_true(col.inx_[0] == jj + start);
            index_t sup = std::min(cnt-jj, rank+1);
            real * dst = mat + ( jj + ldd * jj ) * S_BLOCK_SIZE;
            col[0].addto_lower(dst, ldd, alpha);
            for ( index_t n = 1; n < col.nbb_; ++n )
            {
                // assuming lower triangle is stored:
                assert_true(col.inx_[n] > jj + start);
                auto ij = col.inx_[n] - ( jj + start );
                if ( ij < sup )
                {
                    //fprintf(stderr, "SMSB %4i %4i % .4f\n", ii, jj, a);
                    col[n].addto(dst+ij*S_BLOCK_SIZE, ldd, alpha);
                    //col[n].addto_trans(dst+ij*ldd*S_BLOCK_SIZE, ldd, alpha);
                } else break; // assuming the column is sorted
            }
        }
    }
}


void SparMatSymBlk::addDiagonalTrace(real alpha, real* mat, index_t ldd,
                                     const index_t start, const index_t cnt,
                                     const index_t mul, const index_t rank, const bool sym) const
{
    assert_true( mul == S_BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );

    for ( index_t jj = 0; jj < cnt; ++jj )
    {
        Column & col = column_[jj+start];
        if ( col.nbb_ > 0 )
        {
            assert_true(col.inx_[0] == jj + start);
            index_t sup = std::min(cnt-jj, rank+1);
            // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
            real * dst = mat + ( jj + ldd * jj );
            dst[0] += alpha * col[0].trace();  // diagonal term
            for ( index_t n = 1; n < col.nbb_; ++n )
            {
                // assuming lower triangle is stored:
                assert_true(col.inx_[n] > jj + start);
                auto ij = col.inx_[n] - ( jj + start );
                if ( ij < sup )
                {
                    real a = alpha * col[n].trace();
                    //fprintf(stderr, "SMSB %4lu %4lu : %.4f\n", i, j, a);
                    dst[ij] += a;
                    if ( sym ) dst[ldd*ij] += a;
                } else break; // assuming the column is sorted
            }
        }
    }
}



int SparMatSymBlk::bad() const
{
    if ( rsize_ <= 0 ) return 1;
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Column & col = column_[jj];
        for ( index_t n = 0 ; n < col.nbb_ ; ++n )
        {
            if ( col.inx_[n] >= rsize_ ) return 2;
            if ( col.inx_[n] <= jj )    return 3;
        }
    }
    return 0;
}


/** all allocated elements are counted, even if zero */
size_t SparMatSymBlk::nbElements(index_t start, index_t stop, size_t& alc) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    alc = 0;
    index_t cnt = 0;
    for ( index_t i = start; i < stop; ++i )
    {
        cnt += column_[i].nbb_;
        alc += column_[i].alo_;
    }
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark -


std::string SparMatSymBlk::what() const
{
    size_t alc = 0;
    size_t cnt = nbElements(0, rsize_, alc);
    std::ostringstream msg;
#if SMSB_USES_AVX && REAL_IS_DOUBLE
    msg << "SMSBx ";
#elif SMSB_USES_SSE && REAL_IS_DOUBLE
    msg << "SMSBe ";
#elif SMSB_USES_SSE
    msg << "SMSBf ";
#else
    msg << "SMSB ";
#endif
    msg << Block::what() << "*" << cnt << " (" << alc << ")";
    return msg.str();
}


static void printSparseBlock(std::ostream& os, real epsilon, SparMatSymBlk::Block const& B, index_t ii, index_t jj)
{
    index_t d = ( ii == jj );
    for ( index_t x = 0; x < B.dimension(); ++x )
    for ( index_t y = d*x; y < B.dimension(); ++y )
    {
        real v = B(y, x);
        if ( abs_real(v) > epsilon )
            os << std::setw(6) << ii+y << " " << std::setw(6) << jj+x << " " << std::setw(16) << v << "\n";
    }
}


void SparMatSymBlk::printSparse(std::ostream& os, real epsilon, index_t start, index_t stop) const
{
    os << "% SparMatSymBlk size " << rsize_ << ":\n";
    stop = std::min(stop, rsize_);
    std::streamsize p = os.precision(8);
    if ( ! column_ )
        return;
    for ( index_t jj = start; jj < stop; ++jj )
    {
        Column & col = column_[jj];
        if ( col.notEmpty() )
            os << "% column " << jj << "\n";
        for ( index_t n = 0 ; n < col.nbb_ ; ++n )
            printSparseBlock(os, epsilon, col.blk_[n], col.inx_[n], jj);
    }
    os.precision(p);
}


void SparMatSymBlk::printSummary(std::ostream& os, index_t start, index_t stop)
{
    stop = std::min(stop, rsize_);
    os << "SMSB size " << rsize_ << ":";
    for ( index_t j = start; j < stop; ++j )
        if ( column_[j].notEmpty() )
        {
            os << "\n   " << j << "   " << column_[j].nbb_;
            os << " index " << colidx_[j];
        }
    std::endl(os);
}


void SparMatSymBlk::Column::printBlocks(std::ostream& os) const
{
    for ( index_t n = 0; n < nbb_; ++n )
        os << " " << inx_[n] << " " << blk_[n];
}


void SparMatSymBlk::printBlocks(std::ostream& os) const
{
    for ( index_t j = 0; j < rsize_; ++j )
    {
        os << "\nSMSB  col " << j;
        column_[j].printBlocks(os);
    }
    std::endl(os);
}


//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication


/// A block element of the sparse matrix suitable for qsort()
class SparMatSymBlk::Element
{
public:
    /// block elements
    real blk[S_BLOCK_SIZE*S_BLOCK_SIZE];

    /// index
    index_t inx;
};


/// qsort function comparing line indices
static int compareSMSBElement(const void * A, const void * B)
{
    index_t a = static_cast<SparMatSymBlk::Element const*>(A)->inx;
    index_t b = static_cast<SparMatSymBlk::Element const*>(B)->inx;
    
    return ( a > b ) - ( b > a );
}

/**
 This copies the data to the provided temporary array
 */
void SparMatSymBlk::Column::sortElements(Element tmp[], index_t tmp_size)
{
    assert_true( nbb_ <= tmp_size );
    for ( index_t i = 1; i < nbb_; ++i )
    {
        blk_[i].store(tmp[i].blk);
        tmp[i].inx = inx_[i];
    }
    
    //std::clog << "sizeof(SparMatSymBlk::Element) " << sizeof(Element) << "\n";
    qsort(tmp+1, nbb_-1, sizeof(Element), &compareSMSBElement);
    
    for ( index_t i = 1; i < nbb_; ++i )
    {
        blk_[i].load(tmp[i].blk);
        inx_[i] = tmp[i].inx;
    }
}


index_t SparMatSymBlk::newElements(Element*& ptr, index_t cnt)
{
    constexpr index_t chunk = 16;
    index_t all = ( cnt + chunk - 1 ) & ~( chunk - 1 );
    free(ptr);  // Element has no destructor
    void* tmp = nullptr;
    if ( posix_memalign(&tmp, 32, all*sizeof(Element)) )
        throw std::bad_alloc();
    ptr = new(tmp) Element[all];
    return all;
}


void SparMatSymBlk::sortElements()
{
    //index_t cnt = 0;
    index_t tmp_size = 0;
    Element * tmp = nullptr;
    
    for ( index_t j = colidx_[0]; j < rsize_; j = colidx_[j+1] )
    {
        assert_true( j < rsize_ );
        Column & col = column_[j];
        assert_true( col.nbb_ > 0 );
        //std::clog << "SMSB column " << j << " has " << col.nbb_ << " elements\n";
        
        // order the elements within the column:
        if ( col.nbb_ > 2 )
        {
            if ( tmp_size < col.nbb_ )
                tmp_size = newElements(tmp, col.nbb_);
            col.sortElements(tmp, tmp_size);
        }
        
        //++cnt;
        
        // diagonal element should be first:
        assert_true( col.inx_[0] == j );
        col.blk_[0].copy_lower();
#ifndef NDEBUG
        for ( index_t n = 1 ; n < col.nbb_ ; ++n )
        {
            const index_t i = col.inx_[n];
            assert_true( i < rsize_ );
            assert_true( i > j );
        }
#endif
    }
    
    free(tmp);
    //std::clog << "SparMatSymBlk " << rsize_ << " with " << cnt << " non-empty columns\n";
}


bool SparMatSymBlk::prepareForMultiply(int)
{
    colidx_[rsize_] = rsize_;
    if ( rsize_ > 0 )
    {
        index_t inx = rsize_;
        index_t nxt = rsize_;
        while ( inx-- > 0 )
        {
            if ( column_[inx].notEmpty() )
                nxt = inx;
            else
                column_[inx].deallocate();
            colidx_[inx] = nxt;
        }
    }

    // check if matrix is empty:
    if ( colidx_[0] == rsize_ )
        return false;
    
    sortElements();
    //printSummary(std::cout);
    return true;
}


//------------------------------------------------------------------------------
#pragma mark - Column Vector Multiplication


#if ( S_BLOCK_SIZE == 1 )
void SparMatSymBlk::Column::vecMulAdd1D(const real* X, real* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    const real X0 = X[jj];
    real D = blk_[0].value();
    real Y0 = Y[jj] + D * X0;
    assert_true(inx_[0]==jj);
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = inx_[n];
        const real M = blk_[n].value();
        Y[ii] += M * X0;
        Y0 += M * X[ii];
    }
    Y[jj] = Y0;
}
#endif


#if ( S_BLOCK_SIZE == 2 )
void SparMatSymBlk::Column::vecMulAdd2D(const real* X, real* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    const Vector2 xx(X+jj);
    assert_true(2*inx_[0]==jj);
    assert_small(blk_[0].asymmetry());
    Vector2 yy = blk_[0].vecmul(xx);
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 2 * inx_[n];
        Block const& M = blk_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.add_to(Y+jj);
}
#endif

#if ( S_BLOCK_SIZE == 3 )
void SparMatSymBlk::Column::vecMulAdd3D(const real* X, real* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    const Vector3 xxx(X+jj);
    assert_true(3*inx_[0]==jj);
    assert_small(blk_[0].asymmetry());
    Vector3 yyy = blk_[0].vecmul(xxx);
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 3 * inx_[n];
        Block const& M = blk_[n];
        M.vecmul(xxx).add_to(Y+ii);
        yyy += M.trans_vecmul(X+ii);
    }
    yyy.add_to(Y+jj);
}
#endif


#if ( S_BLOCK_SIZE == 4 )
void SparMatSymBlk::Column::vecMulAdd4D(const real* X, real* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    const Vector4 xxxx(X+jj);
    assert_true(4*inx_[0]==jj);
    assert_small(blk_[0].asymmetry());
    Vector4 yyyy = blk_[0].vecmul(xxxx);
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 4 * inx_[n];
        Block const& M = blk_[n];
        M.vecmul(xxxx).add_to(Y+ii);
        yyyy += M.trans_vecmul4_(X+ii);
    }
    yyyy.add_to(Y+jj);
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Single precision Optimized Vector Multiplication

#if ( S_BLOCK_SIZE == 3 ) && SMSB_USES_SSE && !REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd3D_SSE(const float* X, float* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(3*inx_[0] == jj);
    // load 3x3 matrix diagonal element into 3 vectors:
    Block const& D = blk_[0];
    
    //multiply with the diagonal block, assuming it is symmetric:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    const vec4f xxx = loadu4f(X+jj);
    vec4f s0 = mul4f(D.data0(), xxx);
    vec4f s1 = mul4f(D.data1(), xxx);
    vec4f s2 = mul4f(D.data2(), xxx);
    const vec4f x0 = clear4th(broadcastXf(xxx));
    const vec4f x1 = clear4th(broadcastYf(xxx));
    const vec4f x2 = clear4th(broadcastZf(xxx));

    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 3 * inx_[n];
        const vec4f m012 = blk_[n].data0();
        const vec4f m345 = blk_[n].data1();
        const vec4f m678 = blk_[n].data2();
        // multiply with the full block:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec4f z = fmadd4f(m012, x0, loadu4f(Y+ii));
        z = fmadd4f(m345, x1, z);
        z = fmadd4f(m678, x2, z);
        storeu4f(Y+ii, z);
        
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
        s0 = fmadd4f(m012, xyz, s0);
        s1 = fmadd4f(m345, xyz, s1);
        s2 = fmadd4f(m678, xyz, s2);
    }
    /* finally sum horizontally:
     s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
     to { Y0+Y0+Y0, Y1+Y1+Y1, Y2+Y2+Y2, 0 }
     */
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4f s3 = setzero4f();
    s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
    s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
    s0 = add4f(catshift2f(s0, s2), blend22f(s0, s2));
    assert_true(s0[3] == 0);
    storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
}
#endif


#if ( S_BLOCK_SIZE == 3 ) && SMSB_USES_SSE && !REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd3D_SSEU(const float* X, float* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(3*inx_[0] == jj);
    //std::cout << blk_[0].to_string(7,1); printf(" MSSB %lu : %lu\n", jj, nbb_);
    // load 3x3 matrix diagonal element into 3 vectors:
    Block const& D = blk_[0];
    
    //multiply with the diagonal block, assuming it has been symmetrized:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    const vec4f xxx = loadu4f(X+jj);
    vec4f s0 = mul4f(D.data0(), xxx);
    vec4f s1 = mul4f(D.data1(), xxx);
    vec4f s2 = mul4f(D.data2(), xxx);
    
    if ( nbb_ > 1 )
    {
        const vec4f x0 = clear4th(broadcastXf(xxx));
        const vec4f x1 = clear4th(broadcastYf(xxx));
        const vec4f x2 = clear4th(broadcastZf(xxx));

        index_t n = 1;
        {
            const index_t end = 1 + ( nbb_ & ~1 );
            // process 2 by 2
            #pragma nounroll
            for ( ; n < end; n += 2 )
            {
                const index_t ii = 3 * inx_[n];
                const vec4f m012 = blk_[n].data0();
                const vec4f m345 = blk_[n].data1();
                const vec4f m678 = blk_[n].data2();
                const index_t kk = 3 * inx_[n+1];
                const vec4f p012 = blk_[n+1].data0();
                const vec4f p345 = blk_[n+1].data1();
                const vec4f p678 = blk_[n+1].data2();
                assert_true( ii < kk );
                vec4f z = loadu4f(Y+ii);
                vec4f t = loadu4f(Y+kk);
                // multiply with the full block:
                z = fmadd4f(m012, x0, z);
                t = fmadd4f(p012, x0, t);
                vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
                vec4f tuv = loadu4f(X+kk);  // xyz = { X0 X1 X2 - }
                z = fmadd4f(m345, x1, z);
                t = fmadd4f(p345, x1, t);
                s0 = fmadd4f(m012, xyz, s0);
                s1 = fmadd4f(m345, xyz, s1);
                s2 = fmadd4f(m678, xyz, s2);
                z = fmadd4f(m678, x2, z);
                t = fmadd4f(p678, x2, t);
                s0 = fmadd4f(p012, tuv, s0);
                s1 = fmadd4f(p345, tuv, s1);
                s2 = fmadd4f(p678, tuv, s2);
                storeu4f(Y+ii, z);
                storeu4f(Y+kk, t);
            }
        }
        
        // process remaining blocks
#pragma nounroll
        for ( ; n < nbb_; ++n )
        {
            const index_t ii = 3 * inx_[n];
            const vec4f m012 = blk_[n].data0();
            const vec4f m345 = blk_[n].data1();
            const vec4f m678 = blk_[n].data2();
            // multiply with the full block:
            //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
            //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
            //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
            vec4f z = fmadd4f(m012, x0, loadu4f(Y+ii));
            z = fmadd4f(m345, x1, z);
            z = fmadd4f(m678, x2, z);
            storeu4f(Y+ii, z);
            
            // multiply with the transposed block:
            //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
            //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
            //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
            vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
            s0 = fmadd4f(m012, xyz, s0);
            s1 = fmadd4f(m345, xyz, s1);
            s2 = fmadd4f(m678, xyz, s2);
        }
    }
    /* finally sum horizontally:
     s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
     to { Y0+Y0+Y0, Y1+Y1+Y1, Y2+Y2+Y2, 0 }
     */
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4f s3 = setzero4f();
    s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
    s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
    s0 = add4f(catshift2f(s0, s2), blend22f(s0, s2));
    assert_true(s0[3] == 0);
    storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Double precision Optimized Vector Multiplication

#if ( S_BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSB_USES_SSE
void SparMatSymBlk::Column::vecMulAdd2D_SSE(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    vec2 s1 = load2(X+jj);
    vec2 x0 = unpacklo2(s1, s1);
    vec2 x1 = unpackhi2(s1, s1);
    //const real X0 = X[jj  ];
    //const real X1 = X[jj+1];

    // load 2x2 matrix element into 2 vectors:
    double const* D = blk_->data();
    //assume the block is already symmetrized:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1;
    // Y1 = Y[jj+1] + M[1] * X0 + M[3] * X1;
    vec2 s0 = mul2(load2(D), x0);
    s1 = mul2(load2(D+2), x1);
    
    Block const* blk = blk_ + 1;
    Block const* end = blk_ + nbb_;
    auto const* inx = inx_ + 1;

#if ( 1 )
    Block const* stop = blk + ( (end-blk) & ~1 );
    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( ; blk < stop; inx += 2, blk += 2 )
    {
        const auto i0 = 2 * inx[0];
        const auto i1 = 2 * inx[1];
        // load 2x2 matrix element into 2 vectors:
        double const* M = blk[0].data();
        double const* P = blk[1].data();
        vec2 m01 = load2(M);
        vec2 m23 = load2(M+2);
        vec2 p01 = load2(P);
        vec2 p23 = load2(P+2);
        vec2 xx = load2(X+i0);
        vec2 tt = load2(X+i1);

        // multiply with the full block:
        //Y[ii  ] += M[0] * X0 + M[2] * X1;
        //Y[ii+1] += M[1] * X0 + M[3] * X1;
        store2(Y+i0, fmadd2(m23, x1, fmadd2(m01, x0, load2(Y+i0))));
        store2(Y+i1, fmadd2(p23, x1, fmadd2(p01, x0, load2(Y+i1))));

        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
        //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
        s0 = fmadd2(m01, xx, s0);
        s1 = fmadd2(m23, xx, s1);
        s0 = fmadd2(p01, tt, s0);
        s1 = fmadd2(p23, tt, s1);
    }
#endif
    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( ; blk < end; ++inx, ++blk )
    {
        const auto ii = 2 * inx[0];
        // load 2x2 matrix element into 2 vectors:
        double const* M = blk[0].data();
        vec2 m0 = load2(M);
        vec2 m2 = load2(M+2);
        vec2 xx = load2(X+ii);

        // multiply with the full block:
        //Y[ii  ] += M[0] * X0 + M[2] * X1;
        //Y[ii+1] += M[1] * X0 + M[3] * X1;
        store2(Y+ii, fmadd2(m2, x1, fmadd2(m0, x0, load2(Y+ii))));

        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
        //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
        s0 = fmadd2(m0, xx, s0);
        s1 = fmadd2(m2, xx, s1);
    }
    //Y[jj  ] = Y0;
    //Y[jj+1] = Y1;
    s0 = add2(unpacklo2(s0, s1), unpackhi2(s0, s1));
    store2(Y+jj, add2(s0, load2(Y+jj)));
}
#endif

#if ( S_BLOCK_SIZE == 2 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd2D_AVX(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(2*inx_[0] == jj);
    // xy = { X0 X1 X0 X1 }
    vec4 xy = broadcast2(X+jj);
    //multiply with full block, assuming it is symmetric:
    // Y0 = M[0] * X0 + M[1] * X1;
    // Y1 = M[1] * X0 + M[3] * X1;
    
    // yyyy = { Y0 Y0 Y1 Y1 }
    // load 2x2 matrix element into 2 vectors:
    vec4 ss = mul4(blk_[0].data0(), xy);

    //const real X0 = X[jj  ];
    //const real X1 = X[jj+1];
    // xxyy = { X0 X0 X1 X1 }
    const vec4 xxyy = duplohi4(xy);

    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 2 * inx_[n];
        vec4 mat = blk_[n].data0();     // load 2x2 matrix
        vec4 yy = load2Z(Y+ii);         // yy = { Y0 Y1 0 0 }
        vec4 xx = broadcast2(X+ii);     // xx = { X0 X1 X0 X1 }

        // multiply with the full block:
        //Y[ii  ] += M[0] * X0 + M[2] * X1;
        //Y[ii+1] += M[1] * X0 + M[3] * X1;
        vec4 u = fmadd4(mat, xxyy, yy);
        store2(Y+ii, add2(getlo(u), gethi(u)));
        
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
        //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
        ss = fmadd4(mat, xx, ss);
    }
    // need to collapse yyyy = { S0 S0 S1 S1 }
    // Y[jj  ] += yyyy[0] + yyyy[1];
    // Y[jj+1] += yyyy[2] + yyyy[3];
    vec2 yy = load2(Y+jj);
    vec2 h = gethi(ss);
    store2(Y+jj, add2(yy, add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h))));
}
#endif


#if ( S_BLOCK_SIZE == 2 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
static inline void multiply2D(double const* X, double* Y, index_t ii, vec4 const& mat, vec4 const& xxxx, vec4& ss)
{
    vec4 xx = broadcast2(X+ii);
    vec4 u = fmadd4(mat, xxxx, load2Z(Y+ii));
    store2(Y+ii, add2(getlo(u), gethi(u)));
    ss = fmadd4(mat, xx, ss);
}
#endif


#if ( S_BLOCK_SIZE == 2 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd2D_AVXU(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(2*inx_[0] == jj);
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(blk_[0].data0(), xyxy);
    const vec4 xxyy = duplohi4(xyxy);
    vec4 s1 = setzero4();

    Block const* blk = blk_ + 1;
    Block const* end = blk_ + 1 + ( (nbb_-1) & ~1 );
    auto const* inx = inx_ + 1;
    // process 2 by 2:
    #pragma nounroll
    for ( ; blk < end; blk += 2, inx += 2 )
    {
#if ( 0 )
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        multiply2D(X, Y, 2*inx[0], blk[0].data0(), xxyy, ss);
        multiply2D(X, Y, 2*inx[1], blk[1].data0(), xxyy, s1);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const auto i0 = 2 * inx[0];
        const auto i1 = 2 * inx[1];
        assert_true( i0 < i1 );
        vec4 mat0 = blk[0].data0();
        vec4 mat1 = blk[1].data0();
        vec4 u0 = fmadd4(mat0, xxyy, load2Z(Y+i0));
        vec4 u1 = fmadd4(mat1, xxyy, load2Z(Y+i1));
        ss = fmadd4(mat0, broadcast2(X+i0), ss);
        s1 = fmadd4(mat1, broadcast2(X+i1), s1);
        store2(Y+i0, add2(getlo(u0), gethi(u0)));
        store2(Y+i1, add2(getlo(u1), gethi(u1)));
#endif
    }
    // collapse 'ss'
    ss = add4(ss, s1);
    // process remaining blocks:
    end = blk_ + nbb_;
    //#pragma nounroll
    if ( blk < end ) //for ( ; blk < end; ++blk, ++inx )
        multiply2D(X, Y, 2*inx[0], blk[0].data0(), xxyy, ss);
    /* finally horizontally sum ss = { SX SX SY SY } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


#if ( S_BLOCK_SIZE == 2 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd2D_AVXUU(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(2*inx_[0] == jj);
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(blk_[0].data0(), xyxy);
    const vec4 xxyy = duplohi4(xyxy);
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();

    Block const* blk = blk_ + 1;
    Block const* end = blk_ + 1 + ( (nbb_-1) & ~3 );
    auto const* inx = inx_ + 1;
    // process 4 by 4:
    #pragma nounroll
    for ( ; blk < end; blk += 4, inx += 4 )
    {
#if ( 0 )
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        multiply2D(X, Y, 2*inx[0], blk[0].data0(), xxyy, ss);
        multiply2D(X, Y, 2*inx[1], blk[1].data0(), xxyy, s1);
        multiply2D(X, Y, 2*inx[2], blk[2].data0(), xxyy, s2);
        multiply2D(X, Y, 2*inx[3], blk[3].data0(), xxyy, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        assert_true( inx[0] < inx[1] );
        assert_true( inx[1] < inx[2] );
        assert_true( inx[2] < inx[3] );
        const auto i0 = 2 * inx[0];
        const auto i1 = 2 * inx[1];
        const auto i2 = 2 * inx[2];
        const auto i3 = 2 * inx[3];
        vec4 mat0 = blk[0].data0();
        vec4 mat1 = blk[1].data0();
        vec4 mat2 = blk[2].data0();
        vec4 mat3 = blk[3].data0();
        vec4 u0 = fmadd4(mat0, xxyy, load2Z(Y+i0));
        vec4 u1 = fmadd4(mat1, xxyy, load2Z(Y+i1));
        vec4 u2 = fmadd4(mat2, xxyy, load2Z(Y+i2));
        vec4 u3 = fmadd4(mat3, xxyy, load2Z(Y+i3));
        ss = fmadd4(mat0, broadcast2(X+i0), ss);
        s1 = fmadd4(mat1, broadcast2(X+i1), s1);
        s2 = fmadd4(mat2, broadcast2(X+i2), s2);
        s3 = fmadd4(mat3, broadcast2(X+i3), s3);
        store2(Y+i0, add2(getlo(u0), gethi(u0)));
        store2(Y+i1, add2(getlo(u1), gethi(u1)));
        store2(Y+i2, add2(getlo(u2), gethi(u2)));
        store2(Y+i3, add2(getlo(u3), gethi(u3)));
#endif
    }
    // collapse 'ss'
    ss = add4(add4(ss,s1), add4(s2,s3));
    // process remaining blocks:
    end = blk_ + nbb_;
    //#pragma nounroll
    if ( blk < end ) //for ( ; blk < end; ++blk, ++inx )
        multiply2D(X, Y, 2*inx[0], blk[0].data0(), xxyy, ss);
    /* finally sum ss = { S0 S0 S1 S1 } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


#if ( S_BLOCK_SIZE == 3 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd3D_AVX(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(3*inx_[0] == jj);
    // load 3x3 matrix diagonal element into 3 vectors:
    Block const& D = blk_[0];
    
    //multiply with the diagonal block, assuming it is symmetric:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    vec4 s0, s1, s2;
    vec4 x0, x1, x2;
    {
        vec4 xxx = load3Z(X+jj);
        s0 = mul4(D.data0(), xxx);
        s1 = mul4(D.data1(), xxx);
        s2 = mul4(D.data2(), xxx);
        // sum non-diagonal elements:
#if ( 0 )
        x0 = broadcast1(X+jj);
        x1 = broadcast1(X+jj+1);
        x2 = broadcast1(X+jj+2);
#else
        x0 = swap2f128(xxx);
        x1 = blend22(xxx, x0);
        x2 = blend22(x0, xxx);
        x0 = duplo4(x1);
        x1 = duphi4(x1);
        x2 = duplo4(x2);
#endif
        // zero out the 4th terms:
        x0 = clear4th(x0);
        x1 = clear4th(x1);
        x2 = clear4th(x2);
    }
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 3 * inx_[n];
        const vec4 m012 = blk_[n].data0();
        const vec4 m345 = blk_[n].data1();
        const vec4 m678 = blk_[n].data2();
        // multiply with the full block:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec4 z = fmadd4(m012, x0, loadu4(Y+ii));
        z = fmadd4(m345, x1, z);
        z = fmadd4(m678, x2, z);
        storeu4(Y+ii, z);
        
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        vec4 xyz = load3Z(X+ii);  // xyz = { X0 X1 X2 0 }
        s0 = fmadd4(m012, xyz, s0);
        s1 = fmadd4(m345, xyz, s1);
        s2 = fmadd4(m678, xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
#if ( 0 )
    Y[jj  ] += s0[0] + s0[1] + s0[2];
    Y[jj+1] += s1[0] + s1[1] + s1[2];
    Y[jj+2] += s2[0] + s2[1] + s2[2];
#else
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4 s3 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    s1 = add4(catshift2(s0, s2), blend22(s0, s2));
    assert_true(s1[3] == 0);
    storeu4(Y+jj, add4(loadu4(Y+jj), s1));
#endif
}
#endif


#if ( S_BLOCK_SIZE == 3 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd3D_AVXU(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(3*inx_[0] == jj);
    Block const* blk = blk_;

    vec4 s0, s1, s2;
    vec4 x0, x1, x2;
    vec4 t0 = setzero4();
    vec4 t1 = setzero4();
    vec4 t2 = setzero4();
    // load 3x3 matrix element into 3 vectors:
    {
        vec4 xxx = load3Z(X+jj);
        // multiply by diagonal elements:
        s0 = mul4(blk->data0(), xxx);
        s1 = mul4(blk->data1(), xxx);
        s2 = mul4(blk->data2(), xxx);
        // prepare broadcasted vectors:
        x0 = swap2f128(xxx);
        x1 = blend22(xxx, x0);
        x2 = blend22(x0, xxx);
        // zero out the 4th terms:
        x0 = clear4th(duplo4(x1));
        x1 = clear4th(duphi4(x1));
        x2 = clear4th(duplo4(x2));
    }
    ++blk;
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    Block const* end = blk_ + 1 + ( (nbb_-1) & ~1 );
    auto const* inx = inx_ + 1;
    /*
     Unrolling will reduce the dependency chain, which may be limiting the
     throughput here. However the number of registers (16 for AVX CPU) limits
     the level of unrolling that can be done.
     */
    //process 2 by 2:
    #pragma nounroll
    for ( ; blk < end; blk += 2, inx += 2 )
    {
        const index_t ii = 3 * inx[0];
        const index_t kk = 3 * inx[1];
        assert_true( ii < kk );
        //printf("--- %4i %4i\n", ii, kk);
        vec4 M0 = blk[0].data0();
        vec4 P0 = blk[1].data0();
        vec4 z = fmadd4(M0, x0, loadu4(Y+ii));
        vec4 t = fmadd4(P0, x0, loadu4(Y+kk));
        vec4 xyz = loadu4(X+ii);
        vec4 tuv = loadu4(X+kk);
        s0 = fmadd4(M0, xyz, s0);
        t0 = fmadd4(P0, tuv, t0);
        // multiply with the full block:
        vec4 M1 = blk[0].data1();
        vec4 P1 = blk[1].data1();
        z = fmadd4(M1, x1, z);
        t = fmadd4(P1, x1, t);
        s1 = fmadd4(M1, xyz, s1);
        t1 = fmadd4(P1, tuv, t1);
        vec4 M2 = blk[0].data2();
        vec4 P2 = blk[1].data2();
        z = fmadd4(M2, x2, z);
        t = fmadd4(P2, x2, t);
        s2 = fmadd4(M2, xyz, s2);
        t2 = fmadd4(P2, tuv, t2);
        /*
         Attention: the 4th elements of the vectors z0 and z1 would be correct,
         because only zero was added to the value loaded from 'Y'. However, in the
         case where the indices ii and kk are consecutive and reverted (kk < ii),
         the value stored in z0 would not have been updated giving a wrong results.
         The solution is to either use a 'store3(Y+kk, z1)', or to make sure that
         indices are non-consecutive or ordered in the column in increasing order.
         This affects performance since 'store3' is slower than 'storeu4'
         */
        assert_true(z[3]==Y[ii+3]);
        assert_true(t[3]==Y[kk+3]);
        storeu4(Y+ii, z);
        storeu4(Y+kk, t);
    }
    s0 = add4(s0, t0);
    s1 = add4(s1, t1);
    s2 = add4(s2, t2);
    
    // process remaining blocks:
    #pragma nounroll
    for ( end = blk_ + nbb_; blk < end; ++blk )
    {
        const index_t ii = 3 * inx[0];
        ++inx;
        //printf("--- %4i\n", ii);
        vec4 ma = blk->data0();
        vec4 z = fmadd4(ma, x0, loadu4(Y+ii));
        vec4 xyz = loadu4(X+ii);
        s0 = fmadd4(ma, xyz, s0);
        
        vec4 mb = blk->data1();
        z = fmadd4(mb, x1, z);
        s1 = fmadd4(mb, xyz, s1);
        
        vec4 mc = blk->data2();
        z = fmadd4(mc, x2, z);
        s2 = fmadd4(mc, xyz, s2);
        storeu4(Y+ii, z);
    }
    // finally sum s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    x0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s1 = add4(unpacklo4(s2, x0), unpackhi4(s2, x0));
    s0 = add4(catshift2(s0, s1), blend22(s0, s1));
    assert_true(s0[3] == 0);
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
}
#endif


#if ( S_BLOCK_SIZE == 4 ) && SMSB_USES_AVX && REAL_IS_DOUBLE
void SparMatSymBlk::Column::vecMulAdd4D_AVX(const double* X, double* Y, index_t jj) const
{
    assert_true(nbb_ > 0);
    assert_true(3*inx_[0] == jj);
    Block const* blk = blk_;
    //multiply with the diagonal block, assuming it is symmetric:
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    vec4 s0, s1, s2, s3;
    vec4 x0, x1, x2, x3;
    {
        vec4 tt = load4(X+jj);
        s0 = mul4(blk->data0(), tt);
        s1 = mul4(blk->data1(), tt);
        s2 = mul4(blk->data2(), tt);
        s3 = mul4(blk->data3(), tt);
        // sum non-diagonal elements:
#if ( 0 )
        x0 = broadcast1(X+jj);
        x1 = broadcast1(X+jj+1);
        x2 = broadcast1(X+jj+2);
        x3 = broadcast1(X+jj+3);
#else
        x1 = duplo2f128(tt);
        x3 = duphi2f128(tt);
        x0 = duplo4(x1);
        x1 = duphi4(x1);
        x2 = duplo4(x3);
        x3 = duphi4(x3);
#endif
    }
    ++blk;
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( index_t n = 1; n < nbb_; ++n )
    {
        const index_t ii = 4 * inx_[n];
        const vec4 yy = load4(Y+ii);
        const vec4 xx = load4(X+ii);  // xx = { X0 X1 X2 X3 }
        const vec4 m0 = blk->data0();
        vec4 z = fmadd4(m0, x0, yy);
        s0 = fmadd4(m0, xx, s0);
        
        const vec4 m1 = blk->data1();
        z  = fmadd4(m1, x1, z);
        s1 = fmadd4(m1, xx, s1);

        const vec4 m2 = blk->data2();
        z  = fmadd4(m2, x2, z);
        s2 = fmadd4(m2, xx, s2);

        const vec4 m3 = blk->data3();
        z  = fmadd4(m3, x3, z);
        s3 = fmadd4(m3, xx, s3);
        store4(Y+ii, z);
    }
    // finally sum s0 = { Y0 Y0 Y0 Y0 }, s1 = { Y1 Y1 Y1 Y1 }, s2 = { Y2 Y2 Y2 Y2 }
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    s1 = add4(catshift2(s0, s2), blend22(s0, s2));
    store4(Y+jj, add4(load4(Y+jj), s1));
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector Add-multiply


// multiplication of a vector: Y = Y + M * X
void SparMatSymBlk::vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    for ( index_t j = start; j < stop; ++j )
    if ( column_[j].notEmpty() )
    {
        //std::clog << "SparMatSymBlk column " << j << "  " << rsize_ << " \n";
#if ( S_BLOCK_SIZE == 1 )
        column_[j].vecMulAdd1D(X, Y, j);
#elif ( S_BLOCK_SIZE == 2 )
        column_[j].vecMulAdd2D(X, Y, 2*j);
#elif ( S_BLOCK_SIZE == 3 )
        column_[j].vecMulAdd3D(X, Y, 3*j);
#elif ( S_BLOCK_SIZE == 4 )
        column_[j].vecMulAdd4D(X, Y, 4*j);
#endif
    }
}


#if SMSB_USES_AVX && REAL_IS_DOUBLE
#   define VECMULADD2D vecMulAdd2D_AVXU
#   define VECMULADD3D vecMulAdd3D_AVXU
#   define VECMULADD4D vecMulAdd4D_AVX
#elif SMSB_USES_SSE && REAL_IS_DOUBLE
#   define VECMULADD2D vecMulAdd2D_SSE
#   define VECMULADD3D vecMulAdd3D
#   define VECMULADD4D vecMulAdd4D
#elif SMSB_USES_SSE
#   define VECMULADD2D vecMulAdd2D
#   define VECMULADD3D vecMulAdd3D_SSE
#   define VECMULADD4D vecMulAdd4D
#else
#   define VECMULADD2D vecMulAdd2D
#   define VECMULADD3D vecMulAdd3D
#   define VECMULADD4D vecMulAdd4D
#endif


// multiplication of a vector: Y = Y + M * X
void SparMatSymBlk::vecMulAdd(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    for ( index_t j = colidx_[start]; j < stop; j = colidx_[j+1] )
    {
        //std::clog << "SparMatSymBlk column " << j << "  " << rsize_ << " \n";
#if ( S_BLOCK_SIZE == 1 )
        column_[j].vecMulAdd1D(X, Y, j);
#elif ( S_BLOCK_SIZE == 2 )
        column_[j].VECMULADD2D(X, Y, 2*j);
#elif ( S_BLOCK_SIZE == 3 )
        column_[j].VECMULADD3D(X, Y, 3*j);
#elif ( S_BLOCK_SIZE == 4 )
        column_[j].VECMULADD4D(X, Y, 4*j);
#endif
    }
}

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication

void SparMatSymBlk::vecMul(const real* X, real* Y) const
{
    zero_real(S_BLOCK_SIZE*rsize_, Y);
    vecMulAdd(X, Y, 0, rsize_);
}
