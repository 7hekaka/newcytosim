// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include <cmath>
#include "sparmatsymblkdiag.h"
#include "assert_macro.h"
#include "vector2.h"
#include "vector3.h"
#include <sstream>

// Flags to enable SIMD implementation
#if defined(__AVX__)
#  include "simd.h"
#  include "simd_float.h"
#  define SMSBD_USES_AVX 1
#  define SMSBD_USES_SSE 1
#elif defined(__SSE3__)
#  include "simd.h"
#  include "simd_float.h"
#  define SMSBD_USES_AVX 0
#  define SMSBD_USES_SSE 1
#else
#  define SMSBD_USES_AVX 0
#  define SMSBD_USES_SSE 0
#endif


SparMatSymBlkDiag::SparMatSymBlkDiag()
{
    size_  = 0;
    alloc_ = 0;
    pilar_ = nullptr;
    colix_ = new size_t[2];
    colix_[0] = 0;
}


void SparMatSymBlkDiag::allocate(size_t alc)
{
    if ( alc > alloc_*BLOCK_SIZE )
    {
        /*
         'chunk' can be adjusted to tune performance: by increasing 'chunk',
          more memory will be used, but reallocation will be less frequent
        */
        constexpr size_t chunk = 8;
        alc = ( alc/BLOCK_SIZE + chunk - 1 ) & ~( chunk -1 );
        
        //fprintf(stderr, "SMSBD allocates %u\n", alc);
        Column * ptr = new Column[alc];
        
        if ( pilar_ )
        {
            for (size_t n = 0; n < alloc_; ++n )
                ptr[n] = pilar_[n];
            delete[] pilar_;
        }
        
        pilar_ = ptr;
        alloc_ = alc;
        
        delete[] colix_;
        colix_ = new size_t[alc+1];
        for ( size_t n = 0; n <= alc; ++n )
            colix_[n] = n;
    }
}


void SparMatSymBlkDiag::deallocate()
{
    delete[] pilar_;
    delete[] colix_;
    pilar_ = nullptr;
    colix_ = nullptr;
    alloc_ = 0;
}


//------------------------------------------------------------------------------
#pragma mark - Column

SparMatSymBlkDiag::Column::Column()
{
    size_ = 0;
    allo_ = 0;
    dia_.reset();
    inx_ = nullptr;
    blk_ = nullptr;
}


/*
\todo Columns should use partitions of a single memory pool allocated by Matrix
This may require some smart allocation scheme.
*/
void SparMatSymBlkDiag::Column::allocate(size_t alc)
{
    if ( alc > allo_ )
    {
        //if ( inx_ ) fprintf(stderr, "SMSBD reallocates column %lu for %lu\n", inx_[0], alc);
        //else fprintf(stderr, "SMSBD allocates column for %lu\n", alc);
        /*
         'chunk' can be adjusted to tune performance: by increasing 'chunk',
          more memory will be used, but reallocation will be less frequent
        */
        constexpr size_t chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk - 1 );
        
        // use aligned memory:
        void * ptr = new_real(alc*SB);
        Block * blk_new  = new(ptr) Block[alc];

        if ( posix_memalign(&ptr, 32, alc*sizeof(size_t)) )
            throw std::bad_alloc();
        size_t * inx_new = (size_t*)ptr;

        if ( inx_ )
        {
            for ( size_t n = 0; n < size_; ++n )
                inx_new[n] = inx_[n];
            free(inx_);
        }

        if ( blk_ )
        {
            for ( size_t n = 0; n < size_; ++n )
                blk_new[n] = blk_[n];
            free_real(blk_);
        }
        inx_  = inx_new;
        blk_  = blk_new;
        allo_ = alc;
        
        //std::clog << "Column " << this << "  " << alc << ": ";
        //std::clog << " alignment " << ((uintptr_t)elem_ & 63) << "\n";
    }
}


void SparMatSymBlkDiag::Column::deallocate()
{
    //if ( inx_ ) fprintf(stderr, "SMSBD deallocates column %lu of size %lu\n", inx_[0], allo_);
    free(inx_);
    free_real(blk_);
    inx_ = nullptr;
    blk_ = nullptr;
    allo_ = 0;
    size_ = 0;
}


void SparMatSymBlkDiag::Column::operator =(SparMatSymBlkDiag::Column & col)
{
    //if ( inx_ ) fprintf(stderr, "SMSBD transfers column %u\n", inx_[0]);
    free(inx_);
    free_real(blk_);

    size_ = col.size_;
    allo_ = col.allo_;
    dia_ = col.dia_;
    inx_ = col.inx_;
    blk_ = col.blk_;
    
    col.size_ = 0;
    col.allo_ = 0;
    col.dia_.reset();
    col.inx_ = nullptr;
    col.blk_ = nullptr;
}

/**
 This allocates to be able to hold the matrix element if necessary
 */
SparMatSymBlkDiag::Block& SparMatSymBlkDiag::Column::block(size_t ii, size_t jj)
{
#if 0
    // this is unnecessary
    if ( ii == jj )
        return dia_;
#endif
    
    assert_true( ii > jj );
    if ( size_ > 0 )
    {
        /* This is a silly search that could be optimized */
        for ( size_t n = 0; n < size_; ++n )
            if ( inx_[n] == ii )
                return blk_[n];
    }
    
    // add the requested term:
    size_t n = size_;
    
    // allocate space for new Element if necessary:
    if ( n >= allo_ )
        allocate(n+1);
    
    assert_true( n < allo_ );
    inx_[n] = ii;
    blk_[n].reset();
    size_ = n + 1;
    
    //printColumn(jj);
    return blk_[n];
}


void SparMatSymBlkDiag::Column::reset()
{
    size_ = 0;
    dia_.reset();
}


real& SparMatSymBlkDiag::operator()(size_t iii, size_t jjj)
{
    // branchless code to address lower triangle
    size_t ii = std::max(iii, jjj);
    size_t jj = std::min(iii, jjj);
#if ( BLOCK_SIZE == 1 )
    return pilar_[jj/BLOCK_SIZE].block(ii, jj).value();
#else
    size_t i = ii % BLOCK_SIZE;
    size_t j = jj / BLOCK_SIZE;
    size_t jb = j * BLOCK_SIZE;
    if ( ii-i == j*BLOCK_SIZE )
        return pilar_[j].dia_(i, jj-jb);
    return pilar_[j].block(ii-i, jb)(i, jj-jb);
#endif
}


real* SparMatSymBlkDiag::addr(size_t iii, size_t jjj) const
{
    // branchless code to address lower triangle
    size_t ii = std::max(iii, jjj);
    size_t jj = std::min(iii, jjj);
#if ( BLOCK_SIZE == 1 )
    return &pilar_[jj/BLOCK_SIZE].block(ii, jj).value();
#else
    size_t i = ii % BLOCK_SIZE;
    size_t j = jj % BLOCK_SIZE;
    size_t jb = j * BLOCK_SIZE;
    if ( ii-i == j*BLOCK_SIZE )
        return pilar_[j].dia_.addr(i, jj-jb);
    return pilar_[j].block(ii-i, jb).addr(i, jj-jb);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Stuff

void SparMatSymBlkDiag::reset()
{
    for ( size_t n = 0; n < alloc_; ++n )
        pilar_[n].reset();
}


bool SparMatSymBlkDiag::isNotZero() const
{
    //check for any non-zero sparse term:
    for ( size_t jj = 0; jj < size_/BLOCK_SIZE; ++jj )
    {
        Column & col = pilar_[jj];
        if ( col.dia_ != 0.0 )
            return true;
        for ( size_t n = 0 ; n < col.size_ ; ++n )
            if ( col[n] != 0.0 )
                return true;
    }
    //if here, the matrix is empty
    return false;
}


void SparMatSymBlkDiag::scale(const real alpha)
{
    for ( size_t jj = 0; jj < size_/BLOCK_SIZE; ++jj )
    {
        Column & col = pilar_[jj];
        col.dia_.scale(alpha);
        for ( size_t n = 0 ; n < col.size_ ; ++n )
            col[n].scale(alpha);
    }
}


void SparMatSymBlkDiag::addDiagonalBlock(real* mat, size_t ldd,
                                         const size_t start, const size_t cnt) const
{
    assert_false( start % BLOCK_SIZE );
    assert_false( cnt % BLOCK_SIZE );

    size_t end = start + cnt;
    size_t off = start + ldd * start;
    assert_true( end <= size_ );
    
    for ( size_t jj = start; jj < end; jj += BLOCK_SIZE )
    {
        Column & col = pilar_[jj/BLOCK_SIZE];
        col.dia_.addto_symm(mat+(1+ldd)*jj-off, ldd);
        for ( size_t n = 0; n < col.size_; ++n )
        {
            size_t ii = col.inx_[n];
            // assuming lower triangle is stored:
            assert_true(ii>jj);
            if ( ii < end )
            {
                //fprintf(stderr, "SMSBD %4lu %4lu\n", ii, jj); col[n].print(stderr);
                col[n].addto(mat+(ii+ldd*jj)-off, ldd);
                col[n].addto_trans(mat+(jj+ldd*ii)-off, ldd);
            }
        }
    }
}


void SparMatSymBlkDiag::addLowerBand(real alpha, real* mat, size_t ldd,
                                     const size_t start, const size_t cnt, size_t rank) const
{
    assert_false( start % BLOCK_SIZE );
    assert_false( cnt % BLOCK_SIZE );

    size_t end = start + cnt;
    size_t off = start + ldd * start;
    assert_true( end <= size_ );
    
    for ( size_t jj = start; jj < end; jj += BLOCK_SIZE )
    {
        Column & col = pilar_[jj/BLOCK_SIZE];
        col.dia_.addto_lower(mat+(1+ldd)*jj-off, ldd, alpha);
        for ( size_t n = 0; n < col.size_; ++n )
        {
            size_t ii = col.inx_[n];
            // assuming lower triangle is stored:
            assert_true(ii>jj);
            if ((ii <= jj+rank) & (ii < end))
            {
                //fprintf(stderr, "SMSBD %4lu %4lu\n", ii, jj); col[n].print(stderr);
                col[n].addto(mat+(ii+ldd*jj)-off, ldd, alpha);
                //col[n].addto_trans(mat+(jj+ldd*ii)-off, ldd, alpha);
            }
        }
    }
}

/*
addresses `mat' using lower banded storage for a symmetric matrix
mat(i, j) is stored in mat[i-j+ldd*j]
*/
void SparMatSymBlkDiag::addDiagonalTrace(real alpha, real* mat, size_t ldd,
                                         const size_t start, const size_t cnt,
                                         const size_t rank, bool sym) const
{
    assert_false( start % BLOCK_SIZE );
    assert_false( cnt % BLOCK_SIZE );

    size_t end = start + cnt;
    assert_true( end <= size_ );

    for ( size_t jj = start; jj < end; jj += BLOCK_SIZE )
    {
        Column & col = pilar_[jj/BLOCK_SIZE];
        size_t j = ( jj - start ) / BLOCK_SIZE;
        // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
        mat[j+ldd*j] += alpha * col.dia_.trace();  // diagonal term
        for ( size_t n = 0; n < col.size_; ++n )
        {
            size_t ii = col.inx_[n];
            // assuming lower triangle is stored:
            assert_false( ii % BLOCK_SIZE );
            if (( ii < end ) & ( ii <= jj+rank ))
            {
                size_t i = ( ii - start ) / BLOCK_SIZE;
                assert_true( i > j );
                real a = alpha * col[n].trace();
                //fprintf(stderr, "SMSBD %4lu %4lu : %.4f\n", i, j, a);
                // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
                mat[i+ldd*j] += a;
                if ( sym ) mat[j+ldd*i] += a;
            }
        }
    }
}



int SparMatSymBlkDiag::bad() const
{
    size_t stop = size_ / BLOCK_SIZE;
    for ( size_t j = 0; j < stop; ++j )
    {
        Column & col = pilar_[j];
        for ( size_t n = 0 ; n < col.size_ ; ++n )
        {
            if ( col.inx_[n] >= size_ ) return 2;
            if ( col.inx_[n] <= j ) return 3;
        }
    }
    return 0;
}


/** all elements are counted, even if zero */
size_t SparMatSymBlkDiag::nbElements(size_t start, size_t stop, size_t& alc) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_) / BLOCK_SIZE;
    alc = 0;
    size_t cnt = stop - start; // counting diagonal elements
    for ( size_t i = start/BLOCK_SIZE; i < stop; ++i )
    {
        cnt += pilar_[i].size_;
        alc += pilar_[i].allo_;
    }
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark - I/O


std::string SparMatSymBlkDiag::what() const
{
    size_t alc = 0;
    size_t cnt = nbElements(0, size_, alc);
    std::ostringstream msg;
#if SMSBD_USES_AVX && REAL_IS_DOUBLE
    msg << "SMSBDx ";
#elif SMSBD_USES_SSE && REAL_IS_DOUBLE
    msg << "SMSBDe ";
#else
    msg << "SMSBD ";
#endif
    msg << Block::what() << "*" << cnt << " (" << alc*SB << ")";
    return msg.str();
}


void SparMatSymBlkDiag::printSparse(std::ostream& os, real inf, size_t start, size_t stop) const
{
    stop = std::min(stop, size_) / BLOCK_SIZE;

    char str[256];
    std::streamsize p = os.precision();
    os.precision(8);
    if ( ! pilar_ )
        return;
    for ( size_t j = start/BLOCK_SIZE; j < stop; ++j )
    {
        size_t jj = j * BLOCK_SIZE;
        Column & col = pilar_[j];
        os << "% column " << jj << "\n";
        Block D = col.dia_;
        for ( size_t y = 0; y < BLOCK_SIZE; ++y )
        for ( size_t x = y; x < BLOCK_SIZE; ++x )
        {
            real v = D(y, x);
            if ( abs_real(v) >= inf )
            {
                snprintf(str, sizeof(str), "%6lu %6lu %16.6f\n", jj+y, jj+x, v);
                os << str;
            }
        }
        for ( size_t n = 0 ; n < col.size_ ; ++n )
        {
            size_t ii = col.inx_[n];
            Block B = col.blk_[n];
            for ( size_t y = 0; y < BLOCK_SIZE; ++y )
            for ( size_t x = 0; x < BLOCK_SIZE; ++x )
            {
                real v = B(y, x);
                if ( abs_real(v) >= inf )
                {
                    snprintf(str, sizeof(str), "%6lu %6lu %16.6f\n", ii+y, jj+x, v);
                    os << str;
                }
            }
        }
    }
    os.precision(p);
}


void SparMatSymBlkDiag::printColumns(std::ostream& os, size_t start, size_t stop)
{
    stop = std::min(stop, size_)/BLOCK_SIZE;
    os << "\nSMSBD size " << size_ << ":";
    for ( size_t j = start/BLOCK_SIZE; j < stop; ++j )
    {
        if ( pilar_[j].size_ > 0 )
        {
            os << "\n   " << j << "   " << pilar_[j].size_;
            os << " ---> " << colix_[j];
        }
    }
    std::endl(os);
}


void SparMatSymBlkDiag::Column::print(std::ostream& os) const
{
    for ( size_t n = 0; n < size_; ++n )
        os << "\n" << inx_[n] << " : " << blk_[n] << "\n";
    std::endl(os);
}


//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication


/// A block element of the sparse matrix suitable for qsort()
class alignas(4*sizeof(real)) SparMatSymBlkDiag::Element
{
public:
    /// block element
    real blk[BLOCK_SIZE*BLOCK_SIZE];

    /// index
    size_t inx;
};


/// function for qsort, comparing line indices
static int compareSMSBDElement(const void * A, const void * B)
{
    size_t a = static_cast<SparMatSymBlkDiag::Element const*>(A)->inx;
    size_t b = static_cast<SparMatSymBlkDiag::Element const*>(B)->inx;
    
    return ( a > b ) - ( b > a );
}

/**
 This copies the data to the provided temporary array
 */
void SparMatSymBlkDiag::Column::sortElements(Element tmp[], size_t tmp_size)
{
    assert_true( size_ <= tmp_size );
    for ( size_t i = 0; i < size_; ++i )
    {
        blk_[i].store(tmp[i].blk);
        tmp[i].inx = inx_[i];
    }
    
    //std::clog << "sizeof(SparMatSymBlkDiag::Element) " << sizeof(Element) << "\n";
    qsort(tmp, size_, sizeof(Element), &compareSMSBDElement);
    
    for ( size_t i = 0; i < size_; ++i )
    {
         blk_[i].load(tmp[i].blk);
         inx_[i] = tmp[i].inx;
    }
}


size_t SparMatSymBlkDiag::newElements(Element*& ptr, size_t cnt)
{
    constexpr size_t chunk = 16;
    size_t all = ( cnt + chunk - 1 ) & ~( chunk - 1 );
    free(ptr);  // Element has no destructor
    void* tmp = nullptr;
    if ( posix_memalign(&tmp, 32, all*sizeof(Element)) )
        throw std::bad_alloc();
    ptr = new(tmp) Element[all];
    return all;
}


void SparMatSymBlkDiag::sortElements()
{
    size_t tmp_size = 0;
    Element * tmp = nullptr;
    
    size_t last = size_/BLOCK_SIZE;
    for ( size_t j = 0; j < last; ++j )
    {
        Column & col = pilar_[j];
        //std::clog << "SMSBD column " << j << " has 1+" << col.size_ << " elements\n";

        if ( col.size_ > 1 )
        {
            // order the elements within the column:
            if ( tmp_size < col.size_ )
                tmp_size = newElements(tmp, col.size_);
            col.sortElements(tmp, tmp_size);
        }

#ifndef NDEBUG
        for ( size_t n = 0 ; n < col.size_ ; ++n )
        {
            const size_t i = col.inx_[n];
            assert_true( i < size_ );
            assert_true( i != j );
        }
#endif
    }
    // release memory:
    free(tmp);
}


bool SparMatSymBlkDiag::prepareForMultiply(int)
{
    size_t last = size_/BLOCK_SIZE;
    if ( size_ > 0 )
    {
        size_t inx = last;
        size_t nxt = last;
        while ( inx-- > 0 )
        {
            if ( pilar_[inx].size_ > 0 )
                nxt = inx;
            colix_[inx] = nxt;
        }
    }
    colix_[last] = last;

    for ( size_t j = 0; j < last; ++j )
    {
        pilar_[j].dia_.copy_lower();
    }
    
    sortElements();

    //printColumns(std::cout, 0, size_);
    return true;
}


//------------------------------------------------------------------------------
#pragma mark - Column Vector Multiplication


void SparMatSymBlkDiag::Column::vecMulAdd1D(const real* X, real* Y, size_t jj) const
{
#if ( BLOCK_SIZE == 1 )
    const real X0 = X[jj];
    real D = dia_.value();
    real Y0 = Y[jj] + D * X0;
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        const real M = blk_[n].value();
        Y[ii] += M * X0;
        Y0 += M * X[ii];
    }
    Y[jj] = Y0;
#endif
}


void SparMatSymBlkDiag::Column::vecMulAdd2D(const real* X, real* Y, size_t jj) const
{
#if ( BLOCK_SIZE == 2 )
    const Vector2 xx(X+jj);
    assert_small(dia_.asymmetry());
    Vector2 yy = dia_.vecmul(xx);
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        Block const& M = blk_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.add_to(Y+jj);
#endif
}

void SparMatSymBlkDiag::Column::vecMulAdd3D(const real* X, real* Y, size_t jj) const
{
#if ( BLOCK_SIZE == 3 )
    const Vector3 xxx(X+jj);
    assert_small(dia_.asymmetry());
    Vector3 yyy = dia_.vecmul(xxx);
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        Block const& M = blk_[n];
        M.vecmul(xxx).add_to(Y+ii);
        yyy += M.trans_vecmul(X+ii);
    }
    yyy.add_to(Y+jj);
#endif
}


void SparMatSymBlkDiag::Column::vecMulAdd4D(const real* X, real* Y, size_t jj) const
{
#if ( BLOCK_SIZE == 4 )
    const vec4 xxxx = load4(X+jj);
    assert_small(dia_.asymmetry());
    vec4 yyyy = dia_.vecmul4(xxxx);
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        Block const& M = blk_[n];
        store4(Y+ii, add4(load4(Y+ii), M.vecmul4(xxxx)));
        yyyy += M.trans_vecmul(X+ii);
    }
    store4(Y+jj, add4(yyyy, load4(Y+jj)));
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Single Precision Optimized Vector Multiplication

#if ( BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && defined(__SSE3__)
void SparMatSymBlkDiag::Column::vecMulAdd3D_SSE(const float* X, float* Y, size_t jj) const
{
    // load 3x3 matrix diagonal element into 3 vectors:
    float const* D = dia_;
    
    //multiply with the symmetrized block, assuming it has been symmetrized:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    const vec4f tt = loadu4f(X+jj);
# if ( BLD == 4 )
    vec4f s0 = mul4f(streamload4f(D  ), tt);
    vec4f s1 = mul4f(streamload4f(D+4), tt);
    vec4f s2 = mul4f(streamload4f(D+8), tt);
# else
    vec4f s0 = mul4f(load3fZ(D      ), tt);
    vec4f s1 = mul4f(load3fZ(D+BLD  ), tt);
    vec4f s2 = mul4f(load3fZ(D+BLD*2), tt);
# endif
    const vec4f x0 = permute4f(tt, 0x00);
    const vec4f x1 = permute4f(tt, 0x55);
    const vec4f x2 = permute4f(tt, 0xAA);
    
    Block  const* blk = blk_;
    size_t const* inx = inx_;

    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( size_t n = 0; n < size_; ++n )
    {
        float const* M = *blk++;
        const size_t ii = *inx++;
# if ( BLD == 4 )
        const vec4f M012 = streamload4f(M  );
        const vec4f M345 = streamload4f(M+4);
        const vec4f M678 = streamload4f(M+8);
# else
        const vec4f M012 = load3fZ(M      );
        const vec4f M345 = load3fZ(M+BLD  );
        const vec4f M678 = load3fZ(M+BLD*2);
# endif
        // multiply with the full block:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec4f z = fmadd4f(M012, x0, loadu4f(Y+ii));
        z = fmadd4f(M345, x1, z);
        z = fmadd4f(M678, x2, z);
        storeu4f(Y+ii, z);
        
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
        s0 = fmadd4f(M012, xyz, s0);
        s1 = fmadd4f(M345, xyz, s1);
        s2 = fmadd4f(M678, xyz, s2);
    }
    /* finally sum horizontally:
     s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
     to { Y0+Y0+Y0, Y1+Y1+Y1, Y2+Y2+Y2, 0 }
     */
    vec4f s3 = setzero4f();
    s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
    s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
    s0 = add4f(shuffle4f(s0, s2, 0x4E), shuffle4f(s0, s2, 0xE4));
    storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
}
#endif


#if ( BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && defined(__SSE3__)
void SparMatSymBlkDiag::Column::vecMulAdd3D_SSEU(const float* X, float* Y, size_t jj) const
{
    assert_small(dia_.asymmetry());
    //std::cout << dia_.to_string(7,1); printf(" MSSB %lu : %lu\n", jj, size_);
    // load 3x3 matrix diagonal element into 3 vectors:
    float const* D = dia_;

    //multiply with the diagonal block, assuming it has been symmetrized:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    const vec4f tt = loadu4f(X+jj);
    const vec4f x0 = permute4f(tt, 0x00);
    const vec4f x1 = permute4f(tt, 0x55);
    const vec4f x2 = permute4f(tt, 0xAA);

    Block  const* blk = blk_;
    size_t const* inx = inx_;

    if ( size_ > 0 )
    {
# if ( BLD == 4 )
        vec4f s0 = mul4f(streamload4f(D  ), tt);
        vec4f s1 = mul4f(streamload4f(D+4), tt);
        vec4f s2 = mul4f(streamload4f(D+8), tt);
# else
        vec4f s0 = mul4f(load3fZ(D      ), tt);
        vec4f s1 = mul4f(load3fZ(D+BLD  ), tt);
        vec4f s2 = mul4f(load3fZ(D+BLD*2), tt);
# endif
        size_t n = 0;
        {
            const size_t end = 2 * (size_/2);
            // process 2 by 2
            #pragma nounroll
            for ( ; n < end; n += 2 )
            {
                float const* M = *blk++;
                float const* P = *blk++;
                const size_t ii = *inx++;
                const size_t kk = *inx++;
                assert_true( ii < kk );
# if ( BLD == 4 )
                const vec4f M012 = streamload4f(M  );
                const vec4f M345 = streamload4f(M+4);
                const vec4f M678 = streamload4f(M+8);
                const vec4f P012 = streamload4f(P  );
                const vec4f P345 = streamload4f(P+4);
                const vec4f P678 = streamload4f(P+8);
# else
                const vec4f M012 = load3fZ(M      );
                const vec4f M345 = load3fZ(M+BLD  );
                const vec4f M678 = load3fZ(M+BLD*2);
                const vec4f P012 = load3fZ(P      );
                const vec4f P345 = load3fZ(P+BLD  );
                const vec4f P678 = load3fZ(P+BLD*2);
# endif
                // multiply with the full block:
                vec4f z = fmadd4f(M012, x0, loadu4f(Y+ii));
                vec4f t = fmadd4f(P012, x0, loadu4f(Y+kk));
                vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
                vec4f tuv = loadu4f(X+kk);  // xyz = { X0 X1 X2 - }
                z = fmadd4f(M345, x1, z);
                t = fmadd4f(P345, x1, t);
                s0 = fmadd4f(M012, xyz, s0);
                s1 = fmadd4f(M345, xyz, s1);
                s2 = fmadd4f(M678, xyz, s2);
                z = fmadd4f(M678, x2, z);
                t = fmadd4f(P678, x2, t);
                s0 = fmadd4f(P012, tuv, s0);
                s1 = fmadd4f(P345, tuv, s1);
                s2 = fmadd4f(P678, tuv, s2);
                storeu4f(Y+ii, z);
                storeu4f(Y+kk, t);
            }
        }
        // process remaining blocks
        #pragma nounroll
        for ( ; n < size_; ++n )
        {
            float const* M = *blk++;
            const size_t ii = *inx++;
# if ( BLD == 4 )
            const vec4f M012 = streamload4f(M  );
            const vec4f M345 = streamload4f(M+4);
            const vec4f M678 = streamload4f(M+8);
# else
            const vec4f M012 = load3fZ(M      );
            const vec4f M345 = load3fZ(M+BLD  );
            const vec4f M678 = load3fZ(M+BLD*2);
# endif
            // multiply with the full block:
            //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
            //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
            //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
            vec4f z = fmadd4f(M012, x0, loadu4f(Y+ii));
            z = fmadd4f(M345, x1, z);
            z = fmadd4f(M678, x2, z);
            storeu4f(Y+ii, z);
            
            // multiply with the transposed block:
            //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
            //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
            //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
            vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
            s0 = fmadd4f(M012, xyz, s0);
            s1 = fmadd4f(M345, xyz, s1);
            s2 = fmadd4f(M678, xyz, s2);
        }
        
        /* finally sum horizontally:
         s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
         to { Y0+Y0+Y0, Y1+Y1+Y1, Y2+Y2+Y2, 0 }
         */
        vec4f s3 = setzero4f();
        s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
        s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
        s0 = add4f(shuffle4f(s0, s2, 0x4E), shuffle4f(s0, s2, 0xE4));
        storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
    }
    else
    {
# if ( BLD == 4 )
        vec4f s0 = mul4f(streamload4f(D  ), x0);
        vec4f s1 = mul4f(streamload4f(D+4), x1);
        vec4f s2 = mul4f(streamload4f(D+8), x2);
# else
        vec4f s0 = mul4f(load3fZ(D      ), x0);
        vec4f s1 = mul4f(load3fZ(D+BLD  ), x1);
        vec4f s2 = mul4f(load3fZ(D+BLD*2), x2);
# endif
        storeu4f(Y+jj, add4f(add4f(loadu4f(Y+jj), s0), add4f(s1, s2)));
    }
}
#endif

/**
 Only process off-diagonal terms!
 */
#if ( BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && defined(__SSE3__)
void SparMatSymBlkDiag::Column::vecMulAddTriangle3D_SSE(const float* X, float* Y, size_t jj) const
{
    assert_true(size_ > 0);
    vec4f tt = loadu4f(X+jj);
    
    const vec4f x0 = permute4f(tt, 0x00);
    const vec4f x1 = permute4f(tt, 0x55);
    const vec4f x2 = permute4f(tt, 0xAA);
    
    vec4f s0 = setzero4f();
    vec4f s1 = setzero4f();
    vec4f s2 = setzero4f();

    Block  const* blk = blk_;
    size_t const* inx = inx_;
    size_t const* end = inx_ + size_;
    
    //size_t n = 0;
    {
        // process 2 by 2
        size_t const* stop = end - 1; // inx+1 <= end-1 is inx < end-1
        #pragma nounroll
        for ( ; inx < stop; inx += 2 )
        {
            float const* M = *blk++;
            float const* P = *blk++;
            const size_t ii = inx[0];
            const size_t kk = inx[1];
            assert_true( ii < kk );
# if ( BLD == 4 )
            const vec4f M012 = streamload4f(M  );
            const vec4f M345 = streamload4f(M+4);
            const vec4f M678 = streamload4f(M+8);
            const vec4f P012 = streamload4f(P  );
            const vec4f P345 = streamload4f(P+4);
            const vec4f P678 = streamload4f(P+8);
# else
            const vec4f M012 = load3fZ(M      );
            const vec4f M345 = load3fZ(M+BLD  );
            const vec4f M678 = load3fZ(M+BLD*2);
            const vec4f P012 = load3fZ(P      );
            const vec4f P345 = load3fZ(P+BLD  );
            const vec4f P678 = load3fZ(P+BLD*2);
# endif
            // multiply with the full block:
            vec4f z = fmadd4f(M012, x0, loadu4f(Y+ii));
            vec4f t = fmadd4f(P012, x0, loadu4f(Y+kk));
            vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
            vec4f tuv = loadu4f(X+kk);  // xyz = { X0 X1 X2 - }
            z = fmadd4f(M345, x1, z);
            t = fmadd4f(P345, x1, t);
            s0 = fmadd4f(M012, xyz, s0);
            s1 = fmadd4f(M345, xyz, s1);
            s2 = fmadd4f(M678, xyz, s2);
            z = fmadd4f(M678, x2, z);
            t = fmadd4f(P678, x2, t);
            s0 = fmadd4f(P012, tuv, s0);
            s1 = fmadd4f(P345, tuv, s1);
            s2 = fmadd4f(P678, tuv, s2);
            storeu4f(Y+ii, z);
            storeu4f(Y+kk, t);
        }
    }
    //#pragma nounroll // for ( ; n < size_; ++n )
    // process last remaining block
    if ( inx < end )
    {
        const size_t ii = *inx;
        assert_true( blk == blk_+size_-1 );
        float const* L = *blk; // last block
# if ( BLD == 4 )
        const vec4f L012 = streamload4f(L  );
        const vec4f L345 = streamload4f(L+4);
        const vec4f L678 = streamload4f(L+8);
# else
        const vec4f L012 = load3fZ(L      );
        const vec4f L345 = load3fZ(L+BLD  );
        const vec4f L678 = load3fZ(L+BLD*2);
# endif
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        vec4f xyz = loadu4f(X+ii);  // xyz = { X0 X1 X2 - }
        s0 = fmadd4f(L012, xyz, s0);
        s1 = fmadd4f(L345, xyz, s1);
        s2 = fmadd4f(L678, xyz, s2);

        // multiply with the full block:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec4f z = fmadd4f(L012, x0, loadu4f(Y+ii));
        z = fmadd4f(L345, x1, z);
        z = fmadd4f(L678, x2, z);
        storeu4f(Y+ii, z);
    }
    
    /* finally sum horizontally:
     s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
     to { Y0+Y0+Y0, Y1+Y1+Y1, Y2+Y2+Y2, 0 }
     */
    tt = setzero4f();
    s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
    s2 = add4f(unpacklo4f(s2, tt), unpackhi4f(s2, tt));
    s0 = add4f(shuffle4f(s0, s2, 0x4E), shuffle4f(s0, s2, 0xE4));
    storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Double Precision Optimized Vector Multiplication

#if ( BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && defined(__SSE3__)
void SparMatSymBlkDiag::Column::vecMulAdd2D_SSE(const double* X, double* Y, size_t jj) const
{
    vec2 x0, x1;
    vec2 yy = load2(Y+jj);
    {
        //const real X0 = X[jj  ];
        //const real X1 = X[jj+1];
        vec2 xx = load2(X+jj);
        x0 = unpacklo2(xx, xx);
        x1 = unpackhi2(xx, xx);
        
        // load 2x2 matrix element into 2 vectors:
        double const* M = dia_;
        //assume the block is already symmetrized:
        // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1;
        // Y1 = Y[jj+1] + M[1] * X0 + M[3] * X1;
        xx = add2(mul2(load2(M  ), x0), yy);
        yy = add2(mul2(load2(M+2), x1), xx);
    }
    
    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        vec2 xx = load2(X+ii);
        
        // load 2x2 matrix element into 2 vectors:
        double const* M = blk_[n];
        vec2 m01 = load2(M);
        vec2 m23 = load2(M+2);
        
        // multiply with the full block:
        //Y[ii  ] += M[0] * X0 + M[2] * X1;
        //Y[ii+1] += M[1] * X0 + M[3] * X1;
        vec2 mx0 = add2(mul2(m01, x0), load2(Y+ii));
        mx0 = add2(mul2(m23, x1), mx0);
        store2(Y+ii, mx0);

        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
        //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
        vec2 mxx = mul2(m01, xx);
        vec2 myy = mul2(m23, xx);
        yy = add2(add2(unpacklo2(mxx, myy), unpackhi2(mxx, myy)), yy);
    }
    //Y[jj  ] = Y0;
    //Y[jj+1] = Y1;
    store2(Y+jj, yy);
}
#endif


#if ( BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAdd2D_AVX(const double* X, double* Y, size_t jj) const
{
    // xy = { X0 X1 X0 X1 }
    vec4 xy = broadcast2(X+jj);
    //multiply with full block, assuming it is symmetric:
    // Y0 = M[0] * X0 + M[1] * X1;
    // Y1 = M[1] * X0 + M[3] * X1;
    
    // yyyy = { Y0 Y0 Y1 Y1 }
    // load 2x2 matrix element into 2 vectors:
    vec4 ss = mul4(streamload4(dia_), xy);

    //const real X0 = X[jj  ];
    //const real X1 = X[jj+1];
    // xxyy = { X0 X0 X1 X1 }
    const vec4 xxyy = permute4(xy, 0b1100);

    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        vec4 mat = streamload4(blk_[n]); // load 2x2 matrix
        vec4 yy = load2Z(Y+ii);          // yy = { Y0 Y1 0 0 }
        vec4 xx = broadcast2(X+ii);      // xx = { X0 X1 X0 X1 }

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


#if ( BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
inline static void multiply2D(double const* X, double* Y, size_t ii, vec4 const& mat, vec4 const& xxxx, vec4& ss)
{
    vec4 xx = broadcast2(X+ii);
    vec4 u = fmadd4(mat, xxxx, load2Z(Y+ii));
    store2(Y+ii, add2(getlo(u), gethi(u)));
    ss = fmadd4(mat, xx, ss);
}
#endif


#if ( BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAdd2D_AVXU(const double* X, double* Y, size_t jj) const
{
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(streamload4(dia_), xyxy);
    const vec4 xxyy = permute4(xyxy, 0b1100);
    vec4 s1 = setzero4();
    
    Block  const* blk = blk_;
    size_t const* inx = inx_;

    size_t n = 0;
    const size_t end = 2 * (size_/2);
    // process 2 by 2:
    #pragma nounroll
    for ( ; n < end; n += 2 )
    {
#if ( 0 )
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, ss);
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, s1);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const size_t i0 = *inx++;
        const size_t i1 = *inx++;
        assert_true( i0 < i1 );
        vec4 mat0 = streamload4(*blk++);
        vec4 mat1 = streamload4(*blk++);
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
    #pragma nounroll
    for ( ; n < size_; ++n )
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, ss);
    /* finally horizontally sum ss = { SX SX SY SY } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


#if ( BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAdd2D_AVXUU(const double* X, double* Y, size_t jj) const
{
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(streamload4(dia_), xyxy);
    const vec4 xxyy = permute4(xyxy, 0b1100);
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();

    Block  const* blk = blk_;
    size_t const* inx = inx_;

    size_t n = 0;
    const size_t end = 4 * (size_/4);
    // process 4 by 4:
    #pragma nounroll
    for ( ; n < end; n += 4 )
    {
#if ( 0 )
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, ss);
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, s1);
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, s2);
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const size_t i0 = *inx++;
        const size_t i1 = *inx++;
        const size_t i2 = *inx++;
        const size_t i3 = *inx++;
        assert_true( i0 < i1 );
        assert_true( i1 < i2 );
        assert_true( i2 < i3 );
        vec4 mat0 = streamload4(*blk++);
        vec4 mat1 = streamload4(*blk++);
        vec4 mat2 = streamload4(*blk++);
        vec4 mat3 = streamload4(*blk++);
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
    #pragma nounroll
    for ( ; n < size_; ++n )
        multiply2D(X, Y, *inx++, streamload4(*blk++), xxyy, ss);
    /* finally sum ss = { S0 S0 S1 S1 } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


#if ( BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAdd3D_AVX(const double* X, double* Y, size_t jj) const
{
    // load 3x3 matrix diagonal element into 3 vectors:
    double const* D = dia_;
    
    //multiply with the symmetrized block, assuming it has been symmetrized:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    vec4 s0, s1, s2;
    vec4 x0, x1, x2;
    {
        vec4 tt = loadu4(X+jj);
# if ( BLD == 4 )
        s0 = mul4(streamload4(D  ), tt);
        s1 = mul4(streamload4(D+4), tt);
        s2 = mul4(streamload4(D+8), tt);
# else
        s0 = mul4(load3(D      ), tt);
        s1 = mul4(load3(D+BLD  ), tt);
        s2 = mul4(load3(D+BLD*2), tt);
# endif
        // prepare broacasted vectors
        vec4 p = swap2f128(tt);
        vec4 l = blend22(tt, p);
        vec4 u = blend22(p, tt);
        x0 = duplo4(l);
        x1 = duphi4(l);
        x2 = duplo4(u);
    }
#if 0
    // alternative strategy:
    const vec4 x0 = broadcast1(X+jj);
    const vec4 x1 = broadcast1(X+jj+1);
    const vec4 x2 = broadcast1(X+jj+2);
#endif
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        double const* M = blk_[n];
# if ( BLD == 4 )
        const vec4 M012 = streamload4(M  );
        const vec4 M345 = streamload4(M+4);
        const vec4 M678 = streamload4(M+8);
# else
        const vec4 M012 = load3(M      );
        const vec4 M345 = load3(M+BLD  );
        const vec4 M678 = load3(M+BLD*2);
# endif
        // multiply with the full block:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec4 z = fmadd4(M012, x0, loadu4(Y+ii));
        z = fmadd4(M345, x1, z);
        z = fmadd4(M678, x2, z);
        storeu4(Y+ii, z);
        
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        vec4 xyz = loadu4(X+ii);  // xyz = { X0 X1 X2 - }
        s0 = fmadd4(M012, xyz, s0);
        s1 = fmadd4(M345, xyz, s1);
        s2 = fmadd4(M678, xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
#if ( 0 )
    Y[jj  ] += s0[0] + s0[1] + s0[2];
    Y[jj+1] += s1[0] + s1[1] + s1[2];
    Y[jj+2] += s2[0] + s2[1] + s2[2];
#else
    vec4 s3 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    s0 = add4(catshift2(s0, s2), blend22(s0, s2));
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
#endif
}
#endif


#if ( BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAdd3D_AVXU(const double* X, double* Y, size_t jj) const
{
    vec4 s0, s1, s2;
    vec4 x0, x1, x2;
    // load 3x3 matrix element into 3 vectors:
    {
        const real* D = dia_;
        assert_small(dia_.asymmetry());
        vec4 tt = loadu4(X+jj);
        // multiply by diagonal elements:
        s0 = mul4(streamload4(D  ), tt);
        s1 = mul4(streamload4(D+4), tt);
        s2 = mul4(streamload4(D+8), tt);
        // prepare broadcasted vectors:
        vec4 p = swap2f128(tt);
        vec4 l = blend22(tt, p);
        vec4 u = blend22(p, tt);
        x0 = duplo4(l);
        x1 = duphi4(l);
        x2 = duplo4(u);
    }
    if ( size_ > 0 )
    {
        vec4 t0 = setzero4();
        vec4 t1 = setzero4();
        vec4 t2 = setzero4();
        // There is a dependency in the loop for 's0', 's1' and 's2'.
        const real* M = blk_[0];
        const size_t* inx = inx_;
        const real* end = blk_[2*(size_/2)];
        /*
         Unrolling will reduce the dependency chain, which may be limiting the
         throughput here. However the number of registers (16 for AVX CPU) limits
         the level of unrolling that can be done.
         */
        //process 2 by 2:
#pragma nounroll
        for ( ; M < end; M += 2*SB )
        {
            const size_t i0 = inx[0];
            const size_t i1 = inx[1];
            assert_true( i0 < i1 );
            inx += 2;
            //printf("--- %4i %4i\n", i0, i1);
            vec4 ma0 = streamload4(M);
            vec4 ma1 = streamload4(M+SB);
            vec4 z0 = fmadd4(ma0, x0, loadu4(Y+i0));
            vec4 z1 = fmadd4(ma1, x0, loadu4(Y+i1));
            vec4 xyz0 = loadu4(X+i0);
            vec4 xyz1 = loadu4(X+i1);
            s0 = fmadd4(ma0, xyz0, s0);
            t0 = fmadd4(ma1, xyz1, t0);
            // multiply with the full block:
            vec4 mb0 = streamload4(M+4);
            vec4 mb1 = streamload4(M+(SB+4));
            z0 = fmadd4(mb0, x1, z0);
            z1 = fmadd4(mb1, x1, z1);
            s1 = fmadd4(mb0, xyz0, s1);
            t1 = fmadd4(mb1, xyz1, t1);
            vec4 mc0 = streamload4(M+8);
            vec4 mc1 = streamload4(M+(SB+8));
            z0 = fmadd4(mc0, x2, z0);
            z1 = fmadd4(mc1, x2, z1);
            s2 = fmadd4(mc0, xyz0, s2);
            t2 = fmadd4(mc1, xyz1, t2);
            /*
             Attention: the 4th elements of the vectors z0 and z1 would be correct,
             because only zero was added to the value loaded from 'Y'. However, in the
             case where the indices i0 and i1 are consecutive and reverted (i1 < i0),
             the value stored in z0 would not have been updated giving a wrong results.
             The solution is to either use a 'store3(Y+i1, z1)', or to make sure that
             indices are non-consecutive or ordered in the column in increasing order.
             This affects performance since 'store3' is slower than 'storeu4'
             */
            storeu4(Y+i0, z0);
            storeu4(Y+i1, z1);
        }
        s0 = add4(s0, t0);
        s1 = add4(s1, t1);
        s2 = add4(s2, t2);
        
        // process remaining blocks:
        end = blk_[size_];
        #pragma nounroll
        for ( ; M < end; M += SB )
        {
            const size_t ii = inx[0];
            ++inx;
            //printf("--- %4i\n", ii);
            vec4 ma = streamload4(M);
            vec4 z = fmadd4(ma, x0, loadu4(Y+ii));
            vec4 xyz = loadu4(X+ii);
            s0 = fmadd4(ma, xyz, s0);
            
            vec4 mb = streamload4(M+4);
            z = fmadd4(mb, x1, z);
            s1 = fmadd4(mb, xyz, s1);
            
            vec4 mc = streamload4(M+8);
            z = fmadd4(mc, x2, z);
            s2 = fmadd4(mc, xyz, s2);
            storeu4(Y+ii, z);
        }
    }
    // finally sum s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
    x0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, x0), unpackhi4(s2, x0));
    s0 = add4(catshift2(s0, s2), blend22(s0, s2));
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
}
#endif



#if ( BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAddTriangle3D_AVX(const double* X, double* Y, size_t jj) const
{
    vec4 x0, x1, x2;
    // load 3x3 matrix element into 3 vectors:
    {
        vec4 tt = loadu4(X+jj);
        // prepare broadcasted vectors:
        vec4 p = swap2f128(tt);
        vec4 l = blend22(tt, p);
        vec4 u = blend22(p, tt);
        x0 = duplo4(l);
        x1 = duphi4(l);
        x2 = duplo4(u);
    }
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    if ( size_ > 0 )
    {
        vec4 t0 = setzero4();
        vec4 t1 = setzero4();
        vec4 t2 = setzero4();
        // There is a dependency in the loop for 's0', 's1' and 's2'.
        const real* M = blk_[0];
        const size_t* inx = inx_;
        const real* end = blk_[2*(size_/2)];
        /*
         Unrolling will reduce the dependency chain, which may be limiting the
         throughput here. However the number of registers (16 for AVX CPU) limits
         the level of unrolling that can be done.
         */
        //process 2 by 2:
#pragma nounroll
        for ( ; M < end; M += 2*SB )
        {
            const size_t i0 = inx[0];
            const size_t i1 = inx[1];
            assert_true( i0 < i1 );
            inx += 2;
            //printf("--- %4i %4i\n", i0, i1);
            vec4 ma0 = streamload4(M);
            vec4 ma1 = streamload4(M+SB);
            vec4 z0 = fmadd4(ma0, x0, loadu4(Y+i0));
            vec4 z1 = fmadd4(ma1, x0, loadu4(Y+i1));
            vec4 xyz0 = loadu4(X+i0);
            vec4 xyz1 = loadu4(X+i1);
            s0 = fmadd4(ma0, xyz0, s0);
            t0 = fmadd4(ma1, xyz1, t0);
            // multiply with the full block:
            vec4 mb0 = streamload4(M+4);
            vec4 mb1 = streamload4(M+(SB+4));
            z0 = fmadd4(mb0, x1, z0);
            z1 = fmadd4(mb1, x1, z1);
            s1 = fmadd4(mb0, xyz0, s1);
            t1 = fmadd4(mb1, xyz1, t1);
            vec4 mc0 = streamload4(M+8);
            vec4 mc1 = streamload4(M+(SB+8));
            z0 = fmadd4(mc0, x2, z0);
            z1 = fmadd4(mc1, x2, z1);
            s2 = fmadd4(mc0, xyz0, s2);
            t2 = fmadd4(mc1, xyz1, t2);
            /*
             Attention: the 4th elements of the vectors z0 and z1 would be correct,
             because only zero was added to the value loaded from 'Y'. However, in the
             case where the indices i0 and i1 are consecutive and reverted (i1 < i0),
             the value stored in z0 would not have been updated giving a wrong results.
             The solution is to either use a 'store3(Y+i1, z1)', or to make sure that
             indices are non-consecutive or ordered in the column in increasing order.
             This affects performance since 'store3' is slower than 'storeu4'
             */
            storeu4(Y+i0, z0);
            storeu4(Y+i1, z1);
        }
        s0 = add4(s0, t0);
        s1 = add4(s1, t1);
        s2 = add4(s2, t2);
        
        // process remaining block:
        end = blk_[size_];
        #pragma nounroll
        for ( ; M < end; M += SB )
        {
            const size_t ii = inx[0];
            ++inx;
            //printf("--- %4i\n", ii);
            vec4 ma = streamload4(M);
            vec4 z = fmadd4(ma, x0, loadu4(Y+ii));
            vec4 xyz = loadu4(X+ii);
            s0 = fmadd4(ma, xyz, s0);
            
            vec4 mb = streamload4(M+4);
            z = fmadd4(mb, x1, z);
            s1 = fmadd4(mb, xyz, s1);
            
            vec4 mc = streamload4(M+8);
            z = fmadd4(mc, x2, z);
            s2 = fmadd4(mc, xyz, s2);
            storeu4(Y+ii, z);
        }
    }
    // finally sum s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
    x0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, x0), unpackhi4(s2, x0));
    s0 = add4(catshift2(s0, s2), blend22(s0, s2));
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
}
#endif


#if ( BLOCK_SIZE == 4 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Column::vecMulAdd4D_AVX(const double* X, double* Y, size_t jj) const
{
    double const* D = dia_;
    //multiply with the symmetrized block, assuming it has been symmetrized:
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    vec4 s0, s1, s2;
    {
        vec4 tt = load4(X+jj);
        s0 = mul4(streamload4(D   ), tt);
        s1 = mul4(streamload4(D+4 ), tt);
        s2 = mul4(streamload4(D+8 ), tt);
        s3 = mul4(streamload4(D+12), tt);
    }
    // sum non-diagonal elements:
#if ( 0 )
    const vec4 x0 = broadcast1(X+jj);
    const vec4 x1 = broadcast1(X+jj+1);
    const vec4 x2 = broadcast1(X+jj+2);
    const vec4 x3 = broadcast1(X+jj+3);
#else
    vec4 x0, x1, x2, x3;
    {
        x1 = duplo2f128(tt);
        x3 = duphi2f128(tt);
        x0 = duplo4(x1);
        x1 = duphi4(x1);
        x2 = duplo4(x3);
        x3 = duphi4(x3);
    }
#endif
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( size_t n = 0; n < size_; ++n )
    {
        const size_t ii = inx_[n];
        double const* M = blk_[n];
        const vec4 yy = load4(Y+ii);
        const vec4 xyzt = load4(X+ii);  // xyzt = { X0 X1 X2 X3 }
        const vec4 m0 = streamload4(M);
        vec4 z = fmadd4(m0, x0, yy);
        s0 = fmadd4(m0, xyzt, s0);
        
        const vec4 m1 = streamload4(M+4);
        z  = fmadd4(m1, x1, z);
        s1 = fmadd4(m1, xyzt, s1);

        const vec4 m2 = streamload4(M+8);
        z  = fmadd4(m2, x2, z);
        s2 = fmadd4(m2, xyzt, s2);

        const vec4 m3 = streamload4(M+12);
        z  = fmadd4(m3, x3, z);
        s3 = fmadd4(m3, xyzt, s3);
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

#if SMSBD_USES_AVX && REAL_IS_DOUBLE
#   define VECMULADD2D vecMulAdd2D_AVXU
#   define VECMULADD3D vecMulAdd3D_AVXU
#   define VECMULADD4D vecMulAdd4D_AVX
#elif SMSBD_USES_SSE && REAL_IS_DOUBLE
#   define VECMULADD2D vecMulAdd2D_SSE
#   define VECMULADD3D vecMulAdd3D
#   define VECMULADD4D vecMulAdd4D
#elif SMSBD_USES_SSE
#   define VECMULADD2D vecMulAdd2D
#   define VECMULADD3D vecMulAdd3D_SSEU
#   define VECMULADD4D vecMulAdd4D
#else
#   define VECMULADD2D vecMulAdd2D
#   define VECMULADD3D vecMulAdd3D
#   define VECMULADD4D vecMulAdd4D
#endif


// multiplication of a vector: Y = Y + M * X
void SparMatSymBlkDiag::vecMulAdd(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_) / BLOCK_SIZE;
    for ( size_t jj = start/BLOCK_SIZE; jj < stop; ++jj )
    {
        //std::clog << "SparMatSymBlkDiag column " << jj << "  " << size_ << " \n";
#if ( BLOCK_SIZE == 1 )
        pilar_[jj].vecMulAdd1D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 2 )
        pilar_[jj].VECMULADD2D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 3 )
        pilar_[jj].VECMULADD3D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 4 )
        pilar_[jj].VECMULADD4D(X, Y, jj*BLOCK_SIZE);
#endif
    }
}


// multiplication of a vector: Y = Y + M * X
void SparMatSymBlkDiag::vecMulAdd_ALT(const real* X, real* Y) const
{
    size_t stop = size_ / BLOCK_SIZE;
    for ( size_t jj = 0; jj < stop; ++jj )
    {
        //std::clog << "SparMatSymBlkDiag column " << jj << "  " << size_ << " \n";
#if ( BLOCK_SIZE == 1 )
        pilar_[jj].vecMulAdd1D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 2 )
        pilar_[jj].vecMulAdd2D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 3 )
        pilar_[jj].vecMulAdd3D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 4 )
        pilar_[jj].vecMulAdd4D(X, Y, jj*BLOCK_SIZE);
#endif
    }
}


// multiplication of a vector: Y = Y + M * X
void SparMatSymBlkDiag::vecMulAdd_TIME(const real* X, real* Y) const
{
    size_t stop = size_ / BLOCK_SIZE;
    size_t cnt = 0, col = 0;
    //auto rdt = __rdtscd();
    for ( size_t jj = 0; jj < stop; ++jj )
    {
        col++;
        cnt += pilar_[jj].size_;
        //std::clog << "SparMatSymBlkDiag column " << jj << "  " << size_ << " \n";
#if ( BLOCK_SIZE == 1 )
        pilar_[jj].vecMulAdd1D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 2 )
        pilar_[jj].vecMulAdd2D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 3 )
        pilar_[jj].vecMulAdd3D(X, Y, jj*BLOCK_SIZE);
#elif ( BLOCK_SIZE == 4 )
        pilar_[jj].vecMulAdd4D(X, Y, jj*BLOCK_SIZE);
#endif
    }
    /*
    if ( cnt > 0 )
        fprintf(stderr, "SMSBD %6lu rows %6lu blocks  cycles/block: %5.2f\n",\
                col, cnt, real(__rdtscd()-rdt)/cnt);
     */
}

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication

#if ( BLOCK_SIZE == 4 ) && REAL_IS_DOUBLE && defined(__AVX__)
void SparMatSymBlkDiag::vecMulDiagonal3D(const double* src, double* dst) const
{
    const size_t stop = size_ / BLOCK_SIZE;
    #pragma ivdep unroll (4)
    #pragma clang loop unroll_count(4)
    for ( size_t j = 0; j < stop; ++j )
    {
        double const* M = pilar_[j].dia_;
#if 1
        vec4 x0 = loadu4(src);
        vec4 x2 = swap2f128(x0);
        vec4 x1 = blend22(x0, x2);
        x2 = blend22(x2, x0);
        x0 = duplo4(x1);
        x1 = duphi4(x1);
        x2 = duplo4(x2);
#else
        vec4 x2 = broadcast2(src);
        vec4 x0 = unpacklo4(x2, x2);
        vec4 x1 = unpackhi4(x2, x2);
        x2 = broadcast1(src+2);
#endif
        src += BLOCK_SIZE;
        //multiply with the diagonal block:
        // Y0 = M[0] * X0 + M[3] * X1 + M[6] * X2;
        // Y1 = M[1] * X0 + M[4] * X1 + M[7] * X2;
        // Y2 = M[2] * X0 + M[5] * X1 + M[8] * X2;
# if ( BLD == 4 )
        x0 =   mul4(streamload4(M), x0);
        x0 = fmadd4(streamload4(M+4), x1, x0);
        x0 = fmadd4(streamload4(M+8), x2, x0);
# else
        x0 =   mul4(load3(M), x0);
        x0 = fmadd4(load3(M+BLD  ), x1, x0);
        x0 = fmadd4(load3(M+BLD*2), x2, x0);
# endif
        storeu4(dst, x0);
        dst += BLOCK_SIZE;
    }
}
#endif



#if ( BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
void SparMatSymBlkDiag::vecMulDiagonal3D(const double* src, double* dst) const
{
    size_t j = 0;
    if ( size_ & 1 )
    {
        double const* M = pilar_[0].dia_;
        vec4 x0 = broadcast1(src  );
        vec4 x1 = broadcast1(src+1);
        vec4 x2 = broadcast1(src+2);
        src += BLOCK_SIZE;
        //multiply with the diagonal block:
        // Y0 = M[0] * X0 + M[3] * X1 + M[6] * X2;
        // Y1 = M[1] * X0 + M[4] * X1 + M[7] * X2;
        // Y2 = M[2] * X0 + M[5] * X1 + M[8] * X2;
# if ( BLD == 4 )
        x0 = mul4(streamload4(M), x0);
        x0 = fmadd4(streamload4(M+4), x1, x0);
        x0 = fmadd4(streamload4(M+8), x2, x0);
# else
        x0 = mul4(load3(M), x0);
        x0 = fmadd4(load3(M+BLD  ), x1, x0);
        x0 = fmadd4(load3(M+BLD*2), x2, x0);
# endif
        storeu4(dst, x0);
        dst += BLOCK_SIZE;
        ++j;
    }

    const size_t stop = size_ / BLOCK_SIZE;
    #pragma ivdep
    #pragma clang loop unroll_count(2)
    while ( j < stop )
    {
        double const* M = pilar_[j++].dia_;
        double const* N = pilar_[j++].dia_;
        //broadcast the source vectors:
        vec4 x0 = broadcast1(src);
        vec4 x1 = broadcast1(src+1);
        vec4 x2 = broadcast1(src+2);
        vec4 x3 = broadcast1(src+3);
        vec4 x4 = broadcast1(src+4);
        vec4 x5 = broadcast1(src+5);
        src += 2*BLOCK_SIZE;
        //multiply with the matrix block:
        x0 = mul4(streamload4(M  ), x0);
        x1 = mul4(streamload4(M+4), x1);
        x2 = mul4(streamload4(M+8), x2);
        x3 = mul4(streamload4(N  ), x3);
        x4 = mul4(streamload4(N+4), x4);
        x5 = mul4(streamload4(N+8), x5);
        storeu4(dst,   add4(add4(x0, x1), x2));
        storeu4(dst+3, add4(add4(x3, x4), x5));
        dst += 2*BLOCK_SIZE;
    }
}
#endif


#if ( BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && defined(__SSE3__)
void SparMatSymBlkDiag::vecMulDiagonal3D(const float* X, float* Y) const
{
    size_t stop = size_ / BLOCK_SIZE;
    #pragma unroll (4)
    for ( size_t j = 0; j < stop; ++j )
    {
        float const* D = pilar_[j].dia_;
        //multiply with the diagonal block, assuming it has been symmetrized:
        // Y0 = M[0] * X0 + M[1] * X1 + M[2] * X2;
        // Y1 = M[1] * X0 + M[4] * X1 + M[5] * X2;
        // Y2 = M[2] * X0 + M[5] * X1 + M[8] * X2;
        const vec4f tt = loadu4f(X+3*j);
# if ( BLD == 4 )
        vec4f s0 = mul4f(streamload4f(D  ), permute4f(tt, 0x00));
        vec4f s1 = mul4f(streamload4f(D+4), permute4f(tt, 0x55));
        vec4f s2 = mul4f(streamload4f(D+8), permute4f(tt, 0xAA));
# else
        vec4f s0 = mul4f(load3fZ(D      ), permute4f(tt, 0x00));
        vec4f s1 = mul4f(load3fZ(D+BLD  ), permute4f(tt, 0x55));
        vec4f s2 = mul4f(load3fZ(D+BLD*2), permute4f(tt, 0xAA));
# endif
        storeu4f(Y+3*j, add4f(s2, add4f(s0, s1)));
    }
}
#endif


void SparMatSymBlkDiag::vecMul(const real* X, real* Y) const
{
#if ( BLOCK_SIZE == 3 ) && SMSBD_USES_AVX && REAL_IS_DOUBLE
    
    // process diagonal:
    vecMulDiagonal3D(X, Y);
    
    // process off-diagonal elements:
    const size_t stop = size_ / 3;
    for ( size_t j = colix_[0]; j < stop; j = colix_[j+1] )
        pilar_[j].vecMulAddTriangle3D_AVX(X, Y, j*3);

#elif ( BLOCK_SIZE == 3 ) && SMSBD_USES_SSE && !REAL_IS_DOUBLE
    
    // process diagonal:
    vecMulDiagonal3D(X, Y);
    
    // process off-diagonal elements:
    const size_t stop = size_ / 3;
    for ( size_t j = colix_[0]; j < stop; j = colix_[j+1] )
        pilar_[j].vecMulAddTriangle3D_SSE(X, Y, j*3);
    
#else
    zero_real(size_, Y);
    vecMulAdd(X, Y, 0, size_);
#endif
}
