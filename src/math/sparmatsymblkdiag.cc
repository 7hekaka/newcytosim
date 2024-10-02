// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include <cmath>
#include <sstream>

#include "sparmatsymblkdiag.h"
#include "assert_macro.h"
#include "vector2.h"
#include "vector3.h"
#include "simd.h"

// Flags to enable SIMD implementation
#if defined(__AVX__)
#  include "simd_float.h"
#  define SMSBD_USES_AVX 1
#  define SMSBD_USES_SSE 1
#elif USE_SIMD
#  include "simd_float.h"
#  define SMSBD_USES_AVX 0
#  define SMSBD_USES_SSE 1
#else
#  define SMSBD_USES_AVX 0
#  define SMSBD_USES_SSE 0
#endif


SparMatSymBlkDiag::SparMatSymBlkDiag()
: rsize_(0), alloc_(0), pilar_(nullptr)
{
    colix_ = new unsigned[2]();
}


void SparMatSymBlkDiag::allocate(index_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be adjusted to tune performance: by increasing 'chunk',
          more memory will be used, but reallocation will be less frequent
        */
        constexpr index_t chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        
        //fprintf(stderr, "SMSBD allocates %u\n", alc);
        Pilar * ptr = new Pilar[alc];
        
        if ( pilar_ )
        {
            for (index_t n = 0; n < alloc_; ++n )
                ptr[n] = pilar_[n];
            delete[] pilar_;
        }
        
        pilar_ = ptr;
        alloc_ = alc;
        
        delete[] colix_;
        colix_ = new unsigned[alc+2];
        for ( unsigned n = 0; n <= alc; ++n )
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
#pragma mark - Pilar

SparMatSymBlkDiag::Pilar::Pilar()
{
    noff_ = 0;
    allo_ = 0;
    dia_.reset();
    inx_ = nullptr;
    blk_ = nullptr;
}


/*
This may require some smart allocation scheme.
*/
void SparMatSymBlkDiag::Pilar::allocate(index_t alc)
{
    if ( alc > allo_ )
    {
        /*
         'chunk' can be adjusted to tune performance: by increasing 'chunk',
          more memory will be used, but reallocation will be less frequent
        */
        constexpr index_t chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk - 1 );
        
        //if ( inx_ ) fprintf(stderr, "SMSBD reallocates column %u: %lu\n", inx_[0], alc);
        //else fprintf(stderr, "SMSBD allocates column: %lu\n", alc);

        // use aligned memory, with some extra for SIMD burr:
        void * ptr = new_real(alc*sizeof(Block)/sizeof(real)+chunk);
        Block * blk_new  = new(ptr) Block[alc];

        if ( posix_memalign(&ptr, 32, alc*sizeof(unsigned)) )
            throw std::bad_alloc();
        auto * inx_new = (unsigned*)ptr;

        if ( inx_ )
        {
            for ( index_t n = 0; n < noff_; ++n )
                inx_new[n] = inx_[n];
            free(inx_);
        }

        if ( blk_ )
        {
            for ( index_t n = 0; n < noff_; ++n )
                blk_new[n] = blk_[n];
            free_real(blk_);
        }
        inx_  = inx_new;
        blk_  = blk_new;
        allo_ = alc;
        
        //std::clog << "Pilar " << this << "  " << alc << ": ";
        //std::clog << " alignment " << ((uintptr_t)elem_ & 63) << "\n";
    }
}


void SparMatSymBlkDiag::Pilar::deallocate()
{
    //if ( inx_ ) fprintf(stderr, "SMSBD deallocates column %u of size %lu\n", inx_[0], allo_);
    free(inx_);
    free_real(blk_);
    inx_ = nullptr;
    blk_ = nullptr;
    allo_ = 0;
    noff_ = 0;
}


void SparMatSymBlkDiag::Pilar::operator = (SparMatSymBlkDiag::Pilar & col)
{
    //if ( inx_ ) fprintf(stderr, "SMSBD transfers column %u\n", inx_[0]);
    free(inx_);
    free_real(blk_);

    noff_ = col.noff_;
    allo_ = col.allo_;
    dia_ = col.dia_;
    inx_ = col.inx_;
    blk_ = col.blk_;
    
    col.noff_ = 0;
    col.allo_ = 0;
    col.dia_.reset();
    col.inx_ = nullptr;
    col.blk_ = nullptr;
}


/* This is a silly search that could be optimized */
SparMatSymBlkDiag::Block* SparMatSymBlkDiag::Pilar::find_block(index_t ii) const
{
    for ( index_t n = 0; n < noff_; ++n )
        if ( inx_[n] == ii )
            return blk_ + n;
    return nullptr;
}

/**
 This allocates to be able to hold the matrix element if necessary
 */
SparMatSymBlkDiag::Block& SparMatSymBlkDiag::Pilar::block(index_t ii)
{
    SparMatSymBlkDiag::Block * B = find_block(ii);
    if ( !B )
    {
        allocate(noff_+1);
        B = blk_ + noff_;
        inx_[noff_] = ii;
        B->reset();
        ++noff_;
    }
    //printColumn(jj);
    return *B;
}


void SparMatSymBlkDiag::Pilar::reset()
{
    noff_ = 0;
    dia_.reset();
}


real& SparMatSymBlkDiag::element(index_t iii, index_t jjj)
{
    // branchless code to address lower triangle
    index_t ii = std::max(iii, jjj);
    index_t jj = std::min(iii, jjj);
#if ( SD_BLOCK_SIZE == 1 )
    if ( ii == jj )
        return pilar_[jj].dia_.value();
    return pilar_[jj].block(ii).value();
#else
    index_t i = ii / SD_BLOCK_SIZE;
    index_t j = jj / SD_BLOCK_SIZE;
    index_t ir = ii % SD_BLOCK_SIZE;
    index_t jr = jj % SD_BLOCK_SIZE;
    if ( i == j )
        return pilar_[j].dia_(ir, jr);
    return pilar_[j].block(i)(ir, jr);
#endif
}


real* SparMatSymBlkDiag::address(index_t iii, index_t jjj) const
{
    // branchless code to address lower triangle
    index_t ii = std::max(iii, jjj);
    index_t jj = std::min(iii, jjj);
#if ( SD_BLOCK_SIZE == 1 )
    if ( ii == jj )
        return pilar_[jj].dia_.data();
    return pilar_[jj].block(ii).data();
#else
    index_t i = ii / SD_BLOCK_SIZE;
    index_t j = jj / SD_BLOCK_SIZE;
    index_t ir = ii % SD_BLOCK_SIZE;
    index_t jr = jj % SD_BLOCK_SIZE;
    if ( i == j )
        return pilar_[j].dia_.addr(ir, jr);
    Block * B = pilar_[j].find_block(i);
    if ( B )
        return B->addr(ir, jr);
    return nullptr;
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Stuff

void SparMatSymBlkDiag::reset()
{
    for ( index_t n = 0; n < alloc_; ++n )
        pilar_[n].reset();
}


bool SparMatSymBlkDiag::notZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Pilar & col = pilar_[jj];
        if ( col.dia_ != 0.0 )
            return true;
        for ( index_t n = 0 ; n < col.noff_ ; ++n )
            if ( col[n] != 0.0 )
                return true;
    }
    //if here, the matrix is empty
    return false;
}


void SparMatSymBlkDiag::scale(const real alpha)
{
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Pilar & col = pilar_[jj];
        col.dia_.scale(alpha);
        for ( index_t n = 0 ; n < col.noff_ ; ++n )
            col[n].scale(alpha);
    }
}


void SparMatSymBlkDiag::addDiagonalBlock(real* mat, index_t ldd, const index_t start, const index_t cnt, index_t mul) const
{
    assert_true( mul == SD_BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );
    
    for ( index_t jj = 0; jj < cnt; ++jj )
    {
        Pilar & col = pilar_[jj+start];
        real * dst = mat + ( jj + ldd * jj ) * SD_BLOCK_SIZE;
        col.dia_.addto_symm(dst, ldd);
        for ( index_t n = 0; n < col.noff_; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col.inx_[n] > jj + start );
            auto ij = col.inx_[n] - ( jj + start );
            if ( jj+ij < cnt )
            {
                //fprintf(stderr, "SMSBD %4lu %4lu\n", ii, jj); col[n].print(stderr);
                col[n].addto(dst+ij*SD_BLOCK_SIZE, ldd);
                col[n].addto_trans(dst+ij*ldd*SD_BLOCK_SIZE, ldd);
            }
        }
    }
}


void SparMatSymBlkDiag::addLowerBand(real alpha, real* mat, index_t ldd, const index_t start, const index_t cnt,
                                     const index_t mul, const index_t rank) const
{
    assert_true( mul == SD_BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );

    for ( index_t jj = 0; jj < cnt; ++jj )
    {
        Pilar & col = pilar_[jj+start];
        index_t sup = std::min(cnt-jj, rank+1);
        real * dst = mat + ( jj * ldd + jj ) * SD_BLOCK_SIZE;
        col.dia_.addto_lower(dst, ldd, alpha);
        for ( index_t n = 0; n < col.noff_; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col.inx_[n] > jj + start );
            auto ij = col.inx_[n] - ( jj + start );
            if ( ij < sup )
            {
                //fprintf(stderr, "SMSBD %4lu %4lu\n", ii, jj); col[n].print(stderr);
                col[n].addto(dst+ij*SD_BLOCK_SIZE, ldd, alpha);
                //col[n].addto_trans(dst+ij*ldd*SD_BLOCK_SIZE, ldd, alpha);
            } else break; // assuming the column is sorted
        }
    }
}


void SparMatSymBlkDiag::addDiagonalTrace(real alpha, real* mat, index_t ldd,
                                         const index_t start, const index_t cnt,
                                         const index_t mul, const index_t rank, const bool sym) const
{
    assert_true( mul == SD_BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );

    for ( index_t jj = 0; jj < cnt; ++jj )
    {
        Pilar & col = pilar_[jj+start];
        index_t sup = std::min(cnt-jj, rank+1);
        real * dst = mat + ( jj * ldd + jj );
        dst[0] += alpha * col.dia_.trace();   // diagonal term
        for ( index_t n = 0; n < col.noff_; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col.inx_[n] > jj + start );
            auto ij = col.inx_[n] - ( jj + start );
            if ( ij < sup )
            {
                real a = alpha * col[n].trace();
                //fprintf(stderr, "SMSBD %4lu %4lu : %.4f\n", i, j, a);
                dst[ij] += a;
                if ( sym ) dst[ldd*ij] += a;
            } else break; // assuming the column is sorted
        }
    }
}



int SparMatSymBlkDiag::bad() const
{
    for ( index_t j = 0; j < rsize_; ++j )
    {
        Pilar & col = pilar_[j];
        for ( index_t n = 0 ; n < col.noff_ ; ++n )
        {
            if ( col.inx_[n] >= rsize_ ) return 2;
            if ( col.inx_[n] <= j ) return 3;
        }
    }
    return 0;
}


/** all elements are counted, even if zero */
size_t SparMatSymBlkDiag::nbElements(index_t start, index_t stop, size_t& alc) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    alc = 0;
    //index_t cnt = stop - start; // counting diagonal elements
    index_t cnt = 0; // not counting diagonal elements
    for ( index_t i = start; i < stop; ++i )
    {
        cnt += pilar_[i].noff_;
        alc += pilar_[i].allo_;
    }
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark - I/O


std::string SparMatSymBlkDiag::what() const
{
    size_t alc = 0;
    size_t cnt = nbElements(0, rsize_, alc);
    std::ostringstream msg;
#if SMSBD_USES_AVX && REAL_IS_DOUBLE
    msg << "SMSBDx ";
#elif SMSBD_USES_SSE && REAL_IS_DOUBLE
    msg << "SMSBDe ";
#elif SMSBD_USES_SSE && !REAL_IS_DOUBLE
    msg << "SMSBDf ";
#else
    msg << "SMSBD ";
#endif
    msg << Block::what() << "*" << cnt << " (" << alc << ")";
    return msg.str();
}


static void printSparseBlock(std::ostream& os, real inf, SparMatSymBlkDiag::Block const& B, index_t ii, index_t jj)
{
    index_t d = ( ii == jj );
    for ( index_t x = 0; x < B.dimension(); ++x )
    for ( index_t y = d*x; y < B.dimension(); ++y )
    {
        real v = B(y, x);
        if ( abs_real(v) > inf )
            os << std::setw(6) << ii+y << " " << std::setw(6) << jj+x << " " << std::setw(16) << v << "\n";
    }
}


void SparMatSymBlkDiag::printSparse(std::ostream& os, real inf, index_t start, index_t stop) const
{
    os << "% SparMatSymBlkDiag size " << rsize_ << ":\n";
    stop = std::min(stop, rsize_);
    std::streamsize p = os.precision(8);
    if ( ! pilar_ )
        return;
    for ( index_t j = start; j < stop; ++j )
    {
        Pilar & col = pilar_[j];
        os << "% column " << j << "\n";
        index_t jj = j *SD_BLOCK_SIZE;
        printSparseBlock(os, inf, col.dia_, jj, jj);
        for ( index_t n = 0 ; n < col.noff_ ; ++n )
            printSparseBlock(os, inf, col.blk_[n], col.inx_[n]*SD_BLOCK_SIZE, jj);
    }
    os.precision(p);
}


void SparMatSymBlkDiag::printSummary(std::ostream& os, index_t start, index_t stop)
{
    stop = std::min(stop, rsize_);
    os << "\nSMSBD size " << rsize_ << ":";
    for ( index_t j = start; j < stop; ++j )
    {
        if ( pilar_[j].noff_ > 0 )
        {
            os << "\n   " << j << "   " << pilar_[j].noff_;
            os << " ---> " << colix_[j];
        }
    }
    std::endl(os);
}


void SparMatSymBlkDiag::Pilar::printBlocks(std::ostream& os) const
{
    os << " " << dia_;
    for ( index_t n = 0; n < noff_; ++n )
        os << " " << inx_[n] << " " << blk_[n];
}


void SparMatSymBlkDiag::printBlocks(std::ostream& os) const
{
    for ( index_t j = 0; j < rsize_; ++j )
    {
        os << "\nSMSBD  col " << j;
        pilar_[j].printBlocks(os);
    }
    std::endl(os);
}

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication


/// A block element of the sparse matrix suitable for qsort()
class alignas(4*sizeof(real)) SparMatSymBlkDiag::Element
{
public:
    /// block element
    real blk[SD_BLOCK_SIZE*SD_BLOCK_SIZE];

    /// index
    index_t inx;
};


/// qsort function comparing line indices
static int compareSMSBDElement(const void * A, const void * B)
{
    index_t a = static_cast<SparMatSymBlkDiag::Element const*>(A)->inx;
    index_t b = static_cast<SparMatSymBlkDiag::Element const*>(B)->inx;
    
    return ( a > b ) - ( b > a );
}

/**
 This copies the data to the provided temporary array
 */
void SparMatSymBlkDiag::Pilar::sortElements(Element tmp[], index_t tmp_size)
{
    assert_true( noff_ <= tmp_size );
    for ( index_t i = 0; i < noff_; ++i )
    {
        blk_[i].store(tmp[i].blk);
        tmp[i].inx = inx_[i];
    }
    
    //std::clog << "sizeof(SparMatSymBlkDiag::Element) " << sizeof(Element) << "\n";
    qsort(tmp, noff_, sizeof(Element), &compareSMSBDElement);
    
    for ( index_t i = 0; i < noff_; ++i )
    {
         blk_[i].load(tmp[i].blk);
         inx_[i] = tmp[i].inx;
    }
}


index_t SparMatSymBlkDiag::newElements(Element*& ptr, index_t cnt)
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


void SparMatSymBlkDiag::sortElements()
{
    index_t tmp_size = 0;
    Element * tmp = nullptr;
    
    for ( index_t j = 0; j < rsize_; ++j )
    {
        Pilar & col = pilar_[j];
        //std::clog << "SMSBD column " << j << " has 1+" << col.noff_ << " elements\n";

        if ( col.noff_ > 1 )
        {
            // order the elements within the column:
            if ( tmp_size < col.noff_ )
                tmp_size = newElements(tmp, col.noff_);
            col.sortElements(tmp, tmp_size);
        }

#ifndef NDEBUG
        for ( index_t n = 0 ; n < col.noff_ ; ++n )
        {
            const auto i = col.inx_[n];
            assert_true( i < rsize_ );
            assert_true( i != j );
        }
#endif
    }
    // release memory:
    free(tmp);
}


bool SparMatSymBlkDiag::prepareForMultiply(int)
{
    index_t end = rsize_;
    if ( rsize_ > 0 )
    {
        index_t inx = end;
        index_t nxt = end;
        while ( inx-- > 0 )
        {
            if ( pilar_[inx].noff_ > 0 )
                nxt = inx;
            else if ( pilar_[inx].allo_ > 32 )
                pilar_[inx].deallocate();
            colix_[inx] = nxt;
        }
    }
    colix_[end] = end;

    bool res = false;
    for ( index_t j = 0; j < end; ++j )
    {
        pilar_[j].dia_.copy_lower();
        res |= ( pilar_[j].notEmpty() );
    }
    
    sortElements();
    
    //printSummary(std::cout, 0, rsize_);
#if 0
    std::clog.precision(2);
    for ( index_t j = 0; j < end; ++j )
    {
        Pilar const& pil = pilar_[j];
        for ( index_t i = 0; i < std::min(8UL, pil.noff_); ++i )
            std::clog << std::setw(4) << i << " " << std::setw(9) << pil[i] << "\n";
    }
#endif
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - Column Vector Multiplication, reference implementations


#if ( SD_BLOCK_SIZE == 1 )
void SparMatSymBlkDiag::Pilar::vecMulAdd1D(const real* X, real* Y, index_t jj) const
{
    const real xx = X[jj];
    real D = dia_.value();
    real yy = Y[jj] + D * xx;
    for ( index_t n = 0; n < noff_; ++n )
    {
        const real M = blk_[n].value();
        const auto ii = inx_[n];
        Y[ii] += M * xx;
        yy += M * X[ii];
    }
    Y[jj] = yy;
}

void SparMatSymBlkDiag::Pilar::vecMulAddTriangle1D(const real* X, real* Y, index_t jj) const
{
    assert_true(noff_ > 0);
    const real xx = X[jj];
    real yy = Y[jj];
    for ( index_t n = 0; n < noff_; ++n )
    {
        const real M = blk_[n].value();
        const auto ii = inx_[n];
        Y[ii] += M * xx;
        yy += M * X[ii];
    }
    Y[jj] = yy;
}

void SparMatSymBlkDiag::vecMulDiagonal1D(const real* X, real* Y) const
{
    Pilar const* pil = pilar_;
    for ( index_t j = 0; j < rsize_; ++j, ++pil )
    {
        const real M = pil->dia_.value();
        Y[j] = X[j] * M;
    }
}
#endif


#if ( SD_BLOCK_SIZE == 2 )
void SparMatSymBlkDiag::Pilar::vecMulAdd2D(const real* X, real* Y, index_t jj) const
{
    const Vector2 xx(X+jj);
    assert_small(dia_.asymmetry());
    Vector2 yy = dia_.vecmul(xx);
    for ( index_t n = 0; n < noff_; ++n )
    {
        const auto ii = 2 * inx_[n];
        Block const& M = blk_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.add_to(Y+jj);
}

void SparMatSymBlkDiag::Pilar::vecMulAddTriangle2D(const real* X, real* Y, index_t jj) const
{
    assert_true(noff_ > 0);
    Vector2 yy(Y+jj);
    const Vector2 xx(X+jj);
    for ( index_t n = 0; n < noff_; ++n )
    {
        Block const& M = blk_[n];
        const auto ii = 2 * inx_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.store(Y+jj);
}

void SparMatSymBlkDiag::vecMulDiagonal2D(const real* X, real* Y) const
{
    Pilar const* pil = pilar_;
    for ( index_t j = 0; j < rsize_; ++j, ++pil )
    {
        assert_small(pil->dia_.asymmetry());
        Vector2 yy = pil->dia_.vecmul(X+2*j);
        yy.store(Y+2*j);
    }
}
#endif

#if ( SD_BLOCK_SIZE == 3 )
void SparMatSymBlkDiag::Pilar::vecMulAdd3D(const real* X, real* Y, index_t jj) const
{
    const Vector3 xx(X+jj);
    assert_small(dia_.asymmetry());
    Vector3 yy = dia_.vecmul(xx);
    for ( index_t n = 0; n < noff_; ++n )
    {
        Block const& M = blk_[n];
        const auto ii = 3 * inx_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.add_to(Y+jj);
}

void SparMatSymBlkDiag::Pilar::vecMulAddTriangle3D(const real* X, real* Y, index_t jj) const
{
    assert_true(noff_ > 0);
    const Vector3 xx(X+jj);
    Vector3 yy(Y+jj);
    for ( index_t n = 0; n < noff_; ++n )
    {
        Block const& M = blk_[n];
        const auto ii = 3 * inx_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.store(Y+jj);
}

void SparMatSymBlkDiag::vecMulDiagonal3D(const real* X, real* Y) const
{
    Pilar const* pil = pilar_;
    for ( index_t j = 0; j < rsize_; ++j, ++pil )
    {
        assert_small(pil->dia_.asymmetry());
        Vector3 yy = pil->dia_.vecmul(X+3*j);
        yy.store(Y+3*j);
    }
}
#endif


#if ( SD_BLOCK_SIZE == 4 )
void SparMatSymBlkDiag::Pilar::vecMulAdd4D(const real* X, real* Y, index_t jj) const
{
    const vec4 xxxx = load4(X+jj);
    assert_small(dia_.asymmetry());
    vec4 yyyy = dia_.vecmul4(xxxx);
    for ( index_t n = 0; n < noff_; ++n )
    {
        Block const& M = blk_[n];
        const auto ii = 4 * inx_[n];
        store4(Y+ii, add4(load4(Y+ii), M.vecmul4(xxxx)));
        yyyy += M.trans_vecmul(X+ii);
    }
    store4(Y+jj, add4(yyyy, load4(Y+jj)));
}

void SparMatSymBlkDiag::Pilar::vecMulAddTriangle4D(const real* X, real* Y, index_t jj) const
{
    assert_true(noff_ > 0);
    Vector4 yy(Y+jj);
    const Vector4 xx(X+jj);
    for ( index_t n = 0; n < noff_; ++n )
    {
        Block const& M = blk_[n];
        const auto ii = 4 * inx_[n];
        M.vecmul(xx).add_to(Y+ii);
        yy += M.trans_vecmul(X+ii);
    }
    yy.store(Y+jj);
}

void SparMatSymBlkDiag::vecMulDiagonal4D(const real* X, real* Y) const
{
    Pilar const* pil = pilar_;
    for ( index_t j = 0; j < rsize_; j += 4, ++pil )
    {
        assert_small(pil->dia_.asymmetry());
        Vector4 yy = pil->dia_.vecmul(X+4*j);
        yy.store(Y+4*j);
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - 2D Double Precision Optimized Vector Multiplication

#if ( SD_BLOCK_SIZE == 2 ) && SMSBD_USES_SSE && REAL_IS_DOUBLE
void SparMatSymBlkDiag::Pilar::vecMulAdd2D_SSE(const double* X, double* Y, index_t jj) const
{
    vec2 s1 = load2(X+jj);
    //const real X0 = X[jj  ];
    //const real X1 = X[jj+1];
    vec2 x0 = unpacklo2(s1, s1);
    vec2 x1 = unpackhi2(s1, s1);

    assert_small(dia_.asymmetry());
    // load 2x2 matrix element into 2 vectors:
    double const* D = dia_.data();
    //assume the block is already symmetrized:
    //real Y0 = Y[jj  ] + D[0] * X0 + D[1] * X1;
    //real Y1 = Y[jj+1] + D[1] * X0 + D[3] * X1;
    vec2 s0 = mul2(load2(D), s1);
    s1 = mul2(load2(D+2), s1);
    
    Block const* blk = blk_;
    auto const* inx = inx_;

    Block const* end = blk_ + ( noff_ & ~1 );
#if ( 1 )
    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( ; blk < end; inx += 2, blk += 2 )
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
        s0 = fmadd2(m01, xx, s0);
        s1 = fmadd2(m23, xx, s1);
        vec2 t = fmadd2(m01, x0, load2(Y+i0));
        vec2 u = fmadd2(p01, x0, load2(Y+i1));
        // multiply with the transposed block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1];
        //Y1 += M[2] * X[ii] + M[3] * X[ii+1];
        s0 = fmadd2(p01, tt, s0);
        s1 = fmadd2(p23, tt, s1);
        store2(Y+i0, fmadd2(m23, x1, t));
        store2(Y+i1, fmadd2(p23, x1, u));
    }
#endif
    end = blk_ + noff_;
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
        s0 = fmadd2(m0, xx, s0); // s0 += m0 * xx;
        s1 = fmadd2(m2, xx, s1); // s1 += m2 * xx;
    }
    //Y[jj  ] = Y0;
    //Y[jj+1] = Y1;
    s0 = add2(unpacklo2(s0, s1), unpackhi2(s0, s1));
    store2(Y+jj, add2(s0, load2(Y+jj)));
}
#endif


#if ( SD_BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAdd2D_AVX(const double* X, double* Y, index_t jj) const
{
    // xy = { X0 X1 X0 X1 }
    vec4 xy = broadcast2(X+jj);
    //multiply with full block, assuming it is symmetric:
    // Y0 = M[0] * X0 + M[1] * X1;
    // Y1 = M[1] * X0 + M[3] * X1;
    
    // yyyy = { Y0 Y0 Y1 Y1 }
    // load 2x2 matrix element into 2 vectors:
    vec4 ss = mul4(dia_.data0(), xy);

    //const real X0 = X[jj  ];
    //const real X1 = X[jj+1];
    // xxyy = { X0 X0 X1 X1 }
    const vec4 xxyy = duplohi4(xy);

    // while x0 and x1 are constant, there is a dependency in the loop for 'yy'.
    for ( index_t n = 0; n < noff_; ++n )
    {
        const auto ii = 2 * inx_[n];
        vec4 mat = blk_[n].data0(); // load 2x2 matrix
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


#if ( SD_BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
static inline void multiply2D(double const* X, double* Y, index_t ii, vec4 const& mat, vec4 const& xxxx, vec4& ss)
{
    vec4 xx = broadcast2(X+ii);
    vec4 u = fmadd4(mat, xxxx, load2Z(Y+ii));
    store2(Y+ii, add2(getlo(u), gethi(u)));
    ss = fmadd4(mat, xx, ss);
}
#endif


#if ( SD_BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAdd2D_AVXU(const double* X, double* Y, index_t jj) const
{
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(dia_.data0(), xyxy);
    const vec4 xxyy = duplohi4(xyxy);
    vec4 s1 = setzero4();
    
    Block const* blk = blk_;
    Block const* end = blk_ + ( noff_ & ~1 );
    auto const* inx = inx_;
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
    end = blk_ + noff_;
    //#pragma nounroll
    if ( blk < end ) // for ( ; blk < end; ++blk, ++inx )
        multiply2D(X, Y, 2*inx[0], blk[0].data0(), xxyy, ss);
    /* finally horizontally sum ss = { SX SX SY SY } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


#if ( SD_BLOCK_SIZE == 2 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAdd2D_AVXUU(const double* X, double* Y, index_t jj) const
{
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(dia_.data0(), xyxy);
    const vec4 xxyy = duplohi4(xyxy);
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();

    Block  const* blk = blk_;
    Block const* end = blk_ + ( noff_ & ~3 );
    auto const* inx = inx_;

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
        vec4 mat0 = blk[0].data0();
        vec4 mat1 = blk[1].data0();
        vec4 mat2 = blk[2].data0();
        vec4 mat3 = blk[3].data0();
        const auto i0 = 2 * inx[0];
        const auto i1 = 2 * inx[1];
        const auto i2 = 2 * inx[2];
        const auto i3 = 2 * inx[3];
        assert_true( i0 < i1 );
        assert_true( i1 < i2 );
        assert_true( i2 < i3 );
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
    end = blk_ + noff_;
    //#pragma nounroll
    if ( blk < end ) // for ( ; blk < end; ++blk, ++inx )
        multiply2D(X, Y, 2*inx[0], blk[0].data0(), xxyy, ss);
    /* finally sum ss = { S0 S0 S1 S1 } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


//------------------------------------------------------------------------------
#pragma mark - 3D Single Precision Optimized Vector Multiplication

#if ( SD_BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && SMSBD_USES_SSE
void SparMatSymBlkDiag::Pilar::vecMulAdd3D_SSE(const float* X, float* Y, index_t jj) const
{
    //multiply with the diagonal block, assuming it is symmetric:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    const vec4f xxx = loadu4f(X+jj);
    vec4f s0 = mul4f(dia_.data0(), xxx);
    vec4f s1 = mul4f(dia_.data1(), xxx);
    vec4f s2 = mul4f(dia_.data2(), xxx);
    const vec4f x0 = clear4th(broadcastXf(xxx));
    const vec4f x1 = clear4th(broadcastYf(xxx));
    const vec4f x2 = clear4th(broadcastZf(xxx));

    Block const* blk = blk_;
    auto const* inx = inx_;

    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( index_t n = 0; n < noff_; ++n, ++blk, ++inx )
    {
        const auto ii = 3 * inx[0];
        const vec4f M012 = blk->data0();
        const vec4f M345 = blk->data1();
        const vec4f M678 = blk->data2();
        // multiply with the full block:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec4f z = fmadd4f(M012, x0, loadu4f(Y+ii));
        z = fmadd4f(M345, x1, z);
        z = fmadd4f(M678, x2, z);
        assert_true(z[3]==Y[ii+3]);
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
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4f s3 = setzero4f();
    s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
    s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
    s0 = add4f(catshift2f(s0, s2), blend22f(s0, s2));
    assert_true(s0[3]==0);
    storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
}
#endif


#if ( SD_BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && SMSBD_USES_SSE
void SparMatSymBlkDiag::Pilar::vecMulAdd3D_SSEU(const float* X, float* Y, index_t jj) const
{
    assert_small(dia_.asymmetry());
    //std::cout << dia_.to_string(7,1); printf(" MSSB %lu : %lu\n", jj, rsize_);
    //multiply with the diagonal block, assuming it has been symmetrized:
    // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
    // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
    // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    const vec4f xxx = loadu4f(X+jj);
    const vec4f x0 = clear4th(broadcastXf(xxx));
    const vec4f x1 = clear4th(broadcastYf(xxx));
    const vec4f x2 = clear4th(broadcastZf(xxx));

    Block const* blk = blk_;
    auto const* inx = inx_;

    if ( noff_ > 0 )
    {
        Block const* end = blk_ + ( noff_ & ~1 );
        // load 3x3 matrix diagonal element into 3 vectors:
        vec4f s0 = mul4f(dia_.data0(), xxx);
        vec4f s1 = mul4f(dia_.data1(), xxx);
        vec4f s2 = mul4f(dia_.data2(), xxx);
        {
            // process 2 by 2
            #pragma nounroll
            for ( ; blk < end; blk += 2, inx += 2 )
            {
                Block const& M = blk[0];
                Block const& P = blk[1];
                const auto ii = 3 * inx[0];
                const auto kk = 3 * inx[1];
                const vec4f M012 = M.data0();
                const vec4f M345 = M.data1();
                const vec4f M678 = M.data2();
                const vec4f P012 = P.data0();
                const vec4f P345 = P.data1();
                const vec4f P678 = P.data2();
                assert_true( ii < kk );
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
                assert_true(z[3]==Y[ii+3]);
                assert_true(t[3]==Y[kk+3]);
                storeu4f(Y+ii, z);
                storeu4f(Y+kk, t);
            }
        }
        // process remaining blocks
        end = blk_ + noff_;
#pragma nounroll
        for ( ; blk < end; ++inx, ++blk )
        {
            const auto ii = 3 * inx[0];
            const vec4f M012 = blk->data0();
            const vec4f M345 = blk->data1();
            const vec4f M678 = blk->data2();
            // multiply with the full block:
            //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
            //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
            //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
            vec4f z = fmadd4f(M012, x0, loadu4f(Y+ii));
            z = fmadd4f(M345, x1, z);
            z = fmadd4f(M678, x2, z);
            assert_true(z[3]==Y[ii+3]);
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
        s0 = clear4th(s0); // s0[3] is garbage
        s1 = clear4th(s1); // s1[3] is garbage
        s2 = clear4th(s2); // s2[3] is garbage
        vec4f s3 = setzero4f();
        s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
        s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
        s0 = add4f(catshift2f(s0, s2), blend22f(s0, s2));
        assert_true(s0[3]==0);
        storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
    }
    else
    {
        vec4f s0 = mul4f(dia_.data0(), x0);
        vec4f s1 = mul4f(dia_.data1(), x1);
        vec4f s2 = mul4f(dia_.data2(), x2);
        assert_true(s0[3]==0);
        storeu4f(Y+jj, add4f(add4f(loadu4f(Y+jj), s0), add4f(s1, s2)));
    }
}
#endif

/**
 Only process off-diagonal terms!
 */
#if ( SD_BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && SMSBD_USES_SSE
void SparMatSymBlkDiag::Pilar::vecMulAddTriangle3D_SSE(const float* X, float* Y, index_t jj) const
{
    assert_true(noff_ > 0);
    const vec4f xxx = loadu4f(X+jj);
    const vec4f x0 = clear4th(broadcastXf(xxx));
    const vec4f x1 = clear4th(broadcastYf(xxx));
    const vec4f x2 = clear4th(broadcastZf(xxx));

    vec4f s0 = setzero4f();
    vec4f s1 = setzero4f();
    vec4f s2 = setzero4f();

    Block const* blk = blk_;
    Block const* end = blk_ + noff_;
    auto const* inx = inx_;
    
    //index_t n = 0;
    if ( 1 ) {
        // process 2 by 2, blk+1 <= end-1 is blk < end-1
        #pragma nounroll
        for (Block const* stop = end-1; blk < stop; blk += 2, inx += 2 )
        {
            Block const& M = blk[0];
            Block const& P = blk[1];
            const auto ii = 3 * inx[0];
            const vec4f M012 = M.data0();
            const vec4f M345 = M.data1();
            const vec4f M678 = M.data2();
            const auto kk = 3 * inx[1];
            const vec4f P012 = P.data0();
            const vec4f P345 = P.data1();
            const vec4f P678 = P.data2();
            assert_true( ii < kk );
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
            assert_true(z[3]==Y[ii+3]);
            assert_true(t[3]==Y[kk+3]);
            storeu4f(Y+ii, z);
            storeu4f(Y+kk, t);
        }
    }
    //#pragma nounroll
    // process last remaining block
    //assert_true( blk == blk_+noff_-1 );
    for ( ; blk < end; ++blk, ++inx ) //if ( blk < end )
    {
        const auto ii = 3 * inx[0];
        const vec4f L012 = blk->data0();
        const vec4f L345 = blk->data1();
        const vec4f L678 = blk->data2();
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
        assert_true(z[3]==Y[ii+3]);
        storeu4f(Y+ii, z);
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
    assert_true(s0[3]==0);
    storeu4f(Y+jj, add4f(loadu4f(Y+jj), s0));
}
#endif

//------------------------------------------------------------------------------
#pragma mark - 3D Double Precision Optimized Vector Multiplication

#if ( SD_BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_SSE
void SparMatSymBlkDiag::Pilar::vecMulAdd3D_SIMD(const double* X, double* Y, index_t jj) const
{
    vec2 zero = setzero2();
    /* vec4 s0, s1, s2 accumulate lines of the transposed-matrix multiplied by 'xyz' */
    vec2 s0, s1, s2;
    vec2 t0, t2;
    vec2 xx, yy, zz;
    {
        vec2 xy = loadu2(X+jj);
        zz = loaddup2(X+jj+2);
        //multiply with the diagonal block, assuming it is symmetric:
        // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
        // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
        // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
        auto const& D = dia_.data();
        yy = loadu2(Y+jj);
        s0 = fmadd2(loadu2(D  ), xy, unpacklo2(yy, zero));
        s1 = fmadd2(loadu2(D+3), xy, unpackhi2(yy, zero));
        vec2 mat6 = loadu2(D+6);
        s2 = fmadd2(mat6, xy, unpacklo2(load1(Y+jj+2), zero));
        // prepare broadcasted vectors:
        xx = duplo2(xy);
        yy = duphi2(xy);
        // multiply last column into t0 = { X, Y } and t2 = { Z }:
        t0 = mul2(mat6, zz);
        t2 = mul1(load1(D+8), zz); // upper is garbage
    }
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    for ( index_t n = 0; n < noff_; ++n )
    {
        const auto ii = 3 * inx_[n];
        auto const& mat = blk_[n].data();
        // multiply with the full block in vectors { XY, Z0 }:
        //Y[ii  ] +=  M[0] * X0 + M[3] * X1 + M[6] * X2;
        //Y[ii+1] +=  M[1] * X0 + M[4] * X1 + M[7] * X2;
        //Y[ii+2] +=  M[2] * X0 + M[5] * X1 + M[8] * X2;
        vec2 mat0 = loadu2(mat);
        vec2 mat2 = loadu2(mat+2);
        vec2 XY = fmadd2(mat0, xx, loadu2(Y+ii));
        vec2 Z0 = fmadd1(mat2, xx, load1(Y+ii+2));  // upper not used

        vec2 xyi = loadu2(X+ii);
        vec2 mat4 = loadu2(mat+4);
        s0 = fmadd2(mat0, xyi, s0);
        // multiply with the transposed block in lines { s0, s1, s2 }:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];

        vec2 zzi = loaddup2(X+ii+2);
        vec2 mat3 = catshift(mat2, mat4);
        Z0 = fmadd1(unpackhi2(mat4, zero), yy, Z0);
        XY = fmadd2(mat3, yy, XY);
        s1 = fmadd2(mat3, xyi, s1);
        vec2 mat6 = loadu2(mat+6);
        t0 = fmadd2(blend11(mat2, mat4), zzi, t0);

        vec2 mat8 = load1(mat+8);   // upper is garbage
        XY = fmadd2(mat6, zz, XY);
        Z0 = fmadd1(mat8, zz, Z0);  // upper is garbage
        s2 = fmadd2(mat6, xyi, s2);
        t2 = fmadd1(mat8, zzi, t2); // upper is garbage
        storeu2(Y+ii, XY);
        store1(Y+ii+2, Z0);
    }
    s0 = add2(s0, unpacklo2(t0, zero));
    s1 = add2(s1, unpackhi2(t0, zero));
    s2 = add2(s2, unpacklo2(t2, zero));
    // finally sum s0 = { Y0 Y0 }, s1 = { Y1 Y1 }, s2 = { Y2 Y2 }
    s0 = add2(unpacklo2(s0, s1), unpackhi2(s0, s1));
    s2 = add2(unpacklo2(s2, s2), unpackhi2(s2, s2));
    storeu2(Y+jj, s0);
    store1(Y+jj+2, s2);
}
#endif

#if ( SD_BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAdd3D_AVX(const double* X, double* Y, index_t jj) const
{
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    vec4 s0, s1, s2;
    vec4 x0, x1, x2;
    {
        vec4 xxx = load3Z(X+jj);
        //multiply with the diagonal block, assuming it is symmetric:
        // Y0 = Y[jj  ] + M[0] * X0 + M[1] * X1 + M[2] * X2;
        // Y1 = Y[jj+1] + M[1] * X0 + M[4] * X1 + M[5] * X2;
        // Y2 = Y[jj+2] + M[2] * X0 + M[5] * X1 + M[8] * X2;
        s0 = mul4(dia_.data0(), xxx);
        s1 = mul4(dia_.data1(), xxx);
        s2 = mul4(dia_.data2(), xxx);
        // prepare broacasted vectors
#if 0
        const vec4 x0 = broadcast1(X+jj);
        const vec4 x1 = broadcast1(X+jj+1);
        const vec4 x2 = broadcast1(X+jj+2);
#else
        // alternative strategy:
        x0 = swap2f128(xxx);
        x1 = blend22(xxx, x0);
        x2 = blend22(x0, xxx);
        // zero out the 4th terms:
        x0 = clear4th(duplo4(x1));
        x1 = clear4th(duphi4(x1));
        x2 = clear4th(duplo4(x2));
#endif
    }
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( index_t n = 0; n < noff_; ++n )
    {
        const auto ii = 3 * inx_[n];
        const vec4 M012 = blk_[n].data0();
        const vec4 M345 = blk_[n].data1();
        const vec4 M678 = blk_[n].data2();
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
        vec4 xyz = load3Z(X+ii);  // xyz = { X0 X1 X2 0 }
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
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4 s3 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    s0 = add4(catshift2(s0, s2), blend22(s0, s2));
    assert_true(s0[3]==0);
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
#endif
}
#endif


#if ( SD_BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAdd3D_AVXU(const double* X, double* Y, index_t jj) const
{
    vec4 s0, s1, s2;
    vec4 x0, x1, x2;
    // load 3x3 matrix element into 3 vectors:
    {
        vec4 xxx = load3Z(X+jj);
        //multiply with the symmetric diagonal block:
        assert_small(dia_.asymmetry());
        s0 = mul4(dia_.data0(), xxx);
        s1 = mul4(dia_.data1(), xxx);
        s2 = mul4(dia_.data2(), xxx);
        // prepare broadcasted vectors:
        x0 = swap2f128(xxx);
        x1 = blend22(xxx, x0);
        x2 = blend22(x0, xxx);
        // zero out the 4th terms:
        x0 = clear4th(duplo4(x1));
        x1 = clear4th(duphi4(x1));
        x2 = clear4th(duplo4(x2));
    }
    if ( noff_ > 0 )
    {
        vec4 t0 = setzero4();
        vec4 t1 = setzero4();
        vec4 t2 = setzero4();
        // There is a dependency in the loop for 's0', 's1' and 's2'.
        Block const* blk = blk_;
        Block const* end = blk_ + ( noff_ & ~1 );
        auto const* inx = inx_;
        /*
         Unrolling will reduce the dependency chain, which may be limiting the
         throughput here. However the number of registers (16 for AVX CPU) limits
         the level of unrolling that can be done.
         */
        //process 2 by 2:
#pragma nounroll
        for ( ; blk < end; blk += 2, inx += 2 )
        {
            const auto ii = 3 * inx[0];
            const auto kk = 3 * inx[1];
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
        for ( end = blk_ + noff_; blk < end; ++blk, ++inx )
        {
            const auto ii = 3 * inx[0];
            //printf("--- %4i\n", ii);
            vec4 m0 = blk[0].data0();
            vec4 z = fmadd4(m0, x0, loadu4(Y+ii));
            vec4 xyz = loadu4(X+ii);
            s0 = fmadd4(m0, xyz, s0);
            
            vec4 m1 = blk[0].data1();
            z = fmadd4(m1, x1, z);
            s1 = fmadd4(m1, xyz, s1);
            
            vec4 m2 = blk[0].data2();
            z = fmadd4(m2, x2, z);
            s2 = fmadd4(m2, xyz, s2);
            storeu4(Y+ii, z);
        }
    }
    // finally sum s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    x0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, x0), unpackhi4(s2, x0));
    s0 = add4(catshift2(s0, s2), blend22(s0, s2));
    assert_true(s0[3]==0);
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
}
#endif



#if ( SD_BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAddTriangle3D_AVX(const double* X, double* Y, index_t jj) const
{
    vec4 x0, x1, x2;
    // load 3x3 matrix element into 3 vectors:
    {
        vec4 xxx = load3Z(X+jj);
        // prepare broadcasted vectors:
        x0 = swap2f128(xxx);
        x1 = blend22(xxx, x0);
        x2 = blend22(x0, xxx);
        // zero out the 4th terms:
        x0 = clear4th(duplo4(x1));
        x1 = clear4th(duphi4(x1));
        x2 = clear4th(duplo4(x2));
    }
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    if ( noff_ > 0 )
    {
        vec4 t0 = setzero4();
        vec4 t1 = setzero4();
        vec4 t2 = setzero4();
        // There is a dependency in the loop for 's0', 's1' and 's2'.
        Block const* blk = blk_;
        Block const* end = blk_  + ( noff_ & ~1 );
        auto const* inx = inx_;
        /*
         Unrolling will reduce the dependency chain, which may be limiting the
         throughput here. However the number of registers (16 for AVX CPU) limits
         the level of unrolling that can be done.
         */
        //process 2 by 2:
#pragma nounroll
        for ( ; blk < end; blk += 2, inx += 2 )
        {
            vec4 M0 = blk[0].data0();
            vec4 M1 = blk[0].data1();
            vec4 M2 = blk[0].data2();
            vec4 P0 = blk[1].data0();
            vec4 P1 = blk[1].data1();
            vec4 P2 = blk[1].data2();
            const auto ii = 3 * inx[0];
            const auto kk = 3 * inx[1];
            assert_true( ii < kk );
            //printf("--- %4i %4i\n", ii, kk);
            vec4 z = fmadd4(M0, x0, loadu4(Y+ii));
            vec4 t = fmadd4(P0, x0, loadu4(Y+kk));
            vec4 xyz = loadu4(X+ii);
            vec4 tuv = loadu4(X+kk);
            s0 = fmadd4(M0, xyz, s0);
            t0 = fmadd4(P0, tuv, t0);
            // multiply with the full block:
            z = fmadd4(M1, x1, z);
            t = fmadd4(P1, x1, t);
            s1 = fmadd4(M1, xyz, s1);
            t1 = fmadd4(P1, tuv, t1);
            
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
        
        // process remaining block:
        #pragma nounroll
        for ( end = blk_ + noff_; blk < end; ++blk, ++inx )
        {
            vec4 M0 = blk[0].data0();
            vec4 M1 = blk[0].data1();
            vec4 M2 = blk[0].data2();
            // extract index from the matrix data:
            const auto ii = 3 * inx[0];
            //printf("--- %4i\n", ii);
            
            vec4 z = fmadd4(M0, x0, loadu4(Y+ii));
            vec4 xyz = loadu4(X+ii);
            s0 = fmadd4(M0, xyz, s0);
            
            z = fmadd4(M1, x1, z);
            s1 = fmadd4(M1, xyz, s1);
            
            z = fmadd4(M2, x2, z);
            s2 = fmadd4(M2, xyz, s2);
            storeu4(Y+ii, z);
        }
    }
    // finally sum s0 = { Y0 Y0 Y0 0 }, s1 = { Y1 Y1 Y1 0 }, s2 = { Y2 Y2 Y2 0 }
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    x0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, x0), unpackhi4(s2, x0));
    s0 = add4(catshift2(s0, s2), blend22(s0, s2));
    assert_true(s0[3]==0);
    storeu4(Y+jj, add4(loadu4(Y+jj), s0));
}
#endif


#if ( SD_BLOCK_SIZE == 4 ) && REAL_IS_DOUBLE && SMSBD_USES_AVX
void SparMatSymBlkDiag::Pilar::vecMulAdd4D_AVX(const double* X, double* Y, index_t jj) const
{
    //multiply with the diagonal block, assuming it is symmetric:
    /* vec4 s0, s1, s2 add lines of the transposed-matrix multiplied by 'xyz' */
    vec4 s0, s1, s2, s3;
    vec4 x0, x1, x2, x3;
    {
        vec4 tt = load4(X+jj);
        s0 = mul4(dia_.data0(), tt);
        s1 = mul4(dia_.data1(), tt);
        s2 = mul4(dia_.data2(), tt);
        s3 = mul4(dia_.data3(), tt);
        vec4 x1 = duplo2f128(tt);
        vec4 x3 = duphi2f128(tt);
        vec4 x0 = duplo4(x1);
        x1 = duphi4(x1);
        vec4 x2 = duplo4(x3);
        x3 = duphi4(x3);
    }
    // sum non-diagonal elements:
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    #pragma nounroll
    for ( index_t n = 0; n < noff_; ++n )
    {
        const auto ii = 4 * inx_[n];
        const vec4 yy = load4(Y+ii);
        const vec4 xx = load4(X+ii);  // xx = { X0 X1 X2 X3 }
        const vec4 m0 = blk_[n].data0();
        vec4 z = fmadd4(m0, x0, yy);
        s0 = fmadd4(m0, xx, s0);
        
        const vec4 m1 = blk_[n].data1();
        z  = fmadd4(m1, x1, z);
        s1 = fmadd4(m1, xx, s1);

        const vec4 m2 = blk_[n].data2();
        z  = fmadd4(m2, x2, z);
        s2 = fmadd4(m2, xx, s2);

        const vec4 m3 = blk_[n].data3();
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
void SparMatSymBlkDiag::vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    for ( index_t j = start; j < stop; ++j )
    {
        //std::clog << "SparMatSymBlkDiag column " << j << "  " << rsize_ << " \n";
#if ( SD_BLOCK_SIZE == 1 )
        pilar_[j].vecMulAdd1D(X, Y, j);
#elif ( SD_BLOCK_SIZE == 2 )
        pilar_[j].vecMulAdd2D(X, Y, j*2);
#elif ( SD_BLOCK_SIZE == 3 )
        pilar_[j].vecMulAdd3D(X, Y, j*3);
#elif ( SD_BLOCK_SIZE == 4 )
        pilar_[j].vecMulAdd4D(X, Y, j*4);
#endif
    }
}


#if SMSBD_USES_AVX && REAL_IS_DOUBLE
#   define VECMULADD2D vecMulAdd2D_AVXU
#   define VECMULADD3D vecMulAdd3D_AVXU
#   define VECMULADD4D vecMulAdd4D_AVX
#elif SMSBD_USES_SSE && REAL_IS_DOUBLE
#   define VECMULADD2D vecMulAdd2D_SSE
#   define VECMULADD3D vecMulAdd3D_SIMD
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
void SparMatSymBlkDiag::vecMulAdd(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    for ( index_t j = start; j < stop; ++j )
    {
        //std::clog << "SparMatSymBlkDiag column " << jj << "  " << rsize_ << " \n";
#if ( SD_BLOCK_SIZE == 1 )
        pilar_[j].vecMulAdd1D(X, Y, j);
#elif ( SD_BLOCK_SIZE == 2 )
        pilar_[j].VECMULADD2D(X, Y, j*2);
#elif ( SD_BLOCK_SIZE == 3 )
        pilar_[j].VECMULADD3D(X, Y, j*3);
#elif ( SD_BLOCK_SIZE == 4 )
        pilar_[j].VECMULADD4D(X, Y, j*4);
#endif
    }
}


//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication

/*
This code is disabled here in favor of the next version which is unrolled
*/
#if ( SD_BLOCK_SIZE == 0 ) && REAL_IS_DOUBLE && defined(__AVX__)
void SparMatSymBlkDiag::vecMulDiagonal3D_AVX(const double* src, double* dst) const
{
    #pragma ivdep unroll (4)
    #pragma clang loop unroll_count(4)
    for ( index_t j = 0; j < rsize_; ++j )
    {
        Block const& D = pilar_[j].dia_;
        /**
         Since the diagonal block is symmetric, we only need to load 6 scalars
         instead of 9 here. Moreover, the elements of the vector and matrix can
         be handled together to limit the numbers of swaps.
         Particular in SSE code this should be faster than the code below. */
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
        src += 3;
        //multiply with the diagonal block, which is symmetric:
        // Y0 = M[0] * X0 + M[3] * X1 + M[6] * X2;
        // Y1 = M[1] * X0 + M[4] * X1 + M[7] * X2;
        // Y2 = M[2] * X0 + M[5] * X1 + M[8] * X2;
        x0 = mul4(D.data0(), x0);
        x0 = fmadd4(D.data1(), x1, x0);
        x0 = fmadd4(D.data2(), x2, x0);
        storeu4(dst, x0); // the 4-th is overwritten in next iteration
        dst += 3;
    }
    // clear garbage:
    dst[0] = 0;
}
#endif



#if ( SD_BLOCK_SIZE == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
void SparMatSymBlkDiag::vecMulDiagonal3D_AVX(const double* src, double* dst) const
{
    index_t j = 0;
    if ( rsize_ & 1 )
    {
        Block const& D = pilar_[0].dia_;
        vec4 x0 = broadcast1(src  );
        vec4 x1 = broadcast1(src+1);
        vec4 x2 = broadcast1(src+2);
        src += 3;
        //multiply with the diagonal block, which is symmetric:
        // Y0 = M[0] * X0 + M[3] * X1 + M[6] * X2;
        // Y1 = M[1] * X0 + M[4] * X1 + M[7] * X2;
        // Y2 = M[2] * X0 + M[5] * X1 + M[8] * X2;
        x0 = mul4(D.data0(), x0);
        x0 = fmadd4(D.data1(), x1, x0);
        x0 = fmadd4(D.data2(), x2, x0);
        storeu4(dst, x0); // the 4-th is overwritten in next iteration
        dst += 3;
        ++j;
    }

    #pragma ivdep
    #pragma clang loop unroll_count(2)
    while ( j < rsize_ )
    {
        Block const& D = pilar_[j++].dia_;
        Block const& N = pilar_[j++].dia_;
        //broadcast the source vectors:
        vec4 x0 = broadcast1(src);
        vec4 x1 = broadcast1(src+1);
        vec4 x2 = broadcast1(src+2);
        vec4 x3 = broadcast1(src+3);
        vec4 x4 = broadcast1(src+4);
        vec4 x5 = broadcast1(src+5);
        src += 6;
        //multiply with the matrix block:
        x0 = mul4(D.data0(), x0);
        x1 = mul4(D.data1(), x1);
        x2 = mul4(D.data2(), x2);
        x3 = mul4(N.data0(), x3);
        x4 = mul4(N.data1(), x4);
        x5 = mul4(N.data2(), x5);
        storeu4(dst,   add4(add4(x0, x1), x2)); // the 4-th is overwritten in next iteration
        storeu4(dst+3, add4(add4(x3, x4), x5)); // the 4-th is overwritten in next iteration
        dst += 6;
    }
    // clear garbage:
    dst[0] = 0;
}
#endif


#if ( SD_BLOCK_SIZE == 3 ) && !REAL_IS_DOUBLE && SMSBD_USES_SSE
void SparMatSymBlkDiag::vecMulDiagonal3D_SSE(const float* X, float* Y) const
{
    #pragma unroll (4)
    for ( index_t j = 0; j < rsize_; ++j )
    {
        Block const& D = pilar_[j].dia_;
        //multiply with the diagonal block, assuming it has been symmetrized:
        // Y0 = M[0] * X0 + M[1] * X1 + M[2] * X2;
        // Y1 = M[1] * X0 + M[4] * X1 + M[5] * X2;
        // Y2 = M[2] * X0 + M[5] * X1 + M[8] * X2;
        const vec4f xxx = loadu4f(X+3*j); // garbage 4-th term is not used
        vec4f s0 = mul4f(D.data0(), broadcastXf(xxx));
        vec4f s1 = mul4f(D.data1(), broadcastYf(xxx));
        vec4f s2 = mul4f(D.data2(), broadcastZf(xxx));
        // garbage 4-th term will be overwritten
        storeu4f(Y+3*j, add4f(s2, add4f(s0, s1)));
    }
    // clear garbage:
    Y[3*rsize_] = 0;
}
#endif


void SparMatSymBlkDiag::vecMul(const real* X, real* Y) const
{
#if ( SD_BLOCK_SIZE == 3 ) && SMSBD_USES_AVX && REAL_IS_DOUBLE
    
    // process diagonal:
    vecMulDiagonal3D_AVX(X, Y);
    
    // process off-diagonal elements:
    for ( index_t j = colix_[0]; j < rsize_; j = colix_[j+1] )
        pilar_[j].vecMulAddTriangle3D_AVX(X, Y, 3*j);

#elif ( SD_BLOCK_SIZE == 3 ) && SMSBD_USES_SSE && !REAL_IS_DOUBLE
    
    // process diagonal:
    vecMulDiagonal3D_SSE(X, Y);
    
    // process off-diagonal elements:
    for ( index_t j = colix_[0]; j < rsize_; j = colix_[j+1] )
        pilar_[j].vecMulAddTriangle3D_SSE(X, Y, 3*j);
    
#elif ( SD_BLOCK_SIZE == 0 )
    
    // process diagonal:
    vecMulDiagonal3D(X, Y);
    
    // process off-diagonal elements:
    for ( index_t j = colix_[0]; j < rsize_; j = colix_[j+1] )
        pilar_[j].vecMulAddTriangle3D(X, Y, 3*j);
    
#elif ( SD_BLOCK_SIZE == 2 )
    
    // process diagonal:
    vecMulDiagonal2D(X, Y);
    
    // process off-diagonal elements:
    for ( index_t j = colix_[0]; j < rsize_; j = colix_[j+1] )
        pilar_[j].vecMulAddTriangle2D(X, Y, 2*j);
    
#else
    
    zero_real(SD_BLOCK_SIZE*rsize_, Y);
    vecMulAdd(X, Y, 0, rsize_);
    //for ( size_t j = 0; j < rsize_; ++j ) pilar_[j].vecMulAdd3D_SIMD(X, Y, j*3);

#endif
}
