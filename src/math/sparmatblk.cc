// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include "sparmatblk.h"
#include "assert_macro.h"
#include <sstream>

#define TRANSPOSE_2D_BLOCKS 0


SparMatBlk::SparMatBlk()
{
    rsize_  = 0;
    alloc_  = 0;
    row_    = nullptr;
    blocks_ = nullptr;
    colidx_ = new unsigned[2]();
    is_symmetric = false;
}


void SparMatBlk::allocate(index_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be increased to gain performance:
          more memory will be used, but reallocation will be less frequent
        */
        constexpr index_t chunk = 16;
        alc = ( alc + chunk - 1 ) & ~( chunk - 1 );

        //fprintf(stderr, "SMB allocates %lu\n", alc);
        Line * ptr = new Line[alc];
       
        if ( row_ )
        {
            // using the custom '=' for a Line object
            for ( index_t n = 0; n < alloc_; ++n )
                ptr[n] = row_[n];
            delete[] row_;
        }
        
        row_ = ptr;
        alloc_ = alc;
        
        delete[] colidx_;
        colidx_ = new unsigned[alc+2];
        for ( unsigned n = 0; n <= alc; ++n )
            colidx_[n] = n;
    }
}


void SparMatBlk::deallocate()
{
    delete[] row_;
    delete[] colidx_;
    free_real(blocks_);
    row_ = nullptr;
    colidx_ = nullptr;
    blocks_ = nullptr;
    alloc_ = 0;
}


void SparMatBlk::Line::allocate(index_t alc)
{
    if ( alc > allo_ )
    {
        //fprintf(stderr, "SMB reallocates line %i for %u: %p\n", inx_[0], alc);
        //else fprintf(stderr, "SMB allocates line for %u: %p\n", alc);
        /*
         'chunk' can be increased, to possibly gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr index_t chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk - 1 );
        
        // use aligned memory:
        void * ptr = new_real(alc*sizeof(Block)/sizeof(real)+4);
        Block * blk_new = new(ptr) Block[alc];

        if ( posix_memalign(&ptr, 32, alc*sizeof(index_t)) )
            throw std::bad_alloc();
        auto * inx_new = (index_t*)ptr;

        if ( inx_ )
        {
            for ( index_t n = 0; n < rlen_; ++n )
                inx_new[n] = inx_[n];
            free(inx_);
        }

        if ( blk_ )
        {
            for ( index_t n = 0; n < rlen_; ++n )
                blk_new[n] = blk_[n];
            free_real(blk_);
        }
        inx_  = inx_new;
        blk_  = blk_new;
        allo_ = alc;
        
        //std::clog << "Line " << this << "  " << alc << ": ";
        //std::clog << " alignment " << ((uintptr_t)elem_ & 63) << "\n";
    }
}


void SparMatBlk::Line::deallocate()
{
    //if ( inx_ ) fprintf(stderr, "SMB deallocates column %lu of size %lu\n", inx_[0], allo_);
    free(inx_);
    free_real(blk_);
    inx_ = nullptr;
    blk_ = nullptr;
    allo_ = 0;
    rlen_ = 0;
}


void SparMatBlk::Line::operator = (SparMatBlk::Line & row)
{
    //if ( inx_ ) fprintf(stderr, "SMB transfers line %u\n", inx_[0]);
    free(inx_);
    free_real(blk_);

    rlen_ = row.rlen_;
    allo_ = row.allo_;
    inx_ = row.inx_;
    blk_ = row.blk_;
    sbk_ = nullptr;

    row.rlen_ = 0;
    row.allo_ = 0;
    row.inx_ = nullptr;
    row.blk_ = nullptr;
}


SparMatBlk::Block* SparMatBlk::Line::find_block(index_t jj) const
{
    /* This is a silly search that could be optimized */
    for ( index_t n = 0; n < rlen_; ++n )
        if ( inx_[n] == jj )
            return blk_ + n;
    return nullptr;
}

/**
 This allocates to be able to hold the matrix element if necessary
 */
SparMatBlk::Block& SparMatBlk::Line::block(index_t jj)
{
    SparMatBlk::Block * B = find_block(jj);
    if ( !B )
    {
        allocate(rlen_+1);
        //add the requested term:
        B = blk_ + rlen_;
        inx_[rlen_] = jj;
        B->reset();
        ++rlen_;
    }
    return *B;
}


void SparMatBlk::Line::reset()
{
    rlen_ = 0;
}


real& SparMatBlk::element(index_t ii, index_t jj)
{
    assert_true( ii >= jj );
#if ( BLOCK_SIZE == 1 )
    return row_[ii].block(jj).value();
#else
    index_t i = ii / BLOCK_SIZE;
    index_t j = jj / BLOCK_SIZE;
    assert_true( i < rsize_ );
    return row_[i].block(j)(ii%BLOCK_SIZE, jj%BLOCK_SIZE);
#endif
}


real* SparMatBlk::address(index_t ii, index_t jj) const
{
    assert_true( ii >= jj );
#if ( BLOCK_SIZE == 1 )
    return row_[ii].block(jj).data();
#else
    index_t i = ii / BLOCK_SIZE;
    index_t j = jj / BLOCK_SIZE;
    Block * B = row_[i].find_block(j);
    if ( B )
        return B->addr(ii%BLOCK_SIZE, jj%BLOCK_SIZE);
    return nullptr;
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

void SparMatBlk::reset()
{
    is_symmetric = false;
    for ( index_t n = 0; n < rsize_; ++n )
        row_[n].reset();
}


bool SparMatBlk::notZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Line & row = row_[jj];
        for ( index_t n = 0 ; n < row.rlen_ ; ++n )
            if ( row[n] != 0.0 )
                return true;
    }
    //if here, the matrix is empty
    return false;
}


void SparMatBlk::scale(const real alpha)
{
    for ( index_t jj = 0; jj < rsize_; ++jj )
    {
        Line & row = row_[jj];
        for ( index_t n = 0 ; n < row.rlen_ ; ++n )
            row[n].scale(alpha);
    }
}


void SparMatBlk::addDiagonalBlock(real* mat, index_t ldd, const index_t start, const index_t cnt, index_t mul) const
{
    assert_true( mul == BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );
    index_t end = start + cnt;

    for ( index_t ii = start; ii < end; ++ii )
    {
        index_t is = ii - start;
        Line & row = row_[ii];
        for ( index_t n = 0; n < row.rlen_; ++n )
        {
            // not assuming anything to be safe, since the matrix may be symmetric or not
            auto jj = row.inx_[n];
            if (( start <= jj ) & ( jj < end ))
            {
                auto js = jj - start;
                if ( ii == jj )
                    row[n].addto_symm(mat+(is+ldd*js)*BLOCK_SIZE, ldd);
                else
                    row[n].addto(mat+(is+ldd*js)*BLOCK_SIZE, ldd);
            }
        }
    }
}


void SparMatBlk::addLowerBand(real alpha, real* mat, index_t ldd, index_t start, index_t cnt,
                              const index_t mul, const index_t rank) const
{
    assert_true( mul == BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );
    index_t end = start + cnt;

    for ( index_t ii = start; ii < end; ++ii )
    {
        index_t is = ii - start;
        Line & row = row_[ii];
        for ( index_t n = 0; n < row.rlen_; ++n )
        {
            auto jj = row.inx_[n];
            auto diff = ii > jj ? ii - jj : jj - ii;
            if ((start <= jj) & (jj < end) & (diff <= rank))
            {
                auto js = jj - start;
                if ( ii == jj )
                    row[n].addto_lower(mat+(is+ldd*js)*BLOCK_SIZE, ldd, alpha);
                else if ( ii > jj )
                    row[n].addto(mat+(is+ldd*js)*BLOCK_SIZE, ldd, alpha);
            }
        }
    }
}


void SparMatBlk::addDiagonalTrace(real alpha, real* mat, index_t ldd,
                                  const index_t start, const index_t cnt,
                                  const index_t mul, const index_t rank, const bool sym) const
{
    assert_true( mul == BLOCK_SIZE );
    assert_true( start + cnt <= rsize_ );
    index_t end = start + cnt;

    for ( index_t ii = start; ii < end; ++ii )
    {
        index_t is = ii - start;
        Line & row = row_[ii];
        for ( index_t n = 0; n < row.rlen_; ++n )
        {
            // not assuming anything to be safe:
            auto jj = row.inx_[n];
            auto diff = ii > jj ? ii - jj : jj - ii;
            if (( start <= jj ) & ( jj < end ) & ( diff <= rank ))
            {
                index_t js = jj - start;
                real a = alpha * row[n].trace();
                //fprintf(stderr, "SMB %4lu %4lu : %.4f\n", is, js, a);
                if ( ii >= jj || sym )
                    mat[is+ldd*js] += a;
            }
        }
    }
}



int SparMatBlk::bad() const
{
    for ( index_t i = 0; i < rsize_; ++i )
    {
        Line & row = row_[i];
        index_t sup = ( is_symmetric ? rsize_ : i );
        for ( index_t n = 0 ; n < row.rlen_ ; ++n )
        {
            index_t j = row.inx_[n];
            if ( j > sup )
                return 2;
        }
    }
    return 0;
}


size_t SparMatBlk::nbElements(index_t start, index_t stop, size_t& alc) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    alc = 0;
    index_t cnt = 0;
    for ( index_t i = start; i < stop; ++i )
    {
        cnt += row_[i].rlen_;
        alc += row_[i].allo_;
    }
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark -


std::string SparMatBlk::what() const
{
    size_t alc = 0;
    size_t cnt = nbElements(0, rsize_, alc);
    std::ostringstream msg;
#if SPARMATBLK_USES_AVX
    msg << "SMBx ";
#else
    msg << "SMB ";
#endif
    msg << Block::what() << "*" << cnt << " (" << alc << ")";
    return msg.str();
}


static void printSparseBlock(std::ostream& os, real inf, SparMatBlk::Block const& B, index_t ii, index_t jj)
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


void SparMatBlk::printSparse(std::ostream& os, real epsilon, index_t start, index_t stop) const
{
    os << "% SparMatBlk size " << rsize_ << ":\n";
    stop = std::min(stop, rsize_);
    std::streamsize p = os.precision(8);
    if ( ! row_ )
        return;
    for ( index_t jj = start; jj < stop; ++jj )
    {
        Line & row = row_[jj];
        if ( row.notEmpty() )
            os << "% line " << jj << "\n";
        for ( index_t n = 0 ; n < row.rlen_ ; ++n )
            printSparseBlock(os, epsilon, row.blk_[n], row.inx_[n], jj);
    }
    os.precision(p);
}


void SparMatBlk::printSummary(std::ostream& os)
{
    os << "SMB size " << rsize_ << ":";
    for ( index_t i = 0; i < rsize_; ++i )
        if ( row_[i].notEmpty() )
        {
            os << "\n   " << i << "   " << row_[i].rlen_;
            os << " index " << colidx_[i];
        }
    std::endl(os);
}


void SparMatBlk::Line::printBlocks(std::ostream& os) const
{
    for ( index_t n = 0; n < rlen_; ++n )
        os << " " << inx_[n] << " " << blk_[n];
}


void SparMatBlk::printBlocks(std::ostream& os) const
{
    for ( index_t j = 0; j < rsize_; ++j )
    {
        os << "\nSMB  col " << j;
        row_[j].printBlocks(os);
    }
    std::endl(os);
}

//------------------------------------------------------------------------------
#pragma mark - Prepare for Multiplication


/// A block element of the sparse matrix suitable for qsort()
class SparMatBlk::Element
{
public:
    /// block element
    real blk[BLOCK_SIZE*BLOCK_SIZE];

    /// index
    index_t inx;
};


/// qsort function comparing line indices
static int compareSMBElement(const void * A, const void * B)
{
    index_t a = static_cast<SparMatBlk::Element const*>(A)->inx;
    index_t b = static_cast<SparMatBlk::Element const*>(B)->inx;

    return ( a > b ) - ( b > a );
}


index_t SparMatBlk::newElements(SparMatBlk::Element*& ptr, index_t cnt)
{
    constexpr index_t chunk = 16;
    index_t all = ( cnt + chunk - 1 ) & ~( chunk - 1 );
    free(ptr);  // Element has no destructor 
    void* tmp = nullptr;
    if ( posix_memalign(&tmp, 32, all*sizeof(SparMatBlk::Element)) )
        throw std::bad_alloc();
    ptr = new(tmp) SparMatBlk::Element[all];
    return all;
}

//------------------------------------------------------------------------------
#pragma mark - Preparing operations

/**
 This copies the data to the provided temporary array
 */
void SparMatBlk::Line::sortElements(Element tmp[], index_t tmp_size)
{
    assert_true( rlen_ <= tmp_size );
    for ( index_t i = 0; i < rlen_; ++i )
    {
        blk_[i].store(tmp[i].blk);
        tmp[i].inx = inx_[i];
    }
    
    //std::clog << "sizeof(Element) " << sizeof(Element) << "\n";
    qsort(tmp, rlen_, sizeof(Element), &compareSMBElement);
    
    for ( index_t i = 0; i < rlen_; ++i )
    {
         blk_[i].load(tmp[i].blk);
         inx_[i] = tmp[i].inx;
    }
}

/**
 Sort elements in each line in order of increasing column index
 */
void SparMatBlk::sortElements()
{
    //index_t cnt = 0;
    index_t tmp_size = 0;
    Element * tmp = nullptr;
    
    for ( index_t i = colidx_[0]; i < rsize_; i = colidx_[i+1] )
    {
        assert_true( i < rsize_ );
        Line & row = row_[i];
        assert_true( row.rlen_ > 0 );
        //std::clog << "SMB line " << jj << " has " << row.rlen_ << " elements\n";
        
        // order the elements in each line:
        if ( row.rlen_ > 1 )
        {
            if ( tmp_size < row.rlen_ )
                tmp_size = newElements(tmp, row.rlen_);
            row.sortElements(tmp, tmp_size);
        }
        
        //++cnt;
    }
    
    free(tmp);
    //std::clog << "SparMatBlk " << rsize_ << " with " << cnt << " non-empty lines\n";
}


/**
 allocates a single chunk of memory to hold all the lines consecutively
 */
void SparMatBlk::consolidate()
{
    index_t cnt = 0;
    for ( index_t i = colidx_[0]; i < rsize_; i = colidx_[i+1] )
    {
        cnt += row_[i].rlen_;
        //std::cerr << "\nMatrixSparseBlock line " << i << "  " << row.rlen_ << "  " << row.blk_ << "";
    }
    
    //std::cerr << "\nMatrixSparseBlock:consolidate with " << cnt << " blocks";

    free_real(blocks_);
    real * ptr = new_real(cnt*sizeof(Block)/sizeof(real));
    blocks_ = new(ptr) Block[cnt];
    
    Block * B = blocks_;
    for ( index_t i = 0; i < rsize_; ++i )
    {
        Line & row = row_[i];
        row.sbk_ = B;
#if ( BLOCK_SIZE == 2 ) && TRANSPOSE_2D_BLOCKS
        for ( index_t j = 0; j < row.rlen_; ++j )
            *B++ = row.blk_[j].transposed();
#else
        for ( index_t j = 0; j < row.rlen_; ++j )
            *B++ = row.blk_[j];
#endif
    }
}


/**
 Copy data from the lower triangle to the upper triangle
 */
void SparMatBlk::symmetrize()
{
    for ( index_t i = colidx_[0]; i < rsize_; i = colidx_[i+1] )
    {
        Line & row = row_[i];
        //std::clog << "SMB line " << i << " has " << row.rlen_ << " elements\n";
        
        for ( index_t n = 0 ; n < row.rlen_ ; ++n )
        {
            /// we duplicate blocks below the diagonal:
            index_t j = row.inx_[n];
            if ( i == j )
                row.blk_[n].copy_lower();
            else {
                assert_true( i > j );
                //std::cerr << "copying block at " << i << ", " << j << "\n";
                row_[j].block(i) = row.blk_[n].transposed();
            }

        }
    }
    
#if 0
    /// check that indices are in ascending order:
    for ( index_t i = 0; i < rsize_; ++i )
    {
        Line & row = row_[i];
        //std::clog << "SMB line " << i << " has " << row.rlen_ << " elements\n";
        index_t j = 0;
        for ( index_t n = 0 ; n < row.rlen_ ; ++n )
        {
            if ( row.inx_[n] < j )
                std::clog << "SMB line " << i << " is disordered\n";
            j = row.inx_[n];
        }
    }
#endif
}


bool SparMatBlk::prepareForMultiply(int)
{
    assert_false(bad());
    colidx_[rsize_] = rsize_;
    if ( rsize_ > 0 )
    {
        index_t inx = rsize_;
        index_t nxt = rsize_;
        while ( inx-- > 0 )
        {
            if ( row_[inx].notEmpty() )
                nxt = inx;
            else
                row_[inx].deallocate();
            colidx_[inx] = nxt;
        }
    }

    // check if matrix is empty:
    if ( colidx_[0] == rsize_ )
        return false;

    sortElements();

    if ( !is_symmetric )
    {
        //std::cerr << "\nMatrixSparseBlock:symmetrize " << nbElements();
        symmetrize();
        //std::cerr << " -> " << nbElements() << "  ";
        is_symmetric = true;
    }
    
#if 0
    consolidate();
#else
    for ( index_t i = 0; i < rsize_; ++i )
        row_[i].sbk_ = row_[i].blk_;
#endif
    return true;
}


//------------------------------------------------------------------------------
#pragma mark - Basic Vector Multiplication

void SparMatBlk::Line::vecMulLine(const real* X, real* Y) const
{
    Vector vec(Y);
    for ( index_t n = 0; n < rlen_; ++n )
        vec += blk_[n].vecmul(X+BLOCK_SIZE*inx_[n]);
    vec.store(Y);
}


#if ( BLOCK_SIZE == 1 )
real SparMatBlk::Line::vecMul1D(const real* X) const
{
    real res = 0;
    for ( index_t n = 0; n < rlen_; ++n )
        res += blk_[n].value() * X[inx_[n]];
    return res;
}
#endif

//------------------------------------------------------------------------------
#pragma mark - SIMD Optimized Vector Multiplication

#if SPARMATBLK_USES_AVX

#include "simd.h"

#if ( BLOCK_SIZE == 2 )
vec2 SparMatBlk::Line::vecMul2D(const double* X) const
{
    vec4 ss = setzero4();
    Block const* blk = blk_;
    Block const* end = blk_ + rlen_;
    auto const* inx = inx_;
    #pragma nounroll
    for ( ; blk < end; ++blk )
    {
        vec4 xy = broadcast2(X+2*inx[0]);  // xy = { X Y }
        ++inx;
        //SX += M[0] * X + M[2] * Y;
        //SY += M[1] * X + M[3] * Y;
        //ss[0] += M[0] * xy[0];
        //ss[1] += M[1] * xy[0];
        //ss[2] += M[2] * xy[1];
        //ss[3] += M[3] * xy[1];
        ss = fmadd4(blk->data0(), duplohi4(xy), ss);
    }
    // collapse result:
    return add2(getlo(ss), gethi(ss));
}
#endif


#if ( BLOCK_SIZE == 2 )
vec2 SparMatBlk::Line::vecMul2DU(const double* X) const
{
    vec4 ss = setzero4();
    vec4 tt = setzero4();
    vec4 uu = setzero4();
    vec4 vv = setzero4();
    Block const* blk = sbk_;
    Block const* end = sbk_ + ( rlen_ & ~3 );
    auto const* inx = inx_;
    #pragma nounroll
    for ( ; blk < end; blk += 4, inx += 4 )
    {
        assert_true( inx[0] < inx[1] );
        assert_true( inx[1] < inx[2] );
        assert_true( inx[2] < inx[3] );
        vec4 xy0 = broadcast2(X+2*inx[0]);  // xy = { X Y }
        vec4 xy1 = broadcast2(X+2*inx[1]);  // xy = { X Y }
        vec4 xy2 = broadcast2(X+2*inx[2]);  // xy = { X Y }
        vec4 xy3 = broadcast2(X+2*inx[3]);  // xy = { X Y }
#if TRANSPOSE_2D_BLOCKS
        //SX += M[0] * X + M[1] * Y;
        //SY += M[2] * X + M[3] * Y;
        ss = fmadd4(blk[0].data0(), xy0, ss);
        tt = fmadd4(blk[1].data0(), xy1, tt);
        uu = fmadd4(blk[2].data0(), xy2, uu);
        vv = fmadd4(blk[3].data0(), xy3, vv);
#else
        //SX += M[0] * X + M[2] * Y;
        //SY += M[1] * X + M[3] * Y;
        ss = fmadd4(blk[0].data0(), duplohi4(xy0), ss);
        tt = fmadd4(blk[1].data0(), duplohi4(xy1), tt);
        uu = fmadd4(blk[2].data0(), duplohi4(xy2), uu);
        vv = fmadd4(blk[3].data0(), duplohi4(xy3), vv);
#endif
    }
    ss = add4(add4(ss, tt), add4(uu, vv));
    end = sbk_ +rlen_;
    #pragma nounroll
    for ( ; blk < end; ++blk, ++inx )
    {
        vec4 xy = broadcast2(X+2*inx[0]);  // xy = { X Y }
#if TRANSPOSE_2D_BLOCKS
        ss = fmadd4(blk[0].data0(), xy, ss);
#else
        ss = fmadd4(blk[0].data0(), duplohi4(xy), ss);
#endif
    }
    // collapse result:
#if TRANSPOSE_2D_BLOCKS
    vec2 h = gethi(ss);
    return add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
#else
    return add2(getlo(ss), gethi(ss));
#endif
}
#endif


#if ( BLOCK_SIZE == 3 )
vec4 SparMatBlk::Line::vecMul3D(const double* X) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    Block const* blk = blk_;
    Block const* end = blk_ + rlen_;
    auto const* inx = inx_;
    #pragma nounroll
    for ( ; blk < end; ++blk )
    {
        vec4 xyz = loadu4(X+3*inx[0]);  // xyz = { X0 X1 X2 - }
        ++inx;
        // multiply with the block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        s0 = fmadd4(blk->data0(), xyz, s0);
        s1 = fmadd4(blk->data1(), xyz, s1);
        s2 = fmadd4(blk->data2(), xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4 s3 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    return add4(catshift2(s0, s2), blend22(s0, s2));
}
#endif


#if ( BLOCK_SIZE == 3 )
vec4 SparMatBlk::Line::vecMul3DU(const double* X) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 t0 = setzero4();
    vec4 t1 = setzero4();
    vec4 t2 = setzero4();

    Block const* blk = sbk_;
    Block const* end = sbk_ + ( rlen_ & ~1 );
    auto const* inx = inx_;
    {
        /*
         Unrolling will reduce the dependency chain but the number of registers
         (16 for AVX CPU) may limit the level of unrolling that can be done.
         Moreover, the bottleneck here is the high number of loads needed
         to advance the calculation. AVX-512 loads & muls would work well here.
         */
        // process blocks 2 by 2:
        #pragma nounroll
        for ( ; blk < end; blk += 2, inx += 2 )
        {
            assert_true( inx[0] < inx[1] );
            vec4 A = loadu4(X+3*inx[0]);
            vec4 B = loadu4(X+3*inx[1]);
            // multiply each line of the two blocks:
            s0 = fmadd4(blk[0].data0(), A, s0);
            s1 = fmadd4(blk[0].data1(), A, s1);
            s2 = fmadd4(blk[0].data2(), A, s2);
            t0 = fmadd4(blk[1].data0(), B, t0);
            t1 = fmadd4(blk[1].data1(), B, t1);
            t2 = fmadd4(blk[1].data2(), B, t2);
        }
        s0 = add4(s0, t0);
        s1 = add4(s1, t1);
        s2 = add4(s2, t2);
    }
    // process remaining blocks:
    #pragma nounroll
    for ( end = sbk_ + rlen_ ; blk < end; ++blk )
    {
        vec4 xyz = loadu4(X+3*inx[0]);  // xyz = { X0 X1 X2 - }
        ++inx;
        // multiply with the block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        s0 = fmadd4(blk->data0(), xyz, s0);
        s1 = fmadd4(blk->data1(), xyz, s1);
        s2 = fmadd4(blk->data2(), xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    vec4 s3 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    return add4(catshift2(s0, s2), blend22(s0, s2));
}
#endif


#if ( BLOCK_SIZE == 3 )
vec4 SparMatBlk::Line::vecMul3DUU(const double* X) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 t0 = setzero4();
    vec4 t1 = setzero4();
    vec4 t2 = setzero4();
    vec4 u0 = setzero4();
    vec4 u1 = setzero4();
    vec4 u2 = setzero4();

    Block const* blk = sbk_;
    Block const* end = sbk_ + (rlen_-rlen_%3);
    auto const* inx = inx_;
    {
        /*
         Unrolling will reduce the dependency chain but the number of registers
         (16 for AVX CPU) may limit the level of unrolling that can be done.
         Moreover, the bottleneck here is the high number of loads needed
         to advance the calculation. AVX-512 loads & muls would work well here.
         */
        // process blocks 3 by 3:
        #pragma nounroll
        for ( ; blk < end; blk += 3 )
        {
            assert_true( inx[0] < inx[1] );
            assert_true( inx[1] < inx[2] );
            vec4 A = loadu4(X+3*inx[0]);
            vec4 B = loadu4(X+3*inx[1]);
            vec4 C = loadu4(X+3*inx[2]);
            inx += 3;
            // multiply each line of the two blocks:
            s0 = fmadd4(blk[0].data0(), A, s0);
            s1 = fmadd4(blk[0].data1(), A, s1);
            s2 = fmadd4(blk[0].data2(), A, s2);
            t0 = fmadd4(blk[1].data0(), B, t0);
            t1 = fmadd4(blk[1].data1(), B, t1);
            t2 = fmadd4(blk[1].data2(), B, t2);
            u0 = fmadd4(blk[2].data0(), C, u0);
            u1 = fmadd4(blk[2].data1(), C, u1);
            u2 = fmadd4(blk[2].data2(), C, u2);
        }
        s0 = add4(s0, add4(t0, u0));
        s1 = add4(s1, add4(t1, u1));
        s2 = add4(s2, add4(t2, u2));
    }
    // process remaining blocks:
    #pragma nounroll
    for ( end = sbk_ + rlen_; blk < end; ++blk )
    {
        vec4 xyz = loadu4(X+3*inx[0]);  // xyz = { X0 X1 X2 - }
        ++inx;
        // multiply with the block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        s0 = fmadd4(blk->data0(), xyz, s0);
        s1 = fmadd4(blk->data1(), xyz, s1);
        s2 = fmadd4(blk->data2(), xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
    s0 = clear4th(s0); // s0[3] is garbage
    s1 = clear4th(s1); // s1[3] is garbage
    s2 = clear4th(s2); // s2[3] is garbage
    t0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, t0), unpackhi4(s2, t0));
    return add4(catshift2(s0, s2), blend22(s0, s2));
}
#endif


#if ( BLOCK_SIZE == 4 )
vec4 SparMatBlk::Line::vecMul4D(const double* X) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();

    Block const* blk = blk_;
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    for ( index_t n = 0; n < rlen_; ++n, ++blk )
    {
        const vec4 xx = load4(X+4*inx_[n]);  // xyzt = { X0 X1 X2 X3 }
        s0 = fmadd4(blk->data0(), xx, s0);
        s1 = fmadd4(blk->data1(), xx, s1);
        s2 = fmadd4(blk->data2(), xx, s2);
        s3 = fmadd4(blk->data3(), xx, s3);
    }
    // finally sum s0 = { Y0 Y0 Y0 Y0 }, s1 = { Y1 Y1 Y1 Y1 }, s2 = { Y2 Y2 Y2 Y2 }
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    return add4(catshift2(s0, s2), blend22(s0, s2));
}
#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication

// multiplication of a vector: Y = Y + M * X
void SparMatBlk::vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    for ( index_t i = start; i < stop; ++i )
        row_[i].vecMulLine(X, Y+BLOCK_SIZE*i);
}


// multiplication of a vector: Y = Y + M * X
void SparMatBlk::vecMulAdd(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= rsize_ );
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    for ( index_t i = colidx_[start]; i < stop; i = colidx_[i+1] )
    {
#if ( BLOCK_SIZE == 1 )
        Y[i] += row_[i].vecMul1D(X);
#elif ( BLOCK_SIZE == 2 ) && SPARMATBLK_USES_AVX
        store2(Y+2*i, add2(load2(Y+2*i), row_[i].vecMul2DU(X)));
#elif ( BLOCK_SIZE == 3 ) && SPARMATBLK_USES_AVX
        // we need to use store3 only for the last line, if multithreaded
        store3(Y+3*i, add4(loadu4(Y+3*i), row_[i].vecMul3DU(X)));
#else
        row_[i].vecMulLine(X, Y+BLOCK_SIZE*i);
#endif
    }
}


// multiplication of a vector: Y = M * X
void SparMatBlk::vecMul(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, rsize_);
    //printf("msb %6i %6i : %p\n", start, stop, pthread_self());
    
    /** All values need to be reset since as the matrix is sparse,
     not every line will be addressed below */
    zero_real(BLOCK_SIZE*(stop-start), Y+BLOCK_SIZE*start);
    
    for ( index_t i = colidx_[start]; i < stop; i = colidx_[i+1] )
    {
#if ( BLOCK_SIZE == 1 )
        Y[i] = row_[i].vecMul1D(X);
#elif ( BLOCK_SIZE == 2 ) && SPARMATBLK_USES_AVX
        store2(Y+2*i, row_[i].vecMul2DU(X));
#elif ( BLOCK_SIZE == 3 ) && SPARMATBLK_USES_AVX
        // we need to use store3 only for the last line, if multithreaded
        store3(Y+3*i, row_[i].vecMul3DUU(X));
#else
        row_[i].vecMulLine(X, Y+BLOCK_SIZE*i);
#endif
    }
}
