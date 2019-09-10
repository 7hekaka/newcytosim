// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsparseblk.h"
#include "assert_macro.h"
#include <sstream>

#ifndef DIM
#  define DIM BLOCK_SIZE
#endif

#include "vector.h"

// Flag to enable AVX implementation
#ifdef __AVX2__
#  define MATRIXSB_USES_AVX REAL_IS_DOUBLE
#else
#  define MATRIXSB_USES_AVX 0
#endif


MatrixSparseBlock::MatrixSparseBlock()
{
    allocated_ = 0;
    row_       = nullptr;
    blocks_    = nullptr;
    
    next_      = new index_t[1];
    next_[0]   = 0;
}


void MatrixSparseBlock::allocate(size_t alc)
{
    if ( alc > allocated_ )
    {
        /*
         'chunk' can be increased to gain performance:
          more memory will be used, but reallocation will be less frequent
        */
        constexpr size_t chunk = 16;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "MSB allocates %u\n", alc);
        Line * row_new = new Line[alc];
       
        if ( row_ )
        {
            // using the specialy defined '=' for a Line object
            for (size_t n = 0; n < allocated_; ++n )
                row_new[n] = row_[n];
            delete[] row_;
        }
        
        row_       = row_new;
        allocated_ = alc;
        
        delete[] next_;
        next_ = new index_t[allocated_+1];
    }
}


void MatrixSparseBlock::deallocate()
{
    delete[] row_;
    delete[] next_;
    row_ = nullptr;
    next_ = nullptr;
    free_real(blocks_);
    allocated_ = 0;
}


void MatrixSparseBlock::Line::allocate(size_t alc)
{
    if ( alc > allo_ )
    {
        //fprintf(stderr, "MSB reallocates column %i for %u: %p\n", inx_[0], alc);
        //else fprintf(stderr, "MSB allocates column for %u: %p\n", alc);
        /*
         'chunk' can be increased, to possibly gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr size_t chunk = 32;
        alc = ( alc + chunk - 1 ) & ~( chunk - 1 );
        
        // use aligned memory:
        void * ptr = new_real(alc*sizeof(SubBlock)/sizeof(real));
        SubBlock * blk_new  = new(ptr) SubBlock[alc];

        posix_memalign(&ptr, 32, alc*sizeof(index_t));
        index_t * inx_new = (index_t*)ptr;

        if ( inx_ )
        {
            for ( index_t n = 0; n < size_; ++n )
                inx_new[n] = inx_[n];
            free(inx_);
        }

        if ( blk_ )
        {
            for ( index_t n = 0; n < size_; ++n )
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


void MatrixSparseBlock::Line::deallocate()
{
    //if ( inx_ ) fprintf(stderr, "MSB deallocates column %i\n", inx_[0]);
    free(inx_);
    free_real(blk_);
    inx_ = nullptr;
    blk_ = nullptr;
}


void MatrixSparseBlock::Line::operator =(MatrixSparseBlock::Line & row)
{
    //if ( inx_ ) fprintf(stderr, "MSB transfers column %u\n", inx_[0]);
    free(inx_);
    free_real(blk_);

    size_ = row.size_;
    allo_ = row.allo_;
    inx_ = row.inx_;
    blk_ = row.blk_;
    sbk_ = nullptr;

    row.size_ = 0;
    row.allo_ = 0;
    row.inx_ = nullptr;
    row.blk_ = nullptr;
}

/**
 This allocate to be able to hold the matrix element if necessary
 */
SubBlock& MatrixSparseBlock::Line::block(index_t jj)
{
    index_t n = 0;
    /* This is a silly search that could be optimized */
    for ( n = 0; n < size_; ++n )
        if ( inx_[n] == jj )
            return blk_[n];

    allocate(n+1);
    
    //add the requested term:
    inx_[n] = jj;
    blk_[n].reset();
    size_ = n + 1;
    return blk_[n];
}


void MatrixSparseBlock::Line::reset()
{
    size_ = 0;
}


real& MatrixSparseBlock::operator()(index_t ii, index_t jj)
{
    assert_true( ii >= jj );
#if ( BLOCK_SIZE == 1 )
    return row_[ii].block(jj).value();
#else
    index_t i = ii % BLOCK_SIZE;
    index_t j = jj % BLOCK_SIZE;
    return row_[ii-i].block(jj-j)(i, j);
#endif
}


real* MatrixSparseBlock::addr(index_t ii, index_t jj) const
{
    assert_true( ii >= jj );
#if ( BLOCK_SIZE == 1 )
    return &row_[ii].block(jj).value();
#else
    index_t i = ii % BLOCK_SIZE;
    index_t j = jj % BLOCK_SIZE;
    return row_[ii-i].block(jj-j).addr(i, j);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

void MatrixSparseBlock::reset()
{
    already_symmetric = false;
    for ( index_t n = 0; n < size_; ++n )
        row_[n].reset();
}


bool MatrixSparseBlock::nonZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        Line & row = row_[jj];
        for ( unsigned n = 0 ; n < row.size_ ; ++n )
            if ( row[n] != 0.0 )
                return true;
    }
    //if here, the matrix is empty
    return false;
}


void MatrixSparseBlock::scale(const real alpha)
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        Line & row = row_[jj];
        for ( unsigned n = 0 ; n < row.size_ ; ++n )
            row[n].scale(alpha);
    }
}


void MatrixSparseBlock::addTriangularBlock(real* mat, const unsigned ldd,
                                           const index_t start,
                                           const unsigned cnt,
                                           const unsigned dim) const
{
    assert_false( start % BLOCK_SIZE );
    assert_false( cnt % BLOCK_SIZE );

    index_t end = start + cnt;
    index_t off = start + ldd * start;
    assert_true( end <= size_ );
    
    for ( index_t i = start; i < end; ++i )
    {
        Line & row = row_[i];
        for ( index_t n = 1; n < row.size_; ++n )
        {
            index_t j = row.inx_[n];
            if ( start <= j && j < end )
                row[n].addto(mat+(i+ldd*j)-off, ldd);
        }
    }
}


void MatrixSparseBlock::addDiagonalBlock(real* mat, unsigned ldd,
                                         const index_t start,
                                         const unsigned cnt) const
{
    assert_false( start % BLOCK_SIZE );
    assert_false( cnt % BLOCK_SIZE );

    index_t end = start + cnt;
    index_t off = start + ldd * start;
    assert_true( end <= size_ );
    
    for ( index_t i = start; i < end; ++i )
    {
        Line & row = row_[i];
        for ( index_t n = 0; n < row.size_; ++n )
        {
            index_t j = row.inx_[n];
            if ( start <= j && j < end )
                row[n].addto(mat+(i+ldd*j)-off, ldd);
        }
    }
}


int MatrixSparseBlock::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        Line & row = row_[jj];
        for ( unsigned n = 0 ; n < row.size_ ; ++n )
        {
            if ( row.inx_[n] >= size_ ) return 2;
            if ( row.inx_[n] <= jj )    return 3;
        }
    }
    return 0;
}


/** all allocated elements are counted, even if zero */
size_t MatrixSparseBlock::nbElements() const
{
    size_t cnt = 0;
    for ( index_t jj = 0; jj < size_; ++jj )
        cnt += row_[jj].size_;
    return cnt;
}


//------------------------------------------------------------------------------
#pragma mark -


std::string MatrixSparseBlock::what() const
{
    std::ostringstream msg;
#if MATRIXSB_USES_AVX
    msg << "MSBx " << SubBlock::what() << "*" << nbElements();
#elif defined(__SSE3__) &&  REAL_IS_DOUBLE
    msg << "MSBe " << SubBlock::what() << "*" << nbElements();
#else
    msg << "MSB " << SubBlock::what() << "*" << nbElements();
#endif
    return msg.str();
}


void MatrixSparseBlock::printSparse(std::ostream& os) const
{
    char str[256];
    std::streamsize p = os.precision();
    os.precision(8);
    if ( ! row_ )
        return;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        Line & row = row_[jj];
        if ( row.size_ > 0 )
            os << "% column " << jj << "\n";
        for ( unsigned n = 0 ; n < row.size_ ; ++n )
        {
            index_t ii = row.inx_[n];
            SubBlock blk = row.blk_[n];
            int d = ( ii == jj );
            for ( int x = 0  ; x < BLOCK_SIZE; ++x )
            for ( int y = x*d; y < BLOCK_SIZE; ++y )
            {
                real v = blk(y, x);
                if ( v != 0 )
                {
                    snprintf(str, sizeof(str), "%6i %6i %16.6f\n", ii+y, jj+x, blk(y, x));
                    os << str;
                }
            }
        }
    }
    os.precision(p);
}


void MatrixSparseBlock::printLines(std::ostream& os)
{
    os << "MSB size " << size_ << ":";
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        os << "\n   " << jj << "   " << row_[jj].size_;
        os << " " << next_[jj];
    }
    std::endl(os);
}


void MatrixSparseBlock::Line::print(std::ostream& os) const
{
    for ( index_t n = 0; n < size_; ++n )
        os << "\n" << inx_[n] << " : " << blk_[n] << "\n";
    std::endl(os);
}


//------------------------------------------------------------------------------
#pragma mark - Prepare for Multiplication


/// A block element of the sparse matrix suitable for qsort()
class MatrixSparseBlock::Element
{
public:
    /// block element
    SubBlock blk;
    
    /// index
    index_t inx;
};


/// function for qsort, comparing line indices
int compareMSBElement(const void * p, const void * q)
{
    MatrixSparseBlock::Element const* a = (MatrixSparseBlock::Element const*)(p);
    MatrixSparseBlock::Element const* b = (MatrixSparseBlock::Element const*)(q);
    
    if ( a->inx > b->inx ) return  1;
    if ( a->inx < b->inx ) return -1;
    return 0;
}


size_t MatrixSparseBlock::newElements(MatrixSparseBlock::Element*& ptr, size_t size)
{
    constexpr size_t chunk = 16;
    size_t all = ( size + chunk - 1 ) & ~( chunk - 1 );
    free(ptr);  // Element has no destructor 
    void* tmp = nullptr;
    if ( size > 0 )
    {
        posix_memalign(&tmp, 32, all * sizeof(MatrixSparseBlock::Element));
        ptr = new(tmp) MatrixSparseBlock::Element[all];
    }
    else
        ptr = nullptr;
    return all;
}

/**
 This copies the data to the provided temporary array
 */
void MatrixSparseBlock::Line::sort(Element*& tmp, size_t tmp_size)
{
    assert_true( size_ <= tmp_size );
    for ( unsigned i = 0; i < size_; ++i )
    {
        tmp[i].blk = blk_[i];
        tmp[i].inx = inx_[i];
    }
    
    //std::clog << "sizeof(Element) " << sizeof(Element) << "\n";
    qsort(tmp, size_, sizeof(Element), &compareMSBElement);
    
    for ( unsigned i = 0; i < size_; ++i )
    {
         blk_[i] = tmp[i].blk;
         inx_[i] = tmp[i].inx;
    }
}


void MatrixSparseBlock::sortElements()
{
    //unsigned cnt = 0;
    size_t tmp_size = 0;
    Element * tmp = nullptr;
    
    for ( index_t i = next_[0]; i < size_; i = next_[i+1] )
    {
        assert_true( i < size_ );
        Line & row = row_[i];
        assert_true( row.size_ > 0 );
        //std::clog << "MSB line " << jj << " has " << row.size_ << " elements\n";
        
        // order the elements in each line:
        if ( row.size_ > 1 )
        {
            if ( tmp_size < row.size_ )
                tmp_size = newElements(tmp, row.size_);
            row.sort(tmp, tmp_size);
        }
        
        //++cnt;
    }
    
    newElements(tmp, 0);
    //std::clog << "MatrixSparseBlock " << size_ << " with " << cnt << " non-empty columns\n";
}


void MatrixSparseBlock::consolidate()
{
    size_t cnt = 0;
    
    for ( index_t i = next_[0]; i < size_; i = next_[i+1] )
    {
        cnt += row_[i].size_;
        //std::cerr << "\nMatrixSparseBlock line " << i << "  " << row.size_ << "  " << row.blk_ << "";
    }
    
    //std::cerr << "\nMatrixSparseBlock:consolidate with " << cnt << " blocks";

    free_real(blocks_);
    real * ptr = new_real(cnt*sizeof(SubBlock)/sizeof(real));
    blocks_ = new(ptr) SubBlock[cnt];
    
    cnt = 0;
    for ( index_t i = 0; i < size_; ++i )
    {
        Line & row = row_[i];
        row.sbk_ = blocks_ + cnt;
        cnt += row.size_;
        for ( size_t j = 0; j < row.size_; ++j )
            row.sbk_[j] = row.blk_[j];
    }
}


/**
 Copy data from the lower triangle to the upper triangle
 */
void MatrixSparseBlock::symmetrize()
{
    for ( index_t i = next_[0]; i < size_; i = next_[i+1] )
    {
        Line & row = row_[i];
        //std::clog << "MSB line " << i << " has " << row.size_ << " elements\n";
        
        for ( unsigned n = 0 ; n < row.size_ ; ++n )
        {
            /// we duplicate blocks from the lower triangle
            index_t j = row.inx_[n];
            assert_true( j <= i );
            if ( j < i )
            {
                //std::cerr << "copying block at " << i << ", " << j << "\n";
                assert_true( row_[j].size_ > 0 );
                row_[j].block(i) = row.blk_[n].transposed();
            }
            else if ( i == j )
                row.blk_[n].copy_lower();
        }
    }
    
#if 0
    /// check that indices are in ascending order:
    for ( index_t i = 0; i < size_; ++i )
    {
        Line & row = row_[i];
        //std::clog << "MSB line " << i << " has " << row.size_ << " elements\n";
        index_t j = 0;
        for ( unsigned n = 0 ; n < row.size_ ; ++n )
        {
            if ( row.inx_[n] < j )
                std::clog << "MSB line " << i << " is disordered\n";
            j = row.inx_[n];
        }
    }
#endif
}


void MatrixSparseBlock::prepareForMultiply(int)
{
    next_[size_] = size_;
    
    if ( size_ > 0 )
    {
        index_t inx = size_;
        index_t nxt = size_;
        while ( --inx > 0 )
        {
            if ( row_[inx].size_ > 0 )
                nxt = inx;
            next_[inx] = nxt;
        }
        if ( row_[0].size_ > 0 )
            next_[0] = 0;
        else
            next_[0] = nxt;
    }
    
    sortElements();
    if ( !already_symmetric )
    {
        //std::cerr << "\nMatrixSparseBlock:symmetrize " << nbElements();
        symmetrize();
        //std::cerr << " -> " << nbElements() << "  ";
        already_symmetric = true;
    }
    
    consolidate();
}


//------------------------------------------------------------------------------
#pragma mark - Basic Vector Multiplication


void MatrixSparseBlock::Line::vecMulAdd(const real* X, real* Y) const
{
    Vector res(0,0,0);
    for ( index_t n = 0; n < size_; ++n )
        res += blk_[n] * Vector(X+inx_[n]);
    res.add_to(Y);
}


// multiplication of a vector: Y = Y + M * X
void MatrixSparseBlock::vecMulAdd_SCAL(const real* X, real* Y, index_t start, index_t end) const
{
    assert_true( start <= end );
    assert_true( end <= size_ );
    for ( index_t i = next_[start]; i < end; i = next_[i+1] )
        row_[i].vecMulAdd(X, Y+i);
}

//------------------------------------------------------------------------------
#pragma mark - Manually Optimized Vector Multiplication

#if MATRIXSB_USES_AVX

#include "simd.h"

#if ( BLOCK_SIZE == -2 )
void MatrixSparseBlock::Line::vecMulAdd2D(const real* X, real* Y) const
{
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(load4(blk_[0]), xyxy);
    const vec4 xxyy = permute4(xyxy, 0b1100);
    vec4 s1 = setzero4();

    unsigned n = 1;
    const unsigned stop = 1 + 2 * ( ( size_ - 1 ) / 2 );
    // process 4 by 4:
    for ( ; n < stop; n += 2 )
    {
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        const index_t i0 = inx_[n  ];
        const index_t i1 = inx_[n+1];
        vec4 mat0 = load4(blk_[n  ]);
        vec4 mat1 = load4(blk_[n+1]);
        vec4 u0 = fmadd4(mat0, xxyy, cast4(load2(Y+i0)));
        vec4 u1 = fmadd4(mat1, xxyy, cast4(load2(Y+i1)));
        ss = fmadd4(mat0, broadcast2(X+i0), ss);
        s1 = fmadd4(mat1, broadcast2(X+i1), s1);
        store2(Y+i0, add2(getlo(u0), gethi(u0)));
        store2(Y+i1, add2(getlo(u1), gethi(u1)));
    }
    // collapse 'ss'
    ss = add4(ss, s1);
    // process remaining blocks:
    for ( ; n < size_; ++n )
        multiply2D(X, Y, inx_[n], load4(blk_[n]), xxyy, ss);
    /* finally horizontally sum ss = { SX SX SY SY } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    return add2(load2(Y+jj), h);
}
#endif


#if ( BLOCK_SIZE == -2 )
void MatrixSparseBlock::Line::vecMulAdd2DU(const real* X, real* Y) const
{
    vec4 xyxy = broadcast2(X+jj);
    vec4 ss = mul4(load4(blk_[0]), xyxy);
    const vec4 xxyy = permute4(xyxy, 0b1100);
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();

    unsigned n = 1;
    const unsigned stop = 1 + 4 * ( ( size_ - 1 ) / 4 );
    // process 4 by 4:
    for ( ; n < stop; n += 4 )
    {
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        const index_t i0 = inx_[n  ];
        const index_t i1 = inx_[n+1];
        const index_t i2 = inx_[n+2];
        const index_t i3 = inx_[n+3];
        vec4 mat0 = load4(blk_[n  ]);
        vec4 mat1 = load4(blk_[n+1]);
        vec4 mat2 = load4(blk_[n+2]);
        vec4 mat3 = load4(blk_[n+3]);
        vec4 u0 = fmadd4(mat0, xxyy, cast4(load2(Y+i0)));
        vec4 u1 = fmadd4(mat1, xxyy, cast4(load2(Y+i1)));
        vec4 u2 = fmadd4(mat2, xxyy, cast4(load2(Y+i2)));
        vec4 u3 = fmadd4(mat3, xxyy, cast4(load2(Y+i3)));
        ss = fmadd4(mat0, broadcast2(X+i0), ss);
        s1 = fmadd4(mat1, broadcast2(X+i1), s1);
        s2 = fmadd4(mat2, broadcast2(X+i2), s2);
        s3 = fmadd4(mat3, broadcast2(X+i3), s3);
        store2(Y+i0, add2(getlo(u0), gethi(u0)));
        store2(Y+i1, add2(getlo(u1), gethi(u1)));
        store2(Y+i2, add2(getlo(u2), gethi(u2)));
        store2(Y+i3, add2(getlo(u3), gethi(u3)));
    }
    // collapse 'ss'
    ss = add4(add4(ss,s1), add4(s2,s3));
    // process remaining blocks:
    for ( ; n < size_; ++n )
        multiply2D(X, Y, inx_[n], load4(blk_[n]), xxyy, ss);
    /* finally sum ss = { S0 S0 S1 S1 } */
    vec2 h = gethi(ss);
    h = add2(unpacklo2(getlo(ss), h), unpackhi2(getlo(ss), h));
    store2(Y+jj, add2(load2(Y+jj), h));
}
#endif


#if ( BLOCK_SIZE == 3 )
void MatrixSparseBlock::Line::vecMulAdd3D(const real* X, real* Y) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    // There is a dependency in the loop for 's0', 's1' and 's2'.
    const real* M = blk_[0];
    const real* stop = blk_[size_];
    const index_t * inx = inx_;
    #pragma nounroll
    for ( ; M < stop; M += 12 )
    {
        vec4 xyz = loadu4(X+inx[0]);  // xyz = { X0 X1 X2 - }
        ++inx;
        // multiply with the block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        s0 = fmadd4(load4(M  ), xyz, s0);
        s1 = fmadd4(load4(M+4), xyz, s1);
        s2 = fmadd4(load4(M+8), xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
    vec4 s3 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s1 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    s0 = add4(permute2f128(s0, s1, 0x20), permute2f128(s0, s1, 0x31));
    storeu4(Y, add4(loadu4(Y), s0));
}
#endif


//#include "iacaMarks.h"

#if ( BLOCK_SIZE == 3 )
void MatrixSparseBlock::Line::vecMulAdd3DU(const real* X, real* Y) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 t0 = setzero4();
    vec4 t1 = setzero4();
    vec4 t2 = setzero4();

    assert_true( sizeof(SubBlock) == 12 * sizeof(real) );
    const real* M = sbk_[0];
    const real* stop = sbk_[size_-size_%2];
    const index_t * inx = inx_;
    {
        /*
         Unrolling will reduce the dependency chain but the number of registers
         (16 for AVX CPU) may limit the level of unrolling that can be done.
         Moreover, the bottleneck here is the high number of loads needed
         to advance the calculation. AVX-512 loads & muls would work well here.
         */
        // process blocks 2 by 2:
        #pragma nounroll
        for ( ; M < stop; M += 24 )
        {
            //IACA_START
            vec4 A = loadu4(X+inx[0]);
            vec4 B = loadu4(X+inx[1]);
            inx += 2;
            // multiply each line of the two blocks:
            s0 = fmadd4(load4(M   ), A, s0);
            s1 = fmadd4(load4(M+4 ), A, s1);
            s2 = fmadd4(load4(M+8 ), A, s2);
            t0 = fmadd4(load4(M+12), B, t0);
            t1 = fmadd4(load4(M+16), B, t1);
            t2 = fmadd4(load4(M+20), B, t2);
        }
        //IACA_END
        s0 = add4(s0, t0);
        s1 = add4(s1, t1);
        s2 = add4(s2, t2);
    }
    // process remaining blocks:
    stop = sbk_[size_];
    #pragma nounroll
    for ( ; M < stop; M += 12 )
    {
        vec4 xyz = loadu4(X+inx[0]);  // xyz = { X0 X1 X2 - }
        ++inx;
        // multiply with the block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        s0 = fmadd4(load4(M  ), xyz, s0);
        s1 = fmadd4(load4(M+4), xyz, s1);
        s2 = fmadd4(load4(M+8), xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
    t0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s1 = add4(unpacklo4(s2, t0), unpackhi4(s2, t0));
    s0 = add4(permute2f128(s0, s1, 0x20), permute2f128(s0, s1, 0x31));
    storeu4(Y, add4(loadu4(Y), s0));
}
#endif


#if ( BLOCK_SIZE == 3 )
void MatrixSparseBlock::Line::vecMulAdd3DU4(const real* X, real* Y) const
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

    const real* M = sbk_[0];
    const real* stop = sbk_[size_-size_%3];
    const index_t * inx = inx_;
    {
        /*
         Unrolling will reduce the dependency chain but the number of registers
         (16 for AVX CPU) may limit the level of unrolling that can be done.
         Moreover, the bottleneck here is the high number of loads needed
         to advance the calculation. AVX-512 loads & muls would work well here.
         */
        // process blocks 3 by 3:
        #pragma nounroll
        for ( ; M < stop; M += 36 )
        {
            vec4 A = loadu4(X+inx[0]);
            vec4 B = loadu4(X+inx[1]);
            vec4 C = loadu4(X+inx[2]);
            inx += 3;
            // multiply each line of the two blocks:
            s0 = fmadd4(load4(M   ), A, s0);
            s1 = fmadd4(load4(M+4 ), A, s1);
            s2 = fmadd4(load4(M+8 ), A, s2);
            t0 = fmadd4(load4(M+12), B, t0);
            t1 = fmadd4(load4(M+16), B, t1);
            t2 = fmadd4(load4(M+20), B, t2);
            u0 = fmadd4(load4(M+24), C, u0);
            u1 = fmadd4(load4(M+28), C, u1);
            u2 = fmadd4(load4(M+32), C, u2);
        }
        s0 = add4(s0, add4(t0, u0));
        s1 = add4(s1, add4(t1, u1));
        s2 = add4(s2, add4(t2, u2));
    }
    // process remaining blocks:
    stop = sbk_[size_];
    #pragma nounroll
    for ( ; M < stop; M += 12 )
    {
        vec4 xyz = loadu4(X+inx[0]);  // xyz = { X0 X1 X2 - }
        ++inx;
        // multiply with the block:
        //Y0 += M[0] * X[ii] + M[1] * X[ii+1] + M[2] * X[ii+2];
        //Y1 += M[3] * X[ii] + M[4] * X[ii+1] + M[5] * X[ii+2];
        //Y2 += M[6] * X[ii] + M[7] * X[ii+1] + M[8] * X[ii+2];
        s0 = fmadd4(load4(M  ), xyz, s0);
        s1 = fmadd4(load4(M+4), xyz, s1);
        s2 = fmadd4(load4(M+8), xyz, s2);
    }
    // finally sum s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
    t0 = setzero4();
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s1 = add4(unpacklo4(s2, t0), unpackhi4(s2, t0));
    s0 = add4(permute2f128(s0, s1, 0x20), permute2f128(s0, s1, 0x31));
    storeu4(Y, add4(loadu4(Y), s0));
}
#endif


#if ( BLOCK_SIZE == 4 )
void MatrixSparseBlock::Line::vecMulAdd4D(const real* X, real* Y) const
{
    vec4 s0 = setzero4();
    vec4 s1 = setzero4();
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();

    // There is a dependency in the loop for 's0', 's1' and 's2'.
    for ( index_t n = 0; n < size_; ++n )
    {
        real const* M = blk_[n];
        const vec4 xyz = load4(X+inx_[n]);  // xyzt = { X0 X1 X2 X3 }
        s0 = fmadd4(load4(M   ), xyz, s0);
        s1 = fmadd4(load4(M+ 4), xyz, s1);
        s2 = fmadd4(load4(M+ 8), xyz, s2);
        s3 = fmadd4(load4(M+12), xyz, s3);
    }
    // finally sum s0 = { Y0 Y0 Y0 Y0 }, s1 = { Y1 Y1 Y1 Y1 }, s2 = { Y2 Y2 Y2 Y2 }
    s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
    s1 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
    s0 = add4(permute2f128(s0, s1, 0x20), permute2f128(s0, s1, 0x31));
    store4(Y, add4(load4(Y), s0));
}
#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication

// multiplication of a vector: Y = Y + M * X
void MatrixSparseBlock::vecMulAdd(const real* X, real* Y, index_t start, index_t end) const
{
    for ( index_t i = next_[start]; i < end; i = next_[i+1] )
    {
#if MATRIXSB_USES_AVX
#if ( DIM == 1 )
        row_[i].vecMulAdd(X, Y+i);
#elif ( DIM == 2 )
        row_[i].vecMulAdd2D(X, Y+i);
#elif ( DIM == 3 )
        row_[i].vecMulAdd3D(X, Y+i);
#endif
#else
        row_[i].vecMulAdd(X, Y+i);
#endif
    }
}


// multiplication of a vector: Y = Y + M * X
void MatrixSparseBlock::vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t end) const
{
    for ( index_t i = next_[start]; i < end; i = next_[i+1] )
    {
#if MATRIXSB_USES_AVX
#if ( DIM == 1 )
        row_[i].vecMulAdd(X, Y+i);
#elif ( DIM == 2 )
        row_[i].vecMulAdd2D(X, Y+i);
#elif ( DIM == 3 )
        row_[i].vecMulAdd3DU4(X, Y+i);
#endif
#else
        row_[i].vecMulAdd(X, Y+i);
#endif
    }
}

// multiplication of a vector: Y = Y + M * X
void MatrixSparseBlock::vecMulAdd_TIME(const real* X, real* Y, index_t start, index_t end) const
{
    unsigned long cnt = 0, row = 0;
    unsigned long long time = __rdtsc();
    for ( index_t i = next_[start]; i < end; i = next_[i+1] )
    {
        row++;
        cnt += row_[i].size_;
        row_[i].vecMulAdd(X, Y+i);
    }
    if ( cnt > 0 )
        fprintf(stderr, "MSB %6u with %6lu blocks/column  cycles/block: %6llu\n", size_, cnt/row, (__rdtsc()-time)/cnt);
}

