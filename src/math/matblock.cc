// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matblock.h"
#include "cblas.h"

#include <iomanip>
#include <sstream>

//------------------------------------------------------------------------------
MatrixOfBlocks::MatrixOfBlocks()
{
    allocated_  = 0;
    block_     = nullptr;
    block_alc  = nullptr;
    block_cnt   = 0;
    block_size = nullptr;
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::deallocate()
{
    if ( block_ == nullptr ) return;
    
    for ( size_t ii=0; ii < block_cnt; ++ii )
        delete[] block_[ii];
    
    delete[] block_;
    delete[] block_size;
    delete[] block_alc;
    allocated_ = 0;
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::allocate(size_t nbb)
{
    assert_true( nbb > 0 );
    block_cnt = nbb;
    if ( nbb > allocated_ )
    {
        constexpr size_t chunk = 32;
        size_t nba = ( nbb + chunk - 1 ) & ~( chunk -1 );

        //printf("new block-matrix sz %i\n", nbblock );
        real** block_new = new real*[nba];
        
        size_t * block_alc_new  = new size_t[nba];
        size_t * block_size_new = new size_t[nba];
        
        size_t ii = 0;
        
        if ( block_ )
        {
            for ( ii = 0; ii < allocated_; ++ii )
            {
                block_new[ii]      = block_[ii];
                block_alc_new[ii]  = block_alc[ii];
                block_size_new[ii] = block_size[ii];
            }
            delete[] block_;
            delete[] block_alc;
            delete[] block_size;
         }
        
        for ( ; ii < nba; ++ii )
        {
            block_new[ii]      = nullptr;
            block_alc_new[ii]  = 0;
            block_size_new[ii] = 0;
        }
        
        block_      = block_new;
        block_alc   = block_alc_new;
        block_size  = block_size_new;
        allocated_  = nba;
    }
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::setBlockSize( const size_t inx, size_t arg )
{
    assert_true( block_ != nullptr );
    assert_true( inx < block_cnt );
    
    block_size[inx] = arg;
    
    if ( block_alc[inx] < arg )
    {
        constexpr size_t chunk = 32;
        arg = ( arg + chunk - 1 ) & ~( chunk -1 );

        //printf("MatrixOfBlocks::new block %i size %i\n", bb, sz );
        block_alc[inx] = arg;
        free_real(block_[inx]);
        block_[inx] = new_real(arg*arg);
    }
}

//------------------------------------------------------------------------------
size_t MatrixOfBlocks::calculateSize()
{
    size_ = 0;
    for ( size_t ii = 0; ii < block_cnt; ++ii )
        size_ += block_size[ii];
    return size_;
}


//------------------------------------------------------------------------------
real* MatrixOfBlocks::addr(size_t ii, size_t jj) const
{
    size_t bx = 0; //block index on x
    size_t by = 0; //block index on y
    
    while ( ii >= block_size[bx] ) ii -= block_size[bx++];
    while ( jj >= block_size[by] ) jj -= block_size[by++];
    
    if ( bx != by ) // the element is not on a diagonal block
        return nullptr;
        
    if ( block_[bx] == nullptr )  //the block in not allocated
        return nullptr;
    
    assert_true( ii < block_size[bx] );
    assert_true( jj < block_size[bx] );
    assert_true( block_size[bx] <= block_alc[bx] );
    
    return & block_[bx][ii + block_size[bx] * jj];
}

//------------------------------------------------------------------------------
real& MatrixOfBlocks::operator()( size_t x, size_t y )
{
    return * addr( x, y );
}

//------------------------------------------------------------------------------
size_t MatrixOfBlocks::nbElements() const
{
    size_t result=0;
    for ( size_t bb = 0; bb < block_cnt; ++bb )
        result += block_size[bb] * block_size[bb];
    return result;
}

//------------------------------------------------------------------------------
size_t MatrixOfBlocks::maxBlockSize()  const
{
    size_t result = 0;
    for ( size_t bb=0; bb < block_cnt; ++bb )
        if ( result < block_size[bb] ) result = block_size[bb];
    return result;
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::setBlockToZero( const size_t bb )
{
    real* BS = block_[bb];
    if ( BS )
    {
        for ( size_t kk = 0; kk < block_size[bb] * block_size[bb]; ++kk )
            BS[kk] = 0;
    }
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::reset()
{
    for ( size_t ii = 0; ii < block_cnt ; ++ii )
        setBlockToZero( ii );
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::scaleBlock( const size_t bb, const real a )
{
    real* BS = block_[bb];
    if ( BS ) {
        for ( size_t kk = 0; kk < block_size[bb] * block_size[bb]; ++kk )
            BS[kk] = a * BS[kk];
    }
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::scale(const real a)
{
    for ( size_t bb = 0; bb < block_cnt ; ++bb )
        scaleBlock( bb, a );
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::vecMulAdd( const real* X, real* Y )  const
{
    size_t xx = 0;
    for ( size_t bb = 0; bb < block_cnt; ++bb )
    {
        assert_true( block_[bb] );
        blas::xgemv('N', block_size[bb], block_size[bb], 1.0, 
                   block_[bb], block_size[bb], X+xx, 1, 1.0, Y+xx, 1);
        xx += block_size[bb];
    }
}

//------------------------------------------------------------------------------
void MatrixOfBlocks::vecMul( const real* X, real* Y )  const
{
    size_t xx = 0;
    for ( size_t bb = 0; bb < block_cnt; ++bb )
    {
        assert_true( block_[bb] );
        blas::xgemv('N', block_size[bb], block_size[bb], 1.0, 
                   block_[bb], block_size[bb], X+xx, 1, 0.0, Y+xx, 1);
        xx += block_size[bb];
    }
}

//------------------------------------------------------------------------------
real MatrixOfBlocks::norm_inf() const
{
    real mx = 0;
    for ( size_t bb=0; bb < block_cnt; ++bb )
    {
        real* X = block_[bb];
        size_t m = block_size[bb];
        for ( size_t ii = 0; ii < m*m; ++ii )
            mx = std::max(mx, X[ii]);
    }
    return mx;
}


//------------------------------------------------------------------------------
/// printf debug function in sparse mode: i, j : value
void MatrixOfBlocks::printSparse(std::ostream& os) const
{
    size_t offset=0;
    os.precision(8);
    for ( size_t bb=0; bb < block_cnt; ++bb )
    {
        const size_t bsize = block_size[bb];
        for ( size_t ii=0; ii < bsize; ++ii )
            for ( size_t jj=0; jj < bsize; ++jj )
            {
                os << ii+offset << " " << jj+offset << " ";
                os << std::setw(16) << block_[bb][ii + bsize * jj] << std::endl;
            }
        offset += bsize;
    }
}

//------------------------------------------------------------------------------
std::string MatrixOfBlocks::what() const
{
    std::ostringstream msg;
    msg << "mB " << block_cnt;
    return msg.str();
}

