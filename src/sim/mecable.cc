// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "mecable.h"
#include "exceptions.h"
#include "blas.h"
#include "iowrapper.h"
#include "organizer.h"
#include "space.h"


//------------------------------------------------------------------------------
/**
clear pointers
 */
Mecable::Mecable()
{
    pAllocated = 0;
    nPoints    = 0;
    pBlock     = nullptr;
    pBlockAlc  = 0;
    pBlockType = 0;
    pPos       = nullptr;
    pForce     = nullptr;
    pIndex     = ~0;
#if ADD_PROJECTION_DIFF
    useProjectionDiff = false;
#endif
}


void Mecable::setNbPoints(index_t n)
{
    if ( n != nPoints )
    {
        allocateMecable(n);
        nPoints = (SIZE_T)n;
        assert_true(nPoints==n);
        // invalidate data that depend on number of points:
        pForce = nullptr;
        pBlockType = 0;
    }
}


Mecable::Mecable(const Mecable & o) : Mecable()
{
    setNbPoints(o.nPoints);
    copy_real(DIM*nPoints, o.pPos, pPos);
}


Mecable& Mecable::operator = (const Mecable& o)
{
    setNbPoints(o.nPoints);
    copy_real(DIM*nPoints, o.pPos, pPos);
    return *this;
}

//------------------------------------------------------------------------------

/**
Set block size to 'bks' and allocate as necessary to hold 'alc' reals, and 'pivot' integers
 */
void Mecable::blockSize(index_t bks, index_t alc, index_t pivot)
{
    assert_true( bks <= DIM * nPoints );
    // add enough to cover 'pivot' integers and a bit more:
    alc += 4 + ( pivot * sizeof(int) ) / sizeof(real);
    
    if ( alc > pBlockAlc )
    {
        free_real(pBlock);
        pBlockAlc = chunk_real(alc);
        assert_true( pBlockAlc == chunk_real(alc) );
        // add 4 slots to allow for some SIMD instruction burr:
        pBlock = new_real(pBlockAlc);
        //zero_real(pBlockAlc, pBlock);
        //std::clog << reference() << " allocateBlock(" << bks << " " << pivot << " " << alc << ")\n";
    }
    // use aligned memory:
    if ( pivot > 0 )
        pPivot = reinterpret_cast<int*>(pBlock+pBlockAlc) - (( 1 + pivot ) & ~1);
    else
        pPivot = nullptr;
}


/**
 allocateMecable(N) ensures that the object can hold `nbp` vertices
 @returns pointer to new memory allocated, or nullptr if no allocation was necessary
 some extra space is allowed in 3D to allow for AVX overspill
 */
real* Mecable::allocateMemory(const index_t nbp, index_t add)
{
    if ( nbp + (DIM==3) > pAllocated )
    {
        // in 3D, we want one extra vector for SIMD burr
        index_t all = chunk_real(nbp+(DIM==3));
        // std::clog << "mecable(" << reference() << ") allocates " << all << '\n';
        
        // allocate memory for vertices + requested extra:
        real * mem = new_real(all*(DIM+add));
        // reset the point for a clean start:
        zero_real(all*(DIM+add), mem);

        // transfer existing data:
        if ( pPos )
        {
            assert_true(nPoints < all);
            // copy vertex coordinates:
            copy_real(nPoints*DIM, pPos, mem);
            // copy additional chunks of data:
            copy_real(nPoints*add, pPos+pAllocated*DIM, mem+all*DIM);
            free_real(pPos);
        }
        pPos = mem;
        pAllocated = all;
        return mem + all * DIM;
    }
    return nullptr;
}


void Mecable::release()
{
    free_real(pBlock);
    pBlock = nullptr;
    pBlockAlc = 0;
    
    if ( pAllocated ) free_real(pPos);
    pPos = nullptr;
    
    pForce = nullptr;
    pAllocated = 0;
    nPoints = 0;
}


//------------------------------------------------------------------------------
#pragma mark - Modifying points

index_t Mecable::addPoint(Vector const& vec)
{
    allocateMecable(nPoints+1);
    index_t i = nPoints++;
    //std::clog << "mecable " << reference() << " point" << i+1 << " = " << vec << "\n";
    vec.store(pPos+DIM*i);
    return i;
}


void Mecable::removePoints(const index_t inx, const index_t nbp)
{
    assert_true( inx + nbp <= nPoints );
    
    nPoints -= nbp;
    
    //move part of the array down, to erase 'nbp' points from index 'inx'
    for ( index_t i = DIM*inx; i < DIM*nPoints; ++i )
        pPos[i] = pPos[i+DIM*nbp];
}


void Mecable::shiftPoints(const index_t inx, const index_t nbp)
{
    allocateMecable(nPoints+nbp);
    
    //move part of the array up, making space for 'nbp' points from index 'inx'
    for ( index_t i = DIM*inx; i < DIM*nPoints; ++i )
        pPos[i+DIM*nbp] = pPos[i];
    
    nPoints += nbp;
}

//------------------------------------------------------------------------------
/**
 shifts array to keep only points within [p, last]
 */
void Mecable::truncateM(const index_t p)
{
    assert_true( p < nPoints - 1 );
    
    index_t np = nPoints - p;
    
    for ( index_t i = 0; i < DIM*np; ++i )
        pPos[i] = pPos[i+DIM*p];
    
    nPoints = np;
}

/**
 erase higher indices of array to keep [0, p]
 */
void Mecable::truncateP(const index_t p)
{
    assert_true( p < nPoints );
    assert_true( p > 0 );
    
    nPoints = p+1;
}

//------------------------------------------------------------------------------

void Mecable::resetPoints()
{
    for ( index_t i = 0; i < DIM*pAllocated; ++i )
        pPos[i] = 0;
}


void Mecable::addNoise(const real mag)
{
    for ( index_t i = 0; i < DIM*nPoints; ++i )
        pPos[i] += mag * RNG.sreal();
}


void Mecable::translate(Vector const& T)
{
    for ( index_t i = 0; i < nPoints; ++i )
        T.add_to(pPos+DIM*i);
}


void Mecable::rotate(Rotation const& T)
{
    for ( index_t i = 0; i < nPoints; ++i)
        T.vecmul(pPos+DIM*i).store(pPos+DIM*i);
}


//------------------------------------------------------------------------------
#pragma mark - Export/Inport

/** Assuming that pts[] is rightfully allocated! */
void Mecable::putPoints(real * pts) const
{
    copy_real(DIM*nPoints, pPos, pts);
}


/** Assuming that pts[] is rightfully allocated! */
void Mecable::getPoints(const real * pts)
{
    copy_real(DIM*nPoints, pts, pPos);
}


void Mecable::setPoints(const real pts[], const index_t nbp)
{
    setNbPoints(nbp);
    copy_real(DIM*nbp, pts, pPos);
}


/**
Copy vertex coordinates to given array, converting to single precision.
*/
void Mecable::putPoints(float ptr[], index_t sup) const
{
    sup = std::min(DIM * nbPoints(), sup);
    for ( index_t i = 0; i < sup; ++i )
        ptr[i] = pPos[i];
}


Vector Mecable::netForce(const index_t p) const
{
    if ( pForce )
        return Vector(pForce+DIM*p);
    else
        return Vector(0,0,0);
}

//------------------------------------------------------------------------------
/**
 Returns the center of gravity of all points
 */
Vector Mecable::position() const
{
    Vector sum = posP(0);
    for ( index_t i = 1; i < nPoints; ++i )
        sum += posP(i);
    return sum / real(nPoints);
}


Vector Mecable::interpolatePoints(index_t ref, real const coef[], index_t rank) const
{
    assert_true( rank > 0 );
    assert_true( ref < nPoints );
    index_t end = nbPoints() - ref;
    if ( rank < end )
        end = rank;
    Vector res = coef[0] * posP(ref);
    for ( index_t i = 1; i < end; ++i )
        res += coef[i] * posP(ref+i);
    return res;
}


/**
 Calculate first and second moment of vertex coordinates:
 - avg = sum( P ) / num_points
 - dev = sum( P .* P ) / num_points - square( avg );
 .
 */
void Mecable::calculateMomentum(Vector& avg, Vector& dev)
{
    avg.reset();
    dev.reset();
    
    for ( index_t i = 0; i < nPoints; ++i )
    {
        Vector x = posPoint(i);
        avg += x;
        dev += x.e_squared();
    }
    
    if ( nPoints > 1 )
    {
        avg /= nPoints;
        dev /= nPoints;
    }
    
    dev -= avg.e_squared();
}


void Mecable::foldPosition(Modulo const* m)
{
    Vector off = m->offset(position());
    if ( off.is_not_zero() )
        translate(-off);
}


bool Mecable::allPointsInside(Space const* spc) const
{
    for ( index_t i = 0; i < nPoints; ++i )
    {
        if ( spc->outside(posP(i)) )
            return false;
    }
    return true;
}


//------------------------------------------------------------------------------
#pragma mark - Read/write


void Mecable::write(Outputter& out) const
{
    out.writeUInt16(nPoints);
    for ( index_t i = 0; i < nPoints ; ++i )
        out.writeFloats(pPos+DIM*i, DIM, '\n');
}


void Mecable::read(Inputter& in, Simul&, ObjectTag)
{
    index_t nb = in.readUInt16();
    setNbPoints(nb);
    in.readFloats(nb, pPos, DIM);
}


void Mecable::print(std::ostream& os, real const* ptr) const
{
    os << "new mecable " << reference() << "\n{\n";
    os << " nb_points = " << nPoints << '\n';
    for ( index_t i = 0; i < nPoints ; ++i )
    {
        os << " point" << i+1 << " = " << Vector(ptr+DIM*i) << '\n';
    }
    os << "}\n";
}


std::ostream& operator << (std::ostream& os, Mecable const& arg)
{
    arg.print(os, arg.addrPoints());
    return os;
}


index_t Mecable::point_index(std::string const& str) const
{
    const index_t sup = nbPoints();
    if ( str.size() > 5  &&  str.compare(0,5,"point") == 0 )
    {
        unsigned long i = 0;
        try {
             i = std::stoul(str.substr(5), nullptr, 10);
        }
        catch ( ... ) {
            throw InvalidParameter("a point index must be specified, eg. `point1`");
        }
        if ( i < 1 ) throw InvalidParameter("a point index must must be >= 1");
        if ( i > sup ) throw InvalidParameter("point index is out of range");
        return (index_t)(i - 1);
    }
    throw InvalidParameter("expected a point specification eg. `point1'");
    return 0;
}


int Mecable::invalid() const
{
    for ( index_t i = 0; i < DIM * nPoints ; ++i )
        if ( pPos[i] != pPos[i] )
            return 1;
    return 0;
}
