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
    pPivot     = nullptr;
    pBlockAlc  = 0;
    pPivotAlc  = 0;
#if EXPERIMENTAL_PRECONDITIONNERS
    pBlockAge  = 0;
#endif
    pBlockSize = 0;
    pBlockType = 0;
    pPos       = nullptr;
    pForce     = nullptr;
    pForceMax  = 0;
    pIndex     = 0;
}


Mecable::Mecable(const Mecable & o) : Mecable()
{
    setNbPoints(o.nPoints);
    copy_real(DIM*nPoints, o.pPos, pPos);
}


Mecable& Mecable::operator =(const Mecable& o)
{
    setNbPoints(o.nPoints);
    copy_real(DIM*nPoints, o.pPos, pPos);
    return *this;
}

//------------------------------------------------------------------------------

/**
Set block size to 'bks' and allocate as necessary to hold 'alc' scalars
 */
void Mecable::blockSize(size_t bks, size_t alc, size_t pivot)
{
    assert_true( bks <= DIM * nPoints );
    pBlockSize = bks;
    
    if ( alc > pBlockAlc )
    {
        free_real(pBlock);
        pBlockAlc = static_cast<unsigned>(chunk_real(alc));
        assert_true( pBlockAlc == chunk_real(alc) );
        // add 4 slots to allow for some SIMD instruction burr:
        pBlock = new_real(pBlockAlc+4);
        //zero_real(pBlockAlc+4, pBlock);
        //std::clog << reference() << " allocateBlock " << bks << " " << pBlockAlc << "\n";
    }
    
    if ( pivot > pPivotAlc )
    {
        delete[] pPivot;
        pPivotAlc = chunk_real(pivot);
        pPivot = new int[pPivotAlc];
    }
}


/**
 allocateMecable(size) ensures that the set can hold `size` points
 returns: size of new memory allocated, or 0
 some extra space is allowed in 3D to allow for AVX overspill
 */
size_t Mecable::allocateMecable(const size_t nbp)
{
    pForce = nullptr;
    if ( pAllocated < nbp )
    {
        size_t all = chunk_real(nbp+(DIM==3));
        // std::clog << "mecable(" << reference() << ") allocates " << all << '\n';
        
        // allocate memory:
        real * mem = new_real(DIM*all);

        // retain existing data:
        if ( pPos )
        {
            assert_true(nPoints < all);
            copy_real(DIM*nPoints, pPos, mem);
            free_real(pPos);
        }
        pPos = mem;
        pAllocated = all-(DIM==3);
        return all;
    }
    
    return 0;
}


void Mecable::release()
{
    free_real(pBlock);
    pBlock = nullptr;
    pBlockSize = 0;
    pBlockAlc = 0;

    delete[] pPivot;
    pPivot = nullptr;
    
    free_real(pPos);
    pPos = nullptr;
    
    pForce = nullptr;
    pAllocated = 0;
    nPoints = 0;
}


//------------------------------------------------------------------------------
#pragma mark - Modifying points

size_t Mecable::addPoint(Vector const& vec)
{
    allocateMecable(nPoints+1);
    size_t i = nPoints++;
    //std::clog << "mecable " << reference() << " point" << i+1 << " = " << vec << "\n";
    vec.store(pPos+DIM*i);
    return i;
}


void Mecable::removePoints(const size_t inx, const size_t nbp)
{
    assert_true( inx + nbp <= nPoints );
    
    nPoints -= nbp;
    
    //move part of the array down, to erase 'nbp' points from index 'inx'
    for ( size_t i = DIM*inx; i < DIM*nPoints; ++i )
        pPos[i] = pPos[i+DIM*nbp];
}


void Mecable::shiftPoints(const size_t inx, const size_t nbp)
{
    allocateMecable(nPoints+nbp);
    
    //move part of the array up, making space for 'nbp' points from index 'inx'
    for ( size_t i = DIM*inx; i < DIM*nPoints; ++i )
        pPos[i+DIM*nbp] = pPos[i];
    
    nPoints += nbp;
}

//------------------------------------------------------------------------------
/**
 shifts array to keep only points within [p, last]
 */
void Mecable::truncateM(const size_t p)
{
    assert_true( p < nPoints - 1 );
    
    size_t np = nPoints - p;
    
    for ( size_t i = 0; i < DIM*np; ++i )
        pPos[i] = pPos[i+DIM*p];
    
    nPoints = np;
}

/**
 erase higher indices of array to keep [0, p]
 */
void Mecable::truncateP(const size_t p)
{
    assert_true( p < nPoints );
    assert_true( p > 0 );
    
    nPoints = p+1;
}

//------------------------------------------------------------------------------

void Mecable::resetPoints()
{
    for ( size_t i = 0; i < DIM*pAllocated; ++i )
        pPos[i] = 0;
}


void Mecable::addNoise(const real mag)
{
    for ( size_t i = 0; i < DIM*nPoints; ++i )
        pPos[i] += mag * RNG.sreal();
}


void Mecable::translate(Vector const& T)
{
    for ( size_t i = 0; i < nPoints; ++i )
        T.add_to(pPos+DIM*i);
}


void Mecable::rotate(Rotation const& T)
{
    for ( size_t i = 0; i < nPoints; ++i)
        ( T.vecmul(pPos+DIM*i) ).store(pPos+DIM*i);
}


//------------------------------------------------------------------------------
#pragma mark - Export/Inport

/** Assuming that ptr[] is rightfully allocated! */
void Mecable::putPoints(real * pts) const
{
    copy_real(DIM*nPoints, pPos, pts);
}


/** Assuming that ptr[] is rightfully allocated! */
void Mecable::getPoints(const real * pts)
{
    copy_real(DIM*nPoints, pts, pPos);
}


void Mecable::setPoints(const real pts[], const size_t nbp)
{
    setNbPoints(nbp);
    copy_real(DIM*nbp, pts, pPos);
}


/**
Copy the coordinates of the points of Object to array `ptr[]`, which has been
allocated to hold `cnt` coordinates. Thus in 3D, `cnt` should be >= 3*Object::nbPoints()
Exactly 3*Object::nbPoints() values are set at most. The Z component is set to zero in 2D mode.
The data is converted to single precision.

@return error code: 0 = no error, 1 = insufficient allocation
*/
int Mecable::putPoints(float ptr[], size_t cnt) const
{
#if ( DIM == 2 )
    size_t sup = std::min(DIM*(size_t)nPoints, cnt);
    for ( size_t i = 0; i < sup; ++sup )
        ptr[i] = pPos[i];
#else
    size_t sup = std::min((size_t)nPoints, cnt/3);
    for ( size_t i = 0; i < sup; ++sup )
    {
        for ( size_t d = 0; d < DIM; ++d )
            ptr[3*i+d] = pPos[3*i+d];
        for ( size_t d = DIM; d < 3; ++d )
            ptr[3*i+d] = 0;
    }
#endif
    return ( cnt < 3*nPoints );
}


Vector Mecable::netForce(const size_t p) const
{
    assert_true( !pForce || nPoints==pForceMax );
    
    if (( pForce != nullptr ) & ( p < pForceMax ))
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
    for ( size_t i = 1; i < nPoints; ++i )
        sum += posP(i);
    return sum / real(nPoints);
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
    
    for ( size_t i = 0; i < nPoints; ++i )
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


bool Mecable::allInside(Space const* spc) const
{
    for ( size_t i = 0; i < nPoints; ++i )
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
    for ( size_t i = 0; i < nPoints ; ++i )
        out.writeFloats(pPos+DIM*i, DIM, '\n');
}


void Mecable::read(Inputter& in, Simul&, ObjectTag)
{
    try
    {
        size_t nb = in.readUInt16();
        if ( allocateMecable(nb) )
            resetPoints();    //we reset the point for a clean start
        nPoints = nb;
#if !REAL_IS_DOUBLE
        in.readFloats(nb, pPos, DIM);
#else
        for ( size_t i = 0; i < nb ; ++i )
            in.readFloats(pPos+DIM*i, DIM);
#endif
    }
    catch( Exception & e )
    {
        clearPoints();
        throw;
    }
}


void Mecable::print(std::ostream& os, real const* ptr) const
{
    os << "new mecable " << reference() << "\n{\n";
    os << " nb_points = " << nPoints << '\n';
    for ( size_t i = 0; i < nPoints ; ++i )
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


size_t Mecable::point_index(std::string const& str) const
{
    const size_t sup = nbPoints();
    if ( str.size() > 5  &&  str.compare(0,5,"point") == 0 )
    {
        errno = 0;
        unsigned long i = strtoul(str.c_str()+5, nullptr, 10);
        if ( errno ) throw InvalidParameter("a point index must be specified, eg. `point1`");
        if ( i < 1 ) throw InvalidParameter("a point index must must be >= 1");
        if ( i > sup ) throw InvalidParameter("point index is out of range");
        return i - 1;
    }
    throw InvalidParameter("expected a point specification eg. `point1'");
    return 0;
}

