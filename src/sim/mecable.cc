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
void Mecable::clearMecable()
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
    pIndex     = 0;
}


Mecable::Mecable(const Mecable & o)
{
    clearMecable();
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
Set block size and allocate if necessary to hold 'alc' scalars
 */
void Mecable::blockSize(size_t bks, size_t block, size_t pivot)
{
    assert_true( bks <= DIM * nPoints );
    pBlockSize = bks;
    
    if ( block > pBlockAlc )
    {
        free_real(pBlock);
        pBlockAlc = chunk_real(block);
        //std::clog <<reference()<<" allocateBlock " << pBlockAlc << "\n";
        pBlock = new_real(pBlockAlc);
        //zero_real(pBlockAlc, pBlock);
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
    delete[] pPivot;
    pPivot = nullptr;
    
    pBlockAlc  = 0;
    pBlockSize = 0;
    
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
    size_t p = nPoints++;
    //std::clog << "mecable " << reference() << " point" << p+1 << " = " << vec << "\n";
    vec.store(pPos+DIM*p);
    return p;
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
    
    for ( size_t ii = 0; ii < DIM*np; ++ii )
        pPos[ii] = pPos[ii+DIM*p];
    
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
    if ( pPos )
    {
        for ( size_t p = 0; p < DIM*pAllocated; ++p )
            pPos[p] = 0;
    }
}


void Mecable::addNoise(const real amount)
{
    for ( size_t p = 0; p < DIM*nPoints; ++p )
        pPos[p] += amount * RNG.sreal();
}


void Mecable::translate(Vector const& T)
{
    for ( size_t p = 0; p < nPoints; ++p )
        T.add_to(pPos+DIM*p);
}


void Mecable::rotate(Rotation const& T)
{
    for ( size_t p = 0; p < nPoints; ++p)
        ( T.vecmul(pPos+DIM*p) ).store(pPos+DIM*p);
}


//------------------------------------------------------------------------------
#pragma mark - Export/Inport

void Mecable::putPoints(real * ptr) const
{
    copy_real(DIM*nPoints, pPos, ptr);
}


void Mecable::getPoints(const real * ptr)
{
    copy_real(DIM*nPoints, ptr, pPos);
}


void Mecable::setPoints(const real pts[], const size_t nbp)
{
    setNbPoints(nbp);
    copy_real(DIM*nbp, pts, pPos);
}


/**
Copy the coordinates of the points of Object to array `ptr[]`, previously
allocated to hold `ptr_n` coordinates. Thus in 3D, `ptr_n` should be >= 3*Object::nbPoints()
Exactly 3*Object::nbPoints() values are set at most. The Z component is set to zero in 2D mode.

@return error code: 0 = no error, 1 = insufficient allocation
*/
int Mecable::putPoints(float ptr[], size_t ptr_n) const
{
#if ( DIM == 2 )
    size_t sup = std::min(DIM*nPoints, ptr_n);
    for ( size_t i = 0; i < sup; ++sup )
        ptr[i] = pPos[i];
#else
    size_t sup = std::min(nPoints, ptr_n/3);
    for ( size_t i = 0; i < sup; ++sup )
    {
        for ( size_t d = 0; d < DIM; ++d )
            ptr[3*i+d] = pPos[3*i+d];
        for ( size_t d = DIM; d < 3; ++d )
            ptr[3*i+d] = 0;
    }
#endif
    return ( ptr_n < 3*nPoints );
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
    for ( size_t p = 1; p < nPoints; ++p )
        sum += posP(p);
    return sum / real(nPoints);
}


/**
 Calculate first and second moment of point distribution:
 - avg = sum( P ) / nb_points
 - sec = sum( P * P );
 .
 if 'sub = true', the average is substracted from 'sec'
 */
void Mecable::calculateMomentum(Vector& avg, Vector& sec, bool sub)
{
    avg.reset();
    sec.reset();
    
    // calculate first and second moments:
    for ( size_t p = 0; p < nPoints; ++p )
    {
        avg += posP(p);
        sec += posP(p).e_squared();
        /*
         real const* pp = pPos + DIM*p;
         avg.XX += pp[0];
         sec.XX += pp[0] * pp[0];
         #if ( DIM > 1 )
         avg.YY += pp[1];
         sec.YY += pp[1] * pp[1];
         #endif
         #if ( DIM > 2 )
         avg.ZZ += pp[2];
         sec.ZZ += pp[2] * pp[2];
         #endif
         */
    }
    
    if ( nPoints > 1 )
    {
        avg /= nPoints;
        sec /= nPoints;
    }
    
    if ( sub )
        sec -= avg.e_squared();
}


void Mecable::foldPosition(Modulo const* s)
{
    Vector off = s->offset(position());
    if ( off.is_not_zero() )
        translate(-off);
}


bool Mecable::allInside(Space const* spc) const
{
    for ( size_t ii = 0; ii < nPoints; ++ii )
    {
        if ( spc->outside(posP(ii)) )
            return false;
    }
    return true;
}

//------------------------------------------------------------------------------
#pragma mark - Read/write


void Mecable::write(Outputter& out) const
{
    out.writeUInt16(nPoints);
    for ( size_t p = 0; p < nPoints ; ++p )
        out.writeFloats(pPos+DIM*p, DIM, '\n');
}


void Mecable::read(Inputter& in, Simul&, ObjectTag)
{
    try
    {
        size_t nb = in.readUInt16();
        allocateMecable(nb);
        //we reset the point for a clean start:
        resetPoints();
        nPoints = nb;
#if ( 1 )
        for ( size_t p = 0; p < nb ; ++p )
            in.readFloats(pPos+DIM*p, DIM);
#else
        in.readFloats(pPos, nb, DIM);
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

