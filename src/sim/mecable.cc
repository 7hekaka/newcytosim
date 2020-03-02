// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "mecable.h"
#include "exceptions.h"
#include "cblas.h"
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
    pBlockUse  = false;
    pBlockSize = 0;
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
 
 The first time, we allocate exactly what is demanded, but if new allocation 
 is required again, we allocate with some margin, because it means that this
 object is probably growing.

 */
void Mecable::allocateBlock()
{
    pBlockSize = DIM * nPoints;
    
    if ( pBlockSize > pBlockAlc )
    {
        free_real(pBlock);
        delete[] pPivot;
        size_t bum = chunk_real(pBlockSize);
        //std::clog << "Mecable("<<reference()<<")::allocateBlock " << bum << "\n";
   
        pBlock = new_real(bum*bum);
        pPivot = new int[bum];
        pBlockAlc = bum;
        
        //zero_real(bum*bum, pBlock);
    }
}


/**
 allocateMecable(size) ensures that the set can hold `size` points
 it returns the size if new memory was allocated
 */
size_t Mecable::allocateMecable(const size_t nbp)
{
    pForce = nullptr;
    if ( pAllocated < nbp )
    {
        size_t all = chunk_real(nbp);
        // std::clog << "mecable(" << reference() << ") allocates " << all << '\n';
        
        // allocate memory:
        real * mem = new_real(DIM*all);

        // retain existing data:
        if ( pPos )
        {
            copy_real(DIM*nPoints, pPos, mem);
            free_real(pPos);
        }
        pPos = mem;
        pAllocated = all;
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

void Mecable::copyPoints(const size_t nbp, const real pts[])
{
    setNbPoints(nbp);
    copy_real(DIM*nbp, pts, pPos);
}


void Mecable::putPoints(real * ptr) const
{
    copy_real(DIM*nPoints, pPos, ptr);
}


void Mecable::getPoints(const real * ptr)
{
    copy_real(DIM*nPoints, ptr, pPos);
}


Vector Mecable::netForce(const size_t p) const
{
    assert_true( !pForce || nPoints==pForceMax );
    
    if ( pForce && p < pForceMax )
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
        out.writeFloatVector(pPos+DIM*p, DIM, '\n');
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
            in.readFloatVector(pPos+DIM*p, DIM);
#else
        in.readFloatVector(pPos, nb, DIM);
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


std::ostream& operator << (std::ostream& os, Mecable const& obj)
{
    obj.print(os, obj.data());
    return os;
}


size_t Mecable::point_index(std::string const& str, size_t X)
{
    if ( str.size() > 5  &&  str.compare(0,5,"point") == 0 )
    {
        errno = 0;
        unsigned long i = strtoul(str.c_str()+5, nullptr, 10);
        if ( errno ) throw InvalidParameter("a point index must be specified, eg. `point1`");
        if ( i < 1 ) throw InvalidParameter("a point index must must be >= 1");
        if ( i > X ) throw InvalidParameter("point index is out of range");
        return i - 1;
    }
    throw InvalidParameter("expected a point specification eg. `point1'");
    return 0;
}

