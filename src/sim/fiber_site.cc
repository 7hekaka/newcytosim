// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "fiber_site.h"
#include "iowrapper.h"
#include "simul.h"
#include "cymdef.h"


FiberSite::FiberSite(Fiber* f, real a)
: hFiber(f), hAbs(a)
{
    assert_true(f);
#if FIBER_HAS_LATTICE
    hLattice = nullptr;
    hSite = 0;
#endif
    inter_ = 0;
    segix_ = 0;
    //reinterpolate();
}


void FiberSite::clear()
{
    hAbs = 0;
    inter_ = 0;
    segix_ = 0;
    hFiber = nullptr;
#if FIBER_HAS_LATTICE
    hLattice = nullptr;
#endif
}

void FiberSite::relocateM()
{
    assert_true(hFiber);
    hAbs = hFiber->abscissaM();
    reinterpolate(hFiber->interpolateEndM());
#if FIBER_HAS_LATTICE
    assert_true(!hLattice);
#endif
}


void FiberSite::relocateP()
{
    assert_true(hFiber);
    hAbs = hFiber->abscissaP();
    reinterpolate(hFiber->interpolateEndP());
#if FIBER_HAS_LATTICE
    assert_true(!hLattice);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


#if FIBER_HAS_FAMILY
Vector FiberSite::outerPos() const
{
    if ( hFiber->family_ != hFiber )
    {
        real a = hAbs - hFiber->abscissaM();
#if DIM == 3
        // using the two flanking protofilaments to set the outer direction
        assert_true( hFiber->brother_ != hFiber->sister_ );
        Vector b = hFiber->brother_->posM(a);
        Vector s = hFiber->sister_->posM(a);
        return hFiber->posM(a) + cross(b-s, hFiber->dirM(a));
#else
        // using a central backbone to set the outer direction
        Vector p = hFiber->posM(a);
        Vector c = hFiber->family_->posM(a); //centerline
        return p + 0.5 * ( p - c );
#endif
        //return 2 * pos - hFiber->family_->posM(ab);
    }
    return pos();
}
#endif


FiberEnd FiberSite::nearestEnd() const
{
    assert_true(hFiber);
    if ( hAbs > hFiber->abscissaC() )
        return PLUS_END;
    else
        return MINUS_END;
}


real FiberSite::distanceToEnd(FiberEnd end) const
{
    assert_true(hFiber);
    if ( end == PLUS_END )
        return hFiber->abscissaP() - hAbs;
    else
    {
        assert_true(end == MINUS_END);
        return hAbs - hFiber->abscissaM();
    }
}


real  FiberSite::abscissaFrom(const FiberEnd ref) const
{
    assert_true(hFiber);
    switch( ref )
    {
        case MINUS_END:  return hAbs - hFiber->abscissaM();
        case PLUS_END:   return hFiber->abscissaP() - hAbs;
        case ORIGIN:     return hAbs;
        case CENTER:     return hAbs - hFiber->abscissaC();
        default:         ABORT_NOW("invalid argument value");
    }
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark - I/O


void FiberSite::write(Outputter& out) const
{
    if ( hFiber )
    {
        checkAbscissa();
#if FIBER_HAS_LATTICE
        if ( hLattice )
        {
            Object::writeReference(out, Fiber::TAG_LATTICE, hFiber->identity());
            // in older format, `hAbs` was written here
            out.writeInt32(hSite);
        }
        else
#endif
        {
#if !NEW_COMPACT_STORAGE
            // normal way
            Object::writeReference(out, Fiber::TAG, hFiber->identity());
            out.writeFloat(hAbs);
#else
            // compact format created on 23/06/2021
            // sacrificing precision to save a bit of space (29/06/2021)
            Object::writeReference(out, Fiber::TAG_ALT, hFiber->identity());
            // calculate relative position on fiber, which should be in [0, 1]:
            real x = ( hAbs - hFiber->abscissaM() ) / ( hFiber->length() );
            out.writeFixed(x); // 2 bytes
#endif
        }
    }
    else
    {
        Object::writeNullReference(out);
    }
}


void FiberSite::read(Inputter& in, Simul& sim)
{
    ObjectTag tag = 0;
    hFiber = sim.readFiberReference(in, tag);
    
    if ( hFiber )
    {
        //std::clog << "FiberSite::read() " << (char)tag << '\n';
        if ( tag == Fiber::TAG )
        {
            hAbs = in.readFloat();
#if FIBER_HAS_LATTICE
            if ( hLattice ) // set site to closest integral position
                hSite = hLattice->index(hAbs);
#endif
        }
        else if ( tag == Fiber::TAG_ALT )
        {
            real x = in.readFixed();
            hAbs = x * hFiber->length() + hFiber->abscissaM();
#if FIBER_HAS_LATTICE
            if ( hLattice ) // set site to closest integral position
                hSite = hLattice->index(hAbs);
#endif
        }
        else if ( tag == Fiber::TAG_LATTICE )
        {
#if BACKWARD_COMPATIBILITY < 49
            if ( in.formatID() < 49 )
                hAbs = in.readFloat();
#endif
#if FIBER_HAS_LATTICE
            hSite = in.readInt32();
            hLattice = hFiber->lattice();
            // put in the middle of the site:
            // the abscissa will be adjusted in Fiber::resetLattice()
            hAbs = ( hSite + 0.5 ) * hLattice->unit();
#else
            int32_t t = in.readInt32();
            // relying on the lattice_unit being correct at this stage:
            hAbs = ( t + 0.5 ) * hFiber->prop->lattice_unit;
            //throw InvalidIO("Cannot import Digit without fiber's lattice");
#endif
        }
        else
        {
            ///\todo: we should allow binder to refer to any Mecable
            throw InvalidIO("unexpected class in FiberSite");
        }

        // reinterpolate() will be called in updateFiber();
        //checkAbscissa();
    }
}

void FiberSite::print(std::ostream& os) const
{
    if ( fiber() )
    {
#if FIBER_HAS_LATTICE
        if ( hLattice )
            os << "[" << fiber()->reference() << " " << hSite << "]";
        else
#endif
        {
            std::streamsize p = os.precision();
            os.precision(3);
            os << "(" << fiber()->reference() << " " << std::fixed << abscissa() << ")";
            os.precision(p);
        }
    } else
        os << "[null]";
}

std::ostream& operator << (std::ostream& os, FiberSite const& arg)
{
    arg.print(os);
    return os;
}


//------------------------------------------------------------------------------
#pragma mark -


int FiberSite::checkAbscissa() const
{
    assert_true(hFiber);
    
    real a = hFiber->abscissaM() - hAbs;
    if ( a > real(1e-3) )
    {
        std::cerr << "FiberSite:abscissa < fiber:abscissa(MINUS_END) by " << a << '\n';
        return 2;
    }
    
    real b = hAbs - hFiber->abscissaP();
    if ( b > real(1e-3) )
    {
        std::cerr << "FiberSite:abscissa > fiber:abscissa(PLUS_END) by " << b << '\n';
        return 1;
    }
    return 0;
}


int FiberSite::bad() const
{
    if ( hFiber->betweenMP(hAbs) )
    {
        // the abscissa of the interpolated point:
        real a = hFiber->abscissaPoint(real(segix_)+inter_);

        constexpr real MAG = 1000;
        const real e = MAG * ( hAbs - a );
        
        //std::clog << "Interpolation " << std::scientific << e << '\n';
        if ( abs_real(e) > 1 )
        {
            Interpolation pi = hFiber->interpolate(hAbs);
            real b = hFiber->abscissaPoint(pi.point1() + pi.coef1());
            std::cerr << "FiberSite::Interpolation error " << e << " nm in abscissa:\n";
            std::cerr << "    binder       " << MAG * hAbs << "\n";
            std::cerr << "    interpolated " << MAG * a << "\n";
            std::cerr << "    updated      " << MAG * b << "\n";
            return 8;
        }
    }
    return 0;
}


