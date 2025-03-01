// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "fiber_site.h"
#include "iowrapper.h"
#include "simul.h"
#include "cymdef.h"


FiberSite::FiberSite(Fiber const* f, real a)
: hFiber(f), hAbs(a)
{
    assert_true(f);
    inter_ = 0;
    segix_ = 0;
#if FIBER_HAS_LATTICE
    hSite = 0;
    hLattice = f->lattice();
#endif
    //reinterpolate();
}


void FiberSite::clear()
{
    hAbs = 0;
    inter_ = 0;
    segix_ = 0;
    hFiber = nullptr;
#if FIBER_HAS_LATTICE
    hSite = 0;
    hLattice = nullptr;
#endif
}

void FiberSite::relocateM()
{
    assert_true(hFiber);
    hAbs = hFiber->abscissaM();
    reinterpolate(hFiber->interpolateEndM());
}


void FiberSite::relocateP()
{
    assert_true(hFiber);
    hAbs = hFiber->abscissaP();
    reinterpolate(hFiber->interpolateEndP());
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
    if ( hAbs - hFiber->abscissaC() < 0 )
        return MINUS_END;
    else
        return PLUS_END;
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


/// this will return a negative value if the abscissa is outside the fiber's ends
real FiberSite::distanceToNearestEnd() const
{
    assert_true(hFiber);
    return std::min(hAbs - hFiber->abscissaM(), hFiber->abscissaP() - hAbs);
}


real FiberSite::abscissaFrom(const FiberEnd ref) const
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


void FiberSite::writeFiberSite(Outputter& out) const
{
    if ( hFiber )
    {
        checkAbscissa();
#if FIBER_HAS_LATTICE
        if ( hLattice )
        {
            Object::writeReference(out, Fiber::LATTICE_TAG, hFiber->identity());
            // in older format, `hAbs` was written here
            out.writeInt32(hSite);
        }
        else
#endif
        {
            Object::writeReference(out, Fiber::TAG, hFiber->identity());
            out.writeFloat(hAbs);
        }
    }
    else
    {
        Object::writeReference(out, Object::NULL_TAG, 0);
    }
}


ObjectID FiberSite::readFiberSite(Inputter& in, Simul& sim)
{
    ObjectID id = 0;
    ObjectTag tag = 0;
    hFiber = sim.readFiberReference(in, tag, id);
    
    if ( hFiber )
    {
        //std::clog << "FiberSite::readFiberSite() " << (char)tag << '\n';
        if ( tag == Fiber::TAG )
        {
            hAbs = in.readFloat();
#if FIBER_HAS_LATTICE
            if ( hLattice ) // set site to closest integral position
                hSite = hLattice->index(hAbs);
#endif
        }
        else if ( tag == Fiber::LATTICE_TAG )
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
            int t = in.readInt32();
            // relying on the lattice_unit being correct at this stage:
            hAbs = ( t + 0.5 ) * hFiber->prop->lattice_unit;
            //throw InvalidIO("Cannot import Digit without fiber's lattice");
#endif
        }
        else
        {
            ///\todo: we could allow binder to refer to any Mecable
            throw InvalidIO("unexpected class in FiberSite");
        }

        // reinterpolate() will be called in updateFiber();
        //checkAbscissa();
    }
    return id;
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
            std::streamsize p = os.precision(3);
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
    if ( hFiber && hFiber->betweenMP(hAbs) )
    {
        // the abscissa of the interpolated point:
        real a = hFiber->abscissaPoint(real(segix_)+inter_);

        constexpr real MAG = 1000;
        const real e = MAG * ( hAbs - a );
        
        //std::clog << "Interpolation " << std::scientific << e << '\n';
        if ( abs_real(e) > 1 )
        {
            Interpolation pi = hFiber->interpolateAbs(hAbs);
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


