// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "fiber_site.h"
#include "iowrapper.h"
#include "simul.h"
#include "sim.h"


FiberSite::FiberSite(Fiber* f, real a)
: hFiber(f), hAbs(a)
{
    assert_true(f);
#if FIBER_HAS_LATTICE
    hLattice = nullptr;
    hSite = 0;
#endif
    hTerp = f->interpolate(a);
}


void FiberSite::relocateM()
{
    assert_true(hFiber);
    hAbs = hFiber->abscissaM();
    hTerp = hFiber->interpolateEndM();
#if FIBER_HAS_LATTICE
    assert_true(!hLattice);
#endif
}


void FiberSite::relocateP()
{
    assert_true(hFiber);
    hAbs = hFiber->abscissaP();
    hTerp = hFiber->interpolateEndP();
#if FIBER_HAS_LATTICE
    assert_true(!hLattice);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


#if FIBER_HAS_FAMILY
Vector FiberSite::outerPos() const
{
    if ( hFiber->family_ )
    {
        real ab = hAbs - hFiber->abscissaM();
        Vector pos = hFiber->posM(ab);
        Vector cen = hFiber->family_->posM(ab);
        return pos + 0.5 * ( pos - cen );
        //return 2 * pos - hFiber->family_->posM(ab);
    }
    return hTerp.pos();
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
    out.writeSoftSpace();
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
            Object::writeReference(out, hFiber);
            out.writeFloat(hAbs);
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
    Object * w = sim.readReference(in, tag);
    hFiber = static_cast<Fiber*>(w);

    if ( w )
    {
        //std::clog << "FiberSite::read() " << (char)tag << '\n';
        if ( tag == Fiber::TAG )
        {
            hAbs = in.readFloat();
#if FIBER_HAS_LATTICE
            // set mSite to closest integral position
            if ( hLattice )
                hSite = hLattice->index(hAbs);
#endif
        }
        else if ( tag == Fiber::TAG_LATTICE )
        {
#ifdef BACKWARD_COMPATIBILITY
            if ( in.formatID() < 49 )
                hAbs = in.readFloat();
#endif
#if FIBER_HAS_LATTICE
            hSite = in.readInt32();
            hLattice = &hFiber->lattice();
            // put in the middle of the site:
            // the abscissa will be adjusted in Fiber::resetLattice()
            hAbs = ( hSite + 0.5 ) * hLattice->unit();
#else
            lati_t t = in.readInt32();
            // relying on the lattice_unit being correct at this stage:
            hAbs = ( t + 0.5 ) * hFiber->prop->lattice_unit;
            //throw InvalidIO("Cannot import Digit without fiber's lattice");
#endif
        }
#ifdef BACKWARD_COMPATIBILITY
        else if ( tag == 'm' )
        {
            hAbs = in.readFloat();
        } 
#endif
        else
        {
            ///\todo: we should allow binder to refer to any Mecable
            throw InvalidIO("unexpected class in FiberSite");
        }

        update();
        checkAbscissa();
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
        os << "(null)";
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
    if ( a > 1e-3 )
    {
        std::cerr << "FiberSite:abscissa < fiber:abscissa(MINUS_END) by " << a << '\n';
        return 2;
    }
    
    real b = hAbs - hFiber->abscissaP();
    if ( b > 1e-3 )
    {
        std::cerr << "FiberSite:abscissa > fiber:abscissa(PLUS_END) by " << b << '\n';
        return 1;
    }
    return 0;
}


int FiberSite::bad() const
{
    if ( hFiber != hTerp.mecable() )
    {
        std::cerr << "Interpolation mismatch " << hFiber << " " << hTerp.mecable() << '\n';
        return 7;
    }
    
    if ( hFiber->betweenMP(hAbs) )
    {
        const real e = hAbs - abscissaInterp();
        
        //std::clog << "Interpolation " << std::scientific << e << '\n';
        if ( abs_real(e) > 1e-3 )
        {
            std::cerr << "FiberSite::Interpolation error " << std::scientific << e << "\n";
            std::cerr << " abscissa:\n";
            std::cerr << "    binder       " << hAbs << "\n";
            std::cerr << "    interpolated " << abscissaInterp() << "\n";
            Interpolation pi = hFiber->interpolate(hAbs);
            std::cerr << "    updated      " << hFiber->abscissaPoint(pi.point1()+pi.coef1()) << "\n";
            return 8;
        }
    }
    return 0;
}


