// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_site.h"
#include "iowrapper.h"
#include "simul.h"
#include "sim.h"


/// specifies the position occupied within the Lattice site
/**
 `SUBSITE_POS` should be in [0, 1]:
 - with `0.0`, the attachment position is at the start of the site
 - with `1.0`, the attachment position is at the end of the site
 - with `0.5`, the attachment is exactly midway
 */
const real SUBSITE_POS = 0.5;


FiberSite::FiberSite(Fiber* f, real a)
: fbFiber(f), fbAbs(a)
{
    assert_true(f);
#if FIBER_HAS_LATTICE
    fbLattice = nullptr;
    fbSite = 0;
#endif
    inter = f->interpolate(a);
}


#if FIBER_HAS_LATTICE
void FiberSite::engageLattice()
{
    assert_true( fbFiber );
    fbLattice = &fbFiber->lattice();
        
    if ( !fbLattice->ready() )
        throw InvalidParameter("Fiber:lattice was not initialized");
    
    // attach to site closest from given abscissa:
    fbSite = fbLattice->index(fbAbs);
    // adjust abscissa
    fbAbs = fbLattice->abscissa(fbSite+SUBSITE_POS);
}


void FiberSite::hop(site_t s)
{
    dec();
    fbSite = s;
    inc();
    fbAbs = fbLattice->abscissa(fbSite+SUBSITE_POS);
    assert_true(fiber()->betweenMP(fbAbs));
    update();
}
#endif


void FiberSite::relocateM()
{
    assert_true(fbFiber);
    fbAbs = fbFiber->abscissaM();
    inter = fbFiber->interpolateEndM();
#if FIBER_HAS_LATTICE
    assert_true(!fbLattice);
#endif
}


void FiberSite::relocateP()
{
    assert_true(fbFiber);
    fbAbs = fbFiber->abscissaP();
    inter = fbFiber->interpolateEndP();
#if FIBER_HAS_LATTICE
    assert_true(!fbLattice);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


FiberEnd FiberSite::nearestEnd() const
{
    assert_true(fbFiber);
    if ( fbAbs > fbFiber->abscissaC() )
        return PLUS_END;
    else
        return MINUS_END;
}


real FiberSite::distanceToEnd(FiberEnd end) const
{
    assert_true(fbFiber);
    if ( end == PLUS_END )
        return fbFiber->abscissaP() - fbAbs;
    else
    {
        assert_true(end == MINUS_END);
        return fbAbs - fbFiber->abscissaM();
    }
}


real  FiberSite::abscissaFrom(const FiberEnd ref) const
{
    assert_true(fbFiber);
    switch( ref )
    {
        case MINUS_END:  return fbAbs - fbFiber->abscissaM();
        case PLUS_END:   return fbFiber->abscissaP() - fbAbs;
        case ORIGIN:     return fbAbs;
        case CENTER:     return fbAbs - fbFiber->abscissaC();
        default:         ABORT_NOW("invalid argument value");
    }
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark - I/O

void FiberSite::write(Outputter& out) const
{
    out.writeSoftSpace();
    if ( fbFiber )
    {
        checkAbscissa();
#if FIBER_HAS_LATTICE
        if ( fbLattice )
        {
            fbFiber->writeReference(out, Fiber::TAG_LATTICE);
            // in older format, `fbAbs` was written here
            out.writeInt32(fbSite);
        }
        else
#endif
        {
            fbFiber->writeReference(out);
            out.writeFloat(fbAbs);
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
    fbFiber = static_cast<Fiber*>(w);

    if ( w )
    {
        //std::clog << "FiberSite::read() " << (char)tag << std::endl;
        if ( tag == Fiber::TAG )
        {
            fbAbs  = in.readFloat();
#if FIBER_HAS_LATTICE
            // set mSite to closest integral position
            if ( fbLattice )
                fbSite = fbLattice->index(fbAbs);
#endif
        }
        else if ( tag == Fiber::TAG_LATTICE )
        {
#ifdef BACKWARD_COMPATIBILITY
            if ( in.formatID() < 49 )
                fbAbs = in.readFloat();
#endif
#if FIBER_HAS_LATTICE
            fbSite = in.readUInt32();
            fbLattice = &fbFiber->lattice();
            fbAbs = fbLattice->abscissa(fbSite);
#else
            in.readUInt32();
            throw InvalidIO("Cannot import Digit without fiber's lattice");
#endif
        }
#ifdef BACKWARD_COMPATIBILITY
        else if ( tag == 'm' )
        {
            fbAbs  = in.readFloat();
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
        if ( fbLattice )
            os << "(" << fiber()->reference() << " site " << fbSite << ")";
        else
#endif
            os << "(" << fiber()->reference() << " abs " << abscissa() << ")";
    } else
        os << "(null)";
}

std::ostream& operator << (std::ostream& os, FiberSite const& obj)
{
    obj.print(os);
    return os;
}


//------------------------------------------------------------------------------
#pragma mark -


void FiberSite::checkAbscissa() const
{
    assert_true(fbFiber);
    
    real a = fbFiber->abscissaM() - fbAbs;
    if ( a > 1e-3 )
        Cytosim::warn << "FiberSite:abscissa < fiber:abscissa(MINUS_END) :  " << a << std::endl;
    
    real b = fbAbs - fbFiber->abscissaP();
    if ( b > 1e-3 )
        Cytosim::warn << "FiberSite:abscissa > fiber:abscissa(PLUS_END)  :  " << b << std::endl;
}


int FiberSite::bad() const
{
    if ( fbFiber != inter.mecable() )
    {
        std::cerr << "Interpolation mismatch " << fbFiber << " " << inter.mecable() << std::endl;
        return 7;
    }
    
    if ( fbFiber->betweenMP(fbAbs) )
    {
        const real e = fbAbs - abscissaInter();
        
        //std::clog << "Interpolation " << std::scientific << e << std::endl;
        if ( fabs(e) > 1e-5 )
        {
            std::cerr << "Interpolation error is " << std::scientific << e << "\n";
            std::cerr << " abscissa:\n";
            std::cerr << "    binder       " << fbAbs << "\n";
            std::cerr << "    interpolated " << abscissaInter() << "\n";
            Interpolation pi = fbFiber->interpolate(fbAbs);
            std::cerr << "    updated      " << fbFiber->abscissaPoint(pi.point1()+pi.coef1()) << "\n";
            return 8;
        }
    }
    return 0;
}


