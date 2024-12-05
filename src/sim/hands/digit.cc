// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "digit.h"
#include "digit_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "glossary.h"
#include "lattice.h"
#include "simul_part.h"
#include "hand_monitor.h"


Digit::Digit(DigitProp const* p, HandMonitor* h)
: Hand(p,h)
{
}


/**
 Digit::attachmentAllowed modifies its argument 'sit' in a meaningful way:
 It adjusts in particular the 'site' index as a function of abscissa,
 */
bool Digit::attachmentAllowed(FiberSite& sit) const
{
    if ( Hand::attachmentAllowed(sit) )
    {
#if FIBER_HAS_LATTICE
        FiberLattice const* lat = sit.lattice();
        
        if ( !lat->ready() )
            throw InvalidParameter("a lattice was not defined for `"+sit.fiber()->prop->name()+"'");
        
        // index to site containing given abscissa:
        lati_t s = lat->index(sit.abscissa());
        
        if ( s < lat->entry() )
        {
            if ( prop()->bind_also_end & MINUS_END )
                s = lat->entry();
            else
                return false;
        }
        else if ( lat->fence() < s )
        {
            if ( prop()->bind_also_end & PLUS_END )
                s = lat->fence();
            else
                return false;
        }
        
        if ( valLattice(lat, s) )
            return false;
        
        // adjust to match selected lattice site:
        sit.engageLattice(s, prop()->site_shift);
#endif
        return true;
    }
    return false;
}


/**
 Digit::attachmentAllowed() should have been called before,
 such that the `sit` already points to a valid lattice site
 */
void Digit::attach(FiberSite const& sit)
{
#if FIBER_HAS_LATTICE
    hSite = sit.site();
    hLattice = sit.lattice();
    //std::clog << "offset " << sit.abscissa() - hLattice->unit()*hSite - prop()->site_shift << "\n";
#endif
    Hand::attach(sit);
}


void Digit::detach()
{
    Hand::detach();
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Try to move `n` sites towards the plus end, stopping before any already occupied site.
 */
void Digit::crawlP(const int n)
{
    lati_t s = site();
    lati_t e = s + n;
    
    while ( s < e )
    {
        if ( !valLattice(s+1) )
            ++s;
        else
            break;
    }
    if ( s != site() )
        hopLattice(s);
}


/**
 Try to move `n` sites towards the minus end, stopping before any already occupied site.
 */
void Digit::crawlM(const int n)
{
    lati_t s = site();
    lati_t e = s - n;
    
    while ( s > e )
    {
        if ( !valLattice(s-1) )
            --s;
        else
            break;
    }
    if ( s != site() )
        hopLattice(s);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 The Digit normally does not move by itself
 */
void Digit::handleDisassemblyM()
{
    assert_true( attached() );
    
    if ( RNG.test(prop()->hold_shrinking_end[1]) )
    {
        // do not detach:
        jumpToEndM();
        if ( site() < lattice()->entry() )
            detach();
    }
    else
        detach();
}


/**
 The Digit normally does not move by itself
 */
void Digit::handleDisassemblyP()
{
    assert_true( attached() );
    
    if ( RNG.test(prop()->hold_shrinking_end[0]) )
    {
        // do not detach:
        jumpToEndP();
        if ( site() > lattice()->fence() )
            detach();
    }
    else
        detach();
}


/**
tests detachment
 */
void Digit::stepUnloaded()
{
    assert_true( attached() );
}


/**
(see @ref Stochastic)
 */
void Digit::stepLoaded(Vector const& force)
{
    assert_true( attached() );
}

//------------------------------------------------------------------------------
#pragma mark -

void Digit::promote()
{
#if FIBER_HAS_LATTICE
    if ( hFiber && !hLattice )
    {
        FiberSite sit(*this);
        // attach Hand to Lattice if possible, and detach otherwise
        if ( attachmentAllowed(sit) )
        {
            hLattice = sit.lattice();
            hSite = sit.site();
            static_cast<Digit*>(this)->incLattice();
        }
        else
            detachHand();
    }
#endif
}


std::ostream& operator << (std::ostream& os, Digit const& arg)
{
    os << arg.property()->name() << "(" << arg.fiber()->reference() << ", " << arg.abscissa();
    os << ", " << arg.site() << ")";
    return os;
}


#if FIBER_HAS_LATTICE
bool Fiber::resetLattice(bool lazy)
{
    if ( !fLattice.data() )
    {
        if ( lazy )
            return false;
        fLattice.setRange(abscissaM(), abscissaP());
    }

    fLattice.clear();
    real uni = fLattice.unit();
    
    for ( Hand * h = fHands.front(); h; h = h->next() )
    {
        if ( h->lattice() == lattice() )
        {
            Digit* i = static_cast<Digit*>(h);
            FiberSite::lati_t s = i->site();
            if ( !i->valLattice(s) )
            {
                i->incLattice();
                i->setAbscissa(uni * s + i->prop()->site_shift);
            }
            else
            {
                std::cerr << "doubly occupied site? " << reference() << " @ " << s << '\n';
            }
        }
    }
    return true;
}
#endif

