// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "digit.h"
#include "digit_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "glossary.h"
#include "lattice.h"
#include "simul_part.h"


Digit::Digit(DigitProp const* p, HandMonitor* h)
: Hand(p,h)
{
}

/**
 This sets 'sit.site()' which is thus modified in a meaningful way!
 */
bool Digit::attachmentAllowed(FiberSite& sit) const
{
    if ( Hand::attachmentAllowed(sit) )
    {
#if FIBER_HAS_LATTICE
        FiberLattice* lat = sit.fiber()->lattice();
        
        if ( !lat->ready() )
            throw InvalidParameter("a lattice was not defined for `"+sit.fiber()->prop->name()+"'");
        
        // index to site containing given abscissa:
        lati_t s = lat->index(sit.abscissa());
        
        if ( s < lat->entry() )
        {
            if ( prop->bind_also_end & MINUS_END )
                s = lat->entry();
            else
                return false;
        }
        else if ( lat->fence() < s )
        {
            if ( prop->bind_also_end & PLUS_END )
                s = lat->fence();
            else
                return false;
        }
        
        if ( occupied(lat, s) )
            return false;
        
        // adjust to match selected lattice site:
        sit.engageLattice(lat, s, prop->site_shift);
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
    Hand::attach(sit);
    inc();
}


void Digit::detach()
{
    dec();
    Hand::detach();
}

//------------------------------------------------------------------------------
#pragma mark -

void Digit::hop(lati_t s)
{
    assert_true(attached());
#if FIBER_HAS_LATTICE
    assert_true(hLattice);
    //std::clog << this << ":hop " << hSite << " ---> " << s << "\n";
    dec();
    hSite = s;
    inc();
    hAbs = s * lattice()->unit() + prop()->site_shift;
#else
    hAbs = s * prop()->step_size + prop()->site_shift;
#endif
}

//------------------------------------------------------------------------------

/**
 Try to move `n` sites in the PLUS_END direction,
 stopping if any intermediate position is already occupied.
 */
void Digit::crawlP(const int n)
{
    lati_t s = site();
    lati_t e = s + n;
    
    while ( s < e )
    {
        if ( vacant(s+1) )
            ++s;
        else
            break;
    }
    if ( s != site() )
        hop(s);
}


/**
 Try to move `n` sites in the MINUS_END direction,
 stopping if any intermediate position is already occupied.

 */
void Digit::crawlM(const int n)
{
    lati_t s = site();
    lati_t e = s - n;
    
    while ( s > e )
    {
        if ( vacant(s-1) )
            --s;
        else
            break;
    }
    if ( s != site() )
        hop(s);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 The Digit normally does not move by itself
 */
void Digit::handleDisassemblyM()
{
    assert_true( attached() );
    
    if ( RNG.test(prop()->hold_shrinking_end) )
    {
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
    
    if ( prop()->hold_shrinking_end > 0 )
    {
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
            static_cast<Digit*>(this)->inc();
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
void Fiber::resetLattice()
{
    if ( fLattice.data() )
    {
        fLattice.clear();
        real unit = fLattice.unit();
    
        for ( Hand * ha = fHands.front(); ha; ha = ha->next() )
        {
            if ( ha->lattice() == lattice() )
            {
                Digit* i = static_cast<Digit*>(ha);
                i->inc();
                i->moveTo(unit * i->site() + i->prop()->site_shift);
            }
        }
    }
}
#endif

