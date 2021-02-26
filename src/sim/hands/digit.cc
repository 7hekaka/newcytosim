// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "digit.h"
#include "digit_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "glossary.h"
#include "lattice.h"
#include "simul.h"


Digit::Digit(DigitProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
}

bool Digit::attachmentAllowed(FiberSite& sit) const
{
    if ( Hand::attachmentAllowed(sit) )
    {
#if FIBER_HAS_LATTICE
        FiberLattice& lat = sit.fiber()->lattice();
        
        if ( !lat.ready() )
            throw InvalidParameter("a lattice was not defined for `"+sit.fiber()->prop->name()+"'");
        
        // index to site containing given abscissa:
        lati_t s = lat.index(sit.abscissa());

        if ( lat.outsideMP(s) || occupied(lat, s) )
            return false;
        
        // adjust to match selected lattice site:
        sit.engageLattice(&lat, s, prop->site_shift);
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
    assert_true( attached() );
#if FIBER_HAS_LATTICE
    assert_true( lattice() );
    dec();
    hSite = s;
    inc();
    hAbs = s * lattice()->unit() + prop->site_shift;
#else
    hAbs = s * prop->step_size + prop->site_shift;
#endif
    update();
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
    
    if ( RNG.test(prop->hold_shrinking_end) )
    {
        jumpToEndM();
        if ( site() < lattice()->first() )
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
    
    if ( prop->hold_shrinking_end > 0 )
    {
        jumpToEndP();
        if ( site() >= lattice()->fence() )
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
    testDetachment();
}


/**
(see @ref Stochastic)
 */
void Digit::stepLoaded(Vector const& force, real force_norm)
{
    assert_true( attached() );
    assert_true( nextDetach >= 0 );

    testKramersDetachment(force_norm);
}

//------------------------------------------------------------------------------
#pragma mark -

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
        
        for ( Hand * ha = fHands.front(); ha; ha = ha->next() )
        {
            if ( ha->isDigit() && ha->lattice() == &fLattice )
            {
                Digit* dig = static_cast<Digit*>(ha);
                dig->inc();
                dig->moveTo(fLattice.unit() * dig->site()+ dig->prop->site_shift);
            }
        }
    }
}
#endif

