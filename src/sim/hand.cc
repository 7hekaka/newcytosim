// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "hand.h"
#include "hand_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "fiber_prop.h"
#include "digit.h"
#include "simul.h"
#include "cymdef.h"


/// this HandMonitor does nothing
HandMonitor Hand::dummyMonitor;

//------------------------------------------------------------------------------

Hand::Hand(HandProp const* p, HandMonitor* m)
 : hNext(nullptr), hPrev(nullptr), hMonitor(m), prop(p)
{
    if ( !m )
        hMonitor = &dummyMonitor;
    // initialize in unattached state:
    nextDetach = 0;
}


Hand::~Hand()
{
    // the Hands should be detached in ~Couple and ~Single
    assert_true(!hFiber);
    prop = nullptr;
}


Hand * Hand::otherHand() const
{
    return hMonitor->otherHand(this);
}


Vector Hand::linkBase() const
{
    return hMonitor->linkBase(this);
}


real Hand::linkStiffness() const
{
    return hMonitor->linkStiffness();
}


void Hand::resetTimers()
{
    // initialize the Gillespie counters:
    if ( attached() )
    {
        nextDetach = RNG.exponential();
    }
    else
    {
        nextDetach = 0;
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Hand::relocate(Fiber* f)
{
    assert_true(f);
    if ( hFiber )
    {
        hFiber->removeHand(this);
        hFiber = f;
#if FIBER_HAS_LATTICE
        if ( hLattice )
            hLattice = &f->lattice();
#endif
    }
    f->addHand(this);
    update();
}


void Hand::relocate(Fiber* f, const real a)
{
    assert_true(f);
    if ( f != hFiber )
    {
        if ( hFiber )
        {
            hFiber->removeHand(this);
            hFiber = f;
#if FIBER_HAS_LATTICE
            if ( hLattice )
                hLattice = &f->lattice();
#endif
        }
        f->addHand(this);
    }
    hAbs = a;
    update();
}


void Hand::moveToEnd(const FiberEnd end)
{
    assert_true(hFiber);
    assert_true(end==PLUS_END || end==MINUS_END);
    
    if ( end == PLUS_END )
        FiberSite::relocateP();
    else
        FiberSite::relocateM();
}

//------------------------------------------------------------------------------
#pragma mark -

// only checks the Monitor's permission
bool Hand::monitorAllowsAttachment(FiberSite& sit) const
{
    return hMonitor->allowAttachment(sit, this);
}


/**
Checks that all the conditions required for attachment are met
 */
bool Hand::attachmentAllowed(FiberSite& sit) const
{
    assert_true( sit.attached() );
    
#if NEW_BINDING_LIMITS
    if ( sit.abscissa() < prop->binding_limits[0] )
        return false;
    if ( sit.abscissa() > prop->binding_limits[1] )
        return false;
#endif
    
    // check end-on binding:
    if ( sit.abscissaFromM() < 0 )
    {
        if ( prop->bind_also_end & MINUS_END )
            sit.relocateM();
        else
            return false;
    }
    else if ( sit.abscissaFromP() < 0 )
    {
        if ( prop->bind_also_end & PLUS_END )
            sit.relocateP();
        else
            return false;
    }
    
    FiberEnd end = NO_END;

    switch ( prop->bind_only_end )
    {
        case NO_END:
            break;
        case MINUS_END:
            if ( sit.abscissaFromM() > prop->bind_end_range )
                return false;       // too far from fiber end
            end = MINUS_END;
            break;
        case PLUS_END:
            if ( sit.abscissaFromP() > prop->bind_end_range )
                return false;       // too far from fiber end
            end = PLUS_END;
            break;
        case BOTH_ENDS:
        {
            if ( sit.abscissaFromM() > prop->bind_end_range )
            {
                // too far from MINUS_END
                if ( sit.abscissaFromP() > prop->bind_end_range )
                    return false;       // too far from PLUS_END
                end = PLUS_END;
            }
            else
            {
                // close from MINUS_END
                if ( sit.abscissaFromP() > prop->bind_end_range )
                    end = MINUS_END;    // too far from PLUS_END
                else
                    end = RNG.choice(MINUS_END, PLUS_END);
            }
        } break;
        default:
            throw Exception("Illegal value of hand:bind_only_end");
    }

#if NEW_BIND_ONLY_FREE_END
    // check occupancy near the end (should be done with FiberLattice)
    if ( end != NO_END && prop->bind_only_free_end )
    {
        if ( 0 < sit.fiber()->nbHandsNearEnd(prop->bind_end_range, end) )
            return false;
    }
#endif
    
    // also check the Monitor's permission:
    return hMonitor->allowAttachment(sit, this);
}


void Hand::locate(Fiber* f, real a)
{
    assert_true(f);
    assert_true(!hFiber);
    //assert_true(f->abscissaM() <= a + REAL_EPSILON);
    //assert_true(a <= f->abscissaP() + REAL_EPSILON);

    hAbs = a;
    hFiber = f;
    f->addHand(this);
    update();
    hMonitor->afterAttachment(this);
    nextDetach = RNG.exponential();
}


void Hand::attach(FiberSite const& s)
{
    assert_true(s.attached());
    assert_true(!hFiber);

    locate(s.fiber(), s.abscissa());
#if FIBER_HAS_LATTICE
    hLattice = s.lattice();
    hSite = s.site();
#endif
}


void Hand::detachHand()
{
    assert_true( attached() );
    hFiber->removeHand(this);
    hFiber = nullptr;
#if FIBER_HAS_LATTICE
    hLattice = nullptr;
#endif
}


void Hand::detach()
{
    assert_true( attached() );
    hMonitor->beforeDetachment(this);
    hFiber->removeHand(this);
    hFiber = nullptr;
#if FIBER_HAS_LATTICE
    hLattice = nullptr;
#endif
}

//------------------------------------------------------------------------------
#pragma mark -


void Hand::checkFiberRange()
{
    assert_true( attached() );
    
    if ( hAbs < hFiber->abscissaM() )
        handleDisassemblyM();
    else if ( hAbs > hFiber->abscissaP() )
        handleDisassemblyP();
}


void Hand::handleDisassemblyM()
{
    if ( RNG.test(prop->hold_shrinking_end) )
        relocateM();
    else
        detach();
}

void Hand::handleDisassemblyP()
{
    if ( RNG.test(prop->hold_shrinking_end) )
        relocateP();
    else
        detach();
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Test for attachment to nearby Fibers
 */
void Hand::stepUnattached(Simul& sim, Vector const& pos)
{
    assert_true( unattached() );

    // test for attachment
#if POOL_HAND_ATTACHMENT
    if ( 0 == sim.skipAttach )
#endif
        sim.fiberGrid.tryToAttach(pos, *this);
}


/**
 By default, the unloaded Hand does nothing
 */
void Hand::stepUnloaded()
{
    assert_true( attached() );
}


/**
 By default, the loaded Hand does nothing
 */
void Hand::stepLoaded(Vector const& force)
{
    assert_true( attached() );
}


//------------------------------------------------------------------------------
#pragma mark -


void Hand::write(Outputter& out) const
{
    /*
     it is not necessary to write the property number here,
     since it is set when the Hand is created in class Single or Couple.
     */
    FiberSite::write(out);
}


void Hand::read(Inputter& in, Simul& sim)
{
    Fiber * fib = hFiber;
    FiberSite::read(in, sim);
    resetTimers();
    
    // update fiber's lists:
    if ( fib != hFiber )
    {
        if ( fib )
            fib->removeHand(this);
        if ( hFiber )
            hFiber->addHand(this);
    }
    
#if FIBER_HAS_LATTICE
    if ( hFiber && !hLattice && isDigit() )
    {
        /* Promote a Digit class that is not bound to the lattice */
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


std::ostream& operator << (std::ostream& os, Hand const& arg)
{
    os << "hand(" << arg.fiber()->reference() << " " << arg.abscissa() << ")";
    return os;
}
