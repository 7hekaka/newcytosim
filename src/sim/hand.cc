// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "hand.h"
#include "hand_prop.h"
#include "hand_monitor.h"
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
: next_(nullptr), prev_(nullptr), prop(p), hMonitor(m), nextDetach(0), nextAct(0)
{
    if ( !m )
        hMonitor = &dummyMonitor;
}


Hand::~Hand()
{
    // the Hands should be detached in ~Couple and ~Single
    assert_true(!hFiber);
    prop = nullptr;
}


Hand const* Hand::otherHand() const
{
    return hMonitor->otherHand(this);
}


Vector Hand::linkFoot() const
{
    return hMonitor->linkFoot(this);
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

Vector Hand::unbindingPosition() const
{
    // needed if the interpolation is not up-to-date
    //reinterpolate();
#if ( DIM < 2 )
    /*
     Relocate the Couple unbound position vector to where it is attached.
     This ensures that the diffusion process starts from the correct location
     */
    return pos();
#else
    /*
     Set position near the attachment point, but offset in the perpendicular
     direction at a random distance within the range of attachment of the Hand
     
     This is necessary to achieve detailed balance, which in particular implies
     that rounds of binding/unbinding should not get the Couples closer to
     the Filaments, even after successive rounds of binding / unbinding.
     */
    if ( distanceToNearestEnd() < 0.008 )
    {
        // if direct binding to the end is permitted, we need to reverse this:
        FiberEnd end = nearestEnd();
        if ( prop->bind_also_end & end )
        {
            Vector dir = hFiber->dirEnd(end);
            Vector off = Vector::randB(prop->binding_range);
            real x = dot(dir, off);
            off = off + dir * ( abs_real(x) - x );
            return pos() + off;
        }
    }
    
    return pos() + dirFiber().randOrthoB(prop->binding_range);
#endif
}


void Hand::relocate(Fiber const* f, const real a)
{
    assert_true(f);
    if ( hFiber != f )
    {
        if ( hFiber ) hFiber->removeHand(this);
        f->addHand(this);
        hFiber = f;
    }
#if FIBER_HAS_LATTICE
    if ( hLattice )
        hLattice = f->lattice();
#endif
    hAbs = a;
    reinterpolate(f->interpolateAbs(a));
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

/**
Checks that all the conditions required for attachment are met
 */
bool Hand::attachmentAllowed(FiberSite& sit) const
{
    assert_true( sit.attached() );
    
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
    
    [[maybe_unused]] FiberEnd end = NO_END;

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
                // too far from minus end
                if ( sit.abscissaFromP() > prop->bind_end_range )
                    return false;       // too far from PLUS_END
                end = PLUS_END;
            }
            else
            {
                // close from minus end
                if ( sit.abscissaFromP() > prop->bind_end_range )
                    end = MINUS_END;    // too far from plus end
                else
                    end = RNG.choice(MINUS_END, PLUS_END);
            }
        } break;
        case CENTER:
            if ( abs_real(sit.abscissaFromC()) > 0.5 * prop->bind_end_range )
                return false;       // too far from fiber center
            break;
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
    
    // finally check the Monitor's permission:
    return hMonitor->permitAttachment(sit, this);
}

/**
 This is where member variables are really updated,
 normally all derived Hand::attach() should call this.
 */
void Hand::do_attach(Fiber const* f, real a)
{
    assert_true(f);
    assert_true(!hFiber);
    assert_true(f->abscissaM() <= a + REAL_EPSILON);
    assert_true(a <= f->abscissaP() + REAL_EPSILON);

    hAbs = a;
    hFiber = f;
    f->addHand(this);
    reinterpolate(f->interpolateAbs(a));
    hMonitor->afterAttachment(this);
    nextDetach = RNG.exponential();
    assert_true(nextDetach > 0);
#if FIBER_HAS_LATTICE
    if ( hLattice )
        incLattice();
#endif
}


void Hand::attach(FiberSite const& s)
{
    assert_true(!hFiber);

    do_attach(s.fiber(), s.abscissa());
#if 0
    Hand const* h = otherHand();
    if ( h && h->attached() )
    {
        real x = dot( h->pos() - s.pos(), s.dirFiber());
        std::clog << "attach " << s << " at " << x << " um\n";
    }
#endif
}


void Hand::detachHand()
{
    assert_true( attached() );
    hFiber->removeHand(this);
    hFiber = nullptr;
#if FIBER_HAS_LATTICE
    if ( hLattice ) {
        decLattice();
        hLattice = nullptr;
    }
#endif
}


void Hand::detach()
{
    assert_true( attached() );
    hMonitor->beforeDetachment(this);
    hFiber->removeHand(this);
#if FIBER_HAS_LATTICE
    if ( hLattice ) {
        decLattice();
        hLattice = nullptr;
    }
#endif
    hFiber = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -


void Hand::checkFiberRange(real absM, real absP)
{
    assert_true( attached() );
    
    if ( hAbs < absM )
        handleDisassemblyM();
    else if ( hAbs > absP )
        handleDisassemblyP();
}


void Hand::handleDisassemblyM()
{
    if ( RNG.test(prop->hold_shrinking_end[1]) )
        relocateM();
    else
        detach();
}

void Hand::handleDisassemblyP()
{
    if ( RNG.test(prop->hold_shrinking_end[0]) )
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
#pragma mark - I/O


ObjectID Hand::readHand(Inputter& in, Simul& sim)
{
    Fiber const* fib = hFiber;
    ObjectID id = readFiberSite(in, sim);
    resetTimers();
    
    // update fiber's lists:
    if ( fib != hFiber )
    {
        if ( fib )
            fib->removeHand(this);
        if ( hFiber )
            hFiber->addHand(this);
    }
    return id;
}


std::ostream& operator << (std::ostream& os, Hand const& arg)
{
    os << "hand(" << arg.fiber()->reference() << " " << arg.abscissa() << ")";
    return os;
}
