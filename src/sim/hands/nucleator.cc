// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "nucleator.h"
#include "nucleator_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "fiber_prop.h"
#include "fiber_set.h"
#include "simul.h"


//------------------------------------------------------------------------------

Nucleator::Nucleator(NucleatorProp const* p, HandMonitor* h)
: Hand(p,h)
{
    nextAct = RNG.exponential();
}

//------------------------------------------------------------------------------

ObjectList Nucleator::makeFiber(Simul& sim, Vector pos, std::string const& fiber_type, Glossary& opt)
{
    ObjectList objs;
    ObjectMark mk = 0;
    FiberProp * fip = sim.findProperty<FiberProp>("fiber", fiber_type);

    // determine direction of nucleation:
    Vector dir;
    Hand const* h = hMonitor->otherHand(this);
    if ( h && h->attached() )
    {
        // nucleating on the side of a 'mother' fiber:
        dir = h->dirFiber();
        // equalize marks to highlight amplification:
        mk = h->fiber()->mark();
        if ( mk == 0 )
        {
            mk = h->fiber()->identity();
            h->fiber()->mark(mk);
        }
        // remove key to avoid unused warning:
        opt.clear("direction");
    }
    else
    {
        std::string str;
        if ( opt.set(str, "direction") )
        {
            // nucleating in the bulk:
            std::istringstream ss(str);
            dir = Movable::readDirection(ss, pos, fip->confine_space_ptr);
        }
        else
        {
            // nucleating from a Mecable:
            dir = hMonitor->linkDir(this);
        }
    }
    // flip if nucleator is at the PLUS_END
    if ( prop()->hold_end == PLUS_END )
        dir = -dir;
    
    Fiber * fib = sim.fibers.newFiber(objs, fip, opt);
    // select rotation to align with direction of nucleation:
    Rotation rot = Rotation::randomRotationToVector(dir);
    
    const real A = prop()->nucleation_angle;
    const real L = hMonitor->linkRestingLength();
#if ( DIM == 2 )
    // shift position by the length of the interaction:
    real F = RNG.sflip();
    pos += rot * Vector(0, L*F, 0);
    rot = rot * Rotation::rotation(std::cos(A), std::sin(A)*F);
#elif ( DIM == 3 )
    // shift position by the length of the interaction:
    pos += rot * Vector(0, L, 0);
    rot = rot * Rotation::rotationAroundZ(A);
#endif

    // mark fiber to highlight mode of nucleation:
    opt.set(mk, "mark");
    Simul::mark(objs, mk);

    ObjectSet::rotateObjects(objs, rot);
    /*
     We translate Fiber to match the Nucleator's position,
     and if prop()->hold_end, the Hand is attached to the new fiber
     */
    if ( prop()->hold_end == MINUS_END )
    {
        attachEnd(fib, MINUS_END);
        pos -= fib->posEndM();
    }
    else if ( prop()->hold_end == PLUS_END )
    {
        attachEnd(fib, PLUS_END);
        pos -= fib->posEndP();
    }
    else
        pos -= fib->position();
    
    assert_true(pos.valid());
    ObjectSet::translateObjects(objs, pos);
    //std::clog << "nucleated fiber in direction " << fib->dirEndM() << "\n";

    opt.print_warnings(std::cerr, 1, "nucleator:spec\n");
    assert_true(fib->valid());
    return objs;
}


//------------------------------------------------------------------------------
/**
 Does not attach nearby Fiber, but can nucleate.
 the argument `pos` is the position of the other Hand
 */
void Nucleator::stepUnattached(Simul& sim, Vector const& pos)
{
    assert_false( attached() );
    
    nextAct -= prop()->rate_dt;
    
    if ( nextAct < 0 )
    {
        nextAct = RNG.exponential();
        try {
            Glossary opt(prop()->fiber_spec);
            ObjectList objs = makeFiber(sim, pos, prop()->fiber_type, opt);
            sim.add(objs);
        }
        catch( Exception & e )
        {
            e << "\nException occurred while executing nucleator:code";
            throw;
        }
    }
}


void Nucleator::stepUnloaded()
{
    assert_true( attached() );
    
    /// OPTION 1: delete entire fiber
    if ( prop()->addictive == 2 )
    {
        delete(fiber());
        return;
    }
    
    // may track the end of the Fiber:
    if ( prop()->track_end == MINUS_END )
        relocateM();
    else if ( prop()->track_end == PLUS_END )
        relocateP();
}


void Nucleator::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // may track the end of the Fiber:
    if ( prop()->track_end == MINUS_END )
        relocateM();
    else if ( prop()->track_end == PLUS_END )
        relocateP();
}


void Nucleator::detach()
{
    // if `addictive`, give a poisonous goodbye-kiss to the fiber
    if ( prop()->detached_end_state != STATE_BLACK )
        fiber()->setEndState(nearestEnd(), prop()->detached_end_state);
    else if ( prop()->addictive )
        fiber()->setEndState(nearestEnd(), STATE_RED);
    
    Hand::detach();
}

