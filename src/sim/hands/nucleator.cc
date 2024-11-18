// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "nucleator.h"
#include "nucleator_prop.h"
#include "hand_monitor.h"
#include "primitives.h"
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

ObjectList Nucleator::createFiber(Simul& sim, Vector pos, FiberProp const* fip, Glossary& opt)
{
    // determine direction of nucleation:
    Vector dir(1, 0, 0);
    Hand const* h = otherHand();
    if ( h && h->attached() )
    {
        // nucleating on the side of a 'mother' fiber:
        dir = h->dirFiber();
        // remove key to avoid unused warning:
        opt.clear("direction");
    }
    else
    {
        std::string str;
        if ( opt.set(str, "direction") )
        {
            // nucleating in the bulk:
            dir = Cytosim::readDirection(str, pos, fip->confine_space);
        }
        else
        {
            // nucleating radially out from a Mecable, or randomly:
            dir = hMonitor->linkDir(this);
        }
    }
    // flip direction if nucleator will stay at the plus end:
    if ( prop()->hold_end == PLUS_END )
        dir.negate();
    
    if ( prop()->specificity == NucleatorProp::NUCLEATE_MOSTLY_PARALLEL )
    {
        if ( RNG.flip_8th() )
            dir.negate();
    }

    // specified angle between extant and nucleated filament:
    const real A = prop()->nucleation_angle;
    const real L = hMonitor->linkRestingLength();
    // can flip the side in 2D or when nucleating 'in-plane'
    const real F = RNG.sflip();

    ObjectList objs;
    Fiber * fib = sim.fibers.newFiber(objs, fip, opt);
    // select rotation to align with direction of nucleation:
    Rotation rot(0, 1);
#if ( DIM >= 3 )
    if ( prop()->nucleate_in_plane )
    {
        Space const* spc = fip->confine_space;
        Vector out = spc->normalToEdge(pos);
        // make 'dir' tangent to the Space's edge
        dir = out.orthogonal(dir, 1.0);
        rot = Rotation::rotationAroundZ(std::cos(A), F*std::sin(A));
        rot = Rotation::rotationToVectors(dir, out) * rot;
        // shift position sideways by the length of the interaction:
        pos += rot * Vector(0, F*L, 0);
        rot = rot * Rotation::rotationAroundX(M_PI*RNG.sreal());
     }
    else
#endif
    {
        rot = Rotation::randomRotationToVector(dir);
#if ( DIM == 2 )
        rot = rot * Rotation::rotation(std::cos(A), F*std::sin(A));
#elif ( DIM == 3 )
        rot = rot * Rotation::rotationAroundZ(A);
#endif
        // shift position sideways by the length of the interaction:
        pos += rot * Vector(0, F*L, 0);
    }

    // mark fiber to highlight mode of nucleation:
    ObjectMark mk = 0;
    if ( opt.value_is("mark", 0, "random") )
        mk = RNG.pint32();
    else opt.set(mk, "mark");
    fib->mark(mk);

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
#if 0
    real a = std::acos(dot(fib->dirEndM(), dir));
    std::clog << "nucleated with angle " << std::setw(8) << a << " along " << fib->dirEndM() << "\n";
#endif
    opt.print_warnings(stderr, 1, " in nucleator:spec\n");
    assert_false(fib->invalid());
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
    FiberProp const* fip = prop()->fiber_class;
    
    float damp = 1.f - float(fip->nbFibers()) * prop()->nucleation_limit;

    float R = prop()->nucleation_rate_dt * max_float(0, damp);
    nextAct -= R;
    
    if ( nextAct < 0 )
    {
        nextAct = RNG.exponential();
        try {
            Glossary opt(prop()->fiber_spec);
            ObjectList objs = createFiber(sim, pos, fip, opt);
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
    
    if ( prop()->stabilize != 0 )
    {
        Fiber * fib = modifiableFiber();
        fib->stabilize(nearestEnd(), prop()->stabilize);
    }
}


void Nucleator::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    // may track the end of the Fiber:
    if ( prop()->track_end == MINUS_END )
        relocateM();
    else if ( prop()->track_end == PLUS_END )
        relocateP();
    
    if ( prop()->stabilize != 0 )
    {
        Fiber * fib = modifiableFiber();
        fib->stabilize(nearestEnd(), prop()->stabilize);
    }
}


void Nucleator::detach()
{
    // if `addictive`, give a poisonous goodbye kiss to the fiber
    if ( prop()->addictive )
    {
        Fiber * fib = modifiableFiber();
        fib->setEndState(nearestEnd(), prop()->addictive_state);
    }
    Hand::detach();
}

