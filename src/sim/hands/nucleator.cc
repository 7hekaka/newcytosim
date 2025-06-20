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
    ObjectMark mk = 0;
    Vector dir(1, 0, 0);
    Hand const* h = otherHand();
    
    const real L = hMonitor->linkRestingLength();
    // specified angle between extant and nucleated filament:
    real A = 0;
    // can flip the side in 2D or when nucleating 'in-plane'
    const real F = RNG.sflip();

    // determine direction of nucleation:
    if ( h && h->attached() )
    {
        // deviation angle from 'dir'
        A = prop()->branch_angle;
        // get mark from mother fiber:
        mk = h->fiber()->mark();
        // nucleating on the side of a 'mother' fiber:
        switch ( prop()->branch_direction )
        {
            case NucleatorProp::BRANCH_PARALLEL:
                dir = h->dirFiber(); break;
            case NucleatorProp::BRANCH_MOSTLY_PARALLEL:
                dir = h->dirFiber() * ( RNG.flip_8th() ? -1: +1 ); break;
            case NucleatorProp::BRANCH_ANTI_PARALLEL:
                dir = -h->dirFiber(); break;
            case NucleatorProp::BRANCH_RANDOM:
                dir = h->dirFiber(); A = M_PI * RNG.preal(); break;
            case NucleatorProp::BRANCH_SPECIFIED:
                dir = Cytosim::readDirection(opt.value("direction", 0), pos, fip->confine_space); break;
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
            dir = Cytosim::readDirection(str, pos, fip->confine_space);
        }
        else
        {
            // nucleating radially out from a Mecable, or randomly:
            dir = hMonitor->linkDir(this);
        }
        // remove key to avoid unused warning:
        opt.clear("branch_direction");
    }
    // flip direction if nucleator will stay at the plus end:
    if ( prop()->hold_end == PLUS_END )
        dir.negate();

    ObjectList objs;
    Rotation rot(0, 1);
    Fiber * fib = sim.fibers.newFiber(objs, fip, opt);
    
    // select rotation to align with direction of nucleation:
#if ( DIM >= 3 )
    Space const* spc = prop()->nucleate_space;
    if ( spc )
    {
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
        if ( dir.normSqr() > 0.01 )
            rot = Rotation::randomRotationToVector(dir);
        else
            rot = Rotation::randomRotation();
#if ( DIM == 2 )
        rot = rot * Rotation::rotation(std::cos(A), F*std::sin(A));
#elif ( DIM == 3 )
        rot = rot * Rotation::rotationAroundZ(A);
#endif
        // shift position sideways by the length of the interaction:
        pos += rot * Vector(0, F*L, 0);
    }

    // mark fiber to highlight mode of nucleation:
    if ( opt.value_is("mark", 0, "random") )
        mk = RNG.pint32();
    else if ( opt.value_is("mark", 0, "add") )
        ++mk;
    else opt.set(mk, "mark");
    fib->mark(mk);

    // set signature of fiber:
    ObjectSignature S = 0;
    if ( opt.value_is("signature", 0, "base") )
        fib->signature(hMonitor->baseSignature());
    else if ( opt.set(S, "signature") )
        fib->signature(S);

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
    
    // factor to limit the number of fibers nucleated
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

