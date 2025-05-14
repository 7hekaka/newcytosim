// Cytosim was created by Francois Nedelec. Copyright 2025 NC State University.

#include "dim.h"
#include "assert_macro.h"
#include "bundle.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "interpolation.h"
#include "fiber_prop.h"
#include "glossary.h"
#include "simul.h"
#include "meca.h"


void Bundle::step()
{
    Simul & sim = simul();
    
    for ( size_t ii = 0; ii < nbOrganized(); ++ii )
    {
        if ( !organized(ii)  &&  RNG.test(prop->fiber_prob) )
        {
            ObjectList objs;
            FiberProp const* fip = sim.findProperty<FiberProp>("fiber", prop->fiber_type);
            Fiber * F = sim.fibers.newFiber(objs, fip, prop->fiber_spec);
            F->adjustLength(prop->overlap, prop->pole);
            ///\todo: we should orient the new Fiber in bundle direction
            sim.add(objs);
            grasp(F, ii);
        }
    }
}


Bundle::~Bundle()
{
    prop = nullptr;
}


/*
 Parallel connection near the 'pole' end of the fibers
*/
void Bundle::linkParallelFibers(Meca& meca, Fiber * mt1, Fiber * mt2) const
{
    const real stiff = prop->stiffness;
    const real dis = prop->overlap;
    
    const FiberEnd pole = prop->pole;
    meca.addLink(mt1->interpolateFrom(dis, pole), mt2->interpolateFrom(dis, pole), stiff);
    meca.addLink(mt1->exactEnd(pole), mt2->exactEnd(pole), stiff);
}


/**
 Antiparallel connection near the 'pole' end of the fibers
*/
void Bundle::linkAntiparallelFibers(Meca& meca, Fiber * mt1, Fiber * mt2) const
{
    const real stiff = prop->stiffness;
    const real dis = prop->overlap;

    const FiberEnd pole = prop->pole;
    if ( dis < REAL_EPSILON )
        meca.addLink(mt1->exactEnd(pole), mt2->exactEnd(pole), stiff+stiff);
    else {
        meca.addLink(mt2->exactEnd(pole), mt1->interpolateFrom(dis, pole), stiff);
        meca.addLink(mt1->exactEnd(pole), mt2->interpolateFrom(dis, pole), stiff);
    }
}


/**
 Connect the fibers near their ends, to form a ring:
 1. connect fibers with their neighbors,
 2. close the ring by connecting first and last fibers.
 */
void Bundle::setInteractions(Meca& meca) const
{
    if ( nbOrganized() > 1 )
    {
        const int S = ( prop->bipolar ? -1 : 1 );
        Fiber * mt0 = Fiber::toFiber(organized(0));
        Fiber * mt1 = mt0;
        
        int osi = 1, side = 1;
        for ( size_t i = 1; i < nbOrganized(); ++i )
        {
            Fiber * mt2 = Fiber::toFiber(organized(i));
            side *= S;
            if ( mt2 )
            {
                if ( side == osi )
                    linkParallelFibers(meca, mt1, mt2);
                else
                    linkAntiparallelFibers(meca, mt1, mt2);
                mt1 = mt2;
                osi = side;
            }
        }
        
        // connect first and last fibers:
        if ( mt1 && mt0 )
        {
            if ( side == osi )
                linkParallelFibers(meca, mt0, mt1);
            else
                linkAntiparallelFibers(meca, mt0, mt1);
        }
    }
}


//------------------------------------------------------------------------------

Vector Bundle::position() const
{
    Vector res(0,0,0);
    for ( size_t i = 1 ; i < nbOrganized(); ++i )
    {
        Fiber const* fib = Fiber::toFiber(organized(i));
        res += fib->posEnd(prop->pole);
    }
    return res / (real)nbOrganized();
}


/**
 It is possible to specify the lengths of individual fibers:

 new bundle bundle
 {
    length = 3.0, 4.2
 }

 */
ObjectList Bundle::build(Glossary& opt, Simul& sim)
{
    ObjectList objs;
    assert_true(prop);
    size_t cnt = 0;
    std::string type, spec;
    opt.set(cnt,  "fibers");
    opt.set(type, "fibers", 1);
    opt.set(spec, "fibers", 2);
    
    FiberProp const* fip = sim.findProperty<FiberProp>("fiber", type);

    if ( cnt <= 0 )
        throw InvalidParameter(prop->name()+":fibers[0] (number of fibers) must be specified and >= 1");
    
    nbOrganized(cnt);
    
    const FiberEnd pole = prop->pole;
    for ( size_t inx = 0; inx < cnt; ++inx )
    {
        ObjectList list;
        Fiber * fib = sim.fibers.newFiber(list, fip, spec);
        objs.append(list);
        
        // rotate odd fibers by 180 degrees to make an anti-parallel overlap:
        if ( inx & 1 )
            ObjectSet::rotateObjects(list, Rotation::rotation180());
        
        // translate to adjust the overlap:
        ObjectSet::translateObjects(list, fib->posMiddle()-fib->posFrom(0.5*prop->overlap, pole));
        
        real len;
        if ( opt.set(len, "length", inx) )
            fib->adjustLength(len, pole== PLUS_END?MINUS_END:PLUS_END);
        
        grasp(fib, inx);
    }
    objs.push_back(this);
    return objs;
}

void Bundle::write(Outputter& out) const
{
    writeMarker(out, Organizer::BUNDLE_TAG);
    writeOrganized(out);
}

