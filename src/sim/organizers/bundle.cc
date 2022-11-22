// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
            Glossary opt(prop->fiber_spec);
            Fiber * fib = sim.fibers.newFiber(objs, prop->fiber_type, opt);
            fib->adjustLength(prop->overlap, prop->joint);
            ///\todo: we should orient the new Fiber in bundle direction
            sim.add(objs);
            grasp(fib, ii);
        }
    }
}


Bundle::~Bundle()
{
    prop = nullptr;
}


/*
 Parallel connection near the 'prop->joint' end of the fibers
*/
void Bundle::linkParallel(Meca& meca, Fiber * mt1, Fiber * mt2) const
{
    const real stiff = prop->stiffness;
    const real dis = prop->overlap;
    
    meca.addLink(mt1->interpolateFrom(dis, prop->joint), mt2->interpolateFrom(dis, prop->joint), stiff);
    meca.addLink(mt1->exactEnd(prop->joint), mt2->exactEnd(prop->joint), stiff);
}


/**
 Antiparallel connection near the 'prop->joint' end of the fibers
*/
void Bundle::linkAntiparallel(Meca& meca, Fiber * mt1, Fiber * mt2) const
{
    const real stiff = prop->stiffness;
    const real dis = prop->overlap;

    if ( dis < REAL_EPSILON )
        meca.addLink(mt1->exactEnd(prop->joint), mt2->exactEnd(prop->joint), stiff+stiff);
    else {
        meca.addLink(mt2->exactEnd(prop->joint), mt1->interpolateFrom(dis, prop->joint), stiff);
        meca.addLink(mt1->exactEnd(prop->joint), mt2->interpolateFrom(dis, prop->joint), stiff);
    }
}


/**
 Connect the fibers near their ends, to form a ring:
 1. connect fibers with their neighbors,
 2. close the ring by connecting first and last fibers.
 */
void Bundle::setInteractions(Meca& meca) const
{
    if ( nbOrganized() > 0 )
    {
        Fiber * mt0 = Fiber::toFiber(organized(0));
        Fiber * mt1 = mt0, * mt2 = nullptr;
        
        for ( size_t i = 1 ; i < nbOrganized(); ++i )
        {
            mt2 = Fiber::toFiber(organized(i));
            if ( mt1 && mt2 )
                linkAntiparallel(meca, mt1, mt2);
            mt1 = mt2;
        }
        
        // connect first and last fibers:
        mt1 = mt0;
        
        if ( mt1 && mt2 )
        {
            if ( nbOrganized() & 1 )
                linkParallel(meca, mt1, mt2);
            else
                linkAntiparallel(meca, mt1, mt2);
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
        res += fib->posEnd(prop->joint);
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
    
    if ( cnt <= 0 )
        throw InvalidParameter(prop->name()+":fibers[0] (number of fibers) must be specified and >= 1");
    
    nbOrganized(cnt);
    
    for ( size_t inx = 0; inx < cnt; ++inx )
    {
        ObjectList list;
        Glossary fiber_opt(spec);
        Fiber * fib = sim.fibers.newFiber(list, type, fiber_opt);
        objs.append(list);
        
        // rotate odd fibers by 180 degrees to make an anti-parallel overlap:
        if ( inx & 1 )
            ObjectSet::rotateObjects(list, Rotation::rotation180());
        
        // translate to adjust the overlap:
        ObjectSet::translateObjects(list, fib->posMiddle()-fib->posFrom(0.5*prop->overlap, prop->joint));
        
        real len;
        if ( opt.set(len, "length", inx) )
            fib->adjustLength(len, prop->joint== PLUS_END?MINUS_END:PLUS_END);
        
        grasp(fib, inx);
    }
    objs.push_back(this);
    return objs;
}

void Bundle::write(Outputter& out) const
{
    writeHeader(out, Organizer::TAG_BUNDLE);
    writeOrganized(out);
}

