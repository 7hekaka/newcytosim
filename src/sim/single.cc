// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "dim.h"
#include "cymdef.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "glossary.h"
#include "iowrapper.h"
#include "single.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
Single::Single(SingleProp const* p, Vector const& w)
: sPos(w), sHand(nullptr), prop(p)
{
    assert_true(w==w);
    assert_true(prop->hand_prop);
    sHand = prop->hand_prop->newHand(this);
    assert_true(sHand);
}


Single::~Single()
{
    //std::clog << " delete single " << reference() << "\n";
    if ( linked() )
        objset()->remove(this);

    if ( sHand )
    {
        sHand->detachHand();
        delete(sHand);
        sHand = nullptr;
    }
    
    prop = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - HandMonitor

void Single::afterAttachment(Hand const*)
{
    assert_true( attached() );
    // link into correct SingleSet sublist:
    SingleSet * set = static_cast<SingleSet*>(objset());
    if ( set )
        set->relinkA(this);
}


/**
When a Single transitions into the unboud diffusing state, we set its
position near the current location on the fiber, but offset in the perpendicular
direction by a random distance within the range of attachment of the Hand.

This is necessary to achieve detailed balance, which in particular implies
that rounds of binding/unbinding should not get the Couples any closer to
the Filaments.
*/
void Single::beforeDetachment(Hand const* h)
{
    assert_true( h == sHand );
    
    SingleSet * set = static_cast<SingleSet*>(objset());
    if ( set )
    {
        sHand->reinterpolate();
        sPos = h->posSide();
        
        // link into correct SingleSet sublist:
        set->relinkD(this);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Functions


Vector Single::position() const
{
    if ( sHand->attached() )
        return sHand->pos();
    return sPos;
}

void Single::foldPosition(Modulo const* m)
{
    m->fold(sPos);
}


void Single::randomizePosition()
{
    if ( prop->diffusion > 0 || prop->fast_diffusion )
    {
        if ( prop->confine == CONFINE_ON )
            sPos = prop->confine_pointer->placeOnEdge(1.0);
        else if ( prop->confine == CONFINE_INSIDE )
            sPos = prop->confine_pointer->place();
        else if ( prop->confine != CONFINE_OFF )
            throw InvalidParameter("`confine` is incompatible with `fast_diffusion`");
    }
}


void Single::stepF()
{
    assert_false( sHand->attached() );

#if NEW_MOBILE_SINGLE
    // translation:
    sPos += prop->speed_dt;
#endif

    // diffusion:
    Vector pos = sPos + Vector::randS(prop->diffusion_dt);

    // confinement:
    if ( prop->confine == CONFINE_INSIDE )
    {
        sPos = prop->confine_pointer->bounce(pos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        sPos = prop->confine_pointer->project(pos);
    }
    else
    {
        sPos = pos;
    }
    
    sHand->stepUnattached(simul(), sPos);
}


/**
 The default Single has no force
 */
void Single::stepA()
{
    assert_true( sHand->attached() );
    assert_true( !hasLink() );

#if NEW_MOBILE_SINGLE
    // translation:
    sPos += prop->speed_dt;
#endif
    
    if ( sHand->checkDetachment() )
        sHand->detach();
    else
        sHand->stepUnloaded();
}

/**
 Add confinement force to the bound fiber
 */
void Single::setInteractions(Meca& meca) const
{
    assert_true( sHand->attached() );
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        Space const* spc = prop->confine_pointer;
        spc->setConfinement(sHand->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - I/O

void Single::write(Outputter& out) const
{
    writeMarker(out, TAG);
    sHand->writeHand(out);
    out.writeFloats(sPos, DIM);
}


/**
 To speedup reading, we could implement readF(), and readA()
 Since Single are stored seperatetly, depending of their state
 */
void Single::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    sHand->readHand(in, sim);
    in.readFloats(sPos, DIM);
}

