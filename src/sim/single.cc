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
    assert_true(prop->hand_prop);
    sHand = prop->hand_prop->newHand(this);
    assert_true(sHand);
}


Single::~Single()
{
    if ( linked() )
        objset()->remove(this);

    if ( sHand )
    {
        if ( sHand->attached() )
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


void Single::beforeDetachment(Hand const* h)
{
    assert_true( h == sHand );
    
#if ( DIM < 2 )
    /*
     Relocate Single to the position where it is attached.
     This is necessary to start the diffusion process from the correct location
     */
    sPos = h->posHand();
#else
    /*
     Set position near the attachment point, but offset in the perpendicular
     direction at a random distance within the range of attachment of the Hand.
     
     This is necessary to achieve detailed balance, which in particular implies
     that rounds of binding/unbinding should not get the Singles closer to
     the Filaments to which they bind.
     */
    sPos = h->posHand() + h->dirFiber().randOrthoB(h->prop->binding_range);
#endif

    // link into correct SingleSet sublist:
    SingleSet * set = static_cast<SingleSet*>(objset());
    if ( set )
        set->relinkD(this);
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
    if ( prop->confine == CONFINE_ON )
        sPos = prop->confine_space_ptr->randomPlaceOnEdge(1.0);
    else if ( prop->confine == CONFINE_INSIDE )
        sPos = prop->confine_space_ptr->randomPlace();
    else if ( prop->confine != CONFINE_OFF )
        throw InvalidParameter("`confine` is incompatible with `fast_diffusion`");
}


void Single::stepF()
{
    assert_false( sHand->attached() );

#if NEW_MOBILE_SINGLE
    // translation:
    sPos += prop->speed_dt;
#endif

    // diffusion:
    sPos.addRand(prop->diffusion_dt);
    
    // confinement:
    if ( prop->confine == CONFINE_INSIDE )
    {
        prop->confine_space_ptr->bounce(sPos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        sPos = prop->confine_space_ptr->project(sPos);
    }
    
    sHand->stepUnattached(simul(), sPos);
}


void Single::stepA()
{
    assert_true( sHand->attached() );
    assert_true( !hasLink() );

#if NEW_MOBILE_SINGLE
    // translation:
    sPos += prop->speed_dt;
#endif
    
    if ( sHand->testDetachment() )
        sHand->stepUnloaded();
    else
        sHand->detach();
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
        Space const* spc = prop->confine_space_ptr;
        spc->setConfinement(sHand->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - I/O

void Single::write(Outputter& out) const
{
    sHand->write(out);
    out.writeFloats(sPos, DIM);
}


void Single::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    sHand->read(in, sim);
    in.readFloats(sPos, DIM);
}

