// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "dim.h"
#include "cymdef.h"
#include "couple.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "hand_prop.h"
#include "meca.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"
#include "aster.h"
#include "aster_prop.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

Couple::Couple(CoupleProp const* p, Vector const& w)
: prop(p), cPos(w), cHand1(nullptr), cHand2(nullptr)
{
    cHand1 = prop->hand1_prop->newHand(this);
    cHand2 = prop->hand2_prop->newHand(this);

    assert_true(cHand1);
    assert_true(cHand2);
}


Couple::~Couple()
{
    if ( cHand1 && attached1() )
        cHand1->detach();
    
    if ( cHand2 && attached2() )
        cHand2->detach();
    
    if ( linked() )
        objset()->remove(this);
    
    delete(cHand1);
    cHand1 = nullptr;
    delete(cHand2);
    cHand2 = nullptr;
    prop = nullptr;
}


//------------------------------------------------------------------------------

/** This will fail if any hand is attached */
void Couple::changeProperty(CoupleProp * p)
{
    assert_true( p );
    assert_true( !cHand1->attached() && !cHand2->attached() );
    prop = p;
    
    if ( cHand1 && cHand1->prop != prop->hand1_prop )
    {
        delete(cHand1);
        cHand1 = prop->hand1_prop->newHand(this);
    }
    
    if ( cHand2 && cHand2->prop != prop->hand2_prop )
    {
        delete(cHand2);
        cHand2 = prop->hand2_prop->newHand(this);
    }
}

/**
 Returns the configuration of a crosslink, in discrete categories
     Links on the side of the filaments:
         0 - Parallel if cos(filament1, filament2) > 0.5
         1 - Antiparallel if cos(filament1, filament2) < -0.5
         2 - X = none of the above
     Links at the ends of the filaments:
         3 - T-plus
         4 - V-plus
         5 - T-minus
         6 - V-minus
 by Jamie Li Rickman, ~2017
 */
int Couple::configuration(real len, real max_cos) const
{
    int p = (cHand1->abscissaFrom(PLUS_END) < len) + (cHand2->abscissaFrom(PLUS_END) < len);
    int m = (cHand1->abscissaFrom(MINUS_END) < len) + (cHand2->abscissaFrom(MINUS_END) < len);
    
    if ( p > 0 )
        return 2+p;  // T-plus and V-plus
    
    if ( m > 0 )
        return 4+m;  // T-minus and V-minus

    real c = cosAngle();
    if ( c > max_cos ) // angle < PI/3
        return 0; // P
    if ( c < -max_cos ) // angle > 2PI/3
        return 1; // A
    
    return 2; // X
}

//------------------------------------------------------------------------------
#pragma mark -


void Couple::setInteractions(Meca& meca) const
{
    assert_true( attached1() && attached2() );
    
    meca.addLink(cHand1->interpolation(), cHand2->interpolation(), prop->stiffness);
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        Space const* spc = prop->confine_space_ptr;
        spc->setConfinement(cHand1->interpolation(), meca, prop->stiffness, prop->confine);
        spc->setConfinement(cHand2->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}


void Couple::setInteractionsAF(Meca& meca) const
{
    assert_true( attached1() && !attached2() );
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        Space const* spc = prop->confine_space_ptr;
        spc->setConfinement(cHand1->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}


void Couple::setInteractionsFA(Meca& meca) const
{
    assert_true( !attached1() && attached2() );
    
#ifdef NEW_DANGEROUS_CONFINEMENTS
    if ( prop->confine )
    {
        Space const* spc = prop->confine_space_ptr;
        spc->setConfinement(cHand2->interpolation(), meca, prop->stiffness, prop->confine);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark -

void Couple::diffuse()
{
    Vector pos = cPos + Vector::randS(prop->diffusion_dt);
    
    // confinement:
    if ( prop->confine == CONFINE_INSIDE )
    {
        /**
         @todo dirichlet boundary conditions
         Set concentration of molecules at edges of Space by letting molecules
         out, and put some back at a constant rate
         */
        cPos = pos;
        prop->confine_space_ptr->bounce(cPos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        cPos = prop->confine_space_ptr->project(pos);
    }
}

/**
 Simulates:
 - diffusive motion
 - attachment
 .
 */
void Couple::stepFF()
{
    diffuse();

    /*
     For attachment, we select randomly one of the Hand, with equal chances,
     as if the Hands were occupying the two halves of a sphere moving very fast.
     Note that this divides by 2 the effective binding rate of the Hands.
     */
    if ( RNG.flip() )
    {
        cHand1->stepUnattached(simul(), cPos);
    }
    else
    {
        if ( !prop->trans_activated )
            cHand2->stepUnattached(simul(), cPos);
    }
}


/**
 Simulates:
 - attachment of cHand2
 - attached activity of cHand1
 .
 */
void Couple::stepAF()
{
    //we use cHand1->pos() first, because stepUnloaded() may detach cHand1
    cHand2->stepUnattached(simul(), cHand1->outerPos());
    
    if ( cHand1->testDetachment() )
        cHand1->stepUnloaded();
    else
        cHand1->detach();
}


/**
 Simulates:
 - attachment of cHand1
 - attached activity of cHand2
 .
 */
void Couple::stepFA()
{
    //we use cHand2->pos() first, because stepUnloaded() may detach cHand2
    cHand1->stepUnattached(simul(), cHand2->outerPos());
    
    if ( cHand2->testDetachment() )
        cHand2->stepUnloaded();
    else
        cHand2->detach();
}


/**
 Simulates:
 - attached activity of cHand1
 - attached activity of cHand2
 .
 */
void Couple::stepAA()
{
    Vector f = force();
    real fn = f.norm();
    
    if ( cHand1->testKramersDetachment(fn) )
        cHand1->stepLoaded( f);
    else
        cHand1->detach();
    
    if ( cHand2->testKramersDetachment(fn) )
        cHand2->stepLoaded(-f);
    else
        cHand2->detach();
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 @return:
 - True if attachment is possible
 - False if attachment is forbiden
 .
 
 If ( couple:stiff == true ), the two Hands of the Couple will refuse to be attached
 to the same segment, or to two neighboring segments on the same fiber.
 
 We cannot calculate the force of such 'degenerate' links, and they are undesired in
 most cases.
 
 */

bool Couple::allowAttachment(FiberSite const& sit, Hand const* h)
{
    assert_true( h == cHand1 || h == cHand2 );
    FiberSite const* that = ( h == cHand1 ? cHand2 : cHand1 );

    if ( !that->attached() )
        return true;
    
    Fiber const* fib = that->fiber();
    Fiber const* fob = sit.fiber();
    
#if FIBER_HAS_FAMILY
    // prevent binding if that would make a link inside the same family
    if ( fib->family_  &&  fib->family_==fob->family_ )
        return false;
#endif

    // prevent binding to the same fiber at adjacent locations:
    if ( fib==fob  &&  abs_real(sit.abscissa()-that->abscissa()) <= prop->min_loop )
        return false;
    
#if ( 0 )
    /*
     Test here if binding would create a link inside an aster, near the center:
     i.e. a link between two Fibers from the same aster, very close to the center
     of this aster. Such links would be improductive, and would trap the Couples.
     */
    const Buddy * bud = fib->buddy(0);
    if ( bud  &&  bud == fob->buddy(0) )
    {
        real a = that->abscissa();
        real b = sit.abscissa();
        if ( a < 1 && b < 1 )
            return false;
    }
#endif

    /*
     Allow or disallow binding based on the angle made between the two Fibers.
     The threshold on the cosine of the angle are here somewhat arbitrary
     */
    switch( prop->specificity )
    {
        case CoupleProp::BIND_ALWAYS:
            return true;
            
        case CoupleProp::BIND_PARALLEL:
            if ( dot(sit.dirFiber(), that->dirFiber()) < 0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_NOT_PARALLEL:
            if ( dot(sit.dirFiber(), that->dirFiber()) > 0.5 )
                return false;
            break;
  
        case CoupleProp::BIND_ANTIPARALLEL:
            if ( dot(sit.dirFiber(), that->dirFiber()) > -0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_NOT_ANTIPARALLEL:
            if ( dot(sit.dirFiber(), that->dirFiber()) < -0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_ORTHOGONAL:
            if ( abs_real(dot(sit.dirFiber(), that->dirFiber())) > 0.866025 )
                return false;
            break;
            
        default:
            throw InvalidParameter("unknown couple:specificity");
    }

    //attachment is allowed by default:
    return true;
}


void Couple::afterAttachment(Hand const* h)
{
    // link into correct CoupleSet sublist:
    CoupleSet * set = static_cast<CoupleSet*>(objset());
    if ( set )
    {
        if ( h == cHand1 )
            set->relinkA1(this);
        else
            set->relinkA2(this);
    }
}


void Couple::beforeDetachment(Hand const* h)
{
    assert_true(h->attached());
 
#if ( DIM < 2 )
    /*
     Relocate the Couple unbound position vector to where it is attached.
     This ensures that the diffusion process starts from the correct location
     */
    cPos = h->posHand();
#else
    /*
     Set position near the attachment point, but offset in the perpendicular
     direction at a random distance within the range of attachment of the Hand
     
     This is necessary to achieve detailed balance, which in particular implies
     that rounds of binding/unbinding should not get the Couples closer to
     the Filaments to which they bind.
     */
    if ( ! Couple::otherHand(h)->attached() )
        cPos = h->posHand() + h->dirFiber().randOrthoB(h->prop->binding_range);
#endif
    
    // link into correct CoupleSet sublist:
    CoupleSet * set = static_cast<CoupleSet*>(objset());
    if ( set )
    {
        if ( h == cHand1 )
            set->relinkD1(this);
        else
            set->relinkD2(this);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The position is:
 - cPos if the Couple is free,
 - the position of the attached Hand if only one is attached
 - the average position of the two hands if they are both attached
.
 */
Vector Couple::position() const
{
    if ( attached2() )
    {
        if ( attached1() )
            return 0.5 * ( cHand2->pos() + cHand1->pos() );
        return cHand2->pos();
    }
    if ( attached1() )
    {
        return cHand1->pos();
    }
    return cPos;
}


void Couple::foldPosition(Modulo const* m)
{
    m->fold(cPos);
}

void Couple::randomizePosition()
{
    if ( prop->confine == CONFINE_ON )
        cPos = prop->confine_space_ptr->randomPlaceOnEdge(1.0);
    else if ( prop->confine == CONFINE_INSIDE )
        cPos = prop->confine_space_ptr->randomPlace();
    else if ( prop->confine != CONFINE_OFF )
        throw InvalidParameter("`confine` is incompatible `fast_diffusion`");
}

//------------------------------------------------------------------------------
#pragma mark -

Vector Couple::stretch() const
{
    Vector d = cHand2->pos() - cHand1->pos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return d;
}


Vector Couple::force() const
{
    Vector d = cHand2->pos() - cHand1->pos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


Hand * Couple::attachedHand() const
{
    if ( attached1() )
        return cHand1;
    else if ( attached2() )
        return cHand2;
    else
        return nullptr;
}


Hand* Couple::otherHand(Hand const* h) const
{
    if ( h == cHand1 )
        return cHand2;
    else
        return cHand1;
}


Vector Couple::linkBase(Hand const* h) const
{
    if ( h == cHand1 )
    {
        if ( attached2() )
            return cHand2->pos();
        throw Exception("linkBase() called for unattached Hand2");
    }
    else
    {
        if ( attached1() )
            return cHand1->pos();
        throw Exception("linkBase() called for unattached Hand1");
    }
}


Vector Couple::otherDirection(Hand const* h) const
{
    if ( h == cHand1 )
    {
        if ( attached2() )
            return cHand2->dirFiber();
        else
            return Vector::randU();
    }
    else
    {
        if ( attached1() )
            return cHand1->dirFiber();
        else
            return Vector::randU();
    }
}


//------------------------------------------------------------------------------
#pragma mark -


void Couple::write(Outputter& out) const
{
    //std::clog << "- writing " << state() << " at " << out.pos() << '\n';
    cHand1->write(out);
    cHand2->write(out);
    if ( !attached1() && !attached2() )
        out.writeFloats(cPos, DIM);
}


void Couple::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    const bool a1 = cHand1->read(in, sim);
    const bool a2 = cHand2->read(in, sim);
    
    if ( a1 || a2 )
        ;//cPos = position();
    else
        in.readFloats(cPos, DIM);
}

