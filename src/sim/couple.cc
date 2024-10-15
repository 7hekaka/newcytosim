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

//------------------------------------------------------------------------------

Couple::Couple(CoupleProp const* p, Vector const& w)
: prop(p), cPos(w), cHand1(nullptr), cHand2(nullptr)
{
    cHand1 = prop->hand1_prop->newHand(this);
    cHand2 = prop->hand2_prop->newHand(this);

    assert_true(w==w);
    assert_true(cHand1);
    assert_true(cHand2);
}


Couple::~Couple()
{
    if ( cHand1 )
    {
        if ( cHand1->attached() )
            cHand1->detachHand();
        delete(cHand1);
        cHand1 = nullptr;
    }
    
    if ( cHand2 )
    {
        if ( cHand2->attached() )
            cHand2->detachHand();
        delete(cHand2);
        cHand2 = nullptr;
    }
    
    prop = nullptr;
}


//------------------------------------------------------------------------------

/** This will fail if any hand is attached */
void Couple::changeProperty(CoupleProp * p)
{
    assert_true( p );
    assert_true( !cHand1->attached() && !cHand2->attached() );
    prop = p;
    
    if ( cHand1 && cHand1->property() != prop->hand1_prop )
    {
        delete(cHand1);
        cHand1 = prop->hand1_prop->newHand(this);
    }
    
    if ( cHand2 && cHand2->property() != prop->hand2_prop )
    {
        delete(cHand2);
        cHand2 = prop->hand2_prop->newHand(this);
    }
}

/**
 Returns the configuration of a crosslink, in discrete categories within [0, 6]
     Links on the side of the filaments:
         0 - Parallel if cos(filament1, filament2) > 0.5
         1 - Antiparallel if cos(filament1, filament2) < -0.5
         2 - X = none of the above
     Links at the ends of the filaments:
         3 - T-plus
         4 - V-plus
         5 - T-minus
         6 - V-minus
 by Jamie Li Rickman, Francis Crick Institute, London ~2017
 */
int Couple::configuration(real len, real max_cos) const
{
    int P = (cHand1->abscissaFrom(PLUS_END) < len) + (cHand2->abscissaFrom(PLUS_END) < len);
    int M = (cHand1->abscissaFrom(MINUS_END) < len) + (cHand2->abscissaFrom(MINUS_END) < len);
    
    if ( P > 0 )
        return 2+P; // T-plus and V-plus
    
    if ( M > 0 )
        return 4+M; // T-minus and V-minus

    real C = cosAngle();
    if ( C > max_cos ) // angle < PI/3, if max_cos = 0.5
        return 0; // P = parallel
    if ( C < -max_cos ) // angle > 2PI/3, if max_cos = 0.5
        return 1; // A = anti-parallel
    
    return 2; // X
}

//------------------------------------------------------------------------------
#pragma mark -


void Couple::setInteractions(Meca& meca) const
{
    assert_true( attached1() && attached2() );
    
    meca.addLink(cHand1->interpolation(), cHand2->interpolation(), prop->stiffness);
}


void Couple::setInteractionsAF(Meca& meca) const
{
    assert_true( attached1() && !attached2() );
}


void Couple::setInteractionsFA(Meca& meca) const
{
    assert_true( !attached1() && attached2() );
}

//------------------------------------------------------------------------------
#pragma mark -

void Couple::diffuse()
{
    //RNG.add_srand3(pos, cPos, prop->diffusion_dt);
    Vector pos = cPos + Vector::randS(prop->diffusion_dt);
    
    // confinement:
    if ( prop->confine == CONFINE_INSIDE )
    {
        /**
         @todo Dirichlet boundary conditions
         Set concentration of molecules at edges of Space by letting molecules
         out, and put some back at a constant rate
         */
        cPos = prop->confine_space->bounce(pos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        cPos = prop->confine_space->project(pos);
    }
    else
    {
        cPos = pos;
    }
}

/**
 Simulates:
 - diffusive motion
 - attachment of cHand1 and cHand2
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
        cHand1->stepUnattached(simul(), cPos);
    else if ( !prop->trans_activated )
        cHand2->stepUnattached(simul(), cPos);
}


/**
 Simulates:
 - attachment of cHand2
 - detachment of cHand1
 - attached activity of cHand1
 .
 */
void Couple::stepAF()
{
    //we use cHand1->pos() first, because stepUnloaded() may detach cHand1
    cHand2->stepUnattached(simul(), cHand1->outerPos());
    
    if ( cHand1->checkDetachment() )
        cHand1->detach();
    else
        cHand1->stepUnloaded();
}


/**
 Simulates:
 - attachment of cHand1
 - detachment of cHand2
 - attached activity of cHand2
 .
 */
void Couple::stepFA()
{
    //we use cHand2->pos() first, because stepUnloaded() may detach cHand2
    if ( !prop->trans_activated )
        cHand1->stepUnattached(simul(), cHand2->outerPos());
    
    if ( cHand2->checkDetachment() )
        cHand2->detach();
    else
        cHand2->stepUnloaded();
}


/**
 Simulates:
 - detachment of cHand1
 - attached activity of cHand1
 - detachment of cHand2
 - attached activity of cHand2
 .
 */
void Couple::stepAA()
{
    Vector f = Couple::force();
    real mag = f.norm();
    
    if ( cHand1->checkKramersDetachment(mag) )
        cHand1->detach();
    else
        cHand1->stepLoaded( f);
    
    /*
     Note that at this point cHand1 may be detached, but we still call stepLoaded()
      for cHand2 below, to keep the code symmetric for the Couple.
     This should be suitable as long detachment is a rare event, during one time-step
     */
    
    if ( cHand2->checkKramersDetachment(mag) )
        cHand2->detach();
    else
        cHand2->stepLoaded(-f);
}



/**
 Simulates only detachement & activity of cHand1
 */
void Couple::stepHand1()
{
    if ( cHand1->checkDetachment() )
        cHand1->detach();
    else
        cHand1->stepUnloaded();
}


/**
 Simulates only detachement & activity of cHand2
 */
void Couple::stepHand2()
{
    if ( cHand2->checkDetachment() )
        cHand2->detach();
    else
        cHand2->stepUnloaded();
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

bool Couple::permitAttachment(FiberSite const& sit, Hand const* h) const
{
    assert_true( h == cHand1 || h == cHand2 );
    FiberSite const* that = ( h == cHand1 ? cHand2 : cHand1 );
    
    Fiber const* fox = that->fiber();
    if ( ! fox )
        return true;
    
    Fiber const* fib = sit.fiber();
    
#if FIBER_HAS_FAMILY
    // prevent binding if that would make a link inside the same family
    if ( fox->family_  &&  fox->family_==fib->family_ )
        return false;
#endif

    // prevent binding to the same fiber at adjacent locations:
    if ( fib==fox  &&  abs_real(sit.abscissa()-that->abscissa()) <= prop->min_loop )
        return false;
    
#if ( 0 )
    /*
     Test here if binding would create a link inside an aster, near the center:
     i.e. a link between two Fibers from the same aster, very close to the center
     of this aster. Such links would be improductive, and would trap the Couples.
     */
    if ( fib->isBuddy(fox) )
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
            sit.reinterpolate();
            if ( dot(sit.dirFiber(), that->dirFiber()) < 0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_NOT_PARALLEL:
            sit.reinterpolate();
            if ( dot(sit.dirFiber(), that->dirFiber()) > 0.5 )
                return false;
            break;
  
        case CoupleProp::BIND_ANTIPARALLEL:
            sit.reinterpolate();
            if ( dot(sit.dirFiber(), that->dirFiber()) > -0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_NOT_ANTIPARALLEL:
            sit.reinterpolate();
            if ( dot(sit.dirFiber(), that->dirFiber()) < -0.5 )
                return false;
            break;
            
        case CoupleProp::BIND_ORTHOGONAL:
            sit.reinterpolate();
            if ( abs_real(dot(sit.dirFiber(), that->dirFiber())) > M_SQRT1_2 )
                return false;
            break;
            
        case CoupleProp::BIND_BIPOLAR: {
            Vector dir = normalize(that->pos() - fox->posMiddle());
            sit.reinterpolate();
            if ( dot(sit.dirFiber(), dir) < 0.5 )
                return false;
        } break;
            
        case CoupleProp::BIND_ANTIBIPOLAR: {
            Vector dir = normalize(that->pos() - fox->posMiddle());
            sit.reinterpolate();
            if ( dot(sit.dirFiber(), dir) > -0.5 )
                return false;
        } break;

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

/**
 When a Couple transitions into the unboud diffusing state, we set its
 position near the current location on the fiber, but offset in the perpendicular
 direction by a random distance within the range of attachment of the Hand.
 
 This is necessary to achieve detailed balance, which in particular implies
 that rounds of binding/unbinding should not get the Couples any closer to
 the Filaments.
 */
void Couple::beforeDetachment(Hand const* h)
{
    assert_true(h->attached());
    
    CoupleSet * set = static_cast<CoupleSet*>(objset());
    if ( set )
    {
        if ( h == cHand1 )
        {
            // cHand1 will detach
            if ( cHand2->unattached() )
            {
                cHand1->reinterpolate();
                cPos = cHand1->unbindingPosition();
            }
            set->relinkD1(this);
        }
        else
        {
            assert_true( h == cHand2 );
            // cHand2 will detach
            if ( cHand1->unattached() )
            {
                cHand2->reinterpolate();
                cPos = cHand2->unbindingPosition();
            }
            set->relinkD2(this);
        }
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
        cPos = prop->confine_space->placeOnEdge(1.0);
    else
        cPos = prop->confine_space->place();
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


Hand const* Couple::otherHand(Hand const* h) const
{
    if ( h == cHand1 )
        return cHand2;
    else
        return cHand1;
}


Vector Couple::linkFoot(Hand const* h) const
{
    if ( h == cHand1 )
    {
        if ( attached2() )
            return cHand2->pos();
        throw Exception("linkFoot() called for unattached Hand2");
    }
    else
    {
        if ( attached1() )
            return cHand1->pos();
        throw Exception("linkFoot() called for unattached Hand1");
    }
}

//------------------------------------------------------------------------------
#pragma mark -


void Couple::write(Outputter& out) const
{
    writeMarker(out, TAG);
    //std::clog << "- writing " << state() << " at " << out.pos() << '\n';
    cHand1->writeHand(out);
    cHand2->writeHand(out);
    if ( !attached1() && !attached2() )
        out.writeFloats(cPos, DIM);
}


/**
 To speedup reading, we could implement readFF(), readFA() readAF and readAA()
 Since Couple are stored seperatetly, depending of their state
 */
void Couple::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    ObjectID id1 = cHand1->readHand(in, sim);
    ObjectID id2 = cHand2->readHand(in, sim);
    
    if ( id1 || id2 )
    {
#if 0
        // it can be nice to set the position, but not essential
        if ( attached1() ) cHand1->reinterpolate();
        if ( attached2() ) cHand2->reinterpolate();
        cPos = position();
#endif
    }
    else
        in.readFloats(cPos, DIM);
}

