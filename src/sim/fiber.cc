// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "fiber.h"
#include "field.h"
#include "messages.h"
#include "glossary.h"
#include "iowrapper.h"
#include "fiber_segment.h"
#include "fiber_prop.h"
#include "object_set.h"
#include "single.h"
#include "simul.h"
#include "meca.h"
#include "hand.h"
#include "cymdef.h"


void Fiber::step()
{
    assert_small( length1() - length() );
#if FIBER_HAS_GLUE
    if ( prop->glue )
        setGlue(fGlue, PLUS_END, prop->glue);
#endif
    
    //assert_false(hasKink(0));
#if ( 0 )
    // Cut kinked filaments
    size_t p = hasKink(0);
    if ( p )
    {
        LOG_ONCE("SEVER_KINKED_FIBERS\n");
        objset()->add(severPoint(p));
    }
#endif
#if ( 0 )
    // Delete kinked filaments
    if ( hasKink(0) )
    {
        LOG_ONCE("DELETE_KINKED_FIBERS\n");
        delete(this);
        return;
    }
#endif

    // perform the cuts that were registered by sever()
    if ( pendingCuts.size() )
        severNow();
    
    // delete self if shorter than 'FiberProp::min_length'
    if ( length() < prop->min_length && ! prop->persistent )
    {
        delete(this);
        return;
    }
    
    if ( needUpdate )
    {
        adjustSegmentation();
        updateFiber();
        // We may need to update again, as Hand may act in handleDisassemblyM()
        if ( needUpdate )
            updateFiber();
#if FIBER_HAS_FAMILY
        salute();
#endif
    }
    else
        updateHands();
    assert_false(needUpdate);
    
#if FIBER_HAS_MESH
    
    //std::clog << reference() << " mesh " << std::fixed << fMesh->data(0) << '\n';
    //std::clog << reference() << " mesh sum = " << fMesh->sum() << "     ";
    const real tau = time_step(simul());
    
    if ( prop->mesh_binding_rate > 0 || prop->mesh_unbinding_rate > 0 )
    {
        real on  = prop->mesh_binding_rate * tau;
        real off = -std::expm1( -prop->mesh_unbinding_rate * tau );
        equilibrateMesh(fMesh, prop->field_ptr, on, off);
    }
    
    if ( prop->mesh_flux_speed != 0 )
        fluxMesh(fMesh, prop->field_ptr, prop->mesh_flux_speed);
    
    if ( prop->mesh_cut_fiber )
        cutFiberMesh(fMesh);
    
    if ( prop->mesh_aging_rate > 0 )
    {
        /*
         This implements an evolution toward equilibrium,
         starting from a value of 0, and reached for a value of 1,
         with a time-scale given by 'mesh_aging_rate'.
         */
        real cst = prop->mesh_aging_rate * tau;
        evolveMeshValues(fMesh, cst, 1 - cst);
        //std::clog << reference() << " lattice avg = " << fMesh->sum()*fMesh->unit()/length() << '\n';
    }
    //std::clog << fMesh->sum() << '\n';
#endif
}



void Fiber::step(real addM, real addP, bool split)
{
    //std::clog << reference() << " P " << addP << " M " << addM << " len " << length() << "\n";

    // reset increments; they will be set again by growP and growM
    cDeltaM = 0;
    cDeltaP = 0;
    
#if NEW_FIBER_CHEW
    if ( split )
    {
        if ( fChewM > 0 )
            addM -= std::min(fChewM, prop->max_chewing_speed_dt);
        fChewM = 0;
        if ( fChewP > 0 )
            addP -= std::min(fChewP, prop->max_chewing_speed_dt);
        fChewP = 0;
    }
#endif
    
    real inc = addM + addP;
    real len = length() + inc;

    if ( inc < 0 )
    {
        if ( len > prop->min_length )
        {
            if ( addM != 0 ) growM(addM);
            if ( addP != 0 ) growP(addP);
        }
        else
        {
            if ( !prop->persistent )
            {
                delete(this);
                return;
            }
            else
            {
                // set dormant state:
                if ( addM < 0 ) setEndStateM(STATE_WHITE);
                if ( addP < 0 ) setEndStateP(STATE_WHITE);
                
                if ( split && length() > prop->min_length )
                {
                    // the remaining possible shrinkage is distributed to the two ends:
                    inc = ( prop->min_length - length() ) / inc;
                    if ( addM != 0 ) growM(addM*inc);
                    if ( addP != 0 ) growP(addP*inc);
                }
            }
        }
    }
    else if ( inc > 0 )
    {
        if ( len < prop->max_length )
        {
            if ( addM != 0 ) growM(addM);
            if ( addP != 0 ) growP(addP);
        }
        else
        {
            if ( split &&  length() < prop->max_length )
            {
                // the remaining possible growth is distributed to the two ends:
                inc = ( prop->max_length - length() ) / inc;
                if ( addM != 0 ) growM(addM*inc);
                if ( addP != 0 ) growP(addP*inc);
            }
        }
    }
    
    Fiber::step();
}

//------------------------------------------------------------------------------

/**
 The rest of the initialization is done in FiberProp::newFiber(),
 and other newFiber() functions where the initial length is known.
 */
Fiber::Fiber(FiberProp const* p)
: prop(p), disp(nullptr)
{
    if ( prop )
    {
        targetSegmentation(prop->segmentation);
        
        if ( prop->lattice )
        {
#if FIBER_HAS_LATTICE
            if ( prop->lattice_unit < REAL_EPSILON )
                throw InvalidParameter("the Lattice unit (lattice[1]) must be > 0");
            //Cytosim::log << reference() <<  " new Lattice" << '\n';
            fLattice.setUnit(prop->lattice_unit);
#else
            //throw InvalidParameter("Cytosim does not support fiber:lattice");
#endif
        }
        
        if ( prop->mesh )
        {
#if FIBER_HAS_MESH
            if ( prop->mesh_unit < REAL_EPSILON )
                throw InvalidParameter("the Mesh unit (mesh[1]) must be > 0");
            //Cytosim::log << reference() <<  " new Mesh" << '\n';
            fMesh.setUnit(prop->mesh_unit);
#else
            //throw InvalidParameter("Cytosim does not support fiber:mesh");
#endif
        }
    }
#if FIBER_HAS_GLUE
    fGlue = nullptr;
#endif
#if FIBER_HAS_FAMILY
    family_ = this;
    sister_ = nullptr;
    brother_ = nullptr;
#endif
#if NEW_FIBER_CHEW
    fChewM = 0;
    fChewP = 0;
#endif
}

#if FIBER_HAS_FAMILY
/*
 returns the position of the protofilament centerline that is used for display
 Since the protofilaments are set with an effective diameter of 27nm,
 we rescale here by 19/27 = 0.703
*/
Vector Fiber::displayPosM(real ab) const
{
    Vector I = family_->posM(ab);
    Vector O = posM(ab);
    return I + 0.703 * ( O - I );
}

#endif

Fiber::~Fiber()
{
    fHands.detachAll();
    
#if FIBER_HAS_MESH
    if ( fMesh.data() && prop->field_ptr )
        releaseMeshValues(fMesh, prop->field_ptr);
#endif

#if FIBER_HAS_GLUE
    delete(fGlue);
    fGlue = nullptr;
#endif
    if ( disp )
    {
        /*
         Note: the destructor will not be called here, which is OK
         if LineDisp is a trivial type that does not allocate resources
         */
        free(disp);
        disp = nullptr;
    }
    
    prop = nullptr;
}


real Fiber::projectPoint(Vector const& w, real & dis) const
{
    // initialize with the minus-end:
    dis = distanceSqr(w, posP(0));
    real abs = 0, len = segmentation();
    
    // try all segments
    for ( size_t s = 0; s < nbSegments(); ++s )
    {
        //check the segment:
        real d = INFINITY;
        real a = FiberSegment(this, s).projectPoint0(w, d);
        if ( len < a )
        {
            // test exact point
            real e = distanceSqr(w, posP(s+1));
            if ( e < dis ) {
                abs = abscissaPoint(s+1);
                dis = e;
            }
        }
        else if (( 0 <= a ) & ( d < dis ))
        {
            //the projection is the best found so far
            abs = abscissaPoint(s) + a;
            dis = d;
        }
    }
    
    return abs;
}

//------------------------------------------------------------------------------
#pragma mark - Modifying

/**
 Update abscissa of Hands by applying mirror image around the midpoint of the Fiber
 */
void Fiber::flipHandsPolarity()
{
    real mid = abscissaM() + abscissaP();
    Hand * ha = fHands.front();
    while ( ha )
    {
        Hand * nx = ha->next();
        ha->moveTo(mid-ha->abscissa());
        ha = nx;
    }
}


/**
 A portion of size `len` that includes the MINUS_END is removed.
 The Hands bound within the deleted portion are detached.
 */
void Fiber::cutM(real len)
{
    const real abs = abscissaM() + len;
    
    Chain::cutM(len);
    
    Hand * h = fHands.front();
    while ( h )
    {
        Hand * x = h->next();
        if ( h->abscissa() < abs )
            h->detach();
        else
            h->reinterpolate();
        h = x;
    }
}


/**
 A portion of size `len` that includes the PLUS_END is removed.
 The Hands bound within the deleted portion are detached.
 */
void Fiber::cutP(real len)
{
    const real abs = abscissaP() - len;
    
    Chain::cutP(len);
    
    Hand * h = fHands.front();
    while ( h )
    {
        Hand * x = h->next();
        if ( h->abscissa() > abs )
            h->detach();
        else
            h->reinterpolate();
        h = x;
    }
}


/**
 The Fiber is cut at point P
 - A new Fiber is created from the section [ P , PLUS_END ],
 - all FiberSite attached to this section are transferred to the new Fiber,
 - Lattice content is also transferred,
 - a pointer to the new Fiber is returned, which should be added to the Simul
 .
 @return zero, if `pti` is not an internal point
 */
Fiber* Fiber::severPoint(size_t pti)
{
    if ( pti == 0  ||  pti >= lastPoint() )
        return nullptr;
    
    const real abs = abscissaPoint(pti);

    // create a new Fiber of the same class:
    Fiber* fib = prop->newFiber();
    assert_true( fib->prop == prop );
    
    // copy the Chain part of the object:
    *(static_cast<Chain*>(fib)) = *this;
    *(static_cast<Object*>(fib)) = *this;

    // the signature on both pieces should be conserved:
    fib->signature(signature());
    fib->birthTime(birthTime());

    assert_true( fib->abscissaM() == abscissaM() );
    // remove MINUS_END portion on new piece:
    fib->truncateM(pti);
    assert_true(fib->abscissaM() == abs);
    fib->updateRange();

#if FIBER_HAS_MESH
    if ( fMesh.data() )
    {
        assert_true( fMesh.unit() == fib->fMesh.unit() );
        // transfer Lattice values located above the cut
        fib->fMesh.takeP(fMesh, fMesh.index_round(abs));
    }
#endif

    // remove PLUS_END portion on self
    truncateP(pti);
    
    // transfer Hands above point P, at same abscissa
    Hand * h = fHands.front();
    while ( h )
    {
        Hand * x = h->next();
        if ( h->abscissa() > abs )
            h->relocate(fib);
        else
            h->reinterpolate();
        h = x;
    }
    
    resetLattice();
    fib->resetLattice();
    return fib;
}


/**
The Fiber is cut at distance `abs` from its MINUS_END:
 - current Fiber is truncated to keep only the section [ MINUS_END , abs ],
 - A new Fiber is created representing the other section [ abs , PLUS_END ],
 - Hands are transferred to the new Fiber if appropriate,
 - lattice substances are also transferred,
 .
 A pointer to the new Fiber is returned (containing the PLUS_END), but this
 pointer may be zero, if `abs` was not within the valid range of abscissa.
 If a new Fiber was created, it should be added to the FiberSet.
 */
Fiber* Fiber::severM(real abs)
{
    // create a new Fiber of the same class:
    Fiber* fib = prop->newFiber();
    assert_true( fib->prop == prop );

    // copy the Chain part of the object:
    *(static_cast<Chain*>(fib)) = *this;
    *(static_cast<Object*>(fib)) = *this;

    // the signature on both pieces should be conserved:
    assert_true(fib->signature() == signature());
    assert_true(fib->birthTime() == birthTime());
    assert_small(fib->abscissaM() - abscissaM());
    assert_small(fib->abscissaP() - abscissaP());

    // remove MINUS_END portion on new piece
    fib->Chain::cutM(abs);
    // initialize Lattice on new piece
    fib->updateRange();

#if FIBER_HAS_MESH
    if ( fMesh.data() )
    {
        assert_true( fMesh.unit() == fib->fMesh.unit() );
        // transfer Lattice values located above the cut:
        fib->fMesh.takeP(fMesh, fMesh.index_round(abscissaM()+abs));
    }
#endif
    
    // remove PLUS_END portion on self
    Chain::cutP(length()-abs);
    
    assert_small(fib->abscissaM()-abscissaP());

    // transfer all Hands above cut to new piece
    // their abscissa should not change in this transfer
    abs += abscissaM();
    Hand * h = fHands.front();
    while ( h )
    {
        Hand * x = h->next();
        if ( h->abscissa() >= abs )
            h->relocate(fib);
        else
            h->reinterpolate();
        h = x;
    }

    resetLattice();
    fib->resetLattice();

    return fib;
}


Fiber* Fiber::severNow(const real abs)
{
    std::clog << "sever " << reference() << " at " << abs << '\n';
    if ( abscissaM() < abs && abs < abscissaP() )
        return severM(abs-abscissaM());
    return nullptr;
}


/**
 Perform all cuts registered in `pendingCuts`, and clear that list.
 This deletes Fibers that are shorter than FiberProp::min_length
*/
void Fiber::severNow()
{
    /**
     The std::set keeps its objects always in order of descending abscissa,
     which is essential here to avoid data loss:
     cut from high to low abscissa
     */
    for ( CutFacts const& cut : pendingCuts )
    {
        //std::clog << "cut " << cut.abs << " [ " << abscissaM() << " " << abscissaP() << " ]\n";
        if ( cut.abs - abscissaM() <= prop->min_length )
        {
            // we check the range again, since the fiber tip may have changed:
            if ( cut.abs > abscissaM() )
            {
                cutM(cut.abs-abscissaM());
                setEndStateM(cut.stateM);
            }
            /*
             since we have deleted the MINUS_END section,
             the following cuts in the list, which will be of lower abscissa,
             should not be processed.
             */
            break;
        }
        else if ( abscissaP() - cut.abs <= prop->min_length )
        {
            // we check the range again, since the fiber tip may have changed:
            if ( cut.abs < abscissaP() )
            {
                cutP(abscissaP()-cut.abs);
                setEndStateP(cut.stateP);
            }
        }
        else
        {
            Fiber * frag = severNow(cut.abs);
            
            // special case where the PLUS_END section is simply deleted
            if ( cut.stateM == STATE_BLACK )
            {
                delete(frag);
                continue;
            }

            if ( frag )
            {
                //add new fragment to simulation:
                objset()->add(frag);

                // check that ends spatially match:
                assert_small((frag->posEndM() - posEndP()).norm());
                
                try {
                    // old PLUS_END converves its state:
                    frag->setEndStateP(endStateP());
                    
                    // new ends are set as wished:
                    this->setEndStateP(cut.stateP);
                    frag->setEndStateM(cut.stateM);
                }
                catch ( Exception & e )
                {
                    e << "while cutting fiber " << reference();
                    throw;
                }
            
#ifdef LOGGING
                Cytosim::log << "severed " << reference() << " at abscissa " << cut.abs;
                Cytosim::log << "   creating " << frag->reference();
                Cytosim::log << "   position " << frag->posEndM() << '\n';
#endif
                //Cytosim::log << " severed at X = " << frag->posEndM().XX << '\n';
            }
            else
            {
                Cytosim::log << " sever abscissa " << cut.abs << " is out of range";
                Cytosim::log << " [ " << abscissaM() << "   " << abscissaP() << " ]" << '\n';
            }
        }
    }
    pendingCuts.clear();
}


/**
 returns index of first point for which ( std::cos(angle) < max_cosine ),
 or zero
 */
size_t Fiber::hasKink(const real max_cosine) const
{
    size_t end = nPoints - 2;
    for ( size_t p = 0; p < end; ++p )
    {
        if ( dot(diffPoints(p), diffPoints(p+1)) < max_cosine )
            return p+1;
    }
    return 0;
}


void Fiber::planarCut(Vector const& n, const real a, state_t stateP, state_t stateM)
{
    Array<real> cuts;
    
    /*
     The cuts should be processed in order of decreasing abscissa,
     hence we check intersections from PLUS_END to MINUS_END
    */
    for ( size_t s = nbSegments(); s >0 ; --s )
    {
        real abs = planarIntersect(s-1, n, a);
        if ( 0 <= abs  &&  abs < 1 )
            cuts.push_back(abscissaPoint(s-1+abs));
    }
    
    for ( real abs : cuts )
    {
        Fiber * fib = severNow(abs);
        if ( fib )
        {
            // old PLUS_END converves its state:
            fib->setEndStateP(endStateP());
            // dynamic of new ends are set as usual:
            setEndStateP(stateP);
            fib->setEndStateM(stateM);
            objset()->add(fib);
        }
    }
}


/**
 The given `fib` is added past the PLUS_END of `*this`,
 Hands bound to `fib` are transferred to *this.
 The dynamic state of the PLUS_END is also transferred.
 `fib` is enventually deleted
*/
void Fiber::join(Fiber * fib)
{
    assert_true( fib );
    // the two fibers should be of the same class:
    assert_true( prop == fib->prop );
    
    // shift in abscissa must be calculated before joining
    real shift = abscissaP() - fib->abscissaM();

    // join backbones
    Chain::join(fib);
    // extend Lattice range
    updateFiber();

    //transfer dynamic state of PLUS_END:
    setEndStateP(fib->endStateP());

#if FIBER_HAS_MESH
    if ( fMesh.data() )
    {
        assert_true( fMesh.unit() == fib->fMesh.unit() );
        // transfer Lattice values from other fiber
        fMesh.takeP(fib->fMesh, fMesh.indexM());
    }
#endif

    // transfer all Hands
    Hand * ha = fib->fHands.front();
    while ( ha )
    {
        Hand * nx = ha->next();
        ha->relocate(this, ha->abscissa()+shift);
        ha = nx;
    }
    delete(fib);

    resetLattice();
}


//------------------------------------------------------------------------------
#pragma mark - Mobility


#if NEW_ANISOTROPIC_FIBER_DRAG
    constexpr real DRAG = 4;
#else
    constexpr real DRAG = 3;
#endif


/**
 From "Random Walks in Biology" by HC. Berg, Princeton University Press,
 drag coefficients for an ellipsoid are,

     drag_transverse = 2*drag_parallel = 4*PI*length*visc / log(length/radius)

 We should average the mobility coefficients:  speed = mu * f
     mu_X = mu_parallel   = 2 * mu
     mu_Y = mu_transverse = mu
     mu_Z = mu_transverse = mu
 Hence:
     mu_averaged = ( mu + mu + 2*mu ) / 3 = 4/3 * mu.
 drag_averaged = 3*PI*length*viscosity / log(length/radius)

APPROXIMATE FORMULA FOR ELLIPSOIDAL PARTICLE
 > Clift R, Grace JR, Weber ME. Bubbles, drops, and particles
 > Courier Corporation; 2005.

     aspect = length / diameter;
     drag = 3.0 * M_PI * viscosity * diameter * ( 3 + 2 * length/diameter ) / 5.0;

 */
real Fiber::dragCoefficientEllipsoid(const real len, FiberProp const* prop)
{
    // Stokes' law for a sphere having diameter of filament:
    assert_true( prop->drag_radius > 0 );
    real drag = 6 * prop->drag_radius;

    assert_true( len > REAL_EPSILON );
    assert_true( prop->drag_length > REAL_EPSILON );
    // hydrodynamic cut-off on length:
    const real lenc = std::max(std::min(len, prop->drag_length), prop->drag_radius);

    // drag of 3D ellipsoid
    const real drag_ellipsoid = DRAG * len / std::log( lenc / prop->drag_radius );
    
    // length below which ellipsoid formula is not valid:
    const real min_len = std::exp( 1 + std::log(prop->drag_radius) );
    
    if ( len > min_len )
    {
        // use largest drag coefficient
        drag = std::max(drag, drag_ellipsoid);
    }

    assert_true( prop->viscosity > 0 );
    return M_PI * prop->viscosity * drag;
}



/** 
 dragCoefficientCylinder() calculates the mobility for the entire fiber,
 considering that the cylinder is straight and moving in an infinite fluid.
 fiber:drag_length is a hydrodynamic cutoff that makes the
 drag coefficient proportional to length beyond the cutoff.
 
 The drag is determined by the viscosity and the length and diameter of the
 filament. The aspect ratio is defined by:
 
     shape = length / diameter;

 The formula for a cylinder were calculated numerically in:
 > Tirado and de la Torre. J. Chem. Phys 71(6) 1979
 > http://doi.org/10.1063/1.438613
 > Page 2584, Table 1, last column, last line for infinite aspect ratio

 The translational drag coefficient is averaged over all possible configurations:
 
       drag_cylinder = 3 * PI * viscosity * length / ( log(shape) + 0.312 );
 
 If the length is shorter than the diameter, the formula above fails and may
 even give negative result. Hence we also calculate the drag of a sphere with
 the same radius as the cylinder:

       drag_sphere = 6 * PI * viscosity * radius

 We use the maximum value between 'drag_sphere' and 'drag_cylinder'.
 */
/*
 Ct =  0.312 + 0.565/shape - 0.100/(shape*shape);

 The rotational diffusion coefficient is given by:
 > Tirado and de la Torre. J. Chem. Phys 73(4) 1980
 
     Cr = -0.662 + 0.917/aspect - 0.050/(shape*shape);
     drag_rotation = 1/3*M_PI*viscosity*length^3 / ( log(shape) + Cr )

 */
real Fiber::dragCoefficientCylinder(const real len, FiberProp const* prop)
{
    // Stokes' law for a sphere having diameter of filament:
    assert_true( prop->drag_radius > 0 );
    real drag = 6 * prop->drag_radius;

    assert_true( len > REAL_EPSILON );
    assert_true( prop->drag_length > REAL_EPSILON );
    // hydrodynamic cut-off on length:
    const real lenc = std::max(std::min(len, prop->drag_length), prop->drag_radius);
    
    /// drag of a cylinder:
    const real drag_cylinder = DRAG * len / ( std::log(0.5*lenc/prop->drag_radius) + 0.32 );

    // length below which the formula is not valid anymore ( ~ 3.94 * radius )
    // this corresponds to the minimun value of the formula above
    const real min_len = 2 * prop->drag_radius * std::exp(1.0-0.32);
    
    if ( len > min_len )
    {
        // use largest drag coefficient
        drag = std::max(drag, drag_cylinder);
    }
    
    assert_true( prop->viscosity > 0 );
    return M_PI * prop->viscosity * drag;
}


/**
 dragCoefficientSurface() uses a formula calculated by F. Gittes in:
 > Hunt et al. Biophysical Journal (1994) v 67 pp 766-781  
 > http://dx.doi.org/10.1016/S0006-3495(94)80537-5
 
 It applies to a cylinder moving parallel to its axis and near an immobile surface:

       drag_per_unit_length = 2 &pi &eta / acosh(h/r)
 
 With:
 - r = cylinder radius,
 - h = distance between cylinder bottom and surface,
 - &eta = viscosity of the fluid.
 
 If the cylinder is exactly touching the surface, `h=0` and the drag coefficient is infinite.
 
 The drag coefficient for motion perpendicular to the cylinder axis would be twice higher,
 but for gliding assays, the parallel drag coefficient is the appropriate choice.  
 
 Note that this is usually equivalent to the approximate formula:

       drag_per_unit_length = 2 &pi &eta / log(2*h/r)

 because

       acosh(x) = ln[ x + std::sqrt(x^2-1)) ] ~ ln[2x] if x >> 1

 Hunt et al. also credit this reference for the formula:
 > The slow motion of a cylinder next to a plane wall.
 > Jeffrey, D.J. & Onishi, Y. (1981) Quant. J. Mech. Appl. Math. 34, 129-137.
*/
real Fiber::dragCoefficientSurface(const real len, FiberProp const* prop)
{
    if ( prop->drag_gap <= 0 )
        throw InvalidParameter("fiber:drag_model[1] (height above surface) must set and > 0!");
    
    // use the higher drag: perpendicular to the cylinder (factor 2)
    real drag = 2 * len / acosh( 1 + prop->drag_gap/prop->drag_radius );

    return M_PI * prop->viscosity * drag;
}


/**
 Calculate drag coefficient by calling one of these functions:

        dragCoefficientEllipsoid();
        dragCoefficientCylinder();
        dragCoefficientSurface();

 */
void Fiber::setDragCoefficient()
{
    const real len = length();
    assert_true( len > 0 );
    
    real drag = 0;
    
    if ( prop->drag_model )
    {
        drag = dragCoefficientSurface(len, prop);
#if ( 0 )
        real d = dragCoefficientCylinder(len, prop);
        Cytosim::log << "Drag of Fiber near a planar surface amplified by " << drag/d << '\n';
#endif
    }
    else
        drag = dragCoefficientCylinder(len, prop);

    assert_true( drag > 0 );
    // distribute drag equally to all points, to set point's mobility
    iPointMobility = nPoints / drag;
    
#if ( 0 )
    std::ostream& os = std::cerr; //Cytosim::log;
    os << "Fiber " << std::setw(16) << prop->name() << " length " << std::setw(7) << length();
    os << " drag " << std::setw(9) << drag << " point_mobility " << iPointMobility << '\n';
#endif
}


void Fiber::prepareMecable()
{
    setDragCoefficient();
    storeDirections();
#if NEW_UNCONSTRAINED_LENGTH
    constrainLength(prop->constrain_length);
    if ( !unconstrainLength )
#endif
#if NEW_ANISOTROPIC_FIBER_DRAG
    makeProjectionAnisotropic();
#else
    makeProjection();
#endif
    //printProjection(std::clog);

    assert_true( iPointMobility >= 0 );
    
    // the scaling of the bending elasticity depends on the length of the segments
    iRigidity = prop->rigidity / cube(segmentation());
#if NEW_FIBER_LOOP
    iRigidityLoop = prop->loop;
#endif
#if ( 0 )
    real energy = bendingEnergy();
    real euler = M_PI * M_PI * prop->rigidity / ( length() * length() );
    Cytosim::log << "Euler buckling = " << euler << "    ";
    Cytosim::log << "Bending energy = " << energy << '\n';
#endif
}


//------------------------------------------------------------------------------

/** Add interactions between vertices and Space */
void Fiber::setConfinement(Meca& meca, Confinement mode, Space const* spc, real stiff) const
{
    switch ( mode )
    {
        case CONFINE_INSIDE:
            for ( size_t i = 0; i < nPoints; ++i )
            {
                Vector pos = posP(i);
                if ( spc->outside(pos) )
                    spc->setConfinement(pos, Mecapoint(this, i), meca, stiff);
            }
            break;
            
        case CONFINE_OUTSIDE:
            for ( size_t i = 0; i < nPoints; ++i )
            {
                Vector pos = posP(i);
                if ( spc->inside(pos) )
                    spc->setConfinement(pos, Mecapoint(this, i), meca, stiff);
            }
            break;
            
        case CONFINE_ON:
            for ( size_t i = 0; i < nPoints; ++i )
                spc->setConfinement(posP(i), Mecapoint(this, i), meca, stiff);
            break;
            
        case CONFINE_MINUS_END:
            spc->setConfinement(posP(0), Mecapoint(this, 0), meca, stiff);
            break;
            
        case CONFINE_PLUS_END:
        {
            const size_t L = lastPoint();
            spc->setConfinement(posP(L), Mecapoint(this, L), meca, stiff);
        } break;
            
        case CONFINE_BOTH_ENDS:
        {
            const size_t L = lastPoint();
            spc->setConfinement(posP(0), Mecapoint(this, 0), meca, stiff);
            spc->setConfinement(posP(L), Mecapoint(this, L), meca, stiff);
        } break;
            
        case CONFINE_PLUS_OUT:
        {
            const size_t L = lastPoint();
            Vector pos = posP(L);
            if ( spc->inside(pos) )
                spc->setConfinement(pos, Mecapoint(this, L), meca, stiff);
        } break;
#if NEW_FIBER_CONFINE_RANGE
        case CONFINE_RANGE:
        {
            // we use here the MINUS_END as a reference... which maybe problematic
            size_t S = clampedIndexM(prop->confine_range[0]);
            size_t E = clampedIndexM(prop->confine_range[1]);
            for ( size_t i = S; i < E; ++i )
                spc->setConfinement(posP(i), Mecapoint(this, i), meca, stiff);
        } break;
#endif
        default:
            throw InvalidParameter("Invalid fiber:confine");
    }
}


void Fiber::setInteractions(Meca& meca) const
{
#if OLD_SQUEEZE_FORCE
    if ( prop->squeeze == 1 )
    {
        // squeezing force in the YZ-plane:
        const real f = prop->squeeze_force;
        const real r = prop->squeeze_range;
        for ( size_t i = 0; i < nPoints; ++i )
        {
            Vector P = posP(i);
            if ( P.normYZ() < r )
                meca.addLineClamp(Mecapoint(this, i), Vector(P.XX, 0, 0), Vector(1,0,0), f/r);
            else {
                // forces is capped to a maximum magnitude 'f':
#if ( DIM == 3 )
                Vector n = Vector(0, -P.YY, -P.ZZ).normalized(f);
                meca.addForce(Mecapoint(this, i), n);
#elif ( DIM == 2 )
                Vector n(0, std::copysign(f, -P.YY), 0);
                meca.addForce(Mecapoint(this, i), n);
#endif
            }
        }
    }
#endif
#if NEW_END_FORCE
    switch( prop->end_force_mode )
    {
        case MINUS_END:
            meca.addForce(Mecapoint(this, 0), prop->end_force);
            break;
        case PLUS_END:
            meca.addForce(Mecapoint(this, lastPoint()), prop->end_force);
            break;
        case BOTH_ENDS:
            meca.addForce(Mecapoint(this, 0), prop->end_force);
            meca.addForce(Mecapoint(this, lastPoint()), prop->end_force);
            break;
        case CENTER:
            meca.addForce(interpolateCenter(), prop->end_force);
            break;
        case ORIGIN: //this adds a bending torque where the sum of force is zero
            meca.addForce(Mecapoint(this, 0), prop->end_force);
            meca.addForce(Mecapoint(this, lastPoint()), prop->end_force);
            meca.addForce(interpolateCenter(), -2.0 * prop->end_force);
        default:
        break;
    }
#endif

#if NEW_COLINEAR_FORCE
    /*
     Add a length-dependent force acting parallel to the filament.
     A force proportional to the length of the segments is distributed
     to the vertices.
     */
    if ( prop->colinear_force )
    {
        real s = 0.5 * prop->colinear_force * segmentation();
        for ( size_t i = 0; i < nbSegments(); ++i )
        {
            Vector f = s * dirSegment(i);
            meca.addForce(Mecapoint(this, i  ), f);
            meca.addForce(Mecapoint(this, i+1), f);
        }
    }
#endif
    
    if ( prop->confine != CONFINE_OFF )
        setConfinement(meca, prop->confine, prop->confine_space_ptr, prop->confine_stiffness);
    
#if NEW_FIBER_CONFINE2
    /// add another confinement force
    if ( prop->confine2 != CONFINE_OFF )
        setConfinement(meca, prop->confine2, prop->confine2_space_ptr, prop->confine2_stiffness);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Attached Hands

size_t Fiber::nbHandsInRange(real i, real s, const FiberEnd ref) const
{
    // Convert to absolute abscissa:
    i = abscissaFrom(i, ref);
    s = abscissaFrom(s, ref);
    if ( s < i )
        std::swap(i, s);
    return fHands.countInRange(i, s);
}


size_t Fiber::nbHandsNearEnd(const real len, const FiberEnd ref) const
{
    real i = -INFINITY, s = INFINITY;
    
    if ( ref == PLUS_END )
        i = abscissaP() - len;
    else if ( ref == MINUS_END )
        s = abscissaM() + len;
    else
        throw("invalid argument value to nbHandsNearEnd()");
        
    return fHands.countInRange(i, s);
}

//------------------------------------------------------------------------------
#pragma mark - Dynamic ends

state_t Fiber::endState(FiberEnd end) const
{
    if ( end == PLUS_END )
        return endStateP();
    if ( end == MINUS_END )
        return endStateM();
    ABORT_NOW("invalid argument value");
    return 0;
}


void Fiber::setEndState(const FiberEnd end, const state_t s)
{
    if ( end == PLUS_END )
        setEndStateP(s);
    else if ( end == MINUS_END )
        setEndStateM(s);
}


void Fiber::updateRange()
{
#if FIBER_HAS_LATTICE
    // this will allocate the Lattice's site to cover the range of Abscissa:
    if ( fLattice.ready() )
    {
        // this will reallocate the Lattice
        fLattice.setRange(abscissaM(), abscissaP());
    }
#endif
#if FIBER_HAS_MESH
    // this will allocate the Lattice's site to cover the range of Abscissa:
    if ( fMesh.ready() )
    {
        Field* field = prop->field_ptr;
        fMesh.setRange(abscissaM(), abscissaP());
        
        if ( field )
        {
            real sumM;
            // release Lattice substance located outside the valid abscissa range
            fMesh.collectM(sumM);
            field->cell(posEndM()) += sumM;
            //Cytosim::log << " Fiber::MINUS_END releases " << sumM << '\n';
            
            real sumP;
            fMesh.collectP(sumP);
            field->cell(posEndP()) += sumP;
            //Cytosim::log << " Fiber::PLUS_END releases " << sumP << '\n';
        }
    }
#endif
}


/**
 Update all bound Hands, allowing them to detach if they are out of range
 */
void Fiber::updateHands()
{
    real M = abscissaM();
    real P = abscissaP();
    Hand * h = fHands.front();
    while ( h )
    {
        Hand * x = h->next();
        assert_true(h->fiber()==this);
        // this is equivalent to h->reinterpolate():
        h->hTerp = interpolateM(h->abscissa() - M);
        // must iterate ahead, because `checkFiberRange` may lead to detachment:
        h->checkFiberRange(M, P);
        h = x;
    }
}

/**
 Assuming that the length has changed, or that the abscissa of the ends have changed,
 this updates the segmentation of the fiber if needed, the position of the Hands,
 and the boundaries of the Lattice if present.
 */
void Fiber::updateFiber()
{
    needUpdate = false;
#if ( 0 )
    Cytosim::log << reference() << " update [ "  << std::setw(9) << std::left << abscissaM();
    Cytosim::log << " "  << std::setw(9) << std::left << abscissaP() << " ]" << '\n';
#endif
    
    updateRange();
    updateHands();
}

//------------------------------------------------------------------------------
#pragma mark - Lattice


void Fiber::printLattice(std::ostream& os) const
{
#if FIBER_HAS_LATTICE
    FiberLattice const& lat = fLattice;
    using std::setw;
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    os << "Lattice for " << reference() << ":\n";
    os << "    inf  " << inf << "  " << abscissaM() << "\n";
    os << "    sup  " << sup << "  " << abscissaP() << "\n";
    for ( auto h = inf; h < sup; ++h )
        os << setw(8) << h << "  " << setw(10) << lat.abscissa(h) << setw(10) << lat.data(h) << "\n";
    os << "\n";
#endif
}


void Fiber::infoLattice(size_t& cnt, size_t& vac, real& sum, real& mn, real& mx) const
{
#if FIBER_HAS_LATTICE
    FiberLattice const& lat = fLattice;
    if ( lat.data() )
    {
        const auto sup = lat.indexP();
        for ( auto i = lat.indexM(); i <= sup; ++i )
        {
            ++cnt;
            vac += ( lat.data(i) == 0 );
            sum += lat.data(i);
            real x = lat.data(i);
            mn = std::min(mn, x);
            mx = std::max(mx, x);
        }
    }
#endif
}


VisibleLattice const* Fiber::visibleLattice() const
{
#if FIBER_HAS_MESH
    if ( fMesh.data() )
        return &fMesh;
#elif FIBER_HAS_LATTICE
    if ( fLattice.data() )
        return &fLattice;
#endif
    return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Mesh

#if FIBER_HAS_MESH

/**
 */
void Fiber::setMeshValues(Lattice<real>& lat, real density) const
{
    const real uni = lat.unit();
    assert_true( uni > 0 );
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();
    
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        assert_true( lat.abscissa(inf+1) >= abscissaP() );
        site[inf] = density * ( abscissaP() - abscissaM() );
    }
    else
    {
        // the terminal site may be truncated
        site[inf] = density * ( lat.abscissa(inf+1) - abscissaM() );
        
        for ( auto h = inf+1; h < sup; ++h )
        site[h] = density * uni;
        
        // the terminal site may be truncated
        site[sup] = density * ( abscissaP() - lat.abscissa(sup) );
    }
}

/**
Update all Mesh values according to:

    site[i] <- cst + fac * site[i]

 */
void Fiber::evolveMeshValues(Lattice<real>& lat, real cst, real fac) const
{
    assert_false(lat.bad());
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();

    //std::clog << "evolve " << inf << " " << sup << "\n";
    for ( auto h = inf; h <= sup; ++h )
        site[h] = cst + fac * site[h];
}


void Fiber::bindMesh(Lattice<real>& lat, Field * fld, real bind_rate) const
{
    assert_false(lat.bad());
    // we want roughly one point per cell:
    const real spread = fld->cellWidth();
    
    // each point represents a Fiber chunk of length 'spread':
    const real rate = bind_rate * spread / fld->cellVolume();
    
    // fraction of the cell content that will bind in one time step:
    const real frac = -std::expm1( -rate * time_step(simul()) );
    
    real abs = spread * RNG.exponential();
    const real len = length();
    
    // stochastic sampling with distance 'spread' along the Fiber:
    while ( abs < len )
    {
        Vector pos = posM(abs);
        real& cell = fld->cell(pos);
        assert_true( cell >= 0 );
        
        // amount to be transferred:
        real flux = cell * frac;
        
        cell -= flux;
        lat.cell(abs+abscissaM()) += flux;
        
        abs += spread * RNG.exponential();
    }
}


/**
 Release a fraction 'frac' of the Mesh substance into the Field.
 The subtance in each Lattice site is released in a cell
 corresponding to a random position within this site.
 The factor `frac` must be between 0 and 1.
 */
void Fiber::equilibrateMesh(Lattice<real>& lat, Field * fld, real on, real off) const
{
    const real uni = lat.unit();
    assert_true( uni > 0 );
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();

    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        assert_true( lat.abscissa(inf+1) >= abscissaP() );
        real & cell = fld->cell(posM(RNG.preal()*length()));
        real flux = on * cell - off * site[inf];
        cell      -= flux;
        site[inf] += flux;
    }
    else
    {
        // the terminal site may be truncated
        real a = RNG.preal() * ( uni*(inf+1) - abscissaM() );
        real & cell = fld->cell(posM(a));
        real flux = on * cell - off * site[inf];
        cell      -= flux;
        site[inf] += flux;
        
        for ( auto h = inf+1; h < sup; ++h )
        {
            // we select a random position along each site and find corresponding cell:
            cell = fld->cell(pos(uni*(RNG.preal()+h)));
            flux = on * cell - off * site[h];
            cell    -= flux;
            site[h] += flux;
        }
        
        // the terminal site may be truncated
        a = uni*sup + RNG.preal() * ( abscissaP() - uni*sup );
        cell = fld->cell(pos(a));
        flux = on * cell - off * site[sup];
        cell      -= flux;
        site[sup] += flux;
    }
}


void Fiber::fluxMesh(Lattice<real>& lat, Field * fld, real speed) const
{
    assert_false(lat.bad());
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();

    const real fac = speed * time_step(simul()) / lat.unit();
    
    if ( abs_real(fac) > 1 )
        throw InvalidParameter("mesh_flux_speed * time_step / lattice_unit is too high");

    if ( fac < 0 )
    {
        real s = site[inf];
        
        for ( auto h = inf; h < sup; ++h )
            site[h] -= fac * ( site[h+1] - site[h] );
        
        fld->cell(posEndM()) -= fac * s;
        site[sup] += fac * site[sup];
    }
    else
    {
        real s = site[sup];
        
        for ( auto h = sup; h > inf; --h )
            site[h] -= fac * ( site[h] - site[h-1] );
        
        fld->cell(posEndP()) += fac * s;
        site[inf] -= fac * site[inf];
    }
}


/**
 Release all Mesh substance into the Field.
 The subtance in each Lattice site is released in a cell
 corresponding to a random position within this site.
 */
void Fiber::releaseMeshValues(Lattice<real>& lat, Field * fld) const
{
    const real uni = lat.unit();
    assert_true( uni > 0 );
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();
    
    //@todo Handle the terminal site differently since they are truncated
    for ( auto h = inf; h <= sup; ++h )
    {
        fld->cell(pos(uni*(RNG.preal()+h))) += site[h];
        site[h] = 0;
    }
}


void Fiber::cutFiberMesh(Lattice<real>& lat)
{
    const real uni = lat.unit();
    assert_true( uni > 0 );
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();

    const real fac = 1.0 / time_step(simul());
    
    assert_true( inf >= lat.inf() );
    assert_true( sup <  lat.sup() );
    assert_true( lat.abscissa(inf) <= abscissaM() );
    assert_true( lat.abscissa(sup+1) >= abscissaP() );
    
    auto h = inf;
    real val = fac * RNG.exponential();
    real ai  = abscissaM();
    real as  = lat.abscissa(h+1);
    if ( as > abscissaP() )
        as = abscissaP();
    
    while ( h <= sup )
    {
        assert_true( site[h] >= 0 );
        val -= site[h];
        //assert_true( ai >= abscissaM() );
        //assert_true( as <= abscissaP() );
        
        while ( val < 0 )
        {
            /*
             Since val < 0 and 0 <= val+site[h],
             then        -val/site[h] <= 1
             hence   0 < -val/site[h] < 1
             */
            real abs = ai - ( as - ai ) * val / site[h];
            
            assert_true( abs >= ai - REAL_EPSILON );
            assert_true( abs <= as + REAL_EPSILON );
            
            sever(abs, STATE_RED, STATE_GREEN);
            val += fac * RNG.exponential();
        }
        
        ai = as;
        if ( ++h == sup )
            as = abscissaP();
        else
            as += uni;
    }
}

#endif

void Fiber::infoMesh(real& len, size_t& cnt, real& sm, real& mn, real& mx, bool density) const
{
#if FIBER_HAS_MESH
    Lattice<real> const& lat = fMesh;
    if ( lat.data() )
    {
        len += length();
        const real scale = ( density ? 1.0/lat.unit() : 1.0 );
        const auto sup = lat.indexP();
        for ( auto i = lat.indexM(); i <= sup; ++i )
        {
            ++cnt;
            sm += lat.data(i);
            real x = lat.data(i) * scale;
            mn = std::min(mn, x);
            mx = std::max(mx, x);
        }
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Glue

/**
 setGlue1 keeps the Single attached as long as the Fiber tip is outside the Space.
 The Single's hand is managed to always remain at the tip of the Fiber.
 The Single detaches immediately if the Fiber tip is back inside.
 This generates mostly a pushing force from the cortex
 */
void Fiber::setGlue1(Single* glue, const FiberEnd end, Space const* spc)
{
    assert_true(spc);
    if ( spc->inside(posEnd(end)) )
    {
        //detach immediately if the tip is inside the Space
        if ( glue->attached() )
            glue->detach();
    }
    else
    {
        if ( glue->attached() )
        {
            //always keep tracking the tip:
            glue->moveToEnd(end);
        }
        else {
            //reposition the Single base:
            glue->setPosition(spc->project(posEnd(end)));
            //attach to the MT-tip:
            glue->attachEnd(this, end);
        }
    }
}


/**
 setGlue2 keeps the Single always at the tip of the Fiber.
 The Single's hand detaches only spontaneously.
 This creates both pulling and pushing force from the cortex
 */
void Fiber::setGlue2(Single* glue, const FiberEnd end, Space const* spc)
{
    assert_true(spc);
    if ( glue->attached() )
    {
        //keep tracking the tip of the fiber while attached
        glue->moveToEnd(end);
    }
    else
    {
        // Attach a new grafted if MT-tip is outside and growing:
        if ( isGrowing(end) && spc->outside(posEnd(end)) )
        {
            //reposition the Single base:
            glue->setPosition(spc->project(posEnd(end)));
            //attach to the MT-tip:
            glue->attachEnd(this, end);
        }
    }
}


/**
 setGlue3 keeps the Single where the Fiber first crosses the Space's edge.
 This makes an anchor point exactly at the edge, and only controls attachment.
 The Single's Hand behaves and detaches normally.
 */
void Fiber::setGlue3(Single* glue, Space const* spc)
{    
    assert_true(spc);
    /*
     If the glue is not already attached, we first check if the fiber intersects
     the edge of the Space:
     */
    if ( ! glue->attached() )
    {
        bool in = spc->inside(posEndM());
        
        if ( in == spc->inside(posEndP()) )
            return;
        
        // find a vertex that is on the other side of the Space edge:
        for ( size_t i = 1; i < nPoints; ++i )
        {
            if ( spc->inside(posP(i)) != in )
            {
                // the abscissa is interpolated using the distances of P1 and P2 to the edge
                real d1 = spc->distanceToEdge(posP(i-1));
                real d2 = spc->distanceToEdge(posP(i));
                if ( d1 + d2 > REAL_EPSILON )
                {
                    /* we find the abscissa corresponding to the intersection,
                     assuming that the edge is locally straight */
                    FiberSite sit(this, abscissaPoint(i-1+d1/(d2+d1)));
                    glue->hand()->attach(sit);
                    glue->setPosition(glue->posHand());
                    break;
                }
            }
        }
    }
}


/**
 setGlueG keeps the Single attached to the fiber tip if this tip is Growing.
 */
void Fiber::setGlueG(Single* glue, FiberEnd end)
{
    if ( glue->attached() )
    {
        if ( isShrinking(end) )
            glue->detach();
        else
            glue->moveToEnd(end);
    }
    else if ( isGrowing(end) )
    {
        glue->attachEnd(this, end);
        glue->unbase();
    }
}


/**
 setGlueE keeps the Single at the end of the Fiber until it detaches spontaneously,
 and rebind the Single when the Fiber tip resumes growth.
 */
void Fiber::setGlueE(Single* glue, FiberEnd end)
{
    if ( glue->attached() )
    {
        glue->moveToEnd(end);
    }
    else if ( isGrowing(end) )
    {
        glue->attachEnd(this, end);
        glue->unbase();
    }
}


/**
 Search for a glue in the list of bound HandSingle
 this is useful when a simulation is restarted from file
 */
void Fiber::makeGlue(Single*& glue)
{
    SingleSet& set = simul().singles;

    for ( Single * s = set.firstA(); s; s=s->next() )
    {
        if ( s->hand()->fiber() == this  &&  s->mark() == identity() )
        {
            glue = s;
            //Cytosim::log << "found Fiber:glue for " << reference() << '\n';
            return;
        }
    }
    
    // create the Single if needed
    if ( !glue )
    {
        glue = prop->glue_prop->newSingle();
        glue->mark(identity());
        set.add(glue);
    }
}


/**
 setGlue() creates a Single associated with this Fiber, used to implement
 specific effects or interactions associated with the tip of the Fiber:
 - glue type 1, 2 and 3 anchor themselves at the edge of the Space
 - glue type 4 and 5 anchor themselves on Solids
 .
 This was used for example to implement cortical anchors for astral microtubules.
*/
void Fiber::setGlue(Single*& glue, const FiberEnd end, int mode)
{
    if ( !glue )
        makeGlue(glue);
    
    switch( mode )
    {
        case 1: setGlue1(glue, end, prop->confine_space_ptr); break;
        case 2: setGlue2(glue, end, prop->confine_space_ptr); break;
        case 3: setGlue3(glue, prop->confine_space_ptr); break;
        case 4: setGlueG(glue, end); break;
        case 5: setGlueE(glue, end); break;
        default: throw InvalidParameter("invalid value of fiber:glue");
    }
    
#if ( 1 )
    // we keep the Single linked in the simulation only if it is attached:
    if ( glue->attached() )
    {
        if ( !glue->linked() )
        {
            glue->objset(&simul().singles);
            simul().singles.link(glue);
        }
    }
    else if ( glue->linked() )
    {
        glue->objset(nullptr);
        simul().singles.unlink(glue);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - I/O

void Fiber::write(Outputter& out) const
{
#if !NEW_COMPACT_STORAGE
    // normal way
    writeHeader(out, tag());
    Chain::write(out);
#else
    // compact format created on 23/06/2021
    writeHeader(out, TAG_ALT);
    Chain::writeAngles(out);
#endif
    
#if FIBER_HAS_LATTICE
    /*
     We can save the occupancy Lattice here, but this is not necessary
     as it can be recalculated on the fly, so we save space by skipping
     */
    if ( prop->save_lattice && fLattice.data() )
    {
        writeHeader(out, TAG_LATTICE);
        // fLattice.write(out);
        // only write information corresponding to actual Fiber abscissa range:
        fLattice.write(out, fLattice.indexM(), fLattice.indexP()+1);
    }
#endif
#if FIBER_HAS_MESH
    if ( fMesh.data() )
    {
        writeHeader(out, TAG_FIBMESH);
        // fMesh.write(out);
        // only write information corresponding to actual Fiber abscissa range:
        fMesh.write(out, fMesh.indexM(), fMesh.indexP()+1);
    }
#endif
}


void Fiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
#if BACKWARD_COMPATIBILITY < 57
    // before 23/06/2021, TAG_LATTICE was 'l' and TAG_FIBMESH was 'L'
    if ( in.formatID() < 57 )
    {
        if ( tag == 'l' ) tag = TAG_LATTICE;
        if ( tag == 'L' ) tag = TAG_FIBMESH;
    }
#endif
    //std::clog << " Fiber::read(" << tag << ")\n";
    if ( tag == TAG )
    {
        Chain::read(in, sim, tag);
#if FIBER_HAS_LATTICE
        if ( fLattice.data() )
            fLattice.setRange(abscissaM(), abscissaP());
#endif
#if FIBER_HAS_MESH
        if ( fMesh.data() )
            fMesh.setRange(abscissaM(), abscissaP());
#endif
#if BACKWARD_COMPATIBILITY < 50
        if ( in.formatID() > 47 && in.formatID() < 50 ) // 4.7.2018 added birthTime
            birthTime(in.readFloat());
#endif
        if ( length() + 128*FLT_EPSILON < prop->min_length )
        {
            Cytosim::warn << "fiber:length < min_length ( " << length() << " < " << prop->min_length << " )\n";
        }
#if FIBER_HAS_GLUE
        fGlue = nullptr;
#endif
    }
    else if ( tag == TAG_ALT )
    {
        Chain::readAngles(in, sim, tag);
        updateRange();
#if FIBER_HAS_GLUE
        fGlue = nullptr;
#endif
    }
    else if ( tag == TAG_LATTICE )
    {
#if FIBER_HAS_LATTICE
        if ( fLattice.ready() )
            fLattice.setRange(abscissaM(), abscissaP());
        fLattice.read(in);
#else
        FiberLattice dummy;
        dummy.read(in);
        // store unit, to get digits at the right abscissa
        const_cast<FiberProp*>(prop)->lattice_unit = dummy.unit();
#endif
    }
    else if ( tag == TAG_FIBMESH )
    {
#if FIBER_HAS_MESH
        if ( fMesh.ready() )
            fMesh.setRange(abscissaM(), abscissaP());
        fMesh.read(in);
        //real S = fMesh.sum();
        //if ( S ) std::cerr << reference() << " sum(mesh) = " << S << "\n";
#else
        Lattice<real> dummy;
        dummy.read(in);
#endif
    }
#if BACKWARD_COMPATIBILITY < 100
    else if ( tag == TAG_DYNAMIC )
    {
        static bool virgin = true;
        // that is for Fiber class we do not know...
        if ( virgin )
        {
            virgin = false;
            std::cerr << "INCIDENT: skipping end states for `" << prop->name() << "'\n";
        }
        in.readUInt16();
        in.readUInt16();
        // try to recover the next record in file:
        if ( !isalpha(in.peek()) )
        {
            in.readUInt16();
            in.readUInt16();
        }
    }
#endif
    else
        Cytosim::log << "unknown Fiber TAG `" << (char)tag << "'\n";
}


