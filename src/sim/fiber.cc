// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber.h"
#include "field.h"
#include "messages.h"
#include "glossary.h"
#include "iowrapper.h"
#include "fiber_segment.h"
#include "fiber_prop.h"
#include "simul_prop.h"
#include "meca.h"
#include "hand.h"
#include "sim.h"



void Fiber::step()
{
    assert_small( length1() - length() );
#if FIBER_HAS_GLUE
    //add single that act like glue
    if ( prop->glue )
    {
        setGlue(frGlue, PLUS_END, prop->confine_space_ptr);
    }
#endif
#if NEW_FIBER_CHEW
    if ( frChewP > 0 )
    {
        if ( frChewP > prop->max_chewing_speed_dt )
            frChewP = prop->max_chewing_speed_dt;
        
        if ( length() - frChewM - frChewP < prop->min_length )
        {
            delete(this);
            return;
        }
        growP(-frChewP);
        frChewP = 0;
    }
    if ( frChewM > 0 )
    {
        if ( frChewM > prop->max_chewing_speed_dt )
            frChewM = prop->max_chewing_speed_dt;

        if ( length() - frChewM < prop->min_length )
        {
            delete(this);
            return;
        }
        growM(-frChewM);
        frChewM = 0;
    }
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
#if FIBER_HAS_FAMILY
        handshake();
#endif
    }
    
#if FIBER_HAS_MESH
    
    //std::clog << reference() << " mesh " << std::fixed << frMesh->data(0) << std::endl;
    //std::clog << reference() << " mesh sum = " << frMesh->sum() << "     ";
    
    if ( prop->mesh_binding_rate > 0 || prop->mesh_unbinding_rate > 0 )
    {
        real on  = prop->mesh_binding_rate * simul().prop->time_step;
        real off = -std::expm1( -prop->mesh_unbinding_rate * simul().prop->time_step );
        equilibrateMesh(frMesh, prop->field_ptr, on, off);
    }
    
    if ( prop->mesh_flux_speed != 0 )
        fluxMesh(frMesh, prop->field_ptr, prop->mesh_flux_speed);
    
    if ( prop->mesh_cut_fiber )
        cutFiberMesh(frMesh);
    
    if ( prop->mesh_aging_rate > 0 )
    {
        /*
         This implements an evolution toward equilibrium,
         starting from a value of 0, and reached for a value of 1,
         with a time-scale given by 'mesh_aging_rate'.
         */
        real cst = prop->mesh_aging_rate * simul().prop->time_step;
        evolveMeshValues(frMesh, cst, 1.0-cst);
        //std::clog << reference() << " lattice avg = " << frMesh->sum()*frMesh->unit()/length() << std::endl;
    }
    //std::clog << frMesh->sum() << std::endl;
#endif
}


//------------------------------------------------------------------------------

/**
 The rest of the initialization is done in FiberProp::newFiber(),
 and other newFiber() functions where the initial length is known.
 */
Fiber::Fiber(FiberProp const* p)
: handListFront(nullptr), handListBack(nullptr), prop(p), disp(nullptr)
{
    if ( prop )
    {
        segmentation(prop->segmentation);
        
        if ( prop->lattice )
        {
#if FIBER_HAS_LATTICE
            if ( prop->lattice_unit < REAL_EPSILON )
                throw InvalidParameter("the Lattice unit (lattice[1]) must be > 0");
            //Cytosim::log << reference() <<  " new Lattice" << std::endl;
            frLattice.setUnit(prop->lattice_unit);
#else
            //throw InvalidParameter("Cytosim does not support fiber:lattice");
#endif
        }
        
        if ( prop->mesh )
        {
#if FIBER_HAS_MESH
            if ( prop->mesh_unit < REAL_EPSILON )
                throw InvalidParameter("the Mesh unit (mesh[1]) must be > 0");
            //Cytosim::log << reference() <<  " new value Lattice" << std::endl;
            frMesh.setUnit(prop->mesh_unit);
#else
            //throw InvalidParameter("Cytosim does not support fiber:mesh");
#endif
        }
    }
#if FIBER_HAS_GLUE
    frGlue = nullptr;
#endif
#if FIBER_HAS_FAMILY
    family_ = nullptr;
    sister_ = nullptr;
    brother_ = nullptr;
#endif
#if NEW_FIBER_CHEW
    frChewM = 0;
    frChewP = 0;
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
    if ( family_ )
    {
        Vector I = family_->posM(ab);
        Vector O = posM(ab);
        return I + 0.703 * ( O - I );
    }
    return posM(ab);
}

Vector Fiber::radialDir(real ab) const
{
    assert_true( family_ );
    ab -= abscissaM();
    return posM(ab) - family_->posM(ab);
}

#endif

Fiber::~Fiber()
{
    detachHands();
    
#if FIBER_HAS_MESH
    if ( prop->field_ptr )
        releaseMeshValues(frMesh, prop->field_ptr);
#endif

#if FIBER_HAS_GLUE
    delete(frGlue);
    frGlue = nullptr;
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
    for ( size_t ii = 0; ii < nbSegments(); ++ii )
    {
        //check the segment:
        FiberSegment s(this, ii);
        real d = INFINITY;
        real a = s.projectPoint0(w, d);
        if ( len < a )
        {
            // test exact point
            real e = distanceSqr(w, posP(ii+1));
            if ( e < dis ) {
                abs = abscissaPoint(ii+1);
                dis = e;
            }
        }
        else if ( 0 <= a  &&  d < dis )
        {
            //the projection is the best found so far
            abs = abscissaPoint(ii) + a;
            dis = d;
        }
    }
    
    return abs;
}

//------------------------------------------------------------------------------
#pragma mark - Modifying

void Fiber::flipHandsPolarity()
{
    /* update abscissa of Hands to keep them in place:
     new_abs - minus = plus - old_abs
     new_abs = plus + minus - old_abs
     */
    real mid = abscissaM() + abscissaP();
    Hand * ha = handListFront;
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
    
    Hand * ha = handListFront;
    while ( ha )
    {
        Hand * nx = ha->next();
        if ( ha->abscissa() < abs )
            ha->detach();
        ha = nx;
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
    
    Hand * ha = handListFront;
    while ( ha )
    {
        Hand * nx = ha->next();
        if ( ha->abscissa() > abs )
            ha->detach();
        ha = nx;
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
    
    // the signature on both pieces should be conserved:
    fib->signature(signature());
    fib->birthTime(birthTime());

    assert_true( fib->abscissaM() == abscissaM() );
    // remove MINUS_END portion on new piece:
    fib->truncateM(pti);
    assert_true(fib->abscissaM() == abs);
    
#if FIBER_HAS_MESH
    if ( frMesh.ready() )
    {
        assert_true( frMesh.unit() == fib->frMesh.unit() );
        // transfer Lattice values located above the cut
        fib->frMesh.takeP(frMesh, frMesh.index_round(abs));
    }
#endif

    // remove PLUS_END portion on self
    truncateP(pti);
    
    // transfer Hands above point P, at same abscissa
    Hand * ha = handListFront;
    while ( ha )
    {
        Hand * nx = ha->next();
        if ( ha->abscissa() > abs )
            ha->relocate(fib);
        else
            ha->update();
        ha = nx;
    }
    
    resetLattice();
    fib->resetLattice();
    
    return fib;
}


/**
The Fiber is cut at distance `abs` from its MINUS_END:
 - current Fiber is truncated to keep only the section [ MINUS_END , abs ],
 - A new Fiber is created representing the other section [ abs , PLUS_END ],
 - Hands are transfered to the new Fiber if appropriate,
 - lattice substances are also transfered,
 .
 A pointer to the new Fiber is returned (containing the PLUS_END), but this
 pointer may be zero, if `abs` was not within the valid range of abscissa.
 If a new Fiber was created, it should be added to the FiberSet.
 */
Fiber* Fiber::severP(real abs)
{
    if ( abs <= REAL_EPSILON || abs + REAL_EPSILON >= length() )
        return nullptr;
    
    //std::clog << "severP " << reference() << " at " << abscissaM()+abs << "\n";

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
    // remove MINUS_END portion on new piece
    fib->Chain::cutM(abs);

#if FIBER_HAS_MESH
    if ( frMesh.ready() )
    {
        assert_true( frMesh.unit() == fib->frMesh.unit() );
        // ensure valid range:
        fib->frMesh.setRange(abscissaM()+abs, abscissaP());
        
        // transfer Lattice values located above the cut:
        fib->frMesh.takeP(frMesh, frMesh.index_round(abscissaM()+abs));
    }
#endif
    
    assert_small(fib->abscissaM()-abs-abscissaM());
    
    // remove PLUS_END portion on self
    Chain::cutP(length()-abs);
    
    assert_small(fib->abscissaM()-abscissaP());

    // transfer all Hands above cut to new piece
    // their abscissa should not change in this transfer
    abs += abscissaM();
    Hand * ha = handListFront;
    while ( ha )
    {
        Hand * nx = ha->next();
        if ( ha->abscissa() >= abs )
            ha->relocate(fib);
        else
            ha->update();
        ha = nx;
    }

    resetLattice();
    fib->resetLattice();

    return fib;
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
    for ( SeverPos const& cut : pendingCuts )
    {
        if ( cut.abs - abscissaM() <= prop->min_length )
        {
            // we check the range again, since the fiber tip may have changed:
            if ( cut.abs > abscissaM() )
                cutM(cut.abs-abscissaM());
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
                cutP(abscissaP()-cut.abs);
        }
        else
        {
            Fiber * frag = severP(cut.abs-abscissaM());
            
            // special case where the PLUS_END section is simply deleted
            if ( cut.stateM == STATE_BLACK )
            {
                delete(frag);
                continue;
            }

            //add new fragment to simulation:
            objset()->add(frag);

            if ( frag )
            {
                // check that ends spatially match:
                assert_small((frag->posEndM() - posEndP()).norm());
                
                try {
                    // old PLUS_END converves its state:
                    frag->setDynamicStateP(dynamicStateP());
                    
                    // new ends are set as wished:
                    this->setDynamicStateP(cut.stateP);
                    frag->setDynamicStateM(cut.stateM);
                }
                catch ( Exception & e )
                {
                    e << "while cutting fiber " << reference();
                    throw;
                }
            
#ifdef LOGGING
                Cytosim::log << "severed " << reference() << " at abscissa " << s->abs;
                Cytosim::log << "   creating " << frag->reference();
                Cytosim::log << "   position " << frag->posEndM() << std::endl;
#endif
                //Cytosim::log << " severed at X = " << frag->posEndM().XX << std::endl;
            }
            else
            {
                Cytosim::log << " sever abscissa " << cut.abs << " is out of range";
                Cytosim::log << " [ " << abscissaM() << "   " << abscissaP() << " ]" << std::endl;
            }
        }
    }
    pendingCuts.clear();
}


/**
 returns index of first point for which ( cos(angle) < max_cosine ),
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
    
    for ( real s : cuts )
    {
        Fiber * fib = severNow(s);
        if ( fib )
        {
            // old PLUS_END converves its state:
            fib->setDynamicStateP(dynamicStateP());
            // dynamic of new ends are set as usual:
            setDynamicStateP(stateP);
            fib->setDynamicStateM(stateM);
            //assert_true(!fib->linked());
            objset()->add(fib);
        }
    }
}


/**
 The given `fib` is added past the PLUS_END of `*this`,
 Hands bound to `fib` are transfered to *this.
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
    
    //transfer dynamic state of PLUS_END:
    setDynamicStateP(fib->dynamicStateP());

#if FIBER_HAS_MESH
    if ( frMesh.ready() )
    {
        assert_true( frMesh.unit() == fib->frMesh.unit() );
        // ensure valid range:
        frMesh.setRange(abscissaM(), abscissaP());
        
        // transfer Lattice values from other fiber
        frMesh.takeP(fib->frMesh, frMesh.indexM());
    }
#endif

    // transfer all Hands
    Hand * ha = fib->handListFront;
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

/**
 From "Random Walks in Biology" by HC. Berg, Princeton University Press,
 drag coefficients for an ellipsoid are,

     drag_transverse = 2*drag_parallel = 4*PI*L*visc / log(length/radius)

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


/** 
 dragCoefficientVolume() calculates the mobility for the entire fiber,
 considering that the cylinder is straight and moving in an infinite fluid.
 fiber:hydrodynamic_radius[1] is a hydrodynamic cutoff that makes the
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
real Fiber::dragCoefficientVolume()
{
    real len = length();
    assert_true( len > 0 );
    
    // hydrodynamic cut-off on length:
    assert_true( prop->hydrodynamic_radius[1] > REAL_EPSILON );
    real lenc = std::min(len, prop->hydrodynamic_radius[1]);
    lenc = std::max(lenc, prop->hydrodynamic_radius[0]);
    
    //Stoke's for a sphere:
    assert_true( prop->hydrodynamic_radius[0] > 0 );
    real drag_sphere = 6 * prop->hydrodynamic_radius[0];
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    constexpr real pref = 4;
#else
    constexpr real pref = 3;
#endif

#if ( 0 )
    /*
     For an ellipsoid,
     drag_transverse = 2*drag_parallel = 4*PI*L*visc / log(length/radius)
     We should average the mobility coefficients:  speed = mu * f
           mu_X = mu_parallel   = 2 * mu
           mu_Y = mu_transverse = mu
           mu_Z = mu_transverse = mu
     Hence:
           mu_averaged = ( mu + mu + 2*mu ) / 3 = 4/3 * mu.
     drag_averaged = 3*PI*length*viscosity / log(length/radius)
     See for example "Random Walks in Biology" by HC. Berg, Princeton University Press.
     */
    
    // length below which the formula is not valid:
    real min_len = exp( 1 + log(prop->hydrodynamic_radius[0]) );

    real drag_cylinder = pref * len / log( lenc / prop->hydrodynamic_radius[0] );
#else
    /*
     Tirado and de la Torre. J. Chem. Phys 71(6) 1979
     give the averaged translational friction coefficient for a cylinder:
     3*PI*length*viscosity / ( log(length/diameter) + 0.32 )
     (Page 2584, Table 1, last column, last line for infinite aspect ratio)
     */
    
    // length below which the formula is not valid anymore ( ~ 3.94 * radius )
    // this corresponds to the minimun point of `drag_cylinder`
    real min_len = 2 * prop->hydrodynamic_radius[0] * exp(1.0-0.32);
    
    real drag_cylinder = pref * len / ( log(0.5*lenc/prop->hydrodynamic_radius[0]) + 0.32 );
#endif

    real drag = drag_sphere;
    
    if ( len > min_len )
    {
        // use largest drag coefficient
        drag = std::max(drag_cylinder, drag_sphere);
    }
    
    //Cytosim::log("Drag coefficient of Fiber in infinite fluid = %.1e\n", drag);
    //std::cerr << "Fiber with length " << len << " has bulk_drag = " << drag << '\n';

    return M_PI * prop->viscosity * drag;
}


/**
 dragCoefficientSurface() uses a formula calculated by F. Gittes in:\n
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

       acosh(x) = ln[ x + sqrt(x^2-1)) ] ~ ln[2x] if x >> 1

 Hunt et al. also credit this reference for the formula:\n
 <em>
 The slow motion of a cylinder next to a plane wall.
 Jeffrey, D.J. & Onishi, Y. (1981) Quant. J. Mech. Appl. Math. 34, 129-137.
 </em>
*/
real Fiber::dragCoefficientSurface()
{
    real len = length();
    
    if ( prop->cylinder_height <= 0 )
        throw InvalidParameter("fiber:surface_effect[1] (height above surface) must set and > 0!");
    
    // use the higher drag: perpendicular to the cylinder (factor 2)
    real drag = 2 * M_PI * prop->viscosity * len / acosh( 1 + prop->cylinder_height/prop->hydrodynamic_radius[0] );
    
    //Cytosim::log("Drag coefficient of Fiber near a planar surface = %.1e\n", drag);
    //std::cerr << "Fiber with length " << len << " has surface_drag = " << drag << '\n';

    return drag;
}


/**
 Calculate drag coefficient from two possible formulas

     if ( fiber:surface_effect )
        drag = dragCoefficientSurface();
     else
        drag = dragCoefficientVolume();

 */
void Fiber::setDragCoefficient()
{
    real drag;
    
    if ( prop->surface_effect )
    {
        drag = dragCoefficientSurface();
#if ( 0 )
        real d = dragCoefficientVolume();
        Cytosim::log << "Drag coefficient of Fiber near a planar surface amplified by " << drag/d << std::endl;
#endif
    }
    else
        drag = dragCoefficientVolume();

    //the forces are distributed equally on all points, hence we multiply by nPoints
    assert_true( nPoints > 0 );
    rfPointMobility = nPoints / drag;
    
#if ( 0 )
    Cytosim::log << "Fiber L = " << std::setw(7) << length();
    Cytosim::log << " drag = " << drag << " point_mobility " << rfPointMobility << std::endl;
#endif
}


void Fiber::prepareMecable()
{
    setDragCoefficient();
    storeDirections();
    makeProjection();
    //printProjection(std::clog);
    
    assert_true( rfPointMobility > REAL_EPSILON );

    // the scaling of the bending elasticity depends on the length of the segments
    rfRigidity = prop->rigidity / segmentationCube();
#if NEW_FIBER_LOOP
    rfRigidityLoop = prop->loop;
#endif
#if ( 0 )
    real energy = bendingEnergy();
    real euler = M_PI * M_PI * prop->rigidity / ( length() * length() );
    Cytosim::log << "Euler buckling = " << euler << "    ";
    Cytosim::log << "Bending energy = " << energy << std::endl;
#endif
}


//------------------------------------------------------------------------------

void Fiber::setInteractions(Meca& meca) const
{
#if OLD_SQUEEZE_FORCE
    if ( prop->squeeze == 1 )
    {
        // squeezing force in the YZ-plane:
        const real f = prop->squeeze_force;
        const real r = prop->squeeze_range;
        for ( size_t p = 0; p < nPoints; ++p )
        {
            Vector p = posP(pp);
            if ( p.normYZ() < r )
                meca.addLineClamp(Mecapoint(this, p), Vector(p.XX, 0, 0), Vector(1,0,0), f/r);
            else
                meca.addForce(Mecapoint(this, p), f);
        }
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
    {
        Space const* spc = prop->confine_space_ptr;
        switch ( prop->confine )
        {
            case CONFINE_INSIDE:
                for ( size_t i = 0; i < nPoints; ++i )
                {
                    Vector pos = posP(i);
                    if ( spc->outside(pos) )
                        spc->setInteraction(pos, Mecapoint(this, i), meca, prop->confine_stiffness);
                }
                break;
                
            case CONFINE_OUTSIDE:
                for ( size_t i = 0; i < nPoints; ++i )
                {
                    Vector pos = posP(i);
                    if ( spc->inside(pos) )
                        spc->setInteraction(pos, Mecapoint(this, i), meca, prop->confine_stiffness);
                }
                break;
                
            case CONFINE_ON:
                for ( size_t i = 0; i < nPoints; ++i )
                    spc->setInteraction(posP(i), Mecapoint(this, i), meca, prop->confine_stiffness);
                break;
                
            case CONFINE_MINUS_END:
                spc->setInteraction(posP(0), Mecapoint(this, 0), meca, prop->confine_stiffness);
                break;

            case CONFINE_PLUS_END:
            {
                const size_t L = lastPoint();
                spc->setInteraction(posP(L), Mecapoint(this, L), meca, prop->confine_stiffness);
            } break;
                
            case CONFINE_BOTH_ENDS:
            {
                const size_t L = lastPoint();
                spc->setInteraction(posP(0), Mecapoint(this, 0), meca, prop->confine_stiffness);
                spc->setInteraction(posP(L), Mecapoint(this, L), meca, prop->confine_stiffness);
            } break;
                
            case CONFINE_PLUS_OUT:
            {
                const size_t L = lastPoint();
                Vector pos = posP(L);
                if ( spc->inside(pos) )
                    spc->setInteraction(pos, Mecapoint(this, L), meca, prop->confine_stiffness);
            } break;
#if NEW_FIBER_CONFINE_RANGE
            case CONFINE_RANGE:
            {
                // we use here the MINUS_END as a reference... which maybe problematic
                size_t S = clampedIndexM(prop->confine_range[0]);
                size_t E = clampedIndexM(prop->confine_range[1]);
                for ( size_t i = S; i < E; ++i )
                    spc->setInteraction(posP(i), Mecapoint(this, i), meca, prop->confine_stiffness);
            } break;
#endif
            default:
                throw InvalidParameter("Invalid fiber::confine");
        }
    }
    
#if NEW_FIBER_CONFINE2
    /// add another confinement force
    if ( prop->confine2 != CONFINE_OFF )
    {
        Space const* spc = prop->confine2_space_ptr;
        
        switch ( prop->confine2 )
        {
            case CONFINE_INSIDE:
                for ( size_t i = 0; i < nPoints; ++i )
                {
                    Vector pos = posP(i);
                    if ( spc->outside(pos) )
                        spc->setInteraction(pos, Mecapoint(this, i), meca, prop->confine2_stiffness);
                }
                break;
                
            case CONFINE_OUTSIDE:
                for ( size_t i = 0; i < nPoints; ++i )
                {
                    Vector pos = posP(i);
                    if ( spc->inside(pos) )
                        spc->setInteraction(pos, Mecapoint(this, i), meca, prop->confine2_stiffness);
                }
                break;
                
            case CONFINE_ON:
                for ( size_t i = 0; i < nPoints; ++i )
                    spc->setInteraction(posP(i), Mecapoint(this, i), meca, prop->confine2_stiffness);
                break;
                
            default:
                throw InvalidParameter("Invalid fiber::confine2");
        }
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - list of bound Hands

/**
Link the hand at the front of the list.

Keeping track of the bound Hands is needed to run Cytosim if the filaments are
dynamic, so that a Fiber that has changed can update all the Hands bound to it.
The list is also used in reports, or to quickly count Hands bound to the fiber.
*/
void Fiber::addHand(Hand * n) const
{
    n->prev(nullptr);
    n->next(handListFront);
    if ( handListFront )
        handListFront->prev(n);
    else
        handListBack = n;
    handListFront = n;
}


void Fiber::removeHand(Hand * n) const
{
    Hand * x = n->next();
    if ( n->prev() )
        n->prev()->next(x);
    else {
        assert_true( handListFront == n );
        handListFront = x;
    }
    
    if ( x )
        x->prev(n->prev());
    else {
        assert_true( handListBack == n );
        handListBack = n->prev();
    }
}


void Fiber::updateHands() const
{
    for ( Hand * ha = handListFront; ha; ha = ha->next() )
        ha->update();
}


void Fiber::detachHands() const
{
    // we must iterate one step ahead, because detach() will unlink
    Hand * ha = handListFront;
    while ( ha )
    {
        Hand * nx = ha->next();
        ha->detach();
        ha = nx;
    }
}

/**
Sort in ascending order
*/
int compareAbscissa(const void* a, const void* b)
{
    real ax = static_cast<Hand const*>(a)->abscissa();
    real bx = static_cast<Hand const*>(b)->abscissa();
    return ( ax > bx ) - ( bx > ax );
}

/**
 This sorts the Hands in order of increasing abscissa
 Sorting is done by copying to temporary array space, using std::qsort
 */
void Fiber::sortHands() const
{
    size_t cnt = nbHands();
    if ( cnt > 1 )
    {
        Hand ** tmp = new Hand*[cnt];
        
        size_t i = 0;
        Hand * n = handListFront;
        
        while ( n )
        {
            tmp[i++] = n;
            n = n->next();
        }
        
        qsort(tmp, cnt, sizeof(Hand*), compareAbscissa);
        
        n = tmp[0];
        handListFront = n;
        n->prev(nullptr);
        for ( i = 1; i < cnt; ++i )
        {
            n->next(tmp[i]);
            tmp[i]->prev(n);
            n = tmp[i];
        }
        n->next(nullptr);
        handListBack = n;
        
        delete[] tmp;
    }
}


size_t Fiber::nbHands() const
{
    size_t res = 0;
    
    for ( Hand const* ha = handListFront; ha; ha = ha->next() )
        ++res;
    
    return res;
}


int Fiber::nbHands(int (*count)(Hand const*)) const
{
    int res = 0;
    
    for ( Hand const* ha = handListFront; ha; ha = ha->next() )
        res += count(ha);
    
    //printf("nbHands(%p) = %u\n", count, res);
    return res;
}


size_t Fiber::nbHandsInRange(real a, real b, const FiberEnd ref) const
{
    size_t res = 0;
    if ( b > a )
    {
        // Convert to absolute abscissa:
        a = abscissaFrom(a, ref);
        b = abscissaFrom(b, ref);
        if ( b < a )
            std::swap(a, b);
        
        for ( Hand const* ha = handListFront; ha; ha = ha->next() )
            res += ( a <= ha->abscissa()  &&  ha->abscissa() <= b );
    }
    //printf("nbHandsInRange(%8.2f, %8.2f) = %u\n", a, b, res);
    return res;
}


size_t Fiber::nbHandsNearEnd(const real len, const FiberEnd ref) const
{
    size_t res = 0;
    real i = -INFINITY, s = INFINITY;
    
    if ( len > 0 )
    {
        if ( ref == PLUS_END )
            i = abscissaP() - len;
        else if ( ref == MINUS_END )
            s = abscissaM() + len;
        else
            throw("invalid argument value to nbHadsNearEnd()");
        
        for ( Hand const* ha = handListFront; ha; ha = ha->next() )
            res += ( i <= ha->abscissa() && ha->abscissa() <= s );
    }
    //printf("nbHandsNearEnd(%8.2f) = %u\n", len, res);
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Dynamic ends

state_t Fiber::dynamicState(FiberEnd end) const
{
    if ( end == PLUS_END )
        return dynamicStateP();
    if ( end == MINUS_END )
        return dynamicStateM();
    ABORT_NOW("invalid argument value");
    return 0;
}


void Fiber::setDynamicState(const FiberEnd end, const state_t s)
{
    if ( end == PLUS_END )
        setDynamicStateP(s);
    else if ( end == MINUS_END )
        setDynamicStateM(s);
}


real Fiber::freshAssembly(const FiberEnd end) const
{
    if ( end == PLUS_END )
        return freshAssemblyP();
    if ( end == MINUS_END )
        return freshAssemblyM();
    ABORT_NOW("invalid argument value");
    return 0;
}


/**
 Assuming that the length has changed, or that the abscissa of the ends have changed,
 this updates the segmentation of the fiber if needed, the position of the Hands,
 and the boundaries of the Lattice if present.
 */
void Fiber::updateFiber()
{
#if ( 0 )
    Cytosim::log << reference() << " update [ "  << std::setw(9) << std::left << abscissaM();
    Cytosim::log << " "  << std::setw(9) << std::left << abscissaP() << " ]" << std::endl;
#endif
    
    /*
     Update all bound Hands.
     Some Hands may be updated more than once in a time-step,
     but that is only a small performance penalty.

     We must iterate one step ahead, because `checkFiberRange()` may lead to detachment.
     The loop could be unrolled, or parallelized
     */
    Hand * ha = handListFront;
    while ( ha )
    {
        Hand * nx = ha->next();
        assert_true(ha->fiber()==this);
        ha->update();
        ha->checkFiberRange();
        ha = nx;
    }
    
#if FIBER_HAS_LATTICE
    // this will allocate the Lattice's site to cover the range of Abscissa:
    if ( frLattice.ready() )
    {
        frLattice.setRange(abscissaM(), abscissaP());
    }
#endif
#if FIBER_HAS_MESH
    // this will allocate the Lattice's site to cover the range of Abscissa:
    if ( frMesh.ready() )
    {
        frMesh.setRange(abscissaM(), abscissaP());
        
        if ( prop->field_ptr )
        {
            real sumM;
            // release Lattice substance located outside the valid abscissa range
            frMesh.collectM(sumM);
            prop->field_ptr->cell(posEndM()) += sumM;
            //Cytosim::log << " Fiber::MINUS_END releases " << sumM << std::endl;
            
            real sumP;
            frMesh.collectP(sumP);
            prop->field_ptr->cell(posEndP()) += sumP;
            //Cytosim::log << " Fiber::PLUS_END releases " << sumP << std::endl;
        }
    }
#endif
    needUpdate = false;
}

//------------------------------------------------------------------------------
#pragma mark - Lattice


void Fiber::printLattice(std::ostream& os) const
{
#if FIBER_HAS_LATTICE
    FiberLattice const& lat = frLattice;
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


void Fiber::infoLattice(real& len, size_t& cnt, real& sm, real& mn, real& mx) const
{
#if FIBER_HAS_LATTICE
    FiberLattice const& lat = frLattice;
    if ( lat.ready() )
    {
        len += length();
        const auto sup = lat.indexP();
        for ( auto i = lat.indexM(); i <= sup; ++i )
        {
            ++cnt;
            sm += lat.data(i);
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
    return &frMesh;
#elif FIBER_HAS_LATTICE
    return &frLattice;
#else
    return nullptr;
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Mesh

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
    assert_true(lat.ready());
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();

    //std::clog << "evolve " << inf << " " << sup << "\n";
    for ( auto h = inf; h <= sup; ++h )
        site[h] = cst + fac * site[h];
}


void Fiber::bindMesh(Lattice<real>& lat, Field * fld, real binding_rate) const
{
    assert_true(lat.ready());
    // we want roughly one point per cell:
    const real spread = fld->cellWidth();
    
    // each point represents a Fiber chunk of length 'spread':
    const real rate = binding_rate * spread / fld->cellVolume();
    
    // fraction of the cell content that will bind in one time_step:
    const real frac = -std::expm1( -rate * simul().prop->time_step );
    
    real abs = spread * RNG.exponential();
    const real len = length();
    
    // stochastic sampling with distance 'spread' along the Fiber:
    while ( abs < len )
    {
        Vector pos = posM(abs);
        real& cell = fld->cell(pos);
        assert_true( cell >= 0 );
        
        // amount to be transfered:
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
    assert_true( lat.ready() );
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    auto * site = lat.data();

    const real fac = speed * simul().prop->time_step / lat.unit();
    
    if ( fabs(fac) > 1 )
        throw InvalidParameter("mesh_flux_speed * time_step / lattice_unit is too high");

    if ( fac < 0 )
    {
        real s = site[inf];
        
        for ( auto h = inf; h < sup; ++h )
            site[h] -= fac * ( site[h+1] - site[h] );
        
        prop->field_ptr->cell(posEndM()) -= fac * s;
        site[sup] += fac * site[sup];
    }
    else
    {
        real s = site[sup];
        
        for ( auto h = sup; h > inf; --h )
            site[h] -= fac * ( site[h] - site[h-1] );
        
        prop->field_ptr->cell(posEndP()) += fac * s;
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

    const real fac = 1.0 / simul().prop->time_step;
    
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


void Fiber::infoMesh(real& len, size_t& cnt, real& sm, real& mn, real& mx, bool density) const
{
#if FIBER_HAS_MESH
    Lattice<real> const& lat = frMesh;
    if ( lat.ready() )
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
 fiber:glue=1 creates a Single if the tip of the Fiber goes outside the Space.
 The Single's hand is managed to always be attached at the tip of the Fiber.
 The Single detaches if the Fiber tip is pulled inside.
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
 fiber:glue=2
 The Single's hand is managed to always be attached at the tip of the Fiber.
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
 fiber:glue=3 creates a Single at the position where the Fiber crosses the Space's edge.
 This makes an anchor point exactly at the cortex.
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
                    FiberSite fs(this, abscissaPoint(i-1+d1/(d2+d1)));
                    glue->setPosition(fs.pos());
                    glue->attach(fs);
                    break;
                }
            }
        }
    }
}


/**
 Search for a glue in the list of bound HandSingle
 this is useful when a simulation is restarted from file
 */
void Fiber::makeGlue(Single*& glue)
{
    SingleSet& set = simul().singles;

    for ( Single * gh = set.firstA(); gh; gh=gh->next() )
    {
        if ( gh->hand()->fiber() == this  &&  gh->mark() == identity() )
        {
            glue = gh;
            //Cytosim::log << "found Fiber:glue for " << reference() << std::endl;
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
 setGlue() creates Single when MT interact with the edge of the Space
*/
void Fiber::setGlue(Single*& glue, const FiberEnd end, Space const* space)
{
    assert_true(space);
    
    if ( !glue )
        makeGlue(glue);
    
    switch( prop->glue )
    {
        case 1:  setGlue1(glue, end, space);  break;
        case 2:  setGlue2(glue, end, space);  break;
        case 3:  setGlue3(glue, space);       break;
        default: throw InvalidParameter("invalid value of fiber:glue");
    }
    
#if ( 1 )
    // we keep the Single linked in the simulation only if it is attached:
    if ( glue->attached() )
    {
        if ( !glue->linked() )
            simul().singles.link(glue);
    }
    else
    {
        if ( glue->linked() )
            simul().singles.unlink(glue);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - I/O

void Fiber::write(Outputter& out) const
{
    Chain::write(out);

#if FIBER_HAS_LATTICE && 0
    /*
     We can save the occupancy Lattice here, but this is not necessary
     as it can be recalculated on the fly, so we save space by skipping
     */
    if ( frLattice.ready() )
    {
        writeHeader(out, TAG_LATTICE);
        // frLattice.write(out);
        // only write information corresponding to actual Fiber abscissa range:
        frLattice.write(out, frLattice.indexM(), frLattice.indexP()+1);
    }
#endif
#if FIBER_HAS_MESH
    if ( frMesh.ready() )
    {
        writeHeader(out, TAG_FIBMESH);
        // frMesh.write(out);
        // only write information corresponding to actual Fiber abscissa range:
        frMesh.write(out, frMesh.indexM(), frMesh.indexP()+1);
    }
#endif
}


void Fiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    //std::clog << this << " Fiber::read(" << tag << ")\n";
#ifdef BACKWARD_COMPATIBILITY
        
    if ( in.formatID() == 33 )
        mark(in.readUInt32());
    
    if ( tag == 'm' )
    {
        if ( in.formatID()==31 )
        {
            size_t p = in.readUInt16();
            prop = sim.findProperty<FiberProp>("fiber", p);
        }
        //tag = TAG;
    }
    
    if ( in.formatID() < 31 )
    {
        setDynamicStateM(in.readUInt8());
        setDynamicStateP(in.readUInt8());
    }
    
#endif
    
    if ( tag == TAG )
    {
        Chain::read(in, sim, tag);
        
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() > 47 && in.formatID() < 50 ) // 4.7.2018 added birthTime
            fnBirthTime = in.readFloat();
#endif

        if ( length() + REAL_EPSILON < prop->min_length )
        {
            Cytosim::log << "Warning: fiber length < fiber:min_length";
            Cytosim::log << " ( " << length() << " < " << prop->min_length << " )" << std::endl;
        }
#if FIBER_HAS_GLUE
        frGlue = nullptr;
#endif
    }
    else if ( tag == TAG_LATTICE )
    {
        try {
#if FIBER_HAS_LATTICE
            frLattice.read(in);
#else
            FiberLattice dummy;
            dummy.read(in);
            // store unit, to get digits at the right abscissa
            const_cast<FiberProp*>(prop)->lattice_unit = dummy.unit();
#endif
        }
        catch( Exception & e ) {
            e << "reading Lattice for " << reference();
            throw;
        }
    }
    else if ( tag == TAG_FIBMESH )
    {
        try {
#if FIBER_HAS_MESH
            frMesh.read(in);
#else
            Lattice<real> dummy;
            dummy.read(in);
#endif
        }
        catch( Exception & e ) {
            e << "reading value Lattice for " << reference();
            throw;
        }
    }

#ifdef BACKWARD_COMPATIBILITY
    else if ( tag == TAG_DYNAMIC )
    {
        static bool virgin = true;
        // that is for Fiber class we do not know...
        if ( virgin )
        {
            virgin = false;
            std::cerr << "INCIDENT: skipping dynamic states for unknown fiber class\n";
        }
        in.readUInt32();
        in.readUInt32();
    }
    else if ( tag == 'm' )
    {
        Chain::read(in, sim, tag);
        updateFiber();
    }
#endif
    else
    {
        Cytosim::log << "unknown Fiber TAG `" << (char)tag << "'" << std::endl;
    }
}


