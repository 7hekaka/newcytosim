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
#include "random.h"  

#include <unordered_map>

struct ClusterState {
    Vector C;     // center on inner surface
    real   phi=0; // in-plane rotation angle
    real   t_last = -1;
    bool   init=false;
};
static std::unordered_map<const SingleProp*, ClusterState> g_cluster;

static inline Vector radial_xy(const Vector& p) {
    real r = std::sqrt(p.XX*p.XX + p.YY*p.YY);
    if (r > REAL_EPSILON) return Vector(p.XX/r, p.YY/r, 0);
    return Vector(1,0,0);
}
static inline Vector tan_azimuth(const Vector& n) { return Vector(-n.YY, n.XX, 0); } // φ_hat
static inline Vector project_inner_annulus(const Space* S, Vector P); // forward (see below)

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
        sPos = h->unbindingPosition();
        // move to correct SingleSet sublist:
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
    if ( prop->confine == CONFINE_ON )
        sPos = prop->confine_space->placeOnEdge(1.0);
    else
        sPos = prop->confine_space->place();
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
        sPos = prop->confine_space->bounce(pos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        sPos = prop->confine_space->project(pos);
    }
    else
    {
        sPos = pos;
    }
    
    sHand->stepUnattached(simul(), sPos);


auto& sp = *prop;
if ( sp.confine == CONFINE_ON && sp.confine_space ) {
    auto& cs = g_cluster[prop];   // 'prop' is already a const SingleProp*

    if (!cs.init) {
        cs.C     = project_inner_annulus(sp.confine_space, sPos);
        cs.phi   = 0;
        cs.t_last = simul().time();
        cs.init  = true;
    }

    real dt  = simul().prop.time_step;
    if ( simul().time() > cs.t_last + 0.5*dt ) {
        Vector n  = radial_xy(cs.C);
        Vector et = tan_azimuth(n);
        Vector ez = Vector(0,0,1);
        real sigT = std::sqrt(2.0 * sp.cluster_D    * dt);
        real sigR = std::sqrt(2.0 * sp.cluster_Drot * dt);

        cs.C  += sigT * ( RNG.gauss() * et + RNG.gauss() * ez );
        cs.C   = project_inner_annulus(sp.confine_space, cs.C);
        if ( sp.cluster_Drot > 0 ) cs.phi += sigR * RNG.gauss();
        cs.t_last = simul().time();
    }

    if ( !c_off_init_ ) {
        Vector n  = radial_xy(cs.C);
        Vector d  = sPos - cs.C;      // current displacement
        d -= dot(d,n)*n;              // project to tangent plane
        c_off_T_  = d;
        c_off_init_ = true;
    }

    Vector n  = radial_xy(cs.C);
    Vector et = tan_azimuth(n), ez = Vector(0,0,1);
    real c = std::cos(cs.phi), s = std::sin(cs.phi);
    real a = dot(c_off_T_, et);
    real b = dot(c_off_T_, ez);
    Vector off = ( c*a - s*b ) * et + ( s*a + c*b ) * ez;

    Vector pos = cs.C + off;                  // desired position
    sPos = sp.confine_space->project(pos);    // enforce surface now
}
     
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
 */
void Single::setInteractions(Meca& meca) const
{
    assert_true( sHand->attached() );
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

// single.cc bottom (or space_annulus.cc if you prefer)

static inline Vector project_inner_annulus(const Space* S, Vector P)
{
    // For now, just use the generic surface projection.
    // If you later add explicit "inner-only" in the space class, call it here.
    return S->project(P);
}
//------------------------------------------------------------------------------