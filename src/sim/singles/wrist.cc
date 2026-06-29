// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "wrist.h"
#include "meca.h"
#include "modulo.h"
#include "single_set.h"
#include "simul.h"     // needed to access simul().prop.time_step
#include "random.h"    // RNG.gauss() etc.
#include "solid.h"



// Convert a world tangent step to triad-local coordinates and re-set Interpolation4
static inline void rebase_on_triad(Interpolation4& base, const Vector& new_world_pos)
{
    // 1) find the Solid and the triad reference index:
    Mecable const* m = base.mecable();
    const Solid* S   = Solid::toSolid(m);     // safe cast helper exists in this codebase
    if ( !S ) return;

    // We need the triad center and axes. For Solids, a triad is stored at indices [ref .. ref+3]
    // We can retrieve 'ref' via vertex0():
    const index_t ref = base.vertex0().point();   // Interpolation4 exposes vertex0()
    const Vector  C   = S->posP(ref);             // center
    const Vector  X   = S->posP(ref+1) - C;       // triad arm X
    const Vector  Y   = S->posP(ref+2) - C;       // triad arm Y
    const Vector  Z   = S->posP(ref+3) - C;       // triad arm Z

    // 2) build 3×3 from columns [X Y Z] and solve M * v = (new_world_pos - C)
    //    Since triad is constructed orthonormally (Solid::addTriad), a quick projection works:
    const Vector d  = new_world_pos - C;
    const real xL2 = X.normSqr();
    const real yL2 = Y.normSqr();
    const real zL2 = Z.normSqr();
    Vector v(
        xL2 > REAL_EPSILON ? dot(d, X)/xL2 : 0,
        yL2 > REAL_EPSILON ? dot(d, Y)/yL2 : 0,
        zL2 > REAL_EPSILON ? dot(d, Z)/zL2 : 0
    );

    // 3) re-encode the anchor position in triad-local coordinates:
    base.set(m, ref, v);
}

Wrist::Wrist(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Single(sp)
{
    // 'mec' can be Null when reading from file
    rebase(mec, pti);
#if ( 0 )
    if ( p->diffusion > 0 )
        throw InvalidParameter(name()+":diffusion cannot be > 0 if activity=anchored");
#endif
}


Wrist::~Wrist()
{
}


Vector Wrist::stretch() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - sHand->pos();
    
    if ( modulo )
        modulo->fold(d);
    
    return d;
}


Vector Wrist::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - sHand->pos();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


void Wrist::stepF()
{
    assert_false( sHand->attached() );

    // 1) standard detached-hand kinetics at current foot:
    const Vector foot = posFoot();
    sHand->stepUnattached(simul(), foot);

    // 2) 2D Brownian kick along the surface:
    const real dt = simul().prop.time_step;                         // time step
    const real D  = std::max<real>(prop->diffusion, 0);             // <-- temporary knob
    const real sig = std::sqrt(2.0 * D * dt);

    // Get center and normal at foot, build tangent basis:
    const Mecable* m = base_.mecable();
    const Solid*   S = Solid::toSolid(m);
    if ( !S ) return;

    // center and radius from current triad
    const index_t ref = base_.vertex0().point();
    const Vector  C   = S->posP(ref);
    Vector n = base_.normal();
    if ( n.normSqr() < REAL_EPSILON ) n = (foot - C).normalized();

    Vector t1 = Vector::randU();   t1 -= dot(t1, n) * n;
    if ( t1.normSqr() < REAL_EPSILON ) t1 = Vector(1,0,0) - dot(Vector(1,0,0),n)*n;
    t1.normalize();
    const Vector t2 = cross(n, t1);

    // 2D Gaussian in tangent plane:
    const real xi1 = RNG.gauss();
    const real xi2 = RNG.gauss();
    Vector trial = foot + sig*( xi1 * t1 + xi2 * t2 );

    // Project back to the same spherical shell (keep radius constant):
    const real R = (foot - C).norm();
    if ( R > REAL_EPSILON ) {
        Vector u = trial - C;
        real un = u.norm();
        if ( un > REAL_EPSILON ) trial = C + (R/un) * u;
    }

    // 3) Re-encode in Interpolation4 (this actually moves the anchor on the solid)
    rebase_on_triad(base_, trial);
}


void Wrist::stepA()
{
    assert_true( sHand->attached() );
    Vector f = Wrist::force();
    if ( sHand->checkKramersDetachment(f.norm()) )
        sHand->detach();
    else
        sHand->stepLoaded(f);
}


void Wrist::setInteractions(Meca& meca) const
{
    Interpolation i = sHand->interpolation();
    base_.addLink(meca, i, prop->stiffness);
}


void Wrist::write(Outputter& out) const
{
    writeMarker(out, WRIST_TAG);
    sHand->writeHand(out);
    base_.write(out);
}


void Wrist::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    sHand->readHand(in, sim);
    
#if BACKWARD_COMPATIBILITY < 47
    if ( in.formatID() < 47 )
    {
        Mecapoint base;
        base.read(in, sim);
        base_.set(base.mecable(), base.point());
    }
    else
#endif
        base_.read(in, sim);
}

