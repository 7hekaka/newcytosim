// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#include "dim.h"
#include "fiber.h"
#include "tubule.h"
#include "simul.h"
#include "meca.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"


void Tubule::step()
{
}


Tubule::Tubule(TubuleProp * p) : prop(p)
{
    reset();
}


Tubule::~Tubule()
{
    prop = nullptr;
}


void Tubule::reset()
{
    bone_ = nullptr;
    for ( size_t i = 0; i < FILM; ++i )
    {
        fil_[i] = nullptr;
        // set as left-handed 3-start helix: 3 * 4 nm offset in one turn
        offset_[i] = i * ( -0.012 / NFIL );
    }
}


ObjectList Tubule::build(Glossary& opt, Simul& sim)
{
    FiberProp const* fp = sim.findProperty<FiberProp>("fiber", prop->fiber_type);
    real len = 1.0, dev = 0;
    ObjectList res;

    // get the 'bone'
    if ( prop->bone_type.size() > 0 )
    {
        res = sim.fibers.newObjects(prop->bone_type, opt);
        bone_ = static_cast<Fiber*>(res[0]);
        Buddy::connect(bone_);
        len = bone_->length();
        if ( bone_->prop->segmentation != fp->segmentation )
            throw InvalidParameter("Tubule's bone and filament should have equal segmentation");
    }
    else
    {
        opt.set(len, "length");
        if ( opt.set(dev, "length", 1) )
            len += dev * RNG.sreal();
        len = std::max(len, fp->min_length);
        len = std::min(len, fp->max_length);
    }

    // all constitutive fibers should have the same length!
    for ( size_t i = 0; i < NFIL; ++i )
    {
        Fiber * fib = fp->newFiber();
        fib->setProjection(false);
        fib->setOrigin(offset_[i]);
        fib->setStraight(Vector(-0.5*len,0,0), Vector(1,0,0), len);
        fib->updateFiber();
        if ( bone_ )
            fib->setPoints(bone_->addrPoints(), bone_->nbPoints());
        Buddy::connect(fib);
        res.push_back(fib);
        fil_[i] = fib;
    }
    
#if ( DIM >= 3 )
    Vector dir(0,0,0);
    if ( bone_ )
        dir = bone_->avgDirection();
    else
    {
        // find average direction
        for ( size_t i = 0; i < NFIL; ++i )
            dir += fil_[i]->diffPoints(0);
        dir.normalize();
    }
    
    // set orthonormal coordinate system
    Vector E(0,prop->radius,0), F(0,0,prop->radius);
    dir.orthonormal(E, F, prop->radius);
    
    // translate protofilaments to form a tube:
    real a = 0; //M_PI * RNG.sreal();
    real da = 2 * M_PI / NFIL;
    for ( size_t i = 0; i < NFIL; ++i )
    {
        fil_[i]->translate(std::cos(a)*E+std::sin(a)*F);
        a += da;
    }
#endif

    setFamily(bone_?bone_:fil_[0]);
    return res;
}


Vector Tubule::posCenterlineM(real dis)
{
    Vector vec(0, 0, 0);
    for ( size_t i = 0; i < NFIL; ++i )
        vec += fil_[i]->posM(dis);
    return vec / NFIL;
}


//------------------------------------------------------------------------------
#pragma mark - Family


void Tubule::setFamily(Fiber const* fam)
{
    // wrap array values for convenience
    for ( size_t i = NFIL; i < FILM; ++i )
        fil_[i] = fil_[i-NFIL];

#if FIBER_HAS_FAMILY
    for ( size_t i = 0; i < NFIL; ++i )
    {
        fil_[i]->family_ = fam;
        fil_[i]->sister_ = fil_[(i+NFIL-1)%NFIL];
        fil_[i]->brother_ = fil_[(i+1)%NFIL];
    }
    if ( bone_ )
        bone_->family_ = fam;
#else
    std::clog << "WARNING: to use Tubule, please compile with FIBER_HAS_FAMILY\n";
#endif
}


void Tubule::salute(Buddy const* guy)
{
    //std::clog << reference() << " salute(" << guy << ")\n";
    if ( guy == bone_ )
    {
        assert_true(bone_);
        if (( bone_->freshAssemblyP() != 0 ) | ( bone_->freshAssemblyM() != 0 ))
        {
            // copy the growth exhibited from `bone_`
            real aP = bone_->abscissaP();
            real aM = bone_->abscissaM();
            for ( size_t i = 0; i < NFIL; ++i )
            {
                fil_[i]->growP(aP-offset_[i]-fil_[i]->abscissaP());
                fil_[i]->growM(fil_[i]->abscissaM()-aM+offset_[i]);
                fil_[i]->adjustSegmentation();
                fil_[i]->updateFiber();
            }
            //report(std::clog);
        }
    }
}


void Tubule::goodbye(Buddy const* b)
{
    std::cerr << "ERROR: Tubule's filaments cannot be deleted\n";
    if ( b )
    {
        for ( size_t i = 0; i < NFIL+2; ++i )
            if ( fil_[i] == b )
            {
                fil_[i] = nullptr;
                return;
            }
        std::cerr << "Error: unknown cytosim buddy" << b << '\n';
    }
}

//------------------------------------------------------------------------------
#pragma mark - Interactions

/**
This uses only addSideLink() with appropriate directions
Cambridge, 12.2019
*/
void Tubule::setInteractionsB(Meca& meca) const
{
#if ( DIM >= 3 )
    const real stiff = prop->stiffness[0];
    const real ang = M_PI / NFIL;
    const real len = 2 * prop->radius * std::sin(ang);  // distance between protofilaments
    const real c = std::cos(ang), s = std::sin(ang);
    
    assert_true(fil_[0]);
    const size_t end = fil_[0]->nbSegments();
    
    Rotation mat(0,1);

    for ( size_t i = 0; i < end; ++i )
    {
        // get centerline
        Vector cen(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            cen += fil_[n]->posPoint(i);
        cen /= NFIL;
        
        // get average direction of the Tubule at this location:
        Vector dir(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            dir += fil_[n]->diffPoints(i);
        dir.normalize();
        mat = Rotation::rotationAroundAxis(dir, c, s);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            Vector arm = mat.vecmul(( cen - fil_[n]->posPoint(i) ).normalized(len));
            meca.addSideLinkMT(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), arm, stiff);
        }
    }

    // get centerline
    Vector cen(0,0,0);
    for ( size_t n = 0; n < NFIL; ++n )
        cen += fil_[n]->posPoint(end);
    cen /= NFIL;
    
    for ( size_t n = 0; n < NFIL; ++n )
    {
        Vector arm = mat.vecmul(( cen - fil_[n]->posPoint(end) ).normalized(len));
        meca.addSideLinkMT(Interpolation(fil_[n],end-1,end,1), Mecapoint(fil_[n+1],end), arm, stiff);
    }
#endif
}


/**
 This uses a centerline `bone`
 Cambridge, 18.01.2020
 */
void Tubule::setInteractions(Meca& meca) const
{
    if ( !bone_ )
        throw InvalidParameter("tubule:bone must be defined");

#if ( DIM >= 3 )
    const real stiffL = prop->stiffness[0];
    const real stiffR = prop->stiffness[1];
    const real ang = M_PI / NFIL;
    const real len = 2 * prop->radius * std::sin(ang);  // distance between protofilaments
    const real c = std::cos(ang), s = std::sin(ang);
    
    const size_t e = bone_->nbSegments();
    
    // check all filament length:
    for ( size_t n = 0; n < NFIL; ++n )
    {
        if ( fil_[n]->nbSegments() != e )
        {
            Cytosim::warn << "unequal Tubule filaments\n";
            return;
        }
    }
    
    Rotation mat(0,1);
    real alpha = prop->radius / len;
    Vector cen, dir;
    
    for ( size_t i = 0; i < e; ++i )
    {
        // get centerline
        cen = bone_->posPoint(i);
        dir = bone_->dirSegment(i);
        mat = Rotation::rotationAroundAxis(dir, c, s);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            Vector arm = ( cen - fil_[n]->posPoint(i) ).normalized(len);
            // orthoradial beams:
            meca.addSideLinkMT(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), mat.vecmul(arm), stiffL);
            // radial spoke:
            meca.addSideLinkMT(Interpolation(fil_[n],i,i+1,0), Mecapoint(bone_,i), alpha*cross(dir,arm), stiffR);
        }
    }

    // link last point:
    cen = bone_->posPoint(e);
    for ( size_t n = 0; n < NFIL; ++n )
    {
        Vector arm = ( cen - fil_[n]->posPoint(e) ).normalized(len);
        meca.addSideLinkMT(Interpolation(fil_[n],e-1,e,1), Mecapoint(fil_[n+1],e), mat.vecmul(arm), stiffL);
        meca.addSideLinkMT(Interpolation(fil_[n],e-1,e,1), Mecapoint(bone_,e), alpha*cross(dir,arm), stiffR);
    }
#endif
}


/**
This uses addSideLink() and addTorque() with appropriate directions
Cambridge, 12.2019
*/
void Tubule::setInteractionsC(Meca& meca) const
{
#if ( DIM >= 3 )
    const real ang = M_PI / NFIL;
    const real len = 2 * prop->radius * std::sin(ang);  // distance between protofilaments
    const real stiffL = prop->stiffness[0];
    const real stiffA = prop->stiffness[1];
    real co = std::cos(2*ang), si = std::sin(2*ang);
    
    assert_true(fil_[0]);
    const size_t end = fil_[0]->nbSegments();
    
    MatrixBlock mat;
    for ( size_t i = 0; i <= end; ++i )
    {
        // get centerline
        Vector cen(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            cen += fil_[n]->posPoint(i);
        cen /= NFIL;
        
        if ( i < end )
        {
            // get average direction of the Tubule at this location:
            Vector dir(0,0,0);
            for ( size_t n = 0; n < NFIL; ++n )
                dir += fil_[n]->diffPoints(i);
            dir.normalize();
            
            // create rotation matrix for torque:
            mat = Meca::torqueMatrix(stiffA, dir, co, si);
            
            for ( size_t n = 0; n < NFIL; ++n )
            {
                Vector arm = (2*cen - fil_[n]->posPoint(i)- fil_[n+1]->posPoint(i)).normalized(len);
                meca.addSideLinkMT(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), arm, stiffL);
            }
        }
        else
        {
            for ( size_t n = 0; n < NFIL; ++n )
            {
                Vector arm = (2*cen - fil_[n]->posPoint(i) - fil_[n+1]->posPoint(i)).normalized(len);
                meca.addSideLinkMT(Interpolation(fil_[n],i-1,i,1), Mecapoint(fil_[n+1],i), arm, stiffL);
            }
        }
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            meca.addTorque(Mecapoint(fil_[n],i), Mecapoint(fil_[n+1],i),
                           Mecapoint(fil_[n+2],i), mat, stiffA);
        }
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - I/O


void Tubule::write(Outputter& out) const
{
    out.writeUInt16(NFIL+1);
    out.writeSoftNewline();
    Object::writeReference(out, bone_);
    for ( size_t i = 0; i < NFIL; ++i )
    {
        out.writeSoftSpace();
        Object::writeReference(out, fil_[i]);
    }
}


void Tubule::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    size_t n = in.readUInt16();
    ObjectTag g;
    Object * w = sim.readReference(in, g);
    bone_ = Fiber::toFiber(w);
    for ( size_t i = 0; i < n-1; ++i )
    {
        w = sim.readReference(in, g);
        if ( i < NFIL )
            fil_[i] = Fiber::toFiber(w);
    }
    setFamily(bone_?bone_:fil_[0]);
}


void Tubule::report(std::ostream& os)
{
    os << reference() << " " << bone_->segmentation() << " " << bone_->nbPoints() << " " << bone_->length() << "\n";
    for ( size_t i = 0; i < NFIL; ++i )
        os << std::setw(7) << i << " " << fil_[i]->segmentation() << " " << fil_[i]->nbPoints() << " "<< fil_[i]->length() << "\n";
}

