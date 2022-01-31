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


ObjectList Tubule::build(real rad, Glossary& opt, Simul& sim)
{
    FiberProp const* fp = sim.findProperty<FiberProp>("fiber", prop->fiber_type);
    real len = 1.0, dev = 0;
    ObjectList res;

    // get the 'bone'
    if ( prop->bone_type.size() > 0 )
    {
        sim.fibers.newObjects(res, prop->bone_type, opt);
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
        fib->constrainLength(false);
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
    
    // set orthonormal coordinate system (dir, E, F)
    Vector E(0,rad,0), F(0,0,rad);
    dir.orthonormal(E, F, rad);
    
    /*
     translate protofilaments to form a tube, arranging them
     in a counter-clockwise order, to form a right-handed helix
     */
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
        if ( fil_[i] )
        {
            fil_[i]->family_ = fam;
            fil_[i]->sister_ = fil_[(i+NFIL-1)%NFIL];
            fil_[i]->brother_ = fil_[(i+1)%NFIL];
        }
    }
    if ( bone_ )
        bone_->family_ = fam;
#else
    std::cerr << "ERROR: to use Tubule, please compile with FIBER_HAS_FAMILY\n";
    exit(1);
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
void Tubule::setInteractionsA(Meca& meca) const
{
#if ( DIM >= 3 )
    const real stiff = prop->stiffness[0];
    const real ang = M_PI / NFIL; // this is half the angle of each sector!!
    const real len = prop->radius * 2 * std::sin(ang);  // distance between protofilaments
    const real C = std::cos(ang), S = std::sin(ang);
    
    assert_true(fil_[0]);
    const size_t e = fil_[0]->nbSegments();
    
    Rotation mat(0,1);

    for ( size_t i = 0; i < e; ++i )
    {
        // compute centerline
        Vector cen(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            cen += fil_[n]->posPoint(i);
        cen /= NFIL;
        
        // get average direction of the Tubule at this location:
        Vector dir(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            dir += fil_[n]->diffPoints(i);
        dir.normalize();
        mat = Rotation::rotationAroundAxis(dir, C, S);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            real alpha = len * fil_[n]->segmentationInv();
            Vector leg = mat.vecmul(( cen - fil_[n]->posPoint(i) ).normalized(alpha));
            meca.addSideLink(FiberSegment(fil_[n],i), 0, Mecapoint(fil_[n+1],i), leg, stiff);
        }
    }

    // get centerline
    Vector cen(0,0,0);
    for ( size_t n = 0; n < NFIL; ++n )
        cen += fil_[n]->posPoint(e);
    cen /= NFIL;
    
    for ( size_t n = 0; n < NFIL; ++n )
    {
        real alpha = len * fil_[n]->segmentationInv();
        Vector leg = mat.vecmul(( cen - fil_[n]->posPoint(e) ).normalized(alpha));
        meca.addSideLink(FiberSegment(fil_[n],e-1), 1, Mecapoint(fil_[n+1],e), leg, stiff);
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
    const real C = std::cos(ang), S = std::sin(ang);
    
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
    real beta = prop->radius / len;
    Vector cen, dir;
    
    for ( size_t i = 0; i < e; ++i )
    {
        // get centerline
        cen = bone_->posPoint(i);
        dir = bone_->dirSegment(i);
        mat = Rotation::rotationAroundAxis(dir, C, S);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            real alpha = len * fil_[n]->segmentationInv();
            Vector leg = ( cen - fil_[n]->posPoint(i) ).normalized(alpha);
            // orthoradial beams:
            meca.addSideLink(FiberSegment(fil_[n],i), 0, Mecapoint(fil_[n+1],i), mat.vecmul(leg), stiffL);
            // radial spoke:
            meca.addSideLink(FiberSegment(fil_[n],i), 0, Mecapoint(bone_,i), beta*cross(dir,leg), stiffR);
            // twist stiffness
            meca.addTorque(Mecapoint(fil_[n],i), Mecapoint(fil_[n+1],i), stiffR);
        }
    }

    // link last point:
    cen = bone_->posPoint(e);
    for ( size_t n = 0; n < NFIL; ++n )
    {
        real alpha = len * fil_[n]->segmentationInv();
        Vector leg = ( cen - fil_[n]->posPoint(e) ).normalized(alpha);
        meca.addSideLink(FiberSegment(fil_[n],e-1), 1, Mecapoint(fil_[n+1],e), mat.vecmul(leg), stiffL);
        meca.addSideLink(FiberSegment(fil_[n],e-1), 1, Mecapoint(bone_,e), beta*cross(dir,leg), stiffR);
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
    const real theta = 2 * M_PI / NFIL;
    const real len = 2 * prop->radius * std::sin(0.5*theta);  // distance between protofilaments
    const real stiffL = prop->stiffness[0];
    const real stiffA = prop->stiffness[1];
    const real stiffT = prop->stiffness[1];
    Vector2 ang(std::cos(theta), std::sin(theta));
    
    assert_true(fil_[0]);
    const size_t e = fil_[0]->nbSegments();
    
    MatrixBlock mat;
    for ( size_t i = 0; i < e; ++i )
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
        
        // create rotation matrix for torque:
        mat = Meca::torqueMatrix(stiffA, dir, ang);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            real alpha = len * fil_[n]->segmentationInv();
            Vector leg = ( 2*cen - fil_[n]->posPoint(i) - fil_[n+1]->posPoint(i)).normalized(alpha);
            meca.addSideLink(FiberSegment(fil_[n],i), 0, Mecapoint(fil_[n+1],i), leg, stiffL);
            meca.addTorque(Mecapoint(fil_[n],i), Mecapoint(fil_[n+1],i), Mecapoint(fil_[n+2],i), mat, stiffA);
            meca.addTorque(Mecapoint(fil_[n],i), Mecapoint(fil_[n+1],i), stiffT);
        }
    }
    
    // last segment:
    Vector cen(0,0,0);
    for ( size_t n = 0; n < NFIL; ++n )
        cen += fil_[n]->posPoint(e);
    cen /= NFIL;

    for ( size_t n = 0; n < NFIL; ++n )
    {
        real alpha = len * fil_[n]->segmentationInv();
        Vector leg = (2*cen - fil_[n]->posPoint(e) - fil_[n+1]->posPoint(e)).normalized(alpha);
        meca.addSideLink(FiberSegment(fil_[n],e-1), 1, Mecapoint(fil_[n+1],e), leg, stiffL);
    }
#endif
}

/**
 Using addSideLink() but this causes a serious problem, since Torque are disbalanced...
 Strasbourg, 13.05.2021
 */
void Tubule::setInteractionsD(Meca& meca) const
{
#if ( DIM >= 3 )
    const real ang = M_PI / NFIL;  // half of the sector angle!
    const real len = 2 * prop->radius * std::sin(ang);  // distance between protofilaments
    const real stiffL = prop->stiffness[0];
    const real stiffR = prop->stiffness[1];
    const real C = std::sin(ang), S = std::cos(ang); // rotation of angle PI/2 - ang
    Rotation mat;

    for ( size_t n = 0; n < NFIL; ++n )
    {
        Fiber * fib = fil_[n+1];
        const size_t e = fil_[n]->nbSegments();
        const real alpha = len * fib->segmentationInv();

        if ( fib->nbSegments() != e )
        {
            Cytosim::warn << "unequal Tubule filaments\n";
            return;
        }
        
        for ( size_t i = 0; i < e; ++i )
        {
            mat = Rotation::rotationAroundAxis(fib->dirSegment(i), C, S);
            // rotate direction of previous protofilament pair:
            Vector leg = mat.vecmul(fib->posPoint(i) - fil_[n]->posPoint(i)).normalized(alpha);
            meca.addSideLink(FiberSegment(fib,i), 0, Mecapoint(fil_[n+2],i), leg, stiffL);
            // orthoradial links with bending stiffness
            //meca.addTorque(Mecapoint(fil_[n],i), Mecapoint(fil_[n+1],i), Mecapoint(fil_[n+2],i), 1.0, stiffL);
            // twist stiffness
            meca.addTorque(Mecapoint(fib,i), Mecapoint(fil_[n+2],i), stiffR);
        }

        // use the same rotation matrix for the last point:
        //mat = Rotation::rotationAroundAxis(fib->dirSegment(e-1), C, S);
        Vector leg = mat.vecmul(fib->posPoint(e) - fil_[n]->posPoint(e)).normalized(alpha);
        meca.addSideLink(FiberSegment(fib,e-1), 1, Mecapoint(fil_[n+2],e), leg, stiffL);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - I/O


void Tubule::write(Outputter& out) const
{
    out.writeUInt16(NFIL+1);
    Object::writeReference(out, bone_);
    for ( size_t i = 0; i < NFIL; ++i )
    {
        Object::writeReference(out, fil_[i]);
    }
}


void Tubule::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    
    ObjectTag g;
    size_t n = in.readUInt16()-1;
    Object * w = sim.readReference(in, g);
    bone_ = Fiber::toFiber(w);
    size_t i = 0;
    for (; i < n && i < NFIL; ++i )
    {
        w = sim.readReference(in, g);
        fil_[i] = Fiber::toFiber(w);
    }
    for (; i < NFIL; ++i )
        fil_[i] = nullptr;
    setFamily(bone_?bone_:fil_[0]);
}


void Tubule::report(std::ostream& os) const
{
    os << reference() << " " << bone_->segmentation() << " " << bone_->nbPoints() << " " << bone_->length() << "\n";
    for ( size_t i = 0; i < NFIL; ++i )
    {
        Fiber const* F = fil_[i];
        if ( F )
            os << std::setw(7) << i << " " << F->segmentation() << " " << F->nbPoints() << " "<< F->length() << "\n";
    }
}

