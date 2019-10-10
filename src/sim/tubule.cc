// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#include "dim.h"
#include "tubule.h"
#include "fiber.h"
#include "simul.h"
#include "meca.h"


void Tubule::reset()
{
    for ( size_t i = 0; i < NFIL+2; ++i )
        fil_[i] = nullptr;
}


ObjectList Tubule::build(Glossary& opt, Simul& sim)
{
    ObjectList res;
    for ( size_t i = 0; i < NFIL; ++i )
    {
        ObjectList objs = sim.fibers.newObjects(prop->fiber_type, opt);
        res.append(objs);
        fil_[i] = Fiber::toFiber(objs[0]);
    }
    // wrap array values for convenience
    fil_[NFIL  ] = fil_[0];
    fil_[NFIL+1] = fil_[1];
    
#if ( DIM >= 3 )
    Vector E(0,tube_radius,0), F(0,0,tube_radius);
    {
        // find average direction
        Vector dir(0,0,0);
        for ( size_t i = 0; i < NFIL; ++i )
            dir += fil_[i]->diffPoints(0);
        dir.normalize();
        dir.orthonormal(E, F);
        E *= tube_radius;
        F *= tube_radius;
    }
    // adjust protofilaments to form a tube:
    real a = M_PI * RNG.sreal();
    real da = 2 * M_PI / NFIL;
    for ( size_t i = 0; i < NFIL; ++i )
    {
        fil_[i]->translate(cos(a)*E+sin(a)*F);
        a += da;
    }
#endif
    
    // set as left-handed helix:
    assert_true(signature());
    for ( size_t i = 0; i < NFIL; ++i )
    {
        fil_[i]->setOrigin(i*(-0.012/NFIL));
#if FIBER_HAS_FAMILY
        fil_[i]->family = signature();
#endif
    }

    return res;
}


Tubule::Tubule(TubuleProp * p) : prop(p)
{
    reset();
}


Tubule::~Tubule()
{
    prop = nullptr;
}


void Tubule::step(Simul&)
{
}

    
///
void Tubule::setInteractions(Meca& meca)
{
    const real stiffL = prop->stiffness[0];
    const real stiffA = prop->stiffness[1];
    const real len = fil_offset;  // distance between protofilaments
#if ( DIM >= 3 )
    const real ang = 2 * M_PI / NFIL;
    real co = cos(ang), si = sin(ang);
    
    assert_true(fil_[0]);
    const size_t end = fil_[0]->nbPoints() - 1;
    
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
                meca.addSideLink3D(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), arm, stiffL);
            }
        }
        else
        {
            for ( size_t n = 0; n < NFIL; ++n )
            {
                Vector arm = (2*cen - fil_[n]->posPoint(i) - fil_[n+1]->posPoint(i)).normalized(len);
                meca.addSideLink3D(Interpolation(fil_[n],i-1,i,1), Mecapoint(fil_[n+1],i), arm, stiffL);
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


void Tubule::write(Outputter& out) const
{
}


void Tubule::read(Inputter& in, Simul& sim, ObjectTag tag)
{
}

