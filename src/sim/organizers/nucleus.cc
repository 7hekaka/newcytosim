// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "nucleus.h"
#include "exceptions.h"
#include "sphere_prop.h"
#include "bundle_prop.h"
#include "mecapoint.h"
#include "fiber_set.h"
#include "glossary.h"
#include "bundle.h"
#include "simul.h"
#include "meca.h"


void Nucleus::step()
{
}


void Nucleus::setInteractions(Meca& meca) const
{
    Sphere const* sph = sphere();
    
    if ( sph )
    {
        for ( size_t i = Sphere::nbRefPoints; i < nuSphere->nbPoints(); ++i )
        {
            Fiber const* fib = fiber(i-Sphere::nbRefPoints);
            if ( fib )
                meca.addLink(Mecapoint(sph, i), fib->exactEnd(MINUS_END), prop->stiffness);
        }
    }
}


//------------------------------------------------------------------------------
void Nucleus::build(ObjectList& res, Glossary& opt, Simul& sim)
{
    std::string str, spec;
    assert_true(prop);
    size_t cnt = 0;
    
    real rad = -1;
    if ( !opt.set(rad, "radius" ) || rad <= 0 )
        throw InvalidParameter("parameter `radius` should be specified and > 0");
   
    if ( !opt.set(str, "sphere") )
        throw InvalidParameter("parameter `sphere` should be specified");

    SphereProp * sp = sim.findProperty<SphereProp>("sphere", str);
    nuSphere = new Sphere(sp, rad);
    res.push_back(nuSphere);
    
    // get the center of the sphere
    Vector c = nuSphere->posP(0);
    
    if ( opt.set(cnt, "fibers") && cnt > 0 )
    {
        if ( !opt.set(str, "fibers", 1) )
            throw InvalidParameter("fiber type (fiber[1]) should be specified");
        opt.set(spec, "fibers", 2);

        // create points and clamps and add fiber attached to them
        for ( size_t ii = 0; ii < cnt; ++ii )
        {
            Glossary fiber_opt(spec);
            Fiber * fib = sim.fibers.newFiber(res, str, fiber_opt);
            Vector pos = c + Vector::randU(rad);
            Vector dir = Vector::randU();
            fib->setStraight(pos, dir, fib->length());
            nuSphere->addPoint(pos);
            grasp(fib);
        }
    }
    
    if ( opt.set(cnt, "bundles") && cnt > 0 )
    {
        if ( !opt.set(str, "bundles", 1) )
            throw InvalidParameter("bundle type (bundles[1]) should be specified");
        opt.set(spec, "bundles", 2);

        BundleProp * bp = sim.findProperty<BundleProp>("bundle", str);

        Rotation rot;
        // add bundles
        ObjectList objs;
        const real len = 0.5 * bp->overlap;
        for ( size_t ii = 0; ii < cnt; ++ii  )
        {
            Glossary bundle_opt(spec);
            rot = Rotation::randomRotation();
            //a random position on the sphere:
            Vector pos = rot * Vector(0, rad, 0);
            //a vector tangent to the sphere:
            Vector dir = rot * Vector(len, 0, 0);
            
            objs.clear();
            Bundle * bun = new Bundle(bp);
            bun->build(objs, bundle_opt, sim);
            res.append(objs);
            res.push_back(bun);

            //position the bundle (initially aligned with X) tangentially:
            ObjectSet::moveObjects(objs, Isometry(pos, rot));
            
            nuSphere->addPoint( c + (pos-dir).normalized(rad) );
            grasp(bun->organized(0));
            
            nuSphere->addPoint( c + (pos+dir).normalized(rad) );
            grasp(bun->organized(1));
        }
    }
}


Nucleus::~Nucleus()
{
    //Cytosim::log("destroying %c%lu\n", TAG, identity());
    nuSphere = nullptr;
    prop = nullptr;
}

//------------------------------------------------------------------------------

void Nucleus::write(Outputter& out) const
{
    writeHeader(out, Organizer::TAG_NUCLEUS);
    Object::writeReference(out, nuSphere);
    writeOrganized(out);
}


void Nucleus::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    ObjectTag g;
#if BACKWARD_COMPATIBILITY < 53
    if ( in.formatID() < 53 )
    {
        size_t n = in.readUInt16();
        nuSphere = Sphere::toSphere(sim.readReference(in, g));
        Organizer::readOrganized(in, sim, n-1);
    }
    else
#endif
    {
        nuSphere = Sphere::toSphere(sim.readReference(in, g));
        Organizer::read(in, sim, tag);
    }
    assert_true( nbOrganized() > 0 );
}


//------------------------------------------------------------------------------
#pragma mark - Display

/**
 This sets the ends of the link number `inx`
 or returns zero if the link does not exist
 */
bool Nucleus::getLink(size_t inx, Vector& pos1, Vector& pos2) const
{
    size_t i = inx + Sphere::nbRefPoints;
    if ( sphere() && i < sphere()->nbPoints() )
    {
        pos1 = sphere()->posP(i);
        
        Fiber const* fib = fiber(inx);
        if ( fib )
            pos2 = fib->posEnd(MINUS_END);
        else
            pos2 = pos1;
        return true;
    }
    return false;
}

