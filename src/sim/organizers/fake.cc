// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "glossary.h"
#include "messages.h"
#include "mecapoint.h"
#include "simul.h"
#include "fake.h"
#include "aster.h"
#include "meca.h"


void Fake::step()
{
}


void Fake::setInteractions(Meca& meca) const
{
    assert_true( asterPoints.size() == solidPoints.size() );
    
    for ( size_t n = 0; n < asterPoints.size(); ++n )
        meca.addLink(asterPoints[n], solidPoints[n], prop->stiffness);
}


ObjectList Fake::build(Glossary& opt, Simul& sim)
{
    ObjectList objs;
    real rad = 0;
    if ( ! opt.set(rad, "radius") ||  rad <= 0 )
        throw InvalidParameter("fake:radius must be specified and > 0");

    assert_true(prop);
    
    Solid * as = nullptr, * bs = nullptr;
    std::string str;
    
    // find the Aster specified:
    if ( opt.set(str, "aster1") )
    {
        Aster * a = sim.organizers.pickAster(str);
        if ( !a)
            throw InvalidParameter("could not find Aster `"+str+"'");
        as = a->solid();
    }
    else if ( opt.set(str, "solid1") )
    {
        as = sim.pickSolid(str);
        if ( !as )
            throw InvalidParameter("could not find Solid `"+str+"'");
    }
    else
        throw InvalidParameter("fake:solid1 must be specified");

    // find the Aster specified:
    if ( opt.set(str, "aster2") )
    {
        Aster * a = sim.organizers.pickAster(str);
        if ( !a )
            throw InvalidParameter("could not find Aster `"+str+"'");
        bs = a->solid();
    }
    else if ( opt.set(str, "solid2") )
    {
        bs = sim.pickSolid(str);
        if ( !bs )
            throw InvalidParameter("could not find Solid `"+str+"'");
    }
    else
        throw InvalidParameter("fake:solid2 must be specified");
    
    if ( as == bs )
        throw InvalidParameter("fake cannot be made from only one Solid");
    
    Solid * so = new Solid(as->prop);
    
    Vector apos = as->posP(0);
    Vector bpos = bs->posP(0);
    
    // define two orthogonal directions:
    Vector dir1, dir2;
#if ( DIM >= 3 )
    Vector dir0 = normalize( apos - bpos );
    dir0.orthonormal(dir1, dir2);
#elif ( DIM == 2 )
    dir1 = ( apos - bpos ).orthogonal(1);
    dir2 = dir1;
#else
    throw InvalidParameter("the Fake class is not valid in 1D");
#endif
    
    asterPoints.clear();
    solidPoints.clear();

    for ( int d = 1; d < DIM; ++d )
    {
        for ( int s = -1; s < 2; s += 2 )
        {
            Vector x = s * ( d == 2 ? dir2 : dir1 );
            solidPoints.push_back(Mecapoint(so, so->addSphere(apos+x, rad)));
            asterPoints.push_back(Mecapoint(as, as->addPoint(apos+x)));
        
            solidPoints.push_back(Mecapoint(so, so->addSphere(bpos+x, rad)));
            asterPoints.push_back(Mecapoint(bs, bs->addPoint(bpos+x)));
        }
    }
    
    so->fixShape();
    as->fixShape();
    bs->fixShape();

    /*
     // print for debugging:
     so->write(std::clog, true);
     as->write(std::clog, true);
     bs->write(std::clog, true);
    */
    
    fkSolid = so;
    objs.push_back(so);
    objs.push_back(this);
    return objs;
}


Fake::~Fake()
{
    //Cytosim::log("destroying ", TAG, identity(), "\n");
    fkSolid = nullptr;
    prop = nullptr;
}


void Fake::write(Outputter& out) const
{
    writeMarker(out, Organizer::FAKE_TAG);
    Object::writeReference(out, fkSolid);
    writeOrganized(out);
}


void Fake::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    ObjectTag g;
#if BACKWARD_COMPATIBILITY < 53
    if ( in.formatID() < 53 )
    {
        size_t n = in.readUInt16();
        fkSolid = Solid::toSolid(sim.readReference(in, g));
        if ( n > 0 )
            Organizer::readOrganized(in, sim, n-1);
    }
    else
#endif
    {
        fkSolid = Solid::toSolid(sim.readReference(in, g));
        Organizer::read(in, sim, tag);
    }
}


/**
 This sets the ends of the link number `inx`
 or returns zero if the link does not exist
 */
bool Fake::getLink(index_t inx, Vector& pos1, Vector& pos2) const
{
    if ( inx < asterPoints.size() )
    {
        pos1 = asterPoints[inx].pos();
        pos2 = solidPoints[inx].pos();
        return true;
    }
    return false;
}

