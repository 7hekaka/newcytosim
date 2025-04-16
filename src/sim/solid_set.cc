// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "solid_set.h"
#include "solid_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


void SolidSet::steps()
{
#if ( 0 )
    Solid * obj = first();
    while ( obj )
    {
        Solid * nxt = obj->next();
        obj->step();
        obj = nxt;
    }
#endif
    if ( size() > 1 ) shuffle();
}



/**
 This performs a systematic search over all Solids, to return the closest
 possible match.
 
 //@todo: could pick one of the matching Sphere randomly
 */
Solid* SolidSet::insideSphere(Vector const& pos, real range, size_t& inx, SolidProp const* sel) const
{
    real best = INFINITY;
    Solid* res = nullptr;
    
    for ( Solid* S = first(); S; S = S->next() )
        if ( S->prop == sel )
        {
            for ( index_t p = 0; p < S->nbPoints(); ++p )
            {
                const real rad = S->radius(p);
                if ( rad > 0 )
                {
                    real dd = distanceSqr(S->posPoint(p), pos);
                    if (( dd < best ) && ( dd < square(rad+range)))
                    {
                        best = dd;
                        res = S;
                        inx = p;
                    }
                }
            }
        }
    
    return res;
}


//------------------------------------------------------------------------------

Property* SolidSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "solid" )
        return new SolidProp(cat, nom);
    return nullptr;
}


Object * SolidSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Solid::TAG )
    {
        Property * p = simul_.properties.find("solid", pid);
#if BACKWARD_COMPATIBILITY < 47
        // prior to 04.2016, "bead" and "solid" were used interchangeably
        if ( !p )
             p = simul_.properties.find("bead", pid);
#endif
        if ( !p )
            throw InvalidIO("undefined `solid' class with ID "+std::to_string(pid));
        return new Solid(static_cast<SolidProp*>(p));
    }
    throw InvalidIO("Warning: unknown Solid tag `"+std::to_string(tag)+"'");
    return nullptr;
}


/**
@ref Solid::build
 */
ObjectList SolidSet::newObjects(Property const* p, Glossary& opt)
{
    SolidProp const* pp = static_cast<SolidProp const*>(p);
    Solid * obj = new Solid(pp);
    ObjectList const&& list = obj->build(opt, simul_);
    obj->fixShape();
    return list;
}


void SolidSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}


void SolidSet::defrostMore()
{
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        //std::clog << "delete " << i->reference() << "\n";
        inventory_.unassign(i);
        i->objset(nullptr);
        simul_.singles.deleteWrists(i);
        delete(i);
    }
}


void SolidSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul_.singles.deleteWrists(obj);
}


void SolidSet::foldPositions(Modulo const* m) const
{
    for ( Solid * o=first(); o; o=o->next() )
        o->foldPosition(m);
}

