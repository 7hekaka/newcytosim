// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#include "solid_set.h"
#include "solid_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


void SolidSet::step()
{
#if ( 0 )
    for ( Solid * B = first(); B; B = B->next() )
        B->step();
#endif
#if NEW_SOLID_MAKE_COUPLE
    static float nextCreation = RNG.exponential();
    for ( Solid * B = first(); B; B = B->next() )
    {
        nextCreation -= B->prop->source_rate_dt;
        while ( nextCreation <= 0 )
        {
            nextCreation += RNG.exponential();
            // we only consider creation on the first point!
            Vector pos = B->posPoint(0) + Vector::randB(B->radius(0));
            if ( B->prop->source_couple )
            {
                Couple * C = B->prop->source_couple->newCouple(pos);
                simul_.couples.add(C);
                C->activate();
            }
            else
            {
                Single * S = B->prop->source_single->newSingle(pos);
                simul_.singles.add(S);
            }
        }
    }
#endif
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
            for ( size_t p = 0; p < S->nbPoints(); ++p )
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
    return obj->build(opt, simul_);
}


void SolidSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writeObjects(out, pool_);
#if NEW_SOLID_HAS_TWIN
        for ( Object const* n=pool_.front(); n; n=n->next() )
            static_cast<Solid const*>(n)->writeTwin(out);
#endif
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

