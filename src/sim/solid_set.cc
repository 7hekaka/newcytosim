// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "solid_set.h"
#include "solid_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


#if ( 0 )
void SolidSet::step()
{
    for ( Solid * o = first(); o; o=o->next() )
        o->step();
}
#endif


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
void SolidSet::newObjects(ObjectList& res, const std::string& name, Glossary& opt)
{
    SolidProp * p = simul_.findProperty<SolidProp>("solid", name);
    Solid * obj = new Solid(p);
    
    obj->build(res, opt, simul_);
    res.push_back(obj);
    obj->fixShape();
}


void SolidSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.writeLine("\n#section "+title());
        writeObjects(out, pool_);
    }
}

//------------------------------------------------------------------------------


void SolidSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul_.singles.removeWrists(obj);
}


void SolidSet::foldPositions(Modulo const* m) const
{
    for ( Solid * o=first(); o; o=o->next() )
        o->foldPosition(m);
}

