// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "solid_set.h"
#include "solid_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "wrist.h"


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


Object * SolidSet::newObject(const ObjectTag tag, size_t num)
{
    if ( tag == Solid::TAG )
    {
        Property * p = simul_.properties.find("solid", num);
#ifdef BACKWARD_COMPATIBILITY
        // prior to 04.2016, "bead" and "solid" were used interchangeably
        if ( !p )
             p = simul_.properties.find("bead", num);
#endif
        if ( !p )
            throw InvalidIO("could not find `solid' class with ID "+std::to_string(num));
        return new Solid(static_cast<SolidProp*>(p));
   }
    std::cerr << "Warning: unknown Solid tag `"+std::string(1,tag)+"' requested\n";
    return nullptr;
}


/**
@ref Solid::build
 */
ObjectList SolidSet::newObjects(const std::string& name, Glossary& opt)
{
    SolidProp * p = simul_.findProperty<SolidProp>("solid", name);
    Solid * obj = new Solid(p);
    
    ObjectList res;
    res.push_back(obj);
    res.append(obj->build(opt, simul_));
    obj->fixShape();
    return res;
}


void SolidSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.put_line("\n#section "+title(), out.binary());
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

