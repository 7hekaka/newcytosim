// Cytosim was created by Francois Nedelec. Copyright 2024

#include "fiber_set.h"
#include "fiber_segment.h"
#include "tokenizer.h"
#include "iowrapper.h"
#include "messages.h"
#include "glossary.h"
#include "fiber_prop.h"
#include "growing_fiber_prop.h"
#include "dynamic_fiber_prop.h"
#include "classic_fiber_prop.h"
#include "treadmilling_fiber_prop.h"
#include "lapack.h"
#include "simul.h"
#include "cymdef.h"

//------------------------------------------------------------------------------

/**
 @defgroup FiberGroup Fiber and related
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief The Fiber has fixed length, but derived classes can change length.

 A fiber is a filament of constant length.
 Derived classes are available, where different models of how length may change
 have been implemented.
 
 List of classes accessible by specifying `fiber:activity`.
 
 `activity`        | Class               | Parameters                 |
 ------------------|---------------------|-----------------------------
 `none` (default)  | Fiber               | @ref FiberPar
 `grow`            | GrowingFiber        | @ref GrowingFiberPar
 `classic`         | ClassicFiber        | @ref ClassicFiberPar
 `dynamic`         | DynamicFiber        | @ref DynamicFiberPar
 `treadmill`       | TreadmillingFiber   | @ref TreadmillingFiberPar
 
 */
Property* FiberSet::newProperty(const std::string& cat, const std::string& nom, Glossary& opt) const
{
    if ( cat == "fiber" )
    {
        std::string a;
        if ( opt.peek(a, "activity") )
        {
            if ( a == "classic" )
                return new ClassicFiberProp(nom);
            if ( a == "grow" )
                return new GrowingFiberProp(nom);
            if ( a == "dynamic" )
                return new DynamicFiberProp(nom);
            if ( a == "treadmill" )
                return new TreadmillingFiberProp(nom);
            if ( a == "none" )
                return new FiberProp(nom);

            std::cerr << "INCIDENT: substituting generic Fiber for unknown `"+a+"'\n";
            return new FiberProp(nom);
            //throw InvalidParameter("unknown fiber:activity `"+a+"'");
        }
        return new FiberProp(nom);
    }
    return nullptr;
}


/**
 Many options depend on the type of fiber: Fiber, DynamicFiber, ClassicFiber, etc.
 Here only the common options are described.

 <hr>
 
 You may directly attach Single or Couple to the fiber, in different ways:
 
     new filament
     {
        attach1 = [NUMBER] NAME
     }
 
 `NAME` should designate the Single or Couple that will be attached to the Fiber.
 `NUMBER` will specify how many Single/Couple will be attached (by default: 1).
 Note that for a Couple, the first Hand is attached to the fiber (and not the second).
 In this case, the Single/Couple are anchored at random position distributed 
 uniformly along the Fiber.

     new filament
     {
        attach1 = [NUMBER] NAME, ABSCISSA, REFERENCE [, MODIFIER [, POSITION]]
     }
 
 If `ABSCISSA` is specified (in micrometers), the Single/Couple will be attached
 at the specified distance from the `REFERENCE = { minus_end, plus_end, center }
 (default is `minus_end`). The distance is counted towards the other end.
 
 A `MODIFIER = { uniform, exponential, regular }` can be specified (default: `uniform`):
     - with `uniform` Single/Couple are uniformly distributed over distance `[0, DISTANCE]`
       from the `REFERENCE`.
     - with `regular`, they are distributed regularly.
     .
 Finally, `POSITION` can be specified. This is mostly relevant for Singles with
 activity `fixed`.
 
 Multiple attachement instructions can be specified as `attach1`, `attach2`, etc.
 For example, this attaches one `simplex` at each end of the filaments:
 
     new filament
     {
        attach1 = simplex, 0.0, minus_end
        attach2 = simplex, 0.0, plus_end
     }

 @}
 */
Fiber * FiberSet::newFiber(ObjectList& objs, FiberProp const* fip, Glossary& opt) const
{
    Fiber * fib = fip->newFiber(opt);
    fib->birthTime(simul_.time());
    objs.push_back(fib);
    
#if FIBER_HAS_FAMILY
    std::string str;
    if ( opt.set(str, "family") )
        fib->family_ = simul_.pickFiber(str);
#endif
 
    size_t inp = 1;
    std::string spe, var = "attach1";
#if BACKWARD_COMPATIBILITY <= 50
    if ( opt.has_key("attach") )
    {
        var = "attach";
        inp = 0;
    }
#endif
    //can add Singles or Couples to the Fiber:
    while ( opt.set(spe, var) )
    {
        size_t cnt = 1;
        Tokenizer::split_integer(cnt, spe);
        fib->makeAttachedHands(objs, spe, cnt, opt, var, simul_);
        var = "attach" + std::to_string(++inp);
    }
    return fib;
}


Fiber * FiberSet::newFiber(ObjectList& objs, FiberProp const* fip, std::string const& spec) const
{
    Glossary opt(spec);
    Fiber * F = newFiber(objs, fip, opt);
    std::string war;
    if ( opt.has_warning(war) )
    {
        //print_yellow(stderr, war);
        std::cerr << war << " in `" << spec << "'\n";
    }
    return F;
}


/**
 The returned object is not initialized, since this is used for file input
 */
Object * FiberSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Fiber::TAG || tag == Fiber::COMPACT_TAG )
    {
        FiberProp const* fp = simul_.findProperty<FiberProp>("fiber", pid);
        Fiber * obj = fp->newFiber();
        obj->birthTime(simul_.time());
        return obj;
    }
    throw InvalidIO("Warning: unknown Fiber tag `"+std::to_string(tag)+"'");
    return nullptr;
}


ObjectList FiberSet::newObjects(Property const* p, Glossary& opt)
{
    ObjectList res(4);
    newFiber(res, static_cast<FiberProp const*>(p), opt);
    return res;
}


void FiberSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}


void FiberSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << '\n' << title();
        PropertyList plist = simul_.properties.find_all(title());
        for ( Property const* i : plist )
        {
            FiberProp const* p = static_cast<FiberProp const*>(i);
            size_t cnt = count(p);
            os << '\n' << std::setw(10) << cnt << ' ' << p->name();
            os << " ( " << p->activity << " )";
        }
        if ( plist.size() > 1 )
            os << '\n' << std::setw(10) << size() << " total";
    }
}


//------------------------------------------------------------------------------
#pragma mark - Step

/**
 Calculate the free monomer concentration. 
 Calls step() once for every Fiber.
 */

void FiberSet::steps()
{
    PropertyList plist = simul_.properties.find_all("fiber");
    
    // calculate the total length used for each kind of Fiber:
    for ( Property * i : plist )
    {
        static_cast<FiberProp*>(i)->used_polymer = 0;
        static_cast<FiberProp*>(i)->fiber_count = 0;
    }

    for ( Fiber const* fib=first(); fib; fib=fib->next() )
    {
        assert_false(fib->bad());
        fib->prop->used_polymer += fib->length();
        ++fib->prop->fiber_count;
    }
    
    // calculate the ratio of free polymer for each class of Fiber:
    for ( Property * i : plist )
    {
        FiberProp * fp = static_cast<FiberProp*>(i);

        // update the normalized monomer concentration:
        fp->free_polymer = 1.0 - fp->used_polymer / fp->total_polymer;
        
        if ( fp->free_polymer < 0 )
        {
            // this may happen with a fast grow_speed / large time_step
            Cytosim::warn("Inconsistent monomer pool for ", fp->name(),
                          ": total_polymer=", fp->total_polymer,
                          ", used_polymer= ", fp->used_polymer, "\n");
            fp->free_polymer = 0;
        }
    }
    /*
     We call step() here exactly once for every Fiber.
     New Fiber may be created, for instance by Fiber::severSoon(), but they should
     be linked at the start of the list, and thus not considered here.
     */
    Fiber * obj = first();

    while ( obj )
    {
        Fiber * nxt = obj->next();
        obj->step();
        // delete object that have been flagged
        if ( ! obj->prop ) eraseObject(obj);
        obj = nxt;
    }
    if ( size() > 1 ) shuffle();
}


/**
 Cut all Fibers along the plane defined by `n.pos + a = 0`.
 - new plus ends are set to state `stateP`
 - new plus ends are set to state `stateM`
 - any fragment shorter than `min_length` is deleted
 */
void FiberSet::planarCut(Vector const& n, const real a,
                         state_t stateP, state_t stateM, real min_len)
{
    /*
     We must ensure here that each Fiber is processed only once.
     This code works if newly created Fiber are linked at the head of the list
     */
    Fiber * obj = first();

    while ( obj )
    {
        Fiber * nxt = obj->next();
        obj->planarCut(n, a, stateP, stateM, min_len);
        if ( obj->prop == nullptr )
            eraseObject(obj);
        obj = nxt;
    }
}

/**
 Cut all Fibers in list `objs` along the plane defined by `n.pos + a = 0`.
 - new plus ends are set to state `stateP`
 - new plus ends are set to state `stateM`
 - any fragment shorter than `min_length` is deleted
 */
void FiberSet::planarCut(ObjectList& objs, Vector const& n, const real a,
                         state_t stateP, state_t stateM, real min_len)
{
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            fib->planarCut(n, a, stateP, stateM, min_len);
            if ( i->property() == nullptr )
                eraseObject(i);
        }
    }
}


void FiberSet::foldPositions(Modulo const* m) const
{
    for ( Fiber * o=first(); o; o=o->next() )
        o->foldPosition(m);
}


/**
 Calculate intersection between all fibers,
 and report the corresponding abscissa in arrays 'res1' and 'res2'.
 */
void FiberSet::allIntersections0(FiberSiteList& res1, FiberSiteList& res2,
                                 const real max_distance) const
{
    const real sup = square(max_distance);
    res1.clear();
    res2.clear();

    for ( Fiber * fib1 = first(); fib1; fib1 = fib1->next() )
    {
        for ( index_t s1 = 0; s1 < fib1->nbSegments(); ++s1 )
        {
            FiberSegment seg(fib1, s1);
            // check against other segments of this fiber
            for ( index_t s2 = s1+2; s2 < fib1->nbSegments(); ++s2 )
            {
                FiberSegment soc(fib1, s2);
                real abs1, abs2;
                real dis2 = seg.shortestDistanceSqr(soc, abs1, abs2);
                if ((dis2 < sup) && seg.within(abs1) && soc.within(abs2))
                {
                    res1.emplace(fib1, abs1+fib1->abscissaPoint(s1));
                    res2.emplace(fib1, abs2+fib1->abscissaPoint(s2));
                }
            }
            // check against other fibers:
            for ( Fiber * fib2 = fib1->next(); fib2; fib2 = fib2->next() )
            {
                for ( index_t s2 = 0; s2 < fib2->nbSegments(); ++s2 )
                {
                    FiberSegment soc(fib2, s2);
                    real abs1, abs2;
                    real dis2 = seg.shortestDistanceSqr(soc, abs1, abs2);
                    if ((dis2 < sup) && seg.within(abs1) && soc.within(abs2))
                    {
                        res1.emplace(fib1, abs1+fib1->abscissaPoint(s1));
                        res2.emplace(fib2, abs2+fib2->abscissaPoint(s2));
                    }
                }
            }
        }
    }
}


/**
 Calculate intersection between all fibers,
 and report the corresponding abscissa in arrays 'res1' and 'res2'.
 This version is using the FiberGrid for faster results
 */
void FiberSet::allIntersections(FiberSiteList& res1, FiberSiteList& res2,
                                const real max_distance) const
{
    FiberGrid grid;
    Space const* spc = simul_.spaces.master();

    if ( !spc )
        return allIntersections0(res1, res2, max_distance);

    real grid_step = 1;
    Simul::setFiberGrid(grid, spc, grid_step);

#if ( 0 )
    // check what other method gives:
    allIntersections0(res1, res2, max_distance);
    std::clog << "FiberSet::allIntersections0() found " << res1.size() << " intersections\n";
    for ( size_t i = 0; i < res1.size(); ++i )
        std::clog << res1[i] << " " << res2[i] << "\n";
#endif
    
    // find largest fiber:segmentation
    real len = 0;
    for ( Fiber * fib = first(); fib; fib = fib->next() )
        len = std::max(len, fib->segmentation());
    
    const real sup = square(max_distance);
    res1.clear();
    res2.clear();
    
    // distribute segments
    grid.paintGrid(first(), nullptr, std::sqrt(square(len)+sup));
    
    FiberGrid::SegmentList list;
    for ( Fiber const* fib = first(); fib; fib = fib->next() )
    {
        //std::clog << fib->reference() << ":\n";
        for ( index_t s = 0; s < fib->nbSegments(); ++s )
        {
            FiberSegment seg(fib, s);
            list = grid.cellTargets(seg.middle());
            //std::clog << seg << ":";
            for ( FiberSegment const& soc : list )
            {
                Fiber const* bif = soc.fiber();
                if ( fib < bif )
                {
                    //std::clog << "   " << can;
                    real abs1, abs2;
                    real dis2 = seg.shortestDistanceSqr(soc, abs1, abs2);
                    if ((dis2 < sup) && seg.within(abs1) && soc.within(abs2))
                    {
                        res1.emplace(fib, abs1+fib->abscissaPoint(s));
                        res2.emplace(bif, abs2+bif->abscissaPoint(soc.point()));
                    }
                }
            }
            //std::clog << '\n';
        }
    }
#if ( 0 )
    // detailed debug output
    std::clog << "FiberSet::allIntersections()  found " << res1.size() << " intersections \n";
    for ( size_t i = 0; i < res1.size(); ++i )
        std::clog << res1[i] << " " << res2[i] << "\n";
#endif
}


/**
 Set a list of Locations on the fibers, chosen randomly with uniform sampling.
 The number of sites returned on a section of length `L` is  `L / gap`.
 `gap` is thus the average distance between sites.
 
 Because the list of fiber is regularly shuffled, the sites will consider fibers
 in a random order. However, the sites on one fiber will be listed in the order
 of increasing abscissa.
 
 Condition: ( gap > 0 )
 */
void FiberSet::uniFiberSites(FiberSiteList& res, const real gap) const
{
    assert_true( gap > 0 );

    res.clear();
    Fiber * fib = first();
    real abs = gap * RNG.exponential();
    while ( fib )
    {
        real len = fib->length();
        while ( abs < len )
        {
            res.emplace(fib, abs+fib->abscissaM());
            abs += gap * RNG.exponential();
        }
        abs -= len;
        fib = fib->next();
    }
}


/// a random site on the fiber, equidistributed over length
/**
 This method is unefficient if multiple sites are desired:
 It requires two passes to first add the lengths of the fibers
 In principle, one pass would be sufficient, using one random number per fiber
 */
FiberSite FiberSet::randomSite() const
{
    real abs = 0;
    for ( Fiber const* fib=first(); fib; fib=fib->next() )
        abs += fib->length();

    assert_true( abs > 0 );
    
    abs *= RNG.preal();

    for ( Fiber* fib=first(); fib; fib=fib->next() )
    {
        real len = fib->length();
        if ( abs <= len )
            return FiberSite(fib, fib->abscissaM()+abs);
        abs -= len;
    }
    
    ABORT_NOW("unexpected abscissa overrun");
    return FiberSite(first(), 0);
}


/// a random site on the fibers of class 'prop'
/**
 This method is unefficient if multiple sites are desired
 It requires two passes to first add the lengths of the fibers
 */
FiberSite FiberSet::randomSite(FiberProp const* sel) const
{
    real abs = 0;
    for ( Fiber const* fib=first(); fib; fib=fib->next() )
        if ( fib->property() == sel )
            abs += fib->length();

    if ( abs == 0 )
        throw InvalidParameter("found no fibers of requested class");

    abs *= RNG.preal();
    
    for ( Fiber* fib=first(); fib; fib=fib->next() )
        if ( fib->property() == sel )
        {
            real len = fib->length();
            if ( abs <= len )
                return FiberSite(fib, fib->abscissaM()+abs);
            abs -= len;
        }
    
    ABORT_NOW("unexpected abscissa overrun");
    return FiberSite(first(), 0);
}


/**
 Returns a Fiber location corresponding to what is specified in opt[var]:
 
       attach = FIBER, ABSCISSA, REFERENCE
 
 with
 
       FIBER = microtubule1, fiber1, fiber2, etc.
       ABSCISSA = a distance
       REFERENCE = [ plus_end, minus_end, center ]
 
 */
FiberSite FiberSet::someSite(std::string const& key, Glossary& opt) const
{
    std::string str;
    if ( opt.set(str, key) )
    {
        if ( opt.num_values(key) == 1 )
        {
            // just 'fiber' will designate all fibers:
            if ( str == title() )
                return randomSite();
            
            // check property name, designating all fibers with that name:
            Property const* p = simul_.findProperty(title(), str);
            if ( p )
                return randomSite(static_cast<FiberProp const*>(p));
        }
        // check if some individual fiber was requested:
        Fiber* fib = Fiber::toFiber(pickObject(title(), str));
        
        if ( fib )
        {
            // variables defining an abscissa:
            int mod = 7;
            real abs = 0;
            FiberEnd ref = ORIGIN;
            if ( opt.set(abs, key, 1) )
                mod = 0;
            opt.set(ref, key, 2, {{"plus_end", PLUS_END}, {"minus_end", MINUS_END}, {"center", CENTER}});
            opt.set(mod, key, 3, {{"off", 0}, {"uniform", 1}, {"exponential", 2}});
            
            return FiberSite(fib, fib->someAbscissa(abs, ref, mod, 1.0));
        }
    }
    throw InvalidParameter("could not find fiber `"+str+"'");
    return FiberSite();
}

/**
 Set a list of Locations on fibers, on sections that were recently assembled at
 the plus end. This relies on Fiber::freshAssemblyP() returning the length of
 polymer made in the last time step.
 The number of locations returned will be proportional to the total length of
 new polymer recently made, and thus proportional to simul:timestep.
 
 Because the list of fiber is regularly shuffled, the sites will consider fibers
 in a random order. However, the sites on one fiber will be listed in the order
 of increasing abscissa.

 This is for the plus end
 */
void FiberSet::newFiberSitesP(FiberSiteList& res, const real gap) const
{
    assert_true( gap > 0 );
    
    res.clear();
    Fiber * fib = first();
    real abs = gap * RNG.exponential();
    while ( fib )
    {
        real len = fib->freshAssemblyP();
        while ( abs < len )
        {
            res.emplace(fib, fib->abscissaP()-abs);
            abs += gap * RNG.exponential();
        }
        abs -= len;
        fib = fib->next();
    }
}


/**
 Set a list of Locations on fibers, on sections that were recently assembled at
 the minus end. This relies on Fiber::freshAssemblyM() returning the length of
 polymer made in the last time step.
 The number of locations returned will be proportional to the total length of
 new polymer recently made, and thus proportional to simul:timestep.
 
 Because the list of fiber is regularly shuffled, the sites will consider fibers
 in a random order. However, the sites on one fiber will be listed in the order
 of increasing abscissa.

 This is for the minus end
 */
void FiberSet::newFiberSitesM(FiberSiteList& res, const real gap) const
{
    assert_true( gap > 0 );
    
    res.clear();
    Fiber * fib = first();
    real abs = gap * RNG.exponential();
    while ( fib )
    {
        real a = fib->freshAssemblyM();
        while ( abs < a )
        {
            res.emplace(fib, fib->abscissaM()+abs);
            abs += gap * RNG.exponential();
        }
        abs -= a;
        fib = fib->next();
    }
}


void FiberSet::flipFiberPolarity(FiberProp * sel)
{
    for ( Fiber* fib=first(); fib; fib=fib->next() )
    {
        if ( fib->prop == sel )
            fib->flipPolarity();
    }
}


/**
Update Hands after reading all fibers from file, and reset Lattice values
*/
void FiberSet::updateFibers() const
{
    for ( Fiber* fib=first(); fib; fib=fib->next() )
    {
        assert_false(fib->bad());
        fib->updateRange(nullptr);
        fib->resetLattice(1);
        fib->updateHands();
    }
}


/**
 ad-hoc function to cut fibers on each side of the system
 */
void FiberSet::shortenSpindle(real dL, real dR) const
{
    Fiber * L = nullptr;
    Fiber * R = nullptr;
    real xL = +INFINITY;
    real xR = -INFINITY;
    FiberEnd eL = NO_END;
    FiberEnd eR = NO_END;
    for ( Fiber* fib=first(); fib; fib=fib->next() )
    {
        real xM = fib->posEndM().XX;
        real xP = fib->posEndP().XX;
        if ( xM < xL )
        {
            eL = MINUS_END;
            xL = xM; L = fib;
        }
        if ( xP < xL )
        {
            eL = PLUS_END;
            xL = xP; L = fib;
        }
        if ( xM > xR )
        {
            eR = MINUS_END;
            xR = xM; R = fib;
        }
        if ( xP > xR )
        {
            eR = PLUS_END;
            xR = xP; R = fib;
        }
    }
    if ( L )
    {
        std::clog << "shortenSpindle L " << L->reference() << "  " << eL << "\n";
        if ( eL == PLUS_END )
            L->cutP(dL);
        else
            L->cutM(dL);
    }
    if ( R )
    {
        std::clog << "shortenSpindle R " << R->reference() << "  " << eR << "\n";
        if ( eR == PLUS_END )
            R->cutP(dR);
        else
            R->cutM(dR);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Measurments/Quantification


real FiberSet::totalLength() const
{
    real res = 0;
    
    for ( Fiber const* fib=first(); fib; fib=fib->next() )
        res += fib->length();
    
    return res;
}


real FiberSet::totalLength(FiberProp const* sel) const
{
    real res = 0;
    
    for ( Fiber const* fib=first(); fib; fib=fib->next() )
        if ( fib->prop == sel )
            res += fib->length();
    
    return res;
}


void FiberSet::infoLength(ObjectList const& objs, size_t& cnt,
                          real& avg, real& var, real& mn, real& mx, real& off)
{
    cnt = 0;
    avg = 0;
    var = 0;
    mn = INFINITY;
    mx = 0;
    off = 0;

    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            ++cnt;
            real x = fib->length();
            avg += x;
            var += x * x;
            mn = std::min(mn, x);
            mx = std::max(mx, x);
            off += fib->abscissaM();
        }
    }
    
    if ( cnt > 0 )
    {
        avg /= cnt;
        var -= square(avg)*cnt;
    }
    if ( cnt > 1 )
        var /= real(cnt-1);
}


void FiberSet::infoBirthtime(ObjectList const& objs, size_t& cnt,
                             real& avg, real& var, real& mn, real& mx)
{
    cnt = 0;
    avg = 0;
    var = 0;
    mn = INFINITY;
    mx = 0;
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            ++cnt;
            real x = fib->birthTime();
            avg += x;
            var += x * x;
            mn = std::min(mn, x);
            mx = std::max(mx, x);
        }
    }
    
    if ( cnt > 0 )
    {
        avg /= cnt;
        var -= square(avg)*cnt;
    }
    if ( cnt > 1 )
        var /= real(cnt-1);
}


void FiberSet::infoSegments(ObjectList const& objs, size_t& cnt, size_t& points,
                            real& mn, real& mx, real& dv)
{
    cnt = 0;
    points = 0;
    mn = INFINITY;
    mx = 0;
    dv = 0;
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            ++cnt;
            real n, x;
            points += fib->nbPoints();
            fib->segmentMinMax(n, x);
            mn = std::min(mn, n);
            mx = std::max(mx, x);
            n = abs_real(1.0 - n/fib->segmentation());
            x = abs_real(1.0 - x/fib->segmentation());
            dv = std::max(dv, std::max(n, x));
        }
    }
}


size_t FiberSet::nbKinks(ObjectList const& objs)
{
    size_t cnt = 0;
    
    for ( Object * i : objs )
        cnt += Fiber::toFiber(i)->nbKinks();

    return cnt;
}


/**
 Each Fiber segment is weigthed by its length.
 
 @set M = averaged minus ends
 @set G = average center of gravity
 @set P = averaged plus ends

 @return S = sum of length
 
 An average direction can be obtained from ( P - M ) / S.
 */
real FiberSet::infoPosition(ObjectList const& objs, Vector& M, Vector& G, Vector& P)
{
    real S = 0;
    G.reset();
    P.reset();
    M.reset();
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            const real w = fib->length();
            S += w;
            M += w * fib->posEndM();
            P += w * fib->posEndP();
 
            Vector G1 = 0.5 * ( fib->posEndM() + fib->posEndP() );
            for ( index_t n = 1; n < fib->nbSegments(); ++n )
                G1 += fib->posP(n);
            G += G1 * ( w / fib->nbSegments() );
        }
    }
    
    if ( S > 0 )
    {
        G /= S;
        M /= S;
        P /= S;
    }
    return S;
}


static real computeNematicEigenvectors(real res[9], real sum, real M[9])
{
    // rescale matrix, to ensure eigenvalue = 1 in perfect order
    const real beta = ( DIM >= 3 ) ? 0.5 : 1.0;
    sum = beta * DIM / sum;
    for ( int d = 0; d < 9; ++d )
        M[d] = sum * M[d];
    //std::clog << "trace = " << M[0] + M[4] + M[8] << " (should be DIM/2)\n";
    // subtract trace:
    M[0] -= beta;
    if ( DIM > 1 ) M[4] -= beta;
    if ( DIM > 2 ) M[8] -= beta;

    int nbv = 1;
    real vec[9] = { 0 };
    real val[3] = { 0 };
    real work[32];
    int iwork[16];
    int ifail[4];
    int info = 0;

    // calculate two largest eigenvalues in 3D, one in 2D:
    lapack::xsyevx('V','I','L', DIM, M, 3, 0, 0, 2, DIM, REAL_EPSILON,
                   &nbv, val, vec, 3, work, 32, iwork, ifail, &info);

#if ( DIM > 2 )
    //std::clog << "EigenVector1 " << Vector3(vec) << " (" << val[0] << ")\n";
    //std::clog << "EigenVector2 " << Vector3(vec+3) << " (" << val[1] << ")\n";
    // order the 2 vectors in decreasing eigenvalues (reverse order from LAPACK).
    res[0] = vec[3];
    res[1] = vec[4];
    res[2] = vec[5];
    res[3] = vec[0];
    res[4] = vec[1];
    res[5] = vec[2];
    // calculate third vector as vector product of first two:
    res[6] = res[1]*res[5] - res[2]*res[4];
    res[7] = res[2]*res[3] - res[0]*res[5];
    res[8] = res[0]*res[4] - res[1]*res[3];
#else
    //std::clog << "EigenVector1 " << Vector2(vec) << " (" << val[0] << ")\n";
    res[0] = vec[0];
    res[1] = vec[1];
    res[2] = 0;
    // second vector is orthogonal:
    res[3] =-vec[1];
    res[4] = vec[0];
    res[5] = 0;
    // third vector set in Z-direction
    res[6] = 0;
    res[7] = 0;
    res[8] = 1;
#endif

    // return highest eigenvalue, which is the scalar order parameter
    return val[nbv-1];
}


/**
 Each Fiber segment is weigthed by its length.
 The Nematic direction is an eigenvector of the second rank tensor order parameter.
 
 @sets `res`, a 9-elements matrix containing the two principal eigenvectors
 @return the scalar nematic order parameter S, the eigenvalue of direction 1
 
 if DIM == 2:
     direction 1 is { res[0], res[1] }
 if DIM == 3:
     direction 1 is { res[0], res[1], res[2] }
     direction 2 is { res[3], res[4], res[5] }
 */
real FiberSet::infoNematic(ObjectList const& objs, real res[9])
{
    real sum = 0;
    real M[9] = { 0 };
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            const index_t cnt = fib->nbPoints();
            real XX = 0, XY = 0, XZ = 0, YY = 0, YZ = 0, ZZ = 0;
            Vector Q = fib->posP(0);
            for ( index_t n = 1; n < cnt; ++n )
            {
                Vector P = fib->posP(n);
                Vector d = P - Q;
                XX += d.XX * d.XX;
#if ( DIM > 1 )
                XY += d.YY * d.XX;
                YY += d.YY * d.YY;
#endif
#if ( DIM > 2 )
                XZ += d.ZZ * d.XX;
                YZ += d.ZZ * d.YY;
                ZZ += d.ZZ * d.ZZ;
#endif
                Q = P;
            }
            // we should normalize by 1/L^2, but to weight by length, we use 1/L
            real w = fib->segmentation();
            sum += ( cnt - 1 ) * w; // length of fiber
            w = 1.0 / w;
            // update lower triangle of 3x3 second-rank traceless tensor:
            M[0] += w * XX;
            M[1] += w * XY;
            M[2] += w * XZ;
            M[4] += w * YY;
            M[5] += w * YZ;
            M[8] += w * ZZ;
        }
    }
    
    if ( sum > 0 )
        return computeNematicEigenvectors(res, sum, M);
    return 0;
}


/**
This computes a nematic order parameter, and principal directions, like infoNematic,
 but using a vector that is orthogonal and tangent to the Space's edge.
 This is meaningfull only in 3D
 */
real FiberSet::infoOrthoNematic(ObjectList const& objs, real res[9], Space const* spc)
{
    real sum = 0;
    real M[9] = { 0 };
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            const index_t cnt = fib->nbSegments();
            real XX = 0, XY = 0, XZ = 0, YY = 0, YZ = 0, ZZ = 0;
            for ( index_t n = 0; n < cnt; ++n )
            {
                Vector pos = fib->posPoint(n);
                Vector dir = fib->dirSegment(n);
                Vector nor = spc->normalToEdge(pos);
#if ( DIM > 2 )
                Vector p = cross(nor, dir);
                XX += p.XX * p.XX;
                XY += p.YY * p.XX;
                YY += p.YY * p.YY;
                XZ += p.ZZ * p.XX;
                YZ += p.ZZ * p.YY;
                ZZ += p.ZZ * p.ZZ;
#endif
            }
            // we weight each segment by its length
            real w = fib->segmentation();
            sum += cnt * w; // length of fiber
            // update lower triangle of 3x3 second-rank traceless tensor:
            M[0] += w * XX;
            M[1] += w * XY;
            M[2] += w * XZ;
            M[4] += w * YY;
            M[5] += w * YZ;
            M[8] += w * ZZ;
        }
    }
    
    if ( sum > 0 )
        return computeNematicEigenvectors(res, sum, M);
    return 0;
}


/**
 Calculates the principal component directions of the cloud of vertices.
 Each fiber is weighted by its length.
 
 @set avg[] = average center of gravity
 @set mom[], a 9-elements matrix containing the moments in its lower part:
 - mom[0] = sum( ( X - mean(X) ) * ( X - mean(X) ) ) / S
 - mom[1] = sum( ( X - mean(X) ) * ( Y - mean(Y) ) ) / S
 - mom[2] = sum( ( X - mean(X) ) * ( Z - mean(Z) ) ) / S
 - mom[4] = sum( ( Y - mean(Y) ) * ( Y - mean(Y) ) ) / S
 - mom[5] = sum( ( Y - mean(Y) ) * ( Z - mean(Z) ) ) / S
 - mom[8] = sum( ( Z - mean(Z) ) * ( Z - mean(Z) ) ) / S
 .
 
 @set res[], a 9-elements matrix containing the first two principal component vectors
 
 if DIM == 2:
   Component 1 is { res[0], res[1] }
 if DIM == 3:
   Component 1 is { res[0], res[1], res[2] }
   Component 2 is { res[3], res[4], res[5] }
 
 @return 0 if everything proceeded without error
 */
int FiberSet::infoComponents(ObjectList const& objs,
                             real& sum, real avg[3], real mom[9], real res[9])
{
    sum = 0;
    avg[0] = 0.0;
    avg[1] = 0.0;
    avg[2] = 0.0;
    real M[9] = { 0 };
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            const real w = fib->length() / fib->nbPoints();
            for ( index_t n = 0; n < fib->nbPoints(); ++n )
            {
                Vector p = fib->posP(n);
                avg[0] += w * p.XX;
                M[0]   += w * p.XX * p.XX;
#if ( DIM > 1 )
                avg[1] += w * p.YY;
                M[1]   += w * p.YY * p.XX;
                M[4]   += w * p.YY * p.YY;
#endif
#if ( DIM > 2 )
                avg[2] += w * p.ZZ;
                M[2]   += w * p.ZZ * p.XX;
                M[5]   += w * p.ZZ * p.YY;
                M[8]   += w * p.ZZ * p.ZZ;
#endif
            }
            sum += w * fib->nbPoints();
        }
    }
    
    if ( sum == 0 )
        return 1;
    
    /**
     Remove the mean:
       (x-a)*(x-a) = x*x - 2x*a + a*a, hence <(x-a)^2> = <x^2> - a^2
       (x-a)*(y-b) = x*y - x*b - y*a + a*b
     */
    
    avg[0] /= sum;
    M[0] = M[0]/sum - avg[0] * avg[0];
#if ( DIM > 1 )
    avg[1] /= sum;
    M[1] = M[1]/sum - avg[1] * avg[0];
    M[4] = M[4]/sum - avg[1] * avg[1];
#endif
#if ( DIM > 2 )
    avg[2] /= sum;
    M[2] = M[2]/sum - avg[2] * avg[0];
    M[5] = M[5]/sum - avg[2] * avg[1];
    M[8] = M[8]/sum - avg[2] * avg[2];
#endif
    
    // copy moments:
    for ( int i = 0; i < 9; ++i )
        mom[i] = M[i];
    
    int nbv = 1;
    real vec[9] = { 0 };
    real val[3] = { 0 };
    real work[32];
    int iwork[16];
    int ifail[4];
    int info = 0;

    // calculate two largest eigenvalues in 3D, one in 2D:
    lapack::xsyevx('V','I','L', DIM, M, 3, 0, 0, 2, DIM, REAL_EPSILON,
                   &nbv, val, vec, 3, work, 32, iwork, ifail, &info);

#if ( DIM > 2 )
    real u = sign_real(vec[3]);
    real v = sign_real(vec[0]);
    // order the 2 vectors in decreasing eigenvalues.
    res[0] = u * vec[3];
    res[1] = u * vec[4];
    res[2] = u * vec[5];
    res[3] = v * vec[0];
    res[4] = v * vec[1];
    res[5] = v * vec[2];
    // calculate third vector as vector product for first two:
    res[6] = res[1]*res[5] - res[2]*res[4];
    res[7] = res[2]*res[3] - res[0]*res[5];
    res[8] = res[0]*res[4] - res[1]*res[3];
#else
    real u = sign_real(vec[0]);
    res[0] =  u * vec[0];
    res[1] =  u * vec[1];
    res[2] =  0;
    // second vector is orthogonal:
    res[3] = -u * vec[1];
    res[4] =  u * vec[0];
    res[5] =  0;
    res[6] =  0;
    res[7] =  0;
    res[8] =  1;
#endif
    if ( nbv != DIM-1 )
        return 2;
    return info;
}


/**
 Counts the number of fiber intersecting the plane defined by <em> n.pos + a = 0 </em>
 in two categories, depending on the direction with which they cross the plane:
 - `np` = number of parallel segments ( the scalar product dir.n is strictly positive )
 - `na` = number of anti-parallel segments ( dir.n < 0 )
 .
 */
void FiberSet::infoPlane(int& np, int& na, Vector const& n, real a) const
{
    np = 0;
    na = 0;
    for ( Fiber const* fib=first(); fib; fib=fib->next() )
    {
        for ( index_t s = 0; s < fib->nbSegments(); ++s )
        {
            real abs = fib->planarIntersect(s, n, a);
            if (( 0 <= abs ) & ( abs < 1 ))
            {
                real sec = dot(n, fib->dirSegment(s));
                np += ( sec > 0 );
                na += ( sec < 0 );
            }
        }
    }
}


/**
 Calculate two indices characterizing the organization of the fibers along the axis `n`.
 - `ixa` = average { ( o - i ) }
 - `ixp` = average { ( r - l ) }
 .
 where:
 - `o` = number of fiber pointing outward (away from the mid-plane),
 - `i` = number of fiber pointing inward (toward the mid-plane),
 - `r` = number of fiber pointing right (ie. in the direction of `n`),
 - `l` = number of fiber pointing left.
 .
 
 The indices are averaged over planar sections taken every `dm` units of space,
 and the values for each planar section are weighted by the number of fibers.
 The central symmetry plane is defined by `n.x+a=0`, and the edges correspond to `n.x+a=+/-m`.
 
 The results characterize broadly the type of fiber organization:
 - `ixa =  1, ixp = 0`:   aster,
 - `ixa = -1, ixp = 0`:   anti-aster,
 - `ixa =  0, ixp = 1`:   parallel overlap,
 - `ixa =  0, ixp = 0`:   anti-parallel overlap (50/50).
 .
 */
void FiberSet::infoSpindle(real& ixa, real& ixp, Vector const& dir, real a, real m, real dm) const
{
    ixa = 0;
    ixp = 0;
    int no, ni, nio;
    int sum = 0;
    for ( real p = dm/2 ; p < m ; p += dm )
    {
        // left side
        infoPlane(ni, no, dir, a+p);
        nio = ni + no;
        if ( nio )
        {
            ixa += ( no - ni );
            ixp += ( ni - no );
            sum += nio;
        }
    
        // right side: arguments are swapped!
        infoPlane(no, ni, dir, a-p);
        nio = ni + no;
        if ( nio )
        {
            ixa += ( no - ni );
            ixp += ( no - ni );
            sum += nio;
        }
    }
    if ( sum )
    {
        ixa /= sum;
        ixp /= sum;
    }
}


/**
 Sum elastic bending energy of all the fibers `fib` for which func(fib, arg) == true
 */
void FiberSet::infoBendingEnergy(ObjectList const& objs, size_t& cnt,
                                 real& avg, real& var)
{
    cnt = 0;
    avg = 0;
    var = 0;
    
    for ( Object * i : objs )
    {
        Fiber * fib = Fiber::toFiber(i);
        if ( fib )
        {
            ++cnt;
            real x = fib->bendingEnergy();
            avg += x;
            var += x * x;
        }
    }
    
    if ( cnt > 0 )
    {
        avg /= cnt;
        var -= square(avg)*cnt;
    }
    if ( cnt > 1 )
        var /= real(cnt-1);
}

/**
 Sum tensions of all fiber segments intersecting the plane defined by
      <em> n.pos + a = 0 </em>
 (hence vector `n` defines the direction orthogonal to the plane)
 
 The intersecting segments are determined by testing all Fibers.
 The tension dipole along a segment is obtained from the Lagrange multiplier 
 associated with the length of this segment. It is positive if the segment is stretched.
 The magnitude of the dipole is multiplied by the cosine of the angle measured between 
 the segment and the plane normal, yielding axial components that can be summed.
 
 @return cnt = number of segments intersecting the plane
 @return sum = sum of tension in these segments
 @return inf = minimum of tension
 @return sup = maximum of tension
*/
void FiberSet::infoTension(size_t& cnt, real& sum, real& inf, real& sup, Vector const& n, real a) const
{
    cnt = 0;
    sum = 0;
    inf = INFINITY;
    sup = -INFINITY;

    Vector dir = normalize(n);
    for ( Fiber const* fib=first(); fib; fib=fib->next() )
    {
        for ( index_t s = 0; s < fib->nbSegments(); ++s )
        {
            real abs = fib->planarIntersect(s, n, a);
            if (( 0 <= abs ) & ( abs < 1 ))
            {
                real t = fib->tension(s);
                real h = t * abs_real(dot(dir, fib->dirSegment(s)));
                sum += h;
                inf = std::min(inf, h);
                sup = std::max(sup, h);
                ++cnt;
            }
        }
    }
}


/**
 Sum tension of all the segments
 
 @return cnt = total number of segments
 @return ten = sum of tension
 @return max = maximum of tension
 */
void FiberSet::infoTension(size_t& cnt, real& sum, real& inf, real& sup) const
{
    cnt = 0;
    sum = 0;
    inf = INFINITY;
    sup = -INFINITY;

    for ( Fiber const* fib=first(); fib; fib=fib->next() )
    {
        for ( index_t s = 0; s < fib->nbSegments(); ++s )
        {
            real h = fib->tension(s);
            sum += h;
            inf = std::min(inf, h);
            sup = std::max(sup, h);
            ++cnt;
        }
    }
}


void FiberSet::infoRadius(size_t& cnt, real& rad) const
{
    real r = 0;
    cnt = 0;
    
    for ( Fiber const* f=first(); f; f=f->next() )
    {
        for ( index_t p = 0; p < f->nbPoints() ; ++p )
        {
            r += f->posP(p).norm();
            ++cnt;
        }
    }
    if ( cnt )
        rad = r / (real)cnt;
}


void FiberSet::infoRadius(size_t& cnt, real& rad, FiberEnd end) const
{
    real r = 0;
    cnt = 0;
    
    for ( Fiber const* f=first(); f; f=f->next() )
    {
        r += f->posEnd(end).norm();
        ++cnt;
    }
    if ( cnt )
        rad = r / (real)cnt;
}

