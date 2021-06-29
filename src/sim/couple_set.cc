// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul_prop.h"
#include "couple_set.h"
#include "couple_prop.h"
#include "fork_prop.h"
#include "messages.h"
#include "property_list.h"
#include "crosslink_prop.h"
#include "shackle_prop.h"
#include "bridge_prop.h"
#include "duo_prop.h"
#include "glossary.h"
#include "simul.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

void CoupleSet::prepare(PropertyList const& properties)
{
    uniEnabled = uniPrepare(properties);
}


void CoupleSet::step()
{
    /*
     ATTENTION: We ensure here that step() is called exactly once for each object.
     The Couples are stored in multiple lists, and are automatically transferred
     from one list to another one if their Hands bind or unbind. The code relies
     on the fact that a Couple will be moved to the start of the list to which it
     is transferred, using 'push_front'. By proceeding always from the node, which
     was first before any transfer could occur, we process each Couple only once.
     Moreover, we get the 'next' in the list always before calling 'step()', because
     'step()' may transfer the node to another list, changing the value of 'next()'
     */
    
    /*
    Cytosim::log("CoupleSet::step : FF %5i AF %5i FA %5i AA %5i\n",
                 ffList.size(), afList.size(), faList.size(), aaList.size());
    */
    
    Couple *const ffHead = firstFF();
    Couple *const afHead = firstAF();
    Couple *const faHead = firstFA();
    
    bool const faOdd = faList.size() & 1;
    bool const afOdd = afList.size() & 1;
    bool const ffOdd = ffList.size() & 1;

    Couple * obj, * nxt;
    
    obj = firstAA();
    // this loop is unrolled, processing objects 2 by 2:
    if ( aaList.size() & 1 )
    {
        nxt = obj->next();
        obj->stepAA();
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->stepAA();
        obj = nxt->next();
        nxt->stepAA();
    }
    
    obj = faHead;
    // this loop is unrolled, processing objects 2 by 2:
    if ( faOdd )
    {
        nxt = obj->next();
        obj->stepFA();
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->stepFA();
        obj = nxt->next();
        nxt->stepFA();
    }

    obj = afHead;
    // this loop is unrolled, processing objects 2 by 2:
    if ( afOdd )
    {
        nxt = obj->next();
        obj->stepAF();
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        obj->stepAF();
        obj = nxt->next();
        nxt->stepAF();
    }
    
    // use alternative attachment strategy:
    if ( uniEnabled )
    {
        obj = uniCollect(ffHead);
        uniAttach(simul_.fibers);
        while ( obj )
        {
            nxt = obj->next();
            obj->stepFF();
            obj = nxt;
        }
    }
    else
    {
        //std::clog << "CoupleSet::step : FF " << ffList.size() << " head " << ffHead << '\n';
        // this loop is unrolled, processing objects 2 by 2:
        obj = ffHead;
        if ( ffOdd )
        {
            nxt = obj->next();
            obj->stepFF();
            obj = nxt;
        }
        while ( obj )
        {
            nxt = obj->next();
            obj->stepFF();
            obj = nxt->next();
            nxt->stepFF();
        }
    }

    //printf("  : %lu couples [ %u %u ]\n", size(), inventory.first_identity(), inventory.last_identity());
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 @defgroup CoupleGroup Couple and related
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Couple contains two Hand, and can thus crosslink two Fibers.

 The plain Couple may crosslink two Fiber irrespective of their configuration.
 Derived classes implement specificity, angular stiffness, etc.
 
 List of classes accessible by specifying `couple:activity`.

 `activity`          | Classes                 | Parameters         | Property     |
 --------------------|-------------------------|--------------------|---------------
 `diffuse` (default) | Couple CoupleLong       | @ref CouplePar     | CoupleProp
 `crosslink`         | Crosslink CrosslinkLong | @ref CrosslinkPar  | CrosslinkProp
 `bridge`            | Bridge                  | @ref BridgePar     | BridgeProp
 `duo`               | Duo  DuoLong            | @ref DuoPar        | DuoProp
 `slide`             | Shackle ShackleLong     | @ref ShacklePar    | ShackleProp
 `fork`              | Fork                    | @ref ForkPar       | ForkProp

 Example:

     set couple complex
     {
       hand1 = kinesin
       hand2 = kinesin
       stiffness = 100
       diffusion = 10
       activity = crosslink
       length = 0.02
     }

 */

Property* CoupleSet::newProperty(const std::string& cat, const std::string& nom, Glossary& opt) const
{
    if ( cat == "couple" )
    {
        std::string a;
        if ( opt.peek(a, "activity") )
        {
            if ( a == "fork" )
                return new ForkProp(nom);
            if ( a == "crosslink" )
                return new CrosslinkProp(nom);
            if ( a == "bridge" )
                return new BridgeProp(nom);
            if ( a == "duo" )
                return new DuoProp(nom);
            if ( a == "slide" )
                return new ShackleProp(nom);
            if ( a == "diffuse" )
                return new CoupleProp(nom);
#if ( 0 )
            throw InvalidParameter("unknown couple:activity `"+a+"'");
#else
        // try to proceed anyhow:
        std::cerr << "WARNING: unknown couple:activity `"+a+"'\n";
#endif
        }
        return new CoupleProp(nom);
    }
    return nullptr;
}


Object * CoupleSet::newObject(const ObjectTag tag, size_t num)
{
    if ( tag == Couple::TAG )
    {
        CoupleProp * p = simul_.findProperty<CoupleProp>("couple", num);
        return p->newCouple();
    }
    std::cerr << "Warning: unknown Couple tag `"+std::string(1,tag)+"' requested\n";
    return nullptr;
}


/**
 @addtogroup CoupleGroup

 You can attach the hands of a Couple:
 
     new complex
     {
        attach1 = FIBER, REAL, REFERENCE
        attach2 = FIBER, REAL, REFERENCE
     }
 
 where:
 - FIBER designates the fiber:
     - `fiber1` of `fiber2` correspond to fibers directly
     - `first` or `last` to the oldest and youngest fiber
     - `last-1` the penultimate, etc.
     .
 - REAL is the abscissa of the attachment point.
   If the abscissa is not specified, and random position along
   along the fiber will be selected.
 - REFERENCE can be `minus_end`, `center` or `plus_end` (default = `origin`).
   This defines from which position the abscissa is measured.
 .
 
 */
ObjectList CoupleSet::newObjects(const std::string& name, Glossary& opt)
{
    CoupleProp * p = simul_.findProperty<CoupleProp>("couple", name);
    Couple * obj = p->newCouple(&opt);
    
    ObjectList res;
    res.push_back(obj);
        
    // Allow user to attach hand1:
    if ( opt.has_key("attach1") )
        obj->attach1(simul_.fibers.someSite("attach1", opt));

    // Allow user to attach hand2:
    if ( opt.has_key("attach2") )
        obj->attach2(simul_.fibers.someSite("attach2", opt));

    /* It would be possible to create Couple with custom hand type, and the
    syntax below to attach the Hands could be better used for this */
    
    // Allow user to attach hand1:
    if ( opt.has_key("site1") )
        obj->attach1(simul_.fibers.someSite("site1", opt));
    
    // Allow user to attach hand2:
    if ( opt.has_key("site2") )
        obj->attach2(simul_.fibers.someSite("site2", opt));

    return res;
}

//------------------------------------------------------------------------------
#pragma mark -

void CoupleSet::relinkA1(Couple * obj)
{
    assert_true( obj->attached1() );

    if ( obj->attached2() )
    {
        faList.pop(obj);
        aaList.push_front(obj);
    }
    else
    {
        ffList.pop(obj);
        afList.push_front(obj);
    }
}


void CoupleSet::relinkD1(Couple * obj)
{
    assert_true( obj->attached1() );
    
    if ( obj->attached2() )
    {
        aaList.pop(obj);
        faList.push_front(obj);
    }
    else
    {
        afList.pop(obj);
        ffList.push_front(obj);
    }
}


void CoupleSet::relinkA2(Couple * obj)
{
    assert_true( obj->attached2() );

    if ( obj->attached1() )
    {
        afList.pop(obj);
        aaList.push_front(obj);
    }
    else
    {
        ffList.pop(obj);
        faList.push_front(obj);
    }
}


void CoupleSet::relinkD2(Couple * obj)
{
    assert_true( obj->attached2() );

    if ( obj->attached1() )
    {
        aaList.pop(obj);
        afList.push_front(obj);
    }
    else
    {
        faList.pop(obj);
        ffList.push_front(obj);
    }
}


void CoupleSet::link(Object * obj)
{
    assert_true( obj->tag() == Couple::TAG );
    
    Couple * c = static_cast<Couple*>(obj);
    sublist(c->attached1(), c->attached2()).push_back(obj);
    
    //std::clog << "CoupleSet has " << ffList.size() << "  " << afList.size() << "  " << faList.size() << "  " << aaList.size() << '\n';
}


/**
 This will also detach both Hands
 */
void CoupleSet::unlink(Object * obj)
{
    Couple * c = static_cast<Couple*>(obj);

    sublist(c->attached1(), c->attached2()).pop(obj);
    if ( c->attached1() ) c->hand1()->detach();
    if ( c->attached2() ) c->hand2()->detach();
}


//------------------------------------------------------------------------------
#pragma mark -

void CoupleSet::foldPositions(Modulo const* m) const
{
    Couple * cx;
    for ( cx=firstAA(); cx; cx=cx->next() )  cx->foldPosition(m);
    for ( cx=firstFA(); cx; cx=cx->next() )  cx->foldPosition(m);
    for ( cx=firstAF(); cx; cx=cx->next() )  cx->foldPosition(m);
    for ( cx=firstFF(); cx; cx=cx->next() )  cx->foldPosition(m);
}


void CoupleSet::shuffle()
{
    ffList.shuffle();
    afList.shuffle();
    faList.shuffle();
    aaList.shuffle();
}


void CoupleSet::erase()
{
    relax();
    ObjectSet::erasePool(aaList);
    ObjectSet::erasePool(faList);
    ObjectSet::erasePool(afList);
    ObjectSet::erasePool(ffList);
    inventory_.clear();
}


void CoupleSet::freeze()
{
    relax();
    assert_true(ice_.empty());
    ice_.append(aaList);
    ice_.append(faList);
    ice_.append(afList);
    ice_.append(ffList);
}


void CoupleSet::pruneDetach()
{
    /* After reading from file, the Hands should not
     update any Fiber, Single or Couple as they will be deleted */
    Object * i = ice_.pop_front();
    while ( i )
    {
        Couple* o = static_cast<Couple*>(i);
        i = ice_.pop_front();
        if ( o->attached1() ) o->hand1()->detachHand();
        if ( o->attached2() ) o->hand2()->detachHand();
        link(o);
    }
}


/* After reading from file, the Hands should not
 update any Fiber, Single or Couple as they will be deleted */
void CoupleSet::pruneDelete()
{
    Object * i = ice_.pop_front();
    while ( i )
    {
        Couple* o = static_cast<Couple*>(i);
        i = ice_.pop_front();
        if ( o->attached1() ) o->hand1()->detachHand();
        if ( o->attached2() ) o->hand2()->detachHand();
        inventory_.unassign(o);
        o->objset(nullptr);
        delete(o);
    }
}


void CoupleSet::prune()
{
    if ( prune_mode )
        pruneDetach();
    else
        pruneDelete();
}


//------------------------------------------------------------------------------
#pragma mark -


void CoupleSet::writeAA(Outputter& out) const
{
    out.writeLine("\n#section couple AA");
    writeObjects(out, aaList);
}

void CoupleSet::writeAF(Outputter& out) const
{
    out.writeLine("\n#section couple AF");
    writeObjects(out, afList);
}

void CoupleSet::writeFA(Outputter& out) const
{
    out.writeLine("\n#section couple FA");
    writeObjects(out, faList);
}

void CoupleSet::writeFF(Outputter& out) const
{
    out.writeLine("\n#section couple FF 0");
    writeObjects(out, ffList);
}

void CoupleSet::write(Outputter& out) const
{
    if ( sizeAA() > 0 )
        writeAA(out);
    if ( sizeAF() > 0 )
        writeAF(out);
    if ( sizeFA() > 0 )
        writeFA(out);
    if ( sizeFF() > 0 )
    {
        int skip = simul_.prop->skip_free_couple || prune_mode;
        if ( skip & skip_now )
            out.writeLine("\n#section couple FF 1");
        else
        {
            writeFF(out);
            skip_now = skip;
        }
    }
}


void CoupleSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << '\n' << title();
        PropertyList plist = simul_.properties.find_all(title());
        for ( Property const* i : plist )
        {
            CoupleProp const* p = static_cast<CoupleProp const*>(i);
            size_t cnt = count(match_property, p);
            os << '\n' << std::setw(10) << cnt << " " << p->name();
            os << " ( " << p->hand1 << " | " << p->hand2 << " )";
        }
        if ( plist.size() > 1 )
            os << '\n' << std::setw(10) << size() << " total";
    }
}


ObjectList CoupleSet::collect() const
{
    ObjectList res = ObjectSet::collect(ffList);
    res.append( ObjectSet::collect(afList) );
    res.append( ObjectSet::collect(faList) );
    res.append( ObjectSet::collect(aaList) );
    return res;
}


ObjectList CoupleSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    ObjectList res = ObjectSet::collect(ffList, func, arg);
    res.append( ObjectSet::collect(afList, func, arg) );
    res.append( ObjectSet::collect(faList, func, arg) );
    res.append( ObjectSet::collect(aaList, func, arg) );
    return res;
}


size_t CoupleSet::count(bool (*func)(Object const*, void const*), void const* arg) const
{
    size_t ff = ObjectSet::count(ffList, func, arg);
    size_t af = ObjectSet::count(afList, func, arg);
    size_t fa = ObjectSet::count(faList, func, arg);
    size_t aa = ObjectSet::count(aaList, func, arg);
    return ff + af + fa + aa;
}

/**
 Sum tensions of all Couples stretching accross the plane defined by `n.pos + a = 0`

 The tension is normally positive for stretching
 */
void CoupleSet::infoTension(size_t& cnt, real& sum, real& inf, real& sup, Vector const& n, real a) const
{
    cnt = 0;
    sum = 0;
    inf = INFINITY;
    sup = -INFINITY;

    Vector dir = normalize(n);
    for ( Couple * c = firstAA(); c; c = c->next() )
    {
        Vector h1 = c->posHand1();
        Vector h2 = c->posHand2();
        if ( modulo )
        {
            // calculate image that is closest to plane:
            Vector cen = n * ( -a / n.normSqr() );
            modulo->fold(h1, cen);
            modulo->fold(h2, h1);
        }
        real x = dot(n, h1) + a;
        real y = dot(n, h2) + a;
        if ( x * y < 0 )
        {
            real h = dot(dir, c->force()) * sign_real(y-x);
            inf = std::min(inf, h);
            sup = std::max(sup, h);
            sum += h;
            ++cnt;
        }
    }
}


int CoupleSet::bad() const
{
    int code = 0;
    Couple * obj;
    code |= ffList.bad();
    for ( obj=firstFF(); obj ; obj = obj->next() )
    {
        if ( obj->attached1() || obj->attached2() )
            code |= 8;
    }
    
    code |= afList.bad();
    for ( obj=firstAF(); obj ; obj = obj->next() )
    {
        if ( !obj->attached1() || obj->attached2() )
            code |= 16;
    }
    
    code |= faList.bad();
    for ( obj=firstFA(); obj ; obj = obj->next() )
    {
        if ( obj->attached1() || !obj->attached2() )
            code |= 32;
    }
    
    code |= aaList.bad();
    for ( obj=firstAA(); obj ; obj = obj->next() )
    {
        if ( !obj->attached1() || !obj->attached2() )
            code |= 64;
    }
    return code;
}


//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion


void CoupleSet::uniRefill(CoupleList& can, size_t cnt, CoupleProp const* p)
{
    for ( size_t i = can.size(); i < cnt; ++i )
    {
        Couple* c = p->newCouple();
        inventory_.assign(c);
        can.push_back(c);
    }
}


/**
 Attach Hand1 of exactly one Couple from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void CoupleSet::uniAttach1(Array<FiberSite>& loc, CoupleList& can)
{
    // crop list to match available number of candidates:
    if ( can.size() < loc.size() )
    {
        loc.shuffle();
        loc.truncate(can.size());
    }

    for ( FiberSite & i : loc )
    {
        Couple * c = can.back();
        Hand * h = c->hand1();
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            can.pop_back();
            h->attach(i);
            link(c);
        }
    }
}


/**
 Attach Hand2 of exactly one Couple from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void CoupleSet::uniAttach2(Array<FiberSite>& loc, CoupleList& can)
{
    // crop list to match available number of candidates:
    if ( can.size() < loc.size() )
    {
        loc.shuffle();
        loc.truncate(can.size());
    }

    for ( FiberSite & i : loc )
    {
        Couple * c = can.back();
        Hand * h = c->hand2();
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            can.pop_back();
            h->attach(i);
            link(c);
        }
    }
}


/**
 Distribute up to `nb` Couples from `can` by attaching them to locations specified
 by `loc1` and `loc2`. These positions correspond fibers crossings, calculated by
 FiberSet::allIntersections().
 */
void CoupleSet::uniAttach12(Array<FiberSite>& loc1, Array<FiberSite>& loc2,
                            CoupleList& can, size_t nb)
{
    assert_true(loc1.size() == loc2.size());
    
    const size_t nbc = loc1.size();
    const size_t sup = std::min(nb, can.size());

    for ( size_t n = 0; n < sup; ++n )
    {
        Couple * c = can.back();
        can.pop_back();
        // pick a random point to attach:
        size_t p = RNG.pint32(nbc);
        c->attach1(loc1[p]);
        c->attach2(loc2[p]);
        link(c);
    }
}



/**
 Implements a Monte-Carlo approach for attachments of free Couple, assumming that
 diffusion is sufficiently fast to maintain a uniform spatial distribution,
 and that the distribution of fibers is more-or-less uniform such that the
 attachments are distributed randomly along the fibers.
 
 Diffusing (free) Couple are removed from the standard list, and thus the
 random walk that is used for simulating diffusion will be skipped,
 as well as the detection of neighboring fibers done for attachments.
 The attachment of already attached Couple is unchanged.
 
 Algorithm:
 - Remove diffusing Single from the simulation, transfering them to a 'reserve'.
 - Estimate the distance between binding sites occuring in one time-step, from:
    - the total length of fibers,
    - the volume of the Space,
    - the binding parameters of the relevant Hand.
    .
 - Attach Singles from the reserve, at random positions along the Fibers
 .
 
 Note: there is a similar feature for Single
 */
void CoupleSet::uniAttach(FiberSet const& fibers)
{
    // preallocate array:
    Array<FiberSite> loc(1024, 1024);
    
#if ( 0 )
    // this performs a basic verification of fibers.uniFiberSites()
    size_t rep = 1<<10;
    double avg = 0, var = 0;
    for ( size_t i = 0; i < rep; ++i )
    {
        fibers.uniFiberSites(loc, 1.0);
        real s = loc.size();
        avg += s;
        var += s*s;
    }
    avg /= rep;
    var = var/(rep-1) - avg * avg;
    printf("UNI-FIBER-SITES(1)  avg = %9.2f   var = %9.2f\n", avg, var);
#endif
    
    // uniform attachment for reserved couples:
    for ( CoupleReserve & reserve : uniReserves )
    {
        CoupleProp const * p = reserve.first;
        if ( !p )
            continue;

        CoupleList& can = reserve.second;
        
        // assuming (or not) a fixed number of diffusing molecules
        bool fixed = ( p->fast_diffusion_nb > 0 );
        size_t cnt = ( fixed ? p->fast_diffusion_nb : can.size());

        if ( cnt > 0 )
        {
            const real alpha = 2 * p->spaceVolume() / cnt;
            
            if ( p->fast_diffusion & 2 )
            {
                real dis = alpha / p->hand1_prop->bindingSectionRate();
                fibers.newFiberSitesP(loc, dis);
            }
            else
            {
                real dis = alpha / p->hand1_prop->bindingSectionProb();
                fibers.uniFiberSites(loc, dis);
            }
            
            if ( fixed ) // create enough candidates for all sites
                uniRefill(can, loc.size(), p);
            uniAttach1(loc, can);
            
            if ( p->trans_activated )
                continue;
            
            // if ( couple:trans_activated == true ), Hand2 cannot bind
            if ( p->fast_diffusion & 2 )
            {
                real dis = alpha / p->hand2_prop->bindingSectionRate();
                fibers.newFiberSitesP(loc, dis);
            }
            else
            {
                real dis = alpha / p->hand2_prop->bindingSectionProb();
                fibers.uniFiberSites(loc, dis);
            }
            
            if ( fixed ) // create enough candidates for all sites
                uniRefill(can, loc.size(), p);
            uniAttach2(loc, can);
        }
    }
}


/**
 
 Return true if at least one couple:fast_diffusion is true,
 and in this case allocate uniReserves.
 
 The Volume of the Space is assumed to remain constant until the next uniPrepare() 
 */
bool CoupleSet::uniPrepare(PropertyList const& properties)
{
    bool res = false;
    PropertyID last = 0;
    
    PropertyList allprop = properties.find_all("couple");

    for ( Property const* i : allprop )
    {
        CoupleProp const * p = static_cast<CoupleProp const*>(i);
        last = std::max(last, p->number());
        res |= p->fast_diffusion;
    }
    
    if ( res )
    {
        uniReserves.resize(last+1);
        for ( CoupleReserve & reserve : uniReserves )
            reserve.first = nullptr;
        for ( Property const* i : allprop )
        {
            CoupleProp const* p = static_cast<CoupleProp const*>(i);
            if ( p->fast_diffusion )
                uniReserves[p->number()].first = p;
        }
    }
    
    return res;
}


/**
 Transfer free Couple with `fast_diffusion` to the reserves, starting from `obj`.
 Return first Couple that was not transferred.
*/
Couple* CoupleSet::uniCollect(Couple * obj)
{
    Couple * res = nullptr;
    Couple * nxt;
    while ( obj )
    {
        nxt = obj->next();
        CoupleProp const* p = static_cast<CoupleProp const*>(obj->property());
        if ( p->fast_diffusion )
        {
            ffList.pop(obj);
            assert_true(p->number() < uniReserves.size());
            uniReserves[p->number()].second.push_back(obj);
        }
        else if ( !res )
            res = obj;
        obj = nxt;
    }
    return res;
}


/**
 Release all Couples from the reserves
 */
void CoupleSet::uniRelax()
{
    for ( CoupleReserve & reserve : uniReserves )
    {
        for ( Couple * c : reserve.second )
        {
            assert_true(!c->attached1() && !c->attached2());
            c->randomizePosition();
            link(c);
        }
        reserve.second.clear();
    }
}


//------------------------------------------------------------------------------
# pragma mark - Equilibrate


void CoupleSet::equilibrateSym(FiberSet const& fibers, CoupleList& can, CoupleProp const* cop)
{
    Array<FiberSite> loc1(1024, 1024), loc2(1024, 1024);
    
    if ( cop->hand1_prop != cop->hand2_prop )
        throw InvalidParameter("Cannot equilibrate heterogeneous Couple");
    
    if ( cop->trans_activated )
        throw InvalidParameter("Cannot equilibrate trans_activated Couple");

    const real space_volume = cop->spaceVolume();
    const real total_length = fibers.totalLength();

    if ( space_volume <= 0 )
        throw InvalidParameter("Cannot equilibrate as Space:volume == 0");
    
    const real bind_rate = cop->hand1_prop->binding_rate;
    const real bind_range = cop->hand1_prop->binding_range;
    const real unbind_rate = cop->hand1_prop->unbinding_rate;
    
    // get all crosspoints:
    fibers.allIntersections(loc1, loc2, bind_range);
    const size_t num_crossings = loc1.size();
    //const real num_crossings = square(total_length) / ( M_PI * space_volume );

    const real ratio_fibs = 2 * total_length * bind_range / space_volume;
    const real ratio_cros = 4 * M_PI * num_crossings * square(bind_range) / space_volume;
    
    /*
     The different states are defined in Belmonte et al. 2017, supplementary:
     Free, Bridge, Attached in location that cannot bridge, G=Attached near crosspoint
     */
    real bind = bind_rate / unbind_rate;
    real BsG = bind / 2;
    real AsF = ( ratio_fibs - ratio_cros ) * bind;
    real GsF = ratio_cros * bind;
    
    real popF = can.size() / ( 1 + AsF + GsF + BsG * GsF );
    real popA = AsF * popF;
    real popG = GsF * popF;
    real popB = BsG * popG;
    
#if ( 0 )
    printf("Couple::equilibrate %s (sym):\n", cop->name_str());
    printf("     total %lu\n", reserve.size());
    const real num_fibers = fibers.size();
    const real fiber_length = total_length / num_fibers;
    const real nbc = num_fibers * ( num_fibers - 1 ) * square(fiber_length) / ( M_PI * space_volume );
    //const real nbc = square(total_length) / ( M_PI * space_volume );
    printf("     num_crossings predicted  %9.2f   true %9i\n", nbc, num_crossings);
    printf("     F %9.2f A %9.2f G %9.2f B %9.2f\n", popF, popA, popG, popB);
#endif
    
    if ( !can.empty() )
    {
        // distribute Couples at filament's crosspoints:
        uniAttach12(loc1, loc2, can, RNG.poisson(popB));

        real dis = 2 * total_length / ( popA + popG );

        fibers.uniFiberSites(loc1, dis);
        uniAttach1(loc1, can);
    
        fibers.uniFiberSites(loc2, dis);
        uniAttach2(loc2, can);
    }
}


/**
 This attempts to create a configuration of Couple that is close to the equilibrium
 that would be reached, after a sufficient time is given for binding and unbinding.
 This assumes that the configuration of filaments does not change, and also that
 it is random, in particular without bundles. The motion of the motor is also ignored.
 */
void CoupleSet::equilibrate(FiberSet const& fibers, CoupleList& can, CoupleProp const* cop)
{
    Array<FiberSite> loc1(1024, 1024), loc2(1024, 1024);
    if ( cop->trans_activated )
        throw InvalidParameter("Cannot equilibrate trans_activated Couple");
    
    const real space_volume = cop->spaceVolume();
    const real total_length = fibers.totalLength();
    
    if ( space_volume <= 0 )
        throw InvalidParameter("Cannot equilibrate as Space:volume == 0");

    const real bind_rate1 = cop->hand1_prop->binding_rate;
    const real bind_range1 = cop->hand1_prop->binding_range;
    const real unbind_rate1 = cop->hand1_prop->unbinding_rate;

    const real bind_rate2 = cop->hand2_prop->binding_rate;
    const real bind_range2 = cop->hand2_prop->binding_range;
    const real unbind_rate2 = cop->hand2_prop->unbinding_rate;

    // get all crosspoints:
    fibers.allIntersections(loc1, loc2, std::max(bind_range1, bind_range2));
    const size_t num_crossings = loc1.size();
    
    const real ratio_fibs1 = 2 * total_length * bind_range1 / space_volume;
    const real ratio_fibs2 = 2 * total_length * bind_range2 / space_volume;
    const real ratio_cros1 = 4 * M_PI * num_crossings * square(bind_range1) / space_volume;
    const real ratio_cros2 = 4 * M_PI * num_crossings * square(bind_range2) / space_volume;
    
    /*
     The different states are defined in Belmonte et al. 2017, supplementary:
     Free, Bridge, Attached in location that cannot bridge, G=Attached near crosspoint
     */
    real BsG1 = bind_rate1 / unbind_rate1;
    real BsG2 = bind_rate2 / unbind_rate2;
    real A1sF = ( ratio_fibs1 - ratio_cros1 ) * BsG1 / 2;
    real A2sF = ( ratio_fibs2 - ratio_cros2 ) * BsG2 / 2;
    real G1sF = ratio_cros1 * BsG1 / 2;
    real G2sF = ratio_cros2 * BsG2 / 2;
    
    // the two should be equal
    real BsF = 0.5 * ( BsG1 * G1sF + BsG2 * G2sF );

    real popF = can.size() / ( 1.0 + A1sF + A2sF + G1sF + G2sF + BsF );
    real popA1 = A1sF * popF;
    real popA2 = A2sF * popF;
    real popG1 = G1sF * popF;
    real popG2 = G2sF * popF;
    real popB = BsF * popF;

#if ( 0 ) && ( DIM == 2 )
    printf("Couple::equilibrate %s:\n", cop->name_str());
    printf("     total %lu\n", reserve.size());
    const real num_fibers = fibers.size();
    const real fiber_length = total_length / num_fibers;
    const real nbc = num_fibers * ( num_fibers - 1 ) * square(fiber_length) / ( M_PI * space_volume );
    //const real nbc = square(total_length) / ( M_PI * space_volume );
    printf("     num_crossings predicted  %9.2f   true %9i\n", nbc, num_crossings);
    printf("     F %9.2f A %9.2f G %9.2f B %9.2f\n", popF, popA1+popA2, popG1+popG2, popB);
#endif
    
    if ( !can.empty() )
    {
        // distribute Couples at filament's crosspoints:
        uniAttach12(loc1, loc2, can, RNG.poisson(popB));

        const real dis1 = total_length / ( popA1 + popG1 );
        fibers.uniFiberSites(loc1, dis1);
        uniAttach1(loc1, can);

        const real dis2 = total_length / ( popA2 + popG2 );
        fibers.uniFiberSites(loc2, dis2);
        uniAttach2(loc2, can);
    }
}

/**
Distributes Couples for which `trans_activated!=true` on the filaments
*/
void CoupleSet::equilibrate(FiberSet const& fibers, PropertyList const& properties)
{
    for ( Property * i : properties.find_all("couple") )
    {
        CoupleProp * cop = static_cast<CoupleProp *>(i);
        cop->complete(simul_);
        
        if ( !cop->trans_activated )
        {
            CoupleList can;
            
            // collect all Couple of this kind:
            Couple * c = firstFF(), * nxt;
            while ( c )
            {
                nxt = c->next();
                if ( c->property() == cop )
                {
                    ffList.pop(c);
                    can.push_back(c);
                }
                c = nxt;
            }
            if ( can.size() > 0 )
            {
                equilibrate(fibers, can, cop);
                
                // release all collected Couple
                for ( Couple * cx : can )
                    link(cx);
                can.clear();
            }
        }
    }
    printf("Couple::equilibrate    FF %lu FA %lu AF %lu AA %lu\n", sizeFF(), sizeFA(), sizeAF(), sizeAA());
}


void CoupleSet::equilibrate()
{
    equilibrate(simul_.fibers, simul_.properties);
}


/**
 This takes all the Free Couple and attach them at the intersection points of the network of filaments
 */
void CoupleSet::bindToIntersections(FiberSet const& fibers, CoupleList& can, PropertyList const& properties)
{
    // calculate maximum range of Hands
    real range = 0;
    for ( Property const* i : properties.find_all("couple") )
    {
        CoupleProp const* cop = static_cast<CoupleProp const*>(i);
        range = std::max(range, cop->hand1_prop->binding_range);
        range = std::max(range, cop->hand2_prop->binding_range);
    }
    
    if ( range <= 0 )
        throw InvalidParameter("cannot connect network!");
    
    // get all crosspoints within this range:
    Array<FiberSite> loc1(1024, 1024), loc2(1024, 1024);
    fibers.allIntersections(loc1, loc2, range);
    const size_t nbc = loc1.size();
    assert_true(nbc == loc2.size());
    
    Cytosim::log << "Connect on " << nbc << " intersections within range " << range << "\n";
    
    uniAttach12(loc1, loc2, can, can.size());
}


/** applies to all couples of given class */
void CoupleSet::bindToIntersections(CoupleProp const* cop)
{
    CoupleList can;
    Couple * c = firstFF(), * nxt;
    while ( c )
    {
        nxt = c->next();
        if ( cop == c->property() )
        {
            ffList.pop(c);
            can.push_back(c);
        }
        c = nxt;
    }

    bindToIntersections(simul_.fibers, can, simul_.properties);
}


/** applies to all free couples */
void CoupleSet::bindToIntersections()
{
    CoupleList can;
    Couple * c = firstFF(), * nxt;
    while ( c )
    {
        nxt = c->next();
        ffList.pop(c);
        can.push_back(c);
        c = nxt;
    }

    bindToIntersections(simul_.fibers, can, simul_.properties);
}
