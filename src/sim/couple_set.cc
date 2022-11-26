// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

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


/// templated member function pointer...
/**
 In the loops we get the 'next' in the list always before calling 'FUNC', since
 'step()' may transfer the node to another list, changing the value of 'next()'
 */
template < void (Couple::*FUNC)() >
static inline void step_couples(Couple * obj, bool odd)
{
    Couple * nxt;
    if ( odd )
    {
        nxt = obj->next();
        (obj->*FUNC)();
        obj = nxt;
    }
    // this loop is unrolled, processing objects 2 by 2:
    while ( obj )
    {
        nxt = obj->next();
        (obj->*FUNC)();
        obj = nxt->next();
        (nxt->*FUNC)();
    }
}


/**
 Transfer free Couple with `fast_diffusion` to the reserves, starting from `obj`,
 and execute FUNC for all other Couples
*/
template < void (Couple::*FUNC)() >
void CoupleSet::step_collect(Couple * obj)
{
    Couple * nxt;
    while ( obj )
    {
        nxt = obj->next();
        CoupleProp const* p = static_cast<CoupleProp const*>(obj->property());
        if ( p->fast_diffusion )
        {
            ffList.pop(obj);
            uniReserves[p->number()-1].second.push_back(obj);
        }
        else
            (obj->*FUNC)();
        obj = nxt;
    }
}


/**
 CoupleSet::step() must call the appropriate Couple::step() exactly once
 for each Couple: either stepFF(), stepFA(), stepAF() or stepAA().
 
 The Couples are stored in multiple lists, and are automatically transferred
 from one list to another one if their Hands bind or unbind. The code relies
 on the fact that a Couple will be moved to the start of the list to which it
 is transferred, by 'push_front'. By starting always from the node that
 was first before any transfer could occur, we process each Couple only once.
 */
void CoupleSet::step()
{
    /*
    Cytosim::log("CoupleSet::step : FF %5i AF %5i FA %5i AA %5i\n",
                 ffList.size(), afList.size(), faList.size(), aaList.size());
    */
    
    Couple *const ffHead = firstFF();
    Couple *const afHead = firstAF();
    Couple *const faHead = firstFA();
    
    bool const aaOdd = aaList.size() & 1;
    bool const faOdd = faList.size() & 1;
    bool const afOdd = afList.size() & 1;
    bool const ffOdd = ffList.size() & 1;
    
    step_couples<&Couple::stepAA>(firstAA(), aaOdd);
    if ( RNG.flip() )
    {
        step_couples<&Couple::stepFA>(faHead, faOdd);
        step_couples<&Couple::stepAF>(afHead, afOdd);
    }
    else
    {
        step_couples<&Couple::stepAF>(afHead, afOdd);
        step_couples<&Couple::stepFA>(faHead, faOdd);
    }
    
    // use alternative attachment strategy:
    if ( uniEnabled )
    {
        step_collect<&Couple::stepFF>(ffHead);
        uniAttach(simul_.fibers);
    }
    else
    {
        //std::clog << "CoupleSet::step : FF " << ffList.size() << " head " << ffHead << '\n';
        // this loop is unrolled, processing objects 2 by 2:
        step_couples<&Couple::stepFF>(ffHead, ffOdd);
    }

    //printf("  : %lu couples [ %u %u ]\n", size(), inventory_.lowest(), inventory_.highest());
}


/**
 This version does not simulate the attachment of free Hand, and hence calls
 specialized versions of Couple::step() that do not include Hand::stepUnattached():
 either stepHand1(), stepHand2() or stepAA().
 Couple::stepFF() is never called.
 
 This is only used if POOL_UNATTACHED > 1
 */
void CoupleSet::stepSkipUnattached()
{
    /*
    Cytosim::log("CoupleSet::stepSkipUnattached : FF %5i AF %5i FA %5i AA %5i\n",
                 ffList.size(), afList.size(), faList.size(), aaList.size());
    */
    
    Couple *const ffHead = firstFF();
    Couple *const afHead = firstAF();
    Couple *const faHead = firstFA();
    
    bool const aaOdd = aaList.size() & 1;
    bool const faOdd = faList.size() & 1;
    bool const afOdd = afList.size() & 1;
    
    step_couples<&Couple::stepAA>(firstAA(), aaOdd);
    if ( RNG.flip() )
    {
        step_couples<&Couple::stepHand2>(faHead, faOdd);
        step_couples<&Couple::stepHand1>(afHead, afOdd);
    }
    else
    {
        step_couples<&Couple::stepHand1>(afHead, afOdd);
        step_couples<&Couple::stepHand2>(faHead, faOdd);
    }

    // use alternative attachment strategy:
    if ( uniEnabled )
    {
        step_collect<&Couple::stepBlank>(ffHead);
        uniAttach(simul_.fibers);
    }
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


Object * CoupleSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Couple::TAG )
    {
        CoupleProp * p = simul_.findProperty<CoupleProp>("couple", pid);
        return p->newCouple();
    }
    throw InvalidIO("Warning: unknown Couple tag `"+std::to_string(tag)+"'");
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
ObjectList CoupleSet::newObjects(const Property* p, Glossary& opt)
{
    CoupleProp const* pp = static_cast<CoupleProp const*>(p);
    Couple * obj = pp->newCouple(&opt);
        
    // Allow user to attach hand1:
    if ( opt.has_key("attach1") )
        obj->attach1(simul_.fibers.someSite("attach1", opt));

    // Allow user to attach hand2:
    if ( opt.has_key("attach2") )
        obj->attach2(simul_.fibers.someSite("attach2", opt));

#if BACKWARD_COMPATIBILITY < 100
    /* It would be possible to create Couple with custom hand type, and the
    syntax below to attach the Hands could be better used for this */
    
    // Allow user to attach hand1:
    if ( opt.has_key("site1") )
        obj->attach1(simul_.fibers.someSite("site1", opt));
    
    // Allow user to attach hand2:
    if ( opt.has_key("site2") )
        obj->attach2(simul_.fibers.someSite("site2", opt));
#endif
    
    return ObjectList(obj);
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
    assert_true( obj->objset() == this );

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
    assert_true(obj->objset() == nullptr);
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
    if ( ffList.size() > 1 ) ffList.shuffle();
    if ( afList.size() > 1 ) afList.shuffle();
    if ( faList.size() > 1 ) faList.shuffle();
    if ( aaList.size() > 1 ) aaList.shuffle();
    
    ObjectID id = RNG.pint32(inventory_.lowest(), inventory_.highest());
    Couple * c = static_cast<Couple*>(inventory_.get(id));
    if ( c )
    switch( c->state() )
    {
        case 0: if ( !uniEnabled ) ffList.permute(c);break;
        case 1: afList.permute(c); break;
        case 2: faList.permute(c); break;
        case 3: aaList.permute(c); break;
    }
}


void CoupleSet::erase()
{
    for ( CoupleReserve & reserve : uniReserves )
    {
        for ( Couple * c : reserve.second )
        {
            c->objset(nullptr);
            delete(c);
        }
        reserve.second.clear();
    }
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


void CoupleSet::reheat()
{
    //std::clog << "Couple::reheat " << ice_.size() << "\n";
    Object * i = ice_.pop_front();
    while ( i )
    {
        Couple* o = static_cast<Couple*>(i);
        i = ice_.pop_front();
        o->hand1()->detachHand();
        o->hand2()->detachHand();
        o->randomizePosition();
        linkFF(o);
    }
}


/* After reading from file, the Hands should not
 update any Fiber, Single or Couple as they will be deleted */
void CoupleSet::prune()
{
    Object * i = ice_.pop_front();
    while ( i )
    {
        Couple* o = static_cast<Couple*>(i);
        i = ice_.pop_front();
        inventory_.unassign(o);
        o->objset(nullptr);
        delete(o);
    }
}


//------------------------------------------------------------------------------
#pragma mark -


void CoupleSet::writeAA(Outputter& out) const
{
    out.write("\n#section couple AA");
    writeObjects(out, aaList);
}

void CoupleSet::writeAF(Outputter& out) const
{
    out.write("\n#section couple AF");
    writeObjects(out, afList);
}

void CoupleSet::writeFA(Outputter& out) const
{
    out.write("\n#section couple FA");
    writeObjects(out, faList);
}

void CoupleSet::writeFF(Outputter& out) const
{
    out.write("\n#section couple FF");
    writeObjects(out, ffList);
}

void CoupleSet::writeFF_skip(Outputter& out) const
{
    out.write("\n#section couple FF");
    // record number of objects:
    size_t cnt = ffList.size();
    size_t sup = inventory_.highest();
    out.write("\n#record "+std::to_string(cnt)+" "+std::to_string(sup));
    if ( out.binary() ) out.put_char('\n');
    
    for ( Couple const* n=firstFF(); n; n=n->next() )
    {
        if ( ! n->prop->fast_diffusion )
            n->write(out);
    }
    out.write("\n#section couple reheat");
}

void CoupleSet::writeSet(Outputter& out) const
{
    if ( sizeAA() > 0 )
        writeAA(out);
    if ( sizeAF() > 0 )
        writeAF(out);
    if ( sizeFA() > 0 )
        writeFA(out);
    if ( sizeFF() > 0 )
    {
        if ( simul_.prop.skip_free_couple && skip_now )
            writeFF_skip(out);
        else
        {
            writeFF(out);
            skip_now = 1;
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
    int code = ffList.bad() | afList.bad() | faList.bad() | aaList.bad();
    if ( code )
        return code;

    Couple * obj;
    for ( obj=firstFF(); obj ; obj = obj->next() )
    {
        if ( obj->attached1() || obj->attached2() )
            code |= 8;
    }
    
    for ( obj=firstAF(); obj ; obj = obj->next() )
    {
        if ( !obj->attached1() || obj->attached2() )
            code |= 16;
        if ( obj->hand1()->bad() )
            code |= 128;
    }
    
    for ( obj=firstFA(); obj ; obj = obj->next() )
    {
        if ( obj->attached1() || !obj->attached2() )
            code |= 32;
        if ( obj->hand2()->bad() )
            code |= 128;
    }
    
    for ( obj=firstAA(); obj ; obj = obj->next() )
    {
        if ( !obj->attached1() || !obj->attached2() )
            code |= 64;
        if ( obj->hand1()->bad() )
            code |= 256;
        if ( obj->hand2()->bad() )
            code |= 256;
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
        c->objset(this);
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
            linkFF(c);
            h->attach(i);
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
            linkFF(c);
            h->attach(i);
        }
    }
}


/**
 Distribute up to `nb` Couples from `can` by attaching them to locations specified
 by `loc1` and `loc2`. These positions correspond to fibers crossings, found by
 FiberSet::allIntersections().
 The Couples are distributed randomly on the crosspoints
 */
void CoupleSet::uniAttach12(Array<FiberSite>& loc1, Array<FiberSite>& loc2,
                            CoupleList& can, size_t nb)
{
    const size_t nbc = loc1.size();
    assert_true(nbc == loc2.size());
    
    if ( nbc < 1 )
        return;
    
    const size_t sup = std::min(nb, can.size());

    for ( size_t n = 0; n < sup; ++n )
    {
        Couple * c = can.back();
        can.pop_back();
        linkFF(c);
        // pick randomly with replacement:
        size_t p = RNG.pint32(nbc);
        c->attach1(loc1[p]);
        c->attach2(loc2[p]);
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
    Array<FiberSite> loc(128, 128);
    
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
        bool fixed = ( p->fast_reservoir > 0 );
        size_t cnt = ( fixed ? p->fast_reservoir : can.size());

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
        assert_true(p->number() > 0);
        res |= p->fast_diffusion;
    }
    
    if ( res )
    {
        uniReserves.resize(last);
        for ( CoupleReserve & reserve : uniReserves )
            reserve.first = nullptr;
        for ( Property const* i : allprop )
        {
            CoupleProp const* p = static_cast<CoupleProp const*>(i);
            if ( p->fast_diffusion )
                uniReserves[p->number()-1].first = p;
        }
    }
    
    return res;
}


/**
 Move all Couples from the reserves to the list of unbound Couples
 */
void CoupleSet::uniRelax()
{
    for ( CoupleReserve & reserve : uniReserves )
    {
        for ( Couple * c : reserve.second )
        {
            assert_true(!c->attached1() && !c->attached2());
            c->randomizePosition();
            linkFF(c);
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
                for ( Couple * x : can )
                    linkFF(x);
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
    
    Cytosim::log << "Connect on " << loc1.size() << " intersections within range " << range << "\n";
    
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
