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

//------------------------------------------------------------------------------

void CoupleSet::prepare()
{
    uniPrepare(simul_.properties);
}


/// templated member function pointer...
/**
 In the loops we get the 'next' in the list always before calling 'FUNC', since
 'Couple::step()' may transfer the node to another list, changing the value of 'next()'
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
 Call stepFF() for all free Couple starting from `obj`, except that if
 `fast_diffusion==true`, the Couple is transferred to corresponding reserve list
*/
void CoupleSet::uniStepCollect(Couple * obj)
{
    Couple * nxt;
    while ( obj )
    {
        nxt = obj->next();
        CoupleProp const* P = obj->prop;
        if ( P->fast_diffusion > 0 )
        {
            ffList.pop(obj);
            inventory_.unassign(obj);
            P->stocks.push(obj);
            ++P->uni_counts;
        }
        else
            obj->stepFF();
        obj = nxt;
    }
}


/**
 CoupleSet::steps() must call the appropriate Couple::step() exactly once
 for each Couple: either stepFF(), stepFA(), stepAF() or stepAA().
 
 The Couples are stored in multiple lists, and are automatically transferred
 from one list to another one if their Hands bind or unbind. The code relies
 on the fact that a Couple will be moved to the start of the list to which it
 is transferred, by 'push_front'. By starting always from the node that
 was first before any transfer could occur, we process each Couple only once.
 */
void CoupleSet::steps()
{
    /*
    Cytosim::log.print("CoupleSet::step : FF %5i AF %5i FA %5i AA %5i\n",
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
    if ( uniCouples.size() )
    {
        uniStepCollect(ffHead);
        uniAttach(simul_.fibers);
    }
    else
    {
        //std::clog << "CoupleSet::step : FF " << ffList.size() << " head " << ffHead << '\n';
        // this loop is unrolled, processing objects 2 by 2:
        step_couples<&Couple::stepFF>(ffHead, ffOdd);
    }
    
#if 0
    ObjectID h = inventory_.highest();
    if ( h > 4096 && h > 2 * ( size() + countStocks() ) )
    {
        uniRelax();
        inventory_.reassign();
        Cytosim::log("Single::reassign(", h, " ---> ", inventory_.highest(), ")\n");
    }
#endif
    //printf(" couples: %lu %lu [ %u %u ]\n", size(), inventory_.count(), inventory_.lowest(), inventory_.highest());
    if ( size() > 1 ) shuffle();
}


/**
 This version does not simulate the attachment of free Hand, and hence calls
 Couple::stepHand1() and stepHand2() that do not perform attachment.
 
 This is only used if POOL_UNATTACHED > 1
 */
void CoupleSet::stepsSkippingUnattached()
{
    /*
    Cytosim::log.print("CoupleSet::stepsSkippingUnattached : FF %5i AF %5i FA %5i AA %5i\n",
                        ffList.size(), afList.size(), faList.size(), aaList.size());
    */
    
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
    if ( size() > 1 ) shuffle();
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


void CoupleSet::addFreeCouple(Couple * obj)
{
    assert_true(!obj->attached());
    assert_true(obj->objset()==nullptr || obj->objset()==this);
    obj->objset(this);
    ffList.push_back(obj);
    inventory_.assign(obj);
    //std::clog << "addFreeCouple(" << obj->reference() << ")\n";
}


Couple * CoupleSet::addFreeCouple(CoupleProp const* P, Vector const& pos)
{
    Couple * C = makeCouple(P);
    C->setPosition(pos);
    addFreeCouple(C);
    return C;
}


Object * CoupleSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Couple::TAG || tag == Couple::DUO_TAG )
    {
        CoupleProp * P = simul_.findProperty<CoupleProp>("couple", pid);
        return makeCouple(P);
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
ObjectList CoupleSet::newObjects(Property const* p, Glossary& opt)
{
    CoupleProp const* pp = static_cast<CoupleProp const*>(p);
    Couple * obj = makeCouple(pp);
    
    // possibly activate
    int a = 0;
    if ( opt.set(a, "active") && a )
        obj->activate();

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


/// pick from reserves if possible
Couple * CoupleSet::makeCouple(CoupleProp const* P)
{
    Couple * C = P->stocks.head();
    if ( C )
    {
        //std::clog << "stock -> " << C->reference() << " ";
        P->stocks.pop();
    } else {
        C = P->newCouple();
        //std::clog << "new " << C->reference() << " ";
    }
    return C;
}


void CoupleSet::makeCouples(CoupleProp const* P, size_t cnt)
{
    while ( cnt-- > 0 )
    {
        Couple * C = makeCouple(P);
        C->randomizePosition();
        addFreeCouple(C);
    }
}


void CoupleSet::makeCouples(size_t cnt[], PropertyID n_cnt)
{
    for ( PropertyID i = 1; i < n_cnt; ++i )
    {
        if ( cnt[i] > 0 )
        {
            Property * cp = simul_.properties.find("couple", i);
            if ( cp )
            {
                CoupleProp * P = static_cast<CoupleProp*>(cp);
                // renew pointers to 'confine_space'
                P->complete(simul_);
                makeCouples(P, cnt[i]);
            }
        }
    }
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
    assert_true( obj->tag() == Couple::TAG || obj->tag() == Couple::DUO_TAG );
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
    if ( c->attached1() ) c->hand1()->detach();
    if ( c->attached2() ) c->hand2()->detach();
    ffList.pop(obj);
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
    ObjectID id = RNG.pint32(inventory_.lowest(), inventory_.highest());
    Couple * c = static_cast<Couple*>(inventory_.get(id));
    if ( c )
    switch( c->state() )
    {
        case 0: if ( c->prop->fast_diffusion <= 0 ) ffList.shuffle(c); break;
        case 1: afList.shuffle(c); break;
        case 2: faList.shuffle(c); break;
        case 3: aaList.shuffle(c); break;
    }
    
    id = RNG.pint32(inventory_.lowest(), inventory_.highest());
    c = static_cast<Couple*>(inventory_.get(id));
    if ( c )
    switch( c->state() )
    {
        case 0: if ( c->prop->fast_diffusion <= 0 ) ffList.permute(c); break;
        case 1: afList.permute(c); break;
        case 2: faList.permute(c); break;
        case 3: aaList.permute(c); break;
    }
}


void CoupleSet::erase()
{
    for ( Property const* i : simul_.properties.find_all("couple") )
    {
        CoupleProp const * P = static_cast<CoupleProp const*>(i);
        P->stocks.erase();
        P->uni_counts = 0;
    }
    ObjectSet::erasePool(aaList);
    ObjectSet::erasePool(faList);
    ObjectSet::erasePool(afList);
    ObjectSet::erasePool(ffList);
    inventory_.clear();
}


void CoupleSet::detachAll()
{
    for ( Couple * C=firstAA(); C; C=C->next() )
    {
        C->hand1()->detachHand();
        C->hand2()->detachHand();
    }
    ffList.grab(aaList);
    assert_true(aaList.empty());
    
    for ( Couple * C=firstAF(); C; C=C->next() )
        C->hand1()->detachHand();
    ffList.grab(afList);
    assert_true(afList.empty());
    
    for ( Couple * C=firstFA(); C; C=C->next() )
        C->hand2()->detachHand();
    ffList.grab(faList);
    assert_true(faList.empty());
}


void CoupleSet::defrostStore()
{
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        inventory_.unassign(i);
        Couple * C = static_cast<Couple*>(i);
        if ( C->hand1()->attached() )
            C->hand1()->detachHand();
        if ( C->hand2()->attached() )
            C->hand2()->detachHand();
        C->prop->stocks.push(C);
    }
    //infoStocks(std::clog);
}

//------------------------------------------------------------------------------
#pragma mark -


void CoupleSet::freeze()
{
    assert_true(ice_.empty());
    ice_.grab(aaList);
    ice_.grab(faList);
    ice_.grab(afList);
    ice_.grab(ffList);
}


/** cnt[i] is the number of Couples of type `i` to be released */
void CoupleSet::reheat(size_t cnt[], PropertyID n_cnt)
{
#if 0
    std::clog << "Couple::reheat";
    for ( size_t n = 0; n < n_cnt; ++n )
        std::clog << " " << cnt[n];
    std::clog << "\n";
#endif
    //std::clog << "Couple::reheat " << ice_.size() << "\n";
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        Couple* C = static_cast<Couple*>(i);
        // we want to skip the 'beforeDetachment' here:
        if ( C->hand1()->attached() )
            C->hand1()->detachHand();
        if ( C->hand2()->attached() )
            C->hand2()->detachHand();
        PropertyID id = C->prop->number();
        if ( id < n_cnt && 0 < cnt[id] )
        {
            --cnt[id];
            if ( C->prop->fast_diffusion > 0 )
                C->randomizePosition();
            addFreeCouple(C);
        }
        else
        {
            // place C on reserve
            inventory_.unassign(C);
            C->prop->stocks.push(C);
        }
    }
}


void CoupleSet::reheat()
{
    size_t sup = inventory_.capacity();
    size_t cnt[16] = { 0 };
    for ( int i = 0; i < 16; ++i ) cnt[i] = sup;
    reheat(cnt, 16);
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 This will save some of the objects in the normal way, using `write()`
 and will also save a count of the objects that were skipped
 */
void CoupleSet::writeSomeObjects(Outputter& out) const
{
    writeRecords(out, ffList.size(), inventory_.highest());
    
    std::map<PropertyID, size_t> cnt;
    // count all the elements that are virtually present:
    for ( CoupleProp const* P : uniCouples )
        cnt[P->number()] = P->uni_counts;
    // write Couple if `save_unbound > 0`:
    for ( Couple const* n=firstFF(); n; n=n->next() )
    {
        if ( n->prop->save_unbound )
            n->write(out);
        else
            ++cnt[n->prop->number()];
    }
    if ( !cnt.empty() )
    {
        // get highest prop ID that was skipped:
        PropertyID sup = cnt.rbegin()->first;
        if ( sup > 0 )
        {
            // write counts for each class of unwritten Couple:
            out.write("\n#section couple reheat");
            for ( PropertyID i = 0; i <= sup; ++i )
                out.writeUInt(cnt[i], ' ');
        }
    }
    // decrement `save_unbound`:
    for ( Property * i : simul_.properties.find_all("couple") )
    {
        CoupleProp * P = static_cast<CoupleProp *>(i);
        P->save_unbound -= ( P->save_unbound > 0 );
    }
}


void CoupleSet::writeSet(Outputter& out) const
{
    if ( sizeAA() > 0 )
    {
        out.write("\n#section couple AA");
        writePool(out, aaList);
    }
    if ( sizeAF() > 0 )
    {
        out.write("\n#section couple AF");
        writePool(out, afList);
    }
    if ( sizeFA() > 0 )
    {
        out.write("\n#section couple FA");
        writePool(out, faList);
    }
    if ( sizeFF() > 0 )
    {
        out.write("\n#section couple FF");
        writeSomeObjects(out);
        //writePool(out, ffList);
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
    size_t ff = ffList.count(func, arg);
    size_t af = afList.count(func, arg);
    size_t fa = faList.count(func, arg);
    size_t aa = aaList.count(func, arg);
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
    int err = 0;
    Couple * obj;
    size_t cnt = sizeFF();
    for ( obj=firstFF(); obj ; obj=obj->next() )
    {
        if ( obj->attached1() || obj->attached2() )
            err |= 8;
        if ( cnt-- == 0 )
            return 1;
    }
    if ( cnt ) return 1;
    
    cnt = sizeAF();
    for ( obj=firstAF(); obj ; obj=obj->next() )
    {
        if ( !obj->attached1() || obj->attached2() )
            err |= 16;
        if ( simul_.fibers.badIdentity(obj->fiber1()) )
            err |= 512;
        if ( obj->hand1()->bad() )
            err |= 128;
        if ( cnt-- == 0 )
            return 2;
    }
    if ( cnt ) return 2;

    cnt = sizeFA();
    for ( obj=firstFA(); obj ; obj=obj->next() )
    {
        if ( obj->attached1() || !obj->attached2() )
            err |= 32;
        if ( simul_.fibers.badIdentity(obj->fiber2()) )
            err |= 512;
        if ( obj->hand2()->bad() )
            err |= 128;
        if ( cnt-- == 0 )
            return 3;
    }
    if ( cnt ) return 3;

    cnt = sizeAA();
    for ( obj=firstAA(); obj ; obj=obj->next() )
    {
        if ( !obj->attached1() || !obj->attached2() )
            err |= 64;
        if ( simul_.fibers.badIdentity(obj->fiber1()) )
            err |= 512;
        if ( simul_.fibers.badIdentity(obj->fiber2()) )
            err |= 512;
        if ( obj->hand1()->bad() )
            err |= 256;
        if ( obj->hand2()->bad() )
            err |= 256;
        if ( cnt-- == 0 )
            return 4;
    }
    if ( cnt ) return 4;

    return err;
}


//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion

size_t CoupleSet::countStocks() const
{
    size_t res = 0;
    for ( Property const* i : simul_.properties.find_all("couple") )
    {
        CoupleProp const * P = static_cast<CoupleProp const*>(i);
        res += P->stocks.size();
    }
    return res;
}


void CoupleSet::infoStocks(std::ostream& os) const
{
    os << "  Couple:stocks";
    for ( Property const* i : simul_.properties.find_all("couple") )
    {
        CoupleProp const * P = static_cast<CoupleProp const*>(i);
        size_t cnt = count(match_property, P);
        os << " " << P->number() << ": " << cnt << " ( " << P->stocks.size() << " " << P->uni_counts << " )";
    }
    os << "\n";
}


void CoupleSet::uniRefill(CoupleProp const* cop, size_t cnt)
{
    CoupleStock & can = cop->stocks;
    for ( size_t i = can.size(); i < cnt; ++i )
        can.push(cop->newCouple());
}


/**
 Attach Hand1 of exactly one Couple from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void CoupleSet::uniAttach1(FiberSiteList& loc, CoupleStock& can)
{
    // crop list to match available number of candidates:
    loc.shuffle_truncate(can.size());

    for ( FiberSite & i : loc )
    {
        Couple * C = can.head();
        Hand * h = C->hand1();
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            can.pop();
            addFreeCouple(C);
            h->attach(i);
        }
    }
}


/**
 Attach Hand2 of exactly one Couple from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void CoupleSet::uniAttach2(FiberSiteList& loc, CoupleStock& can)
{
    // crop list to match available number of candidates:
    loc.shuffle_truncate(can.size());

    for ( FiberSite & i : loc )
    {
        Couple * C = can.head();
        Hand * h = C->hand2();
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            can.pop();
            addFreeCouple(C);
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
void CoupleSet::uniAttach12(FiberSiteList& loc1, FiberSiteList& loc2,
                            CoupleStock& can, unsigned sup)
{
    const unsigned nbc = loc1.size();
    assert_true(nbc == loc2.size());
    
    if ( nbc < 1 )
        return;
    
    sup = std::min(sup, (unsigned)can.size());

    for ( unsigned n = 0; n < sup; ++n )
    {
        Couple * C = can.head();
        can.pop();
        addFreeCouple(C);
        // pick randomly with replacement:
        unsigned i = RNG.pint32(nbc);
        C->attach1(loc1[i]);
        C->attach2(loc2[i]);
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
    FiberSiteList loc(128);
    
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
    for ( CoupleProp const* P : uniCouples )
    {
        CoupleStock& can = P->stocks;
        // assuming (or not) a given number of diffusing molecules
        bool fixed = ( P->fast_reservoir > 0 );
        size_t cnt = ( fixed ? P->fast_reservoir : P->uni_counts );
        
        if ( cnt > 0 )
        {
            const real alpha = 2 * P->spaceVolume() / cnt;
            
            if ( P->fast_diffusion == 2 )
            {
                real dis = alpha / P->hand1_prop->bindingSectionRate();
                fibers.newFiberSitesP(loc, dis);
            }
            else
            {
                real dis = alpha / P->hand1_prop->bindingSectionProb();
                fibers.uniFiberSites(loc, dis);
            }
            
            if ( fixed ) // create enough candidates for all sites
                uniRefill(P, loc.size());
            
            size_t total = size();
            uniAttach1(loc, can);
            
            if ( !P->trans_activated )
            {
                // if ( couple:trans_activated == true ), Hand2 cannot bind
                if ( P->fast_diffusion == 2 )
                {
                    real dis = alpha / P->hand2_prop->bindingSectionRate();
                    fibers.newFiberSitesP(loc, dis);
                }
                else
                {
                    real dis = alpha / P->hand2_prop->bindingSectionProb();
                    fibers.uniFiberSites(loc, dis);
                }
                
                if ( fixed ) // create enough candidates for all sites
                    uniRefill(P, loc.size());
                
                uniAttach2(loc, can);
            }
            if ( !fixed )
                P->uni_counts -= size() - total;
        }
    }
}


/**
 Sets the list `uniCouples` with couple class concerned with `fast_diffusion`
*/
void CoupleSet::uniPrepare(PropertyList const& properties)
{
    uniCouples.clear();
    for ( Property const* i : simul_.properties.find_all("couple") )
    {
        CoupleProp const * P = static_cast<CoupleProp const*>(i);
        //assert_true( P->stocks.size() == P->stocks.recount() );
        if ( P->fast_diffusion > 0 )
            uniCouples.push_back(P);
    }
}


/**
 Move all Couples from the reserves to the list of unbound Couples
 */
void CoupleSet::uniRelax()
{
    for ( CoupleProp const* P : uniCouples )
    {
        makeCouples(P, P->uni_counts);
        P->uni_counts = 0;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Equilibration


void CoupleSet::equilibrateSym(FiberSet const& fibers, CoupleProp const* cop, size_t total)
{
    Array<FiberSite, 0> loc1(16), loc2(16);
    CoupleStock & can = cop->stocks;
    
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
    
    real popF = total / ( 1 + AsF + GsF + BsG * GsF );
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
    
    // distribute Couples at filament's crosspoints:
    uniAttach12(loc1, loc2, can, RNG.poisson(popB));
    
    real dis = 2 * total_length / ( popA + popG );
    
    fibers.uniFiberSites(loc1, dis);
    uniAttach1(loc1, can);
    
    fibers.uniFiberSites(loc2, dis);
    uniAttach2(loc2, can);
}


/**
 This attempts to create a configuration of Couple that is close to the equilibrium
 that would be reached, after a sufficient time is given for binding and unbinding.
 This assumes that the configuration of filaments does not change, and also that
 it is random, in particular without bundles. The motion of the motor is also ignored.
 */
void CoupleSet::equilibrate(FiberSet const& fibers, CoupleProp const* cop, CoupleStock& can, size_t total)
{
    FiberSiteList loc1(1024), loc2(1024);
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

    real popF = total / ( 1.0 + A1sF + A2sF + G1sF + G2sF + BsF );
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
    
    // distribute Couples at filament's crosspoints:
    uniAttach12(loc1, loc2, can, RNG.poisson(popB));
    
    const real dis1 = total_length / ( popA1 + popG1 );
    fibers.uniFiberSites(loc1, dis1);
    uniAttach1(loc1, can);
    
    const real dis2 = total_length / ( popA2 + popG2 );
    fibers.uniFiberSites(loc2, dis2);
    uniAttach2(loc2, can);
}

/**
Distributes Couples for which `trans_activated!=true` on the filaments
*/
void CoupleSet::equilibrate(FiberSet const& fibers, CoupleProp const* P)
{
    if ( !P->trans_activated )
    {
        CoupleStock can;
        
        // collect all Couple of this kind:
        Couple * C = firstFF(), * nxt;
        while ( C )
        {
            nxt = C->next();
            if ( C->property() == P )
            {
                ffList.pop(C);
                can.push(C);
            }
            C = nxt;
        }
        if ( can.head() )
        {
            equilibrate(fibers, P, P->stocks, P->stocks.size());
            
            // release all collected Couple
            while (( C = can.head() ))
            {
                can.pop();
                addFreeCouple(C);
            }
        }
    }
}


void CoupleSet::equilibrate()
{
    for ( Property * P : simul_.properties.find_all("couple") )
    {
        P->complete(simul_);
        equilibrate(simul_.fibers, static_cast<CoupleProp*>(P));
    }
    printf("Couple::equilibrate    FF %lu FA %lu AF %lu AA %lu\n", sizeFF(), sizeFA(), sizeAF(), sizeAA());
}


/**
 Calculate maximum range of Hands from all know Couple types
 */
real CoupleSet::maxBindingRange() const
{
    real res = 0;
    for ( Property const* i : simul_.properties.find_all("couple") )
    {
        CoupleProp const* P = static_cast<CoupleProp const*>(i);
        res = std::max(res, P->hand1_prop->binding_range);
        res = std::max(res, P->hand2_prop->binding_range);
    }
    return res;
}


/**
 Attach couples in given stock at the intersection points of given fiber set
 */
void CoupleSet::bindToIntersections(FiberSet const& fibers, CoupleStock& can, real rge)
{
    if ( rge <= 0 )
        throw InvalidParameter("cannot connect fibers, with null range!");
    
    // get all crosspoints within this range:
    FiberSiteList loc1(1024), loc2(1024);
    fibers.allIntersections(loc1, loc2, rge);
    
    Cytosim::log("Connect ", loc1.size(), " intersections within range ", rge, "\n");
    
    uniAttach12(loc1, loc2, can, can.size());
}


/** applies to all couples of given class */
void CoupleSet::bindToIntersections(CoupleProp const* cop)
{
    CoupleStock can;
    Couple * C = firstFF(), * nxt;
    while ( C )
    {
        nxt = C->next();
        if ( cop == C->property() )
        {
            ffList.pop(C);
            can.push(C);
        }
        C = nxt;
    }

    const real rge = maxBindingRange();
    bindToIntersections(simul_.fibers, can, rge);
}


/** applies to all free couples */
void CoupleSet::bindToIntersections()
{
    CoupleStock can;
    Couple * C = firstFF(), * nxt;
    while ( C )
    {
        nxt = C->next();
        ffList.pop(C);
        can.push(C);
        C = nxt;
    }
    
    const real rge = maxBindingRange();
    bindToIntersections(simul_.fibers, can, rge);
}
