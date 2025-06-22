// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#include "single_set.h"
#include "single_prop.h"
#include "glossary.h"
#include "iowrapper.h"
#include "messages.h"

#include "simul.h"
#include "simul_prop.h"
#include "property_list.h"
#include "wrist.h"
#include "wrist_long.h"



//------------------------------------------------------------------------------

void SingleSet::prepare()
{
    uniPrepare(simul_.properties);
}


/// templated member function pointer...
/**
In the loops we get the 'next' in the list always before calling 'FUNC', since
'steps()' may transfer the node to another list, changing the value of 'next()'
*/
template < void (Single::*FUNC)() >
static inline void step_singles(Single * obj, bool odd)
{
    Single * nxt;
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
 Call stepF() for all free Single starting from `obj`, except that if
 `fast_diffusion > 0`, the Single is transferred to corresponding reserve list
*/
void SingleSet::uniStepCollect(Single * obj)
{
    Single * nxt;
    while ( obj )
    {
        nxt = obj->next();
        SingleProp const* P = obj->prop;
        if ( P->fast_diffusion > 0 && !obj->base() )
        {
            fList.pop(obj);
            inventory_.unassign(obj);
            P->stocks.push(obj);
            ++P->uni_counts;
        }
        else
            obj->stepF();
        obj = nxt;
    }
}


/**
SingleSet::steps() must call the appropriate Single::step() exactly once
for each Single: either stepF() or stepA().

The Singles are stored in two lists, and are automatically transferred
from one list to the other if their Hands bind or unbind. The code relies
on the fact that a Single will be moved to the start of the list to which it
is transferred, by 'push_front'. By starting always from the node that
was first before any transfer could occur, we process each Couple only once.
*/
void SingleSet::steps()
{
    //Cytosim::log.print("SingleSet: F %5lu A %5lu\n", sizeF(), sizeA());
    
    Single *const fHead = firstF();
    bool fOdd = sizeF() & 1;
    
    step_singles<&Single::stepA>(firstA(), sizeA() & 1);
    
    // use alternative attachment strategy:
    if ( uniSingles.size() )
    {
        uniStepCollect(fHead);
        uniAttach(simul_.fibers);
    }
    else
    {
        step_singles<&Single::stepF>(fHead, fOdd);
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
    if ( aList.size() > 1 ) aList.shuffle();
    if ( fList.size() > 1 ) fList.shuffle();
}


/**
 This version does not simulate the attachment of free Hand, and hence skips
 `Single::stepF()` that performs attachment

 This is only used if POOL_UNATTACHED > 1
*/
void SingleSet::stepsSkippingUnattached()
{
    //Cytosim::log.print("SingleSet::stepsSkippingUnattached : F %5i A %5i\n", fList.size(), aList.size());
    
    step_singles<&Single::stepA>(firstA(), sizeA() & 1);
    if ( aList.size() > 1 ) aList.shuffle();
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 @copydetails SingleGroup
 */
Property* SingleSet::newProperty(const std::string& cat, const std::string& nom, Glossary& opt) const
{
    if ( cat == "single" )
        return new SingleProp(nom);
    else
        return nullptr;
}


void SingleSet::addFreeSingle(Single * obj)
{
    assert_true(!obj->attached());
    assert_true(obj->objset()==nullptr || obj->objset()==this);
    obj->objset(this);
    fList.push_back(obj);
    inventory_.assign(obj);
    //std::clog << "addFreeSingle(" << obj->reference() << ")\n";
}


Single * SingleSet::addFreeSingle(SingleProp const* P, Vector const& pos)
{
    Single * S = makeSingle(P);
    S->setPosition(pos);
    addFreeSingle(S);
    return S;
}


Object * SingleSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Single::TAG )
    {
        SingleProp * P = simul_.findProperty<SingleProp>("single", pid);
        return makeSingle(P);
    }
    else if ( tag == Single::WRIST_TAG )
    {
        SingleProp * P = simul_.findProperty<SingleProp>("single", pid);
        return P->newWrist(nullptr, 0);
    }
    throw InvalidIO("unknown Single tag `"+std::to_string(tag)+"'");
    return nullptr;
}

/**
 @addtogroup SingleGroup
 
 A newly created Single can be anchored to a Mecable:
 
     new NAME {
       base = OBJECT, POINT
     }
 
 where:
 - OBJECT is the concatenation of the class name with the serial number of the object:
     - 'bead1' for the first bead
     - 'bead2' for the second...
     - 'bead0' designates the last bead made,
     - 'bead-1' is the penultimate one, etc.
     .
 - POINT designates a point on this object:
     - point1 = first point
     - point2 = second point...
     .
 .

 You can attach a Single to a fiber:
 
     new simplex
     {
        attach = FIBER, REAL, REFERENCE
     }
 
 where:
 - FIBER designates the fiber:
     - `fiber1` of `fiber2` correspond to fibers directly
     - `first` or `last` to the oldest and youngest fiber
     - `fiber-1` the penultimate, etc.
     .
 - REAL is the abscissa of the attachment point.
   If the abscissa is not specified, and random position along
   along the fiber will be selected.
 - REFERENCE can be `minus_end`, `center` or `plus_end` (default = `origin`).
   This defines from which position the abscissa is measured.
 .
 */
ObjectList SingleSet::newObjects(Property const* p, Glossary& opt)
{
    SingleProp const* pp = static_cast<SingleProp const*>(p);
    Single * obj = nullptr;
    std::string str;
    if ( opt.set(str, "base") )
    {
        Mecable * mec = simul_.pickMecable(str);
        if ( !mec )
            throw InvalidParameter("could not find Mecable specified in single:base `"+str+"'");
        // get index of point in second argument
        index_t ip = 0;
        if ( opt.set(str, "base", 1) )
            ip = mec->point_index(str);
         
        obj = pp->newWrist(mec, ip);
    }
    else
        obj = makeSingle(pp);

    // Allow user to attach Hand to an existing fiber
    if ( opt.has_key("attach") )
        obj->attach(simul_.fibers.someSite("attach", opt));
    
    // Allow user to attach Hand to an existing fiber
    if ( opt.has_key("site") )
        obj->attach(simul_.fibers.someSite("site", opt));

    return ObjectList(obj);
}


//------------------------------------------------------------------------------
#pragma mark -


/// pick from reserves if possible
Single * SingleSet::makeSingle(SingleProp const* P)
{
    Single * S = P->stocks.head();
    if ( S )
        P->stocks.pop();
    else
        S = P->newSingle();
    //std::clog << "makeSingle(" << S->reference() << ") ";
    return S;
}


void SingleSet::makeSingles(SingleProp const* P, size_t cnt)
{
    reserve(inventory_.highest()+cnt);
    while ( cnt-- > 0 )
    {
        Single * S = makeSingle(P);
        S->randomizePosition();
        addFreeSingle(S);
    }
}


void SingleSet::makeSingles(size_t cnt[], PropertyID n_cnt)
{
    //std::clog << "makeSingles " << cnt[0] << " " << cnt[1] << " " << cnt[2] << "\n";
    // note that id=0 is invalid
    for ( PropertyID i = 1; i < n_cnt; ++i )
    {
        Property * sp = simul_.properties.find("single", i);
        if ( sp )
        {
            SingleProp * P = static_cast<SingleProp*>(sp);
            // renew pointers to 'confine_space'
            P->complete(simul_);
            makeSingles(P, cnt[i]);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -


void SingleSet::relinkA(Single * obj)
{
    fList.pop(obj);
    aList.push_front(obj);
}


void SingleSet::relinkD(Single * obj)
{
    aList.pop(obj);
    fList.push_front(obj);
}


void SingleSet::link(Object * obj)
{
    assert_true( obj->tag()==Single::TAG || obj->tag()==Single::WRIST_TAG );
    assert_true( obj->objset() == this );
    
    Single * S = static_cast<Single*>(obj);
    
    if ( S->attached() )
        aList.push_front(obj);
    else
        fList.push_front(obj);
    
    //std::clog << "SingleSet has " << fList.size() << "  " << aList.size() << '\n';
}

/**
 This will also detach the Hand
 */
void SingleSet::unlink(Object * obj)
{
    Single * s = static_cast<Single*>(obj);
    if ( s->attached() ) s->detach();
    fList.pop(obj);
}

//------------------------------------------------------------------------------
#pragma mark -


void SingleSet::foldPositions(Modulo const* m) const
{
    //std::cerr << "SingleSet::foldPositions()\n";
    Single * obj;
    for ( obj=firstF(); obj; obj=obj->next() ) obj->foldPosition(m);
    for ( obj=firstA(); obj; obj=obj->next() ) obj->foldPosition(m);
}


void SingleSet::shuffle()
{
    if ( aList.size() > 1 ) aList.shuffle();
    if ( fList.size() > 1 ) fList.shuffle();
}


void SingleSet::erase()
{
    for ( Property const* i : simul_.properties.find_all("single") )
    {
        SingleProp const * P = static_cast<SingleProp const*>(i);
        P->stocks.erase();
        P->uni_counts = 0;
    }

    ObjectSet::erasePool(fList);
    ObjectSet::erasePool(aList);
    inventory_.clear();
}


void SingleSet::detachAll()
{
    for ( Single * S=firstA(); S; S=S->next() )
        S->hand()->detachHand();
    fList.grab(aList);
    assert_true(aList.empty());
}


void SingleSet::defrostStore()
{
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        inventory_.unassign(i);
        Single * S = static_cast<Single*>(i);
        if ( S->hand()->attached() )
            S->hand()->detachHand();
        S->prop->stocks.push(S);
    }
    //infoStocks(std::clog);
}

//------------------------------------------------------------------------------
#pragma mark -


void SingleSet::freeze()
{
    //std::clog << "Single::freeze " << fList.size() << " " << aList.size() << "\n";
    assert_true(ice_.empty());
    ice_.grab(aList);
    ice_.grab(fList);
}


/** cnt[i] is the number of Singles of type `i` to be released */
void SingleSet::reheat(size_t cnt[], PropertyID n_cnt)
{
    //std::clog << "Single::reheat " << ice_.size() << " " << cnt[1] << " " << cnt[2] << "\n";
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        Single* S = static_cast<Single*>(i);
        // we want to skip the 'beforeDetachment' here:
        if ( S->attached() )
            S->hand()->detachHand();
        PropertyID id = S->prop->number();
        if ( id < n_cnt && 0 < cnt[id] )
        {
            --cnt[id];
            if ( S->prop->fast_diffusion > 0 )
                S->randomizePosition();
            addFreeSingle(S);
        }
        else
        {
            inventory_.unassign(S);
            S->prop->stocks.push(S);
        }
    }
}


void SingleSet::reheat()
{
    size_t sup = inventory_.capacity();
    size_t cnt[16] = { 0 };
    for ( int i = 0; i < 16; ++i ) cnt[i] = sup;
    reheat(cnt, 16);
}


/**
 This will save some of the objects in the normal way, using `write()`
 and will also save a count of the objects that were skipped
 */
void SingleSet::writeSomeObjects(Outputter& out) const
{
    writeRecords(out, fList.size(), inventory_.highest());
    
    std::map<PropertyID, size_t> cnt;
    // count all the elements that are virtually present:
    for ( SingleProp const* P : uniSingles )
        cnt[P->number()] = P->uni_counts;
    // write Single if `save_unbound > 0`:
    for ( Single const* n=firstF(); n; n=n->next() )
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
            // write counts for each class of unwritten Single:
            out.write("\n#section single reheat");
            for ( PropertyID i = 0; i <= sup; ++i )
                out.writeUInt(cnt[i], ' ');
        }
    }
    // decrement `save_unbound`:
    for ( Property * i : simul_.properties.find_all("single") )
    {
        SingleProp * P = static_cast<SingleProp *>(i);
        //printf("%s %u %lu\n", P->name().c_str(), P->save_unbound, cnt[P->number()]);
        P->save_unbound -= ( P->save_unbound > 0 );
    }
}


void SingleSet::writeSet(Outputter& out) const
{
    if ( sizeA() > 0 )
    {
        out.write("\n#section single A");
        writePool(out, aList);
    }
    if ( sizeF() > 0 )
    {
        out.write("\n#section single F");
        writeSomeObjects(out);
        //writePool(out, fList);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


void SingleSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << '\n' << title();
        PropertyList plist = simul_.properties.find_all(title());
        for ( Property const* i : plist )
        {
            SingleProp const* P = static_cast<SingleProp const*>(i);
            size_t cnt = count(match_property, P);
            os << '\n' << std::setw(10) << cnt << ' ' << P->name();
            os << " ( " << P->hand << " )";
        }
        if ( plist.size() > 1 )
            os << '\n' << std::setw(10) << size() << " total";
    }
}


ObjectList SingleSet::collect() const
{
    ObjectList res = ObjectSet::collect(fList);
    res.append( ObjectSet::collect(aList) );
    return res;
}


ObjectList SingleSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    ObjectList res = ObjectSet::collect(fList, func, arg);
    res.append( ObjectSet::collect(aList, func, arg) );
    return res;
}


size_t SingleSet::count(bool (*func)(Object const*, void const*), void const* arg) const
{
    size_t f = fList.count(func, arg);
    size_t a = aList.count(func, arg);
    return f + a;
}


int SingleSet::bad() const
{
    int err = 0;
    Single * obj;
    size_t cnt = sizeF();
    for ( obj = firstF(); obj ; obj=obj->next() )
    {
        if ( obj->attached() )
            err |= 8;
        if ( cnt-- == 0 )
            return 1;
    }
    if ( cnt ) return 1;
    
    cnt = sizeA();
    for ( obj = firstA();  obj ; obj=obj->next() )
    {
        if ( !obj->attached() )
            err |= 16;
        if ( simul_.fibers.badIdentity(obj->fiber()) )
            err |= 64;
        if ( obj->hand()->bad() )
            err |= 32;
        if ( cnt-- == 0 )
            return 2;
    }
    if ( cnt ) return 2;
    return err;
}

//------------------------------------------------------------------------------
#pragma mark - Wrists


/** Create Wrists anchored, and distributed to beads `name` */
void SingleSet::distributeWrists(ObjectList& objs, SingleProp const* sp,
                                 size_t cnt, std::string const& name) const
{
    BeadProp * bip = simul_.findProperty<BeadProp>("bead", name);
    if ( !bip )
        throw InvalidParameter("could not find bead type `"+name+"'");
    ObjectList list = simul_.beads.collect(bip);
    if ( list.empty() )
        throw InvalidParameter("could not find any bead of type `"+name+"'");
    list.shuffle_truncate(cnt);
    // create one Single on each of 'cnt' Beads:
    for ( Object const* i : list )
        objs.push_back(sp->newWrist(static_cast<Bead const*>(i), 0));
}


/**
 This will create Wrists with `mec` as Base, following the specifications given in `arg`.
 These Wrists will be anchored on points `[fip, fip+nbp[` of `mec`.
 
 The syntax understood for `arg` is as follows:

     [INTEGER] NAME_OF_SINGLE [each]

 The first optional integer specifies the number of Singles to be attached.
 Then follows the name of the Single to be created.
 If 'each' is specified, this number is multiplied by the number of point `nbp`,
 and every point receives the same number of Singles.
 
 This is used to attach Single to Bead, Solid and Sphere
 */
void SingleSet::makeWrists(ObjectList& objs, Mecable const* mec, index_t fip, index_t nbp, std::string const& arg)
{
    size_t num = 1;
    std::istringstream iss(arg);
    iss >> num;
    if ( iss.fail() )
    {
        num = 1;
        iss.clear();
    }
    if ( num == 0 || nbp == 0 )
        return;
    
    std::string str, mod;
    iss >> str >> mod;
    
    SingleProp * sip = simul_.findProperty<SingleProp>("single", str);

    if ( mod == "each" )
    {
        for ( index_t u = 0; u < num; ++u )
        {
            for ( index_t i = 0; i < nbp; ++i )
                objs.push_back(sip->newWrist(mec, fip+i));
        }
    }
    else
    {
        for ( index_t u = 0; u < num; ++u )
            objs.push_back(sip->newWrist(mec, fip+RNG.pint32(nbp)));
    }
}


SingleList SingleSet::collectWrists(Object const* arg) const
{
    SingleList res;
    
    for ( Single * s=firstF(); s; s=s->next() )
        if ( s->base() == arg )
            res.push_back(s);
    
    for ( Single * s=firstA(); s; s=s->next() )
        if ( s->base() == arg )
            res.push_back(s);
    
    return res;
}


void SingleSet::detachWrists(Object const* arg)
{
    Single * obj = firstA();
    while ( obj )
    {
        Single * nxt = obj->next();
        if ( obj->base() == arg )
            obj->detach();
        obj = nxt;
    }
}


void SingleSet::deleteWrists(Object const* arg)
{
    Single *nxt, *obj;
    
    obj = firstF();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->base() == arg )
            eraseObject(obj);
        obj = nxt;
    }

    obj = firstA();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->base() == arg )
            eraseObject(obj);
        obj = nxt;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion


size_t SingleSet::countStocks() const
{
    size_t res = 0;
    for ( Property const* i : simul_.properties.find_all("single") )
        res += static_cast<SingleProp const*>(i)->stocks.size();
    return res;
}


void SingleSet::infoStocks(std::ostream& os) const
{
    os << "  Single:stocks";
    for ( Property const* i : simul_.properties.find_all("single") )
    {
        SingleProp const * P = static_cast<SingleProp const*>(i);
        size_t cnt = count(match_property, P);
        os << " " << P->number() << ": " << cnt << " ( " << P->stocks.size() << " " << P->uni_counts << " )";
    }
    os << "\n";
}


void SingleSet::uniRefill(SingleProp const* sip, size_t cnt)
{
    SingleStock & can = sip->stocks;
    for ( size_t i = can.size(); i < cnt; ++i )
        can.push(sip->newSingle());
}


/**
 Attach exactly one Single from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void SingleSet::uniAttach(FiberSiteList& loc, SingleStock& can)
{
    // crop list to match available number of candidates:
    loc.shuffle_truncate(can.size());

    for ( FiberSite & i : loc )
    {
        Single * S = can.head();
        Hand const* h = S->hand();
        
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            i.reinterpolate();
            Vector pos = i.pos();

            Space const* spc = S->confineSpace();
            if ( spc )
            {
                // Only attach if position is near the edge of the Space:
                Vector prj = spc->project(pos);
                if ( distanceSqr(pos, prj) > square(h->property()->binding_range) )
                    continue;
                // Single will be placed on the edge of the Space:
                pos = prj;
            }
            else
            {
#if ( DIM > 1 )
                /*
                 Place the Single in the line perpendicular to the attachment point,
                 at a random distance within the range of attachment of the Hand.
                 This simulates a uniform spatial distribution of Single.
                 
                 This is important for certain Single such as Picket and PicketLong,
                 since they act as link between the anchoring position and the Hand.
                 */
                pos += i.dirFiber().randOrthoB(h->property()->binding_range);
#endif
            }

            can.pop();
            addFreeSingle(S);
            S->setPosition(pos);
            S->attach(i);
        }
    }
}


/**
 Implements a Monte-Carlo approach for attachments of free Single, assumming that
 diffusion is sufficiently fast to maintain a uniform spatial distribution,
 and that the distribution of fibers is more-or-less uniform such that the
 attachments are distributed randomly along the fibers.
 
 Diffusing (free) Single are removed from the standard list, and thus the
 random walk that is used for simulating diffusion will be skipped,
 as well as the detection of neighboring fibers done for attachments.
 
 Algorithm:
 - Remove diffusing Single from the simulation, transfering them to a 'reserve'.
 - Estimate the distance between binding sites occuring in one time-step, from:
    - the total length of fibers,
    - the volume of the Space,
    - the binding parameters of the relevant Hand.
    .
 - Attach Singles from the reserve, at random positions along the Fibers
 .
 
 Note: there is a similar feature for Couple
 */
void SingleSet::uniAttach(FiberSet const& fibers)
{
    FiberSiteList loc(128);
    
    // uniform attachment for selected classes:
    for ( SingleProp const* P : uniSingles )
    {
        SingleStock& can = P->stocks;
        
        // assuming (or not) a fixed number of diffusing molecules
        bool fixed = ( P->fast_reservoir > 0 );
        size_t cnt = ( fixed ? P->fast_reservoir : P->uni_counts );

        if ( cnt > 0 )
        {
            const real vol = P->spaceVolume();
            if ( P->fast_diffusion == 2 )
            {
                real dis = vol / ( cnt * P->hand_prop->bindingSectionRate() );
                fibers.newFiberSitesP(loc, dis);
            }
            else
            {
                real dis = vol / ( cnt * P->hand_prop->bindingSectionProb() );
                fibers.uniFiberSites(loc, dis);
            }
            
            if ( fixed ) // create enough candidates for all sites
                uniRefill(P, loc.size());

            size_t total = size();
            uniAttach(loc, can);
            if ( !fixed )
                P->uni_counts -= size() - total;
        }
    }
}


/**
 Sets the list `uniSingles` with single class concerned with `fast_diffusion`
*/
void SingleSet::uniPrepare(PropertyList const& properties)
{
    uniSingles.clear();
    for ( Property const* i : properties.find_all("single") )
    {
        SingleProp const* P = static_cast<SingleProp const*>(i);
        //std::clog << i->name() << "  " << P->fast_diffusion << "\n";
        if ( P->fast_diffusion > 0 )
            uniSingles.push_back(P);
    }
}


/**
 Revive all singles from the reserves, setting them in the free state
*/
void SingleSet::uniRelax()
{
    for ( SingleProp const* P : uniSingles )
    {
        makeSingles(P, P->uni_counts);
        P->uni_counts = 0;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Equilibration

/**
 Attach all unbound Single to nearby Fibers, to approximate an equilibrated state,
 as defined by the ratio of binding to unbinding rates of the Hand
 */
void SingleSet::equilibrate(FiberSet const& fibers, SingleProp const* sip)
{
    Space const* spc = sip->confine_space;
    if ( !spc )
        throw InvalidSyntax(sip->name()+":space must be defined first!");
    HandProp const* HP = sip->hand_prop;
    const real rge = HP->binding_range;
    const real sup = square(rge);
    // probability of finding hand in the bound state at equilibrium:
    const real prob = HP->binding_rate / ( HP->binding_rate + HP->unbinding_rate );

    FiberGrid grid;
    real grid_step = 1;
    Simul::setFiberGrid(grid, spc, grid_step);
    grid.paintGrid(fibers.first(), nullptr, rge);

    // equilibrate unattached Singles:
    Single *nxt, *obj = firstF();
    
    while ( obj )
    {
        nxt = obj->next();
        if ( !sip || obj->property() == sip )
        {
            Vector pos = obj->position();
            Hand * ha = obj->hand();
            
            // check all segments within grid cell
            for ( FiberSegment const& seg : grid.nearbySegments(pos) )
            {
                real dis = INFINITY;
                real abs = seg.projectPoint(pos, dis);
                if ( dis < sup && RNG.test(prob) )
                {
                    // ATTENTION: convert `abs` relative to the segment to Fiber's abscissa
                    FiberSite sit(seg.fiber(), seg.abscissa1()+abs);
                    if ( ha->attachmentAllowed(sit) )
                    {
                        ha->attach(sit);
                        //std::clog << "   bind " << sit << " at " << 1000*std::sqrt(dis) << " nm\n";
                        break;
                    }
                }
            }
        }
        obj = nxt;
    }
}


void SingleSet::equilibrate(SingleProp const* P)
{
    equilibrate(simul_.fibers, P);
}


void SingleSet::equilibrate()
{
    for ( Property const* i : simul_.properties.find_all("single") )
    {
        SingleProp const* P = static_cast<SingleProp const*>(i);
        equilibrate(simul_.fibers, P);
    }
}

