// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "single_set.h"
#include "single_prop.h"
#include "glossary.h"
#include "iowrapper.h"

#include "simul.h"
#include "simul_prop.h"
#include "property_list.h"
#include "wrist.h"
#include "wrist_long.h"



//------------------------------------------------------------------------------

void SingleSet::prepare(PropertyList const& properties)
{
    uniEnabled = uniPrepare(properties);
}


/// templated member function pointer...
/**
In the loops we get the 'next' in the list always before calling 'FUNC', since
'step()' may transfer the node to another list, changing the value of 'next()'
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
SingleSet::step() must call the appropriate Single::step() exactly once
for each Single: either stepF() or stepAA().

The Singles are stored in two lists, and are automatically transferred
from one list to the other if their Hands bind or unbind. The code relies
on the fact that a Single will be moved to the start of the list to which it
is transferred, by 'push_front'. By starting always from the node that
was first before any transfer could occur, we process each Couple only once.
*/
void SingleSet::step()
{
    //Cytosim::log("SingleSet::step entry : F %5i A %5i\n", fList.size(), aList.size());
    
    Single *const fHead = firstF();
    bool fOdd = sizeF() & 1;
    
    step_singles<&Single::stepA>(firstA(), sizeA() & 1);
    
    // use alternative attachment strategy:
    if ( uniEnabled )
    {
        Single * obj = uniCollect(fHead);
        uniAttach(simul_.fibers);
        // handle Single for which 'fast_diffusion = false':
        {
            while ( obj )
            {
                Single * nxt = obj->next();
                obj->stepF();
                obj = nxt;
            }
        }
    }
    else
    {
        step_singles<&Single::stepF>(fHead, fOdd);
    }
}


/**
 This version does not simulate the attachment of free Hand, and hence skips
 Single::stepF() that calls Hand::stepUnattached().

 This is only used if POOL_UNATTACHED > 1
*/
void SingleSet::stepSkipUnattached()
{
    //Cytosim::log("SingleSet::stepSkipUnattached : F %5i A %5i\n", fList.size(), aList.size());
    
    Single *const fHead = firstF();
    
    step_singles<&Single::stepA>(firstA(), sizeA() & 1);
    
    // use alternative attachment strategy:
    if ( uniEnabled )
    {
        uniCollect(fHead);
        uniAttach(simul_.fibers);
    }
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


Object * SingleSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Single::TAG )
    {
        SingleProp * p = simul_.findProperty<SingleProp>("single", pid);
        return p->newSingle();
    }
    else if ( tag == Single::TAG_WRIST )
    {
        SingleProp * p = simul_.findProperty<SingleProp>("single", pid);
        return p->newWrist(nullptr, 0);
    }
    throw InvalidIO("Warning: unknown Single tag `"+std::to_string(tag)+"'");
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
void SingleSet::newObjects(ObjectList& res, const std::string& name, Glossary& opt)
{
    SingleProp * p = simul_.findProperty<SingleProp>("single", name);
    
    Single * obj = nullptr;
    std::string str;
    if ( opt.set(str, "base") )
    {
        Mecable * mec = simul_.findMecable(str);
        if ( !mec )
            throw InvalidParameter("could not find Mecable specified in single:base `"+str+"'");
        // get index of point in second argument
        size_t ip = 0;
        if ( opt.set(str, "base", 1) )
            ip = mec->point_index(str);
         
        obj = p->newWrist(mec, ip);
    }
    else
        obj = p->newSingle();

    // Allow user to attach Hand to an existing fiber
    if ( opt.has_key("attach") )
        obj->attach(simul_.fibers.someSite("attach", opt));
    
    // Allow user to attach Hand to an existing fiber
    if ( opt.has_key("site") )
        obj->attach(simul_.fibers.someSite("site", opt));

    res.push_back(obj);
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
    assert_true( obj->tag()==Single::TAG || obj->tag()==Single::TAG_WRIST );
    assert_true( obj->objset() == this );

    if ( static_cast<Single*>(obj)->attached() )
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
    assert_true(obj->objset() == nullptr);

    if ( s->attached() )
    {
        aList.pop(obj);
        s->detach();
    }
    else
    {
        fList.pop(obj);
    }
}


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
    relax();
    ObjectSet::erasePool(fList);
    ObjectSet::erasePool(aList);
    inventory_.clear();
}


void SingleSet::freeze()
{
    relax();
    assert_true(ice_.empty());
    ice_.append(aList);
    ice_.append(fList);
}


void SingleSet::pruneDetach()
{
    /* After reading from file, the Hands should not update
     any Fiber, Single or Couple as they will be deleted */
    Object * i = ice_.pop_front();
    while ( i )
    {
        Single* o = static_cast<Single*>(i);
        i = ice_.pop_front();
        o->hand()->detachHand();
        o->randomizePosition();
        linkF(o);
    }
}

/** After reading from file, the Hands should not update
 any Fiber, Single or Couple as they will be deleted */
void SingleSet::pruneDelete()
{
    Object * i = ice_.pop_front();
    while ( i )
    {
        Single* o = static_cast<Single*>(i);
        i = ice_.pop_front();
        inventory_.unassign(o);
        o->objset(nullptr);
        delete(o);
    }
}

void SingleSet::prune()
{
    if ( prune_mode )
        pruneDetach();
    else
        pruneDelete();
}


void SingleSet::write(Outputter& out) const
{
    if ( sizeA() > 0 )
    {
        out.writeLine("\n#section single A");
        writeObjects(out, aList);
    }
    if ( sizeF() > 0 )
    {
        int skip = simul_.prop.skip_free_single || prune_mode;
        if ( skip & skip_now )
            out.writeLine("\n#section single F 1");
        else
        {
            out.writeLine("\n#section single F 0");
            writeObjects(out, fList);
            skip_now = skip;
        }
    }
}


int SingleSet::bad() const
{
    int code = fList.bad() | aList.bad();

    if ( code )
        return code;
    
    Single * obj;
    for ( obj = firstF(); obj ; obj=obj->next() )
    {
        if ( obj->attached() )
            code |= 8;
    }
    
    for ( obj = firstA();  obj ; obj=obj->next() )
    {
        if ( !obj->attached() )
            code |= 16;
        if ( obj->hand()->bad() )
            code |= 32;
    }
    return code;
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
            SingleProp const* p = static_cast<SingleProp const*>(i);
            size_t cnt = count(match_property, p);
            os << '\n' << std::setw(10) << cnt << ' ' << p->name();
            os << " ( " << p->hand << " )";
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
    size_t f = ObjectSet::count(fList, func, arg);
    size_t a = ObjectSet::count(aList, func, arg);
    return f + a;
}

//------------------------------------------------------------------------------
#pragma mark - Wrists

/**
 This will create Wrists with `obj` as Base, following the specifications given in `arg`.
 These Wrists will be anchored on points `fip` to `fip+nbp-1` of `obj`.
 
 The syntax understood for `arg` is as follows:

     [INTEGER] NAME_OF_SINGLE [each]

 The first optional integer specifies the number of Singles to be attached.
 Then follows the name of the Single to be created.
 If 'each' is specified, this number is multiplied by the number of point `nbp`,
 and every point receives the same number of Singles.
 
 This is used to attach Single to Bead, Solid and Sphere
 */
ObjectList SingleSet::makeWrists(Mecable const* obj, size_t fip, size_t nbp, std::string& arg)
{
    ObjectList res;
    size_t num = 1;

    std::istringstream iss(arg);
    iss >> num;
    
    if ( iss.fail() )
    {
        num = 1;
        iss.clear();
    }
    
    if ( num == 0 || nbp == 0 )
        return res;
    
    std::string str, mod;
    iss >> str >> mod;
    
    SingleProp * sip = simul_.findProperty<SingleProp>("single", str);

    if ( mod == "each" )
    {
        for ( size_t u = 0; u < num; ++u )
        {
            for ( size_t i = 0; i < nbp; ++i )
                res.push_back(sip->newWrist(obj, fip+i));
        }
    }
    else
    {
        for ( size_t u = 0; u < num; ++u )
        {
            res.push_back(sip->newWrist(obj, fip+RNG.pint32(nbp)));
        }
    }
    
    return res;
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


void SingleSet::removeWrists(Object const* arg)
{
    Single *nxt, *obj;
    
    obj = firstF();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->base() == arg )
            delete(obj);
        obj = nxt;
    }

    obj = firstA();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->base() == arg )
            delete(obj);
        obj = nxt;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion


void SingleSet::uniRefill(SingleList& can, size_t cnt, SingleProp const* p)
{
    for ( size_t i = can.size(); i < cnt; ++i )
    {
        Single* s = p->newSingle();
        inventory_.assign(s);
        s->objset(this);
        can.push_back(s);
    }
}


/**
 Attach exactly one Single from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void SingleSet::uniAttach(Array<FiberSite>& loc, SingleList& can)
{
    // crop list to match available number of candidates:
    if ( can.size() < loc.size() )
    {
        loc.shuffle();
        loc.truncate(can.size());
    }
    
    for ( FiberSite & i : loc )
    {
        Single * s = can.back();
        Hand const* h = s->hand();
        
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            Vector pos = i.pos();
            Space const* spc = i.fiber()->prop->confine_space_ptr;

            // Only attach if position is within the confining Space:
            if ( spc && spc->outside(pos) )
                continue;

            if ( s->prop->fast_diffusion & 8 )
            {
                if ( ! spc )
                    continue;
                // Only attach if position is near the edge of the Space:
                Vector prj = spc->project(pos);
                if ( distanceSqr(pos, prj) >= square(h->prop->binding_range) )
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
                pos += i.dirFiber().randOrthoB(h->prop->binding_range);
#endif
            }

            can.pop_back();
            linkF(s);
            s->setPosition(pos);
            s->attach(i);
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
    Array<FiberSite> loc(128, 128);
    
    // uniform attachment for the reserves:
    for ( SingleReserve & reserve : uniReserves )
    {
        SingleProp const* p = reserve.first;
        if ( !p )
            continue;
        
        const real vol = p->spaceVolume();
        SingleList& can = reserve.second;
        
        // assuming (or not) a fixed number of diffusing molecules
        bool fixed = ( p->fast_diffusion_nb > 0 );
        size_t cnt = ( fixed ? p->fast_diffusion_nb : can.size());

        if ( cnt > 0 )
        {
            if ( p->fast_diffusion & 2 )
            {
                real dis = vol / ( cnt * p->hand_prop->bindingSectionRate() );
                fibers.newFiberSitesP(loc, dis);
            }
            else
            {
                real dis = vol / ( cnt * p->hand_prop->bindingSectionProb() );
                fibers.uniFiberSites(loc, dis);
            }
            
            if ( fixed ) // create enough candidates for all sites
                uniRefill(can, loc.size(), p);

            uniAttach(loc, can);
        }
    }
}


/**
 
 Return true if at least one single:fast_diffusion is true,
 and in this case allocate uniReserves.
 
 The Volume of the Space is assumed to remain constant until the next uniPrepare()
 */
bool SingleSet::uniPrepare(PropertyList const& properties)
{
    bool res = false;
    PropertyID last = 0;

    PropertyList allprop = properties.find_all("single");
    
    for ( Property const* i : allprop )
    {
        SingleProp const* p = static_cast<SingleProp const*>(i);
        last = std::max(last, p->number());
        assert_true(p->number() > 0);
        res |= p->fast_diffusion;
    }
    
    if ( res )
    {
        uniReserves.resize(last);
        for ( SingleReserve & reserve : uniReserves )
            reserve.first = nullptr;
        for ( Property const* i : allprop )
        {
            SingleProp const* p = static_cast<SingleProp const*>(i);
            if ( p->fast_diffusion )
                uniReserves[p->number()-1].first = p;
        }
    }
    
    return res;
}


/**
 Transfer free Single with `fast_diffusion` to the reserves, starting from `obj`.
 Return first Single that was not transferred
*/
Single* SingleSet::uniCollect(Single * obj)
{
    Single * res = nullptr;
    Single * nxt;
    while ( obj )
    {
        nxt = obj->next();
        SingleProp const* p = obj->prop;
        if ( p->fast_diffusion )
        {
            fList.pop(obj);
            uniReserves[p->number()-1].second.push_back(obj);
        }
        else if ( !res )
            res = obj;
        obj = nxt;
    }
    return res;
}


/**
Release all Singles from the reserves
*/
void SingleSet::uniRelax()
{
    for ( SingleReserve & reserve : uniReserves )
    {
        for ( Single * s : reserve.second )
        {
            assert_true(!s->attached());
            s->randomizePosition();
            linkF(s);
        }
        reserve.second.clear();
    }
}

