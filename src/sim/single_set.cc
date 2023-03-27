// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
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
 Transfer free Single with `fast_diffusion` to the reserves, starting from `obj`.
 Return first Single that was not transferred
*/
void SingleSet::uniStepCollect(Single * obj)
{
    Single * nxt;
    while ( obj )
    {
        nxt = obj->next();
        SingleProp const* P = obj->prop;
        if ( P->fast_diffusion && !obj->base() )
        {
            fList.pop(obj);
            uniReserves[P->number()].push(obj);
        }
        else
            obj->stepF();
        obj = nxt;
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
    ObjectID h = inventory_.highest();
    if ( h > 4096 && h > 2 * ( size() + uniCounts() ) )
    {
        inventory_.reassign();
        std::clog << "Single::reassign(" << h << " ---> " << inventory_.highest() << ")\n";
    }
    //Cytosim::log("SingleSet::step entry : F %5i A %5i\n", fList.size(), aList.size());
    
    Single *const fHead = firstF();
    bool fOdd = sizeF() & 1;
    
    step_singles<&Single::stepA>(firstA(), sizeA() & 1);
    
    // use alternative attachment strategy:
    if ( uniEnabled )
    {
        uniStepCollect(fHead);
        uniAttach(simul_.fibers);
    }
    else
    {
        step_singles<&Single::stepF>(fHead, fOdd);
    }
}


/**
 This version does not simulate the attachment of free Hand, and hence skips
 Single::stepF() that performs attachment

 This is only used if POOL_UNATTACHED > 1
*/
void SingleSet::stepSkipUnattached()
{
    //Cytosim::log("SingleSet::stepSkipUnattached : F %5i A %5i\n", fList.size(), aList.size());
    
    step_singles<&Single::stepA>(firstA(), sizeA() & 1);
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
        SingleProp * P = simul_.findProperty<SingleProp>("single", pid);
        return P->newSingle();
    }
    else if ( tag == Single::TAG_WRIST )
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
        size_t ip = 0;
        if ( opt.set(str, "base", 1) )
            ip = mec->point_index(str);
         
        obj = pp->newWrist(mec, ip);
    }
    else
        obj = pp->newSingle();

    // Allow user to attach Hand to an existing fiber
    if ( opt.has_key("attach") )
        obj->attach(simul_.fibers.someSite("attach", opt));
    
    // Allow user to attach Hand to an existing fiber
    if ( opt.has_key("site") )
        obj->attach(simul_.fibers.someSite("site", opt));

    return ObjectList(obj);
}


ObjectList SingleSet::distributeWrists(SingleProp const* sp, size_t cnt,
                                       std::string const& name) const
{
    ObjectList objs;
    BeadProp * bip = simul_.findProperty<BeadProp>("bead", name);
    if ( !bip )
        throw InvalidParameter("could not find Bead type `"+name+"'");
    ObjectList list = simul_.beads.collect(bip);
    if ( list.empty() )
        throw InvalidParameter("could not find Bead of type `"+name+"'");
    if ( cnt < list.size() )
    {
        list.shuffle();
        list.truncate(cnt);
    }
    // create one Single on each of 'cnt' Beads:
    for ( Object const* i : list )
        objs.push_back(sp->newWrist(static_cast<Bead const*>(i), 0));
    return objs;
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
    for ( SingleReserve & can : uniReserves )
        can.erase();
    ObjectSet::erasePool(fList);
    ObjectSet::erasePool(aList);
    inventory_.clear();
}


void SingleSet::freeze()
{
    relax();
#ifdef MORE_ROBUST_READING
    flag(aList, 7);
    flag(fList, 7);
#endif
    assert_true(ice_.empty());
    ice_.grab(aList);
    ice_.grab(fList);
}


void SingleSet::reheat()
{
    //std::clog << "Single::reheat " << ice_.size() << "\n";
    Object * i = ice_.front();
    while ( i )
    {
        Single* o = static_cast<Single*>(i);
        i = i->next();
        if ( o->prop->fast_diffusion )
        {
            ice_.pop(o);
            // we want to skip the 'beforeDetachment' here:
            o->hand()->detachHand();
            o->randomizePosition();
            linkF(o);
        }
    }
}


void SingleSet::writeF_skip(Outputter& out) const
{
    // record number of objects:
    size_t cnt = fList.size();
    size_t sup = inventory_.highest();
    out.write("\n#record "+std::to_string(cnt)+" "+std::to_string(sup));
    if ( out.binary() ) out.put_char('\n');
    
    for ( Single const* n=firstF(); n; n=n->next() )
    {
        if ( ! n->prop->fast_diffusion )
            n->write(out);
    }
    out.write("\n#section single reheat");
}


void SingleSet::writeSet(Outputter& out) const
{
    if ( sizeA() > 0 )
    {
        out.write("\n#section single A");
        writeObjects(out, aList);
    }
    if ( sizeF() > 0 )
    {
        out.write("\n#section single F");
        if ( simul_.prop.skip_free_single && skip_now )
            writeF_skip(out);
        else
        {
            writeObjects(out, fList);
            skip_now = 1;
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
void SingleSet::makeWrists(ObjectList& objs, Mecable const* mec, size_t fip, size_t nbp, std::string& arg)
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
        for ( size_t u = 0; u < num; ++u )
        {
            for ( size_t i = 0; i < nbp; ++i )
                objs.push_back(sip->newWrist(mec, fip+i));
        }
    }
    else
    {
        for ( size_t u = 0; u < num; ++u )
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


void SingleSet::deleteInvalidWrists()
{
    Single *nxt, *obj;

    obj = firstF();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->invalid() )
        {
            std::clog << "WARNING: invalid base for " << obj->reference() << "\n";
            delete(obj);
        }
        obj = nxt;
    }

    obj = firstA();
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->invalid() )
        {
            std::clog << "WARNING: invalid base for " << obj->reference() << "\n";
            delete(obj);
        }
        obj = nxt;
    }
}

//------------------------------------------------------------------------------
#pragma mark - Fast Diffusion


size_t SingleSet::uniCounts() const
{
    size_t res = 0;
    for ( SingleReserve const& can : uniReserves )
        res += can.size();
    return res;
}


void SingleSet::uniRefill(SingleReserve& can, size_t cnt, SingleProp const* sip)
{
    for ( size_t i = can.size(); i < cnt; ++i )
    {
        Single* S = sip->newSingle();
        inventory_.assign(S);
        S->objset(this);
        can.push(S);
    }
}


/**
 Attach exactly one Single from `can` to each site in `loc`.
 If `can` is not large enough, a subset of `loc` is selected.
 */
void SingleSet::uniAttach(Array<FiberSite>& loc, SingleReserve& can)
{
    // crop list to match available number of candidates:
    if ( can.size() < loc.size() )
    {
        loc.shuffle();
        loc.truncate(can.size());
    }
    
    for ( FiberSite & i : loc )
    {
        Single * s = can.head();
        Hand const* h = s->hand();
        
        if ( h->keyMatch(i.fiber()) &&  h->attachmentAllowed(i) )
        {
            i.reinterpolate();
            Vector pos = i.pos();

            if ( s->prop->confine == CONFINE_ON )
            {
                // Only attach if position is near the edge of the Space:
                Vector prj = s->confineSpace()->project(pos);
                if ( distanceSqr(pos, prj) >= square(h->property()->binding_range) )
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
    for ( SingleReserve & can : uniReserves )
    {
        SingleProp const* P = can.property();
        if ( !P )
            continue;
        
        const real vol = P->spaceVolume();
        
        // assuming (or not) a fixed number of diffusing molecules
        bool fixed = ( P->fast_reservoir > 0 );
        size_t cnt = ( fixed ? P->fast_reservoir : can.size());

        if ( cnt > 0 )
        {
            if ( P->fast_diffusion & 2 )
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
                uniRefill(can, loc.size(), P);

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
        SingleProp const* P = static_cast<SingleProp const*>(i);
        last = std::max(last, P->number());
        assert_true(P->number() > 0);
        res |= P->fast_diffusion;
    }
    
    if ( res )
    {
        uniReserves.resize(last+1);
        for ( SingleReserve & can : uniReserves )
            can.property(nullptr);
        for ( Property const* i : allprop )
        {
            SingleProp const* P = static_cast<SingleProp const*>(i);
            if ( P->fast_diffusion )
                uniReserves[P->number()].property(P);
        }
    }
    
    return res;
}


/**
Release all Singles from the reserves
*/
void SingleSet::uniRelax()
{
    for ( SingleReserve & can : uniReserves )
    {
        Single * o = can.head();
        while ( o )
        {
            assert_true(!o->attached());
            can.pop();
            o->randomizePosition();
            linkF(o);
            o = can.head();
        }
    }
}

