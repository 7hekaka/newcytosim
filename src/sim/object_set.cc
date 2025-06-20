// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University


#include "object_set.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "tokenizer.h"
#include "glossary.h"
#include "modulo.h"
#include "space.h"
#include "property_list.h"
#include "simul.h"
#include <errno.h>

//------------------------------------------------------------------------------

/**
 The object is added at the front of the list
 */
void ObjectSet::link(Object * obj)
{
    assert_true( obj->objset() == this );
    pool_.push_back(obj);
    
    //std::clog << "ObjectSet has " << pool_.size() << '\n';
}


void ObjectSet::unlink(Object * obj)
{
    //assert_true( obj->objset() == this );
    pool_.pop(obj);
}

//------------------------------------------------------------------------------
#pragma mark -


void ObjectSet::flagObjects(ObjectList const& objs, ObjectFlag f)
{
    for ( Object * obj : objs )
        obj->flag(f);
}


void ObjectSet::markObjects(ObjectList const& objs, ObjectMark k)
{
    for ( Object * obj : objs )
        obj->mark(k);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Translate all listed movable objects ( Object::mobile() & 1 ) by `vec`
 */
void ObjectSet::translateObjects(ObjectList const& objs, Vector const& vec)
{
    for ( Object * obj : objs )
        if ( obj->mobile() & 1 )
            obj->translate(vec);
}

/**
 Apply Rotation around the origin to all movable objects in list
 */
void ObjectSet::rotateObjects(ObjectList const& objs, Rotation const& rot)
{
    for ( Object * obj : objs )
        if ( obj->mobile() & 2 )
            obj->rotate(rot);
}

/**
 Apply isometry to all objects
 */
void ObjectSet::moveObjects(ObjectList const& objs, Isometry const& iso)
{
    //std::clog << "moving " << objs.size() << " objects" << '\n';
    for ( Object * obj : objs )
        obj->move(iso);
}


/**
 Translate movable objects in list if ( obj->flag() != f )
 */
void ObjectSet::translateObjects(ObjectList const& objs, Vector const& vec, ObjectFlag f)
{
    for ( Object * obj : objs )
    {
        if ( obj->mobile() & 1 && obj->flag() != f )
        {
            obj->translate(vec);
            obj->flag(f);
        }
    }
}

/**
 Apply Rotation around the origin to objects in list if ( obj->flag() != f )
 */
void ObjectSet::rotateObjects(ObjectList const& objs, Rotation const& rot, ObjectFlag f)
{
    for ( Object * obj : objs )
    {
        if ( obj->mobile() & 2 && obj->flag() != f )
        {
            obj->rotate(rot);
            obj->flag(f);
        }
    }
}

/**
Apply isometry to objects in list if ( obj->flag() != f )
 */
void ObjectSet::moveObjects(ObjectList const& objs, Isometry const& iso, ObjectFlag f)
{
    //std::clog << "moving " << objs.size() << " objects" << '\n';
    for ( Object * obj : objs )
    {
        if ( obj->flag() != f )
        {
            //std::clog << "    moving " << obj->reference() << '\n';
            obj->move(iso);
            obj->flag(f);
        }
        //else std::clog << "    already moved " << obj->reference() << '\n';
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void ObjectSet::add(Object * obj)
{
    //std::clog << "ObjectSet::add " << obj->reference() << '\n';
    if ( !obj->linked() )
    {
        assert_true( !obj->objset() || obj->objset() == this );
        inventory_.assign(obj);
        obj->objset(this);
        link(obj);
        //std::clog << "ObjectSet::add(" << obj->reference() << ") " << obj->identity() << "\n";
    }
    else
    {
        std::clog << "Warning: attempted to re-link "+obj->reference()+" \n";
    }
}


void ObjectSet::add(ObjectList const& list)
{
    for ( Object * obj : list )
        add(obj);
}


void ObjectSet::remove(Object * obj)
{
    //std::clog << "ObjectSet::remove " <<  obj->reference() << '\n';
    assert_true( obj->objset() == this );
    inventory_.unassign(obj);
    obj->objset(nullptr);
    unlink(obj);
}


void ObjectSet::remove(ObjectList const& list)
{
    for ( Object * obj : list )
        remove(obj);
}


void ObjectSet::eraseObject(Object * obj)
{
    //std::clog << "ObjectSet::erase " << obj->reference() << '\n';
    remove(obj);
    delete(obj);
}


void ObjectSet::erasePool(ObjectPool & list)
{
    Object * i;
    while (( i = list.front() ))
    {
        list.pop_front();
        //i->objset(nullptr);
        delete(i);
    }
}


void ObjectSet::erase()
{
    erasePool(pool_);
    inventory_.clear();
}


void ObjectSet::eraseObjects(ObjectList const& objs)
{
    //std::clog << " ObjectSet::erase(" << objs.size() << " objects):\n";
    for ( Object * obj : objs )
    {
        assert_true( obj );
        //std::clog << "   erase " << obj->reference() << '\n';
        assert_true( obj->objset() == this );
        inventory_.unassign(obj);
        unlink(obj);
        delete(obj);
    }
}


/// return string that will identify the object in `pickObject`
std::string ObjectSet::nameObject(Object const* obj) const
{
    if ( !obj )
        return "none";
    Property const* prp = obj->property();
    ObjectID end = obj->identity();
    ObjectID cnt = 1;

    for ( ObjectID id = inventory_.lowest(); id < end; ++id )
    {
        Inventoried * i = inventory_[id];
        if ( i )
            cnt += ( static_cast<Object*>(i)->property() == prp );
    }
    return prp->name() + std::to_string(cnt);
}


Object* ObjectSet::findObject(Property const* pp, long num) const
{
    Inventoried* inv = nullptr;
    if ( num > 0 )
    {
        // 'microtubule1' would return the first created microtubule
        // std::clog << "findObject -> highest pick `" << spec << num << "'\n";
        inv = inventory_.first();
        while ( inv )
        {
            num -= ( static_cast<Object*>(inv)->property() == pp );
            if ( num <= 0 )
                break;
            inv = inventory_.next(inv);
        }
    }
    else
    {
        // 'microtubule0' would return the last created microtubule
        //std::clog << "findObject -> highest pick `" << spec << num << "'\n";
        inv = inventory_.last();
        while ( inv )
        {
            num += ( static_cast<Object*>(inv)->property() == pp );
            if ( num >= 0 )
                break;
            inv = inventory_.previous(inv);
        }
    }
    return static_cast<Object*>(inv);
}


Object* ObjectSet::findObject(const std::string& cat, std::string spec, long num) const
{
    //std::clog << "findObject(" << cat << "|" << num << "|" << spec << ")\n";
    // check for a string starting with the class name (eg. 'fiber'):
    if ( spec == cat )
    {
        Inventoried * inv = nullptr;
        if ( num > 0 )
        {
            inv = inventory_.get(num);
        }
        else
        {
            // start from the end of the list:
            inv = inventory_.last();
            while ( inv  &&  ++num < 0 )
                inv = inventory_.previous(inv);
        }
        return static_cast<Object*>(inv);
    }
    
    // check if string starts with 'first'
    if ( spec == "first" )
    {
        Inventoried* inv = inventory_.first();
        while ( inv  &&  --num >= 0 )
            inv = inventory_.next(inv);
        return static_cast<Object*>(inv);
    }
    
    // check if string starts with 'last'
    if ( spec == "last" )
    {
        Inventoried* inv = inventory_.last();
        while ( inv  &&  ++num <= 0 )
            inv = inventory_.previous(inv);
        return static_cast<Object*>(inv);
    }
    
    if ( spec.empty() )
        return nullptr;

    // finally search for a property name:
    Property * pp = simul_.findProperty(cat, spec);
    if ( pp )
        return findObject(pp, num);
    
    return nullptr;
}


/*
 There are several ways to designate an object.
 For example, if the class name (title) is 'fiber', one may use:
 - `fiber1`  indicates fiber number 1
 - `fiber2`  indicates fiber number 2, etc.
 - `first`   indicates the oldest fiber remaining
 - `first+1` indicates the second oldest fiber remaining
 - `last`    indicates the last fiber created
 - `last-1`  indicates the penultimate fiber created
 - `fiber0`  the last fiber created,
 - `fiber-1` the penultimate fiber, etc.
 .
 */
Object* ObjectSet::pickObject(const std::string& cat, std::string spec) const
{
    //std::clog << "ObjectSet::findObject(" << cat << ", " << spec << ")\n";
    
    if ( spec == "first" )
        return static_cast<Object*>(inventory_.first());
    
    if ( spec == "last" )
        return static_cast<Object*>(inventory_.last());

    // try to split into a word and a number:
    long num = 0;
    if ( Tokenizer::split_polysymbol(spec, num) )
    {
        //std::clog << "pickObject(" << spec << " " << num << ")\n";
        return findObject(cat, spec, num);
    }
    
    // check category name, eg. 'fiber':
    if ( spec == cat )
    {
        ObjectList all = collect();
        //std::clog << "pickObject(" << sel.size() << " " << title << ")\n";
        if ( all.size() > 0 )
            return all.pick_one();
    }
    
    // check property name:
    Property const* P = simul_.findProperty(cat, spec);
    if ( P )
    {
        ObjectList sel = collect(match_property, P);
        //std::clog << "pickObject(" << sel.size() << " " << spec << ")\n";
        if ( sel.size() > 0 )
            return sel.pick_one();
    }

    return nullptr;
}


/**
 return the first object encountered with the given property,
 but it can be any one of them, since the lists are regularly
 shuffled to randomize the order in the list.
 */
Object * ObjectSet::pickObject(Property const* p) const
{
    for ( Object* obj=first(); obj; obj=obj->next() )
        if ( obj->property() == p )
            return obj;
    return nullptr;
}


ObjectList ObjectSet::collect(const ObjectPool & list)
{
    ObjectList res;
    for ( Object* n = list.front(); n; n=n->next() )
        res.push_back(n);
    return res;
}


ObjectList ObjectSet::collect(const ObjectPool & list,
                              bool (*func)(Object const*, void const*), void const* arg)
{
    ObjectList res(list.size());
    Object * n = list.front();
    while ( n )
    {
        if ( func(n, arg) )
            res.push_back(n);
        n = n->next();
    }
    return res;
}


ObjectList ObjectSet::collect() const
{
    return collect(pool_);
}


ObjectList ObjectSet::collect(const size_t cnt) const
{
    ObjectList objs = collect();
    objs.shuffle_truncate(cnt);
    return objs;
}


ObjectList ObjectSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    return collect(pool_, func, arg);
}


ObjectList ObjectSet::collect(bool (*func)(Object const*, void const*), void const* arg, const size_t cnt) const
{
    ObjectList objs = collect(func, arg);
    objs.shuffle_truncate(cnt);
    return objs;
}


size_t ObjectSet::count(bool (*func)(Object const*, void const*), void const* arg) const
{
    return pool_.count(func, arg);
}

//------------------------------------------------------------------------------
#pragma mark - I/O


void ObjectSet::flag(ObjectPool const& list, ObjectFlag f)
{
    for ( Object * n=list.front(); n; n=n->next() )
        n->flag(f);
}


void ObjectSet::freeze()
{
    assert_true(ice_.empty());
    ice_.grab(pool_);
}


void ObjectSet::defrost()
{
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        //std::clog << "delete " << i->reference() << "\n";
        inventory_.unassign(i);
        i->objset(nullptr);
        delete(i);
    }
}


void ObjectSet::thaw()
{
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        link(i);
    }
}

// record number of objects:
void ObjectSet::writeRecords(Outputter& out, size_t tot, size_t sup) const
{
    out.write("\n#record "+std::to_string(tot)+" "+std::to_string(sup));
    if ( out.binary() ) out.put_char('\n');
}

/**
 Write Reference and Object's data, for all Objects in `list`
 */
void ObjectSet::writePool(Outputter& out, ObjectPool const& list) const
{
    writeRecords(out, list.size(), inventory_.highest());
    // record objects:
    for ( Object const* n=list.front(); n; n=n->next() )
    {
        //std::clog << "write " << n->reference() << '\n';
        n->write(out);
    }
}


/**  read Object's mark in binary format. This should match Object::writeMarker() */
static void readMarker(Inputter& in, PropertyID& ix, ObjectID& id)
{
#if 1
    union { uint16_t u; uint8_t c[2]; } u16;
    ix = in.get_char();
    u16.c[0] = in.get_char();
    u16.c[1] = in.get_char();
    id = u16.u;
#else
    ix = in.readUInt8();
    id = in.readUInt16();
#endif
}


/**  read Object's mark in binary format. This should match Object::writeMarker() */
static void readMarkerFat(Inputter& in, PropertyID& ix, ObjectID& id, ObjectMark& mk)
{
    union { uint16_t u; uint8_t c[2]; } u16;
    union { uint32_t u; uint8_t c[4]; } u32;
    
#if BACKWARD_COMPATIBILITY < 58
    if ( in.formatID() < 58 ) // 26.11.2022
    {
        ix = in.readUInt16();
        id = in.readUInt32();
        mk = in.readUInt32();
    }
    else
#endif
    {
        mk = in.get_char();
        u16.c[0] = in.get_char();
        u16.c[1] = in.get_char();
        ix = u16.u;
        u32.c[0] = in.get_char();
        u32.c[1] = in.get_char();
        u32.c[2] = in.get_char();
        u32.c[3] = in.get_char();
        id = u32.u;
    }
}


/**  read Object's mark in text format. This should match Object::writeMarker() */
static void readMarkerASCII(Inputter& in, PropertyID& ix, ObjectID& id, ObjectMark& mk)
{
    unsigned u = 0;
    FILE * f = in.file();
    if ( 1 != fscanf(f, "%u", &u) )
        throw InvalidIO("invalid Object header");
    ix = (PropertyID)u;
    if ( in.get_char() != ':' )
        throw InvalidIO("invalid Object header");
    if ( 1 != fscanf(f, "%u", &u) )
        throw InvalidIO("invalid Object header");
    id = (ObjectID)u;
    int c = in.get_char();
    if ( c == ':' )
    {
        if ( 1 != fscanf(f, "%u", &u) )
            throw InvalidIO("invalid Object header");
        mk = (ObjectMark)u;
        if ( (unsigned)mk != u )
            throw InvalidIO("overflow ObjectMark");
    }
    else
        in.unget_char(c);
}


/**
 Load one object from file
 
 If 'fat==true', read the larger Marker format
 */
void ObjectSet::loadObject(Inputter& in, const ObjectTag tag, int bin)
{
    PropertyID pid = 0;
    ObjectID id = 0;
    ObjectMark mk = 0;
    
    switch ( bin )
    {
        case 0: readMarkerASCII(in, pid, id, mk); break;
        case 1: readMarker(in, pid, id); break;
        default: readMarkerFat(in, pid, id, mk); break;
    }
    
#if BACKWARD_COMPATIBILITY < 45
    // before 18/09/2015, Property IDs started at zero:
    pid += ( in.formatID() < 45 );
#endif

    //std::clog << "- loading " << Object::make_reference(tag, pid, id) << '\n';

    Object * obj = identifyObject(id);
    
    /*
     A lowercase TAG indicates a new object, while uppercase describes some
     information that is associated to an already constructed object.
     LATTICE_TAG = 'l' for backward compatibility with format 56 (before 23/06/2021)
     */
#if BACKWARD_COMPATIBILITY <= 56
    bool primary = ( islower(tag) && tag != 'l' );
#else
    bool primary = islower(tag);
#endif
    
    if ( obj && primary )
    {
        assert_true(obj->property());
        /* check that property index has not changed: this can happen however, if
         objects are created/destroyed since the identity numbers are recycled */
        if ( obj->property()->number() != pid )
        {
#if 0
            Property const* P = obj->property();
            std::clog << "Remark: transmuting " << P->category() << P->number() << " `" << P->name();
            std::clog << "' to load object with property #" << pid << '\n';
#endif
            // the orphan Object remains on the 'ice_' to be deleted during pruning:
            inventory_.unassign(obj);
            obj->setIdentity(0);
            obj = nullptr;
        }
        else if ( obj->objset() )
            ice_.pop(obj);
        else
            obj->objset(this);
    }
    
    if ( !obj )
    {
        if ( id == 0 )
            throw InvalidIO("Invalid ObjectID referenced in file");
        assert_true(primary);
        assert_true(isprint(tag));
        //std::clog << "- new " << Object::reference(tag, pid, id) << '\n';
        // create new object of required class, identified by property-id
        obj = newObject(tag, pid);
        if ( !obj )
            throw InvalidIO("unknown Object referenced in file");
        obj->objset(this);
        obj->setIdentity(id);
        //inventory_.get(id);
        inventory_.assign(obj);
    }
    assert_true( obj->identity() == id );
    
    try {
        //std::clog << "- read " << obj->reference() << '\n';
        // read object data:
        obj->read(in, simul_, tag);
    }
    catch( Exception & e )
    {
        std::cerr << e.brief() << e.info() << " while loading " << obj->reference() << '\n';
        if ( primary ) link(obj);
        throw;
    }
    
    if ( primary )
    {
        if ( obj->invalid() )
        {
            inventory_.unassign(obj);
            delete(obj);
        }
        else
        {
            link(obj);
            if ( mk ) obj->mark(mk);
        }
    }
}


//------------------------------------------------------------------------------


void ObjectSet::writeReport(std::ostream& os, const std::string& title) const
{
    if ( size() > 0 )
    {
        os << '\n' << title;
        PropertyList plist = simul_.properties.find_all(title);
        if ( plist.size() > 0 )
        {
            for ( Property * p : plist )
            {
                size_t cnt = count(p);
                os << '\n' << std::setw(10) << cnt << " " << p->name();
            }
            if ( plist.size() > 1 )
                os << '\n' << std::setw(10) << size() << " total";
        }
        else
        {
            os << '\n' << std::setw(10) << size() << " " << title;
        }
    }
}

