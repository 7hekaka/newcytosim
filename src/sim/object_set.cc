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

extern Modulo const* modulo;

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
    {
        switch ( obj->mobile() )
        {
            case 1: obj->rotateT(iso); obj->translate(iso); break;
            case 2: obj->rotate(iso); break;
            case 3: obj->rotate(iso); obj->translate(iso); break;
        }
    }
}


void ObjectSet::flagObjects(ObjectList const& objs, ObjectFlag f)
{
    for ( Object * obj : objs )
        obj->flag(f);
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
            switch ( obj->mobile() )
            {
                case 1: obj->rotateT(iso); obj->translate(iso); break;
                case 2: obj->rotate(iso); break;
                case 3: obj->rotate(iso); obj->translate(iso); break;
            }
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
        //std::clog << "ObjectSet::add(" << obj->reference() << ")\n";
    }
    else
    {
        std::cerr << "Warning: attempted to re-link "+obj->reference()+" \n";
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
    obj->objset(nullptr);
    unlink(obj);
    inventory_.unassign(obj);
}


void ObjectSet::remove(ObjectList const& list)
{
    for ( Object * obj : list )
        remove(obj);
}


void ObjectSet::erase(Object * obj)
{
    //std::clog << "ObjectSet::erase " << obj->reference() << '\n';
    remove(obj);
    delete(obj);
}


void ObjectSet::erasePool(ObjectPool & list)
{
    Object * i = list.pop_front();
    while ( i )
    {
        static_cast<Object*>(i)->objset(nullptr);
        delete(i);
        i = list.pop_front();
    }
}


void ObjectSet::erase()
{
    erasePool(pool_);
    inventory_.clear();
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
            while ( inv  &&  ++num <= 0 )
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
    
    // finally search for a property name:
    Property * pp = simul_.findProperty(cat, spec);

    if ( pp )
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
                if ( num > 0 )
                    break;
                inv = inventory_.previous(inv);
            }
        }
        return static_cast<Object*>(inv);
    }
    
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
        return findObject(cat, spec, num);
    
    // check category name, eg. 'fiber':
    if ( spec == cat )
    {
        ObjectList all = collect();
        //std::clog << "findObject -> random pick among " << sel.size() << " " << title << "\n";
        if ( all.size() > 0 )
            return all.pick_one();
    }
    
    // check property name:
    Property const* p = simul_.findProperty(cat, spec);
    if ( p )
    {
        ObjectList sel = collect(match_property, p);
        //std::clog << "findObject -> random pick among " << sel.size() << " " << spec << "\n";
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


size_t ObjectSet::count(const ObjectPool & list,
                        bool (*func)(Object const*, void const*), void const* arg)
{
    size_t res = 0;
    Object const* n = list.front();
    while ( n )
    {
        res += func(n, arg);
        n = n->next();
    }
    return res;
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
    ObjectList res(list.size(), 16);
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
    if ( cnt < objs.size() )
    {
        // limit the list to a random subset
        if ( cnt == 1 )
            objs[0] = objs.pick_one();
        else
            objs.shuffle();
        objs.truncate(cnt);
    }
    return objs;
}


ObjectList ObjectSet::collect(bool (*func)(Object const*, void const*), void const* arg) const
{
    return collect(pool_, func, arg);
}


ObjectList ObjectSet::collect(bool (*func)(Object const*, void const*), void const* arg, const size_t cnt) const
{
    ObjectList objs = collect(pool_, func, arg);
    if ( cnt < objs.size() )
    {
        // limit the list to a random subset
        if ( cnt == 1 )
            objs[0] = objs.pick_one();
        else
            objs.shuffle();
        objs.truncate(cnt);
    }
    return objs;
}


ObjectList ObjectSet::collect(Property const* p) const
{
    return collect(match_property, p);
}


size_t ObjectSet::count(bool (*func)(Object const*, void const*), void const* arg) const
{
    return count(pool_, func, arg);
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
    ice_.append(pool_);
}


void ObjectSet::prune()
{
    Object * i = ice_.pop_front();
    while ( i )
    {
        Object * o = i;
        i = ice_.pop_front();
        inventory_.unassign(o);
        o->objset(nullptr);
        delete(o);
    }
}


void ObjectSet::thaw()
{
    Object * i = ice_.pop_front();
    while ( i )
    {
        link(i);
        i = ice_.pop_front();
    }
}


/**
 Write Reference and Object's data, for all Objects in `list`
 */
void ObjectSet::writeObjects(Outputter& out, ObjectPool const& list) const
{
    // record number of objects:
    size_t cnt = list.size();
    size_t sup = inventory_.highest();
    out.write("\n#record "+std::to_string(cnt)+" "+std::to_string(sup));
    if ( out.binary() ) out.put_char('\n');
    
    // record objects:
    for ( Object const* n=list.front(); n; n=n->next() )
    {
        //std::clog << "write " << n->reference() << '\n';
        n->write(out);
    }
}


/** This should match Object::writeHeader() */
static void readObjectHeader(Inputter& in, bool fat, PropertyID& ix, ObjectID& id, ObjectMark& mk)
{
    if ( in.binary() )
    {
        // read header in binary format
        if ( fat )
        {
            ix = in.readUInt16();
            id = in.readUInt32();
            mk = in.readUInt32();
        }
        else
        {
            ix = in.readUInt8();
            id = in.readUInt16();
        }
    }
    else
    {
        // read header in text format
        FILE * f = in.file();
        if ( 1 != fscanf(f, "%u", &ix) )
            throw InvalidIO("invalid Object header");
        if ( in.get_char() != ':' )
            throw InvalidIO("invalid Object header");
        if ( 1 != fscanf(f, "%u", &id) )
            throw InvalidIO("invalid Object header");
        int c = in.get_char();
        if ( c == ':' )
        {
            if ( 1 != fscanf(f, "%u", &mk) )
            throw InvalidIO("invalid Object header");
        }
        else
            in.unget(c);
    }
}


/**
 Load one object from file
 
 If 'fat==true', read the larger header format
 If 'update==true', the corresponding object is changed, or a new object is created.
 */
void ObjectSet::loadObject(Inputter& in, const ObjectTag tag, bool fat, bool update)
{
    PropertyID pid = 0;
    ObjectID id = 0;
    ObjectMark mk = 0;
    Object * obj = nullptr;
    
    readObjectHeader(in, fat, pid, id, mk);
    
    if ( id == 0 )
        throw InvalidIO("Invalid ObjectID referenced in file");
    //std::clog << "- load " << Object::reference(tag, pid, id) << '\n';

    if ( update )
    {
        obj = findID(id);
        /*
         only 'primary' objects with a lowercase TAG are linked, but we exclude the old
         TAG_LATTICE = 'l' for backward compatibility with format 56 (before 23/06/2021)
         */
        if ( obj && islower(tag) && tag != 'l' )
        {
            ice_.pop(obj);
#if ( 1 )
            // check that property index has not changed:
            if ( obj->property()->number() != pid )
            {
                Property const* P = obj->property();
                std::clog << "Warning: erasing " << P->category() << P->number() << " `" << P->name();
                std::clog << "' to load object with property #" << pid << '\n';
                inventory_.unassign(obj);
                obj->objset(nullptr);
                delete(obj);
                obj = nullptr;
            }
#endif
        }
        else
            update = false;
    }
    
    if ( !obj )
    {
        assert_true(isprint(tag));
        // create new object of required class, indentified by property-id
        obj = newObject(tag, pid);
        assert_true(obj);
        update = true;
        obj->setIdentity(id);
        obj->objset(this);
        inventory_.assign(obj);
        //std::clog << "- new " << Object::reference(tag, pid, id) << '\n';
    }
    assert_true( obj->identity() == id );
    assert_true( obj->property() );
    
    try {
        //std::clog << "- read " << Object::reference(tag, pid, id) << '\n';
        // read object data:
        obj->read(in, simul_, tag);
    }
    catch( Exception & e )
    {
        e << "while loading " << Object::reference(tag, pid, id);
        throw;
    }
    
    if ( update )
    {
        link(obj);
        if ( mk ) obj->mark(mk);
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
                size_t cnt = count(match_property, p);
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

