// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "object.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "property.h"
#include "object_set.h"
#include "cymdef.h"


Object::~Object()
{
    //std::clog << "~Object " << this << '\n';

    if ( set_ )
    {
        //std::clog << "unlinking deleted Object " << this << '\n';
        set_->remove(this);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 This is an accessory function used in ObjectSet::collect()
 */
bool match_all(Object const*, void const*)
{
    return true;
}

/**
 This is an accessory function used in ObjectSet::collect()
 */
bool match_mark(Object const* obj, void const* mrk)
{
    return ( obj->mark() == *((ObjectMark*)mrk) );
}

/**
 This is an accessory function used in ObjectSet::collect()
 */
bool match_property(Object const* obj, void const* prop)
{
    return ( obj->property() == prop );
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The ASCII reference has the format `XP:N` where:
 - X is one ascii character, indicating the category (tag),
 - P > 0 is the index of the property, indicating the type of object,
 - N > 0 is the ID of the object within its category.
 .

 For example 'f0:021' is the fiber number 21 of class 0 (ie. property index = 0),
 For example 'f1:010' is the fiber number 10 of class 1 (ie. property index = 1)
*/
std::string Object::reference(ObjectTag tag, size_t pip, ObjectID id)
{
    assert_true( pip > 0 );
    assert_true( id > 0 );
    char tmp[32];
    snprintf(tmp, sizeof(tmp), "%c%lu:%04u", tag, pip, id);
    return std::string(tmp);
}


/**
 The ASCII reference has the format `XP:N` where:
 - X is one ascii character, indicating the category (ObjectSet::tag()),
 - P > 0 is the index of the property, indicating the type of object,
 - N > 0 is the ID of the object within its category.
 .
 
 For example 'f0:21' is the fiber of property 0, number 21
 For example 'f1:10:1' is the fiber of property 1, number 10, and it is marked as '1'
 */
std::string Object::reference() const
{
    if ( property() )
        return reference(tag(), property()->number(), identity());
    else
        return reference(tag(), 0, identity());
}


/**
 Two binary formats are used:
 - A slim format:
     - 1 byte for the tag()
     - 2 bytes for the identity
     .
 - A fat format:
     - 1 byte containing tag() + 128
     - 4 bytes for the identity
     .
 .
 The ascii-based format is always the same.
 This data is read by Simul::readReference()
 */
void Object::writeReference(Outputter& out, ObjectTag g, ObjectID id)
{
    assert_true( id > 0 );

    if ( out.binary() )
    {
#if 1
        // slim/fat format 50 -- 56 using different size
        if ( id > 65535 )
        {
            // set the highest bit of the byte, which is not used by ASCII codes
            out.writeChar(g|HIGH_BIT);
            out.writeUInt32(id);
        }
        else
        {
            out.writeChar(g);
            out.writeUInt16(id);
        }
#else
        /*
         format tryout on 11.06.2021:
         combine `tag` and 'id', leaving 3 bytes and at most 16777216 objects
         Note that the topmost bit of ASCII is not used
         */
        out.writeUInt32(id|(uint32_t(g)<<24));
#endif
    }
    else
    {
        out.writeChar(' ');
        out.writeUInt(id, g);
    }
}


void Object::writeNullReference(Outputter& out)
{
    if ( out.binary() )
    {
#if 1
        out.writeChar(NULL_TAG);
#else
        /*
         format tryout on 11.06.2021:
         combining `tag` and 'id', leaving 3 bytes and at most 16777216 objects
         Note that the topmost bit of ASCII is not used
         */
        out.writeUInt32(uint32_t(NULL_TAG)<<24);
#endif
    }
    else
    {
        out.writeChar(' ');
        out.writeChar(NULL_TAG);
    }
}


void Object::writeReference(Outputter& out, Object const* i)
{
    if ( i )
        writeReference(out, i->tag(), i->identity());
    else
        writeNullReference(out);
}


/**
Writes the info that is common to all objects to file
 Two binary formats are used:
 - A slim format:
     - 1 byte for the tag()
     - 1 byte for the index of the property
     - 2 bytes for the identity
     .
 - A fat format:
     - 1 byte for the tag() with the highest bit set
     - 2 bytes for the index of the property
     - 4 bytes for the identity
     - 4 bytes for the mark
     .
 .
 The ascii based format is invariant.
 This data is read back in `readObjectHeader()`
 */
void Object::writeHeader(Outputter& out, ObjectTag g) const
{
    if ( out.binary() )
    {
        if ( identity() > 16777216 )
            throw InvalidIO("binary file data format overflow");
        if ( identity() > 65535 || property()->number() > 255 || mark() )
        {
            // set the highest bit of the byte, which is not used by ASCII codes
            out.writeChar(g|HIGH_BIT);
            out.writeUInt16(property()->number());
            out.writeUInt32(identity());
            out.writeUInt32(mark());
        }
        else
        {
            out.writeChar(g);
            out.writeUInt8(property()->number());
            out.writeUInt16(identity());
        }
    }
    else
    {
        out.writeChar('\n');
        out.writeChar(g);
        out.writeUInt(property()->number());
        out.writeUInt(identity(), ':');
        if ( mark() )
            out.writeUInt(mark(), ':');
    }
}



/// print a list of objects
std::ostream& operator << (std::ostream& os, ObjectList const& arg)
{
    os << "ObjectList " << &arg << " of size " << arg.size() << "\n{\n";
    for ( Object * obj : arg )
        os << "   " << obj->property()->name() << " " << obj->reference() << '\n';
    os << "}" << '\n';
    return os;
}

