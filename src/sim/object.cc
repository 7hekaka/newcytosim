// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "object.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "property.h"
#include "object_set.h"
#include "cymdef.h"


Object::~Object()
{
    //std::clog << "~Object " << this << "  " << reference() << '\n';

    if ( set_ )
    {
        //std::clog << "unlinking deleted Object " << this << '\n';
        //set_->remove(this);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 This is a selection function used in ObjectSet::collect()
 */
bool match_all(Object const*, void const*)
{
    return true;
}

/**
 This is a selection function used in ObjectSet::collect()
 */
bool match_mark(Object const* obj, void const* mrk)
{
    return ( obj->mark() == reinterpret_cast<uintptr_t>(mrk) );
}

/**
 This is a selection function used in ObjectSet::collect()
 */
bool match_property(Object const* obj, void const* arg)
{
    return ( arg == nullptr ) || ( obj->property() == arg );
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
std::string Object::make_reference(ObjectTag tag, unsigned pip, ObjectID id)
{
    char tmp[32];
    snprintf(tmp, sizeof(tmp), "%c%u:%04u", tag, pip, id);
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
        return make_reference(tag(), property()->number(), identity());
    else
        return make_reference(tag(), 0, identity());
}


/**
 Two binary formats are used:
 - A slim format:
     - 1 byte for the tag()
     - 2 bytes for the identity
     .
 - A fat format:
     - 1 byte containing tag() + HIGH_BIT
     - 3 bytes for the identity
     .
 .
 The ascii-based format is always the same.
 This data is read by Simul::readReference()
 */
void Object::writeReference(Outputter& out, ObjectTag g, ObjectID id)
{
    assert_true(isalpha(g) || g==NULL_TAG);

    if ( out.binary() )
    {
        if ( id > 65535 )
        {
            /* The fat format is signaled by the highest bit of the byte,
             which is not used by ASCII codes. */
            if ( 1 )
            {
                /* format 58 enabled on 26.11.2022: combining `tag` and 'id',
                 leaving 3 bytes and at most 16 777 216 objects */
                if ( id > 1<<24 ) // ~16 Million objects
                    throw InvalidIO("binary file data format overflow");
                out.writeChar(g|HIGH_BIT);
                out.writeChar((id>>16)&0xFF);
                out.writeChar((id>>8)&0xFF);
                out.writeChar(id&0xFF);
            }
            else
            {
                out.writeChar(g|HIGH_BIT);
                out.writeUInt32(id);
            }
        }
        else
        {
            // slim format
            out.writeChar(g);
            if ( g != NULL_TAG )
                out.writeUInt16(id);
        }
    }
    else
    {
        out.writeChar(' ');
        if ( g != NULL_TAG )
            out.writeUInt(id, g);
        else
            out.writeChar(g);
    }
}


void Object::writeReference(Outputter& out, Object const* i)
{
    if ( i )
        writeReference(out, i->tag(), i->identity());
    else
        writeReference(out, NULL_TAG, 0);
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
     - 1 byte for the mark
     - 2 bytes for the index of the property
     - 4 bytes for the identity
     .
 .
 The ascii based format is invariant.
 This data is read back in `readMarker()`
 */
void Object::writeMarker(Outputter& out, ObjectTag g) const
{
    if ( !out.binary() )
    {
        out.writeChar('\n');
        out.writeChar(g);
        out.writeUInt(property()->number());
        out.writeUInt(identity(), ':');
        if ( mark() )
            out.writeUInt(mark(), ':');
    }
    else
    {
        if ( identity() > 65535 || property()->number() > 255 || mark() )
        {
            // using 8 bytes
            // set the highest bit, which is not used by ASCII codes
            out.writeChar(g|HIGH_BIT);
            out.writeUInt8(255&mark());
            out.writeUInt16(property()->number());
            out.writeUInt32(identity());
            /* Prior 26.11.2022, in format 56, this was using 11 bytes:
             out.writeChar(g|HIGH_BIT);
             out.writeUInt16(property()->number());
             out.writeUInt32(identity());
             out.writeUInt32(mark());
             */
        }
        else
        {
            // using 4 bytes
            out.writeChar(g);
            out.writeUInt8(property()->number());
            out.writeUInt16(identity());
        }
    }
}



/// print a list of objects
std::ostream& operator << (std::ostream& os, ObjectList const& arg)
{
    if ( arg.size() == 0 )
        os << "ObjectList " << &arg << " is empty\n{\n";
    else if ( arg.size() == 1 )
    {
        Object * obj = arg[0];
        os << "ObjectList " << &arg << " contains " << obj->property()->name() << " " << obj->reference() << '\n';
    }
    else
    {
        os << "ObjectList " << &arg << " has " << arg.size() << " objects:\n{\n";
        for ( Object * obj : arg )
            os << "   " << obj->property()->name() << " " << obj->reference() << '\n';
        os << "}" << '\n';
    }
    return os;
}

