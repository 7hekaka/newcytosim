// Cytosim was created by Francois Nedelec. Copyright 2024

#include "field_set.h"
#include "field_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "field.h"

//------------------------------------------------------------------------------

void FieldSet::prepare()
{
    for ( Field * f=first(); f; f=f->next() )
    {
        assert_true( f->hasField() );
        f->prepare();
    }
}


void FieldSet::steps()
{
    Field * obj = first();
    while ( obj )
    {
        Field * nxt = obj->next();
        if ( obj->hasField() )
        {
            LOG_ONCE("!!!! Field is active\n");
            obj->step(simul_.fibers);
        }
        obj = nxt;
    }
    if ( size() > 1 ) shuffle();
}

//------------------------------------------------------------------------------
#pragma mark -

Property* FieldSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "field" )
        return new FieldProp(nom);
    return nullptr;
}


Object * FieldSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Field::TAG )
    {
        FieldProp * p = simul_.findProperty<FieldProp>("field", pid);
        return new Field(p);
    }
    throw InvalidIO("Warning: unknown Field tag `"+std::to_string(tag)+"'");
    return nullptr;
}


/**
 @ingroup NewObject

 Specify the initial value of the Field:
 
     new field NAME
     {
        value = 0
     }
 
 \todo: read the value of the field from a file, at initialization
 */
ObjectList FieldSet::newObjects(Property const* p, Glossary& opt)
{
    FieldProp const* pp = static_cast<FieldProp const*>(p);
    Field * obj = new Field(pp);
        
    // initialize field:
    obj->setField();
        
    // an initial concentration can be specified:
    FieldCell val(0);
    if ( opt.set(val, "value", "initial_value") )
    {
        std::string str;
        if ( opt.set(str, "value", 1) )
        {
            Space const* spc = simul_.findSpace(str);
            if ( !spc )
                spc = obj->prop->field_space_ptr;
            obj->setConcentration(spc, val, 0);
        }
        else
        {
            obj->setConcentration(val);
        }
    }

    return ObjectList(obj);
}


void FieldSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}

